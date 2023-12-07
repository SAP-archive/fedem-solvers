// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifdef FT_HAS_HWAFLS
#include "win_fedem_protocol.h"
#include "win_HWAFLS_Queue.h"
#endif

#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>

#include "FFaLib/FFaOS/FFaFortran.H"

#ifdef FT_HAS_HWAFLS
static HWAFLS_Queue* ourQueue = 0;
static HWAFLS_data*  ourHdata = 0;
#else
static int M = 0;
static int N = 0;
static std::vector<double> wdata;
static std::vector<double> wkin;
#endif
static double depth = 0.0;
static bool working = false;

#ifdef FT_USE_PROFILER
extern "C" float _cdecl CLKSEC (const float* t);
static float   wHWA = 0.0;
static float   wMap = 0.0;
static clock_t tCPU = 0;
static clock_t tMap = 0;
static size_t  nCPU = 0;
static size_t  nMap = 0;
#endif

static std::vector<double>  XZ;
static std::map<int,size_t> nodeMap;


#ifndef FT_HAS_HWAFLS
//! Software implementation used for verification.
//! This is exactly how it is implemented on the HWAFLS card.
void finiteDepthWave (double x, double z, double t,
                      double& h, double& u, double& w, double& du, double& dw)
{
  h = u = w = du = dw = 0.0;

  const double pi = acos(-1.0);
  const double zz = z < 0.0 ? depth+z : depth; // No wheeler stretching
  double A, omega, eps, k, kd, kz, ch, sh, s, c;

  for (int i = 0; i < N; i++)
  {
    A     = wdata[  4*i];
    omega = wdata[1+4*i];
    eps   = wdata[2+4*i];
    k     = wdata[3+4*i];

    kz = k*zz;
    kd = k*depth;
    s  = sin(omega*t+eps-k*x);
    c  = cos(omega*t+eps-k*x);

    if (kd > pi && kz > pi)
      sh = ch = exp(kz-kd);
    else if (kd > pi && kz < pi && kz-kd < -13.0)
      sh = ch = 0.0;
    else
    {
      ch = cosh(kz)/sinh(kd);
      sh = sinh(kz)/sinh(kd);
    }

    h  += A*s;
    u  += A*omega*s*ch;
    w  += A*omega*c*sh;
    du += A*omega*omega*c*ch;
    dw -= A*omega*omega*s*sh;
  }
}
#endif


SUBROUTINE(hwafls_init,HWAFLS_INIT) (const int& nNod, const int& nComp,
                                     const double& D, double* RFUNC)
{
  std::cout <<"HWAFLS_init: nNod = "<< nNod <<", nComp = "<< nComp
            <<", D = "<< D << std::endl;

#ifdef FT_HAS_HWAFLS
  if (ourQueue) return;

  ourQueue = new HWAFLS_Queue();
  ourHdata = new HWAFLS_data(nNod,nComp);

  ourHdata->set_d(depth=D);
  ourHdata->set_A_omega_eps_k(RFUNC);
#else
  M = nNod;
  N = nComp;
  depth = D;
  wdata.resize(4*nComp);
  memcpy(&wdata.front(),RFUNC,4*nComp*sizeof(double));
  wkin.resize(5*nNod);
#endif

  XZ.resize(2*nNod,0.0);

#ifdef FT_USE_PROFILER
  wHWA = wMap = 0.0;
  tCPU = tMap = 0;
  nCPU = nMap = 0;
#endif
}


SUBROUTINE(hwafls_close,HWAFLS_CLOSE) ()
{
#ifdef FT_HAS_HWAFLS
  if (!ourQueue) return;

  ourQueue->closeQ();

  delete ourHdata;
  delete ourQueue;
#else
  wdata.clear();
#endif

#ifdef FT_USE_PROFILER
  if (tCPU > 0)
    std::cout <<"Total time in HWAFLS: "<< tCPU <<" [clocks] "
              << (double)tCPU/CLOCKS_PER_SEC <<" [s], wall="
              << wHWA <<" [s] ("
              << nCPU <<" instances)"<< std::endl;

  if (tMap > 0)
    std::cout <<"Total time in nodemap lookup: "<< tMap <<" [clocks] "
              << (double)tMap/CLOCKS_PER_SEC <<" [s], wall="
              << wMap <<" [s] ("
              << nMap <<" instances)"<< std::endl;
#endif

  XZ.clear();
  nodeMap.clear();
}


SUBROUTINE(hwafls_setnode,HWAFLS_SETNODE) (const int& inod,
					   const double& X, const double& Z)
{
#ifdef FT_USE_PROFILER
  const float zero(0.0);
  clock_t start(clock());
  float tc = CLKSEC(&zero);
#endif

  size_t ip = nodeMap.size();
  std::map<int,size_t>::const_iterator nit = nodeMap.find(inod);
  if (nit == nodeMap.end())
    nodeMap[inod] = ip;
  else
    ip = nit->second;

  if (2*ip > XZ.size())
    std::cerr <<" *** HWAFLS_setNode: Node index "<< ip <<" out of range [0,"
              << XZ.size()/2 <<"] inod="<< inod << std::endl;
  else
  {
    XZ[2*ip]   = X;
    XZ[2*ip+1] = Z;
  }

#ifdef FT_USE_PROFILER
  wMap += CLKSEC(&tc);
  clock_t stop(clock());
  tMap += stop-start;
  nMap ++;
#endif
}


SUBROUTINE(hwafls_evaluate,HWAFLS_EVALUATE) (const double& t)
{
#ifdef FT_USE_PROFILER
  const float zero(0.0);
  clock_t start(clock());
  float tc = CLKSEC(&zero);
#endif
#ifdef FT_HAS_HWAFLS
  ourHdata->set_t(t);
  ourHdata->set_x_z(&XZ.front());
  ourQueue->send(ourHdata);
#else
  // Software emulation of the HWAFLS module (for verification)
  const double* p = &XZ.front();
  double* q = &wkin.front();
  for (int i = 0; i < M; i++, p += 2, q +=5)
    finiteDepthWave(*p,*(p+1),t,*q,*(q+1),*(q+2),*(q+3),*(q+4));
#endif
#ifdef FT_USE_PROFILER
  wHWA += CLKSEC(&tc);
  clock_t stop(clock());
  tCPU += stop-start;
  nCPU ++;
#endif

  working = true;
}


SUBROUTINE(hwafls_getnode,HWAFLS_GETNODE) (const int& inod, double* waveData)
{
  if (working)
  {
#ifdef FT_USE_PROFILER
    const float zero(0.0);
    clock_t start(clock());
    float tc = CLKSEC(&zero);
#endif
#ifdef FT_HAS_HWAFLS
    ourQueue->get(ourHdata);
#endif
#ifdef FT_USE_PROFILER
    wHWA += CLKSEC(&tc);
    clock_t stop(clock());
    tCPU += stop-start;
    nCPU ++;
#endif
  }
  working = false;
  memset(waveData,0,9*sizeof(double));

#ifdef FT_USE_PROFILER
  const float zero(0.0);
  float tc = CLKSEC(&zero);
  clock_t start(clock());
#endif
  std::map<int,size_t>::const_iterator nit = nodeMap.find(inod);
  if (nit == nodeMap.end())
    std::cerr <<" *** HWAFLS_getNode: Node "<< inod <<" not found."<< std::endl;
  else
  {
    double Zp = XZ[2*nit->second+1];
    size_t ip = 5*nit->second;
#ifdef FT_HAS_HWAFLS
    waveData[2] = ourHdata->returndata[ip++];
#else
    waveData[2] = wkin[ip++];
#endif
    if (Zp >= -depth && Zp <= waveData[2]) // check if point is in the water
    {
#ifdef FT_HAS_HWAFLS
      waveData[3] = ourHdata->returndata[ip++];
      waveData[5] = ourHdata->returndata[ip++];
      waveData[6] = ourHdata->returndata[ip++];
      waveData[8] = ourHdata->returndata[ip];
#else
      waveData[3] = wkin[ip++];
      waveData[5] = wkin[ip++];
      waveData[6] = wkin[ip++];
      waveData[8] = wkin[ip];
#endif
    }
  }
#ifdef FT_USE_PROFILER
  wMap += CLKSEC(&tc);
  clock_t stop(clock());
  tMap += stop-start;
  nMap ++;
#endif
}
