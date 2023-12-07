// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testShells.C

  \brief Convergence tests for shell elements.

  \details This file contains a main program for executing a set of convergence
  tests of the shell elements on some standard benchmark test models.

  \author Knut Morten Okstad, SAP SE

  \date 5 Jul 2021
*/

#include <array>
#include <fstream>
#include <cstring>

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFilePath.H"
#include "FFaLib/FFaOS/FFaFortran.H"
#include "Admin/FedemAdmin.H"

//! \cond DO_NOT_DOCUMENT
int createFEModel (int iPart, int nels, int nel2,
                   double = 0.0, double = 0.0, double = 0.0,
                   bool = false, bool solve = true);
SUBROUTINE (solver,SOLVER) (const int&, const int*, double*, int& iret);
SUBROUTINE (setdate,SETDATE) (const char*, const int nchar = 0);
SUBROUTINE (setversion,SETVERSION) (const char*, const int nchar = 0);
#define ADDOPTION FFaCmdLineArg::instance()->addOption
//! \endcond


/*!
  \brief Main program for execution of convergence tests.
*/

int main (int argc, char** argv)
{
  // Initialize command-line parser
  FFaCmdLineArg::init (argc,argv);
  // Only define command-line options relevant for this test program
  ADDOPTION ("help",false,"Print out this help text");
  ADDOPTION ("version",false,"Print out program version");
  ADDOPTION ("terminal",6,"File unit number for terminal output");
  ADDOPTION ("debug",0,"Debug print switch");
  ADDOPTION ("printArray",0,"Debug print switch for enumerated arrays");
  ADDOPTION ("linkId",1,"Link base-ID number");
  ADDOPTION ("ftlout","","Name of link output file in FTL-format");
  ADDOPTION ("samfile","","Name of SAM data file");
  ADDOPTION ("extNodes","","List of external nodes");
  ADDOPTION ("resfile","","Name of result output file");
  ADDOPTION ("frsfile","","Name of results database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database files");
  ADDOPTION ("lumpedmass",false,"Use lumped element mass matrices");
  ADDOPTION ("denseSolver",false,"Use the LAPACK equation solver");
  ADDOPTION ("skylineSolver",false,"Use the skyline equation solver");
  ADDOPTION ("tolFactorize",1.0e-12,"Equation solver singularity criterion");
  ADDOPTION ("useANDESformulation",true,"Shell formulation option");
  ADDOPTION ("useRotCpl",true,"Use rotation coupling in WAVGM's");
  ADDOPTION ("useEccCpl",true,"Use eccentricity coupling in WAVGM's");
  ADDOPTION ("autoStiffMethod",3,"Method for automatic stiffness computations");
  ADDOPTION ("autoStiffScale",1.0e2,"Scale factor for auto-added springs");
  ADDOPTION ("autoMassScale",1.0e-9,"Scale factor for auto-added masses");
  ADDOPTION ("drillingStiff",1.0e-6,"Fictitious drilling DOF stiffness");
  ADDOPTION ("useParabolic",false,"Use parabolic shell elements");
  ADDOPTION ("useTriangles",false,"Use triangular shell elements");
  ADDOPTION ("gvec",DoubleVec(),"Gravity vector in global coordinates");
  ADDOPTION ("n1min",1,"Min number of elements in each/first direction");
  ADDOPTION ("n1max",64,"Max number of elements in each/first direction");
  ADDOPTION ("case","PDiCyl","Benchmark identifier");
  ADDOPTION ("mute",true,"Suppress warnings on undefined options");
  FFaCmdLineArg::instance()->getValue("mute",FFaCmdLineArg::mute);

  // Get the program version and build date
  const char* my_version = FedemAdmin::getVersion();
  const char* build_date = FedemAdmin::getBuildDate();

  if (FFaCmdLineArg::instance()->isOptionSetOnCmdLine("help"))
  {
    // Print the help text and quit
    std::string hlpText;
    FFaCmdLineArg::instance()->composeHelpText (hlpText);
    std::cout <<"\n"<< argv[0] <<" version "<< my_version <<" "<< build_date;
    std::cout <<"\n\nAvailable command-line options:\n"<< hlpText;
    return 0;
  }
  if (FFaCmdLineArg::instance()->isOptionSetOnCmdLine("version"))
  {
    // Print program version and quit
    std::cout << argv[0] <<" version "<< my_version <<" "<< build_date <<"\n";
    return 0;
  }

  std::string caseID, ftlOut;
  int  n1min = 0, n1max = 0, ipart = 0;
  bool parabolic = false;
  FFaCmdLineArg::instance()->getValue("n1min",n1min);
  FFaCmdLineArg::instance()->getValue("n1max",n1max);
  FFaCmdLineArg::instance()->getValue("useParabolic",parabolic);
  FFaCmdLineArg::instance()->getValue("case",caseID);
  if (caseID == "PDiCyl")
    ipart = 8;
  else if (caseID == "OPHSph")
    ipart = 9;
  else if (caseID == "ScoLo")
    ipart = 10;
  else
  {
    std::cerr <<" *** Unknown benchmark identifier "<< caseID << std::endl;
    return 1;
  }

  FFaCmdLineArg::instance()->getValue("ftlout",ftlOut);
  ftlOut = FFaFilePath::getBaseName(ftlOut);

  // Pass the version tag and date to the fortran modules
  F90_NAME(setversion,SETVERSION) (my_version,(int)strlen(my_version));
  F90_NAME(setdate,SETDATE) (build_date,(int)strlen(build_date));

  int ierr = 0;
  std::string suffix("_1.ftl");
  std::map<int,std::array<double,2>> conv;
  if (n1min < 2) n1min = parabolic ? 1 : 2;
  for (int n1 = n1min; n1 <= n1max; n1 *= 2, ++suffix[1])
  {
    // Set output filename to use for this mesh
    if (!ftlOut.empty() && n1min*2 <= n1max)
      FFaCmdLineArg::instance()->setValue("ftlout",ftlOut+suffix);

    // Create the FE mesh with n1xn1 elements
    if ((ierr = createFEModel(ipart,n1,n1)))
      break;

    // Find two nodes to extract displacements at
    int node1 = parabolic ? 2*n1+1 : n1+1;
    std::array<int,2> nodes;
    if (ipart == 9)
      nodes = { 1, node1 };
    else
      nodes = { node1*node1, node1*(node1-1)+1 };

    // Solve the FE problem
    double dis[6];
    F90_NAME(solver,SOLVER) (2,nodes.data(),dis,ierr);
    if (ierr) break;

    std::cout <<"\n   * Deflection at nodes "<< nodes[0] <<" "<< nodes[1] <<":";
    for (int i = 0; i < 6; i++)
      std::cout <<" "<< dis[i];
    std::cout << std::endl;

    switch (ipart) {
    case  8: conv[node1*node1] = { dis[1], dis[3] }; break;
    case  9: conv[node1*node1] = { dis[0], dis[4] }; break;
    case 10: conv[node1*node1] = { dis[1], dis[4] }; break;
    }
  }
  if (conv.size() < 2)
    return ierr;

  std::ofstream os(caseID+".dat");
  os <<"# nnod displacements...";
  for (const std::pair<const int,std::array<double,2>>& c : conv)
    os <<"\n"<< c.first <<"\t"<< c.second[0] <<"\t"<< c.second[1];
  os << std::endl;
  return ierr;
}
