// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaString/FFaTokenizer.H"
#include "FFaLib/FFaOS/FFaFortran.H"
#include <cstdlib>

#if defined(win32) || defined(win64)
#include <windows.h>
#define DLLexport(ret) extern "C" __declspec(dllexport) ret
#else
#define DLLexport(ret) extern "C" ret
#endif

SUBROUTINE (modes,MODES) (const int& nStep, const int& nModes, double* tSteps,
                          int* modeNr, int* processMode, int& iret);
SUBROUTINE (stress,STRESS) (int& iret);
SUBROUTINE (gage,GAGE) (int& iret);
SUBROUTINE (fpp,FPP) (int& iret);

extern void cmdLineArgInit (int argc, char** argv);
extern void readOptionFiles (const char* program, int whichFile = 1);


static bool getModesToExpand (std::vector<double>& tSteps,
                              std::vector<int>& modeNumbs,
                              std::vector<int>& processThisMode)
{
  std::string recoverModes;
  FFaCmdLineArg::instance()->getValue ("recover_modes",recoverModes);
  FFaTokenizer timeTokens(recoverModes,'<','>',',');
  const size_t nStep = timeTokens.size();
  if (nStep < 1)
  {
    std::cerr <<"fedem_modes: No time steps, check option -recover_modes \""
              << recoverModes <<"\"\n";
    return false;
  }

  tSteps.clear();
  tSteps.insert(tSteps.end(),nStep,0.0);
  modeNumbs.clear();
  processThisMode.clear();

  // For each time step, get list of mode numbers to expand
  size_t i, j, k;
  for (i = 0; i < nStep; i++)
  {
    FFaTokenizer modeTokens(timeTokens[i],'<','>',',');
    if (modeTokens.size() < 2)
    {
      std::cerr <<"fedem_modes: No modes at time step "<< i
                <<" check sub-token \""<< timeTokens[i]
                <<"\" in -recover_modes\n";
      return false;
    }
    tSteps[i] = atof(modeTokens.front().c_str());
    for (j = 1; j < modeTokens.size(); j++)
    {
      // Check if this mode number has been specified earlier
      int jMode = atoi(modeTokens[j].c_str());
      for (k = 0; k < modeNumbs.size(); k++)
        if (modeNumbs[k] == jMode) break;

      if (k == modeNumbs.size())
      {
        // This is a new mode number, add it to the list of modes
        modeNumbs.push_back(jMode);
        processThisMode.insert(processThisMode.end(),nStep,0);
      }

      // Flag that mode modeNumbs[k] shall be process at the i'th time step
      processThisMode[nStep*k+i] = 1;
    }
  }

  return true;
}


DLLexport(void) initSolverArgs (int argc, char** argv)
{
  cmdLineArgInit (argc,argv);
}


DLLexport(int) solveModes ()
{
  int ierr = 1;
  readOptionFiles ("fedem_modes");
  std::vector<double> timeSteps;
  std::vector<int> modeNums, processThisMode;
  if (getModesToExpand(timeSteps,modeNums,processThisMode))
    F90_NAME(modes,MODES) ((int)timeSteps.size(),(int)modeNums.size(),
                           timeSteps.data(),modeNums.data(),
                           processThisMode.data(),ierr);
  return ierr;
}


DLLexport(int) solveStress ()
{
  int ierr = 0;
  readOptionFiles ("fedem_stress");
  F90_NAME(stress,STRESS) (ierr);
  return ierr;
}


DLLexport(int) solveGage ()
{
  int ierr = 0;
  readOptionFiles ("fedem_gage");
  F90_NAME(gage,GAGE) (ierr);
  return ierr;
}


DLLexport(int) solveFpp ()
{
  int ierr = 0;
  readOptionFiles ("fedem_fpp");
  F90_NAME(fpp,FPP) (ierr);
  return ierr;
}
