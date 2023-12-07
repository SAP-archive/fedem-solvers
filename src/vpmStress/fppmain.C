// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "vpmCommon/FFCmdLineArgInterface.H"

extern "C" {
  void initSolverArgs(int,char**);
  int  solveFpp();
}


int main (int argc, char** argv)
{
  // Initialize the command-line parser
  initSolverArgs (argc,argv);

  // Define all command-line options here
  ADDOPTION ("Bramsize",-1,"In-core size (MB) of displacement recovery matrix"
             "\n< 0: Use the same as in the reducer (default)"
             "\n= 0: Store full matrix in core");
  ADD_PRIVATE_OPTION ("dmramsize",-1,"Same as -Bramsize"
                      " but in terms of double words");
  ADDOPTION ("linkId",0,"Link base-ID number");
  ADDOPTION ("linkfile","","Name of link input file");
  ADDOPTION ("Bmatfile","","Name of B-matrix file");
  ADDOPTION ("eigfile","","Name of eigenvector file");
  ADDOPTION ("dispfile","","Name of gravitation displacement file");
  ADDOPTION ("resfile","","Name of result output file");
  ADDOPTION ("samfile","","Name of SAM data file");
  ADDOPTION ("fsifile","fedem_solver.fsi","Name of solver input file");
  ADDOPTION ("resStressFile","","Name of residual stress input file");
  ADDOPTION ("resStressSet","","Name of residual stress set");
  ADDOPTION ("frsfile","","Name of solver results database file");
  ADDOPTION ("fppfile","","Name of fpp output file");
  ADDOPTION ("rdbfile","","Name of strain coat results database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database file");
  ADDOPTION ("double",false,"Save results in double precision");
  ADD_PRIVATE_OPTION ("writeHistory",false,"Write history frs-files instead");
  ADD_PRIVATE_OPTION ("oldRange",false,"Use old stress/strain range meassures");
  ADDOPTION ("group","","List of element groups to do calculations for");
  ADDOPTION ("blockSize",2000,"Max number of elements processed together");
  ADDOPTION ("BufSizeInc",20,"Buffer increment size");
  ADDOPTION ("PVXGate",10.0f,"Gate value for the Peak Valley extraction"
             "\n(MPa or microns depending on HistDataType)");
  ADDOPTION ("biAxialGate",10.0,"Gate value for the biaxiality calculation");
  ADDOPTION ("angleBins",541,"Number of bins in search for most popular angle");
  ADDOPTION ("HistXMin",-100.0f,"Histogram min X-value");
  ADDOPTION ("HistXMax", 100.0f,"Histogram max X-value");
  ADDOPTION ("HistYMin",-100.0f,"Histogram min Y-value");
  ADDOPTION ("HistYMax", 100.0f,"Histogram max Y-value");
  ADDOPTION ("HistXBins",64,"Histogram number of X-bins");
  ADDOPTION ("HistYBins",64,"Histogram number of Y-bins");
  ADDOPTION ("HistDataType",0,"Histogram data type\n= 0: None"
             "\n= 1: Signed abs max stress\n= 2: Signed abs max strain");
  ADDOPTION ("surface",0,"Surface selection option"
             "\n= 0: All element surfaces"
             "\n= 1: Bottom shell surfaces only"
             "\n= 2: Middle shell surfaces only"
             "\n= 3: Top shell surfaces only");
  ADDOPTION ("SNfile","","Name of SN-curve definition file");
  ADDOPTION ("stressToMPaScale",1.0e-6,"Stress convertion factor to MPa");
  ADDOPTION ("statm",0.0,"Start time");
  ADDOPTION ("stotm",1.0,"Stop time");
  ADDOPTION ("tinc",0.0,"Time increment (= 0.0: process all time steps)");

  // Now launch fedem_fpp
  return solveFpp();
}
