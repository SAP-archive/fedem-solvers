// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "vpmCommon/FFCmdLineArgInterface.H"

extern "C" {
  void initSolverArgs(int,char**);
  int  solveGage();
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
  ADDOPTION ("resfile","","Name of results output file");
  ADDOPTION ("rdbfile","","Name of strain gage results database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database file");
  ADDOPTION ("samfile","","Name of SAM data file");
  ADDOPTION ("fsifile","fedem_solver.fsi","Name of solver input file");
  ADDOPTION ("frsfile","","Name of solver results database file");
  ADDOPTION ("rosfile","","Name of strain rosette input file");
  ADDOPTION ("writeAsciiFiles",false,"Write rosette results to ASCII files");
  ADDOPTION ("deformation",false,"Save nodal deformations to results database");
  ADDOPTION ("nullify_start_rosettestrains",false,
             "Set start strains to zero for the rosettes");
  ADDOPTION ("statm",0.0,"Start time");
  ADDOPTION ("stotm",1.0,"Stop time");
  ADDOPTION ("tinc",0.0,"Time increment (= 0.0: process all time steps)");
  ADDOPTION ("dac_sampleinc",0.001,"Sampling increment for dac output files");
  ADDOPTION ("flushinc",-1.0,"Time between each database file flush"
             "\n< 0.0: Do not flush results database (let the OS decide)"
             "\n= 0.0: Flush at each time step, no external buffers"
             "\n> 0.0: Flush at specified time interval, use external buffers");
  ADDOPTION ("fatigue",0,"Perform damage calculation on the gage stresses");
  ADDOPTION ("stressToMPaScale",1.0e-6,"Scale factor scaling stresses to MPa");
  ADDOPTION ("gate",25.0,"Stress gate value for the damage calculation [MPa]");
  ADDOPTION ("binSize",10.0,"Bin size for stress cycle counting [MPa]");
  ADDOPTION ("loga1",15.117,"Parameter log(a1) of the S-N curve");
  ADDOPTION ("loga2",17.146,"Parameter log(a2) of the S-N curve");
  ADDOPTION ("m1",4.0,"Parameter m1 of the S-N curve");
#ifdef win32
  const bool defEnd = true;
#else
  const bool defEnd = false;
#endif
  ADDOPTION ("littleEndian",defEnd,"Use Little Endian formatting of DAC files");

  // Now launch fedem_gage
  return solveGage();
}
