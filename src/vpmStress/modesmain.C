// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "vpmCommon/FFCmdLineArgInterface.H"

extern "C" {
  void initSolverArgs(int,char**);
  int  solveModes();
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
  ADDOPTION ("samfile","","Name of SAM data file");
  ADDOPTION ("fsifile","fedem_solver.fsi","Name of solver input file");
  ADDOPTION ("frsfile","","Name of solver results database file");
  ADDOPTION ("rdbfile","","Name of modes results database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database file");
  ADDOPTION ("VTFfile","","Name of VTF output file");
  ADDOPTION ("VTFoffset",0,"VTF result block id offset");
  ADDOPTION ("VTFparts",0,"Number of parts in VTF-file");
  ADDOPTION ("VTFexpress",false,"Write express VTF-files (one file per mode)");
  ADDOPTION ("VTFdscale",1.0,"Deformation scaling factor for VTF output");
  ADDOPTION ("double",false,"Save all results in double precision");
  ADDOPTION ("damped",false,"Complex modes are calculated");
  ADDOPTION ("recover_modes","","List of mode numbers to expand");
  ADDOPTION ("write_nodes",false,"Save results as nodal data");
  ADDOPTION ("write_vector",true,"Save results as vector data");
  ADDOPTION ("energy_density",false,"Save scaled strain energy density");
  ADD_PRIVATE_OPTION("stressForm",0,"General stress formulation option");
  ADD_PRIVATE_OPTION("ffqStressForm",2,"Stress formulation for the FFQ shell");
  ADD_PRIVATE_OPTION("fftStressForm",1,"Stress formulation for the FFT shell");
  ADD_PRIVATE_OPTION("useIncompatibleModes",false,"Linear hexahedron option");

  // Now launch fedem_modes
  return solveModes();
}
