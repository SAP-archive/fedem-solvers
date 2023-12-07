// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "vpmCommon/FFCmdLineArgInterface.H"

extern "C" {
  void initSolverArgs(int,char**);
  int  solveStress();
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
  ADDOPTION ("rdbfile","","Name of stress results database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database file");
  ADDOPTION ("VTFfile","","Name of VTF output file");
  ADDOPTION ("VTFoffset",0,"VTF result block id offset");
  ADDOPTION ("VTFparts",0,"Number of parts in VTF-file");
  ADDOPTION ("VTFavgelm",true,"Write averaged element results to VTF-file");
  ADDOPTION ("VTFinit",false,"Write initial state to VTF-file");
  ADDOPTION ("VTFdscale",1.0,"Deformation scaling factor for VTF output");
  ADDOPTION ("double",false,"Save all results in double precision");
  ADDOPTION ("group","","List of element groups to do calculations for");
  ADDOPTION ("nodalForces",false,"Compute and print nodal forces");
  ADDOPTION ("SR",false,"Save stress resultants to results database");
  ADDOPTION ("stress",false,"Save stress tensors to results database");
  ADDOPTION ("strain",false,"Save strain tensors to results database");
  ADDOPTION ("vmStress",false,"Save von Mises stress to results database");
  ADDOPTION ("vmStrain",false,"Save von Mises strain to results database");
  ADDOPTION ("maxPStress",false,"Save max principal stress to results database");
  ADDOPTION ("maxPStrain",false,"Save max principal strain to results database");
  ADDOPTION ("minPStress",false,"Save min principal stress to results database");
  ADDOPTION ("minPStrain",false,"Save min principal strain to results database");
  ADDOPTION ("maxSStress",false,"Save max shear stress to results database");
  ADDOPTION ("maxSStrain",false,"Save max shear strain to results database");
  ADDOPTION ("deformation",false,"Save deformations to results database");
  ADDOPTION ("dumpDefNas",false,"Save deformations to Nastran bulk data files");
  ADDOPTION ("write_nodes",true,"Save deformations as nodal data");
  ADDOPTION ("write_vector",false,"Save deformations as vector data");
  ADDOPTION ("statm",0.0,"Start time");
  ADDOPTION ("stotm",1.0,"Stop time");
  ADDOPTION ("tinc",0.1,"Time increment (= 0.0: process all time steps)");

  ADD_PRIVATE_OPTION("stressForm",0,"General stress formulation option\n"
                     "= 0: Direct evaluation in nodes\n"
                     "= 1: Volume averaged or mid-point evaluation\n"
                     "= 2: Extrapolation from Gauss integration points");
  ADD_PRIVATE_OPTION("ffqStressForm",2,"Stress formulation for the FFQ shell\n"
                     "= 0: Direct evaluation in nodes\n"
                     "= 1: Extrapolation from mid-point\n"
                     "= 2: Extrapolation from 2x2 gauss points");
  ADD_PRIVATE_OPTION("fftStressForm",1,"Stress formulation for the FFT shell\n"
                     "= 1: Membrane part based on the HLST element\n"
                     "= 2: Membrane part based on the TMRF element");
  ADD_PRIVATE_OPTION("useIncompatibleModes",false,"Linear hexahedron option");

  // Now launch fedem_stress
  return solveStress();
}
