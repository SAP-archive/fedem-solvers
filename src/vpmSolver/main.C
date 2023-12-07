/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file vpmSolver/main.C

  \brief This file contains the C++ main program for the FEDEM dynamics solver.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 2 Dec 2016
*/

#include "solverInterface.h"
#include <stdlib.h>
#include <string.h>

int compareResponse (const char*, const char*, double, int);


/*!
  \brief Main program for the FEDEM dynamics solver.

  \details The main program contains very little logic.
  It uses functions from the solverInterface API to initialize and solve the
  dynamic problem at each time step. It can also invoke a response verification
  before program termination, in which the calculated response is compared with
  some reference data. This is mainly used to set up regression tests.

  \callgraph
*/

int main (int argc, char** argv)
{
  // Check if a file for response verification was specified
  const char* verify = NULL;
  double epsTol = 0.0;
  int skipIL = 0;
  if (argc > 3 && !strcmp(argv[argc-3],"-verify"))
  {
    verify = argv[argc-2];
    epsTol = atof(argv[argc-1]);
    argc -= 3;
  }
  else if (argc > 4 && !strcmp(argv[argc-4],"-verify"))
  {
    verify = argv[argc-3];
    epsTol = atof(argv[argc-2]);
    skipIL = atoi(argv[argc-1]);
    argc -= 4;
  }

  // Read input files, preprocess the model and setup the initial configuration
  int status = solverInit(argc,argv);
  if (status) return status;

  // Time step loop.
  // Invoke the solver step-by-step until specified end time is reached,
  // or an error occur.
  while (solveNext(&status))
    if (status) return status; // Simulation failed, aborting...

  // Simulation finished, terminate by closing down the result database, etc.
  if (status || (status = solverDone()))
    return status;

  // Verify exported curve data against the reference data, if specified
  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-curvePlotFile") && i+1 < argc)
      return compareResponse(argv[i+1],verify,epsTol,skipIL);

  return 0;
}
