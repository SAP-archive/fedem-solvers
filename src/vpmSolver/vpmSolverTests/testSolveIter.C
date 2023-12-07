// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testSolveIter.C

  \brief Test program for Newton iteration manipulation.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 30 Mar 2017
*/

#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "../solverInterface.h"


/*!
  \brief Test program for Newton iteration manipulation.

  \callgraph
*/

int main (int argc, char** argv)
{
  // Check for specification of an object base ID
  // for which to test equation number extraction for
  // and for specification of an external force to insert
  std::vector<double> f;
  int i = 1, baseId = 0;
  while (i < argc)
    if (!strcmp(argv[i],"-bid") && i+1 < argc)
    {
      baseId = atoi(argv[i+1]);
      if (i+2 < argc) memmove(argv+i,argv+i+2,(argc-i-2)*sizeof(char*));
      argc -= 2;
    }
    else if (!strcmp(argv[i],"-force") && i+1 < argc)
    {
      int j = i;
      while (++j < argc && isdigit(argv[j][0])) f.push_back(atof(argv[j]));
      if (j < argc) memmove(argv+i,argv+j,(argc-j)*sizeof(char*));
      argc -= j-i;
    }
    else
      i++;

  // Read input files, preprocess the model and setup the initial configuration
  int status = solverInit(argc,argv);
  if (status != 0) return status; // Simulation failed, aborting...

  if (baseId > 0)
  {
    // Check the equation numbers associated with the object with given baseID
    int meqn[6];
    int ndof = getEquations(baseId,meqn);
    std::cout <<"The object with baseID="<< baseId <<" has the equations:";
    for (i = 0; i < ndof; i++) std::cout <<" "<< meqn[i];
    std::cout << std::endl;
  }

  // Get the dimension of the system matrices and vector
  int ndim = getSystemSize();
  std::cout <<"Dimension of system matrices: "<< ndim << std::endl;

  // Time step loop
  for (int istep = 1; status == 0; istep++)
  {
    // Start a new time step
    bool doContinue = startStep(&status);
    if (status < 0) return status; // Simulation failed, aborting...
    if (!doContinue) break; // Simulation finished, terminate...

    // Iteration loop
    for (int iter = 0; doContinue && status == 0; iter++)
    {
      // Enter here any manipulatation of the current linearized system
      // before trying to solving it...
      if (istep == 1 && iter == 0)
      {
	double* Kmat = new double[ndim*ndim];
	if (!getNewtonMatrix(Kmat)) return -99;

	double* Rvec = new double[ndim];
	if (!getRhsVector(Rvec)) return -99;

	std::cout <<"\nHere is the stiffness matrix:";
	for (i = 1; i <= ndim; i++)
	{
	  std::cout <<"\nEquation "<< i <<":";
	  for (int j = 0; j < ndim; j++)
	    std::cout <<" "<< Kmat[j*ndim+i-1];
	}

	std::cout <<"\n\nHere is the right-hand-side force vector:";
	for (i = 1; i <= ndim; i++)
	  std::cout <<"\nEquation "<< i <<": "<< Rvec[i-1];
	std::cout << std::endl;

        if (!f.empty())
        {
          std::cout <<"\nForce vector to insert:";
          for (size_t j = 0; j < f.size() && j < (size_t)ndim; j++)
          {
            Rvec[j] = f[j];
            std::cout <<" "<< f[j];
          }
          std::cout << std::endl;
          if (!setRhsVector(Rvec)) return -99;
          f.clear();

          for (i = 0; i < ndim; i++) Rvec[i] = i;
          if (!getRhsVector(Rvec)) return -99;
          std::cout <<"Force vector inserted:";
          for (i = 0; i < ndim; i++) std::cout <<" "<< Rvec[i];
          std::cout << std::endl;
        }

        delete[] Kmat;
        delete[] Rvec;
      }

      // Solve the next iteration
      doContinue = solveIteration(&status);
      if (status < 0) return status; // Simulation failed, aborting...

      if (istep == 1)
      {
        double* Rvec = new double[ndim];
        if (!getRhsVector(Rvec)) return -99;
        std::cout <<"Next force vector:";
        for (i = 0; i < ndim; i++) std::cout <<" "<< Rvec[i];
        std::cout << std::endl;
      }
    }
  }

  // Simulation finished, terminate by closing down the result database, etc.
  return solverDone();
}
