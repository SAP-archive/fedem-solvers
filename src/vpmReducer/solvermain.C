// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file solvermain.C

  \brief This file contains the C++ main program for the FEDEM linear FE solver.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 15 Aug 2012
*/

extern "C" {
  int initSolverArgs(int argc, char** argv, bool first = true);
  int solvePart();
}


/*!
  \brief Main program for the FEDEM linear FE part solver.

  \details The main program contains very little logic.
  It uses two functions from the reducerInterface API to initialize the
  command-line arguments handler, and then to solve the specified FE part.

  \callgraph
*/

int main (int argc, char** argv)
{
  return initSolverArgs(argc,argv) == 0 ? solvePart() : -1;
}
