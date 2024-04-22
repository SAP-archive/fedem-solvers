// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file reducermain.C

  \brief This file contains the C++ main program for the FEDEM FE part reducer.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 15 Oct 2000
*/

extern "C" {
  int initSolverArgs(int argc, char** argv, bool first, bool muted);
  int reducePart(bool first = true, bool last = true);
}


/*!
  \brief Main program for the FEDEM FE part reducer.

  \details The main program contains very little logic.
  It uses two functions from the reducerInterface API to initialize the
  command-line arguments handler, and then to reduce the specified FE part.

  \callgraph
*/

int main (int argc, char** argv)
{
  return initSolverArgs(argc,argv,true,false) == 0 ? reducePart() : -1;
}
