// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "../solverInterface.h"

#include <iostream>
#include <cstdlib>
#include <cctype>


int main (int argc, char** argv)
{
  if (argc < 2 || !isdigit(argv[1][0]))
  {
    std::cerr <<"usage: "<< argv[0] <<" <ntrial> [fedem_solver options]\n";
    return 1;
  }

  int status = 0, ntrial = atoi(argv[1]);
  argv[1] = argv[0];

  char answr = 'y';
  std::cout <<"Start ? ";
  std::cin >> answr;
  for (int i = 0; i < ntrial && !status && tolower(answr) != 'n'; i++)
  {
    status = solverInit(argc-1,argv+1,NULL,NULL,0,NULL,0);
    std::cout <<"Initialized "<< i <<": "<< status;
    std::cout <<"\nContinue ? ";
    std::cin >> answr;
    if (!status) status = solverDone(i == ntrial-1 || answr == 'n');
    std::cout <<"Done "<< status;
    std::cout <<"\nContinue ? ";
    std::cin >> answr;
  }

  return status;
}
