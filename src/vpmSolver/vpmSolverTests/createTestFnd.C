// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file createTestFnd.C
  \brief This program is used to create a binary file for unit testing.

  \details To compile on Linux:
  \code
  g++ -std=c++11 createTestFnd.C -Ipublic/commonCore -Lpublic/commonCore/FFaLib/FFaOS -lFFaOS
  \endcode
*/

#include "FFaLib/FFaOS/FFaTag.H"
#include <vector>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cstring>


int main (int argc, char** argv)
{
  // Open binary file for output
  const char* fileName = "motions.fnd";
  FT_FILE fileDes = FT_open(fileName,FT_wb);
  if (fileDes <= (FT_FILE)0)
  {
    perror(fileName);
    return 1;
  }

  // Write file header identifying this file type
  const char* fileTag = "#FEDEM nodal displacements";
  if (FFaTag::write(fileDes,fileTag,strlen(fileTag),12345678,LEN_TAG) < 0 ||
      FT_write(";7.3;",1,5,fileDes) < 5)
  {
    perror(fileName);
    FT_close(fileDes);
    return 2;
  }

  size_t nBytes = FT_tell(fileDes);

  // Lambda function that writes data to the binary file.
  auto&& writeData = [fileDes,&nBytes](const void* data, size_t n, size_t m = 1)
  {
    int nWrote = FT_write(data,n,m,fileDes);
    if (nWrote > 0)
    {
      nBytes += n*nWrote;
      return true;
    }

    std::cerr <<"Write failure."<< std::endl;
    return false;
  };

  // Write some node numbers to the binary file
  std::vector<int> nodes = { 7, 12, 38, 42, 59 };
  std::string cNodes = std::to_string(nodes.size()) + ";\n";
  nBytes += FT_write(cNodes.c_str(),1,cNodes.size(),fileDes);
  if (!writeData(nodes.data(),sizeof(int),nodes.size()))
    return 3;

  std::vector<double> displ(3*nodes.size());
  std::iota(displ.begin(),displ.end(),1.0);

  // Time step loop
  double time = 0.0;
  for (int step = 0; step < 3; step++, time += 0.1)
    if (!writeData(&step,sizeof(int)) || !writeData(&time,sizeof(double)))
      return 4;
    else if (!writeData(displ.data(),sizeof(double),displ.size()))
      return 5;
    else
    {
      std::cout << std::setw(16) << nBytes <<" bytes written."<< std::endl;
      for (double& v : displ) v += 100.0;
    }

  // Close up and quit
  std::cout <<"\nDone."<< std::endl;
  FT_close(fileDes);
  return 0;
}
