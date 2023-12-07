// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file readFSI.C

  \brief Utility for reading read a solver input file into a string variable.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 19 Dec 2016
*/

#include <fstream>
#include <string>
#include <cstring>
#include <cctype>
#include <cstdio>


/*!
  \brief Utility function to read a solver input file into a string variable.
  \param[in] cwd Current working directory
  \param[in] fname Path to the solver input file on disk
  \param[out] chfsi Text string containing the compressed file content

  \details All blank lines and comment lines are ignored, and all white spaces
  are replaced by a single space character, to make the string more compact.
  A dollar-sign ($) is appended at the end of each line instead of the newline
  character, since the latter is not so well handled by the Fortran routines
  that read the solver input.

  \callergraph
*/

bool readFsiFile (const char* cwd, const char* fname, std::string& chfsi)
{
  if (!fname || strlen(fname) < 1) return true;

  std::string fName;
  std::ifstream fsi;
  if (!cwd || strlen(cwd) < 1)
    fsi.open(fname,std::ios::in);
  else if (fname[0] == '\\' || fname[0] == '/')
    fsi.open(fname,std::ios::in); // An absolute path was specified, ignore cwd
  else if (strlen(fname) > 1 && isalpha(fname[0]) && fname[1] == ':')
    fsi.open(fname,std::ios::in); // Absolute path on Windows, ignore cwd
  else
  {
    // A current working directory was specified.
    // The input file path is then assumed to be relative to this path.
#if defined(win32) || defined(win64)
    const char slash = '\\';
#else
    const char slash = '/';
#endif
    fName = cwd;
    if (fName.find_last_of(slash) != fName.size()-1)
      fName.append(1,slash);
    fName.append(fname);
    fname = fName.c_str();
    fsi.open(fname,std::ios::in);
  }

  if (!fsi)
  {
    perror(fname);
    return false;
  }

  char cline[1024];
  while (fsi.getline(cline,1024))
    if (strlen(cline) > 0 && cline[0] != '\'') // ignore comment and blank lines
    {
      size_t nchar = strlen(cline);
#if !defined(win32) && !defined(win64)
      // Remove trailing \r on UNIX systems
      if (cline[nchar-1] == 13) --nchar;
      if (nchar == 0) continue;
#endif
      // Reduce all spaces to a single space
      size_t i = 0, j = 0;
      bool inSpace = false;
      for (j = 0; j < nchar; j++)
        if (cline[j] != ' ')
        {
          cline[i++] = cline[j];
          inSpace = false;
        }
        else if (!inSpace)
        {
          if (j > 0)
            cline[i++] = ' ';
          inSpace = true;
        }

      cline[i] = '$'; // end-of-line marker
      chfsi.append(cline,i+1);
    }

  return true;
}
