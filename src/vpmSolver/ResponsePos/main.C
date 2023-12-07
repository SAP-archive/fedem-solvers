// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFortran.H"
#include "Admin/FedemAdmin.H"
#include <cstring>

#define ADDOPTION FFaCmdLineArg::instance()->addOption

SUBROUTINE(setversion,SETVERSION) (const char* theVer, const int nchar);
SUBROUTINE(setdate,SETDATE) (const char* theDate, const int nchar);
SUBROUTINE(response_pos,RESPONSE_POS) (int& iret);


int main (int argc, char** argv)
{
  // Initialize the command-line parser
  FFaCmdLineArg::init (argc,argv);

  // Initialize only the options that are actually evaluated from this program
  ADDOPTION ("consolemsg",true,"Output error messages to console");
  ADDOPTION ("rdbinc",1,"Increment number for the results database file");
  ADDOPTION ("errorfile","","Name of error message file");
  ADDOPTION ("linkfile","","Name of link input file");
  ADDOPTION ("deformation",false,"Save deformations to results database");
  ADDOPTION ("profile",false,"Print out profiling data at termination",false);

  // Get program version and build date
  const char* fedem_version = FedemAdmin::getVersion();
  const char* build_date = FedemAdmin::getBuildDate();

  // Pass the version tag and date to the fortran modules
  F90_NAME(setversion,SETVERSION) (fedem_version,strlen(fedem_version));
  F90_NAME(setdate,SETDATE) (build_date,strlen(build_date));

  // Lauch the Fortran program
  int iret = 0;
  F90_NAME(response_pos,RESPONSE_POS) (iret);
  FFaCmdLineArg::removeInstance ();
  return iret;
}
