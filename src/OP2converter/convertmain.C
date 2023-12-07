// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFilePath.H"
#include "FFaLib/FFaOS/FFaFortran.H"
#include "Admin/FedemAdmin.H"

#include <iostream>

#if defined(win32) || defined(win64)
#include <direct.h>
#else
#include <unistd.h>
#endif

#if _MSC_VER > 1310
#define chdir _chdir
#endif

#define ADDOPTION FFaCmdLineArg::instance()->addOption

SUBROUTINE(convert_cms,CONVERT_CMS) (int& iret);


int main (int argc, char** argv)
{
  const char* program = "fedem_op2fmx";

  // Initialize the command-line parser
  FFaCmdLineArg::init(argc,argv);

  // Define all command-line options for this program
  ADDOPTION ("cwd","","Change working directory");
  ADDOPTION ("help",false,"Display this help and exit");
  ADDOPTION ("version",false,"Display program version and exit");
  ADDOPTION ("partName","","Fedem FE part name");
  ADDOPTION ("ndim",0,"Total number of DOFs (in ASET)");

  // Get the program version and build date
  const char* fedem_version = FedemAdmin::getVersion();
  const char* build_date = FedemAdmin::getBuildDate();

  // Print the help text, if wanted
  bool wantHelp;
  FFaCmdLineArg::instance()->getValue("help",wantHelp);
  if (wantHelp)
  {
    std::string hlpText;
    FFaCmdLineArg::instance()->composeHelpText(hlpText);
    std::cout <<"\n"<< program <<" version "<< fedem_version <<" "<< build_date;
    if (!FedemAdmin::is64bit()) std::cout <<"   (32bit)";
    std::cout <<"\n\nAvailable command-line options:\n"<< hlpText;
    return 0;
  }

  // Print the program version, if wanted
  FFaCmdLineArg::instance()->getValue("version",wantHelp);
  if (wantHelp)
  {
    std::cout << program <<" version "<< fedem_version <<" "<< build_date;
    if (!FedemAdmin::is64bit()) std::cout <<"   (32bit)";
    std::cout <<"\n";
    return 0;
  }

  std::cout <<"\nFedem OP2 conversion utility "
            << fedem_version <<" "<< build_date;
  if (!FedemAdmin::is64bit()) std::cout <<"   (32bit)";
  std::cout <<"\n\n";

  // Change the working directory, if specified
  std::string cwd;
  FFaCmdLineArg::instance()->getValue("cwd",cwd);
  if (!cwd.empty())
  {
    FFaFilePath::checkName(cwd);
    if (chdir(cwd.c_str()))
    {
      // Invalid directory, abort...
      perror((program + (": " + cwd)).c_str());
      return 1;
    }
  }

  // Now launch convert_cms
  int iret = 0;
  if (FFaCmdLineArg::instance()->isOptionSetOnCmdLine("partName"))
    F90_NAME(convert_cms,CONVERT_CMS) (iret);
  else
    std::cout <<"Execute \'"<< program <<" -help\' to see available options.\n";

  return iret;
}
