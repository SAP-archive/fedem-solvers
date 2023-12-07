// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file cmdLineArgInitStd.C

  \brief This file contains C++ code for setting up the command-line handler.

  \details The file contains some global functions that are invoked by the
  C++ main programs for initializing the command-line argument handler,
  printing of help texts to console, etc.
*/

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaDefinitions/FFaAppInfo.H"
#include "FFaLib/FFaOS/FFaFilePath.H"
#include "Admin/FedemAdmin.H"

#include <iostream>

#if defined(win32) || defined(win64)
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <stdio.h>

#if _MSC_VER > 1310
#define chdir _chdir
#define getcwd _getcwd
#endif

//! \brief Convenience macro.
#define ADDOPTION FFaCmdLineArg::instance()->addOption

//! \brief Previous working directory, to return to after finished.
static char* oldwd = NULL;


/*!
  \brief Resets the current working directory after the solver has finished.
*/

void resetCWD ()
{
  if (!oldwd)
    return;
  else if (chdir(oldwd))
    perror(oldwd);

  free(oldwd);
  oldwd = NULL;
}


/*!
  \brief Initializes the command-line parser.

  \param[in] argc Number of command-line arguments (including program name)
  \param[in] argv List of command-line arguments

  \details This function also defines some general command-line options that
  are common for all solver modules.
*/

void cmdLineArgInitStd (int argc, char** argv)
{
#ifdef FT_DEBUG
  if (argc > 0)
    std::cout <<"\n### Starting execution of "<< argv[0] <<" ###"<< std::endl;
#endif
  FFaAppInfo::init (argv[0]);
  FFaCmdLineArg::init (argc,argv);

  ADDOPTION ("fao","","Read additional options from this file");
  ADDOPTION ("fco","","Read calculation options from this file");
  ADDOPTION ("fop","","Read output options from this file");
  ADDOPTION ("cwd","","Change working directory");
  ADDOPTION ("help",false,"Print out this help text");
  ADDOPTION ("helpAll",false,"Print out this help text"
             "\nincluding the private options, if any",false);
  ADDOPTION ("helpSL",false,"Print out this help text"
             " as tab-separated one-liners",false);
  ADDOPTION ("version",false,"Print out program version");
}


/*!
  \brief Reads the command-line options files.

  \param[in] program Name of the program module

  \details The help text containing all available command-line options
  may also be printed to console if the -help (or -helpAll) is specified.
*/

void readOptionFilesStd (const char* program)
{
  // Change the working directory, if specified
  std::string cwd;
  FFaCmdLineArg::instance()->getValue ("cwd",cwd);
  if (!cwd.empty())
  {
    oldwd = getcwd(NULL,128);
    FFaFilePath::checkName (cwd);
    if (chdir(cwd.c_str()))
    {
      // Invalid directory, abort...
      perror((program + (": " + cwd)).c_str());
      free(oldwd);
      exit(1);
    }
#ifdef FT_DEBUG
    std::cout <<"Executing in "<< cwd << std::endl;
#endif
  }

  // Read the options files, if any
  std::string optfile;
  FFaCmdLineArg::instance()->getValue ("fao",optfile);
  if (!optfile.empty()) FFaFilePath::checkName (optfile);
  FFaCmdLineArg::instance()->readOptionsFile (optfile);

  FFaCmdLineArg::instance()->getValue ("fco",optfile);
  if (!optfile.empty()) FFaFilePath::checkName (optfile);
  FFaCmdLineArg::instance()->readOptionsFile (optfile);

  FFaCmdLineArg::instance()->getValue ("fop",optfile);
  if (!optfile.empty()) FFaFilePath::checkName (optfile);
  FFaCmdLineArg::instance()->readOptionsFile (optfile);

  // Get the program version and build date
  const char* fedem_version = FedemAdmin::getVersion();
  const char* build_date = FedemAdmin::getBuildDate();

  // Print the help text, if wanted
  bool wantHelp = false, wantAllHelp = false;
  FFaCmdLineArg::instance()->getValue ("help",wantHelp);
  FFaCmdLineArg::instance()->getValue ("helpAll",wantAllHelp);
  if (wantHelp || wantAllHelp)
  {
    std::string hlpText;
    FFaCmdLineArg::instance()->composeHelpText (hlpText,wantAllHelp);
    std::cout <<"\n"<< program <<" version "<< fedem_version <<" "<< build_date;
    if (!FedemAdmin::is64bit()) std::cout <<"   (32bit)";
    std::cout <<"\n\nAvailable command-line options:\n"<< hlpText;
    resetCWD();
    exit(0);
  }

  // Print the help text as tab-separated columns, one line per entry
  FFaCmdLineArg::instance()->getValue ("helpSL",wantHelp);
  if (wantHelp)
  {
    std::string hlpTxt;
    FFaCmdLineArg::instance()->composeSingleLineHelpText (hlpTxt);
    std::cout << hlpTxt;
    resetCWD();
    exit(0);
  }

  // Print the program version, if wanted
  FFaCmdLineArg::instance()->getValue ("version",wantHelp);
  if (wantHelp)
  {
    std::cout << program <<" version "<< fedem_version <<" "<< build_date;
    if (!FedemAdmin::is64bit()) std::cout <<"   (32bit)";
    std::cout <<"\n";
    resetCWD();
    exit(0);
  }

  // Check for expiration (for internal releases)
  int expired = FedemAdmin::getExpireAfter();
  if (expired > 0)
  {
    std::cout << program <<" version "<< fedem_version <<" "<< build_date
              <<"\n\n";
    int myAge = FedemAdmin::getDaysSinceBuilt();
    if (myAge >= expired)
    {
      std::cout <<"*** This version expired "<< myAge-expired
                <<" days ago ***\n";
      resetCWD();
      exit(1);
    }
    else
      std::cout <<"*** This is an internal (non-public) version that will "
                <<" expire in "<< expired-myAge <<" days ***"<< std::endl;
  }
}
