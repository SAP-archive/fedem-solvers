// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file cmdLineArgInit.C

  \brief This file contains C++ code for setting up the command-line handler.

  \details The file contains some extensions of the functions located in
  cmdLineArgInitStd.C, to deal with output redirection, program versioning,
  signal handling, etc.
*/

#include "SignalHandler.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFortran.H"
#include "Admin/FedemAdmin.H"
#include <cstring>

//! \cond DO_NOT_DOCUMENT
SUBROUTINE(closelogfiles,CLOSELOGFILES) ();
SUBROUTINE(setversion,SETVERSION) (const char* theVer, const int nchar);
SUBROUTINE(setdate,SETDATE) (const char* theDate, const int nchar);
SUBROUTINE(wrimsg,WRIMSG)
#ifdef _NCHAR_AFTER_CHARARG
  (const int& ifile, const char* message, const int nchar, const int&);
#define WRITE_MESSAGE(f,s,nchar,nl) F90_NAME(wrimsg,WRIMSG) (f,s,nchar,nl)
#else
  (const int& ifile, const char* message, const int&, const int nchar);
#define WRITE_MESSAGE(f,s,nchar,nl) F90_NAME(wrimsg,WRIMSG) (f,s,nl,nchar)
#endif

#define ADDOPTION FFaCmdLineArg::instance()->addOption

extern "C" void closeAllBinaryDB ();

extern void cmdLineArgInitStd (int argc, char** argv);
extern void readOptionFilesStd (const char* program);
//! \endcond


/*!
  \brief Cleans up on disk on emergency exists, exceptions, etc.
  \details Used by signal handlers only.
*/

static void emergencyExit ()
{
  // Enter here all clean-up needed on emergency exits, ctrl-C, etc.
  F90_NAME(closelogfiles,CLOSELOGFILES) (); // Close the res-file
  closeAllBinaryDB (); // Close all binary files (deleting temporary files)
}


/*!
  \brief A sub-class of FFaMsg that passes the list-messages to Fortran.

  \details All other member functions are empty (does nothing).
*/

class F90Msg : public FFaMsg
{
  int ifile; //!< File option (0: write to error log file, 1: write to res-file)

public:
  //! \brief Default constructor.
  F90Msg(int whichFile = 0) { ifile = whichFile; }
  //! \brief Empty destructor.
  virtual ~F90Msg() {}

protected:
  //! \brief Writes the \a message to the Fortran file, line by line.
  virtual void listVt (const std::string& message, bool)
  {
    // Write line by line omitting the '\n' character
    const char* str = message.c_str();
    const size_t n = message.length();
    size_t i, j = 0;
    for (i = 0; i < n; i++)
      if (str[i] == '\n')
      {
        WRITE_MESSAGE (ifile,str+j,(int)(i-j),1);
        j = i+1;
      }

    // Write last line, not terminated with '\n'
    if (j < n) WRITE_MESSAGE (ifile,str+j,(int)(n-j),0);
  }
};


/*!
  \brief Initializes the command-line parser.

  \param[in] argc Number of command-line arguments (including program name)
  \param[in] argv List of command-line arguments
*/

void cmdLineArgInit (int argc, char** argv)
{
  cmdLineArgInitStd (argc,argv);

  // Add command-line options which are common for all fortran modules
  ADDOPTION ("debug",0,"Debug print switch");
  ADDOPTION ("terminal",6,"File unit number for terminal output");
  ADDOPTION ("consolemsg",false,"Output error messages to console");
  ADDOPTION ("profile",false,"Print out profiling data at termination",false);
  ADDOPTION ("noSignalHandling",false,"Disable OS signal handling",false);
}


/*!
  \brief Redirects all error messages from C++ to F90.

  \param[in] whichFile Defines to which file we should write (see below)

  \details
  - \a whichFile &lt; 0: No redirection of C++ messages.
  - \a whichFile = 0: Redirection to Output List only.
  - \a whichFile = 1: Redirection to the res-file only (default).
  - \a whichFile = 2: Redirection to the res-file and Output List.
*/

void redirectOutput2Ftn (int whichFile)
{
  FFaMsg::setMessager (new F90Msg(whichFile));
}


/*!
  \brief Reads the command-line options files.

  \param[in] program Name of the program module
  \param[in] whichFile Option for console output redirection
*/

void readOptionFiles (const char* program, int whichFile)
{
  redirectOutput2Ftn (whichFile);
  readOptionFilesStd (program);

  // Initialize the signal handler
  SignalHandler::init (program,&emergencyExit);

  // Pass the version tag and date to the fortran modules
  std::string version(FedemAdmin::getVersion());
  const char* build_date = FedemAdmin::getBuildDate();
  if (!FedemAdmin::is64bit()) version += " (32bit)";
  F90_NAME(setversion,SETVERSION) (version.c_str(),(int)version.size());
  F90_NAME(setdate,SETDATE) (build_date,(int)strlen(build_date));
}
