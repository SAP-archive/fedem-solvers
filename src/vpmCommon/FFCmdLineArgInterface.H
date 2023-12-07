// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFCmdLineArgInterface.H

  \brief Wrapper for FFaCmdLineArg to facilitate use in shared object libraries.

  \details The purpose of this wrapper is to redefine the ADDOPTION symbol into
  functions, that can be used when the command-line argument handler is built
  into a shared library, while the main program using it is not in that library.
*/

#ifndef FF_CMD_LINE_ARG_INTERFACE_H
#define FF_CMD_LINE_ARG_INTERFACE_H

#ifdef FT_USE_DLL
#include <string>
#include <vector>

typedef std::vector<double> DoubleVec; //!< Convenience type definition

extern "C" {
  void AddIntCmdLineOption   (const char*, int        , const char*, bool);
  void AddFloatCmdLineOption (const char*, float      , const char*, bool);
  void AddDoubleCmdLineOption(const char*, double     , const char*, bool);
  void AddCmdLineOption      (const char*, const char*, const char*, bool);
  void AddBoolCmdLineOption  (const char*, bool       , const char*, bool);
  void AddDblVecCmdLineOption(const char*, const DoubleVec&, const char*, bool);
}


//! \brief Adds an int-valued command-line option.
void ADDOPTION (const char* ident, int defaultVal,
                const char* descr, bool showToAll = true)
{
  AddIntCmdLineOption(ident,defaultVal,descr,showToAll);
}

//! \brief Adds a float-valued command-line option.
void ADDOPTION (const char* ident, float defaultVal,
                const char* descr, bool showToAll = true)
{
  AddFloatCmdLineOption(ident,defaultVal,descr,showToAll);
}

//! \brief Adds a double-valued command-line option.
void ADDOPTION (const char* ident, double defaultVal,
                const char* descr, bool showToAll = true)
{
  AddDoubleCmdLineOption(ident,defaultVal,descr,showToAll);
}

//! \brief Adds a string-valued command-line option.
void ADDOPTION (const char* ident, const char* defaultVal,
                const char* descr, bool showToAll = true)
{
  AddCmdLineOption(ident,defaultVal,descr,showToAll);
}

//! \brief Adds a bool-valued command-line option.
void ADDOPTION (const char* ident, bool defaultVal,
                const char* descr, bool showToAll = true)
{
  AddBoolCmdLineOption(ident,defaultVal,descr,showToAll);
}

//! \brief Adds a vector-of-doubles-valued command-line option.
void ADDOPTION (const char* ident, const DoubleVec& defaultVal,
                const char* descr, bool showToAll = true)
{
  AddDblVecCmdLineOption(ident,defaultVal,descr,showToAll);
}


#else
#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"

//! \brief Convenience macro for adding command-line options.
//! \details Used when both FFaCmdLineArg and the user of it (the main program),
//! reside in the same library (or in two different static libraries).
#define ADDOPTION FFaCmdLineArg::instance()->addOption

#endif

//! \brief Convenience macro, for adding private command-line options.
#define ADD_PRIVATE_OPTION(a,b,c) ADDOPTION(a,b,c,false)

#endif
