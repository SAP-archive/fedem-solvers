// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file FFCmdLineArgInterface.C

  \brief Wrapper for FFaCmdLineArg to facilitate use in shared object libraries.
*/

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"

#if defined(win32) || defined(win64)
#include <windows.h>
//! \brief Export of symbols form shared libraries on Windows
#define DLLexport extern "C" __declspec(dllexport)
#else
//! \brief Export of symbols form shared libraries on Linux
#define DLLexport extern "C"
#endif


//! \brief Adds an int-valued command-line option.
DLLexport void AddIntCmdLineOption (const char* ident, int defaultVal,
                                    const char* descr, bool showToAll)
{
  FFaCmdLineArg::instance()->addOption(ident,defaultVal,descr,showToAll);
}


//! \brief Adds a float-valued command-line option.
DLLexport void AddFloatCmdLineOption (const char* ident, float defaultVal,
                                      const char* descr, bool showToAll)
{
  FFaCmdLineArg::instance()->addOption(ident,defaultVal,descr,showToAll);
}


//! \brief Adds a double-valued command-line option.
DLLexport void AddDoubleCmdLineOption (const char* ident, double defaultVal,
                                       const char* descr, bool showToAll)
{
  FFaCmdLineArg::instance()->addOption(ident,defaultVal,descr,showToAll);
}


//! \brief Adds a string-valued command-line option.
DLLexport void AddCmdLineOption (const char* ident, const char* defaultVal,
                                 const char* descr, bool showToAll)
{
  FFaCmdLineArg::instance()->addOption(ident,std::string(defaultVal),
                                       descr,showToAll);
}


//! \brief Adds a bool-valued command-line option.
DLLexport void AddBoolCmdLineOption (const char* ident, bool defaultVal,
                                     const char* descr, bool showToAll)
{
  FFaCmdLineArg::instance()->addOption(ident,defaultVal,descr,showToAll);
}


//! \brief Adds a vector-of-doubles-valued command-line option.
DLLexport void AddDblVecCmdLineOption (const char* ident,
                                       const DoubleVec& defaultVal,
                                       const char* descr, bool showToAll)
{
  FFaCmdLineArg::instance()->addOption(ident,defaultVal,descr,showToAll);
}
