/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file computerConfig.C
  \brief Global functions for extracting information on the running computer.
*/

#include "FFaLib/FFaOS/FFaFortran.H"

#if defined(win32) || defined(win64)
#include <windows.h>
#else
#include <sys/utsname.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>


/*!
  \brief Static helper to append a string to a fixed-size character buffer.
*/

static void appends (char* dest, const char* src, const int nchar)
{
  int mchar = nchar - strlen(dest) - 2;
  if (mchar > 0)
    strncat(strcat(dest," "),src,mchar);
}


/*!
  \brief Returns a text string identifying the computer the program is run on.
*/

SUBROUTINE(getcomputerconfig,GETCOMPUTERCONFIG) (char* cid, const int nchar)
{
  size_t l;
#if defined(win32) || defined(win64)
  DWORD cchBuff = BUFSIZ;
  TCHAR tchBuffer[BUFSIZ];
  LPTSTR hn = tchBuffer;
  char* os; char* id;
  GetComputerName(hn,&cchBuff);
  strncpy(cid,hn,nchar);
  os = getenv("OS");
  appends(cid,os,nchar);
  id = getenv("PROCESSOR_IDENTIFIER");
  appends(cid,id,nchar);
#else
  struct utsname name;
  uname(&name);
  strncpy(cid,name.nodename,nchar);
  appends(cid,name.sysname,nchar);
  appends(cid,name.release,nchar);
  appends(cid,name.version,nchar);
  appends(cid,name.machine,nchar);
#endif
  l = strlen(cid);
  if ((int)l < nchar) memset(cid+l,' ',nchar-l);
}


/*!
  \brief Returns the user name that runs the program.
*/

SUBROUTINE(getusername,GETUSERNAME) (char* cuser, const int nchar)
{
  size_t l;
#if defined(win32) || defined(win64)
  DWORD cchBuff = BUFSIZ;
  TCHAR tchBuffer[BUFSIZ];
  LPTSTR usr = tchBuffer;
  GetUserName(usr,&cchBuff);
#else
  char* usr = getenv("USER");
#endif
  strncpy(cuser, usr ? usr : "(none)", nchar-1);
  l = usr ? strlen(usr) : 6;
  if ((int)l < nchar) memset(cuser+l,' ',nchar-l);
}
