/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file clksec.c
  \brief System-dependent global functions for time measurement.
*/

#include "FFaLib/FFaOS/FFaFortran.H"

#if defined(win32) || defined(win64)
#include <time.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#include <stdlib.h> // to get the definition of the NULL macro
#endif

static unsigned long years = 1261440000; /*!< Seconds in 40 years (1970-2010) */


/*!
  \brief Returns current wall clock time in seconds.

  This machine dependent function returns as its value the wall clock time in
  seconds since Jan 1 2003 00:00:00, minus the value of the input argument \a *T.
*/

REAL_FUNCTION(clksec,CLKSEC) (const float* T)
{
#if defined(win32) || defined(win64)
  struct _timeb val;
  _ftime(&val);
  return (float)(val.time-years) + 0.001f*val.millitm - *T;
#else
  struct timeval val;
  gettimeofday(&val,NULL);
  return (float)(val.tv_sec-years) + 0.000001f*val.tv_usec - *T;
#endif
}


/*!
  \brief Initializes the static years variable.

  This machine dependent function re-initializes the number of seconds passed
  since Jan 1 1970, and up to the time of its invokation. The function should
  be called before the first invokation of CLKSEC. It ensures that the value
  returned by CLKSEC always is small enough to be represented by a float.
*/

SUBROUTINE(clkini,CLKINI) ()
{
#if defined(win32) || defined(win64)
  struct _timeb val;
  _ftime(&val);
  years = (unsigned long)val.time;
#else
  struct timeval val;
  gettimeofday(&val,NULL);
  years = (unsigned long)val.tv_sec;
#endif
}
