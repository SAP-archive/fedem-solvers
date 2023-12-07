// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFaLib/FFaOS/FFaFortran.H"
#include <vector>
#include <set>
#ifdef FT_DEBUG
#include <iostream>
#include <iomanip>
#endif
#include <cstddef>

static std::vector< std::set<int> > irow;


SUBROUTINE(spr_init,SPR_INIT) (const int& neq)
{
  irow.clear();
  irow.resize(neq);
}


SUBROUTINE(spr_preass,SPR_PREASS) (const int& nedof, const int* meen,
                                   const int* meqn, const int* mpmceq,
                                   const int* mmceq)
{
  int i, j, ip, jp, ieq, jeq, iceq, jceq;

  for (i = 0; i < nedof; i++)
    if ((ieq = meen[i]) > 0)
      for (j = 0; j < nedof; j++)
        if ((jeq = meen[j]) >= ieq) // Upper triangle only (assuming symmtry)
          irow[ieq-1].insert(jeq);

  for (i = 0; i < nedof; i++)
    if ((iceq = -meen[i]) > 0)
      for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
        if (mmceq[ip] > 0)
        {
          ieq = meqn[mmceq[ip]-1];
          for (j = 0; j < nedof; j++)
            if ((jeq = meen[j]) >= ieq) // Upper triangle only
              irow[ieq-1].insert(jeq);
            else if (jeq > 0)
              irow[jeq-1].insert(ieq);
            else if ((jceq = -jeq) > 0)
              for (jp = mpmceq[jceq-1]; jp < mpmceq[jceq]-1; jp++)
                if (mmceq[jp] > 0)
                  if ((jeq = meqn[mmceq[jp]-1]) >= ieq) // Upper triangle only
                    irow[ieq-1].insert(jeq);
        }
}


INTEGER_FUNCTION(spr_getnnz,SPR_GETNNZ) ()
{
  int nnz = 0;
  for (size_t i = 0; i < irow.size(); ++i)
    nnz += irow[i].size();

#ifdef FT_DEBUG
  std::cout <<" Size(NNZ) of sparse system matrix: "<< nnz << std::endl;
#endif
  return nnz;
}


SUBROUTINE(spr_getpattern,SPR_GETPATTERN) (int* ia, int* ja)
{
  size_t i, j = 0;
  std::set<int>::const_iterator it;

  ia[0] = 1;
  for (i = 0; i < irow.size(); ia[++i] = j+1)
    for (it = irow[i].begin(); it != irow[i].end(); ++it, ++j)
      ja[j] = *it;

#if FT_DEBUG > 1
  for (i = 0; i < irow.size(); i++)
  {
    std::cout <<" Eq"<< std::setw(5) << i+1 <<":";
    for (j = ia[i]; j < ia[i+1]; j++)
      std::cout <<" "<< ja[j-1];
    std::cout << std::endl;
  }
#endif
  irow.clear();
}
