!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine GetConnectivityFEData (T, nels, nnpc, xnpc, npc)

  use kind_values, only : is
  use FEData     , only : FEDataInput

  implicit none

  type(FEDataInput), intent(in)  :: T
  integer(is)      , intent(out) :: nels, nnpc
  integer(is)      , pointer     :: xnpc(:), npc(:)

  nels =  T%sam%nel
  nnpc =  T%sam%nmmnpc
  xnpc => T%sam%mpmnpc
  npc  => T%sam%mmnpc

end subroutine GetConnectivityFEData


subroutine GetConstraintsFEData (T, nceq, nnceq, xceq, ceq, tcc)

  use kind_values, only : is, wp
  use FEData     , only : FEDataInput

  implicit none

  type(FEDataInput), intent(in)        :: T
  integer(is)      , intent(out)       :: nceq, nnceq
  integer(is)      , pointer           :: xceq(:), ceq(:)
  real(wp)         , pointer, optional :: tcc(:)

  nceq  =  T%sam%nceq
  nnceq =  T%sam%nmmceq
  xceq  => T%sam%mpmceq
  ceq   => T%sam%mmceq
  if (present(tcc)) tcc => T%sam%ttcc

end subroutine GetConstraintsFEData


subroutine GetPartitionFEData (T, nnod, ndof, madof, msc, minex)

  use kind_values, only : is
  use FEData     , only : FEDataInput

  implicit none

  type(FEDataInput), intent(in)  :: T
  integer(is)      , intent(out) :: nnod, ndof
  integer(is)      , pointer     :: madof(:), msc(:), minex(:)

  nnod  =  T%sam%nnod
  ndof  =  T%sam%ndof
  madof => T%sam%madof
  msc   => T%sam%msc
  minex => T%sam%minex

end subroutine GetPartitionFEData
