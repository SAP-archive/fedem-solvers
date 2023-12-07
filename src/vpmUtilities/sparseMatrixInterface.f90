!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module SparseMatrixUtils

  implicit none

  interface

     subroutine spr_init (neq)
       integer, intent(in) :: neq
     end subroutine spr_init

     subroutine spr_preass (nedof, meen, meqn, mpmceq, mmceq)
       integer, intent(in) :: nedof, meen, meqn, mpmceq, mmceq
     end subroutine spr_preass

     function spr_getnnz ()
       integer :: spr_getnnz
     end function spr_getnnz

     subroutine spr_getpattern (ia,ja)
       integer, intent(out) :: ia, ja
     end subroutine spr_getpattern

  end interface

end module SparseMatrixUtils
