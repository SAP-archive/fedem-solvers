!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file sprExtensionModule.f90
!> @brief Sparse matrix extensions.

!!==============================================================================
!> @brief Sparse matrix extensions.
!> @details This module contains some subroutines for extraction of off-diagonal
!> sub-matrices from system matrices into a sparse matrix structure.
!> @sa sparsematrixmodule.

module sprExtensionModule

  use kindModule, only : dp

  implicit none

  real(dp), parameter, private :: epsZero_p = 1.0e-16_dp !< Zero tolerance

  private :: SPRKC1Sparse


contains

  !!============================================================================
  !> @brief Extracts an off-diagonal sub-matrix from a skyline system matrix.
  !>
  !> @param[in] MSKY Matrix of skyline definitions
  !> @param[in] SM System matrix elements stored in sky-line format
  !> @param[out] S12 Off-diagonal sub-matrix associated with the external DOFs
  !> @param[out] IERR Error flag
  !>
  !> @details This subroutine extracts the rectangular submatrix @a S12
  !> corresponding to 1-status DOFs row-wise and 2-status DOFs column-wise,
  !> from a symmetric system matrix @a SM stored in the 'sky-line' format.
  !> @a S12 is returned as a sparsematrixmodule::sparsematrixtype object.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Mar 2013

  subroutine GETS12sparse ( MSKY, SM, S12, IERR )

    use SparseMatrixModule, only : SparseMatrixType, smSize, smSetValue
    use ReportErrorModule , only : getErrorFile

    integer               , intent(in)  :: MSKY(:)
    real(dp)              , intent(in)  :: SM(:)
    type(SparseMatrixType), intent(out) :: S12
    integer               , intent(out) :: IERR

    !! Local variables
    integer :: I, IP, ISTART, J, NDOF1, NDOF2

    !! --- Logic section ---

    IERR  = 0
    NDOF1 = smSize(S12,1)
    NDOF2 = smSize(S12,2)
    if (NDOF1+NDOF2 < size(MSKY)) then
       call ASMERR (16,NDOF1,NDOF2,getErrorFile(),IERR)
       return
    end if

    do J = NDOF1+1, NDOF1+NDOF2
       ISTART = J - MSKY(J) + MSKY(J-1) + 1
       if (ISTART < 1 .or. ISTART > J) then
          call ASMERR (15,ISTART,J,getErrorFile(),IERR)
          return
       else if (ISTART <= NDOF1) then
          IP = MSKY(J-1) + 1
          do I = ISTART, NDOF1
             if (abs(SM(IP)) > epsZero_p) then ! Only insert the non-zero terms
                call smSetValue (S12,SM(IP),I,J-NDOF1,IERR)
                if (IERR < 0) return
             end if
             IP = IP + 1
          end do
       end if
    end do

  end subroutine GETS12sparse


  !!============================================================================
  !> @brief Extracts an off-diagonal sub-matrix from a sparse system matrix.
  !>
  !> @param[in] MSPAR Matrix of sparse parameters
  !> @param[in] MTREES Matrix of elimination assembly trees
  !> @param[in] MSIFA Matrix of storage information for FA
  !> @param[in] SM System matrix elements stored in sparse format
  !> @param[out] K12 Off-diagonal sub-matrix associated with the external DOFs
  !> @param[out] IERR Error flag
  !>
  !> @details This subroutine extracts the lower triangular sub-matrix
  !> associated with the external dofs from @a SM. If this is done prior to
  !> factorization, the results is the coefficient matrix A12.
  !>
  !> @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 8 Apr 1999

  subroutine SPRK12sparse ( MSPAR, MTREES, MSIFA, SM, K12, IERR )

    use sprKindModule     , only : ik
    use SparseMatrixModule, only : SparseMatrixType
    use ScratchArrayModule, only : getIntegerScratchArray

    integer(ik)           , intent(in)  :: MSPAR(:), MTREES(:), MSIFA(:)
    real(dp)              , intent(in)  :: SM(:)
    type(SparseMatrixType), intent(out) :: K12
    integer               , intent(out) :: IERR

    !! Local variables
    logical     :: LIPERM
    integer(ik) :: NEQ   , NSUPER, NOFSUB, SUPRET
    integer(ik) :: XSUPER, XLINDX, SUPSUP
    integer(ik) :: PERE2I, PEREND, LINDX , XLNZ

    integer, pointer :: RINDEX(:) !< Working array. Contains the row index
    !! for each equation in the system that is external, i.e., the equation.

    !! --- Logic section ---

    !! --- MSPAR
    NEQ    = MSPAR( 8)
    NSUPER = MSPAR(11)
    NOFSUB = MSPAR(15)
    SUPRET = MSPAR(53)
    LIPERM = MSPAR(58) > 0_ik

    !! --- MTREES
    XSUPER = MSPAR(45)
    XLINDX = MSPAR(46)
    SUPSUP = XLINDX + NSUPER + 1_ik

    !! --- MSIFA
    if (LIPERM) then
       PERE2I = 1_ik + NEQ
       PEREND = NEQ + NEQ
    else
       PERE2I = 1_ik
       PEREND = 1_ik
    end if
    LINDX  = MSPAR(47)
    XLNZ   = MSPAR(48)

    RINDEX => getIntegerScratchArray(int(NEQ),IERR)
    if (ierr < 0) return

    !! --- Extract K12.
    call SPRKC1sparse (LIPERM, NEQ, NSUPER, &
         &             MTREES(XSUPER:XSUPER+NSUPER), &
         &             MTREES(XLINDX:XLINDX+NSUPER), &
         &             MTREES(SUPSUP:SUPSUP+NSUPER-1_ik), &
         &             MSIFA(PERE2I:PEREND), &
         &             MSIFA(LINDX:LINDX+NOFSUB-1_ik), &
         &             MSIFA(XLNZ:XLNZ+NSUPER), &
         &             SM,K12,RINDEX,IERR)

  end subroutine SPRK12sparse


  !!============================================================================
  !> @brief Extracts an off-diagonal sub-matrix from a sparse system matrix.

  subroutine SPRKC1sparse ( LIPERM, NEQ, NSUPER, &
       &                    XSUPER, XLINDX, SUPSUP, PERE2I, LINDX, XLNZ, &
       &                    SM, K12, RINDEX, IERR)

    use sprKindModule     , only : ik
    use SparseMatrixModule, only : SparseMatrixType, smSetValue

    logical               , intent(in)  :: LIPERM
    integer(ik)           , intent(in)  :: NEQ, NSUPER
    integer(ik)           , intent(in)  :: XSUPER(:), XLINDX(:), SUPSUP(:)
    integer(ik)           , intent(in)  :: PERE2I(:), LINDX(:), XLNZ(:)
    real(dp)              , intent(in)  :: SM(:)
    type(SparseMatrixType), intent(out) :: K12
    integer               , intent(out) :: RINDEX(:), IERR

    !! Local variables
    integer(ik) :: I, J, JS, IPSM, ISTRT, ISTOP, JSTRT, JSTOP
    integer     :: COLIND, ROWIND

    !! --- Logic section ---

    IERR   = 0
    RINDEX = 0

    ! Compute RINDEX.
    do JS = 1_ik, NSUPER
       if ( SUPSUP(JS) > 0_ik ) then

          ! This supernode is in the retained set and we mark the indices.
          do I = XSUPER(JS), XSUPER(JS+1)-1_ik
             RINDEX(I) = int(I)
          end do

       end if
    end do

    COLIND = 0
    ROWIND = 0
    do J = 1, NEQ
       if ( LIPERM ) then
          I = PERE2I(J)
       else
          I = J
       end if

       if ( RINDEX(I) > 0 ) then

          ! This index belongs to retained set and is a row index in K21.
          ROWIND = ROWIND + 1
          RINDEX(I) = ROWIND

       else

          ! This index belongs to internal set and is a column index in K21.
          COLIND = COLIND + 1
          RINDEX(I) = -COLIND

       end if

    end do

    do JS = 1, NSUPER
       if ( SUPSUP(JS) > 0_ik ) cycle

       ! This supernode is in the internal set
       ! and we need to scatter its coefficients
       ! that belong to K21 to K12.
       ISTRT = XLINDX(JS)
       ISTOP = XLINDX(JS+1)-1_ik
       JSTRT = XSUPER(JS)
       JSTOP = XSUPER(JS+1)-1_ik

       ! Next position in SM.
       IPSM = XLNZ(JS)

       ! Insert the coefficient set for the current supernode.
       do J = JSTRT, JSTOP

          ! Get the column of K21. We know that it is negative
          ! and we need to negate.
          COLIND = -RINDEX(J)

          ! Now, traverse down the column of J while picking
          ! the row indices that belong to K21.
          do I = ISTRT, ISTOP

             ! --------------- Get the row index.
             ROWIND = RINDEX(LINDX(I))
             if ( ROWIND > 0 ) then

                ! This is a row index that belongs to K21,
                ! i.e. we have a column of K12 and COLIND
                ! is the corresponding row.

                if (abs(SM(IPSM)) > epsZero_p) then ! Only insert non-zero terms
                   call smSetValue (K12,SM(IPSM),COLIND,ROWIND,IERR)
                   if (IERR < 0) return
                end if

             end if

             IPSM = IPSM + 1

          end do

          ! Update index pointer.
          ISTRT = ISTRT + 1

       end do

    end do

  end subroutine SPRKC1sparse

end module sprExtensionModule
