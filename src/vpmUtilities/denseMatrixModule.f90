!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file denseMatrixModule.f90
!>
!> @brief Subroutines for solution of dense linear equation systems.

!!==============================================================================
!> @brief Module with subroutines for solution of dense linear equation systems.
!>
!> @details This module contains subroutines for solving dense linear systems of
!> equations, by wrapping subroutines from the LAPACK linear algebra package.

module DenseMatrixModule

  implicit none

  !> @brief Solves an eigenvalue problem.
  interface solveEigenvalues
     module procedure solveEigenvalueA
     module procedure solveEigenvalueB
  end interface


contains

  !!============================================================================
  !> @brief Solves the linear equation system A*x = B.
  !>
  !> @param A Dense coefficient matrix
  !> @param B Right-hand-side vectors on input, solution vectors on output
  !> @param[out] ierr Error flag
  !> @param[in]  ndim Dimension of equation system,
  !> specify only if less than @a size(A,1)
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jan 2004

  subroutine solveAxB (A,B,ierr,ndim)

    use KindModule       , only : dp
    use ReportErrorModule, only : internalError
    use ReportErrorModule, only : reportError, error_p, debugFileOnly_p

    real(dp)         , intent(inout) :: A(:,:), B(:,:)
    integer          , intent(out)   :: ierr
    integer, optional, intent(in)    :: ndim

    !! Local variables
    integer           :: n, m, ipiv(128)
    character(len=32) :: errMsg

    !! --- Logic section ---

    if (present(ndim)) then
       m = ndim
    else
       m = size(A,1)
    end if
    n = min(m,size(A,2))
    if (size(ipiv) < n) then
       ierr = internalError('solveAxB: Local array ipiv is too small')
       return
    end if

    !! Factorize matrix A
    call DGETRF (m,n,A(1,1),size(A,1),ipiv(1),ierr)
    if (ierr /= 0) then
       write(errMsg,"('LAPACK::DGETRF, INFO =',I8)") ierr
       goto 900
    end if

    !! Solve for x = A^-1*B
    call DGETRS ('N',m,size(B,2),A(1,1),size(A,1),ipiv(1),B(1,1),size(B,1),ierr)
    if (ierr /= 0) then
       write(errMsg,"('LAPACK::DGETRS, INFO =',I8)") ierr
       goto 900
    end if

    return

900 call reportError (error_p,'Cannot solve the equation system Ax=B', &
         &            addString=errMsg)
    call reportError (debugFileOnly_p,'solveAxB')

  end subroutine solveAxB


  !!============================================================================
  !> @brief Solves a symmetric-definite eigenvalue problem.
  !>
  !> @param A Dense coefficient matrix
  !> @param[out] lambda Computed eigenvalues
  !> @param[out] Z Computed eigenvectors
  !> @param[in]  n Dimension of the eigenvalue problem
  !> @param[in]  nVal Number of eigenvalues to solve for
  !> @param[in]  nVec Number of eigenvectors to solve for
  !> @param[out] ierr Error flag
  !>
  !> @details The eigenvalues and (optionally) the associated eigenvectors
  !> of the symmetric-definite eigenproblem A*x - &lambda;x = 0
  !> is computed using LAPACK::DSYEVX, where @a A is a dense symmetric matrix.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Jul 2021

  subroutine solveEigenvalueA (A,lambda,Z,n,nVal,nVec,ierr)

    use KindModule       , only : dp, nbi_p, nbd_p
    use AllocationModule , only : logAllocMem, doLogMem
    use ReportErrorModule, only : reportError, error_p, debugFileOnly_p
    use ReportErrorModule, only : AllocationError, InternalError

    integer , intent(in)    :: n, nVal, nVec
    real(dp), intent(inout) :: A(:,:)
    real(dp), intent(out)   :: lambda(:), Z(:,:)
    integer , intent(out)   :: ierr

    !! Local variables
    character(len=1)   :: jZ
    character(len=64)  :: errMsg
    integer            :: i, Lwork, nEvl
    integer , target   :: iDum(1)
    integer , pointer  :: Iwork(:), Ifail(:)
    real(dp)           :: ATOL
    real(dp), target   :: Nwork(1)
    real(dp), pointer  :: work(:)
    real(dp), external :: DLAMCH

    !! --- Logic section ---

    if (nVec > 0) then
       jZ = 'V'
    else
       jZ = 'N'
    end if
    atol = 2.0_dp*DLAMCH('S')
    Lwork = -1
    iDum = 0

    do i = 1, 2

       if (i == 1) then
          work => Nwork
          Iwork => iDum
          iFail => iDum
       else
          Lwork = INT(Nwork(1))
          allocate(work(Lwork),Iwork(5*N),Ifail(N),STAT=ierr)
          if (ierr /= 0) then
             ierr = AllocationError('solveEigenvalues: workspace for DSYEVX')
             return
          else if (doLogMem) then
             call logAllocMem ('solveEigenvalues',0,Lwork,nbd_p)
             call logAllocMem ('solveEigenvalues',0,6*N,nbi_p)
          end if
       end if

       !! Solve ([A] - lambda*[I])*{Z} = 0
       call DSYEVX (jZ,'I','U',N,A(1,1),N,work(1),work(1), &
            &       1,nVal,atol,nEvl,lambda(1),Z(1,1),N,work(1), &
            &       Lwork,Iwork(1),Ifail(1),ierr)
       if (ierr < 0) then
          write(errMsg,"('Invalid value on input argument number',i3)") -ierr
          ierr = internalError('DSYEVX: '//errMsg)
          goto (120,100) i
       else if (ierr > N) then
          write(errMsg,"('The leading minor of order',i6)") ierr-N
          call reportError (error_p,'No eigenvalues can be computed.', &
               trim(errMsg)//' of the matrix is not positive definite.', &
               addString='DSYEVX')
          goto 100
       else if (ierr > 0) then
          write(errMsg,"(i6,' eigenvectors failed to converge.')") ierr
          call reportError (error_p,errMsg,addString='DSYEVX')
       end if

    end do

100 continue
    if (doLogMem) then
       call logAllocMem ('solveEigenvalues',Lwork,0,nbd_p)
       call logAllocMem ('solveEigenvalues',6*N,0,nbi_p)
    end if
    deallocate(work,Iwork,Ifail)
120 continue
    if (ierr /= 0) call reportError (debugFileOnly_p,'solveEigenvalues')

  end subroutine solveEigenvalueA


  !!============================================================================
  !> @brief Solves a generalized symmetric-definite eigenvalue problem.
  !>
  !> @param A Dense stiffness matrix
  !> @param B Dense or diagonal mass matrix
  !> @param[out] lambda Computed eigenvalues
  !> @param[out] Z Computed eigenvectors
  !> @param[in]  n Dimension of the eigenvalue problem
  !> @param[in]  nVal Number of eigenvalues to solve for
  !> @param[in]  nVec Number of eigenvectors to solve for
  !> @param[out] ierr Error flag
  !>
  !> @details The eigenvalues and (optionally) the associated eigenvectors
  !> of the generalized symmetric-definite eigenproblem A*x - &lambda;*B*x = 0
  !> is computed using LAPACK::DSYGVX, where @a A and @a B are dense matrices.
  !> If @a B is a diagonal matrix, the eigenvalue problem is transformed into
  !> the form A*x - &lambda;*I*x = 0 where @a I is the identity matrix, and
  !> is solved using the LAPACK::DSYEVX subroutine instead.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 8 Apr 2011

  subroutine solveEigenvalueB (A,B,lambda,Z,n,nVal,nVec,ierr)

    use KindModule       , only : dp, epsDiv0_p, nbi_p, nbd_p
    use AllocationModule , only : logAllocMem, doLogMem
    use ReportErrorModule, only : reportError, error_p, debugFileOnly_p
    use ReportErrorModule, only : AllocationError, InternalError

    integer , intent(in)          :: n, nVal, nVec
    real(dp), intent(inout)       :: A(:,:), B(:,:)
    real(dp), intent(out)         :: lambda(:)
    real(dp), intent(out), target :: Z(:,:)
    integer , intent(out)         :: ierr

    !! Local variables
    character(len=1)   :: jZ
    character(len=64)  :: errMsg
    integer            :: i, Lwork, nEvl
    integer , target   :: iDum(1)
    integer , pointer  :: Iwork(:), Ifail(:)
    real(dp)           :: ATOL
    real(dp), target   :: Nwork(1), Zwork(1,1)
    real(dp), pointer  :: work(:), tZ(:,:)
    real(dp), external :: DLAMCH

    !! --- Logic section ---

    if (size(B,2) == 1) then
       !! Transform the eigenvalue problem to standard form A*x-lambda*I*x=0
       ierr = 0
       do i = 1, N
          if (B(i,1) > epsDiv0_p) then
             atol = 1.0_dp/sqrt(B(i,1))
             A(i,:) = A(i,:) * atol
             A(:,i) = A(:,i) * atol
          else
             ierr = ierr + 1
          end if
       end do
       if (ierr > 0) then
          write(errMsg,"(i6,' zero or negative diagonal elements.')") ierr
          call reportError (error_p,'The mass matrix is not positive define.', &
               'It contains '//adjustl(errMsg), &
               'Check that you have a proper mass distribution.')
       else
          call solveEigenvalueA (A,lambda,Z,n,nVal,nVec,ierr)
       end if
       if (ierr == 0) then
          do i = 1, nVal
             lambda(i) = lambda(i) * sqrt(B(i,1))
          end do
       end if
       return
    else if (size(B,2) /= size(A,2)) then
       ierr = internalError('solveEigenvalues: Incompatible matrix dimensions')
       return
    end if

    if (nVec > 0) then
       jZ = 'V'
       tZ => Z
    else
       jZ = 'N'
       tZ => Zwork
    end if
    atol = 2.0_dp*DLAMCH('S')
    Lwork = -1
    iDum = 0

    do i = 1, 2

       if (i == 1) then
          work => Nwork
          Iwork => iDum
          iFail => iDum
       else
          Lwork = INT(Nwork(1))
          allocate(work(Lwork),Iwork(5*N),Ifail(N),STAT=ierr)
          if (ierr /= 0) then
             ierr = AllocationError('solveEigenvalues: workspace for DSYGVX')
             return
          else if (doLogMem) then
             call logAllocMem ('solveEigenvalues',0,Lwork,nbd_p)
             call logAllocMem ('solveEigenvalues',0,6*N,nbi_p)
          end if
       end if

       !! Solve ([A] - lambda*[B])*{Z} = 0
       call DSYGVX (1,jZ,'I','U',N,A(1,1),N,B(1,1),N,work(1),work(1), &
            &       1,nVal,atol,nEvl,lambda(1),tZ(1,1),N,work(1), &
            &       Lwork,Iwork(1),Ifail(1),ierr)

       if (ierr < 0) then
          write(errMsg,"('Invalid value on input argument number',i3)") -ierr
          ierr = internalError('DSYGVX: '//errMsg)
          goto (120,100) i
       else if (ierr > N) then
          write(errMsg,"('The leading minor of order',i6)") ierr-N
          call reportError (error_p,'No eigenvalues can be computed.', &
               trim(errMsg)//' of the mass matrix is not positive definite.', &
               'Check that you have a proper mass distribution.', &
               addString='DSYGVX')
          goto 100
       else if (ierr > 0) then
          write(errMsg,"(i6,' eigenvectors failed to converge.')") ierr
          call reportError (error_p,errMsg, &
               'Check that the FE model is consistent.', &
               addString='DSYGVX')
       end if

    end do

100 continue
    if (doLogMem) then
       call logAllocMem ('solveEigenvalues',Lwork,0,nbd_p)
       call logAllocMem ('solveEigenvalues',6*N,0,nbi_p)
    end if
    deallocate(work,Iwork,Ifail)
120 continue
    if (ierr /= 0) call reportError (debugFileOnly_p,'solveEigenvalues')

  end subroutine solveEigenvalueB

end module DenseMatrixModule
