!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file matExtensionModule.f90
!> @brief System matrix extension module.

!!==============================================================================
!> @brief System matrix extension module.
!>
!> @details This module contains subroutines for performing various operations
!> on system matrices, such as matrix-matrix multiplication, matrix-vector
!> multiplications, and extraction of sub-matrices from the system matrix.
!> The module defines an interface between the FEDEM solver applications
!> and the various linear algebra packages that are being used.

module MatExtensionModule

  implicit none

  private

  !> @brief Performs a triple matrix multiplication.
  interface csTransform
     module procedure csTransformMat
     module procedure csTransformVec
  end interface

  !> @brief Performs a matrix-matrix multiplication.
  interface csPremult
     module procedure csPremultMat
     module procedure csPremultVec
     module procedure csPremultMatrix
     module procedure csPremultVector
  end interface

  !> @brief Extracts and off-diagonal submatrix from the system matrix.
  interface csGetSub12
     module procedure csGetSub12dense
     module procedure csGetSub12sparse
  end interface

  public :: csTransform, csPremult, csAddMat, csCopyMat
  public :: csGetSub11, csGetSub12, csGetSub22, csGetDiag


contains

  !!============================================================================
  !> @brief Performs a triple matrix multiplication.
  !>
  !> @param[in] aMat System matrix @b A
  !> @param[in] tMat Transformation matrix @b T
  !> @param[out] bMat Resulting matrix @b B
  !> @param[out] ierr Error flag
  !>
  !> @details The matrix multiplication @b B = @b T <sup>T</sup> @b A @b T
  !> is performed, where @b A is on compact (sparse) form and @b T is a dense
  !> rectangular matrix.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csTransformMat (aMat,tMat,bMat,ierr)

    use SprKindModule       , only : ik, dp
    use SysMatrixTypeModule , only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use FELinearSolver      , only : FETRA
    use manipMatrixModule   , only : diagTransform
    use scratchArrayModule  , only : getRealScratchArray
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p
    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType), intent(inout) :: aMat
    real(dp)           , intent(in)    :: tMat(:,:)
    real(dp)           , intent(out)   :: bMat(:,:)
    integer            , intent(out)   :: ierr

    !! Local variables
    character(len=1)  :: cOpt
    integer           :: sizeA, sizeB
    integer(ik)       :: lerr
    real(dp), pointer :: rWork(:)

    character(len=31), parameter :: prgnam_p = 'MatExtensionModule::csTransform'

    !! --- Logic section ---

    sizeA = size(tMat,1)
    sizeB = size(tMat,2)
    if (sizeA > aMat%dim) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions 1')
       return
    else if (sizeA /= aMat%dim .and. aMat%storageType /= outOfCore_p) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions 2')
       return
    else if (sizeB /= size(bMat,1) .or. sizeB /= size(bMat,2)) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions 3')
       return
    end if

    select case (aMat%storageType)

    case (diagonalMatrix_p)

       call diagTransform (aMat%value,tMat,bMat,ierr)

    case (denseMatrix_p)

       rWork => getRealScratchArray(sizeA*sizeB,ierr)
       if (ierr == 0) then

          call DGEMM ('N','N',sizeA,sizeB,sizeA,1.0_dp,aMat%value(1),aMat%dim, &
               &      tMat(1,1),sizeA,0.0_dp,rWork(1),sizeA)
          call DGEMM ('T','N',sizeB,sizeA,sizeA,1.0_dp,tMat(1,1),sizeA, &
               &      rWork(1),sizeA,0.0_dp,bMat(1,1),sizeB)

       end if

    case (skylineMatrix_p)

       rWork => getRealScratchArray(sizeA,ierr)
       if (ierr == 0) then

          call SKYTRA (aMat%value(1), tMat(1,1), bMat(1,1), &
               &       rWork(1), aMat%skyline%msky(1), &
               &       sizeA, sizeB, 1, getErrorFile(), ierr)

       end if

    case (sparseMatrix_p)

       rWork => getRealScratchArray(sizeA,ierr)
       if (ierr == 0) then

          call SPRTRA (aMat%value(sizeA+1), tMat(1,1), bMat(1,1), rWork(1), &
               &       aMat%sparse%mspar(1), aMat%sparse%mtrees(1), &
               &       aMat%sparse%msifa(1), int(sizeA,ik), int(sizeB,ik), &
               &       1_ik, int(getErrorFile(),ik), lerr)
          ierr = int(lerr)

       end if

    case (outOfCore_p)

       if (sizeA < aMat%dim) then
          cOpt = 'I'
       else
          cOpt = 'A'
       end if
       call ffa_startTimer ('FETRA')
       call FETRA (cOpt, 'P', sizeA, sizeB, aMat%gsf%Q, &
            &      1.0_dp, aMat%gsf%S, tMat, 0.0_dp, bMat, GSFinfo)
       call ffa_stopTimer ('FETRA')
       ierr = checkGSFinfo(GSFinfo)

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    if (ierr < 0) then
       call reportError (error_p,'Failed to transform system '// &
            &            matrixType_p(aMat%storageType),addString=prgnam_p)
    end if

  end subroutine csTransformMat


  !!============================================================================
  !> @brief Performs a vector-matrix-vector multiplication.
  !>
  !> @param[in] aMat System matrix @b A
  !> @param[in] v System vector
  !> @param[out] bVal Resulting value
  !> @param[out] ierr Error flag
  !>
  !> @details The multiplication @a b = @b v <sup>T</sup> @b A @b v
  !> is performed, where @b A is on compact (sparse) form and @b v is a vector.
  !> The result @a b is s scalar quantity.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jul 2003

  subroutine csTransformVec (aMat,v,bVal,ierr)

    use SprKindModule       , only : ik, dp
    use SysMatrixTypeModule , only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use FELinearSolver      , only : FETRA
    use manipMatrixModule   , only : MatrixPointToArray_real
    use scratchArrayModule  , only : getRealScratchArray
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p
    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType), intent(inout) :: aMat
    real(dp), target   , intent(in)    :: v(:)
    real(dp)           , intent(out)   :: bVal
    integer            , intent(out)   :: ierr

    !! Local variables
    character(len=1)  :: cOpt
    integer           :: sizeV
    integer(ik)       :: lerr
    real(dp), pointer :: rWork(:), vec(:), vMat(:,:)
    real(dp)          :: bMat(1,1)

    character(len=31), parameter :: prgnam_p = 'MatExtensionModule::csTransform'

    !! --- Logic section ---

    sizeV = size(v)
    if (sizeV > aMat%dim) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions')
       return
    end if

    select case (aMat%storageType)

    case (diagonalMatrix_p)

       ierr = 0
       bVal = dot_product(v*aMat%value(1:sizeV),v)

    case (denseMatrix_p)

       rWork => getRealScratchArray(sizeV,ierr)
       if (ierr == 0) then

          call DGEMV ('N',sizeV,sizeV,1.0_dp,aMat%value(1),aMat%dim, &
               &      v(1),1,0.0_dp,rWork(1),1)
          bVal = dot_product(v,rWork)

       end if

    case (skylineMatrix_p)

       rWork => getRealScratchArray(sizeV,ierr)
       if (ierr == 0) then

          call SKYTRA (aMat%value(1), v(1), bVal, &
               &       rWork(1), aMat%skyline%msky(1), &
               &       sizeV, 1, 1, getErrorFile(), ierr)

       end if

    case (sparseMatrix_p)

       rWork => getRealScratchArray(2*aMat%dim,ierr)
       if (ierr == 0) then

          if (sizeV < aMat%dim) then
             !! Assume that the equations corresponding to the sizeV DOFs in
             !! input vector v is ordered first and zero out the following eqs.
             vec => rWork(aMat%dim+1:)
             call DCOPY (sizeV,v(1),1,vec(1),1)
             call DCOPY (aMat%dim-sizeV,0.0_dp,0,vec(sizeV+1),1)
          else
             vec => v
          end if
          call SPRTRA (aMat%value(aMat%dim+1), vec(1), bVal, rWork(1), &
               &       aMat%sparse%mspar(1), aMat%sparse%mtrees(1), &
               &       aMat%sparse%msifa(1), int(aMat%dim,ik), 1_ik, &
               &       1_ik, int(getErrorFile(),ik), lerr)
          ierr = int(lerr)

       end if

    case (outOfCore_p)

       if (sizeV < aMat%dim) then
          cOpt = 'I'
       else
          cOpt = 'A'
       end if
       vMat => MatrixPointToArray_real(v,sizeV,1)
       call ffa_startTimer ('FETRA')
       call FETRA (cOpt, 'P', sizeV, 1, aMat%gsf%Q, 1.0_dp, aMat%gsf%S, &
            &      vMat, 0.0_dp, bMat, GSFinfo)
       call ffa_stopTimer ('FETRA')
       ierr = checkGSFinfo(GSFinfo)
       bVal = bMat(1,1)

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    if (ierr < 0) then
       call reportError (error_p,'Failed to transform system '// &
            &            matrixType_p(aMat%storageType),addString=prgnam_p)
    end if

  end subroutine csTransformVec


  !!============================================================================
  !> @brief Performs a matrix-matrix multiplication.
  !>
  !> @param[in] aMat System matrix @b A
  !> @param[in] bMat Rectangular dense matrix @b B
  !> @param[out] cMat Resulting matrix @b C
  !> @param[in] iflag Operation flag (see below)
  !> @param[out] ierr Error flag
  !>
  !> @details The matrix multiplication @b C = @b A @b B is performed,
  !> where @b A is on compact (sparse) form whereas @b C and @b B are dense
  !> rectangular matrices. The value on @a iflag determines the sign on the
  !> resulting matrix, as follows:
  !> - iflag =  0 : @b C = @b A @b B (where @b B and @b C may be the same array)
  !> - iflag = +1 : @b C = @b A @b B (@b B and @b C cannot be the same array)
  !> - iflag = -1 : @b C = -@b A @b B
  !> - iflag = +2 : @b C = @b C + @b A @b B
  !> - iflag = -2 : @b C = @b C - @b A @b B
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csPremultMat (aMat,bMat,cMat,iflag,ierr)

    use SprKindModule       , only : ik, dp
    use SysMatrixTypeModule , only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use FELinearSolver      , only : FEPRM
    use scratchArrayModule  , only : getRealScratchArray, getRealScratchMatrix
    use manipMatrixModule   , only : diagmatMul
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p
    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType), intent(inout) :: aMat
    real(dp)           , intent(in)    :: bMat(:,:)
    real(dp), target   , intent(out)   :: cMat(:,:)
    integer            , intent(in)    :: iflag
    integer            , intent(out)   :: ierr

    !! Local variables
    character(len=1)  :: cOpt
    integer           :: sizeA, sizeB, sizeC, jflag
    integer(ik)       :: lerr
    real(dp)          :: alpha, beta
    real(dp), pointer :: rWork(:), rMat(:,:)

    character(len=29), parameter :: prgnam_p = 'MatExtensionModule::csPremult'

    !! --- Logic section ---

    if (iflag == 0) then
       jflag = 1
    else
       jflag = iflag
    end if
    sizeA = size(cMat,1)
    sizeB = size(bMat,1)
    sizeC = size(cMat,2)

    if (sizeA > aMat%dim .or. sizeB > aMat%dim) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions 1')
       return
    else if ((sizeA /= sizeB .or. sizeA /= aMat%dim) .and. &
         &   (aMat%storageType == skylineMatrix_p .or. &
         &    aMat%storageType == sparseMatrix_p)) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions 2')
       return
    else if (sizeC /= size(bMat,2)) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions 3')
       return
    end if

    select case (aMat%storageType)

    case (diagonalMatrix_p)

       call diagMatMul (aMat%value, bMat, cMat, jflag, ierr)

    case (denseMatrix_p)

       if (iflag == 0) then
          alpha = 1.0_dp
          beta = 0.0_dp
          rMat => getRealScratchMatrix(sizeA,sizeC,ierr)
       else
          alpha = real(sign(1,iflag),dp)
          beta = real(abs(iflag)-1,dp)
          rMat => cMat
       end if
       if (ierr == 0) then

          call DGEMM ('N','N',sizeA,sizeC,sizeB,alpha,aMat%value(1),aMat%dim, &
               &      bMat(1,1),sizeB,beta,rMat(1,1),sizeA)
          if (iflag == 0) then
             call DCOPY (sizeA*sizeC,rMat(1,1),1,cMat(1,1),1)
          end if

       end if

    case (skylineMatrix_p)

       rWork => getRealScratchArray(aMat%dim,ierr)
       if (ierr == 0) then

          call PRESKY (aMat%value(1),  bMat(1,1), cMat(1,1), rWork(1), &
               &       aMat%skyline%msky(1), sizeA, sizeC, 1, jflag, &
               &       getErrorFile(), ierr)

       end if

    case (sparseMatrix_p)

       rWork => getRealScratchArray(aMat%dim,ierr)
       if (ierr == 0) then

          call SPRPRM (aMat%value(sizeA+1), bMat(1,1), cMat(1,1), rWork(1), &
               &       aMat%sparse%mspar(1), aMat%sparse%mtrees(1), &
               &       aMat%sparse%msifa(1), int(sizeA,ik), int(sizeC,ik), &
               &       1_ik, int(jflag,ik), int(getErrorFile(),ik), lerr)
          ierr = int(lerr)

       end if

    case (outOfCore_p)

       if (sizeB < aMat%dim) then
          cOpt = 'I'
       else
          cOpt = 'A'
       end if
       alpha = real(sign(1,jflag))
       beta = real(abs(jflag)-1,dp)
       call ffa_startTimer ('FEPRM')
       call FEPRM (cOpt, 'P', sizeB, sizeC, aMat%gsf%Q, &
            &      alpha, aMat%gsf%S, bMat, beta, cMat, GSFinfo)
       call ffa_stopTimer ('FEPRM')
       ierr = checkGSFinfo(GSFinfo)

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    if (ierr < 0) then
       call reportError (error_p,'Failed to premultiply with system '// &
            &            matrixType_p(aMat%storageType),addString=prgnam_p)
    end if

  end subroutine csPremultMat


  !!============================================================================
  !> @brief Performs a matrix-matrix multiplication.
  !>
  !> @param[in] aMat System matrix @b A
  !> @param bMat Rectangular dense matrix @b B
  !> @param[out] ierr Error flag
  !>
  !> @details The matrix multiplication @b B = @b A @b B is performed,
  !> where @b A is on compact (sparse) form and @b B is a dense
  !> rectangular matrix.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jan 2018

  subroutine csPremultMatrix (aMat,bMat,ierr)

    use SysMatrixTypeModule, only : SysMatrixType, dp

    type(SysMatrixType), intent(inout) :: aMat
    real(dp), target   , intent(inout) :: bMat(:,:)
    integer            , intent(out)   :: ierr

    !! Local variables
    real(dp), pointer :: cMat(:,:)

    !! --- Logic section ---

    cMat => bMat
    call csPremultMat (aMat,bMat,cMat,0,ierr)

  end subroutine csPremultMatrix


  !!============================================================================
  !> @brief Performs a matrix-vector multiplication.
  !>
  !> @param[in] aMat System matrix @b A
  !> @param[in] bVec System vector @b b
  !> @param[out] cVec Resulting matrix @b c
  !> @param[in] iflag Operation flag (see below)
  !> @param[out] ierr Error flag
  !>
  !> @details The matrix-vector multiplication @b c = @b A @b b is performed,
  !> where @b A is on compact (sparse) form whereas @b c and @b b are vectors.
  !> The value on @a iflag determines the sign on the resulting vector,
  !> as follows:
  !> - iflag = +1 : @b c = @b A @b b
  !> - iflag = -1 : @b c = -@b A @b b
  !> - iflag = +2 : @b c = @b c + @b A @b b
  !> - iflag = -2 : @b c = @b c - @b A @b b
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csPremultVec (aMat,bVec,cVec,iflag,ierr)

    use SprKindModule       , only : ik, dp
    use SysMatrixTypeModule , only : SysMatrixType, dp, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use FELinearSolver      , only : FEPRM
    use scratchArrayModule  , only : getRealScratchArray
    use manipMatrixModule   , only : diagmatMul, MatrixPointToArray_real
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p
    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType), intent(inout) :: aMat
    real(dp)           , intent(in)    :: bVec(:)
    real(dp), target   , intent(out)   :: cVec(:)
    integer            , intent(in)    :: iflag
    integer            , intent(out)   :: ierr

    !! Local variables
    character(len=1)  :: cOpt
    integer           :: sizeA, sizeB, jflag
    integer(ik)       :: lerr
    real(dp)          :: alpha, beta
    real(dp), pointer :: rWork(:), bMat(:,:), cMat(:,:)

    character(len=29), parameter :: prgnam_p = 'MatExtensionModule::csPremult'

    !! --- Logic section ---

    ierr  = 0
    jflag = iflag
    sizeA = aMat%dim
    sizeB = size(bVec)
    rWork => cVec
    if (sizeB /= size(cVec) .or. sizeA < sizeB) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions 1')
       return
    else if (aMat%storageType == skylineMatrix_p .or. &
         &   aMat%storageType == sparseMatrix_p) then
       if (sizeB /= sizeA) then
          ierr = internalError(prgnam_p//': Invalid matrix dimensions 2')
          return
       else if (iflag == 1) then
          jflag = 88
       else
          rWork => getRealScratchArray(sizeB,ierr)
       end if
    end if

    select case (aMat%storageType)

    case (diagonalMatrix_p)

       call diagMatMul (aMat%value, bVec, cVec, iflag, ierr)

    case (denseMatrix_p)

       alpha = real(sign(1,iflag),dp)
       beta = real(abs(iflag)-1,dp)
       call DGEMV ('N',sizeB,sizeB,alpha,aMat%value(1),sizeA, &
            &      bVec(1),1,beta,cVec(1),1)

    case (skylineMatrix_p)

       if (ierr == 0) then
          call PRESKY (aMat%value(1), bVec(1), cVec(1), rWork(1), &
               &       aMat%skyline%msky(1), sizeA, 1, 1, jflag, &
               &       getErrorFile(), ierr)
       end if

    case (sparseMatrix_p)

       if (ierr == 0) then
          call SPRPRM (aMat%value(sizeA+1), bVec(1), cVec(1), rWork(1), &
               &       aMat%sparse%mspar(1), aMat%sparse%mtrees(1), &
               &       aMat%sparse%msifa(1), int(sizeA,ik), 1_ik, &
               &       1_ik, int(jflag,ik), int(getErrorFile(),ik), lerr)
          ierr = int(lerr)
       end if

    case (outOfCore_p)

       if (sizeB < sizeA) then
          cOpt = 'I'
       else
          cOpt = 'A'
       end if
       alpha = real(sign(1,iflag),dp)
       beta = real(abs(iflag)-1,dp)
       bMat => MatrixPointToArray_real(bVec,sizeB,1)
       cMat => MatrixPointToArray_real(cVec,sizeB,1)
       call ffa_startTimer ('FEPRM')
       call FEPRM (cOpt, 'P', sizeB, 1, aMat%gsf%Q, alpha, aMat%gsf%S, &
            &      bMat, beta, cMat, GSFinfo)
       call ffa_stopTimer ('FEPRM')
       ierr = checkGSFinfo(GSFinfo)

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    if (ierr < 0) then
       call reportError (error_p,'Failed to premultiply with system '// &
            &            matrixType_p(aMat%storageType),addString=prgnam_p)
    end if

  end subroutine csPremultVec


  !!============================================================================
  !> @brief Performs a matrix-vector multiplication.
  !>
  !> @param[in] aMat System matrix @b A
  !> @param bVec System vector @b b
  !> @param[out] ierr Error flag
  !>
  !> @details The matrix-vector multiplication @b b = @b A @b b is performed,
  !> where @b A is on compact (sparse) form and @b b is a vector.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jan 2018

  subroutine csPremultVector (aMat,bVec,ierr)

    use SysMatrixTypeModule, only : SysMatrixType, dp

    type(SysMatrixType), intent(inout) :: aMat
    real(dp), target   , intent(inout) :: bVec(:)
    integer            , intent(out)   :: ierr

    !! Local variables
    real(dp), pointer :: cVec(:)

    !! --- Logic section ---

    cVec => bVec
    call csPremultVec (aMat,bVec,cVec,1,ierr)

  end subroutine csPremultVector


  !!============================================================================
  !> @brief Performs the matrix operation @b A = @b A + &alpha; @b B.
  !>
  !> @param aMat System matrix @b A on compact (sparse) form
  !> @param[in] bMat System matrix @b B on compact (sparse) form
  !> @param[in] alpha Scaling factor &alpha;
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date May 2009

  subroutine csAddMat (aMat,bMat,alpha,ierr)

    use SysMatrixTypeModule , only : SysMatrixType, dp, matrixType_p
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use reportErrorModule   , only : internalError, reportError, error_p

    type(SysMatrixType), intent(inout) :: aMat
    type(SysMatrixType), intent(in)    :: bMat
    real(dp)           , intent(in)    :: alpha
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i, j

    character(len=29), parameter :: prgnam_p = 'MatExtensionModule::csAddMat'

    !! --- Logic section ---

    if (aMat%dim /= bMat%dim) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions')
       return
    else
       ierr = 0
    end if

    select case (aMat%storageType)

    case (diagonalMatrix_p)

       if (bMat%storageType == diagonalMatrix_p) then
          call DAXPY (aMat%dim,alpha,bMat%value(1),1,aMat%value(1),1)
          return
       end if

    case (denseMatrix_p)

       if (bMat%storageType == denseMatrix_p) then
          call DAXPY (aMat%dim*aMat%dim,alpha,bMat%value(1),1,aMat%value(1),1)
          return
       else if (bMat%storageType == diagonalMatrix_p) then
          call DAXPY (aMat%dim,alpha,bMat%value(1),1,aMat%value(1),1+aMat%dim)
          return
       end if

    case (skylineMatrix_p)

       if (bMat%storageType == skylineMatrix_p) then
          if (associated(bMat%skyline,aMat%skyline)) then
             call DAXPY (size(aMat%value),alpha,bMat%value(1),1,aMat%value(1),1)
             return
          end if
       else if (bMat%storageType == diagonalMatrix_p) then
          do i = 1, aMat%dim
             j = aMat%skyline%msky(j)
             aMat%value(j) = aMat%value(j) + alpha*bMat%value(i)
          end do
          return
       end if

    case (sparseMatrix_p)

       i = 1 + aMat%dim
       if (bMat%storageType == sparseMatrix_p) then
          if (associated(bMat%sparse%msica,aMat%sparse%msica)) then
             j = int(aMat%sparse%mspar(16))
             call DAXPY (j,alpha,bMat%value(i),1,aMat%value(i),1)
             return
          end if
       else if (bMat%storageType == diagonalMatrix_p) then
          call DAXPY (aMat%dim,alpha,bMat%value(1),1,aMat%value(i),1)
          return
       end if

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    ierr = -1
    call reportError (error_p,'Failed to add system '// &
         &            trim(matrixType_p(bMat%storageType))//' to system '// &
         &            matrixType_p(aMat%storageType),addString=prgnam_p)

  end subroutine csAddMat


  !!============================================================================
  !> @brief Performs the matrix operation @b A = @b B.
  !>
  !> @param aMat System matrix @b A on compact (sparse) form
  !> @param[in] bMat System matrix @b B on compact (sparse) form
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date May 2009

  subroutine csCopyMat (aMat,bMat,ierr)

    use SysMatrixTypeModule , only : SysMatrixType, matrixType_p
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use reportErrorModule   , only : internalError, reportError, error_p

    type(SysMatrixType), intent(inout) :: aMat
    type(SysMatrixType), intent(in)    :: bMat
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: imat, nmat

    character(len=29), parameter :: prgnam_p = 'MatExtensionModule::csCopyMat'

    !! --- Logic section ---

    if (aMat%dim /= bMat%dim) then
       ierr = internalError(prgnam_p//': Invalid matrix dimensions')
       return
    else if (aMat%storageType /= bMat%storageType) then
       ierr = internalError(prgnam_p//': Different matrix types')
       return
    end if

    ierr = 0
    imat = 1
    nmat = 0
    select case (aMat%storageType)

    case (diagonalMatrix_p)

       nmat = aMat%dim

    case (denseMatrix_p)

       nmat = aMat%dim*aMat%dim

    case (skylineMatrix_p)

       if (associated(bMat%skyline,aMat%skyline)) nmat = size(aMat%value)

    case (sparseMatrix_p)

       if (associated(bMat%sparse%msica,aMat%sparse%msica)) then
          imat = 1 + aMat%dim
          nmat = int(aMat%sparse%mspar(16))
       end if

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    if (nmat > 0) then
       call DCOPY (nmat,bMat%value(imat),1,aMat%value(imat),1)
    else
       ierr = -1
       call reportError (error_p,'Failed to copy system '// &
            &            matrixType_p(aMat%storageType),addString=prgnam_p)
    end if

  end subroutine csCopyMat


  !!============================================================================
  !> @brief Extracts the 11 sub-matrix from a system matrix.
  !>
  !> @param[in] sysMat The system matrix to extract the sub-matrix from
  !> @param[in] meqn1 Matrix of status 1 equation numbers
  !> @param[out] subMat11 Extracted sub-matrix
  !> @param[out] ierr Error flag
  !>
  !> @details The system matrix @a sysMat is stored on compact (sparse) form,
  !> but is logically composed of four submatrices, like
  !>
  !> @verbatim
  !>    sysMat = | sub11 sub12 |
  !>             | sub21 sub22 |
  !> @endverbatim
  !>
  !> The sub-matrix @a sub11 is returned as a dense square matrix.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Apr 2011

  subroutine csGetSub11 (sysMat,meqn1,subMat11,ierr)

    use SprKindModule       , only : ik, i4, dp
    use SysMatrixTypeModule , only : SysMatrixType
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use scratchArrayModule  , only : getIntegerScratchArray
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p

    type(SysMatrixType), intent(inout) :: sysMat
    integer            , intent(in)    :: meqn1(:)
    real(dp)           , intent(out)   :: subMat11(:,:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer          :: j, ndof1
    integer, pointer :: iWork(:)

    character(len=30), parameter :: prgnam_p = 'MatExtensionModule::csGetSub11'

    !! --- Logic section ---

    ndof1 = size(subMat11,1)
    if (size(subMat11,2) == ndof1) then
       ierr = 0
    else
       ierr = internalError(prgnam_p//': Invalid matrix dimension (no square)')
       return
    end if

    select case (sysMat%storageType)

    case (diagonalMatrix_p)

       if (size(meqn1) >= ndof1) then
          call DGATHR (ndof1,meqn1(1),sysMat%value(1),subMat11(1,1),1)
       else if (size(meqn1) == 0) then
          call DCOPY (ndof1,sysMat%value(1),1,subMat11(1,1),1)
       else
          ierr = internalError(prgnam_p//': Invalid size of array meqn1')
       end if

    case (denseMatrix_p)

       if (size(meqn1) >= ndof1) then
          do j = 1, ndof1
             call DGATHR (ndof1,meqn1(1), &
                  &       sysMat%value(1+sysMat%dim*(meqn1(j)-1)), &
                  &       subMat11(1,j),1)
          end do
       else if (size(meqn1) == 0) then
          do j = 1, ndof1
             call DCOPY (ndof1,sysMat%value(1+sysMat%dim*(j-1)),1, &
                  &      subMat11(1,j),1)
          end do
       else
          ierr = internalError(prgnam_p//': Invalid size of array meqn1')
       end if

    case (skylineMatrix_p)

       call EXSEM (sysMat%value(1), sysMat%skyline%msky(1), &
            &      -sysMat%dim, ndof1, getErrorFile(), subMat11(1,1), ierr)
       if (ierr < 0) goto 999

    case (sparseMatrix_p)

       if (ik > i4) then
          iWork => getIntegerScratchArray(sysMat%dim*2,ierr)
       else
          iWork => getIntegerScratchArray(sysMat%dim,ierr)
       end if
       if (ierr < 0) goto 999

       call SPRK11 (sysMat%sparse%mspar(1), sysMat%sparse%mtrees(1), &
            &       sysMat%sparse%msifa(1), sysMat%value(sysMat%dim+1), 2_ik, &
            &       iWork(1), subMat11(1,1), iWork(1))

    case (outOfCore_p)

       ierr = internalError(prgnam_p//': FEExtractK11 not implemented yet')

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')

    end select

    return

999 continue
    call reportError (error_p,'Failed to extract sub-matrix from system '// &
         &            matrixType_p(sysMat%storageType),addString=prgnam_p)

  end subroutine csGetSub11


  !!============================================================================
  !> @brief Extracts the 12 sub-matrix from a system matrix.
  !>
  !> @param[in] sysMat The system matrix to extract the sub-matrix from
  !> @param[in] meqn1 Matrix of status 1 equation numbers
  !> @param[in] meqn2 Matrix of status 2 equation numbers
  !> @param[out] subMat12 Extracted sub-matrix
  !> @param[out] ierr Error flag
  !>
  !> @details The system matrix @a sysMat is stored on compact (sparse) form,
  !> but is logically composed of four submatrices, like
  !>
  !> @verbatim
  !>    sysMat = | sub11 sub12 |
  !>             | sub21 sub22 |
  !> @endverbatim
  !>
  !> The sub-matrix @a sub12 is returned as a dense rectangular matrix.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csGetSub12dense (sysMat,meqn1,meqn2,subMat12,ierr)

    use SprKindModule       , only : ik, i4, dp
    use SysMatrixTypeModule , only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use FELinearSolver      , only : FEExtractK12
    use scratchArrayModule  , only : getIntegerScratchArray
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p
    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType), intent(inout) :: sysMat
    integer            , intent(in)    :: meqn1(:), meqn2(:)
    real(dp)           , intent(out)   :: subMat12(:,:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer          :: j, neq, ndof1, ndof2
    integer, pointer :: iWork(:)

    character(len=30), parameter :: prgnam_p = 'MatExtensionModule::csGetSub12'

    !! --- Logic section ---

    select case (sysMat%storageType)

    case (diagonalMatrix_p)

       ierr = 0
       call DCOPY (size(subMat12),0.0_dp,0,subMat12(1,1),1)

    case (denseMatrix_p)

       ierr = 0
       ndof1 = size(subMat12,1)
       ndof2 = size(subMat12,2)
       if (ndof1 > size(meqn1) .or. ndof2 > size(meqn2)) then
          ierr = internalError(prgnam_p//': Invalid matrix dimension')
          return
       end if

       neq = sysMat%dim
       do j = 1, ndof2
          call DGATHR (ndof1,meqn1(1),sysMat%value(1+neq*(meqn2(j)-1)), &
               &       subMat12(1,j),1)
       end do

    case (skylineMatrix_p)

       neq   = sysMat%dim
       ndof1 = size(subMat12,1)
       ndof2 = size(subMat12,2)

       call GETS12 (sysMat%value(1), sysmat%skyline%msky(1), &
            &       ndof1, ndof2, neq, getErrorFile(), subMat12(1,1), ierr)

    case (sparseMatrix_p)

       neq = sysMat%dim
       if (ik > i4) then
          iWork => getIntegerScratchArray(neq*2,ierr)
       else
          iWork => getIntegerScratchArray(neq,ierr)
       end if
       if (ierr < 0) goto 999

       call SPRK12 (sysMat%sparse%mspar(1), sysMat%sparse%mtrees(1), &
            &       sysMat%sparse%msifa(1), sysMat%value(neq+1), &
            &       subMat12(1,1), iWork(1))

    case (outOfCore_p)

       iWork => getIntegerScratchArray(sysMat%dim,ierr)
       if (ierr < 0) goto 999

       call ffa_startTimer ('FEExtractK12')
       call FEExtractK12 (sysMat%gsf%Q, sysMat%gsf%S, subMat12, iWork, GSFinfo)
       call ffa_stopTimer ('FEExtractK12')
       ierr = checkGSFinfo(GSFinfo)

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    if (ierr >= 0) return
999 call reportError (error_p,'Failed to extract 12 submatrix from system '// &
         &            matrixType_p(sysMat%storageType),addString=prgnam_p)

  end subroutine csGetSub12dense


  !!============================================================================
  !> @brief Extracts the 12 sub-matrix from a system matrix.
  !>
  !> @param[in] sysMat The system matrix to extract the sub-matrix from
  !> @param[out] subMat12 Extracted sub-matrix
  !> @param[out] ierr Error flag
  !>
  !> @details The system matrix @a sysMat is stored on compact (sparse) form,
  !> but is logically composed of four submatrices, like
  !>
  !> @verbatim
  !>    sysMat = | sub11 sub12 |
  !>             | sub21 sub22 |
  !> @endverbatim
  !>
  !> The sub-matrix @a sub12 is returned as a sparse matrix.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Jan 1999
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csGetSub12sparse (sysMat,subMat12,ierr)

    use SysMatrixTypeModule , only : SysMatrixType, dp, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use SparseMatrixModule  , only : SparseMatrixType, smSetMatrix
    use sprExtensionModule  , only : GETS12sparse, SPRK12sparse
    use FELinearSolver      , only : FEExtractK12
    use scratchArrayModule  , only : getIntegerScratchArray
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p
    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType)   , intent(inout) :: sysMat
    type(SparseMatrixType), intent(out)   :: subMat12
    integer               , intent(out)   :: ierr

    !! Local variables
    real(dp), pointer :: values(:)
    integer , pointer :: iRow(:), jCol(:), iWork(:)

    character(len=30), parameter :: prgnam_p = 'MatExtensionModule::csGetSub12'

    !! --- Logic section ---

    select case (sysMat%storageType)

    case (skylineMatrix_p)

       call GETS12sparse (sysMat%skyline%msky, sysMat%value, subMat12, ierr)

    case (sparseMatrix_p)

       call SPRK12sparse (sysMat%sparse%mspar, sysMat%sparse%mtrees, &
            &             sysMat%sparse%msifa, sysMat%value(sysMat%dim+1:), &
            &             subMat12, ierr)

    case (outOfCore_p)

       iWork => getIntegerScratchArray(sysMat%dim,ierr)
       if (ierr < 0) goto 999

       call ffa_startTimer ('FEExtractK12')
       call FEExtractK12 (sysMat%gsf%Q, sysMat%gsf%S, iRow, jCol, values, &
            &             iWork, GSFinfo)
       call ffa_stopTimer ('FEExtractK12')
       ierr = checkGSFinfo(GSFinfo)
       if (ierr < 0) goto 999

       call smSetMatrix (subMat12, iRow, jCol, values, ierr)

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')
       return

    end select

    if (ierr >= 0) return
999 call reportError (error_p,'Failed to extract sub-matrix from system '// &
         &            matrixType_p(sysMat%storageType),addString=prgnam_p)

  end subroutine csGetSub12sparse


  !!============================================================================
  !> @brief Extracts the 22 sub-matrix from a system matrix.
  !>
  !> @param[in] sysMat The system matrix to extract the sub-matrix from
  !> @param[in] meqn2 Matrix of status 2 equation numbers
  !> @param[out] subMat22 Extracted sub-matrix
  !> @param[out] ierr Error flag
  !>
  !> @details The system matrix @a sysMat is stored on compact (sparse) form,
  !> but is logically composed of four submatrices, like
  !>
  !> @verbatim
  !>    sysMat = | sub11 sub12 |
  !>             | sub21 sub22 |
  !> @endverbatim
  !>
  !> The sub-matrix @a sub22 is returned as a dense square matrix.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Jan 1999
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csGetSub22 (sysMat,meqn2,subMat22,ierr)

    use SprKindModule       , only : ik, i4, dp
    use SysMatrixTypeModule , only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule , only : diagonalMatrix_p, denseMatrix_p
    use SysMatrixTypeModule , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule , only : outOfCore_p, matrixType_p
    use FELinearSolver      , only : FEExtractK22
    use scratchArrayModule  , only : getIntegerScratchArray
    use reportErrorModule   , only : internalError, getErrorFile
    use reportErrorModule   , only : reportError, error_p
    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType), intent(inout) :: sysMat
    integer            , intent(in)    :: meqn2(:)
    real(dp)           , intent(out)   :: subMat22(:,:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer          :: j, ndof2
    integer, pointer :: iWork(:)

    character(len=30), parameter :: prgnam_p = 'MatExtensionModule::csGetSub22'

    !! --- Logic section ---

    ndof2 = size(subMat22,1)
    if (ndof2 > size(meqn2)) then
       ierr = internalError(prgnam_p//': Invalid matrix dimension')
       return
    end if

    select case (sysMat%storageType)

    case (diagonalMatrix_p)

       ierr = 0
       call DGATHR (ndof2,meqn2(1),sysMat%value(1),subMat22(1,1),1)

    case (denseMatrix_p)

       if (ndof2 > size(subMat22,2)) then
          ierr = internalError(prgnam_p//': Invalid matrix dimension')
          return
       end if

       ierr = 0
       do j = 1, ndof2
          call DGATHR (ndof2,meqn2(1),sysMat%value(1+sysMat%dim*(meqn2(j)-1)), &
               &       subMat22(1,j),1)
       end do

    case (skylineMatrix_p)

       call EXSEM (sysMat%value(1), sysMat%skyline%msky(1), &
            &      sysMat%dim, ndof2, getErrorFile(), subMat22(1,1), ierr)
       if (ierr < 0) goto 999

    case (sparseMatrix_p)

       if (ik > i4) then
          iWork => getIntegerScratchArray(sysMat%dim*2,ierr)
       else
          iWork => getIntegerScratchArray(sysMat%dim,ierr)
       end if
       if (ierr < 0) goto 999

       call SPRK22 (sysMat%sparse%mspar(1), sysMat%sparse%mtrees(1), &
            &       sysMat%sparse%msifa(1), sysMat%value(sysMat%dim+1), 2_ik, &
            &       iWork(1), subMat22(1,1), iWork(1))

    case (outOfCore_p)

       call ffa_startTimer ('FEExtractK22')
       call FEExtractK22 ('A', sysMat%gsf%Q, sysMat%gsf%S, subMat22, GSFinfo)
       call ffa_stopTimer ('FEExtractK22')
       ierr = checkGSFinfo(GSFinfo)
       if (ierr < 0) goto 999

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')

    end select

    return

999 continue
    call reportError (error_p,'Failed to extract sub-matrix from system '// &
         &            matrixType_p(sysMat%storageType),addString=prgnam_p)

  end subroutine csGetSub22


  !!============================================================================
  !> @brief Extracts characteristic measures for a system matrix.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] sysMat The system matrix to extract data from
  !> @param[out] dt Min, mean and max translational diagonal coefficients
  !> @param[out] dr Min, mean and max rotational diagonal coefficients
  !> @param[out] ierr Error flag
  !>
  !> @details The minimum, mean and maximum, translational and rotational,
  !> diagonal coefficents are extracted from the given system matrix.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Mar 2004

  subroutine csGetDiag (samData,sysMat,dt,dr,ierr)

    use kindModule         , only : dp, hugeVal_p, epsDiv0_p
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType, isAllocated, extractDiagonal
    use scratchArrayModule , only : getRealScratchArray
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(SamType)      , intent(in)    :: samData
    type(SysMatrixType), intent(inout) :: sysMat
    real(dp)           , intent(out)   :: dt(3), dr(3)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer           :: i, j, ieq, nTraD, nRotD
    real(dp), pointer :: diag(:)

    !! --- Logic section ---

    dt = 0.0_dp
    dr = 0.0_dp
    ierr = 0
    if (.not. isAllocated(sysMat)) return

    diag => getRealScratchArray(sysMat%dim,ierr)
    if (ierr < 0) goto 999

    call extractDiagonal (sysMat,diag,ierr)
    if (ierr < 0) goto 999

    nTraD = 0
    nRotD = 0
    dt(1) = hugeVal_p
    dr(1) = hugeVal_p

    do i = 1, samData%nnod
       do j = samData%madof(i), samData%madof(i+1)-1
          ieq = samData%meqn(j)
          if (ieq <= 0) cycle
          if (diag(ieq) < epsDiv0_p) cycle
          if (j-samData%madof(i) < 3) then ! translational dof
             dt(1) = min(dt(1),diag(ieq))
             dt(2) =     dt(2)+diag(ieq)
             dt(3) = max(dt(3),diag(ieq))
             nTraD = nTraD + 1
          else if (j-samData%madof(i) < 6) then ! rotational dof
             dr(1) = min(dr(1),diag(ieq))
             dr(2) =     dr(2)+diag(ieq)
             dr(3) = max(dr(3),diag(ieq))
             nRotD = nRotD + 1
          end if
       end do
    end do

    if (nTraD > 1) dt(2) = dt(2) / dble(nTraD)
    if (nRotD > 1) dr(2) = dr(2) / dble(nRotD)

    return

999 continue
    call reportError (debugFileOnly_p,'MatExtensionModule::csGetDiag')

  end subroutine csGetDiag

end module MatExtensionModule
