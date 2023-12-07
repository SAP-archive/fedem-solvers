!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module NormTypeModule

  use kindModule, only : dp

  implicit none

  integer, parameter :: oneOf_p = 1, &
       &                allOf_p = 2

  integer, parameter :: nNormTypes_p = 4
  integer, parameter :: iVecNorm_p = 1, &
       &                iInfTra_p  = 2, &
       &                iInfRot_p  = 3, &
       &                iInfGen_p  = 4


  type TestItemType
     character(len=12) :: title ! Name to be used in res-file headings
     integer  :: code      ! ON/OFF switch for this norm (0, oneOf_p, allOf_p)
     integer  :: nIterIncr ! Number of iterations with increasing norm value
     real(dp) :: value     ! The current value of this convergence norm
     real(dp) :: tolerance ! The actual convergence tolerance for this norm
     real(dp) :: absTol    ! Absolute convergence tolerance for this norm
     real(dp) :: relTol    ! Velocity proportional tolerance for this norm
  end type TestItemType

  type TestSetType
     logical            :: doingWell
     logical            :: isActive(3)
     integer            :: startMonitor
     integer , pointer  :: worstDOFs(:,:)
     real(dp), pointer  :: worstEnerg(:)
     type(TestItemType) :: disNorms(nNormTypes_p)
     type(TestItemType) :: velNorms(nNormTypes_p)
     type(TestItemType) :: accNorms(nNormTypes_p)
     type(TestItemType) :: resNorms(nNormTypes_p)
     type(TestItemType) :: energyNorms(2)
  end type TestSetType

  private :: InitTolerance


contains

  subroutine NullifyConvSet (convSet)

    !!==========================================================================
    !! Initialize the TestSetType object.
    !!
    !! Programmer : Knut Morten Okstad                           30 Jun 2004/1.0
    !!==========================================================================

    type(TestSetType), intent(out) :: convSet

    !! --- Logic section ---

    convSet%doingWell = .true.
    convSet%isActive = .false.
    convSet%startMonitor = 0
    nullify(convSet%worstDOFs)
    nullify(convSet%worstEnerg)

  end subroutine NullifyConvSet


  subroutine DeallocateConvSet (convSet)

    !!==========================================================================
    !! Deallocate the TestSetType object.
    !!
    !! Programmer : Knut Morten Okstad                           23 Jan 2017/1.0
    !!==========================================================================

    type(TestSetType), intent(inout) :: convSet

    !! --- Logic section ---

    if (associated(convSet%worstDOFs))  deallocate(convSet%worstDOFs)
    if (associated(convSet%worstEnerg)) deallocate(convSet%worstEnerg)

    call nullifyConvSet (convSet)

  end subroutine DeallocateConvSet


  subroutine InitConvChecks (convSet,maxIt,monWorst,monIter,relTol,err)

    !!==========================================================================
    !! Initialize the TestSetType object with tolerance data from command-line.
    !!
    !! Programmer : Knut Morten Okstad                           30 Jun 2004/1.0
    !!==========================================================================

    use reportErrorModule, only : allocationError, reportError, note_p,warning_p

    type(TestSetType), intent(inout) :: convSet
    integer          , intent(in)    :: maxIt, monWorst, monIter
    real(dp)         , intent(in)    :: relTol
    integer          , intent(out)   :: err

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call initTolerance ('tolDispNorm',' L2(displ)',convSet%disNorms(iVecNorm_p))
    call initTolerance ('tolDispTra','Inf(traDis)',convSet%disNorms(iInfTra_p))
    call initTolerance ('tolDispRot','Inf(angDis)',convSet%disNorms(iInfRot_p))
    call initTolerance ('tolDispGen','Inf(genDis)',convSet%disNorms(iInfGen_p))
    call initTolerance ('tolVelNorm','  L2(veloc)',convSet%velNorms(iVecNorm_p))
    call initTolerance ('tolVelTra' ,'Inf(traVel)',convSet%velNorms(iInfTra_p))
    call initTolerance ('tolVelRot' ,'Inf(angVel)',convSet%velNorms(iInfRot_p))
    call initTolerance ('tolVelGen' ,'Inf(genVel)',convSet%velNorms(iInfGen_p))
    call initTolerance ('tolAccNorm','  L2(accel)',convSet%accNorms(iVecNorm_p))
    call initTolerance ('tolAccTra' ,'Inf(traAcc)',convSet%accNorms(iInfTra_p))
    call initTolerance ('tolAccRot' ,'Inf(angAcc)',convSet%accNorms(iInfRot_p))
    call initTolerance ('tolAccGen' ,'Inf(genAcc)',convSet%accNorms(iInfGen_p))
    call initTolerance ('tolResNorm','    L2(res)',convSet%resNorms(iVecNorm_p))
    call initTolerance ('tolResTra' ,'Inf(traRes)',convSet%resNorms(iInfTra_p))
    call initTolerance ('tolResRot' ,'Inf(angRes)',convSet%resNorms(iInfRot_p))
    call initTolerance ('tolResGen' ,'Inf(genRes)',convSet%resNorms(iInfGen_p))
    call initTolerance ('tolEnerMax','E-norm(max)',convSet%energyNorms(1))
    call initTolerance ('tolEnerSum','E-norm(sum)',convSet%energyNorms(2))

    do i = 1, nNormTypes_p
       if (convSet%disNorms(i)%code > 0) convSet%isActive(1) = .true.
       if (convSet%velNorms(i)%code > 0) convSet%isActive(1) = .true.
       if (convSet%accNorms(i)%code > 0) convSet%isActive(1) = .true.
       if (convSet%resNorms(i)%code > 0) convSet%isActive(2) = .true.
    end do
    if (convSet%energyNorms(1)%code > 0) convSet%isActive(3) = .true.
    if (convSet%energyNorms(2)%code > 0) convSet%isActive(3) = .true.

    if (relTol > 0.0_dp .and. convSet%velNorms(iVecNorm_p)%code > 0) then
       convSet%velNorms(iVecNorm_p)%relTol = relTol
       call reportError (note_p,'Using velocity proportional tolerance '// &
            &            'for the scaled vector norm of velocity corrections')
    else if (.not. any(convSet%isActive)) then
       call reportError (warning_p,'No convergence tests defined')
    end if

    convSet%startMonitor = maxIt - monIter
    allocate(convSet%worstDOFs(monWorst,3), &
         &   convSet%worstEnerg(monWorst), stat=err)
    if (err /= 0) err = allocationError('InitConvChecks')

  end subroutine InitConvChecks


  subroutine InitTolerance (tolType,normName,testItem)

    !!==========================================================================
    !! Initialize the TestItemType object with tolerance data from command-line.
    !!
    !! Programmer : Bjorn Haugen                      date/rev : 28 Mar 2003/1.0
    !!              Knut Morten Okstad                           29 Apr 2003/2.0
    !!==========================================================================

    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble

    character(len=*)  , intent(in)    :: tolType, normName
    type(TestItemType), intent(inout) :: testItem

    !! Local variables
    real(dp) :: value

    !! --- Logic section ---

    testItem%title     = normName
    testItem%nIterIncr = 0
    testItem%code      = -1
    testItem%value     = 0.0_dp
    testItem%tolerance = 0.0_dp
    testItem%absTol    = 0.0_dp
    testItem%relTol    = 0.0_dp

    call ffa_cmdlinearg_getdouble (tolType,value)
    if (value > 0.0_dp) then
       testItem%absTol = value
       testItem%code   = allOf_p
    else if (value < 0.0_dp) then
       testItem%absTol = -value
       testItem%code   = oneOf_p
    end if

  end subroutine InitTolerance


  function HasConverged (convSet,factor_opt)

    !!==========================================================================
    !! Check if the state of the convSet is a converged set.
    !!
    !! Programmer : Bjorn Haugen                      date/rev : 28 Mar 2003/1.0
    !!==========================================================================

    type(TestSetType) , intent(in) :: convSet
    real(dp), optional, intent(in) :: factor_opt
    logical                        :: hasConverged

    !! Local variables
    logical :: testsExist(2), oneOfTest, allOfTest
    integer :: iNorm

    !! --- Logic section ---

    testsExist = .false.
    oneOfTest  = .false.
    allOfTest  = .true.

    do iNorm = 1, nNormTypes_p
       call TestItem (convSet%disNorms(iNorm))
       call TestItem (convSet%velNorms(iNorm))
       call TestItem (convSet%accNorms(iNorm))
       call TestItem (convSet%resNorms(iNorm))
    end do
    call TestItem (convSet%energyNorms(1))
    call TestItem (convSet%energyNorms(2))

    hasConverged = .true.
    if (testsExist(oneOf_p)) hasConverged = hasConverged .and. oneOfTest
    if (testsExist(allOf_p)) hasConverged = hasConverged .and. allOfTest

  contains

    subroutine TestItem (item)
      type(TestItemType), intent(in) :: item
      real(dp) :: tol
      if (present(factor_opt)) then
         tol = item%tolerance*factor_opt
      else
         tol = item%tolerance
      end if
      select case (item%code)
      case(oneOf_p); if (item%value <  tol) oneOfTest = .true.
      case(allOf_p); if (item%value >= tol) allOfTest = .false.
      case default ; return
      end select
      testsExist(item%code) = .true.
    end subroutine TestItem

  end function HasConverged


  subroutine checkDivergence (testItem,value,mayDiverge)

    !!==========================================================================
    !! Check if a given norm is showing sign of possible divergence.
    !! If the norm value is increasing in maxIterIncrBeforeWarning_p sequential
    !! iterations within a time step, a possible divergence is flagged.
    !!
    !! Programmer : Knut Morten Okstad                date/rev : 19 Mar 2004/1.0
    !!==========================================================================

    type(TestItemType), intent(inout) :: testItem
    real(dp)          , intent(in)    :: value
    logical           , intent(inout) :: mayDiverge

    !! Local variables
    integer, parameter :: maxIterIncrBeforeWarning_p = 2

    !! --- Logic section ---

    if (testItem%code <= 0) return ! Only check norms used in convergence checks

    if (value >= testItem%value .and. testItem%value >= testItem%tolerance) then
       testItem%nIterIncr = testItem%nIterIncr + 1
    else
       testItem%nIterIncr = 0
    end if

    if (testItem%nIterIncr >= maxIterIncrBeforeWarning_p) mayDiverge = .true.

  end subroutine checkDivergence

end module NormTypeModule
