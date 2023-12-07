!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file normRoutinesModule.f90
!>
!> @brief Subroutines and functions for solution norm calculation.

!!==============================================================================
!> @brief Module with subroutines and functions for solution norm calculation.

module NormRoutinesModule

  implicit none

  private :: ExpandVector, reportNaNerror


contains

  !!============================================================================
  !> @brief Defines the convergence tolerance to use for a given test item.
  !>
  !> @param item The test item to define convergence tolerance for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] u System vector to use with the relative tolerance @a relTol
  !> @param[in] w Weighting factors for the DOFs of different kind
  !>
  !> @details If a relative tolerance is specified (@a relTol > 0.0),
  !> the actual tolerance is taken as the scaled norm of the given vector @a u
  !> times this value, plus the absolute tolerance @a absTol.
  !> The vector @a u is typically the current velocity state.
  !>
  !> @author Bjorn Haugen                                      @date 28 Mar 2003
  !> @author Knut Morten Okstad                                @date 28 Jan 2004

  subroutine CalcTolerance (item,sam,u,w)

    use SamModule     , only : SamType, dp
    use NormTypeModule, only : TestItemType

    type(TestItemType)     , intent(inout) :: item
    type(SamType), optional, intent(in)    :: sam
    real(dp)     , optional, intent(in)    :: u(:), w(:)

    !! --- Logic section ---

    if (item%relTol > 0.0_dp .and. present(u) .and. present(w)) then

       !! u-proportional addition to the absolute tolerance
       item%tolerance = item%absTol + item%relTol * ScaledNorm(sam,u,w)

    else

       !! Absolute tolerance only
       item%tolerance = item%absTol

    end if

  end subroutine CalcTolerance


  !!============================================================================
  !> @brief Calculates the scaled norm of the given vector.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] vec System vector to calculate the scaled norm for
  !> @param[in] wFac Weighting factors for the DOFs of different kind
  !> @param[in] iPow The power of the L-norm (@a iPow = 2 implies L2-norm)
  !> @param[in] doExpand If .true., @a vec is assumed to be equation-oriented
  !> @param[out] ierr Error flag
  !> @return The norm value
  !>
  !> @callergraph @callgraph
  !>
  !> @author Knut Morten Okstad                                @date 27 Jun 2000
  !> @author Bjorn Haugen                                      @date 28 Mar 2003
  !> @author Knut Morten Okstad                                @date 28 Jan 2004

  function ScaledNorm (sam,vec,wFac,iPow,doExpand,ierr)

    use SamModule        , only : SamType, dp
    use ReportErrorModule, only : reportError, debugFileOnly_p

    type(SamType)    , intent(in)  :: sam
    real(dp)         , intent(in)  :: vec(:), wFac(:)
    integer, optional, intent(in)  :: iPow
    logical, optional, intent(in)  :: doExpand
    integer, optional, intent(out) :: ierr
    real(dp)                       :: ScaledNorm

    !! Local variables
    integer           :: i
    real(dp)          :: dPow, rn, sw, wf
    real(dp), pointer :: temp(:)

    !! --- Logic section ---

    ScaledNorm = -999.999_dp
    temp => ExpandVector(sam,vec,doExpand,ierr=ierr)
    if (present(ierr)) then
       if (ierr < 0) goto 990
    else if (.not.associated(temp)) then
       return
    end if

    if (present(iPow)) then
       dPow = real(max(1,iPow),dp)
    else
       dPow = 2.0_dp ! L2-norm is the default
    end if

    sw = 0.0_dp
    rn = 0.0_dp
    do i = 1, sam%ndof
       if (sam%meqn(i) > 0) then ! Has equation number (free DOF)
          if (isNaN(temp(i))) then
             if (present(ierr)) then
                ScaledNorm = temp(i) ! Norm is NaN
                ierr = ierr - 1
             else
                ScaledNorm = 0.0_dp ! Silently ignore NaN norms
                return
             end if
          end if
          wf = wFac(sam%dofType(i))
          sw = sw + wf**dPow
          rn = rn + abs(wf*temp(i))**dPow
       end if
    end do
    if (present(ierr)) then
       if (ierr < 0) goto 900
    end if

    if (sw > 0.0_dp) then
       ScaledNorm = (rn/sw)**(1.0_dp/dPow)
    else
       ScaledNorm = 0.0_dp
    end if

    return

900 call reportNaNerror (-ierr)
990 call reportError (debugFileOnly_p,'NormRoutinesModule::ScaledNorm')

  end function ScaledNorm


  !!============================================================================
  !> @brief Calculates infinity norms of the given vector.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] vec System vector to calculate infinity norms for
  !> @param[in] findWorst If .true., record the worst DOFs in array @a worsDOFs
  !> @param[out] infTra The maximum translation DOF value
  !> @param[out] infRot The maximum rotation DOF value
  !> @param[out] infGen The maximum generalized DOF value
  !> @param[out] worstDOFs Array to record the worst DOFs (highest values) in
  !> @param[in] doExpand If .true., @a vec is assumed to be equation-oriented
  !> @param[out] ierr Error flag
  !>
  !> @callergraph @callgraph
  !>
  !> @author Bjorn Haugen                                      @date 28 Mar 2003
  !> @author Knut Morten Okstad                                @date 28 Jan 2004

  subroutine InfiniteNorms (sam,vec,findWorst,infTra,infRot,infGen, &
       &                    worstDOFs,doExpand,ierr)

    use SamModule          , only : SamType, dp
    use SearchAndSortModule, only : findMaxRealValues
    use ReportErrorModule  , only : reportError, debugFileOnly_p

    type(SamType)    , intent(in)  :: sam
    real(dp)         , intent(in)  :: vec(:)
    logical          , intent(in)  :: findWorst
    real(dp)         , intent(out) :: infTra, infRot, infGen
    integer, optional, intent(out) :: worstDOFs(:)
    logical, optional, intent(in)  :: doExpand
    integer          , intent(out) :: ierr

    !! Local variables
    integer           :: i, iDof, nWorst
    real(dp)          :: infNorm(3), aVal
    real(dp), pointer :: temp(:)

    !! --- Logic section ---

    if (findWorst .and. present(worstDOFs)) then
       nWorst = size(worstDOFs)
    else
       nWorst = 0
    end if

    temp => ExpandVector(sam,vec,doExpand,nWorst>0,ierr)
    if (ierr < 0 .or. .not.associated(temp)) goto 990

    infNorm = 0.0_dp
    do i = 1, sam%ndof
       if (sam%meqn(i) > 0) then ! Has equation number (free DOF)
          iDof = sam%dofType(i)
          if (isNaN(temp(i))) then
             infNorm(iDof) = temp(i) ! Norm is NaN
             ierr = ierr - 1
          else
             aVal = abs(temp(i))
             if (.not.isNaN(infNorm(iDof)) .and. aVal > infNorm(iDof)) then
                infNorm(iDof) = aVal
             end if
             if (nWorst > 0) temp(i) = aVal
          end if
       else if (nWorst > 0) then
          temp(i) = 0.0_dp
       end if
    end do
    if (ierr < 0) goto 900

    infTra = infNorm(1) ! Translatory DOFs
    infRot = infNorm(2) ! Rotational DOFs
    infGen = infNorm(3) ! Generalized DOFs

    if (nWorst > 0) then
       call findMaxRealValues (temp,worstDOFs,ierr=ierr)
       if (ierr < 0) goto 990
    end if

    return

900 call reportNaNerror (-ierr)
990 call reportError (debugFileOnly_p,'NormRoutinesModule:InfiniteNorms')

  end subroutine InfiniteNorms


  !!============================================================================
  !> @brief Calculates energy norms for given force and displacement vectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] disp System displacement vector to calculate energy norms for
  !> @param[in] force System force vector to calculate energy norms for
  !> @param[in] findWorst If .true., record the worst DOFs in array @a worsDOFs
  !> @param[out] E_absSum The mean energy density among all DOFs
  !> @param[out] E_absMax The maximum energy density among all DOFs
  !> @param[out] worstDOFs Array to record the worst DOFs (highest values) in
  !> @param[out] worstVals Array to record the worst DOF values in
  !> @param[out] ierr Error flag
  !>
  !> @callergraph @callgraph
  !>
  !> @author Bjorn Haugen                                      @date 28 Mar 2003
  !> @author Knut Morten Okstad                                @date 28 Jan 2004

  subroutine EnergyNorms (sam,disp,force,findWorst, &
       &                  E_absSum,E_absMax,worstDOFs,worstVals,ierr)

    use SamModule          , only : SamType, dp
    use SearchAndSortModule, only : findMaxRealValues
    use ReportErrorModule  , only : internalError, reportError, debugFileOnly_p

    type(SamType), intent(in)  :: sam
    real(dp)     , intent(in)  :: disp(:), force(:)
    logical      , intent(in)  :: findWorst
    real(dp)     , intent(out) :: E_absSum, E_absMax
    integer      , intent(out) :: worstDOFs(:)
    real(dp)     , intent(out) :: worstVals(:)
    integer      , intent(out) :: ierr

    !! Local variables
    logical, parameter :: doExpand_p = .true.
    integer            :: i, nDof, nWorst
    real(dp)           :: E_abs
    real(dp), pointer  :: temp(:)

    !! --- Logic section ---

    if (findWorst) then
       nWorst = size(worstDOFs)
    else
       nWorst = 0
    end if

    !! Assume the incomming force vector is in equation order
    temp => ExpandVector(sam,force,doExpand_p,ierr=ierr)
    if (ierr < 0 .or. .not.associated(temp)) goto 990

    !! Assume the incomming displacement vector already is in dof order
    nDof = size(disp)
    ierr = nDof - size(temp)
    if (ierr /= 0) then
       ierr = internalError('NormRoutinesModule::EnergyNorms')
       return
    end if

    E_absSum = 0.0_dp
    E_absMax = 0.0_dp
    do i = 1, nDof
       if (sam%meqn(i) > 0) then ! Has equation number (free DOF)
          !! Use abs here so that no accidental zero sum
          E_abs = abs(disp(i)*temp(i))
          if (isNaN(E_abs)) then ! Norm is NaN
             E_absMax = E_abs
             E_absSum = E_abs
             ierr = ierr - 1
          else if (.not. isNaN(E_absMax)) then
             E_absSum = E_absSum + E_abs
             if (E_abs > E_absMax) E_absMax = E_abs
             if (nWorst > 0) temp(i) = E_abs
          end if
       else if (nWorst > 0) then
          temp(i) = 0.0_dp
       end if
    end do
    if (ierr < 0) goto 900

    !! Use mean so that it is more physical equivalent to infinite norm
    if (nDof > 1) E_absSum = E_absSum / dble(sam%neq)

    if (nWorst > 0) then
       !! Find the worst energy DOFs and store their associated values
       call findMaxRealValues (temp,worstDOFs,worstVals,ierr)
       if (ierr < 0) goto 990
    end if

    return

900 call reportNaNerror (-ierr)
990 call reportError (debugFileOnly_p,'NormRoutinesModule::EnergyNorms')

  end subroutine EnergyNorms


  !!============================================================================
  !> @brief Returns a pointer to the expanded solution vector.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] vec System vector to expand from equation-order to DOF-order
  !> @param[in] doExpand If not present or .false., don't do any expansion
  !> @param[in] forceCopy If .true., copy the vector if it is not expanded
  !> @param[out] ierr Error flag
  !> @return Pointer to expanded or copied vector
  !>
  !> @details The expanded/copied vector is actually stored in the internal
  !> scratch array managed by the scratcharraymodule, so use with utmost care.
  !>
  !> @sa solextensionmodule::csexpand
  !>
  !> @callergraph @callgraph
  !>
  !> @author Bjorn Haugen                                      @date 13 Mar 2003
  !> @author Knut Morten Okstad                                @date 28 Jan 2004

  function ExpandVector (sam,vec,doExpand,forceCopy,ierr) result(temp)

    use SamModule         , only : SamType, dp
    use SolExtensionModule, only : csExpand
    use ScratchArrayModule, only : getRealScratchArray
    use ReportErrorModule , only : internalError, reportError, debugFileOnly_p

    type(SamType)    , intent(in)  :: sam
    real(dp), target , intent(in)  :: vec(:)
    logical, optional, intent(in)  :: doExpand, forceCopy
    integer, optional, intent(out) :: ierr
    real(dp)         , pointer     :: temp(:)

    !! Local variables
    logical :: do_Expand, force_Copy
    integer :: lerr

    !! --- Logic section ---

    if (present(doExpand)) then
       do_Expand = doExpand
    else
       do_Expand = .false.
    end if
    if (present(forceCopy)) then
       force_Copy = forceCopy
    else
       force_Copy = .false.
    end if

    if (do_Expand) then
       lerr = size(vec) - max(1,sam%neq)
    else
       lerr = size(vec) - sam%ndof
    end if
    if (lerr /= 0) then
       nullify(temp)
       lerr = internalError('NormRoutinesModule::ExpandVector')
    else if (do_Expand .or. force_Copy) then
       temp => getRealScratchArray(sam%ndof,lerr)
       if (lerr /= 0) then
          nullify(temp)
          call reportError (debugFileOnly_p,'NormRoutinesModule::ExpandVector')
       else if (do_Expand) then
          call csExpand (sam,vec,temp)
       else
          temp = vec ! Vector is not expanded, but force a copy of it
       end if
    else
       temp => vec ! No expansion needed, return a pointer to the input vector
    end if
    if (present(ierr)) ierr = lerr

  end function ExpandVector


  !> @cond FULL_DOC
  !!============================================================================
  !> @brief Gives an error message on NaN's in the solution vector.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Nov 2008

  subroutine reportNaNerror (nNaN)

    use ReportErrorModule, only : reportError, error_p

    integer, intent(in) :: nNaN !< Number of not-a-number values detected

    !! Local variables
    character(len=64) :: errmsg

    !! --- Logic section ---

    write(errMsg,"(I8,' NaN values detected in the solution vector.')") nNaN
    call reportError (error_p,adjustl(errmsg), &
         &            'The model is highly divergent. Check model consistency!')

  end subroutine reportNaNerror
  !> @endcond

end module NormRoutinesModule
