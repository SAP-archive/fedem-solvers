!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateModesTypeModule.f90
!> @brief Subroutine for initialization of the eigenvalue solver.

!!==============================================================================
!> @brief Module with a subroutine for initialization of the eigenvalue solver.

module initiateModesTypeModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes the eigenvalue solver.
  !>
  !> @param sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param modes Data for eigenmodes
  !> @param[out] err Error flag
  !> @param[in] nModes Number of eigen modes to calculate
  !> @param[in] eigSolver Which eigensolver to use (override default)
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 31 Aug 2000

  subroutine initiateModes (sam,sys,modes,err,nModes,eigSolver)

    use KindModule            , only : dp
    use SamModule             , only : SamType
    use SysMatrixTypeModule   , only : skylineMatrix_p, denseMatrix_p
    use SysMatrixTypeModule   , only : shareMatrixStructure, allocateSysMatrix
    use SystemTypeModule      , only : SystemType
    use ModesTypeModule       , only : ModesType, nullifyModes
    use reportErrorModule     , only : allocationError
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdoubles
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    type(SamType)   , intent(inout) :: sam
    type(SystemType), intent(inout) :: sys
    type(ModesType) , intent(inout) :: modes
    integer         , intent(out)   :: err
    integer,optional, intent(in)    :: nModes, eigSolver

    !! Local variables
    logical :: doInit

    !! --- Logic section ---

    err = 0

    !! Check if eigenvalue analysis is requested
    if (present(nModes)) then
       if (present(eigSolver)) modes%solver = eigSolver
       !! Already initialized, check that arrays are large enough
       doInit = modes%nModes < 1
       if (nModes > modes%nModes) then
          !! Need to reallocate
          if (associated(modes%effMass)) deallocate(modes%effMass)
          if (associated(modes%ReVal))   deallocate(modes%ReVal)
          if (associated(modes%ReVec))   deallocate(modes%ReVec)
          if (associated(modes%dampRat)) deallocate(modes%dampRat)
          if (associated(modes%ImVal))   deallocate(modes%ImVal)
          if (associated(modes%ImVec))   deallocate(modes%ImVec)
          modes%nModes = nModes
       else
          modes%nModes = nModes
          return
       end if
    else
       doInit = .true.
       call nullifyModes (modes)
       call ffa_cmdlinearg_getint ('numEigModes',modes%nModes)
       if (modes%nModes < 1) return
    end if

    if (doInit) then

       call ffa_cmdlinearg_getbool ('stressStiffEig',modes%stressStiffIsOn)

       if (.not.present(eigSolver)) then
          if (sys%Nmat%storageType == skylineMatrix_p) then
             if (ffa_cmdlinearg_isTrue('lancz1')) then
                modes%solver = 1 ! Using LANCZ1 (workaround for smaller systems)
             else
                modes%solver = 2 ! Using LANCZ2 is the default
             end if
          else if (sys%Nmat%storageType == denseMatrix_p) then
             modes%solver = 3 ! Using DGGEVX eigensolver for dense systems
          end if
          if (ffa_cmdlinearg_isTrue('damped')) then
             modes%solver = 4 ! Using DGGEVX eigensolver for damped systems
          else if (ffa_cmdlinearg_isTrue('symmetric')) then
             modes%solver = 5 ! Using DSYGVX eigensolver for symmetric systems
          else if (ffa_cmdlinearg_isTrue('undamped')) then
             modes%solver = 3 ! Using DGGEVX eigensolver for undamped systems
          end if
       end if
       if (sys%Nmat%storageType == skylineMatrix_p) then
          sys%Nmat%skyline%eigenSolver = modes%solver
       end if

       !! Share system matrix datastructures
       call shareMatrixStructure (modes%Kmat,sys%Nmat)
       call shareMatrixStructure (modes%Cmat,sys%Nmat)
       call shareMatrixStructure (modes%Mmat,sys%Nmat)

       !! Allocate system matrices for eigenvalue solver
       call allocateSysMatrix (modes%Kmat,err)
       call allocateSysMatrix (modes%Cmat,err)
       call allocateSysMatrix (modes%Mmat,err)
       if (err < 0) then
          call reportError (debugFileOnly_p,'InitiateModes')
          return
       end if

       call ffa_cmdlinearg_getdoubles ('tolFactorize',modes%tol,2)
       call ffa_cmdlinearg_getdouble ('tolEigval',modes%tol(1))
       call ffa_cmdlinearg_getdouble ('tolEigvector',modes%tol(3))
       call ffa_cmdlinearg_getdouble ('eigenshift',modes%shift)
       call ffa_cmdlinearg_getbool ('factorMass_eigensolver',modes%factorMass)
       call ffa_cmdlinearg_getbool ('addBC_eigensolver',modes%addBC)

    end if

    modes%maxLan = 3*modes%nModes + 12 ! Safe estimate (used by Lanczos only)
    if (modes%addBC) then
       if (modes%nModes > sam%ndof1) modes%nModes = sam%ndof1
       if (modes%maxLan > sam%ndof1) modes%maxLan = sam%ndof1
    else
       if (modes%nModes > sam%neq)   modes%nModes = sam%neq
       if (modes%maxLan > sam%neq)   modes%maxLan = sam%neq
    end if

    if (ffa_cmdlinearg_isTrue('effModalMass')) then
       allocate(modes%effMass(modes%nModes,6), STAT=err)
       if (err /= 0) then
          err = AllocationError('InitiateModes 1')
          return
       end if
       modes%effMass = 0.0_dp
    end if

    !! Allocate real eigenvalues and eigenvectors
    allocate(modes%ReVal(modes%nModes), &
         &   modes%ReVec(sam%ndof,modes%nModes), &
         &   modes%dampRat(modes%nModes), STAT=err)
    if (err /= 0) then
       err = AllocationError('InitiateModes 2')
       return
    end if

    modes%ReVal = 0.0_dp
    modes%ReVec = 0.0_dp
    modes%dampRat = 0.0_dp
    if (modes%solver /= 4) return

    !! Allocate imaginary eigenvalues and eigenvectors
    allocate(modes%ImVal(modes%nModes), &
         &   modes%ImVec(sam%ndof,modes%nModes), STAT=err)
    if (err /= 0) then
       err = AllocationError('InitiateModes 3')
       return
    end if

    modes%ImVal = 0.0_dp
    modes%ImVec = 0.0_dp

  end subroutine initiateModes

end module initiateModesTypeModule
