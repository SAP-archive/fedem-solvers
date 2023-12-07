!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file vpmSolver/environmentTypeModule.f90
!> @brief Environmental data container.

!!==============================================================================
!> @brief Module with environment data containers.
!>
!> @details This module contains some global parameters and data not associated
!> with the mechanism itself, but the environment it resides in. The data are
!> the used to calculate the external loads on the mechanism that depends on
!> the environment, such as wind, waves, current, etc.

module EnvironmentTypeModule

  use FunctionTypeModule, only : FunctionType, EngineType, dp

  implicit none

  !> @brief Data type containing some global environmental data.
  type EnvironmentType

     real(dp)                    :: rhoA       !< Air mass density
     real(dp)                    :: rhoW       !< Water mass densities
     real(dp)                    :: rhoInt     !< Internal fluid mass density
     real(dp)                    :: rhoG       !< Marine growth mass density
     real(dp)                    :: tG         !< Marine growth thickness
     real(dp)                    :: zG(2)      !< Marine growth depth range
     real(dp)                    :: grav0(3)   !< Unscaled gravitation vector
     real(dp)                    :: gravity(3) !< Scaled gravitation vector
     real(dp)                    :: sea0       !< Mean sea water level, constant
     real(dp)                    :: seaLevel   !< Sea level at current time step
     real(dp)                    :: seaScale   !< Sea level scaling factor
     real(dp)                    :: seaDepth   !< Sea depth
     real(dp)                    :: currScale  !< Sea current scaling factor
     real(dp)                    :: hdfScale   !< Hydrodynamic force scaling
     real(dp)                    :: Tsea(3,4)  !< Sea coordinate system
     integer                     :: waveTheory !< Which wave theory model to use
     type(FunctionType), pointer :: waveFunc   !< Sea level function, f(t)
     type(EngineType)  , pointer :: wScaleEng  !< Sea level scaling engine
     type(FunctionType), pointer :: currFunc   !< Sea current function, f(z)
     type(FunctionType), pointer :: cDirFunc   !< Sea current direction, f(z)
     type(EngineType)  , pointer :: cScaleEng  !< Sea current scaling engine
     type(EngineType)  , pointer :: hScaleEng  !< Hydrodynamic force scaling

  end type EnvironmentType

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteEnvironmentType
  end interface

  private :: WriteEnvironmentType


contains

  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] env The environmenttypemodule::environmenttype object to write
  !> @param[in] io File unit number to write to
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Oct 2009

  subroutine WriteEnvironmentType (env,io)

    use IdTypeModule, only : getId

    type(EnvironmentType), intent(in) :: env
    integer              , intent(in) :: io

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'Environment','{'
    write(io,*) 'rho(air)    =', env%rhoa
    write(io,*) 'rho(water)  =', env%rhow
    write(io,*) 'rho(int)    =', env%rhoint
    write(io,*) 'rho(growth) =', env%rhog
    write(io,*) 'tg, zg      =', env%tg, env%zg
    write(io,*) 'gravity     =', env%gravity
    write(io,*) 'sea depth   =', env%seaDepth
    write(io,*) 'sea0(level) =', env%sea0
    write(io,*) 'sea(level)  =', env%seaLevel
    write(io,*) 'sea(scale)  =', env%seaScale
    write(io,*) 'curr(scale) =', env%currScale
    write(io,*) 'hdf(scale)  =', env%hdfScale
    write(io,600) (env%Tsea(i,:),i=1,3)
    write(io,*) 'waveTheory  =', env%waveTheory
    if (associated(env%waveFunc)) then
       write(io,*) 'wave(func)  =', trim(getId(env%waveFunc%id))
    end if
    if (associated(env%wScaleEng)) then
       write(io,*) 'wscale(eng) =', trim(getId(env%wScaleEng%id))
    end if
    if (associated(env%currFunc)) then
       write(io,*) 'curr(func)  =', trim(getId(env%currFunc%id))
    end if
    if (associated(env%cDirFunc)) then
       write(io,*) 'cudir(func) =', trim(getId(env%cDirFunc%id))
    end if
    if (associated(env%cScaleEng)) then
       write(io,*) 'cscale(eng) =', trim(getId(env%cScaleEng%id))
    end if
    if (associated(env%hScaleEng)) then
       write(io,*) 'hscale(eng) =', trim(getId(env%hScaleEng%id))
    end if
    write(io,'(A/)') '}'

600 format(' T(sea)      =',4ES15.6E3/(14X,4ES15.6E3))

  end subroutine WriteEnvironmentType


  !!============================================================================
  !> @brief Initializes the EnvironmentType object.
  !>
  !> @param env The environmenttypemodule::environmenttype object to initialize
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Oct 2009

  subroutine NullifyEnvironment (env)

    type(EnvironmentType), intent(out) :: env

    !! --- Logic section ---

    env%rhoa = 0.0_dp
    env%rhow = 0.0_dp
    env%rhoInt = 0.0_dp
    env%rhoG = 0.0_dp
    env%tG = 0.0_dp
    env%zG = 0.0_dp
    env%grav0 = 0.0_dp
    env%gravity = 0.0_dp
    env%seaDepth = 0.0_dp
    env%sea0 = 0.0_dp
    env%seaLevel = 0.0_dp
    env%seaScale = 1.0_dp
    env%currScale = 1.0_dp
    env%hdfScale = 1.0_dp
    env%Tsea = 0.0_dp
    env%waveTheory = 0

    nullify(env%waveFunc)
    nullify(env%wScaleEng)
    nullify(env%currFunc)
    nullify(env%cDirFunc)
    nullify(env%cScaleEng)
    nullify(env%hScaleEng)

  end subroutine NullifyEnvironment


  !!============================================================================
  !> @brief Initializes the environmental data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] env The environmental data of the model
  !> @param[out] err Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2009

  subroutine ReadEnvironment (infp,env,err)

    use FunctionTypeModule, only : GetPtrToId, nullifyFunction, nullifyEngine
    use IdTypeModule      , only : ReportInputError
    use inputUtilities    , only : iuSetPosAtNextEntry
    use reportErrorModule , only : allocationError
    use progressModule    , only : lterm
    use EnvironmentNamelistModule ! Defines all variables in /ENVIRONMENT/
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer              , intent(in)  :: infp
    type(EnvironmentType), intent(out) :: env
    integer              , intent(out) :: err

    !! --- Logic section ---

    !! Read environment data (if any) from the solver input file
    if (.not. iuSetPosAtNextEntry(infp,'&ENVIRONMENT')) then
       rewind(infp)
       err = 0
       return
    end if

    call read_ENVIRONMENT (infp,err)
    if (err /= 0) then
       call ReportInputError ('ENVIRONMENT')
       return
    end if

    write(lterm,*) 'Read &ENVIRONMENT'
    if (rhoAir      > 0.0_dp) write(lterm,6) 'rhoAir      =',rhoAir
    if (rhoWater    > 0.0_dp) write(lterm,6) 'rhoWater    =',rhoWater
    if (rhoInternal > 0.0_dp) write(lterm,6) 'rhoInternal =',rhoInternal
    if (rhoGrowth   > 0.0_dp) write(lterm,6) 'rhoGrowth   =',rhoGrowth
    if (tGrowth     > 0.0_dp) write(lterm,6) 'tGrowth     =',tGrowth
    if (any(abs(zGrowth) > 0.0_dp)) write(lterm,6) 'zGrowth     =',zGrowth
    if (any(abs(gravity) > 0.0_dp)) write(lterm,6) 'gravity     =',gravity
    if (seaDepth    > 0.0_dp) write(lterm,6) 'seaDepth    =',seaDepth
    if (any(abs(seaCS) > 0.0_dp)) then
       write(lterm,6) 'seaCS       =',seaCS(:,1)
       write(lterm,6) '             ',seaCS(:,2)
       write(lterm,6) '             ',seaCS(:,3)
    end if
  6 format(4X,A,1P4E13.5)

    env%rhoa = rhoAir
    env%rhow = rhoWater
    env%rhoInt = rhoInternal
    env%rhoG = rhoGrowth
    env%tG = tGrowth
    env%zG = zGrowth
    env%grav0 = gravity
    if ( ffa_cmdlinearg_intValue('rampSteps') > 0 .and. &
         ffa_cmdlinearg_isTrue('rampGravity') ) then
       env%gravity = 0.0_dp ! Gravity forces are ramped and will start at zero
    else
       env%gravity = gravity
    end if
    env%seaDepth = seaDepth
    env%sea0 = sea0
    env%Tsea = transpose(seaCS)

    !! Allocate temporary function objects, for storage of their baseID
    !! until the actual function objects are read later
    if (err == 0 .and. waveFunction > 0) then
       allocate(env%waveFunc,STAT=err)
       if (err == 0 .and. waveScaleEngine > 0) allocate(env%wScaleEng,STAT=err)
       if (err == 0 .and. hdfScaleEngine > 0)  allocate(env%hScaleEng,STAT=err)
    end if
    if (err == 0 .and. currFunction > 0) then
       allocate(env%currFunc,STAT=err)
       if (err == 0 .and. currDirFunction > 0) allocate(env%cDirFunc ,STAT=err)
       if (err == 0 .and. currScaleEngine > 0) allocate(env%cScaleEng,STAT=err)
    end if

    if (associated(env%waveFunc)) then
       call nullifyFunction (env%waveFunc)
       env%waveFunc%id%baseId = waveFunction
    end if
    if (associated(env%wScaleEng)) then
       call nullifyEngine (env%wScaleEng)
       env%wScaleEng%id%baseId = waveScaleEngine
    end if
    if (associated(env%currFunc)) then
       call nullifyFunction (env%currFunc)
       env%currFunc%id%baseId = currFunction
    end if
    if (associated(env%cDirFunc)) then
       call nullifyFunction (env%cDirFunc)
       env%cDirFunc%id%baseId = currDirFunction
    end if
    if (associated(env%cScaleEng)) then
       call nullifyEngine (env%cScaleEng)
       env%cScaleEng%id%baseId = currScaleEngine
    end if
    if (associated(env%hScaleEng)) then
       call nullifyEngine (env%hScaleEng)
       env%hScaleEng%id%baseId = hdfScaleEngine
    end if

    if (err /= 0) err = allocationError('ReadEnvironment')

  end subroutine ReadEnvironment


  !!============================================================================
  !> @brief Initializes environmental settings after reading the input file.
  !>
  !> @param env The environmental data of the model
  !> @param[in] functions All function shapes in the model
  !> @param[in] engines All general functions in the model
  !> @param[out] err Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2009

  subroutine InitiateEnvironment (env,functions,engines,err)

    use FunctionTypeModule, only : GetPtrToId
    use IdTypeModule      , only : getId
    use reportErrorModule , only : reportError, note_p, error_p

    type(EnvironmentType)     , intent(out) :: env
    type(FunctionType), target, intent(in)  :: functions(:)
    type(EngineType)  , target, intent(in)  :: engines(:)
    integer                   , intent(out) :: err

    !! Local variables
    integer            :: baseId
    character(len=128) :: MSG

    !! --- Logic section ---

    err = 0

    if (associated(env%waveFunc)) then
       baseId = env%waveFunc%id%baseId
       deallocate(env%waveFunc)
       env%waveFunc => GetPtrToId(functions,baseId)
       if (associated(env%waveFunc)) then
          call reportError (note_p,'Sea water level is governed by Function'// &
               &            getId(env%waveFunc%id))
          if (env%seaDepth > 0.0_dp) then
             write(MSG,610) env%seaDepth
             call reportError (note_p,MSG)
610          format('Using finite water depth (',1PE12.5,' [m])')
          end if
       else
          err = err - 1
       end if
    else
       env%seaLevel = env%sea0 ! Constant water level
    end if

    if (associated(env%wScaleEng)) then
       baseId = env%wScaleEng%id%baseId
       deallocate(env%wScaleEng)
       env%wScaleEng => GetPtrToId(engines,baseId)
       if (associated(env%wScaleEng)) then
          call reportError (note_p,'Sea water level is scaled by Engine'// &
               &            getId(env%wScaleEng%id))
       else
          err = err - 1
       end if
    end if

    if (associated(env%currFunc)) then
       baseId = env%currFunc%id%baseId
       deallocate(env%currFunc)
       env%currFunc => GetPtrToId(functions,baseId)
       if (associated(env%currFunc)) then
          call reportError (note_p,'Sea current is governed by Function'// &
               &            getId(env%currFunc%id))
       else
          err = err - 1
       end if
    end if

    if (associated(env%cDirFunc)) then
       baseId = env%cDirFunc%id%baseId
       deallocate(env%cDirFunc)
       env%cDirFunc => GetPtrToId(functions,baseId)
       if (associated(env%cDirFunc)) then
          call reportError (note_p,'Sea current direction is governed '// &
               &            'by Function'//getId(env%cDirFunc%id))
       else
          err = err - 1
       end if
    end if

    if (associated(env%cScaleEng)) then
       baseId = env%cScaleEng%id%baseId
       deallocate(env%cScaleEng)
       env%cScaleEng => GetPtrToId(engines,baseId)
       if (associated(env%cScaleEng)) then
          call reportError (note_p,'Sea current is scaled by Engine'// &
               &            getId(env%cScaleEng%id))
       else
          err = err - 1
       end if
    end if

    if (associated(env%hScaleEng)) then
       baseId = env%hScaleEng%id%baseId
       deallocate(env%hScaleEng)
       env%hScaleEng => GetPtrToId(engines,baseId)
       if (associated(env%hScaleEng)) then
          call reportError (note_p,'Hydrodynamic forces are scaled by Engine'//&
               &            getId(env%hScaleEng%id))
       else
          err = err - 1
       end if
    end if

    if (err < 0) then
       call reportError (error_p,'Invalid sea function/engine specification', &
            &            addString='InitiateEnvironment')
    end if

  end subroutine InitiateEnvironment


  !!============================================================================
  !> @brief Transforms given global coordinates to the sea coordinate system.
  !>
  !> @param[in] env The environmental data of the model
  !> @param[in] xg Spatial coordinates of a point in global coordinate system
  !> @return Spatial coordinates of a point in the sea coordinate system
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Feb 2013

  function glob2Sea (env,xg) result(xs)

    use manipMatrixModule, only : matmul34, invert34

    type(EnvironmentType), intent(in) :: env
    real(dp)             , intent(in) :: xg(3)

    !! Local variables
    real(dp) :: xs(3)

    !! --- Logic section ---

    xs = matmul34(invert34(env%Tsea),xg)

  end function glob2Sea


  !!============================================================================
  !> @brief Updates the scaled gravitation vector in case of dynamic up-ramping.
  !>
  !> @param env The environmental data of the model
  !> @param[in] time Current simulation time
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Aug 2022

  subroutine updateGravity (env,time,ierr)

    use FunctionTypeModule, only : getRampValue
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(in)    :: time
    integer              , intent(out)   :: ierr

    !! Local variables
    real(dp) :: tmpT, rampVal

    !! --- Logic section ---

    ierr = 0
    tmpT = time
    call getRampValue (tmpT,0,rampVal,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'updateGravity')
    else
       env%gravity = rampVal * env%grav0
    end if

  end subroutine updateGravity

end module EnvironmentTypeModule
