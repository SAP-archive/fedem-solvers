!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateWindTurbineModule.f90
!> @brief Initialization of wind turbine objects from the solver input file.

!!==============================================================================
!> @brief Initialization of wind turbine objects from the solver input file.

module initiateWindTurbineModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes the wind turbine with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param env Environmental data
  !> @param turbine Wind turbine configuration object
  !> @param triads Array of all triads in the model
  !> @param[in] sups Array of all superelements in the model
  !> @param[in] joints Array of all joints in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2009

  subroutine ReadTurbineConfig (infp,env,turbine,triads,sups,joints,err)

    use EnvironmentTypeModule     , only : EnvironmentType
    use WindTurbineTypeModule     , only : TurbineConfig, nullifyTurbine
    use TriadTypeModule           , only : TriadType, TriadPtrType, GetPtrToId
    use SupElTypeModule           , only : SupElType, GetPtrToId, IsBeam
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, GetPtrToId
    use MasterSlaveJointTypeModule, only : REVOLUTE_p
    use kindModule                , only : lfnam_p, dp
    use IdTypeModule              , only : initId, getId, StrId
    use IdTypeModule              , only : ldesc_p, ReportInputError
    use WindTurbineRoutinesModule , only : AeroInput
    use progressModule            , only : lterm
    use inputUtilities            , only : iuSetPosAtNextEntry
    use inputUtilities            , only : iuGetNumberOfEntries
    use reportErrorModule         , only : allocationError, internalError
    use reportErrorModule         , only : reportError, getErrorFile
    use reportErrorModule         , only : debugFileOnly_p
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_isTrue

    integer                   , intent(in)    :: infp
    type(EnvironmentType)     , intent(inout) :: env
    type(TurbineConfig)       , pointer       :: turbine
    type(TriadType), target   , intent(inout) :: triads(:)
    type(SupElType)           , intent(in)    :: sups(:)
    type(MasterSlaveJointType), intent(in)    :: joints(:)
    integer                   , intent(out)   :: err

    !! Local variables
    integer, parameter   :: maxTB_p = 100
    integer              :: i, j, n, lpu, stat, inod, nBelm
    integer, allocatable :: MENC(:,:), MNEC(:,:)
    real(dp)             :: R, DR, a(3), b(3), v(3)

    type(TriadType)           , pointer       :: triad
    type(TriadPtrType)        , allocatable   :: bldTriads(:,:)
    type(MasterSlaveJointType), pointer, save :: generator, joint

    !! Define the TURBINE_CONFIG namelist
    integer  :: id, extId(10)
    character(ldesc_p) :: extDescr
    character(lfnam_p) :: ADFile
    real(dp) :: PtfmRef, HubRad, TipRad, acOffSet(maxTB_p), ADcentre(2,maxTB_p)
    integer  :: hubId, generatorJoint, pitchJoint(10)
    integer  :: towerTriad, nacelleTriad, shaftTriad, azimuthTriad, hubTriad
    integer  :: nBlade, nTB, firstTriadId(10), ADnodes(maxTB_p)
    logical  :: CompAero, CompNoise, SumPrint, UserID

    namelist /TURBINE_CONFIG/ id, extId, extDescr, ADFile, &
         &                    PtfmRef, HubRad, TipRad, acOffSet, ADcentre, &
         &                    hubId, generatorJoint, pitchJoint, nBlade, nTB, &
         &                    firstTriadId, towerTriad, nacelleTriad, &
         &                    shaftTriad, azimuthTriad, hubTriad, ADnodes, &
         &                    CompAero, CompNoise, SumPrint, UserID

    !! Define the WIND_POINT namelist
    real(dp) :: coord(3)

    namelist /WIND_POINT/ coord


    !! --- Logic section ---

    err = 0
    if (.not. iuSetPosAtNextEntry(infp,'&TURBINE_CONFIG')) return
    if (ffa_cmdlinearg_isTrue('ignoreAD')) return ! Skip all aerodynamic input

    allocate(turbine,STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadTurbineConfig 1')
       return
    end if

    call nullifyTurbine (turbine)

    !! Default values
    id = 0; extId = 0; extDescr = ''; ADFile = ''
    PtfmRef = 0.0_dp; HubRad = 0.0_dp; TipRad = 0.0_dp
    acOffset = 0.0_dp; ADcentre = 0.0_dp
    hubId = 0; generatorJoint = 0; pitchJoint = 0; nBlade = 0; nTB = 0
    firstTriadId = 0; towerTriad = 0; nacelleTriad = 0
    shaftTriad = 0; azimuthTriad = 0; hubTriad = 0; ADnodes = 1
    CompAero = .true.; CompNoise = .false.; SumPrint = .false.; UserID = .true.

    read(infp,nml=TURBINE_CONFIG,iostat=stat)
    if (stat == 0) then

       call initId (turbine%id,id,extId,extDescr,stat)
       turbine%PtfmRef = PtfmRef
       turbine%HubRad = HubRad
       turbine%TipRad = TipRad

       if (nBlade < 1 .or. nTB < 2) then
          nBlade = 0
          nBelm = 0
          err = err - 1
          call ReportInputError('TURBINE_CONFIG',id=turbine%id, &
               msg='Invalid turbine configuration')
       else if (nBlade > size(firstTriadId) .or. nBlade > size(pitchJoint)) then
          err = internalError('ReadTurbineConfig: Too many turbine blades'// &
               &              StrId(nBlade))
          return ! Must abort since the firstTriad array has been overwritten
       else if (nTB-1 > maxTB_p) then
          err = internalError('ReadTurbineConfig: Too many triads per blade'// &
               &              StrId(nTB))
          return ! Must abort since the acOffset array has been overwritten
       else
          ADnodes(nTB) = 0 ! Ensure the last blade triad is never an AD-node
          nBelm = count(ADnodes(1:nTB) == 1) ! Number of AD-elements per blade
       end if

       allocate(MENC(size(sups),2), MNEC(size(triads),2), &
            &   turbine%blade(nBlade), turbine%node(nBlade,nBelm), &
            &   turbine%bladeElm(nBelm), bldTriads(nTB,nBlade), STAT=stat)
       if (stat /= 0) then
          err = allocationError('ReadTurbineConfig 2')
          return
       end if

       !! Set up temporary topology arrays for the blades
       MENC = 0
       MNEC = 0
       stat = err
       do i = 1, size(sups)
          if (IsBeam(sups(i))) then
             triad => GetPtrToId(triads,sups(i)%triads(1)%p%id%baseId,MENC(i,1))
             if (.not. associated(triad)) err = err - 1
             triad => GetPtrToId(triads,sups(i)%triads(2)%p%id%baseId,MENC(i,2))
             if (.not. associated(triad)) err = err - 1
             do j = 1, 2
                n = MENC(i,j)
                if (n > 0) then
                   if (MNEC(n,1) == 0) then
                      MNEC(n,1) = i
                   else if (MNEC(n,2) == 0) then
                      MNEC(n,2) = i
                   else
                      !! Probably not a blade element
                      if (MENC(i,1) > 0) MNEC(MENC(i,1),:) = 0
                      if (MENC(i,2) > 0) MNEC(MENC(i,2),:) = 0
                      MENC(i,:) = 0
                   end if
                end if
             end do
          end if
       end do
       if (err < stat) then
          CompAero = .false.
          call ReportInputError('TURBINE_CONFIG',id=turbine%id, &
               msg='Detected'//trim(StrId(stat-err))//' blade topology errors')
          stat = err
       end if

       !! Detect some key turbine configuration triads
       turbine%tower => GetPtrToId(triads,towerTriad,userId=UserID)
       if (.not. associated(turbine%tower)) err = err - 1
       if (nacelleTriad > 0) then
          turbine%nacelle => GetPtrToId(triads,nacelleTriad,userId=UserID)
          if (.not. associated(turbine%nacelle)) err = err - 1
       end if
       if (shaftTriad > 0) then
          turbine%shaft => GetPtrToId(triads,shaftTriad,userId=UserID)
          if (.not. associated(turbine%shaft)) err = err - 1
       end if
       if (azimuthTriad > 0) then
          turbine%azimuth => GetPtrToId(triads,azimuthTriad,userId=UserID)
          if (.not. associated(turbine%azimuth)) err = err - 1
       end if
       turbine%hub => GetPtrToId(triads,hubTriad,userId=UserID)
       if (.not. associated(turbine%hub)) err = err - 1

       !! Traverse the blade topology to detect all blade nodes
       do i = 1, nBlade
          nullify(turbine%blade(i)%cone)
          nullify(turbine%blade(i)%tip)
          nullify(turbine%blade(i)%pitch)
          do j = 1, nBelm
             nullify(turbine%node(i,j)%triad)
             turbine%node(i,j)%aeroForce = 0.0_dp
             turbine%node(i,j)%aeroForcePrev = 0.0_dp
          end do

          id = firstTriadId(i)
          turbine%blade(i)%triad => GetPtrToId(triads,id,n,UserID)
          if (.not. associated(turbine%blade(i)%triad)) then
             err = err - (nBelm+1)
             cycle
          end if

          inod = 0
          do j = 1, nTB
             if (MNEC(n,1) > 0 .and. MNEC(n,2) == 0) then
                id = MNEC(n,1)
                MNEC(n,1) = 0
             else if (MNEC(n,2) > 0 .and. MNEC(n,1) == 0) then
                id = MNEC(n,2)
                MNEC(n,2) = 0
             else
                err = err - (nBelm-inod)
                exit
             end if
             if (MENC(id,1) == n) then
                n = MENC(id,2)
                MENC(id,:) = 0
             else if (MENC(id,2) == n) then
                n = MENC(id,1)
                MENC(id,:) = 0
             else
                err = err - (nBelm-inod)
                exit
             end if
             if (ADnodes(j) == 1) then
                inod = inod + 1
                turbine%node(i,inod)%triad => triads(n)
                triads(n)%aeroForce => turbine%node(i,inod)%aeroForce
                allocate(triads(n)%ur_def(12),STAT=stat)
                if (stat /= 0) then
                   err = allocationError('ReadTurbineConfig 3')
                   return
                end if
                triads(n)%ur_def = 0.0_dp
             end if
             if (MNEC(n,1) == id) then
                MNEC(n,1) = 0
             else if (MNEC(n,2) == id) then
                MNEC(n,2) = 0
             else
                err = err - (nBelm-inod)
                exit
             end if
             if (j == nTB .and. ADnodes(j) == 0) then
                !! For the blade tip we compute deflections only (no AeroDyn)
                turbine%blade(i)%tip => triads(n)
                allocate(triads(n)%ur_def(12),STAT=stat)
                if (stat /= 0) then
                   err = allocationError('ReadTurbineConfig 4')
                   return
                end if
                triads(n)%ur_def = 0.0_dp
             end if
             bldTriads(j,i)%p => triads(n)
          end do
          if (inod /= nBelm .and. err == stat) then
             err = err - abs(nBelm-inod)
          end if
       end do

       if (err < stat) then
          CompAero = .false.
          call ReportInputError('TURBINE_CONFIG',id=turbine%id, &
               msg='Failed to resolve'//trim(StrId(stat-err))//' triads')
          stat = err
       end if

       deallocate(MENC,MNEC)

       !! Define the generator DOFs
       if (UserID) then
          generator => GetPtrToId(joints,generatorJoint,jointType=REVOLUTE_p)
       else
          generator => GetPtrToId(joints,generatorJoint)
       end if
       if (.not. associated(generator)) then
          err = err - 1
       else if (generator%type /= REVOLUTE_p) then
          err = err - 1
       else
          turbine%genDOF => generator%jointDOfs(1)%jvar
          if (shaftTriad < 1) turbine%shaft => generator%JMTriads(1)%triad
          if (azimuthTriad < 1) turbine%azimuth => generator%STriad
       end if

       if (hubId > 0) then
          turbine%hubEl => GetPtrToId(sups,hubId,UserID)
          if (.not. associated(turbine%hubEl)) then
             hubId = 0
             err = err - 1
             call ReportInputError('TURBINE_CONFIG',id=turbine%id, &
                  msg='Failed to resolve hub part')
             stat = err
          end if
       end if

       !! Initialize some blade parameters
       do i = 1, nBlade
          if (UserID) then
             joint => GetPtrToId(joints,pitchJoint(i),jointType=REVOLUTE_p)
          else
             joint => GetPtrToId(joints,pitchJoint(i))
          end if
          if (.not. associated(joint)) then
             err = err - 1
          else if (joint%type /= REVOLUTE_p) then
             err = err - 1
          else
             turbine%blade(i)%cone  => joint%JMTriads(1)%triad
             turbine%blade(i)%pitch => joint%jointDOfs(1)%jvar(1)
          end if
          if (err < 0) then
             turbine%blade(i)%precone = 0.0_dp
          else if (i == 1) then
             !! Compute the precone angle for blade 1
             lpu = getErrorFile()
             turbine%blade(i)%precone = getPreCone(joint,turbine%hub)
             if (err < 0) then
                call ReportInputError('TURBINE_CONFIG',id=turbine%id, &
                     msg='First blade must have zero azimuth initially')
             end if
          else
             !! Assume all blades have the same precone angle
             turbine%blade(i)%precone = turbine%blade(i-1)%precone
          end if
       end do

       if (err < stat) then
          CompAero = .false.
          call ReportInputError('TURBINE_CONFIG',id=turbine%id, &
               msg='Failed to resolve'//trim(StrId(stat-err))//' joints')
          stat = err
       end if

       if (any(abs(ADcentre) > 1.0e-15_dp)) then
          do j = 1, nBelm
             turbine%bladeElm(j)%ACoffset = ADcentre(:,j)
          end do
       else
          do j = 1, nBelm
             turbine%bladeElm(j)%ACoffset(1) = acOffset(j)
             turbine%bladeElm(j)%ACoffset(2) = 0.0_dp
          end do
       end if

       !! Calculate some geometry parameters for each blade element.
       !! AeroDyn requires all blades to have identical configurations.
       R = 0.0_dp
       do i = 1, nBlade
          R = HubRad
          a = turbine%blade(i)%triad%ur(:,4) ! Blade root location
          do j = 1, nBelm
             if (.not. associated(turbine%node(i,j)%triad)) exit
             b  = turbine%node(i,j)%triad%ur(:,4)
             v  = b - a
             DR = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
             R  = R + DR
             DR = 2.0_dp*DR
             if (j > 1) then
                DR = DR - turbine%bladeElm(j-1)%DR
             end if
             if (i == 1) then
                turbine%bladeElm(j)%DR = DR
                turbine%bladeElm(j)%R  = R
             else if (abs(turbine%bladeElm(j)%R-R) > 1.0e-4_dp) then
                err = err - 1
             end if
             a = b
          end do
          R = R + 0.5_dp*turbine%bladeElm(nBelm)%DR
       end do
       if (err < stat) then
          call ReportInputError('TURBINE_CONFIG',id=turbine%id, &
               msg='Blades must have identical aerodynamic center positions')
       end if

       if (abs(turbine%genDOF(2)) > 1.0e-15_dp .and. err == 0) then
          !! Calculate initial velocity at all triads rotating with the turbine
          call initiateTurbineVelocity (turbine%hubEl,turbine%blade,bldTriads, &
               &                        generator%JPosInG,turbine%genDOF(2),err)
          if (err < 0) return
       end if
       deallocate(bldTriads)

       !! Define the tip radius (unless user specified it explicitly)
       if (TipRad <= 0.0_dp) turbine%TipRad = R

       !! Process the AeroDyn input (separate input file)
       call AeroInput (env,turbine,trim(ADFile),nBlade, &
            &          CompAero,CompNoise,SumPrint,stat)
       if (stat /= 0) err = err - 1
       if (.not. CompAero) then
          !! Deactivate AeroDynamic calculation in this run
          deallocate(turbine%node)
          allocate(turbine%node(0,0))
       end if

    else
       err = -1
       call ReportInputError('TURBINE_CONFIG')
    end if

    !! Read in the user-defined wind sampling points, if any

    n = iuGetNumberOfEntries(infp,'&WIND_POINT',stat)
    if (stat /= 0) then
       err = err - abs(stat)
    else if (n > 0) then
       write(lterm,*) 'Number of &WIND_POINT =',n
       allocate(turbine%windPt(3,n),STAT=stat)
       if (stat /= 0) then
          err = allocationError('ReadTurbineConfig 5')
          return
       end if
       do i = 1, n
          if (iuSetPosAtNextEntry(infp,'&WIND_POINT')) then
             coord = 0.0_dp
             read(infp,nml=WIND_POINT,iostat=stat)
             if (stat == 0) then
                turbine%windPt(:,i) = coord
                cycle
             end if
          end if
          err = err - 1
          call ReportInputError('WIND_POINT',i)
       end do
    else
       !! Allocate to zero length, such that the size() operator still works
       allocate(turbine%windPt(0,0))
    end if

    if (err < 0) call reportError (debugFileOnly_p,'ReadTurbineConfig')

  contains

    !> @brief Computes the precone angle for a wind turbine blade.
    function getPreCone (pitchJoint,hubTriad)
      use rotationModule   , only : deltaRot
      use manipMatrixModule, only : writeObject

      type(MasterSlaveJointType), intent(in) :: pitchJoint
      type(TriadType)           , intent(in) :: hubTriad

      real(dp)            :: getPreCone, coneAng(3), hub(3,3)
      real(dp), pointer   :: cone(:,:)
      real(dp), parameter :: eps_p = 1.0e-6_dp

      cone => pitchJoint%JMTriads(1)%triad%ur
      hub(:,1) = hubTriad%ur(:,3)
      hub(:,2) = hubTriad%ur(:,1)
      hub(:,3) = hubTriad%ur(:,2)

      !! Account for possible opposite orientation of the cone and hub triads
      if (dot_product(hub(:,2),cone(:,2)) < 0.0_dp) then
         hub(:,1) = -hub(:,1)
         hub(:,2) = -hub(:,2)
      end if
      coneAng = deltaRot(hub,cone)
      getPreCone = coneAng(2)

      !! Verify that the hubTriad and pitchJoint have common Y-axis
      if (abs(coneAng(1)) < eps_p .and. abs(coneAng(3)) < eps_p) return

      err = err - 1
      call writeObject(hub,lpu,'Triad'//getId(hubTriad%id))
      call writeObject(cone,lpu,'Triad'//getId(pitchJoint%JMTriads(1)%triad%id))
      write(lpu,"('Resulting cone angles =',1P,3E12.5)") coneAng
    end function getPreCone

  end subroutine ReadTurbineConfig


  !!============================================================================
  !> @brief Calculates initial velocity for the wind turbine triads.
  !>
  !> @param[in] hub Superelement representing the hub of the wind turbine rotor.
  !> @param blades Root nodes of thw wind turbine blades
  !> @param[in] bladeTriads Array of pointers to all triads of the wind turbine.
  !> @param[in] urGen Coordinate system for the generator triad
  !> @param[in] omega0 Initial generator angular velocity
  !> @param[out] ierr Error flag
  !>
  !> @details The initial velocity at all triads rotating with the turbine is
  !> calculated from the given initial generator velocity @a omega0.
  !> The corresponding centripetal acceleration is also computed.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Jan 2010

  subroutine initiateTurbineVelocity (hub,blades,bladeTriads,urGen,omega0,ierr)

    use SupElTypeModule      , only : SupElType, dp
    use WindTurbineTypeModule, only : BladeRoot
    use TriadTypeModule      , only : TriadPtrType
    use IdTypeModule         , only : getId
    use ReportErrorModule    , only : allocationError, getErrorFile

    type(SupElType)   , pointer       :: hub
    type(BladeRoot)   , intent(inout) :: blades(:)
    type(TriadPtrType), intent(in)    :: bladeTriads(:,:)
    real(dp)          , intent(in)    :: urGen(:,:), omega0
    integer           , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, n, lpu
    real(dp) :: x(3), u(3)
    type(TriadPtrType), allocatable :: triads(:)

    !! --- Logic section ---

    n = size(blades) + size(bladeTriads)
    if (associated(hub)) n = n + size(hub%triads)
    allocate(triads(n),STAT=ierr)
    if (ierr < 0) then
       ierr = allocationError('initiateTurbineVelocity')
       return
    end if

    !! Find all triads that rotate with the turbine rotor
    n = 0
    do i = 1, size(blades)
       n = n + 1
       triads(n)%p => blades(i)%triad
    end do
    do j = 1, size(bladeTriads,2)
       do i = 1, size(bladeTriads,1)
          n = n + 1
          triads(n)%p => bladeTriads(i,j)%p
       end do
    end do
    if (associated(hub)) then
       do i = 1, size(hub%triads)
          n = n + 1
          triads(n)%p => hub%triads(i)%p
       end do
    end if

    !! Compute triad velocities
    lpu = getErrorFile()
    u(3) = 0.0_dp
    do i = 1, n
       if (triads(i)%p%dependent) continue ! Skip dependent triads

       !! Relative position in generator coordinate system
       x = matmul(triads(i)%p%ur(:,4) - urGen(:,4), urGen(:,1:3))

       !! Translational velocity
       u(1) = -omega0*x(2)
       u(2) =  omega0*x(1)
       triads(i)%p%urd(1:3) = matmul(urGen(:,1:3),u)

       !! Angular velocity
       if (triads(i)%p%nDOFs >= 6) triads(i)%p%urd(4:6) = omega0*urGen(:,3)

       !! Centripetal acceleration
       u(1:2) = -x(1:2); u(3) = 0.0_dp
       triads(i)%p%urdd(1:3) = omega0*omega0*matmul(urGen(:,1:3),u)

       write(lpu,600) trim(getId(triads(i)%p%id))
       write(lpu,610) x, triads(i)%p%urd, triads(i)%p%urdd
    end do

    deallocate(triads)

600 format(/5X,'>>> Computed initial velocity and acceleration for Triad',A)
610 format( 9X,'X =',3F6.2,' urd  = ',1P,6E13.5 / 30X,' urdd = ',6E13.5)

  end subroutine initiateTurbineVelocity

end module initiateWindTurbineModule
