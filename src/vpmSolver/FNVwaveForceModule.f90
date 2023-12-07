!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module FNVwaveForceModule

  !!============================================================================
  !! This module contains data structures and associated subroutines for
  !! administrating higher-order wave force calculations on vertical cylinders
  !! in Fedem, based on the FNV-method.
  !!
  !! The core of the computations is found in FNVmodule.
  !!============================================================================

  use TriadTypeModule    , only : TriadType, dp
  use SupElTypeModule    , only : SupElPtrType
  use FiniteElementModule, only : BeamType
  use FNVmodule          , only : MatrixPtr

  implicit none

  type FNVColumnType
     real(dp)                    :: X(2)      !! Horizontal position of column
     real(dp)                    :: D         !! Column diameter
     type(TriadType)   , pointer :: surface   !! The triad at the water surface
     type(SupElPtrType), pointer :: elm(:)    !! Vertical chain of beam elements
     type(BeamType)    , pointer :: beam      !! Auxilliary pointer
     real(dp)          , pointer :: wave(:,:) !! Wave kinematics data
     real(dp)          , pointer :: F3Tot(:)  !! 3rd-order wave force history
     real(dp)          , pointer :: M4Tot(:)  !! 4th-order wave moment history
     real(dp)          , pointer :: F12t(:,:) !! 1st+2nd-order distributed force
     real(dp)          , pointer :: z(:)      !! Nodal Z-coordinates w.r.t. MSL
     type(MatrixPtr)   , pointer :: seakin(:) !! Sea kinematics at nodal points
  end type FNVColumnType

  type(FNVColumnType), save, pointer :: columns(:) => null()
  real(dp), save, private, pointer :: waveData(:,:) => null()

  real(dp), save, private :: g     !! Gravitation constant
  real(dp), save, private :: LPcut !! Cut-off frequency in LP-filter


contains

  subroutine findColumns (Tsea,beamProp,ierr)

    !!==========================================================================
    !! Detects the structural parts that are subjected to FNV wave forces.
    !! That is, chains of vertical beam elements with hydrodynamics enabled.
    !! It is assumed that the beam elements of each column are numbered from
    !! the bottom and up towards the free surface.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Oct 2011/1.0
    !!==========================================================================

    use FiniteElementModule, only : BeamPropertyType, GetPtrToId, beams
    use IdTypeModule       , only : getId
    use reportErrorModule  , only : allocationError, getErrorFile
    use reportErrorModule  , only : reportError, error_p

    real(dp)              , intent(in)  :: Tsea(3,4)
    type(BeamPropertyType), intent(in)  :: beamProp(:)
    integer               , intent(out) :: ierr

    !! Local variables
    integer             :: i, j, k, lpu, iCol, nCol, nElm, nNod
    character(len=128)  :: errMsg
    real(dp)            :: n(3), dX(3), X2(3)
    real(dp), parameter :: epsVer_p = 0.01_dp
    type(BeamType)        , pointer :: current
    type(BeamPropertyType), pointer :: bProp
    type(TriadType)       , pointer :: p0

    !! --- Logic section ---

    n = Tsea(:,3) ! Sea surface normal vector, defines the vertical direction

    !! Count the vertical beam columns with hydrodynamic loads
    nCol = 0
    current => beams
    do while (associated(current))
       nullify(p0)
       dX = 0.0_dp
       do i = 1, size(current%elm)
          if (associated(current%elm(i)%p%hydyn)) then
             if (.not. associated(p0)) then
                p0 => current%elm(i)%p%triads(1)%p
             end if
             dX = current%elm(i)%p%triads(2)%p%ur(:,4) - p0%ur(:,4)
          else if (associated(p0)) then
             exit
          end if
       end do
       if (associated(p0)) then
          !! Check that this column candidate is vertical and its elements are
          !! ordered from the bottom and upwards, and not the other way around
          if (cross_product(dX,n) < epsVer_p*dot_product(dX,n)) nCol = nCol + 1
       end if
       current => current%next
    end do

    if (associated(columns)) deallocate(columns)
    allocate(columns(nCol),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('findColumns 1')
       return
    else if (nCol == 0) then
       return
    end if

    !! Find the size (in terms of number of beam elements) of each column
    j = 1
    nCol = 0
    current => beams
    do while (associated(current))
       nullify(p0)
       dX = 0.0_dp
       nElm = 0
       do i = 1, size(current%elm)
          if (associated(current%elm(i)%p%hydyn)) then
             if (.not. associated(p0)) then
                p0 => current%elm(i)%p%triads(1)%p
             end if
             dX = current%elm(i)%p%triads(2)%p%ur(:,4) - p0%ur(:,4)
             nElm = nElm + 1
          else if (associated(p0)) then
             exit
          end if
       end do
       if (nElm > 0 .and. cross_product(dX,n) < epsVer_p*dot_product(dX,n)) then
          nCol = nCol + 1
          columns(nCol)%beam => current
          columns(nCol)%X = 0.0_dp
          columns(nCol)%D = 0.0_dp
          nullify(columns(nCol)%surface)
          nullify(columns(nCol)%wave)
          nullify(columns(nCol)%F3Tot)
          nullify(columns(nCol)%M4Tot)
          nullify(columns(nCol)%F12t)
          nullify(columns(nCol)%z)
          nullify(columns(nCol)%seakin)
          allocate(columns(nCol)%elm(nElm),STAT=ierr)
          if (ierr /= 0) then
             ierr = allocationError('findColumns 2')
             return
          end if
       end if
       current => current%next
       j = j + 1
    end do

    lpu = getErrorFile()

    !! Find the column elements and the horizontal position
    do iCol = 1, nCol
       current => columns(iCol)%beam
       do i = 1, size(current%elm)
          if (associated(current%elm(i)%p%hydyn)) then
             dX = glob2loc(Tsea,current%elm(i)%p%triads(1)%p%ur(:,4))
             nNod = 0
             do j = 1, size(columns(iCol)%elm)
                k = i+j-1
                columns(iCol)%elm(j)%p => current%elm(k)%p
                if (columns(iCol)%D <= 0.0_dp) then
                   !! Find column diameter from the beam properties
                   bProp => GetPtrToId(beamProp,-current%elm(k)%p%rigidFlag)
                   if (associated(bProp)) columns(iCol)%D = bProp%Morison(4)
                end if
                X2 = glob2loc(Tsea,current%elm(k)%p%triads(2)%p%ur(:,4))
                dX = dX + X2
                if (X2(3) <= 0.0_dp) then
                   nNod = 1+j ! Number of nodes under MSL for this column
                   columns(iCol)%surface => current%elm(k)%p%triads(2)%p
                else if (.not. associated(columns(iCol)%surface)) then
                   ierr = -iCol
                   write(errMsg,610) iCol,k,X2(3), &
                        getId(current%elm(k)%p%triads(2)%p%id)
610                format('Column',I3,'  Node',I4,'  Z =',1PE13.5,'  Triad',A84)
                   call reportError (error_p,'The FNV columns must be '// &
                        &            'ordered from the bottom and up.',errMsg, &
                        &            addString='findColumns')
                   return
                end if
             end do
             if (nNod == 0) then ! This column is entirely above MSL, ignore
                allocate(columns(iCol)%seakin(0),columns(iCol)%F12t(0,0))
                exit
             end if

             columns(iCol)%X = dX(1:2) / real(j,dp) ! Horizontal column position

             !! Allocate storage for the distributed data
             allocate(columns(iCol)%z(nNod), &
                  &   columns(iCol)%seakin(nNod),STAT=ierr)
             if (ierr /= 0) then
                ierr = allocationError('findColumns 3')
                return
             end if

             !! Calculate z-coordinates for the nodes below MSL
             X2 = glob2loc(Tsea,columns(iCol)%elm(1)%p%triads(1)%p%ur(:,4))
             columns(iCol)%z(1) = X2(3)
             nullify(columns(iCol)%seakin(1)%p)
             do j = 2, nNod
                X2 = glob2loc(Tsea,columns(iCol)%elm(j-1)%p%triads(2)%p%ur(:,4))
                columns(iCol)%z(j) = X2(3)
                nullify(columns(iCol)%seakin(j)%p)
             end do
             exit
          end if
       end do

       if (associated(columns(iCol)%surface)) then
          write(lpu,600) i,columns(i)%X,columns(i)%D, &
               &trim(getId(columns(i)%surface%id)), &
               (trim(getId(columns(i)%elm(j)%p%id)),j=1,size(columns(i)%elm))
       end if
600    format(/5X,'>>> FNV column',I3,' <<<', &
            & /9X,'Position:',2F8.2, &
            & /9X,'Diameter:', F8.2, &
            & /9X,'Surface node: Triad',A, &
            & /9X,'Beam members: Beam',A,/(23X,'Beam',A))

    end do

  contains

    function cross_product (a,b)
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: cross_product, x, y, z
      x = a(2)*b(3) - a(3)*b(2)
      y = a(3)*b(1) - a(1)*b(3)
      z = a(1)*b(2) - a(2)*b(1)
      cross_product = sqrt(x*x + y*y + z*z)
    end function cross_product

    function glob2loc (Tlg,xg) result(xl)
      use manipMatrixModule, only : matmul34, invert34
      real(dp), intent(in) :: Tlg(3,4), xg(3)
      real(dp)             :: xl(3)
      xl = matmul34(invert34(Tlg),xg)
    end function glob2loc

  end subroutine findColumns


  subroutine InitiateWaveSpectrumFNV (sys,env,ierr)

    !!==========================================================================
    !! Initiates the wave spectrum used by the FNV wave force model.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Oct 2011/1.0
    !!==========================================================================

    use SystemTypeModule       , only : SystemType
    use EnvironmentTypeModule  , only : EnvironmentType
    use explicitFunctionsModule, only : SINUSOIDAL_p, WAVE_SINUS_p
    use waveFunctionsModule    , only : twoPi, initFUNC4, waveNumber, checkDepth
    use reportErrorModule      , only : reportError, error_p, note_p
    use reportErrorModule      , only : allocationError, getErrorFile
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_intValue
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble

    type(SystemType)     , intent(in)    :: sys
    type(EnvironmentType), intent(inout) :: env
    integer              , intent(out)   :: ierr

    !! Local variables
    integer           :: i, j, lpu, nFreq
    real(dp)          :: Hs, Tp, Omega0, dOmega, Tsim, realData(5)
    character(len=64) :: cmsg

    !! --- Logic section ---

    if (.not. associated(columns)) then ! No FNV calculations in this run
       ierr = 0
       return
    else if (sys%tEnd <= sys%tStart) then
       ierr = 0
       call reportError (note_p,'FNV calculation requires a finite time domain')
       return
    else if (.not. associated(env%waveFunc)) then
       ierr = -1
       call reportError (error_p,'No wave function specified', &
            &            addString='InitiateWaveSpectrumFNV')
       return
    else if (associated(sys%tIncEngine) .or. sys%varInc > 0) then
       ierr = -2
       call reportError (error_p,'Time step must be constant when using FNV', &
            &            addString='InitiateWaveSpectrumFNV')
       return
    end if

    realData = env%waveFunc%realParameters(1:5)

    call ffa_cmdlinearg_getdouble('FNVlength',Tsim)
    if (Tsim < sys%tEnd-sys%tStart) then
       Tsim = sys%tEnd - sys%tStart ! Simulation length
    else
       write(cmsg,"('Simulation length Tsim =',1PE10.3)") Tsim
       call reportError (note_p,'FNV calculation is based on '//cmsg)
    end if
    Omega0 = 0.0_dp     ! Start frequency
    dOmega = twoPi/Tsim ! Delta angular frequency
    select case (env%waveFunc%type)
    case (WAVE_SINUS_p)
       nFreq = int(Tsim/sys%tInc)/2+1 ! Number of frequency components
       Hs = realData(1) ! Significant wave height
       Tp = realData(2) ! Mean wave period
       deallocate(env%waveFunc%realParameters)
       nullify(env%waveFunc%realParameters) ! Replaced by FNV spectrum below
    case (SINUSOIDAL_p)
       nFreq = 2 ! Regular wave, just one single component + offset
       Hs = 2.0_dp*realData(3) ! Wave amplitude
       Tp = realData(2) ! Wave period
    case default
       ierr = -1
       call reportError (error_p,'Invalid wave function data type', &
            &            addString='InitiateWaveSpectrumFNV')
       return
    end select

    if (associated(waveData)) deallocate(waveData)
    allocate(waveData(4,nFreq),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateWaveSpectrumFNV')
       return
    end if

    g = abs(dot_product(env%Tsea(:,3),env%gravity)) ! Gravitation constant
    call ffa_cmdlinearg_getdouble('FNVcutoff',LPcut)
    if (LPcut <= 0.0_dp) then
       LPcut = sqrt(2.0_dp*g/Hs) ! Cut-off frequency in LP-filter [rad/s]
    end if

    if (nFreq == 2) then ! Regular wave
       waveData(:,1) =  0.0_dp
       waveData(1,1) =  realData(4)
       waveData(1,2) =  realData(3)
       waveData(2,2) =  realData(1)*twoPi
       waveData(3,2) = -realData(2)*twoPi
       waveData(4,2) =  waveNumber(waveData(2,2),g,env%seaDepth)
    else
       !! Generate the Jonswap wave spectrum, uniform frequency distribution
       if (env%waveFunc%intParameters(2) > 4) then
          call initFUNC4 (Hs,Tp,realData(5),-Omega0,dOmega,g, &
               &          -env%seaDepth,waveData,ierr, &
               &          env%waveFunc%intParameters(3))
       else
          call initFUNC4 (3,Hs,Tp,realData(5),Omega0,dOmega,g, &
               &          -env%seaDepth,waveData,ierr, &
               &          env%waveFunc%intParameters(3))
       end if
       if (ierr < 0 .and. ierr > -32) goto 900

       !! Detect wave components that are within the legal wave number range.
       !! They are used in drag force evaluations only.
       j = 0
       do i = 1, nFreq
          if ( waveData(1,i) > 0.0_dp .and. waveData(4,i) > 0.0_dp .and. &
               .not. checkDepth(env%seaDepth,waveData(4,i)) ) then
             if (j == 0) j = i
          else if (j > 0) then
             exit
          end if
       end do
       if (j == 0) then
          ierr = min(0,ierr) - 4
          call reportError (error_p,'Entire wave spectrum is zero', &
               &            addString='InitiateWaveSpectrumFNV')
       else if (i-j < nFreq) then
          lpu = getErrorFile()
          if (j > 1) then
             write(lpu,601) 1,j-1,waveData(1,1),waveData(1,j-1), &
                  &               waveData(4,1),waveData(4,j-1)
          end if
          if (i <= nFreq) then
             write(lpu,601) i,nFreq,waveData(1,i),waveData(1,nFreq), &
                  &                 waveData(4,i),waveData(4,nFreq)
          end if
601       format('  ** Ignoring wave components',I6,' -',I6,' : A = [', &
               & 1PE12.5,',',E12.5,']  k = [',E12.5,',',E12.5,']')
          env%waveFunc%intParameters(3) = i-j
          env%waveFunc%realParameters => arrayPt(waveData(:,j:i-1),4*(i-j))
       end if
    end if

    if (ffa_cmdlinearg_intValue('printFunc') == 1) then
       lpu = getErrorFile()
       write(lpu,610)
       write(lpu,"(I8,1P4E13.5)") (i,waveData(:,i),i=1,nFreq)
       write(lpu,*)
       call flush(lpu)
610    format(/5X,'>>> Wave spectrum for FNV <<<', &
            & /5X,'  i      A_i       omega_i      epsilon_i      k_i')
    end if

    if (ierr >= 0) return

900 call reportError (error_p,'Failed to initialize wave data for FNV', &
         &            addString='InitiateWaveSpectrumFNV')

  contains

    function arrayPt (array,nval)
      integer , intent(in)         :: nval
      real(dp), intent(in), target :: array(nval)
      real(dp), pointer            :: arrayPt(:)
      arrayPt => array
    end function arrayPt

  end subroutine InitiateWaveSpectrumFNV


  subroutine InitiateFNVkinematics (TsimIn,dT,depth,ierr)

    !!==========================================================================
    !! Evaluates the wave kinematics for FNV-calculations on the beam columns.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Oct 2011/1.0
    !!==========================================================================

    use FNVmodule             , only : ring_kin_lin, reg_kin_lin
    use fileUtilitiesModule   , only : getDBGfile
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble

    real(dp), intent(in)    :: TsimIn, dT, depth
    integer , intent(inout) :: ierr

    !! Local variables
    integer  :: i, j, n, lerr, iord, nStep, lpu
    real(dp) :: wAng, Tsim

    !! --- Logic section ---

    ierr = 0
    iord = ffa_cmdlinearg_intValue('FNV')/1000 ! Two first digits
    call ffa_cmdlinearg_getdouble('FNVlength',Tsim)
    if (Tsim < TsimIn) Tsim = TsimIn

    if (associated(WaveData) .and. associated(columns)) then
       lpu = 0
       wAng = 0.0_dp
       do i = 1, size(columns)
          if (.not. associated(columns(i)%z)) continue ! Ignore dry columns

          if (size(WaveData,2) > 2) then
             call ring_kin_lin (columns(i)%X(1),columns(i)%X(2),1.0_dp, &
                  &             wAng,g,depth,Tsim,dT,LPcut,iord, &
                  &             WaveData,columns(i)%Wave,lerr)
          else
             call reg_kin_lin (columns(i)%X(1),columns(i)%X(2),1.0_dp, &
                  &            wAng,g,depth,Tsim,dT, &
                  &            WaveData,columns(i)%Wave,lerr)
          end if
          ierr = ierr + lerr
          do j = 1, size(columns(i)%z)
             if (size(WaveData,2) > 2) then
                call ring_kin_lin (columns(i)%X(1),columns(i)%X(2), &
                     &             columns(i)%z(j),wAng,g,depth,Tsim,dT,LPcut, &
                     &             iord,WaveData,columns(i)%seakin(j)%p,lerr)
             else
                call reg_kin_lin (columns(i)%X(1),columns(i)%X(2), &
                     &            columns(i)%z(j),wAng,g,depth,Tsim,dT, &
                     &            WaveData,columns(i)%seakin(j)%p,lerr)
             end if
             ierr = ierr + lerr
          end do

          if (lpu > 0 .or. ffa_cmdlinearg_intValue('debug') <= 0) continue

          lpu = getDBGfile(22,'FNVkinematics.asc')
          j = size(columns(i)%seakin)
          nStep = size(columns(i)%Wave,1)
          write(lpu,"('#   Time',9X,'elevation',4X,'u',12X,'w')")
          write(lpu,"(1P,4E13.5)") (dT*real(n,dp),columns(i)%wave(n,8), &
               &                    columns(i)%seakin(j)%p(n,1), &
               &                    columns(i)%seakin(j)%p(n,2),n=1,nStep)
          close(lpu)
       end do

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'InitiateFNVkinematics')

  end subroutine InitiateFNVkinematics


  subroutine getWaveKinematics (istep,inod,waterMotion,stat)

    !!==========================================================================
    !! Returns the pre-evaluated wave kinematics at a given node and time step.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 7 Nov 2011/1.0
    !!==========================================================================

    use IdTypeModule     , only : StrId
    use reportErrorModule, only : reportError, warning_p

    integer , intent(in)  :: istep, inod
    real(dp), intent(out) :: waterMotion(:)
    integer , intent(out) :: stat

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    do i = 1, size(columns)
       do j = 1, size(columns(i)%seakin)
          if (columns(i)%elm(j)%p%triads(1)%p%samNodNum == inod) then
             if (istep > 0 .and. istep <= size(columns(i)%seakin(j)%p,1)) then
                stat = 3
                waterMotion(3) = columns(i)%wave(istep,8)        ! eta
                waterMotion(4) = columns(i)%seakin(j)%p(istep,1) ! uu
                waterMotion(6) = columns(i)%seakin(j)%p(istep,2) ! ww
                waterMotion(7) = columns(i)%seakin(j)%p(istep,3) ! dudt
                return
             end if
          end if
       end do
    end do

    stat = 3
    waterMotion = 0.0_dp
    if (istep > 1) return

    call reportError (warning_p,'No wave kinematics for node '//StrId(inod), &
         &            addString='getWaveKinematics')

  end subroutine getWaveKinematics


  subroutine calcFNVforces (TsimIn,rhow,ierr)

    !!==========================================================================
    !! Perform the FNV force calculations on the beam columns.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Oct 2011/1.0
    !!==========================================================================

    use FNVmodule             , only : FNV12order, FNV23order
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use fileUtilitiesModule   , only : getDBGfile
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble

    real(dp), intent(in)  :: TsimIn, rhow
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, lpu, lerr, n, nStep, useFNV
    real(dp) :: dt, Tsim
    real(dp), pointer :: F12Tot(:,:)

    !! --- Logic section ---

    ierr = 0
    useFNV = ffa_cmdlinearg_intValue('FNV')
    if (useFNV > 0 .and. associated(WaveData) .and. associated(columns)) then

       call ffa_cmdlinearg_getdouble('FNVlength',Tsim)
       if (Tsim < TsimIn) Tsim = TsimIn

       lpu = 0
       nullify(F12Tot)
       do i = 1, size(columns)
          if (.not. associated(columns(i)%z)) continue ! Ignore dry columns

          call FNV12order (columns(i)%z,columns(i)%D,rhow,columns(i)%seakin, &
               &           columns(i)%F12t,F12Tot,useFNV,lerr)
          ierr = ierr + lerr
          call FNV23order (columns(i)%D,rhow,g,LPcut,Tsim,WaveData, &
               &           columns(i)%Wave, &
               &           columns(i)%F3Tot,columns(i)%M4Tot,useFNV,lerr)
          ierr = ierr + lerr

          if (lpu > 0 .or. ffa_cmdlinearg_intValue('debug') <= 0) continue

          lpu = getDBGfile(23,'FNV123order.asc')
          nStep = size(columns(i)%Wave,1)
          dt = Tsim/real(nStep,dp)
          write(lpu,"('#   Time',9X,'F1Tot',8X,'F2Tot',8X,'F3Tot',8X,'M4Tot')")
          if (associated(columns(i)%F3Tot)) then
             if (associated(columns(i)%M4Tot)) then
                write(lpu,"(1P,5E13.5)") (dt*real(n,dp),F12Tot(n,:), &
                     &                    columns(i)%F3Tot(n), &
                     &                    columns(i)%M4Tot(n),n=1,nStep)
             else
                write(lpu,"(1P,4E13.5)") (dt*real(n,dp),F12Tot(n,:), &
                     &                    columns(i)%F3Tot(n),n=1,nStep)
             end if
          else
             write(lpu,"(1P,3E13.5)") (dt*real(n,dp),F12Tot(n,:),n=1,nStep)
          end if
          close(lpu)
       end do
       if (associated(F12Tot)) deallocate(F12Tot)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'calculateFNVforces')

  end subroutine calcFNVforces


  subroutine addInFNVwaveForces (Q,RF,Tsea,sam,istep,ierr)

    !!==========================================================================
    !! Adds the FNV wave forces into the system vectors Q and RF.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 7 Nov 2011/1.0
    !!==========================================================================

    use SamModule             , only : SamType, dp
    use AsmExtensionModule    , only : csAddNV
    use reportErrorModule     , only : reportError, debugFileOnly_p

    type(SamType), intent(in)    :: sam
    real(dp)     , intent(in)    :: Tsea(3,4)
    real(dp)     , intent(inout) :: Q(:), RF(:)
    integer      , intent(in)    :: istep
    integer      , intent(inout) :: ierr

    !! Local variables
    integer  :: i, j, inod, err, lerr
    real(dp) :: eV(6)

    !! --- Logic section ---

    if (.not. associated(columns)) return ! No FNV calculations in this run

    lerr = ierr
    do i = 1, size(columns)
       if (istep > 0 .and. istep <= size(columns(i)%F12t,1)) then
          eV(4:6) = 0.0_dp
          do j = 1, size(columns(i)%F12t,2)
             if (j > size(columns(i)%elm)) then
                !! The entire column is submerged, get the end node
                !! (this should be the same as columns(i)%surface%samNodNum)
                inod = columns(i)%elm(j-1)%p%triads(2)%p%samNodNum
             else
                inod = columns(i)%elm(j)%p%triads(1)%p%samNodNum
             end if
             eV(1:3) = columns(i)%F12t(istep,j)*Tsea(:,1) ! Force in X-direction
             call csAddNV (sam, -inod, eV, Q, RF, err)
             if (err /= 0) ierr = ierr - 1
          end do
          if (associated(columns(i)%F3Tot)) then
             eV(1:3) = columns(i)%F3Tot(istep)*Tsea(:,1) ! Force in X-direction
          else
             eV(1:3) = 0.0_dp
          end if
          if (associated(columns(i)%M4Tot)) then
             eV(4:6) = columns(i)%M4Tot(istep)*Tsea(:,2) ! Moment about Y-axis
          end if
          inod = columns(i)%surface%samNodNum
          call csAddNV (sam, -inod, eV, Q, RF, err)
          if (err /= 0) ierr = ierr - 1
       end if
    end do

    if (ierr < lerr) call reportError (debugFileOnly_p,'addInFNVwaveForces')

  end subroutine addInFNVwaveForces

end module FNVwaveForceModule
