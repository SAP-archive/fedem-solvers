!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine INFO( rtype, i, j, rm, output, errflag )

  !!==============================================================
  !! Mimick the functionality of the (obsolete) ADAMS INFO routine
  !!==============================================================

  use JD_tireRoutinesModule, only : markers, dp
  use RotationModule       , only : deltaRot, mat_to_vec
  use EngineRoutinesModule , only : EngineValue
#ifdef FT_DEBUG
  use dbgUnitsModule       , only : dbgTire
#endif

  implicit none

  character(*), intent(in) :: rtype
  integer, intent(in)  :: i, j, rm
  real(dp),intent(out) :: output(*)
  logical, intent(out) :: errflag

  !! Local variables
  real(dp) :: vec(3), rot(3)
  integer  :: ierr

  !! --- Logic section ---

  ierr    = 0
  errflag = .false.

  !TODO,bh: take in account offset here

  if ( rtype(1:3) == 'DIS' ) then
     if ( j > 0 ) then
        vec = markers(i)%pWCinG(:,4) - markers(j)%pWCinG(:,4)
        rot = deltaRot(markers(j)%pWCinG(:,1:3),markers(i)%pWCinG(:,1:3))
     else
        vec = markers(i)%pWCinG(:,4)
        call mat_to_vec(markers(i)%pWCinG(:,1:3),rot)
     end if

     !! If the road is moved around with engines we can "fake" this
     !! by moving the tire in stead
     if ( associated(markers(i)%xTrEng) ) then
        vec(1) = vec(1) - EngineValue(markers(i)%xTrEng,ierr)
     end if
     if ( associated(markers(i)%yTrEng) ) then
        vec(2) = vec(2) - EngineValue(markers(i)%yTrEng,ierr)
     end if
     if ( associated(markers(i)%zTrEng) ) then
        vec(3) = vec(3) - EngineValue(markers(i)%zTrEng,ierr)
     end if
     errflag = ierr < 0

  else if ( rtype(1:3) == 'VEL' ) then
     if ( j > 0 ) then
        !! Should correct for Z offset here
        vec  = markers(i)%pTriad%urd(1:3) - markers(j)%pTriad%urd(1:3)
        rot  = markers(i)%pTriad%urd(4:6) - markers(j)%pTriad%urd(4:6)
     else
        vec  = markers(i)%pTriad%urd(1:3)
        rot  = markers(i)%pTriad%urd(4:6)
     end if

  else  if (rtype(1:3) == 'ACC') then
     if ( j > 0 ) then
        vec  = markers(i)%pTriad%urdd(1:3) - markers(j)%pTriad%urdd(1:3)
        rot  = markers(i)%pTriad%urdd(4:6) - markers(j)%pTriad%urdd(4:6)
     else
        vec  = markers(i)%pTriad%urdd(1:3)
        rot  = markers(i)%pTriad%urdd(4:6)
     end if

  else
     errflag = .true.
  end if

  ! == Get result vector in requested coordinate system

  if ( rm == 0 ) then ! Global system
     output(1:3) = vec
     output(4:6) = rot
  else                ! requested wheel coordinate system
#ifdef FT_DEBUG
     write(dbgTire,"('WCinG directions:'/,1p,3(3e12.3,/))") &
          &   markers(rm)%pWCinG(1,1:3), &
          &   markers(rm)%pWCinG(2,1:3), &
          &   markers(rm)%pWCinG(3,1:3)
#endif
     output(1:3) = matmul(vec,markers(rm)%pWCinG(:,1:3))
     output(4:6) = matmul(rot,markers(rm)%pWCinG(:,1:3))
  end if

end subroutine INFO


subroutine SYSARY( rtype, data, nsize, states, nstates, errflag )

  !!===============================================================
  !! Mimick the functionality of the ADAMS SYSARY routine
  !!===============================================================

  use KindModule, only : dp

  implicit none

  character(*), intent(in) :: rtype
  integer, intent(in)  :: nsize, data(:)
  real(dp),intent(out) :: states(:)
  integer, intent(out) :: nstates
  logical, intent(out) :: errflag

  !! --- Logic section ---

  nstates = 0
  nstates = 1/nstates ! Cause a segmentation fault
  call INFO( rtype, data(1), data(2), data(3), states, errflag )
  nstates = 6

end subroutine SYSARY


subroutine AKISPL(xval,zval,id,iord,array,errflag)

  use KindModule, only : dp

  implicit none

  real(dp), intent(in)  :: xval, zval
  integer,  intent(in)  :: id, iord
  real(dp), intent(inout) :: array(:)
  logical,  intent(out) :: errflag

  !! --- Logic section ---

  errflag = iord /= 0

end subroutine AKISPL


subroutine JDTYRE(NDEV  , ISWTCH, JOBFLG, IDTYRE, TIME  , &
     &            DIS   , TRAMAT, ANGTWC, VEL   , OMEGA , &
     &            OMEGAR, NDEQVR, DEQVAR, NTYPAR, TYPARR, &
     &            NCHTDS, CHTDST, ROAD  , IDROAD, NROPAR, &
     &            ROPAR , NCHRDS, CHRDST, FORCEWC,TORQUEWC, &
     &            DEQINI, DEQDER, TYRMOD, NVARS , VARINF, &
     &            RSTIFF, RDAMP , YSTIFF, YDAMP , IERR  )

  !=============================================================================
  !
  !  Project     : JD-tyre
  !  Version     : FEDEM implementation
  !  Description : Main interface to DADS 9.0, selects appropriate tyre model
  !                Reads in tyre data according to DADS lay-out
  !                Adjusts tyre scaling factors according to data from DADS-GUI
  !
  !=============================================================================

  !... Declaration for DTYRE

  use KindModule           , only : dp, epsDiv0_p
  use TireTypeModule       , only : JD_TIRE_p
  use TireRoutinesModule   , only : rtm_gTires
  use JD_tireRoutinesModule, only : jdt_Allocate_arrays
  use JD_tireRoutinesModule, only : jdt_InitializeMarkerData
  use FileUtilitiesModule  , only : findUnitNumber, printOpenFiles
  use ManipMatrixModule    , only : cross_product
#ifdef FT_DEBUG
  use dbgUnitsModule       , only : dbgTire
#endif
  use ReportErrorModule    , only : reportError, error_p, debugFileOnly_p
  use ReportErrorModule    , only : allocationError
  use JDT_inp              , only : ilat, akl, clat, yse
  use JDT_fdm

  implicit none

  integer  :: NDEV  , ISWTCH, JOBFLG, IDTYRE, NDEQVR, NTYPAR
  integer  :: NCHTDS, IDROAD, NROPAR, NCHRDS, NVARS , IERR

  real(dp) :: TIME, ANGTWC, OMEGAR
  real(dp) :: DIS(3), TRAMAT(3,3), VEL(3), OMEGA(3)
  real(dp) :: DEQVAR(NDEQVR), TYPARR(NTYPAR), ROPAR(NROPAR)
  real(dp) :: FORCEWC(3), TORQUEWC(3)
  real(dp) :: DEQINI(NDEQVR), DEQDER(NDEQVR), VARINF(NVARS)
  real(dp) :: RSTIFF(2), RDAMP, YSTIFF, YDAMP(2)
  character(len=*) :: CHTDST, CHRDST, TYRMOD

  external ROAD

  !! Local variables
  logical  :: iflag, dflag
  integer  :: i, j, npar, idof, id, idIn, maxJDT
  real*8   :: fmag, par(10), results(10)
  real(dp) :: stiff(3), damp(3)
  real(dp) :: forceL(3), torqueL(3), forceG(3), torqueG(3), traJD(3,3)
  real(dp) :: dFL(3), dDispL(3), dVelL(3)
  real(dp) :: absSlipAngle, forceYalpha, yDampExp, radius

  character(len=64) :: errMsg

  logical,  save :: useNumDamp, useMaxDamp, useMinDamp, useRotDamp
  real(dp), save :: alpha, dampFactor, stiffFactor

  !! --- Logic section ---

  ierr  = 0

  if (JOBFLG.eq.1) then

     !... Initialisation calls  (JOBFLG=1, determine array sizes)

     forcePrevL  = 0.0_dp
     torquePrevL = 0.0_dp
     forcePrevLout = 0.0_dp
     torquePrevLout = 0.0_dp
     posPrevG = 0.0_dp
     velPrevG = 0.0_dp
     call ffa_cmdlinearg_getdouble ('JDTalpha',alpha)
     call ffa_cmdlinearg_getdouble ('JDTdampFactor',dampFactor)
     call ffa_cmdlinearg_getdouble ('JDTstiffFactor',stiffFactor)
     call ffa_cmdlinearg_getint ('JDTdamping',id)
     useMaxDamp = mod(id,10) == 1
     useMinDamp = mod(id,10) == 2
     useNumDamp = mod(id,10) >= 3
     useRotDamp = id >= 10
     NDEQVR = 3   !TODO,bh: this should be zero really
     NVARS  = 35
     NTYPAR = 300

     if ( IDTYRE == 1 ) then
        do j = size(rtm_gTires), 1, -1
           if (rtm_gTires(j)%type == JD_TIRE_p) then
              !! Allocate internal arrays for the JD tires
              call jdt_allocate_arrays(j,ierr)
              exit
           end if
        end do
        if (ierr /= 0) ierr = AllocationError('jdt_allocate_arrays')
     end if

  else if (JOBFLG.eq.2) then

     !... Set device for error messages (entry point in DTEMSG routine)

     outFile = ndev

     j = IDTYRE
     maxJDT = size(hasReadTire)
     if (j < 1 .or. j > maxJDT) then
        write(errMsg,100) j, maxJDT
100     format('Invalid tire index',i4,'.  Maximum',i4,' JD-tires are allowed.')
        call reportError (error_p,errMsg,addString='JDTYRE')
        ierr = min(-1,-j)
        return
     end if

     !! Open the tire data file, using scratch file to avoid ctrl-M problems
     hasReadTire(j) = .false.
     tireFiles(j) = findUnitNumber(100+2*j)
     call openDataFile (tireFiles(j),CHTDST,'Tire',IDTYRE,ierr)
     if (ierr /= 0) return

     !! Open the road data file, avoid opening more than once
     if (NCHRDS == 0) then
        roadFiles(j) = 0
        FedemRoadIdin(j) = IDROAD
     else
        roadFiles(j) = findUnitNumber(100+2*j+1)
        FedemRoadIdin(j) = 0
        call openDataFile (roadFiles(j),CHRDST,'Road',IDROAD,ierr)
        if (ierr /= 0) return
     end if

     write(outFile,"('Initializing John Deere tire ',i3)") IDTYRE
     if ( j == size(tireFiles) ) then
        call printOpenFiles(1,150,outfile)
        do i = 1,size(tireFiles)
           write(outFile,"(' Tire file ',i2,' unit number ',i5)") i,tireFiles(i)
        end do
     end if

     call jdt_InitializeMarkerData (j,rtm_gTires(j),rtm_gTires(j)%road,ierr)
     if ( ierr /= 0 ) return

     dflag  = .false.
     iflag  = .false.
     par(1) = 0
     npar   = size(par)

     call SFOSUB(j, time, par(1), npar, dflag, iflag, fmag, ierr)
     if ( ierr /= 0 ) return

     if ( YSE_from_fedem(j) ) then
        yse(j) = DIS(2)
        write(outfile,"('Setting YSE(',i2,') to ',1p,e12.3)") j, yse(j)
     else
        write(outfile,"('Keeping YSE(',i2,') at ',1p,e12.3)") j, yse(j)
     end if

     TYPARR(7) = rollRadius(j)

     write(outFile,"('Finished initializing John Deere tire ',i3)") IDTYRE

  else

     dflag  = .false.
     iflag  = .false.
     par(1) = 0
     npar   = size(par)

     idIn   = IDTYRE

     do iDof = 1,6

        id = idIn + idof*100 !! Code to get correct force component
        call SFOSUB(id, time, par, npar, dflag, iflag, fmag, ierr)
        if ( ierr /= 0 ) return

        if      (iDof == 1)  then; forceL(idof)    = -fmag
        else if (iDof <= 3)  then; forceL(idof)    =  fmag
        else;                      torqueL(idof-3) =  fmag
        end if

     end do

     call JDT_stiff_XZ(idIn,time)

     !! Calculate the coordinate system for the JD routines
     traJD(:,3) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
     traJD(:,2) = tramat(:,2) ! Y-axis equal to Y-axis of WC system
     traJD(:,1) = cross_product(traJD(:,2),traJD(:,3))
     fmag = sqrt(dot_product(traJD(:,1),traJD(:,1)))
     traJD(:,1) = traJD(:,1)/fmag
     traJD(:,2) = cross_product(traJD(:,3),traJD(:,1))

!!$     forceL(2)  = (1.0_dp -alpha)*forceL(2)  + alpha*forcePrevL(2,idIn)
!!$     torqueL(1) = (1.0_dp -alpha)*torqueL(1) + alpha*torquePrevL(1,idIn)
!!$     force(1:3)  = (1.0_dp -alpha)*force(1:3)  + alpha*forcePrev(1:3,idIn)
!!$     torque(1:3) = (1.0_dp -alpha)*torque(1:3) + alpha*torquePrev(1:3,idIn)

     dFL    = forceL(1:3) - forcePrevL(:,idIn)
     dDispL = matmul(DIS - posPrevG(:,idIn),TRAMAT)
     dVelL  = VEL - matmul(velPrevG(:,idIn),TRAMAT)

     if ( ILAT(idIn) == 2 ) then
        forceL(2)   = -sign(forceL(2),VEL(2))
        forceYalpha = (1.0_dp-alpha)*forceL(2) + alpha*forcePrevLout(2,idIn)
        forceYalpha = -sign(forceYalpha,VEL(2))
     end if

     forcePrevL(1:3,idIn)  = forceL(1:3)
     torquePrevL(1:3,idIn) = torqueL(1:3)
     posPrevG(:,idIn)      = DIS
     velPrevG(:,idIn)      = matmul(TRAMAT,VEL)

     do i = 1,3
        if ( abs(dDispL(i)) > epsDiv0_p ) then; stiff(i) = -dFL(i)/dDispL(i)
        else;                                   stiff(i) =  0.0_dp
        end if
        if ( abs(dVelL(i))  > epsDiv0_p ) then; damp(i) = -dFL(i)/dVelL(i)
        else;                                   damp(i) =  0.0_dp
        end if
     end do

     !! Theoretical damping

     if ( abs(VEL(2)) > epsDiv0_p ) then
        if ( abs(yFricSlipAngle(idIn)) > 90.0_dp ) then
           absSlipAngle = abs( 180.0_dp - abs(yFricSlipAngle(idIn)) )
        else
           absSlipAngle = abs( yFricSlipAngle(idIn) )
        end if
        yDampExp = (yFricExp(idIn)*57.29_dp*forceL(3)*yFricCoeff(idIn) &
             &      *exp(-yFricExp(idIn)*absSlipAngle))/abs(VEL(2))
     else
        yDampExp = 0.0_dp
     end if

#ifdef FT_DEBUG
     write(dbgTire,"('===IdIn: ',i2,'  Time :',1p,e12.3,' =================')") idIn, time
     write(dbgTire,"('   Stiff   :',1p,3e12.3)") stiff(:)
     write(dbgTire,"('   Damping :',1p,3e12.3)") damp(:)
     write(dbgTire,"('   Damping Y theoretical:',1p,e12.3,' slipAngle:',e12.3)") &
          &         YDAMP(1), yFricSlipAngle(idIn)
     write(dbgTire,"('   ILAT ',i2)") ILAT(idIn)
#endif

     if ( ILAT(idIn) == 1 ) then
        YSTIFF   = AKL(idIn)
        YDAMP(1) = CLAT(idIn)
     else ! if ( ILAT(idIn) == 2 & 3 ) then

        YSTIFF = abs(stiff(2))*stiffFactor
        if ( useMaxDamp ) then
           YDAMP(1) = max(abs(damp(2)),yDampExp)*dampFactor
        else if ( useMinDamp ) then
           YDAMP(1) = min(abs(damp(2)),yDampExp)*dampFactor
        else if ( useNumDamp ) then
           YDAMP(1) = abs(damp(2))*dampFactor
        else
           YDAMP(1) = yDampExp*dampFactor
        end if

     end if

     radius = VARINF(12)
     if ( useRotDamp ) then
        YDAMP(2) = YDAMP(1)*radius*radius
     else
        YDAMP(2) = 0.0_dp
     end if

#ifdef FT_DEBUG
     write(dbgTire,"('   Using Ystiff, Ydamp, xRdamp :',1p3e12.3)") YSTIFF,YDAMP
#endif

     if ( ILAT(idIn) == 2 ) then
        forceL(2)  = forceYalpha
        torqueL(1) = forceL(2)*radius
     end if

     forcePrevLOut(:,IDTYRE) = forceL
     torquePrevLOut(:,IDTYRE) = torqueL

     !! Finally output in WC coordinate system

     forceG  =  forceL(1)*traJD(:,1) +  forceL(2)*traJD(:,2) +  forceL(3)*traJD(:,3)
     torqueG = torqueL(1)*traJD(:,1) + torqueL(2)*traJD(:,2) + torqueL(3)*traJD(:,3)

     forceWC  = matmul(forceG,tramat)   !! Output force in WC coordinate system
     torqueWC = matmul(torqueG,tramat)  !! Output torque in WC coordinate system


     !! Longitudinal data
     par(1) = 0
     id = idIn + 100
     results = 0.0_dp
     call reqsub(id, time, par, npar, iflag, results)

     VARINF( 2) = results(4) !! Longitudinal slip
     VARINF( 4) = results(8) !! FX
     VARINF(12) = results(6) !! RR, rolling radius
!bh     varinf(13) = DIS(1)     !! XDISP !!TODO,bh: look at this

     !! Lateral data
     id = idIn + 200
     results = 0.0_dp
     call reqsub(id, time, par, npar, iflag, results)

     !! Slip angle
     if (abs(results(4)) > 90.0_dp) then
        VARINF(1) = sign(180.0_dp-abs(results(4)),results(4))
     else
        VARINF(1) = results(4)
     end if

     VARINF(14) = DIS(2)     !! YDISP !!TODO,bh: look at this
     VARINF(14) = results(2) !! YDISP !!TODO,bh: look at this
     VARINF( 5) = results(7) !! FY

     !! Vertical data
     id = idIn + 300
     results = 0.0_dp
     call reqsub(id, time, par, npar, iflag, results)

     VARINF(13) = results(2) !! XDISP
     VARINF(15) = results(3) !! ZDISP !! need to subtract radius (downwards)
     RDAMP      = results(5) !! damping
     VARINF(10) = results(6) !! deflection, Essential that this is here !!!
     VARINF( 6) = results(7) !! FZ

     if ( abs(results(6)) > epsDiv0_p ) then
        RSTIFF(1) = abs(results(1)/results(6))
        RSTIFF(2) = 0.0_dp
     end if

#ifdef FT_DEBUG
     write(dbgTire,"('Radial stiffness:',1p2e12.3)") RSTIFF
     write(dbgTire,"('Time:',1p,e12.3,'  Forces:',6e12.3)") &
          time, forceL(1:3), torqueL(1:3)
     write(dbgTire,"('Time:',1p,e12.3,'  SlipAngle REQSUB, SFOSUB:',2e12.3)") &
          time, VARINF(1), yFricSlipAngle(idIn)
#endif

  endif

contains

  subroutine openDataFile (iunit,fname,objType,id,ierr)
    use IdTypeModule  , only : StrId
    use InputUtilities, only : iuCopyToScratch
    integer         , intent(in)  :: iunit, id
    character(len=*), intent(in)  :: fname, objType
    integer         , intent(out) :: ierr
    open(iunit,STATUS='SCRATCH',IOSTAT=ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open scratch input file for '// &
            &            objType//' #'//StrId(id))
    else
       !! Preprocess the data file, removing DOS-formatting, etc.
       call iuCopyToScratch(iunit,fname,ierr)
       if (ierr == 0) rewind(iunit)
    end if
    if (ierr /= 0) call reportError (debugFileOnly_p,'JDTYRE')
  end subroutine openDataFile

end subroutine JDTYRE


SUBROUTINE TIREF_FDM(J, time, delta, forceX, forceZ)

  ! TIREF USES A DISTRIBUTED CONTACT TIRE MODEL TO PREDICT THE FORCE
  ! ACTING ON A TIRE AS IT TRAVELS OVER ROUGH GROUND

  use jdt_var
  use roadRoutinesModule, only : RoadType, dp
  use tireRoutinesModule, only : rtm_gTires
#ifdef FT_DEBUG
  use dbgUnitsModule    , only : dbgTire
#endif

  IMPLICIT NONE

  integer :: j
  real*8  :: time
  real*8  :: delta(3), forceX, forceZ


  integer :: i, k, itrm1, nstam1, nseg, nsegd2, klast, iseg, nstart
  real*8  :: xtest1, slope, x2, y2, x3, y3, st, xx, yy, dist, xtest2, darea
  real*8  :: absfz, alph, fmag, fcorr
  real*8, external :: table1, table1_grad

  integer  :: iter, nRoad, ierr, dflag
  real(dp) :: tireCenter(3), roadContact(3), dz(2), ddz(3), prsurf
  real(dp) :: dS, dR, zUndefTire, xOnSpringSlope
  type(RoadType), pointer :: pRoad


  ! CALCULATE THE LONGITUDINAL AND VERTICAL TIRE SPRING FORCE COMPONENTS
  ! THE SEGMENTED TIRE MODEL

  NSEG=int(ANSEG(J))
  NSEGD2=NSEG/2

  ! CALCULATE THE LENGTH, TMAG(ISEG), OF EACH SEGMENT SPRING. IF THE SEGMENT
  ! SPRING IS NOT COMPRESSED BY THE TERRAIN, THIS LENGTH IS SET EQUAL TO THE
  ! UNDEFLECTED WHEEL RADIUS. IF THE SEGMENT SPRING IS COMPRESSED, THE
  ! LENGTH IS THE DISTANCE BETWEEN THE WHEEL CENTER AND THE POINT OF INTER-
  ! SECTION OF THE SPRING WITH THE TERRAIN.
  ! FIRST PERFORM THIS DETERMINATION FOR THE SEGMENT SPRINGS LYING BEHIND THE
  ! WHEEL CENTER

  if ( associated(rtm_gTires(j)%road%roadFunc) ) then

     pRoad => rtm_gTires(j)%road

     do i = 1,3
        tireCenter(i) = rtm_gTires(j)%WCinG(i,4) + delta(i)
     end do

#ifdef FT_DEBUG
     write(dbgTire,"('x1(j)        = ',1p,e12.3,'       y1(j)   = ',e12.3)") x1(j), y1(j)
     write(dbgTire,"('tireCenter_x = ',1p,e12.3,'  tireCenter_z = ',e12.3)")  &
          & tireCenter(1), tireCenter(3)
     write(dbgTire,"('deltaX       = ',1p,e12.3,'  deltaY       = ',2e12.3)")  &
          & (delta(i),i=1,3)
#endif

     !! Find lengt of all section springs for the tire omdel
     do iseg = 1,nseg

        !!TODO, bh: handle turning wheels here
        roadContact(1) = tireCenter(1) + wrad(j)*ssega(iseg,j)
        roadContact(2) = tireCenter(2)
        zUndefTire     = tireCenter(3) - wrad(j)*csega(iseg,j)
        dS             = sm(iseg,j)  !! Slope of spring
        dR             = dS

        do iter = 1,20
           call FedemRoadMain (time, roadContact(1), 1, 0, j, pRoad%idIn, &
                &              pRoad%nropar, pRoad%roparr(1),  &
                &              pRoad%nCharRoadDataFileName, pRoad%roadDataFileName, 0, &
                &              nroad, roadContact(3), dz, ddz, dflag, prsurf, ierr)

           if ( roadContact(3) <= zUndefTire ) exit !! No contact with road

           xOnSpringSlope = tireCenter(1) - (tireCenter(3)-roadContact(3))/dS

           if ( abs(xOnSpringSlope-roadContact(1)) < 0.00001d0*wrad(j)) then
              !! Road contact found to be close enough to spring line
              exit
           else
              !! Find improved x coordinate for intersection of spring and road
              if (abs(dR - dS) > 1.0d-5) then
                 if ( dFlag < 1 ) then; dR = 0.0d0
                 else;                  dR = dz(1)
                 end if

                 roadContact(1)  = ( (tireCenter(3)  - dS*tireCenter(1)) &
                      &            - (roadContact(3) - dR*roadContact(1)) ) &
                      &             /( dR - dS )
              else
                 !! Improve x coordinate can not be found
                 roadContact(1) = tireCenter(1) + wrad(j)*ssega(iseg,j)
                 roadContact(3) = zUndefTire
              exit
           end if
        end if

        end do

#ifdef FT_DEBUG
        write(dbgTire,"('nIter for tire = ',i2)")  iter
#endif

        dist=sqrt( (tireCenter(1)-roadContact(1))**2 &
             &    +(tireCenter(3)-roadContact(3))**2 ) !! Y-coordinate assumed the same

        if(dist > wrad(j)) then
           tmag(iseg)=wrad(j)
        else
           tmag(iseg)=dist
        endif
     end do

!!$!! bh Begin temp
!!$
!!$     ! FIND THE NUMBER, ITR, OF THE TERRAIN STATION IMMEDIATELY AHEAD OF THE
!!$     ! CENTER OF WHEEL J
!!$#ifdef FT_DEBUG
!!$     write(dbgTire,"('Getting NTSTA(',i2,') = ',I5)") j, ntsta(j)
!!$#endif
!!$     write(*      ,"('Getting NTSTA(',i2,') = ',I5)") j, ntsta(j)
!!$
!!$     ITRM1=ITR(J)-3
!!$     IF(ITR(J).LE.3) ITRM1=1
!!$     NSTAM1=NTSTA(J)-1
!!$     DO  I=ITRM1,NSTAM1
!!$        IF (X1(J).GE.STA(I,J).AND.X1(J).LE.STA(I+1,J)) then
!!$           ITR(J)=I+1
!!$           exit
!!$        end IF
!!$     end DO
!!$
!!$     KLAST=0
!!$     DO I=1,NSEGD2
!!$        ISEG=NSEGD2+1-I
!!$        DO K=1,500
!!$           XTEST1=X1(J)-STA(ITR(J)-KLAST-K,J)
!!$           IF(DABS(XTEST1).LE.0.01) cycle
!!$           slope=(Y1(J)-HGT(ITR(J)-KLAST-K,J))/XTEST1
!!$           IF(slope.GT.SM(ISEG,J)) cycle
!!$           X2=STA(ITR(J)-KLAST-K,J)
!!$           Y2=HGT(ITR(J)-KLAST-K,J)
!!$           X3=STA(ITR(J)-KLAST-K+1,J)
!!$           Y3=HGT(ITR(J)-KLAST-K+1,J)
!!$           ST=(Y3-Y2)/(X3-X2)
!!$           XX=(Y2-Y1(J)+SM(ISEG,J)*X1(J)-ST*X2)/(SM(ISEG,J)-ST)
!!$           YY=ST*(XX-X2)+Y2
!!$           DIST=DSQRT((X1(J)-XX)**2+(Y1(J)-YY)**2)
!!$           IF(DIST.GE.WRAD(J)) TMAG(ISEG)=WRAD(J)
!!$           IF(DIST.LT.WRAD(J)) TMAG(ISEG)=DIST
!!$           KLAST=KLAST+K-1
!!$           exit
!!$        end DO
!!$     end DO
!!$
!!$     ! NEXT PERFORM THIS DETERMINATION FOR THE SEGMENT SPRINGS LYING
!!$     ! AHEAD OF THE WHEEL CENTER
!!$
!!$     NSTART=NSEGD2+1
!!$     KLAST=-1
!!$
!!$     DO I=NSTART,NSEG
!!$        ISEG=I
!!$        DO K=1,500
!!$           XTEST2=STA(ITR(J)+KLAST+K,J)-X1(J)
!!$           IF(DABS(XTEST2).LE.0.01) cycle
!!$           slope=(Y1(J)-HGT(ITR(J)+KLAST+K,J))/XTEST2
!!$           IF(slope.GT.DABS(SM(ISEG,J))) cycle
!!$           X2=STA(ITR(J)+KLAST+K-1,J)
!!$           Y2=HGT(ITR(J)+KLAST+K-1,J)
!!$           X3=STA(ITR(J)+KLAST+K,J)
!!$           Y3=HGT(ITR(J)+KLAST+K,J)
!!$           ST=(Y3-Y2)/(X3-X2)
!!$           XX=(Y2-Y1(J)+SM(ISEG,J)*X1(J)-ST*X2)/(SM(ISEG,J)-ST)
!!$           YY=ST*(XX-X2)+Y2
!!$           DIST=DSQRT((X1(J)-XX)**2+(Y1(J)-YY)**2)
!!$           IF(DIST.GE.WRAD(J)) TMAG(ISEG)=WRAD(J)
!!$           IF(DIST.LT.WRAD(J)) TMAG(ISEG)=DIST
!!$           KLAST=KLAST+K-1
!!$           exit
!!$        end DO
!!$     end DO
!!$
!!$#ifdef FT_DEBUG
!!$     write(dbgTire,"(' ')")
!!$     write(dbgTire,"('======= Segment deflections new and old road file model',i2)") J
!!$     write(dbgTire,"('delta = ',1p,3e12.3)") (delta(k),k=1,3)
!!$     do iseg = 1,nseg
!!$        write(dbgTire,"('iseg = ',I2,1p,' tmag tmag2 =',2e12.3)") iseg, tmag(iseg),tmag2(iseg)
!!$     end do
!!$#endif
!!$
!!$!! bh End temp

  else

     ! FIND THE NUMBER, ITR, OF THE TERRAIN STATION IMMEDIATELY AHEAD OF THE
     ! CENTER OF WHEEL J

     ITRM1=ITR(J)-3
     IF(ITR(J).LE.3) ITRM1=1
     NSTAM1=NTSTA(J)-1
     DO  I=ITRM1,NSTAM1
        IF (X1(J).GE.STA(I,J).AND.X1(J).LE.STA(I+1,J)) then
           ITR(J)=I+1
           exit
        end IF
     end DO

     KLAST=0
     DO I=1,NSEGD2
        ISEG=NSEGD2+1-I
        DO K=1,500
           XTEST1=X1(J)-STA(ITR(J)-KLAST-K,J)
           IF(DABS(XTEST1).LE.0.01) cycle
           slope=(Y1(J)-HGT(ITR(J)-KLAST-K,J))/XTEST1
           IF(slope.GT.SM(ISEG,J)) cycle
           X2=STA(ITR(J)-KLAST-K,J)
           Y2=HGT(ITR(J)-KLAST-K,J)
           X3=STA(ITR(J)-KLAST-K+1,J)
           Y3=HGT(ITR(J)-KLAST-K+1,J)
           ST=(Y3-Y2)/(X3-X2)
           XX=(Y2-Y1(J)+SM(ISEG,J)*X1(J)-ST*X2)/(SM(ISEG,J)-ST)
           YY=ST*(XX-X2)+Y2
           DIST=DSQRT((X1(J)-XX)**2+(Y1(J)-YY)**2)
           IF(DIST.GE.WRAD(J)) TMAG(ISEG)=WRAD(J)
           IF(DIST.LT.WRAD(J)) TMAG(ISEG)=DIST
           KLAST=KLAST+K-1
           exit
        end DO
     end DO

     ! NEXT PERFORM THIS DETERMINATION FOR THE SEGMENT SPRINGS LYING
     ! AHEAD OF THE WHEEL CENTER

     NSTART=NSEGD2+1
     KLAST=-1

     DO I=NSTART,NSEG
        ISEG=I
        DO K=1,500
           XTEST2=STA(ITR(J)+KLAST+K,J)-X1(J)
           IF(DABS(XTEST2).LE.0.01) cycle
           slope=(Y1(J)-HGT(ITR(J)+KLAST+K,J))/XTEST2
           IF(slope.GT.DABS(SM(ISEG,J))) cycle
           X2=STA(ITR(J)+KLAST+K-1,J)
           Y2=HGT(ITR(J)+KLAST+K-1,J)
           X3=STA(ITR(J)+KLAST+K,J)
           Y3=HGT(ITR(J)+KLAST+K,J)
           ST=(Y3-Y2)/(X3-X2)
           XX=(Y2-Y1(J)+SM(ISEG,J)*X1(J)-ST*X2)/(SM(ISEG,J)-ST)
           YY=ST*(XX-X2)+Y2
           DIST=DSQRT((X1(J)-XX)**2+(Y1(J)-YY)**2)
           IF(DIST.GE.WRAD(J)) TMAG(ISEG)=WRAD(J)
           IF(DIST.LT.WRAD(J)) TMAG(ISEG)=DIST
           KLAST=KLAST+K-1
           exit
        end DO
     end DO

  end if
  ! CALCULATE THE SPRING FORCE ACTING IN EACH SEGMENT

  DAREA=0.0
  DO I=1,NSEG
     TMAG(I)=WRAD(J)-TMAG(I)
     IF(LTIRE.EQ.0) then
        TMAG(I) = TMAG(I)*ASEGK(I,J)
        SpringStiff(i) = asegk(i,j)
     else
        DAREA  = DAREA+(WRAD(J)*TMAG(I)-0.5*TMAG(I)**2) * DTHET(J)
        TMAG(I)= TABLE1(DEFL,FORCE,TMAG(I),2,NLDPT(J),1,J)
        SpringStiff(i) = table1_grad(defl,force,tmag(i),nldpt(j),j)
     end IF
  end DO

  ! SUM THE LONGITUDINAL AND VERTICAL COMPONENTS OF EACH SEGMENT SPRING FORCE
  ! TO FIND THE NET LONGITUDINAL AND VERTICAL TIRE SPRING FORCE COMPONENTS

  forceX=0.0
  forceZ=0.0
!!$  do ii = 1,3
!!$     do jj = 1,3
!!$        stiffM(ii,jj,J) = 0.0d0
!!$     end do
!!$  end do

  DO I=1,NSEG
     forceX=forceX-TMAG(I)*SSEGA(I,J)
     forceZ=forceZ+TMAG(I)*CSEGA(I,J)

!!$     !! Material stiffness
!!$     stiffM(1,1,J) = stiffM(1,1,J) + SpringStiff(i)*ssega(i,j)*ssega(i,j)
!!$     stiffM(3,3,J) = stiffM(3,3,J) + SpringStiff(i)*csega(i,j)*csega(i,j)
!!$     stiffM(1,3,J) = stiffM(1,3,J) - SpringStiff(i)*csega(i,j)*ssega(i,j)
!!$     stiffM(3,1,J) = stiffM(1,3,J)
!!$
!!$     !! Geometric stiffness
!!$     stiffGeo = -TMAG(j)/wrad(j)
!!$     stiffM(1,1,J) = stiffM(1,1,J) + stiffGeo*csega(i,j)*csega(i,j)
!!$     stiffM(3,3,J) = stiffM(3,3,J) + stiffGeo*ssega(i,j)*ssega(i,j)
!!$     stiffM(1,3,J) = stiffM(1,3,J) + stiffGeo*csega(i,j)*ssega(i,j)
!!$     stiffM(3,1,J) = stiffM(1,3,J)

  end DO

  ! IF THE INTERPOLATION METHOD IS TO BE USED, FIND THE NET LONGITUDINAL
  ! AND VERTICAL TIRE SPRING FORCE COMPONENTS

  if ( LTIRE /= 0 .or. forceZ > 0.0 ) then

     fCorr = 1.0_dp/forceZ

     ABSFZ=DABS(forceZ)
     ALPH=DATAN2(forceX,ABSFZ)
     FMAG=TABLE1(AREA,FORCE,DAREA,2,NLDPT(J),1,J)
     forceX=FMAG*DSIN(ALPH)
     forceZ=FMAG*DCOS(ALPH)

     fCorr = fCorr*forceZ
!!$     write(56,"('Deflection to area correction, fCorr = ',1p, e12.3)") fCorr
!!$     do ii = 1,3
!!$        do jj = 1,3
!!$           stiffM(ii,jj,J) = stiffM(ii,jj,j)*fCorr
!!$        end do
!!$     end do
  end if

  ! CHECK TO SEE IF THE VERTICAL TIRE SPRING FORCE IS NEGATIVE INDICATING
  ! THE TIRE HAS LOST CONTACT WITH THE GROUND. IF SO, SET ALL THE TIRE FORCE
  ! COMPONENTS TO ZERO.

  IF(forceZ.LE.0.0) THEN
     forceX=0.0
     forceZ=0.0
!!$     do ii = 1.3
!!$        do jj = 1.3
!!$           stiffM(ii,jj,J) = 0.0_dp
!!$        end do
!!$     end do
  ENDIF

  RETURN

END SUBROUTINE TIREF_FDM


function table1_grad(x,y,xarg,ndim,k) result(grad)
  !
  !  THIS FUNCTION USES LAGRANGIAN INTERPOLATION OF DEGREE IDEG TO
  !         DETERMINE A VALUE OF Y CORRESPONDING TO XARG, BASED UPON
  !         THE NDIM PAIRED VALUES OF X AND Y.
  !  EXAMPLE USE:
  !         TWO ARRAYS, TEMP(I) AND HUMID(I), IN SOME PROGRAM CONTAIN
  !             CONSECUTIVE VALUES OF TWO VARIABLES. THERE ARE NDIM
  !             OF EACH VARIABLE, WITH THOSE HAVING EQUAL SUBSCRIPTS
  !             CORRESPONDING TO ONE ANOTHER. TO OBTAIN THE INTERPOLATED
  !             VALUE, HVALUE, CORRESPONDING TO A SPECIFIC VALUE WITHIN
  !             THE TEMP(I) RANGE, SAY TVALUE, THE FOLLOWING ASSIGNMENT
  !             STATEMENT WOULD BE USED IN THE PROGRAM WHICH HAS THE TWO
  !             ARRAYS STORED:
  !                  HVALUE = TABLE(TEMP,HUMID,TVALUE,IDEG,NDIM,JMIN)
  !
  !  *** CAUTION ***   X AND Y MUST BE DIMENSIONED NDIM IN MAIN PROGRAM
  !
  use jdt_max, only : maxCPT

  IMPLICIT NONE

  integer :: ndim, k
  real*8  :: X(maxCPT,*),Y(maxCPT,*)
  real*8  :: xarg, grad

  integer :: i

  !  Search for the interval containing xarg

  grad = 0.0d0
  do i = 1,ndim-2
     if ( xarg >= x(i,k) ) then
        ! Compute gradient
        grad = (y(i+1,k) - y(i,k))/(x(i+1,k) - x(i,k))
     end if
  end do

end function table1_grad


subroutine JDT_stiff_XZ(J,time)

  use JDT_var, only : wrad, stiffM

  implicit none

  integer J
  real*8 time

  real*8 delta(3), fXp, fZp, fXm, fZm, dd

  dd = wrad(j)*0.0001d0

  stiffM(1:3,1:3,J) = 0.0d0

  !! Perturbation in X coordinate
  delta = 0.0d0
  delta(1) = -dd
  call tiref_fdm(j,time,delta,fXm,fZm)
  fXm = -fXm
  fZm = -fZm

  delta(1) = dd
  call tiref_fdm(j,time,delta,fXp,fZp)
  fXp = -fXp
  fZp = -fZp
  stiffM(1,1,j) = (fXp - fXm)/(2.0d0*dd)
  stiffM(3,1,j) = (fZp - fZm)/(2.0d0*dd)

  !! Perturbation in Z coordinate
  delta = 0.0d0
  delta(3) = -dd
  call tiref_fdm(j,time,delta,fXm,fZm)
  fXm = -fXm
  fZm = -fZm
  delta(3) = dd
  call tiref_fdm(j,time,delta,fXp,fZp)
  fXp = -fXp
  fZp = -fZp
  stiffM(1,3,j) = (fXp - fXm)/(2.0d0*dd)
  stiffM(3,3,j) = (fZp - fZm)/(2.0d0*dd)

!bh  stiffM(1,3,j) = (stiffM(1,3,j) + stiffM(3,1,J))*0.5d0
!bh  stiffM(3,1,j) = stiffM(1,3,j)

end subroutine JDT_stiff_XZ
