!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module TireRoutinesModule

  use TireTypeModule, only : dp, TireType, ScaleToType

  implicit none

  type WCdata
     real(dp) :: time
     real(dp) :: pos(3)
     real(dp) :: trans(3,3)
  end type WCdata

  type(ScaleToType), save :: scaleTo ! Scaling factors from modeling units to SI

  private :: WCdata
  private :: calcTireForces, calcInertiaForces
  private :: STI, initiateSTIstateVar, updateSTIstateVar

  !! Global array of all tires. To be used as a back door into the data
  !! structure for some tire routines (John Deere tire), which have to have a
  !! preset argument list and be F77-callable.
  type(TireType), save, pointer :: rtm_gTires(:)


contains

  subroutine addInTireStiffMat (lDynamics,scaleK,Nmat,tire,sam,ierr,Rhs)

    !!==========================================================================
    !! Calculates the stiffness matrix for a tire and adds to
    !! the system Newton matrix, multiplied by a factor, scaleK.
    !!
    !! Programmer : Bjorn Haugen                          date/rev: Aug 2002/1.0
    !!==========================================================================

    use SamModule             , only : SamType, dp
    use SysMatrixTypeModule   , only : SysMatrixType
    use TireTypeModule        , only : MF_p, JD_TIRE_p
    use TriadTypeModule       , only : transGlobToSys
    use JD_tireRoutinesModule , only : jdt_stiffness
    use AsmExtensionModule    , only : csAddEM
    use rotationModule        , only : EccExpand
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgTire
#endif
    use reportErrorModule     , only : reportError, debugFileOnly_p

    logical            , intent(in)    :: lDynamics
    real(dp)           , intent(in)    :: scaleK
    type(SysMatrixType), intent(inout) :: Nmat
    type(TireType)     , intent(in)    :: tire
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    real(dp) :: s(3), eM(6,6)

    !! --- Logic section ---

    ierr = 0
    if (.not. tire%isInContact) return ! No stiffness when tire is jumping

    !! We've got road contact, or are in static equilibium iterations
    if (tire%type == MF_p .and. tire%rStiff(1) <= 0.0_dp) then
       !! Use STI model file value
       s(3) = tire%sti%typarr(15)*(scaleTo%m/scaleTo%N)*scaleK
    else if (tire%rStiff(1) > 0.0_dp) then
       !! Use value from solver input file
       s(3) = (tire%rStiff(1) + tire%rStiff(2)*tire%defl)*scaleK
    else
       return
    end if

    if (tire%isSliding) then
       s(1:2) = 0.0_dp ! no tangential stiffness
    else if (tire%type == JD_TIRE_p .and. lDynamics) then
       s(1) = 0.0_dp             ! in the dynamics simulation, use
       s(2) = tire%yStiff*scaleK ! stiffness only in Y-direction for JD-tire
    else
       s(1:2) = s(3) ! use radial stiffness in both tangential directions
    end if

    if (tire%type == JD_TIRE_p) then
       call jdt_stiffness (lDynamics,tire,s,scaleK,eM)
    else
       call EccExpand (tire%roadContactInG-tire%RimTriad%ur(:,4),s,eM)
    end if

#ifdef FT_DEBUG
    write(dbgTire,600) tire%id%userId, eM(3,3), scaleK
600 format('addInTireStiffMat:  iTire=',i2,' eM(3,3)=',1pe12.3,' scaleK=',e12.3)
#endif

    call transGlobToSys (tire%RimTriad,eM,1,tire%RimTriad%nDOFs)
    call csAddEM (sam,tire%samElnum,tire%RimTriad%nDOFs,eM,Nmat,ierr,Rhs)
    if (ierr < 0) call reportError (debugFileOnly_p,'addInTireStiffMat')

  end subroutine addInTireStiffMat


  subroutine addInTireDamperMat (scaleC,Nmat,tire,sam,ierr,Rhs)

    !!==========================================================================
    !! Calculates the damping matrix for a tire and adds to
    !! the system Newton matrix, multiplied by a factor, scaleC.
    !!
    !! Programmer : Bjorn Haugen                          date/rev: Aug 2002/1.0
    !!==========================================================================

    use SamModule          , only : SamType, dp
    use SysMatrixTypeModule, only : SysMatrixType
    use TireTypeModule     , only : MF_p, JD_TIRE_p
    use TriadTypeModule    , only : transGlobToSys
    use AsmExtensionModule , only : csAddEM
    use manipmatrixModule  , only : dyadic_product
#ifdef FT_DEBUG
    use dbgUnitsModule     , only : dbgTire
#endif
    use reportErrorModule  , only : reportError, debugFileOnly_p

    real(dp)           , intent(in)    :: scaleC
    type(SysMatrixType), intent(inout) :: Nmat
    type(TireType)     , intent(in)    :: tire
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    real(dp) :: eM(6,6)

    !! --- Logic section ---

    ierr = 0
    if (.not. tire%isInContact) return ! No damping when tire is jumping

    eM = 0.0_dp ! We've got road contact
    if (tire%rDamp <= 0.0_dp) then
       if (tire%type == MF_p) then
          !! Use STI model file value
          eM(3,3) = tire%sti%typarr(16)*(scaleTo%s/scaleTo%kg)*scaleC
       else
          return
       end if
    else if (tire%type == JD_TIRE_p) then
       !! Use values from solver input file
       eM(1:3,1:3) = dyadic_product(tire%WCinG(:,2)) * tire%yDamp(1)*scaleC
       eM(4:6,4:6) = dyadic_product(tire%WCinG(:,1)) * tire%yDamp(2)*scaleC
       eM(3,3) = eM(3,3) + tire%rDamp*scaleC
#ifdef FT_DEBUG
       write(dbgTire,"('Adding damping (scaleC =',1pe12.3'):',3e12.3)") &
            &     scaleC, em(1,1),em(2,2),em(3,3)
#endif
    else
       !! Use value from solver input file
       eM(3,3) = tire%rDamp*scaleC
    end if

#ifdef FT_DEBUG
    write(dbgTire,600) tire%id%userId, eM(3,3), scaleC
600 format('addInTireDamperMat: iTire=',i2,' eM(3,3)=',1pe12.3,' scaleC=',e12.3)
#endif

    call transGlobToSys (tire%RimTriad,eM,1,tire%RimTriad%nDOFs)
    call csAddEM (sam,tire%samElnum,tire%RimTriad%nDOFs,eM,Nmat,ierr,Rhs)
    if (ierr < 0) call reportError (debugFileOnly_p,'addInTireDamperMat')

  end subroutine addInTireDamperMat


  subroutine addInTireMassMat (scaleM,Nmat,tire,sam,ierr,Rhs)

    !!==========================================================================
    !! Calculates the mass matrix for a tire element and adds to
    !! the system Newton matrix, multiplied by a factor, scaleM.
    !!
    !! Programmer : Bjorn Haugen                          date/rev: Aug 2002/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use KindModule         , only : dp, epsDiv0_p
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType
    use TriadTypeModule    , only : transTriadToSys
    use AsmExtensionModule , only : csAddEM
#ifdef FT_DEBUG
    use dbgUnitsModule     , only : dbgTire
#endif
    use reportErrorModule  , only : reportError, debugFileOnly_p

    real(dp)           , intent(in)    :: scaleM
    type(SysMatrixType), intent(inout) :: Nmat
    type(TireType)     , intent(in)    :: tire
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    logical  :: haveMass
    real(dp) :: eM(6,6)

    !! --- Logic section ---

    eM = 0.0_dp
    ierr = 0
    haveMass = .false.
    if (tire%mass > epsDiv0_p) then
       haveMass = .true.
       eM(1,1) = tire%mass*scaleM
       eM(2,2) = eM(1,1)
       eM(3,3) = eM(1,1)
    end if
    if (tire%Iyy > epsDiv0_p) then
       haveMass = .true.
       eM(4,4) = tire%Iyy*scaleM
       eM(5,5) = tire%Iyy*scaleM
    end if
    if (tire%Izz > epsDiv0_p) then
       haveMass = .true.
       eM(6,6) = tire%Izz*scaleM
    end if

    if (.not. haveMass) return

#ifdef FT_DEBUG
    write(dbgTire,600) tire%id%userId, eM(1,1), eM(5,5), eM(6,6), scaleM
600 format('addInTireMassMat: iTire=',i2,' eM=',1p3e12.3,' scaleM=',e12.3)
#endif

    call transTriadToSys (tire%RimTriad,eM,1,tire%RimTriad%nDOFs)
    call csAddEM (sam,tire%samElnum,tire%RimTriad%nDOFs,eM,Nmat,ierr,Rhs)
    if (ierr < 0) call reportError (debugFileOnly_p,'addInTireMassMat')

  end subroutine addInTireMassMat


  subroutine addInTireForces (Q,RF,tire,gravity,sam,ierr)

    !!==========================================================================
    !! Calculates the contributions to the system force vector from a tire
    !! element and adds them into the system vectors Q and RF.
    !!
    !! Programmer : Bjorn Haugen                          date/rev: Aug 2002/1.0
    !!              Knut Morten Okstad                              May 2008/2.0
    !!==========================================================================

    use SamModule         , only : SamType, dp
    use AsmExtensionModule, only : csAddEV
#ifdef FT_DEBUG
    use TireTypeModule    , only : MF_p
    use dbgUnitsModule    , only : dbgTire
#endif
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp)      , intent(inout) :: Q(:), RF(:)
    type(TireType), intent(in)    :: tire
    real(dp)      , intent(in)    :: gravity(3)
    type(SamType) , intent(in)    :: sam
    integer       , intent(inout) :: ierr

    !! Local variables
    integer  :: err
    real(dp) :: eV(6)

    !! --- Logic section ---

    call calcTireForces (tire,gravity,eV)

    call csAddEV (sam,-tire%samElnum,tire%RimTriad%nDOFs,eV,Q,RF,err)
    if (err /= 0) then
       ierr = ierr - 1
       call reportError (debugFileOnly_p,'addInTireForces')
    end if

#ifdef FT_DEBUG
    write(dbgTire,600) tire%id%userId, eV(1:tire%RimTriad%nDOFs)
600 format('addInTireForces: iTire=',i2,' eV=',1p6e12.3)
    if (tire%type == MF_p) then
       write(dbgTire,"(24x,'disp=',1pe12.3)") tire%sti%varinf(10)
    end if
#endif

  end subroutine addInTireForces


  subroutine addInTireInertiaForces (FI,RF,tire,sam,ierr)

    !!==========================================================================
    !! Calculates the contributions to the system inertia force vector from a
    !! tire element and adds them into the system vectors FI and RF.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Oct 2003/1.0
    !!              Knut Morten Okstad                              May 2008/2.0
    !!==========================================================================

    use SamModule         , only : SamType, dp
    use AsmExtensionModule, only : csAddEV
#ifdef FT_DEBUG
    use dbgUnitsModule    , only : dbgTire
#endif
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp)      , intent(inout) :: FI(:), RF(:)
    type(TireType), intent(in)    :: tire
    type(SamType) , intent(in)    :: sam
    integer       , intent(inout) :: ierr

    !! Local variables
    integer  :: err
    real(dp) :: eV(6)

    !! --- Logic section ---

    call calcInertiaForces (tire,eV,err)
    if (err > 0) then
       call csAddEV (sam,tire%samElnum,tire%RimTriad%nDOFs,eV,FI,RF,err)
    end if
    if (err /= 0) then
       ierr = ierr - 1
       call reportError (debugFileOnly_p,'addInTireInertiaForces')
    end if

#ifdef FT_DEBUG
    write(dbgTire,600) tire%id%userId, eV(1:tire%RimTriad%nDOFs)
600 format('addInTireInertiaForces: iTire=',i2,' eV=',1p6e12.3)
#endif

  end subroutine addInTireInertiaForces


  subroutine calcTireForces (tire,gravity,eV)

    !!==========================================================================
    !! Calculates the force vector for a tire element.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: May 2008/1.0
    !!==========================================================================

    use KindModule     , only : dp, epsDiv0_p
    use TriadTypeModule, only : transGlobToSys

    type(TireType), intent(in)  :: tire
    real(dp)      , intent(in)  :: gravity(3)
    real(dp)      , intent(out) :: eV(6)

    !! --- Logic section ---

    eV = tire%forceInG
    if (tire%mass > epsDiv0_p) then
       eV(1:3) = eV(1:3) + gravity*tire%mass
    end if

    call transGlobToSys (tire%RimTriad,eV)

  end subroutine calcTireForces


  subroutine calcInertiaForces (tire,eV,stat)

    !!==========================================================================
    !! Calculates the inertia force vector for a tire element.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: May 2008/1.0
    !!==========================================================================

    use KindModule       , only : dp, epsDiv0_p
    use TriadTypeModule  , only : transTriadToGlob, transGlobToSys
    use manipMatrixModule, only : cross_product

    type(TireType), intent(in)  :: tire
    real(dp)      , intent(out) :: eV(6)
    integer       , intent(out) :: stat

    !! Local variables
    real(dp) :: vel(3), acc(3), II(3,3)

    !! --- Logic section ---

    if (tire%mass > epsDiv0_p) then
       stat = 3
       eV(1:3) = tire%mass * tire%RimTriad%urdd(1:3)
    else
       stat = 0
       eV(1:3) = 0.0_dp
    end if
    if (tire%Iyy > epsDiv0_p .or. tire%Izz > epsDiv0_p) then
       II = 0.0_dp
       II(1,1) = tire%Iyy
       II(2,2) = tire%Iyy
       II(3,3) = tire%Izz
       call transTriadToGlob (tire%RimTriad,II,1,3)
       vel = tire%RimTriad%urd(4:6)
       acc = tire%RimTriad%urdd(4:6)
       eV(4:6) = matmul(II,acc) + cross_product(vel,matmul(II,vel))
       stat = 6
    else
       eV(4:6) = 0.0_dp
    end if

    if (stat > 0) then
       call transGlobToSys (tire%RimTriad,eV(1:stat))
    end if

  end subroutine calcInertiaForces


  subroutine updateTireForces (tires,gravity)

    !!==========================================================================
    !! Updates the triad forces with tire contributions.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: May 2008/1.0
    !!==========================================================================

    use TriadTypeModule, only : updateNodeForce, dp

    type(TireType), intent(inout) :: tires(:)
    real(dp)      , intent(in)    :: gravity(3)

    !! Local variables
    integer  :: i, stat
    real(dp) :: eV(6)

    !! --- Logic section ---

    do i = 1, size(tires)

       call calcTireForces (tires(i), gravity, eV)
       call updateNodeForce (tires(i)%RimTriad, -eV)

       call calcInertiaForces (tires(i), eV, stat)
       call updateNodeForce (tires(i)%RimTriad, eV)

    end do

  end subroutine updateTireForces


  subroutine setTireStateVar (sys,tires)

    !!==========================================================================
    !! Initiates all tire state variables with predictor-step values.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Jan 2003/2.0
    !!==========================================================================

    use SystemTypeModule      , only : SystemType
    use TireTypeModule        , only : STI_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(SystemType), intent(in)    :: sys
    type(TireType)  , intent(inout) :: tires(:)

    !! Local variables
    integer :: i, tireIntegrator

    !! --- Logic section ---

    call ffa_cmdlinearg_getint ('tireIntegrator',tireIntegrator)

    do i = 1, size(tires)
       if (tires(i)%api == STI_p) then
          call initiateSTIstateVar (tireIntegrator, sys%timeStep*scaleTo%s, &
               &                    tires(i)%sti)
       end if
    end do

  end subroutine setTireStateVar


  subroutine updateTires (sys,tires,ierr,updateCTI,staticIterations)

    !!==========================================================================
    !! Updates all tire elements.
    !!
    !! Programmer : Bjorn Haugen                          date/rev: Aug 2002/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use SystemTypeModule      , only : SystemType
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgSolve, dbgTire
#endif
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(SystemType), intent(in)    :: sys
    type(TireType)  , intent(inout) :: tires(:)
    integer         , intent(inout) :: ierr
    integer,optional, intent(in)    :: updateCTI
    logical,optional, intent(in)    :: staticIterations

    !! Local variables
    integer :: i, lerr, tireStaticEquil, tireIntegrator, okCTI

    !! --- Logic section ---

    if (present(updateCTI)) then
       okCTI = updateCTI
    else
       okCTI = -1
    end if

    if (present(staticIterations)) then ! Invoked from staticIncAndUpdate
       if (staticIterations) then
          !! We are within the initial static equilibrium iterations
          call ffa_cmdlinearg_getint ('tireStaticEquil',tireStaticEquil)
          if (tireStaticEquil > 0) then
             okCTI = -(1+tireStaticEquil)
          else
             return ! Ignore tires in the static equilibrium iterations
          end if
       else
          return ! No static equilibrium iterations requested
       end if
    end if

    call ffa_cmdlinearg_getint ('tireIntegrator',tireIntegrator)

#ifdef FT_DEBUG
    if (size(tires) > 0) write(dbgSolve,"(10x,'UpdateTires:',i3)") okCTI
    if (size(tires) > 0) write(dbgTire,600) sys%time,sys%nIterThisStep,okCTI
600 format(/'updateTires: time=',1pe12.5,' iter=',i3,' okCTI=',i3)
#endif

    lerr = ierr
    do i = 1, size(tires)
       call updateTire (tireIntegrator,okCTI,sys,tires(i),ierr)
    end do

    if (ierr < lerr) call reportError (debugFileOnly_p,'updateTires')

  end subroutine updateTires


  subroutine updateCTIresAtConvergence (sys,tires,ierr)

    !!==========================================================================
    !! Updates all CTI tire elements after convergence has been achieved.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Jan 2003/1.0
    !!==========================================================================

    use SystemTypeModule, only : SystemType
    use TireTypeModule  , only : CTI_p
#ifdef FT_DEBUG
    use dbgUnitsModule  , only : dbgTire
#endif

    type(SystemType), intent(in)    :: sys
    type(TireType)  , intent(inout) :: tires(:)
    integer         , intent(inout) :: ierr

    !! Local variables
    integer :: i

    integer, parameter :: okCTI_p = 1 ! System state has been accepted

    !! --- Logic section ---

#ifdef FT_DEBUG
    if (size(tires) > 0) write(dbgTire,600) sys%time,sys%nIterThisStep
600 format(/'updateCTIresAtConvergence: time=',1pe12.5,' iter=',i3)
#endif

    do i = 1, size(tires)
       if (tires(i)%api == CTI_p) call updateTire (0,okCTI_p,sys,tires(i),ierr)
    end do

  end subroutine updateCTIresAtConvergence


  subroutine closeTireModules (tires)

    !!==========================================================================
    !! Close all tire modules upon program completion.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Feb 2003/1.0
    !!==========================================================================

    use TireTypeModule  , only : CTI_p

    type(TireType)  , intent(inout) :: tires(:)

    !! Local variables
    integer :: i, nCTI

    !! --- Logic section ---

    nCTI = 0
    do i = 1, size(tires)
       if (tires(i)%api == CTI_p) nCTI = nCTI + 1
    end do

#ifdef FT_HAS_CTI
    if (nCTI > 0) call CTICLS ! Close the CTI interface
#endif

  end subroutine closeTireModules


  subroutine initiateTire (tire,time,ierr)

    !!==========================================================================
    !! Initiates model data for a tire element.
    !!
    !! Programmer : Bjorn Haugen                          date/rev: Aug 2002/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use KindModule         , only : dp, epsDiv0_p
    use IdTypeModule       , only : getId
    use TireTypeModule     , only : allocateSTIapi, STI_p, CTI_p
    use manipMatrixModule  , only : matmul34, invert34, cross_product
#ifdef FT_DEBUG
    use fileUtilitiesModule, only : getDBGfile
    use dbgUnitsModule     , only : dbgTire
#endif
    use reportErrorModule  , only : allocationError, internalError
    use reportErrorModule  , only : reportError, error_p

    type(TireType), intent(inout) :: tire
    real(dp)      , intent(in)    :: time
    integer       , intent(out)   :: ierr

    !! Local variables
    integer      :: jobflg
    type(WCdata) :: wc

    !! --- Logic section ---

#ifdef FT_DEBUG
    if (dbgTire == 0) dbgTire = getDBGfile(88,'tire.dbg')
#endif

    select case (tire%api)

    case (STI_p)

       !! Calculate the Wheel Carrier (WC) system in global coordinates
       tire%WCinG = getWCsystem(tire%WCtriad%ur,tire%ZoffSet,tire%WCYalongZ)
       if (ierr /= 0) return

       !! Get WC-system relative to the WC triad system
       tire%WCinT = matmul34(invert34(tire%WCtriad%ur),tire%WCinG)

       wc%time  = time
       wc%pos   = tire%WCinG(:,4)*scaleTo%m
       wc%trans = tire%WCinG(:,1:3)

       !! jobflag = 1 : get necessary array sizes and allocate arrays
       !! jobflag = 2 : initialize tire model

       do jobflg = 1, 2

          call STI (jobflg, .false., wc, &
               &    tire%WCtriad%urd(1:3)*scaleTo%m/scaleTo%s, &
               &    tire%WCtriad%urd(4:6)/scaleTo%s, tire, ierr)
          if (ierr /= 0) then
             call reportError (error_p, &
                  'Unable to initialize Tire'//getId(tire%id), &
                  addString='initiateTire')
             exit
          end if

          if (jobflg == 1) then
             call allocateSTIapi (tire%sti,ierr)
             if (ierr /= 0) then
                ierr = allocationError('initiateTire')
                return
             end if
          end if

       end do

       tire%radius = tire%sti%typarr(7) / scaleTo%m ! Unloaded tire radius

    case (CTI_p)

       !! Calculate the rim system in global coordinates
       tire%WCinG = getWCsystem(tire%RimTriad%ur,tire%ZoffSet,tire%WCYalongZ)
       if (ierr /= 0) return

       !! Get rim system relative to the rim triad system
       tire%WCinT = matmul34(invert34(tire%RimTriad%ur),tire%WCinG)

#ifdef FT_HAS_CTI
#ifdef CTI_DEBUG
       call CTIV (tire%idIn,1) ! CTI verbose mode
#endif

       !! Hand over the tire data file to CTI
       call CTILTF (tire%idIn,ierr,trim(tire%tireDataFileName))
       if (ierr /= 0) then
          call reportError (error_p, &
               'Unable to initialize Tire'//getId(tire%id), &
               'Failed to load tire data file '//tire%tireDataFileName, &
               addString='initiateTire')
          return
       end if

       !! Hand over the road data file to CTI
       call CTILRF (tire%idIn,ierr,trim(tire%road%roadDataFileName))
       if (ierr /= 0) then
          call reportError (error_p, &
               'Unable to initialize Tire'//getId(tire%id), &
               'Failed to load road data file '//tire%road%roadDataFileName, &
               addString='initiateTire')
          return
       end if
#else
       ierr = -2
       call reportError (error_p,'The CTI tire interface is not available '// &
            'in this version.',addString='initiateTire')
#endif

    case default

       ierr = -1
       call reportError (error_p,'Invalid tire API for Tire'//getId(tire%id), &
            'Check that consistent tire type and API is used.', &
            addString='initiateTire')

    end select

  contains

    function getWCsystem (WCtriad,ZoffSet,WCYalongZ) result(T)
      real(dp), intent(in) :: WCtriad(:,:), ZoffSet
      logical , intent(in) :: WCYalongZ
      real(dp)             :: T(3,4), lenSqr

      !! The WC position equals position of WC-triad
      !! + offset along the local Z-axis of the WC triad
      T(:,4) = WCtriad(:,4) + ZoffSet*WCtriad(:,3)

      !! Y-axis coincides with the Z-direction of the WC-triad
      !! since we currently are using a revolute joint triad here
      if (WCYalongZ) then
         T(:,2) =  WCtriad(:,3)
      else
         T(:,2) = -WCtriad(:,3)
      end if

      !! X-axis as cross product of the Y-axis and the global Z-axis = <0,0,1>
      T(1,1) =  T(2,2)
      T(2,1) = -T(1,2)
      T(3,1) =  0.0_dp
      lenSqr =  T(1,1)*T(1,1) + T(2,1)*T(2,1)
      if (lenSqr < epsDiv0_p) then
         ierr = 1
         call reportError (error_p, &
              'Unable to initialize Tire'//getId(tire%id), &
              'The wheel rotation axis is parallel to the global Z-axis.', &
              addString='initiateTire')
      else
         ierr = 0
         T(:,1) = T(:,1)/sqrt(lenSqr)
         !! Z-axis as cross product of the X- and Y-axes
         T(:,3) = cross_product(T(:,1),T(:,2))
      end if

    end function getWCsystem

  end subroutine initiateTire


  subroutine updateTire (tireIntegrator,okCTI,sys,tire,err)

    !!==========================================================================
    !! Updates a tire element.
    !!
    !! Note: The interpretation of the 'okCTI' argument depends on the API used:
    !!       for STI: = -3: Static equilibrium iterations, initial velocities.
    !!                = -2: Static equilibrium iterations with zero velocities.
    !!                > -2: Dynamic equilibrium iterations.
    !!       for CTI: <  0: Do nothing.
    !!                >= 0: As defined by the CTI API (see COSIN documentation).
    !!
    !! Programmer : Bjorn Haugen                          date/rev: Aug 2002/1.0
    !!              Knut Morten Okstad                              Jan 2003/2.0
    !!==========================================================================

    use KindModule       , only : dp, epsDiv0_p
    use SystemTypeModule , only : SystemType
    use IdTypeModule     , only : getId
    use TireTypeModule   , only : STI_p, CTI_p
    use TriadTypeModule  , only : transVGlobToTriad
    use manipMatrixModule, only : matmul34, cross_product
    use profilerModule   , only : startTimer, stopTimer, tir_p
#ifdef FT_DEBUG
    use dbgUnitsModule   , only : dbgTire
#endif
    use reportErrorModule, only : internalError, reportError
    use reportErrorModule, only : empty_p, error_p, warning_p, warningFileOnly_p

    integer         , intent(in)    :: tireIntegrator, okCTI
    type(SystemType), intent(in)    :: sys
    type(TireType)  , intent(inout) :: tire
    integer         , intent(inout) :: err

    !! Local variables
    integer  :: ierr, jerr
    real(dp) :: velInG(3), omegaInG(3), eccVecInG(3)
#ifdef FT_HAS_CTI
    real(dp) :: rmax
#endif
    type(WCdata) :: wc

    integer, parameter :: jobflg_p = 0 ! Dynamic analysis (normal mode)

    character(len=32) :: errMsg
    integer, save     :: nWarn = 0

    !! --- Logic section ---

    ierr = 0
    if (okCTI <= -2) then ! Static equilibrium iterations

       !! Set contact to true and sliding to false so that stiffness is always
       !! added to the stiffness matrix, even when not in contact. This is done
       !! in order to avoid "Zero pivot" before contact actually is obtained.
       tire%isInContact = .true.
       tire%isSliding   = .false.

    end if

    select case (tire%api)

    case (STI_p)

       !! Update the WC position
       tire%WCinG = matmul34(tire%WCTriad%ur,tire%WCinT)

       wc%time  = sys%time * scaleTo%s
       wc%pos   = tire%WCinG(:,4) * scaleTo%m
       wc%trans = tire%WCinG(:,1:3)

       !! Velocity of WC
       velInG   = tire%WCTriad%urd(1:3) * scaleTo%m/scaleTo%s
       omegaInG = tire%WCTriad%urd(4:6) / scaleTo%s
       if (tire%hasEccVec) then
          !! Add eccentric effect to velocity
          eccVecInG = matmul(tire%WCTriad%ur(:,1:3),tire%WCinT(:,4)) * scaleTo%m
          velInG  = velInG + cross_product(omegaInG,eccVecInG)
       end if

       !! Invoke tire model through the STI interface
       call startTimer (tir_p)
       call STI (jobflg_p, okCTI == -2, wc, velInG, omegaInG, tire, ierr)
       call stopTimer (tir_p)

       if (ierr == 1) then ! Warning from DTYRE, try continue simulation

          ierr = 0
          if (nWarn < 10) then
             jerr = warning_p
          else
             jerr = warningFileOnly_p
          end if
          nWarn = nWarn + 1
          write(errMsg,"(', time =',1pe12.5,', iter =',i3)") &
               &       sys%time, sys%nIterThisStep
          call reportError (jerr,'Problems updating Tire'// &
               &            trim(getId(tire%id))//errMsg)
          if (nWarn == 10) then
             call reportError (empty_p,'Last message of this kind. '// &
                  'More such messages will be printed to the res-file only.')
          end if

       end if

       !! Integrate the tire state variables
       call updateSTIstateVar (tireIntegrator, sys%timeStep*scaleTo%s, tire%sti)

#ifdef FT_DEBUG
       if (size(tire%sti%deqvar) > 0) then
          write(dbgTire,"('iTire=',i2,' degvar=',1pe12.3,' deqder=',e12.3)") &
               &        tire%id%userId, tire%sti%deqvar(1), tire%sti%deqder(1)
          do jerr = 2, size(tire%sti%deqvar)
             write(dbgTire,"(1pe28.3,e20.3)") &
                  &        tire%sti%deqvar(jerr), tire%sti%deqder(jerr)
          end do
       end if
#endif

       tire%defl           = tire%sti%varinf(10) / scaleTo%m
       tire%radius         = tire%sti%varinf(12) / scaleTo%m
       tire%roadContactInG = tire%sti%varinf(13:15) / scaleTo%m

       if (okCTI > -2) then
          !! Check if tire is in road contact, and if it is sliding or not
          tire%isInContact = tire%defl > epsDiv0_p
          if (tire%isInContact) then
             tire%isSliding = abs(tire%sti%varinf(2)) > 0.5_dp
          end if
       end if

    case (CTI_p)
#ifdef FT_HAS_CTI
       if (okCTI < 0) return ! No action (yet) for static equilibrium iterations

       !! Update the rim position
       tire%WCinG = matmul34(tire%RimTriad%ur,tire%WCinT)

       !! Position of rim
       posInG = tire%WCinG(:,4) * scaleTo%m

       !! Rim velocity
       velInG   = tire%RimTriad%urd(1:3) * scaleTo%m/scaleTo%s
       omegaInG = tire%RimTriad%urd(4:6) / scaleTo%s
       if (tire%hasEccVec) then
          eccVecInG = 0.0_dp !TODO,kmo: Add proper eccentric effect to velocity
       end if

#ifdef FT_DEBUG
       write(dbgTire,"('updateTire CTI:',2i3,1pe13.5)") tire%idIn,okCTI,sys%time
       write(dbgTire,"('      position:',1p6e13.5)") posInG,tire%WCinG(1,1:3)
       write(dbgTire,"('               ',39x,1p3e13.5)")    tire%WCinG(2,1:3)
       write(dbgTire,"('               ',39x,1p3e13.5)")    tire%WCinG(3,1:3)
       write(dbgTire,"('      velocity:',1p6e13.5)") velInG,omegaInG
#endif

       !! Invoke tire model through the CTI interface
       call startTimer (tir_p)
       call CTI (tire%idIn, sys%time*scaleTo%s, &
            &    posInG, tire%WCinG, velInG, omegaInG, okCTI, &
            &    tire%forceInG(1:3), tire%forceInG(4:6), ierr)
       call stopTimer (tir_p)

#ifdef FT_DEBUG
       write(dbgTire,"(9x,'force:',1p6e13.5)") tire%forceInG
#endif

       !! Calculate forces acting on the rim triad in WC coordinates
       tire%forceInWC = transVGlobToTriad(tire%WCtriad,tire%forceInG)

       tire%isInContact = .false. !TODO,kmo: Need to set this based on tire defl

       !! Get some tire properties
       call CTIPTP (tire%idIn,rmax,tire%radius,tire%mass,tire%Izz,tire%Iyy, &
            &       tire%rStiff(1),tire%rStiff(2),jerr)
       if (jerr /= 0) then
          if (ierr == 0) ierr = jerr
          call reportError (error_p,'Tire data can not be used', &
               &            addString='CTIPTP')
       else
          !! Scale from SI to the modelling units used
          tire%radius = tire%radius / scaleTo%m
          tire%mass   = tire%mass   / scaleTo%kg
          tire%Iyy    = tire%Iyy    /(scaleTo%kg*scaleTo%m*scaleTo%m)
          tire%Izz    = tire%Izz    /(scaleTo%kg*scaleTo%m*scaleTo%m)
          tire%rStiff = tire%rStiff * scaleTo%m/scaleTo%N
       end if
#ifdef FT_DEBUG
       write(dbgTire,"('initial radius:',1p2e13.5)") rmax
       write(dbgTire,"('dynamic radius:',1p2e13.5)") tire%radius
       write(dbgTire,"('rim-fixed mass:',1p2e13.5)") tire%mass
       write(dbgTire,"('--||-- inertia:',1p2e13.5)") tire%Iyy, tire%Izz
       write(dbgTire,"('     stiffness:',1p2e13.5)") tire%rStiff
#endif
#else
       err = internalError('updateTire: The CTI interface is not available')
       return
#endif

    case default

       err = internalError('updateTire: Invalid tire API')
       return

    end select

    if (ierr /= 0) then ! Error occurred, abort simulation

       err = err - abs(ierr)
       write(errMsg,"(', time =',1pe12.5,', ierr =',i3)") sys%time,ierr
       call reportError (error_p,'Unable to update Tire'// &
            &            trim(getId(tire%id))//errMsg, &
            &            addString='updateTire')
       return

    end if

    if (tire%hasEccVec) then
       !! Add eccentric effect
       tire%forceInG(4:6) = tire%forceInG(4:6) &
            &             + cross_product(eccVecInG,tire%forceInG(1:3))
    end if

    !! Scale forces back to modeling units
    tire%forceInG (1:3) = tire%forceInG (1:3) / scaleTo%N
    tire%forceInWC(1:3) = tire%forceInWC(1:3) / scaleTo%N

    !! Scale torques back to modeling units
    tire%forceInG (4:6) = tire%forceInG (4:6) / (scaleTo%N*scaleTo%m)
    tire%forceInWC(4:6) = tire%forceInWC(4:6) / (scaleTo%N*scaleTo%m)

  end subroutine updateTire


  subroutine STI (jobflg,zeroVel,WC,velWC,omegaWC,tire,ierr)

    !!==========================================================================
    !! F90-incapsulation of the DTYRE subroutine of the STI interface.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Jan 2003/1.0
    !!==========================================================================

    use KindModule       , only : dp
#ifdef FT_HAS_TNO_TYRE
    use TireTypeModule   , only : MF_p, SWIFT_p
#endif
#ifdef FT_HAS_JDTYRE
    use TireTypeModule   , only : JD_TIRE_p
#endif
    use IdTypeModule     , only : getId
    use reportErrorModule, only : getErrorFile, reportError, error_p

    integer       , intent(in)    :: jobflg
    logical       , intent(in)    :: zeroVel
    type(WCdata)  , intent(in)    :: WC
    real(dp)      , intent(in)    :: velWC(:), omegaWC(:)
    type(TireType), intent(inout) :: tire
    integer       , intent(out)   :: ierr

    !! Local variables
    real(dp) :: velWCinWC(3), omegaWCinWC(3), alphaRimToWC, omegaRim

    external FedemRoadMain

    !! --- Logic section ---

    !! Transform the WC velocity from global to WC coordinate system
    velWCinWC    = matmul(velWC,WC%trans)
    omegaWCinWC  = matmul(omegaWC,WC%trans)

    !! Rotation angle and velocity of the wheel (using revolute joint variables)
    alphaRimToWC = tire%pJvar(1)
    omegaRim     = tire%pJvar(2) / scaleTo%s

    if (zeroVel) then
       !bh: Used by static equilibrium algorithm. Necessary when initial
       !bh: velocities are non-zero. Should ideally be able to use the
       !bh: initial velocities here as well, but this needs some more work.
       velWCinWC   = 0.0_dp
       omegaWCinWC = 0.0_dp
       omegaRim    = 0.0_dp
    end if

    !! Store current tire radius in ropar(1)
    if (tire%road%nropar > 0) tire%road%roparr(1) = tire%radius

    select case (tire%type)
#ifdef FT_HAS_TNO_TYRE
    case (MF_p, SWIFT_p); call DTYRE (getErrorFile(), &
         tire%sti%iswtch, &
         jobflg, &
         tire%idIn, &
         WC%time, &
         WC%pos(1), &           ! input, global WC coordinates
         WC%trans(1,1), &       ! input, transformation from WC-system to global
         alphaRimToWC, &        ! input, angle between rim and WC
         velWCinWC(1), &        ! input, velocity of WC in WC-system
         omegaWCinWC(1), &      ! input, angular velocity of WC in WC-system
         omegaRim, &            ! input, angular velocity of rim relative to WC
         tire%sti%ndeqvar, &    ! input, number of state variables
         tire%sti%deqvar(1), &  ! input, value of state variables
         tire%sti%nTyparr, &    ! input, number of tyre parameters
         tire%sti%typarr(1), &  ! input, array of tyre parameters
         tire%nCharTireDataFileName, &
         tire%tireDataFileName, &
         FedemRoadMain, &       ! subroutine for road contact
         tire%road%idIn, &      ! input, road identifier
         tire%road%nRopar, &    ! input, number of road parameters
         tire%road%roparr(1), & ! input, array of road parameters
         tire%road%nCharRoadDataFileName, &
         tire%road%roadDataFileName, &
         tire%forceInWC(1), &   ! output, force applied to wheel center
         tire%forceInWC(4), &   ! output, torque applied to wheel center
         tire%sti%deqini(1), &  ! input, initial value for state variables
         tire%sti%deqder(1), &  ! output, gradients of state variables
         tire%sti%tyrmod, &     ! output, tyre model description
         tire%sti%nvars, &      ! input, dimension of varinf
         tire%sti%varinf(1), &  ! output, see file "tno_varptr.inc"
         tire%sti%nwork, &
         tire%sti%wrkarr(1), &
         tire%sti%niwork, &
         tire%sti%iwrkarr(1), &
         ierr)
#endif
#ifdef FT_HAS_JDTYRE
    case (JD_TIRE_p); call JDTYRE (getErrorFile(), &
         tire%sti%iswtch, &
         jobflg, &
         tire%idIn, &
         WC%time, &
         WC%pos(1), &           ! input, global WC coordinates
         WC%trans(1,1), &       ! input, transformation from WC-system to global
         alphaRimToWC, &        ! input, angle between rim and WC
         velWCinWC(1), &        ! input, velocity of WC in WC-system
         omegaWCinWC(1), &      ! input, angular velocity of WC in WC-system
         omegaRim, &            ! input, angular velocity of rim relative to WC
         tire%sti%ndeqvar, &    ! input, number of state variables
         tire%sti%deqvar(1), &  ! input, value of state variables
         tire%sti%nTyparr, &    ! input, number of tyre parameters
         tire%sti%typarr(1), &  ! input, array of tyre parameters
         tire%nCharTireDataFileName, &
         tire%tireDataFileName, &
         FedemRoadMain, &       ! subroutine for road contact
         tire%road%idIn, &      ! input, road identifier
         tire%road%nRopar, &    ! input, number of road parameters
         tire%road%roparr(1), & ! input, array of road parameters
         tire%road%nCharRoadDataFileName, &
         tire%road%roadDataFileName, &
         tire%forceInWC(1), &   ! output, force applied to wheel center
         tire%forceInWC(4), &   ! output, torque applied to wheel center
         tire%sti%deqini(1), &  ! input, initial value for state variables
         tire%sti%deqder(1), &  ! output, gradients of state variables
         tire%sti%tyrmod, &     ! output, tyre model description
         tire%sti%nvars, &      ! input, dimension of varinf
         tire%sti%varinf(1), &  ! output, see file "tno_varptr.inc"
         tire%rStiff(1), &
         tire%rDamp, &
         tire%yStiff, &
         tire%yDamp(1), &
         ierr)
#endif

    case default

       ierr = -1
       call reportError (error_p,'Invalid tire type for Tire'//getId(tire%id), &
            &            'Check if type is consistent with the STI tyre API.', &
            &            addString='STI')

    end select

    if (ierr > 1 .or. jobflg > 0) return

    !! Calculate forces acting on the rim triad in global coordinates
    tire%forceInG(1:3) = matmul(WC%trans,tire%forceInWC(1:3))
    tire%forceInG(4:6) = matmul(WC%trans,tire%forceInWC(4:6))

  end subroutine STI


  subroutine initiateSTIstateVar (tireIntegrator,dt,sti)

    !!==========================================================================
    !! Initiate STI tire state variables with predictor-step values.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Jan 2003/1.0
    !!==========================================================================

    use TireTypeModule, only : STIapiType, dp

    integer         , intent(in)    :: tireIntegrator
    real(dp)        , intent(in)    :: dt
    type(STIapiType), intent(inout) :: sti

    !! --- Logic section ---

    sti%deqvar_prev = sti%deqvar

    select case (tireIntegrator)

    case (0:1) ! Forward/Backward Euler
       sti%deqvar      = sti%deqvar_prev + sti%deqder*dt
       sti%deqder_prev = sti%deqder

    case (2)   ! Improved Euler (Heuns method)
       sti%deqvar      = sti%deqvar_prev &
            &          + (sti%deqder + sti%deqder_prev) * 0.5_dp*dt
       sti%deqder_prev = sti%deqder

    case (3)   ! Basic theta-method (theta=0.505)
       sti%deqvar      = sti%deqvar_prev &
            &          + (0.505_dp*sti%deqder_prev + 0.495_dp*sti%deqder) * dt
       sti%deqder_prev = sti%deqder

    case (4)   ! Adams-Moulton 2
       sti%deqvar      = sti%deqvar_prev &
            &          + (5.0_dp*sti%deqder + sti%deqder_prev) * dt/12.0_dp
       sti%deqder_prev =  8.0_dp*sti%deqder - sti%deqder_prev

    end select

  end subroutine initiateSTIstateVar


  subroutine updateSTIstateVar (tireIntegrator,dt,sti)

    !!==========================================================================
    !! Integrate the STI tire state variables in time.
    !!
    !! Programmer : Knut Morten Okstad                    date/rev: Jan 2003/1.0
    !!==========================================================================

    use TireTypeModule, only : STIapiType, dp

    integer         , intent(in)    :: tireIntegrator
    real(dp)        , intent(in)    :: dt
    type(STIapiType), intent(inout) :: sti

    !! --- Logic section ---

    select case (tireIntegrator)

    case (1) ! Backward Euler
       sti%deqvar = sti%deqvar_prev + sti%deqder*dt

    case (2) ! Improved Euler (Heuns method)
       sti%deqvar = sti%deqvar_prev + (sti%deqder_prev + sti%deqder) * 0.5_dp*dt

    case (3) ! Basic theta-method (theta=0.505)
       sti%deqvar = sti%deqvar_prev &
            &     + (0.505_dp*sti%deqder_prev + 0.495_dp*sti%deqder) * dt

    case (4) ! Adams-Moulton 2
       sti%deqvar = sti%deqvar_prev &
            &     + (sti%deqder_prev + 5.0_dp*sti%deqder) * dt/12.0_dp

    end select

  end subroutine updateSTIstateVar

end module TireRoutinesModule
