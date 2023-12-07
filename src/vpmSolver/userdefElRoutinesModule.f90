!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file userdefElRoutinesModule.f90
!> @brief Subroutines for user-defined element calculations.

!!==============================================================================
!> @brief Module with subroutines for user-defined element calculations.
!>
!> @details This module contains a set of subroutines for performing various
!> computation tasks on the userdefeltypemodule::userdefeltype objects
!> in the model during the dynamic or quasi-static simulation.

module UserdefElRoutinesModule

  implicit none

  private :: updateUDE, calcBuoyancyForces, calcMorisonForces


contains

  !!============================================================================
  !> @brief Initializes a user-defined element.
  !>
  !> @param[in] elmId Id of the user-defined element to initialize
  !> @param[in] etype User-defined element type flag
  !> @param[in] nedof Dimension of the element matrices
  !> @param[in] triads The triads connected to this user-defined element
  !> @param iwork Integer work array for this user-defined element
  !> @param rwork Real work array for this user-defined element
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine InitializeUDE (elmId,etype,nedof,triads,iwork,rwork,ierr)

    use IdTypeModule      , only : IdType, getId
    use TriadTypeModule   , only : TriadPtrType, dp
    use FiUserElmInterface, only : Fi_UDE0, Fi_UDE1
    use allocationModule  , only : reAllocate
    use scratchArrayModule, only : getRealScratchArray
    use reportErrorModule , only : reportError, error_p, debugFileOnly_p

    type(IdType)      , intent(in)  :: elmId
    integer           , intent(in)  :: etype, nedof
    type(TriadPtrType), intent(in)  :: triads(:)
    integer           , pointer     :: iwork(:)
    real(dp)          , pointer     :: rwork(:)
    integer           , intent(out) :: ierr

    !! Local variables
    integer           :: i, ipT, nenod, niwork, nrwork, nwork
    real(dp), pointer :: nodeVar(:)

    !! --- Logic section ---

    nenod = size(triads)
    nwork = 3*nenod
    do i = 1, nenod
       if (triads(i)%p%nDOFs >= 6) nwork = nwork + 9
    end do
    if (nwork == 3*nenod) nwork = nwork + 1
    nodeVar => getRealScratchArray(nwork,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'InitializeUDE')
       return
    end if

    !! Extract arrays of nodal variables from the triads

    ipT = 3*nenod + 1
    do i = 1, nenod
       nodeVar(3*i-2:3*i) = triads(i)%p%ur(:,4)
       if (triads(i)%p%nDOFs >= 6) then
          nodeVar(ipT:ipT+8) = reshape(triads(i)%p%ur(:,1:3),(/9/))
          ipT = ipT + 9
       end if
    end do

    !! Re-allocate the work arrays to fit the requested size from the plugin

    call Fi_UDE0 (elmId%baseId,etype,nenod,nedof,niwork,nrwork)
    call reAllocate ('InitializeUDE',iwork,niwork,ierr,.true.)
    call reAllocate ('InitializeUDE',rwork,nrwork,ierr,.true.)
    if (ierr < 0) return

    ierr = min(niwork-3,nrwork-4-iwork(3))
    if (ierr < 0) then
       call reportError (error_p,'Invalid element plugin.', &
            &            'The work array dimensions are too small.', &
            &            addString='InitializeUDE')
       return
    end if

    iwork(1) = niwork
    iwork(2) = nrwork

    !! Initialize the user-defined element

    ipT = 3*nenod + 1
    call Fi_UDE1 (elmId%baseId,etype,nenod,nedof,nodeVar(1),nodeVar(ipT), &
         &        iwork(1),rwork(1),ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failed to initialize User-defined Element'// &
            &            getId(elmId),addString='InitializeUDE')
    end if

  end subroutine InitializeUDE


  !!============================================================================
  !> @brief Updates a user-defined element.
  !>
  !> @param[in] elmId Id of the user-defined element to update
  !> @param[in] etype User-defined element type flag
  !> @param[in] nedof Dimension of the element matrices
  !> @param[in] triads The triads connected to this user-defined element
  !> @param[in] engines The engines connected to this user-defined element
  !> @param iwork Integer work array for this user-defined element
  !> @param rwork Real work array for this user-defined element
  !> @param hydyn Data for hydrodynamic force calculation for this element
  !> @param envir Environmental data
  !> @param[out] K Tangent stiffness matrix
  !> @param[out] C Damping matrix
  !> @param[out] M Mass matrix
  !> @param[out] Fs Internal elastic forces
  !> @param[out] Fd Damping forces
  !> @param[out] Fi Intertia forces
  !> @param[out] Q External forces
  !> @param[out] Tlg Local-to-global transformation matrix for the element
  !> @param[in] t Current simulation time
  !> @param[in] dt Time increment size
  !> @param[in] istep Time increment counter
  !> @param[in] iter Iteration counter
  !> @param[out] ierr Error flag
  !>
  !> @details This is the main work-horse for the user-defined element objects.
  !> It is responsible for calculation of all tangent matrices and associated
  !> right-hand-side force vectors eminating from the linearized dynamic
  !> equilibrium equation. The subroutine is typically invoked once per
  !> Newton iteration for each element.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine UpdateUDE (elmId,etype,nedof,triads,engines, &
       &                iwork,rwork,hydyn,envir, &
       &                K,C,M,Fs,Fd,Fi,Q,Tlg,t,dt,istep,iter,ierr)

    use IdTypeModule         , only : IdType, getId
    use TriadTypeModule      , only : TriadPtrType, dp
    use FunctionTypeModule   , only : EnginePtrType
    use SupElTypeModule      , only : HydroDynType
    use EnvironmentTypeModule, only : EnvironmentType
    use EngineRoutinesModule , only : EngineValue
    use FiUserElmInterface   , only : Fi_UDE2, Fi_UDE3
    use scratchArrayModule   , only : getRealScratchArray
#ifdef FT_DEBUG
    use manipMatrixModule, only : writeObject
    use dbgUnitsModule   , only : dbgUDE
#endif
    use reportErrorModule, only : reportError, error_p, debugFileOnly_p

    type(IdType)         , intent(in)    :: elmId
    integer              , intent(in)    :: etype, nedof
    type(TriadPtrType)   , intent(in)    :: triads(:)
    type(EnginePtrType)  , intent(in)    :: engines(:)
    integer              , intent(inout) :: iwork(:)
    real(dp)             , intent(inout) :: rwork(:)
    type(HydroDynType)   , pointer       :: hydyn
    type(EnvironmentType), intent(inout) :: envir
    real(dp)             , intent(out)   :: K(:,:), C(:,:), M(:,:)
    real(dp)             , intent(out)   :: Fs(:), Fd(:), Fi(:), Q(:), Tlg(3,4)
    real(dp)             , intent(in)    :: t, dt
    integer              , intent(in)    :: istep, iter
    integer              , intent(out)   :: ierr

    !! Local variables
    integer           :: i, ipT, ipvel, ipacc, itmp, nenod, nndof, nwork
    real(dp), pointer :: nodeVar(:)

    !! --- Logic section ---

    nenod = size(triads)
    ipT   = 3*nenod + 1
    ipvel = ipT
    do i = 1, nenod
       if (triads(i)%p%nDOFs >= 6) ipvel = ipvel + 9
    end do
    ipacc = ipvel + nedof
    nwork = ipacc + nedof-1
    nodeVar => getRealScratchArray(nwork,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'UpdateUDE')
       return
    end if

    !! Extract arrays of nodal variables from the triads

    itmp = ipvel
    do i = 1, nenod
       nndof = triads(i)%p%nDOFs
       nodeVar(3*i-2:3*i) = triads(i)%p%ur(:,4)
       if (nndof >= 6) then
          nodeVar(ipT:ipT+8) = reshape(triads(i)%p%ur(:,1:3),(/9/))
          ipT = ipT + 9
       end if
       if (nndof > 0) then
          nodeVar(ipvel:ipvel+nndof-1) = triads(i)%p%urd
          nodeVar(ipacc:ipacc+nndof-1) = triads(i)%p%urdd
          ipvel = ipvel + nndof
          ipacc = ipacc + nndof
       end if
    end do

    !! Update state-dependent variables, if any

    ipT = iwork(2) - iwork(3)
    do i = 1, size(engines)
       if (associated(engines(i)%p)) then
          rwork(ipT+i) = EngineValue(engines(i)%p,ierr)
       end if
    end do
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'UpdateUDE')
       return
    end if

    !! Update element state

    ipT   = 3*nenod + 1
    ipvel = itmp
    ipacc = itmp + nedof
    call Fi_UDE2 (elmId%baseId,etype,nenod,nedof,nodeVar(1),nodeVar(ipT), &
         &        nodeVar(ipVel),nodeVar(ipAcc),iwork(1),rwork(1), &
         &        K(1,1),C(1,1),M(1,1),Fs(1),Fd(1),Fi(1),Q(1), &
         &        t,dt,istep,iter,ierr)

    !! Update element coordinate system

    if (ierr >= 0) then
       call Fi_UDE3 (elmId%baseId,etype,nenod,nodeVar(1),nodeVar(ipT), &
            &        iwork(1),rwork(1),Tlg(1,1),ierr)
    end if

    if (nenod == 2 .and. envir%rhow > 0.0_dp .and. associated(hydyn)) then

       !! Calculate Morison forces and add to the external forces.
       !! Also add associated damping and added mass terms.

       if (dt > 0.0_dp) then
          call calcMorisonForces ('User-defined Element'//getId(elmId), &
               &                  nedof,triads,Tlg(:,1:3), &
               &                  nodeVar(ipVel:ipVel+nedof-1), &
               &                  nodeVar(ipAcc:ipAcc+nedof-1), &
               &                  C,M,Q,hydyn,envir, &
               &                  t,istep,abs(iter),ierr)
       else ! Quasi-static solution, calculate buoyancy only
          call calcBuoyancyForces ('User-defined Element'//getId(elmId), &
               &                   nedof,triads,Tlg(:,1:3), &
               &                   Q,hydyn,envir,t,abs(iter),ierr)
       end if

    end if

#ifdef FT_DEBUG
    write(dbgUDE,600) elmId%userId, iStep, iter
600 format(/'==== In UpdateUDE: iElm =',i6,'  iStep =',i6,' iter =',i4,' ===='/)
    call writeObject (K ,dbgUDE,'Userdef K')
    call writeObject (C ,dbgUDE,'Userdef C')
    call writeObject (M ,dbgUDE,'Userdef M')
    call writeObject (Fs,dbgUDE,'Userdef Fs')
    call writeObject (Fd,dbgUDE,'Userdef Fd')
    call writeObject (Fi,dbgUDE,'Userdef Fi')
    call writeObject (Q ,dbgUDE,'Userdef Q')
    call writeObject (Tlg,dbgUDE,'Userdef Tlg')
#endif

    if (ierr < 0) then
       call reportError (error_p,'Failed to update User-defined Element'// &
            &            getId(elmId),addString='UpdateUDE')
    end if

  end subroutine UpdateUDE


  !!============================================================================
  !> @brief Updates all user-defined elements in the model.
  !>
  !> @param elms Array of all user-defined elements in the model
  !> @param env Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] timeStep Time increment size
  !> @param[in] istep Time increment counter
  !> @param[in] iter Iteration counter
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine UpdateUDEs (elms,env,time,timeStep,istep,iter,ierr)

    use UserDefElTypeModule   , only : UserdefElType, dp
    use EnvironmentTypeModule , only : EnvironmentType
    use TriadTypeModule       , only : UpdateNodeForce
    use reportErrorModule     , only : reportError, error_p, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue

    type(UserdefElType)  , intent(inout) :: elms(:)
    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(in)    :: time, timeStep
    integer              , intent(in)    :: istep, iter
    integer              , intent(inout) :: ierr

    !! Local variables
    integer  :: i, j, n, m, jerr, lerr
    real(dp) :: RelmToSys(3,3)
    real(dp), target  :: dummy(6)
    real(dp), pointer :: f(:,:), fn(:)

    !! --- Logic section ---

    if (ffa_cmdlinearg_intValue('dragForm') < 1) then
       ierr = -1
       call reportError (error_p,'-dragForm=0 is not supported '// &
            &            'for user-defined elements', &
            &            addString='UpdateUDEs')
       return
    end if

    jerr = ierr
    do i = 1, size(elms)

       !! Update the state of the user-defined element

       call updateUDE (elms(i)%id,elms(i)%type,elms(i)%nTotDofs, &
            &          elms(i)%triads,elms(i)%engines, &
            &          elms(i)%iwork,elms(i)%rwork,elms(i)%hydyn,env, &
            &          elms(i)%Kmat,elms(i)%Cmat,elms(i)%Mmat, &
            &          elms(i)%FS,elms(i)%FD,elms(i)%FI,elms(i)%Q,elms(i)%Tlg, &
            &          time,timeStep,istep,iter,lerr)
       if (lerr < 0) ierr = ierr + lerr

       !! Transform element forces to system coordinates and update nodal forces

       do j = 1, size(elms(i)%triads)
          if (elms(i)%triads(j)%p%nDOFs < 3) cycle

          n = elms(i)%triads(j)%firstDOF
          m = n-1 + elms(i)%triads(j)%p%nDOFs
          if (associated(elms(i)%triads(j)%p%sysDirInG) .and. m-n >= 2) then
             RelmToSys = transpose(elms(i)%triads(j)%p%sysDirInG)
             elms(i)%Q (n:n+2) = matmul(RelmToSys,elms(i)%Q (n:n+2))
             elms(i)%FS(n:n+2) = matmul(RelmToSys,elms(i)%FS(n:n+2))
             if (timeStep > 0.0_dp) then
                elms(i)%FD(n:n+2) = matmul(RelmToSys,elms(i)%FD(n:n+2))
                elms(i)%FI(n:n+2) = matmul(RelmToSys,elms(i)%FI(n:n+2))
             end if
             if (m-n >= 5) then
                elms(i)%Q (n+3:n+5) = matmul(RelmToSys,elms(i)%Q (n+3:n+5))
                elms(i)%FS(n+3:n+5) = matmul(RelmToSys,elms(i)%FS(n+3:n+5))
                if (timeStep > 0.0_dp) then
                   elms(i)%FD(n+3:n+5) = matmul(RelmToSys,elms(i)%FD(n+3:n+5))
                   elms(i)%FI(n+3:n+5) = matmul(RelmToSys,elms(i)%FI(n+3:n+5))
                end if
             end if
          end if

          if (associated(elms(i)%triads(j)%p%supForce)) then
             fn => elms(i)%triads(j)%p%supForce
          else
             fn => dummy(1:elms(i)%triads(j)%p%nDOFs)
          end if
          fn = elms(i)%FS(n:m) - elms(i)%Q(n:m)
          if (timeStep > 0.0_dp) then
             fn = fn + elms(i)%FD(n:m) + elms(i)%FI(n:m)
          end if

          call updateNodeForce (elms(i)%triads(j)%p,fn)

          if (associated(elms(i)%triads(j)%p%force)) then
             f => elms(i)%triads(j)%p%force
             f(:,1) = f(:,1) + elms(i)%FS(n:m)
             if (timeStep > 0.0_dp) then
                f(:,2) = f(:,2) + elms(i)%FD(n:m)
                f(:,3) = f(:,3) + elms(i)%FI(n:m)
             end if
             f(:,4) = f(:,4) + elms(i)%Q(n:m)
          end if

       end do

    end do

    if (ierr < jerr) call reportError (debugFileOnly_p,'UpdateUDEs')

  end subroutine UpdateUDEs


  !!============================================================================
  !> @brief Calculates Morison force contributions for a two-noded element.
  !>
  !> @param[in] elmId ID string of the element to calculate for
  !> @param[in] nedof Dimension of the element matrices
  !> @param[in] triads The triads connected to this user-defined element
  !> @param[in] Tlg Local-to-global transformation matrix for the element
  !> @param[in] urd Element velocity vector, in global coordinates
  !> @param[in] urdd Element acceleration vector, in global coordinates
  !> @param C Element damping matrix
  !> @param M Element mass matrix
  !> @param Q External element forces
  !> @param hydyn Data for hydrodynamic force calculation for this element
  !> @param envir Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] istep Time increment counter
  !> @param[in] iter Iteration counter
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Sep 2019

  subroutine calcMorisonForces (elmId,nedof,triads,Tlg,urd,urdd,C,M,Q, &
       &                        hydyn,envir,time,istep,iter,ierr)

    use TriadTypeModule      , only : TriadPtrType, dp
    use SupElTypeModule      , only : HydroDynType
    use EnvironmentTypeModule, only : EnvironmentType
    use HydrodynamicsModule  , only : getMorisonForces
    use reportErrorModule    , only : reportError, debugFileOnly_p

    character(len=*)     , intent(in)    :: elmId
    type(TriadPtrType)   , intent(in)    :: triads(:)
    real(dp)             , intent(in)    :: Tlg(:,:), urd(:), urdd(:), time
    real(dp)             , intent(inout) :: C(:,:), M(:,:), Q(:)
    type(HydroDynType)   , intent(inout) :: hydyn
    type(EnvironmentType), intent(inout) :: envir
    integer              , intent(in)    :: nedof, istep, iter
    integer              , intent(inout) :: ierr

    !! Local variables
    integer  :: i, jerr, lerr
    real(dp) :: vl(nedof), al(nedof), Ql(nedof), tmp(3)
    real(dp) :: Cd(nedof,nedof), Ma(nedof,nedof)

    !! --- Logic section ---

    lerr = ierr

    !! Transform the velocity and acceleration to local element axes
    do i = 1, nedof-2, 3
       vl(i:i+2) = matmul(urd(i:i+2),Tlg)
       al(i:i+2) = matmul(urdd(i:i+2),Tlg)
       Ql(i:i+2) = matmul(Q(i:i+2),Tlg)
    end do

    !! Calculate the hydrodynamic forces on this element
    call getMorisonForces (elmId,triads,Tlg,vl,al,Ql,Ma,Cd, &
         &                 hydyn,envir,time,istep,iter,jerr)
    if (jerr < 0) then
       ierr = ierr + jerr
    else if (jerr > 0) then
       return ! No hydrodynamics (yet)
    end if

    !! Transform added mass and drag damping matrices to global coordinate axes
    do i = 1, nedof-2, 3
       call MATTRA (Tlg(1,1),Ma(1,1),tmp(1),nedof,3,i,0,jerr)
       if (jerr < 0) ierr = ierr + jerr
       call MATTRA (Tlg(1,1),Cd(1,1),tmp(1),nedof,3,i,0,jerr)
       if (jerr < 0) ierr = ierr + jerr
    end do
    if (ierr < lerr) then
       call reportError (debugFileOnly_p,'calcMorisonForces')
       return
    end if

    !! Add to the global element matrices
    C = C + Cd
    M = M + Ma
    do i = 1, nedof-2, 3
       Q(i:i+2) = matmul(Tlg,Ql(i:i+2))
    end do

  end subroutine calcMorisonForces


  !!============================================================================
  !> @brief Calculates buoyancy force contributions for a two-noded element.
  !>
  !> @param[in] elmId ID string of the element to calculate for
  !> @param[in] nedof Dimension of the element matrices
  !> @param[in] triads The triads connected to this user-defined element
  !> @param[in] Tlg Local-to-global transformation matrix for the element
  !> @param Q External element forces
  !> @param hydyn Data for hydrodynamic force calculation for this element
  !> @param envir Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] iter Iteration counter
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Sep 2019

  subroutine calcBuoyancyForces (elmId,nedof,triads,Tlg,Q, &
       &                         hydyn,envir,time,iter,ierr)

    use TriadTypeModule      , only : TriadPtrType, dp
    use SupElTypeModule      , only : HydroDynType
    use EnvironmentTypeModule, only : EnvironmentType
    use HydrodynamicsModule  , only : getBuoyancyForces
    use reportErrorModule    , only : reportError, debugFileOnly_p

    character(len=*)     , intent(in)    :: elmId
    type(TriadPtrType)   , intent(in)    :: triads(:)
    real(dp)             , intent(in)    :: Tlg(:,:), time
    real(dp)             , intent(inout) :: Q(:)
    type(HydroDynType)   , intent(inout) :: hydyn
    type(EnvironmentType), intent(inout) :: envir
    integer              , intent(in)    :: nedof, iter
    integer              , intent(inout) :: ierr

    !! Local variables
    integer  :: i, lerr
    real(dp) :: g, Ql(nedof)

    !! --- Logic section ---

    !! Transform the external load vector to local element axes
    do i = 1, nedof-2, 3
       Ql(i:i+2) = matmul(Q(i:i+2),Tlg)
    end do

    !! Add buoyancy forces this element
    call getBuoyancyForces (elmId,triads,Tlg,hydyn,envir,g,time,iter,Ql,lerr)
    if (lerr < 0) then
       ierr = ierr + lerr
       call reportError (debugFileOnly_p,'calcBuoyancyForces')
       return
    end if

    !! Transform back to global coordinate axes
    do i = 1, nedof-2, 3
       Q(i:i+2) = matmul(Tlg,Ql(i:i+2))
    end do

  end subroutine calcBuoyancyForces


  !!============================================================================
  !> @brief Adds element forces into corresponding system force vectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] elms All user-defined elements in the model
  !> @param FSk Internal stiffness force vector
  !> @param FDk Internal damping force vector
  !> @param FIk Internal inertia force vector
  !> @param Qk  External force vector (gravitation loads)
  !> @param RFk Reaction forces
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine AddInUDEForces (sam,elms,FSk,FDk,FIk,Qk,RFk,ierr)

    use SamModule          , only : SamType
    use UserdefElTypeModule, only : UserdefElType, dp
    use IdTypeModule       , only : getId
    use AsmExtensionModule , only : csAddEV
    use reportErrorModule  , only : reportError, error_p

    type(SamType)      , intent(in)    :: sam
    type(UserdefElType), intent(in)    :: elms(:)
    real(dp)           , intent(inout) :: FSk(:), FDk(:), FIk(:), Qk(:), RFk(:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i, iel, err, nedof

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(elms)
       iel   = elms(i)%samElNum
       nedof = elms(i)%nTotDofs
       call csAddEV (sam, iel, nedof, elms(i)%FS, FSk, RFk, err); ierr=ierr+err
       call csAddEV (sam, iel, nedof, elms(i)%FD, FDk, RFk, err); ierr=ierr+err
       call csAddEV (sam, iel, nedof, elms(i)%FI, FIk, RFk, err); ierr=ierr+err
       call csAddEV (sam,-iel, nedof, elms(i)%Q , Qk,  RFk, err); ierr=ierr+err
       if (ierr /= 0) then
          call reportError (error_p,'Inconsistent control information.', &
               'Encountered during assembly of User-defined Element'// &
               getId(elms(i)%id),addString='AddInUDEForces')
          return
       end if
    end do

  end subroutine AddInUDEForces


  !!============================================================================
  !> @brief Adds static element forces into corresponding system force vectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] elms All user-defined elements in the model
  !> @param FSk Internal stiffness force vector
  !> @param Qk  External force vector (gravitation loads)
  !> @param RFk Reaction forces
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine AddInStaticUDEForces (sam,elms,FSk,Qk,RFk,ierr)

    use SamModule          , only : SamType
    use UserdefElTypeModule, only : UserdefElType, dp
    use IdTypeModule       , only : getId
    use AsmExtensionModule , only : csAddEV
    use reportErrorModule  , only : reportError, error_p

    type(SamType)      , intent(in)    :: sam
    type(UserdefElType), intent(in)    :: elms(:)
    real(dp)           , intent(inout) :: FSk(:)
    real(dp)           , intent(inout) :: Qk(:)
    real(dp)           , intent(inout) :: RFk(:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i, iel, err, nedof

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(elms)
       iel   = elms(i)%samElNum
       nedof = elms(i)%nTotDofs
       call csAddEV (sam, iel, nedof, elms(i)%FS, FSk, RFk, err); ierr=ierr+err
       call csAddEV (sam,-iel, nedof, elms(i)%Q , Qk,  RFk, err); ierr=ierr+err
       if (ierr /= 0) then
          call reportError (error_p,'Inconsistent control information.', &
               'Encountered during assembly of User-defined Element'// &
               getId(elms(i)%id),addString='AddInStaticUDEForces')
          return
       end if
    end do

  end subroutine AddInStaticUDEForces


  !!============================================================================
  !> @brief Adds an element matrix into the equivalent system matrix.
  !>
  !> @param[in] udeMat The element matrix to add into the system
  !> @param sysMat System coefficiant matrix (stiffness, mass, ...)
  !> @param[in] elm The user-defined element the matrix is associated with
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !> @param sysRhs System right-hand-side vector
  !>
  !> @details This subroutine transforms the user-defined element matrix
  !> @a udeMat to the system directions for each DOF associated with it,
  !> and then adds it into the corresponding system matrix @a sysMat.
  !> If any of the element DOFs are prescribed, the associated force
  !> contributions are added into the system force vector @a sysRhs, if present.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine AddInUDEMat (udeMat,sysMat,elm,sam,ierr,sysRhs)

    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType
    use UserdefElTypeModule, only : UserdefElType, dp
    use AsmExtensionModule , only : csAddEM
    use scratchArrayModule , only : getRealScratchMatrix
    use reportErrorModule  , only : getErrorFile, reportError, debugFileOnly_p

    real(dp)           , intent(in)    :: udeMat(:,:)
    type(SysMatrixType), intent(inout) :: sysMat
    type(UserdefElType), intent(in)    :: elm
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: sysRhs(:)

    !! Local variables
    integer           :: i, n, nDim, lpu
    real(dp)          :: RelmToSys(3,3), tmp(3)
    real(dp), pointer :: uMat(:,:)

    !! --- Logic section ---

    lpu = getErrorFile()
    nDim = size(udeMat,1)
    uMat => getRealScratchMatrix(nDim,nDim,ierr)
    if (ierr < 0) goto 900

    uMat = udeMat
    do i = 1, size(elm%triads)
       if (.not. associated(elm%triads(i)%p%sysDirInG)) cycle

       n = elm%triads(i)%firstDOF
       RelmToSys = transpose(elm%triads(i)%p%sysDirInG)

       if (elm%triads(i)%p%nDOFs >= 3) then
          call MATTRA (RelmToSys(1,1),uMat(1,1),tmp(1),nDim,3,n,lpu,ierr)
          if (ierr < 0) goto 900
       end if

       if (elm%triads(i)%p%nDOFs >= 6) then
          call MATTRA (RelmToSys(1,1),uMat(1,1),tmp(1),nDim,3,n+3,lpu,ierr)
          if (ierr < 0) goto 900
       end if

    end do

    call csAddEM (sam,elm%samElNum,elm%nTotDofs,uMat,sysMat,ierr,sysRhs)
    if (ierr == 0) return

900 continue
    call reportError (debugFileOnly_p,'AddInUDEMat')

  end subroutine AddInUDEMat


  !!============================================================================
  !> @brief Computes the user-defined element Newton matrix
  !>
  !> @param[in] scaleM Mass matrix scaling factor
  !> @param[in] scaleC Damping matrix scaling factor
  !> @param[in] scaleK Stiffness matrix scaling factor
  !> @param elm The user-defined element to compute Newton matrix for
  !> @param[out] ierr Error flag
  !>
  !> @details The Newton matrix is a linear combination of the stiffness-,
  !> damping- (if any) and stiffness matrices of the user-defined element.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine BuildUDENewtonMat (scaleM,scaleC,scaleK,elm,ierr)

    use UserdefElTypeModule, only : UserdefElType, dp

    real(dp)           , intent(in)    :: scaleM, scaleC, scaleK
    type(UserdefElType), intent(inout) :: elm
    integer            , intent(out)   :: ierr

    !! --- Logic section ---

    ierr = 0
    if (associated(elm%Cmat)) then
       call BuildNewtonMat (elm%Nmat,elm%Mmat,elm%Cmat,elm%Kmat)
    else
       call BuildUndampedNewtonMat (elm%Nmat,elm%Mmat,elm%Kmat)
    end if

  contains

    !> @brief Calculates the Newton matrix with damping.
    subroutine BuildNewtonMat (N,M,C,K)
      real(dp), intent(out) :: N(:,:)
      real(dp), intent(in)  :: M(:,:), C(:,:), K(:,:)
      N = scaleM*M + scaleC*C + scaleK*K
    end subroutine BuildNewtonMat

    !> @brief Calculates the Newton matrix without damping.
    subroutine BuildUndampedNewtonMat (N,M,K)
      real(dp), intent(out) :: N(:,:)
      real(dp), intent(in)  :: M(:,:), K(:,:)
      N = scaleM*M + scaleK*K
    end subroutine BuildUndampedNewtonMat

  end subroutine BuildUDENewtonMat

end module UserdefElRoutinesModule
