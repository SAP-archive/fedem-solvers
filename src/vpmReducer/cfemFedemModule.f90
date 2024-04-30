!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module cfemFedemModule
#ifdef FT_HAS_CFEM

  implicit none

contains

  subroutine readReducerDataAndSolveCfem (sam,sysK,sysM,tenc,medof, &
       &                      diagMass,lumpedMass,ipsw,lpu,ierr, &
       &                      stiff_f2c,force_disp_f2c, &
       &                      numStiff, numStiffCalculated, &
       &                      additionalCfemInput)

    !!==========================================================================
    !! Administer data input and preprocessing for fedem_reducer.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Sep 2000/1.0
    !!==========================================================================

    use SamModule              , only : SamType, dp
    use SamModule              , only : writeObject, writeSamSize
    use SamReducerModule       , only : initiateSAM
    use SysMatrixTypeModule    , only : SysMatrixType
    use SysMatrixTypeModule    , only : diagonalmatrix_p, skylinematrix_p
    use SysMatrixTypeModule    , only : nullifySysMatrix, shareMatrixStructure
    use SysMatrixTypeModule    , only : check, writeObject
    use AsmExtensionModule     , only : csAllocPointerArrays
    use AsmExtensionModule     , only : csRenumber, csPreassemble
    use ManipMatrixModule      , only : writeObject
    use AllocationModule       , only : reAllocate
    use ProgressModule         , only : writeProgress
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_massprop
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getbool
    use InputReducerModule     , only : SUBPRIN
    use InputReducerModule     , only : checkSingularNodes, checkLinearCouplings

    type(SamType)      , intent(out) :: sam
    type(SysMatrixType), intent(out) :: sysK, sysM
    real(dp)           , pointer     :: tenc(:,:)
    integer            , pointer     :: medof(:)
    logical            , intent(in)  :: diagMass, lumpedMass
    integer            , intent(in)  :: ipsw, lpu
    integer            , intent(out) :: ierr

    real(dp),            pointer     :: stiff_f2c(:), force_disp_f2c(:,:)
    integer,             intent(in)  :: numStiff
    integer,             intent(out) :: numStiffCalculated
    character(len=*),    intent(inout)  :: additionalCfemInput

    !! Local variables
    logical            :: factorMass
    integer            :: i, nenod, neval, ngen, ndim
    real(dp)           :: mass, cg(3), inertia(6)
    real(dp), pointer  :: tncoor(:,:)
    character(len=128) :: errMsg

    !! --- Logic section ---

    ! --- Initialize the SAM datastructure
    nullify(tncoor)
    call initiateSAM (sam,tncoor,.false.,ipsw,lpu,ierr)
    if (ierr /= 0) goto 900

    sam%mpar(40) = numStiff
    call ffa_cmdlinearg_getint ('linkId',sam%mpar(18))
    call ffa_cmdlinearg_getint ('ngen',ngen)
    if (ngen > 0) then
       call ffa_cmdlinearg_getint ('neval',neval)
       call ffa_cmdlinearg_getbool ('factorMass',factorMass)
       sam%mpar(19) = max(neval,ngen)
    else
       factorMass = .false.
    end if

    ! --- Compute the number of external nodes and the reduced matrix dimension
    nenod = 0
    ndim = ngen
    do i = 1, sam%nnod
       if (sam%mnnn(i) > 1) then
          nenod = nenod + 1
          ndim = ndim + sam%madof(i+1)-sam%madof(i)
       end if
    end do

    sam%mpar(22) = ngen
    sam%mpar(24) = ndim
    sam%mpar(29) = nenod

    write(lpu,600) sam%nel,sam%nnod,sam%mpar(29),sam%mpar(22), &
         &         sam%mpar(17),sam%mpar(27),sam%mpar(24)

    ! --- Compute and report mass properties of the unreduced FE model
    call writeProgress (' --> Computing mass properties')
    call ffl_massprop (mass,cg,inertia)
    write(lpu,610) mass,cg,inertia(1),inertia(4:5), &
         &         inertia(2),inertia(6),inertia(3)
    call flush(lpu)

600 format(//4X,'FE MODEL SUMMARY' / 4X,16('-') &
         & //4X,'Number of elements             =',I8 &
         &  /4X,'Number of nodes (total)        =',I8 &
         &  /4X,'Number of external nodes       =',I8 &
         &  /4X,'Number of generalized modes    =',I8 &
         &  /4X,'Number of materials            =',I8 &
         &  /4X,'Number of geometric properties =',I8 &
         &  /4X,'Dimension of reduced matrices  =',I8 / )
610 format(  4X,'Total mass                     =',1PE13.5 &
         &  /4X,'Center of gravity              =',3E13.5 &
         &  /4X,'Moments of inertia (tensor)    =',3E13.5 &
         &  /49X,2E13.5 /38X,'symm.',19X,E13.5 / )

    ! --- Store away nodal coordinates of the external nodes (needed in JCMS)
    call reAllocate ('readReducerData',tenc,3,nenod,ierr)
    call reAllocate ('readReducerData',medof,nenod+1,ierr)
    if (ierr < 0) goto 900

    nenod = 0
    medof(1) = 1
    do i = 1, sam%nnod
       if (sam%mnnn(i) > 1) then
          nenod = nenod + 1
          tenc(:,nenod) = tncoor(i,:)
          medof(nenod+1) = medof(nenod) + sam%madof(i+1)-sam%madof(i)
       end if
    end do

    ! --- Check if any of the nodes lack stiffness or boundary conditions
    call writeProgress (' --> Checking model consistency')
    call CheckSingularNodes (sam%nel,sam%nnod,sam%minex, &
         &                   sam%mpmnpc,sam%mmnpc,sam%madof,sam%msc,ierr)
    if (ierr > 0) then
       if (ipsw > 0) then
          call writeObject(sam%minex,lpu,'MINEX')
          call writeObject(sam%mpmnpc,lpu,'MPMNPC')
          call writeObject(sam%mmnpc,lpu,'MMNPC')
          call writeObject(sam%madof,lpu,'MADOF')
          call writeObject(sam%msc,lpu,'MSC')
       end if
       goto 990
    end if

    if (associated(sam%mmceq)) then

       ! --- Check linear couplings
       call CheckLinearCouplings (sam%mpmceq,sam%mmceq,sam%msc,ierr)
       if (ierr > 0) then
          if (ipsw > 0) then
             call writeObject(sam%mpmceq,lpu,'MPMCEQ')
             call writeObject(sam%mmceq,lpu,'MMCEQ')
             call writeObject(sam%msc,lpu,'MSC')
          end if
          goto 990
       end if
    end if

    ! --- Allocate the pointer structure for the linear equation solver
    call writeProgress (' --> Creating data structures for equation solver')
    call ffa_cmdlinearg_getint ('cachesize',ngen)
    call csAllocPointerArrays (sam,sysK,skylineMatrix_p,lpu,ierr,ipsw)
    if (ierr < 0) goto 900

    ! --- Optimize the equation system
    call csRenumber (sam,sysK,lpu,ierr)
    if (ierr < 0) goto 900

    ! --- Preassembly
    call writeProgress (' --> Pre-assembling the stiffness matrix')
    call csPreassemble (sam,sysK,.true.,lpu,ierr,ipsw)
    if (ierr < 0) goto 900

    ! --- Check for consistency
    call check (sysK,sam%mpar,lpu,ierr)
    if (ierr < 0) goto 900

    if (.not.associated(sam%meqn)) sam%meqn => sysK%meqn

    ! --- Create the mass matrix data structure

    if (.not.diagMass .and. .not.lumpedMass) then

       ! Let the mass- and stiffness matrices share datastructures
       call shareMatrixStructure (sysM,sysK,factorMass,ierr=ierr)
       if (ierr < 0) goto 900

    else if (sam%nmmceq > 0 .and. lumpedMass) then

       ! Lumped mass is requested for a link that has constraint equations.
       ! The system mass matrix will then have a different sparsity pattern
       ! than the system stiffness matrix.

       call csAllocPointerArrays (sam,sysM,-sysK%storageType,lpu,ierr, &
            &                     ipsw,ngen,factorMass,sysK%gsf)
       if (ierr < 0) goto 900

       ! --- Optimize the equation system
       call csRenumber (sam,sysM,lpu,ierr)
       if (ierr < 0) goto 900

       ! --- Preassembly
       call writeProgress (' --> Pre-assembling the lumped mass matrix')
       call csPreassemble (sam,sysM,.false.,lpu,ierr,ipsw)
       if (ierr < 0) goto 900

       ! --- Check for consistency
       call check (sysM,sam%mpar,lpu,ierr)
       if (ierr < 0) goto 900

    else

       ! Lumped mass is requested for a link that has no constraint equations.
       ! The system mass matrix will then be a pure diagonal matrix.
       call shareMatrixStructure (sysM,sysK,factorMass,diagonalMatrix_p)

    end if


    !! Solve link in nonlinear solver (CFEM)
    nullify(stiff_f2c)
    nullify(force_disp_f2c)
    call cfem_model (sam, sysK%skyline%msky, tncoor, &
         &           stiff_f2c, force_disp_f2c, numStiffCalculated, &
         &           additionalCfemInput, ierr)
    if (ierr < 0) goto 900

    if (ipsw > 0) then

       ! --- Print out substructure data to res-file
       call SUBPRIN (sam%mpar(18),sam%mpar,sam%mpmnpc,sam%mmnpc,sam%melcon, &
            &        sam%minex,tncoor(:,1),tncoor(:,2),tncoor(:,3),sam%madof, &
            &        sam%msc,sam%mpmceq,sam%mmceq,sam%ttcc,sam%mnnn,LPU)

    end if

    call writeSAMsize (sam,lpu)

    if (sam%ndof2 < 1) then
       ierr = -1
       call reportError (error_p,'This FE part has no external DOFs.', &
            'Check that it is properly attached to other model members.')
       goto 900
    else if (sam%mpar(22) > sam%ndof1) then
       ierr = -1
       write(errMsg,690) sam%mpar(22), sam%ndof1
690    format('NGEN =',I8,',  NDOF1 =',I8)
       call reportError (error_p,'The number of generalized DOFs (NGEN)', &
            'is larger than the number of internal nodal DOFs (NDOF1).', &
            errMsg,'This is not allowed, reduce NGEN to less than NDOF1.')
       goto 900
    end if

900 continue
    if (ipsw > 1) then
       call writeObject (sam,lpu,4)
       call writeObject (sysK,sam%mpar,lpu,' === System stiffness matrix',10,2)
       call writeObject (sysM,sam%mpar,lpu,' === System mass matrix',10,2)
    end if
    if (ipsw >= -7 .and. ipsw <= -5) call writeObject (sam,lpu,5)
990 call reAllocate ('readReducerData',tncoor)
    if (ierr /= 0) call reportError (debugFileOnly_p,'readReducerData')

  end subroutine readReducerDataAndSolveCfem


  subroutine cfem_model(sam, msky, tncoor, stiff_f2c, force_disp_f2c, &
       &                num_stiff, additonalCfemInput, ierr)

    use SamModule              , only : SamType, dp
    use ProgressModule         , only : lterm
    use ReportErrorModule      , only : allocationError
    use FFlLinkHandlerInterface, only : ffl_getmat, ffl_getthick
    use FFlLinkHandlerInterface, only : ffl_getBeamSection

    real(dp),pointer              :: tncoor(:,:)
    real(dp),pointer              :: stiff_f2c(:)
    real(dp),pointer              :: force_disp_f2c(:,:)
    type(SamType)   , intent(in)  :: sam
    integer         , intent(in)  :: msky(:)
    integer         , intent(out) :: num_stiff, ierr
    character(len=*), intent(inout) :: additonalCfemInput

    !! Local variables
    integer  :: iel, iMat, iProp, nenod, nMat, nProp, bRedirectStdout
    real(dp) :: e, E_prev, nu, nu_prev, rho, prop(14), prop_prev(14)
    real(dp), allocatable :: matData(:,:), propData(:,:)
    integer,  allocatable :: elMat(:), elProp(:)
    real(dp), allocatable :: force_f2c(:), disp_f2c(:)




    !**************f2c****** Materials to CFEM *****************
    nMat = 0
    E_prev  = -1.0_dp
    nu_prev = -1.0_dp

    do iel = 1, sam%nel
       call ffl_getmat (e,nu,rho,iel,ierr)
       if ( e /= E_prev .or. nu /= nu_prev ) then
          nMat = nMat + 1
          E_prev = e
          nu_prev = nu
       end if
    end do

    allocate(elMat(sam%nel),matData(3,nMat),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('cfem_model 1')
       return
    end if

    elMat   = 0
    matData = 0.0_dp
    E_prev  = -1.0_dp
    nu_prev = -1.0_dp

    iMat = 0
    do iel = 1, sam%nel
       call ffl_getmat (e,nu,rho,iel,ierr)
       if ( e /= E_prev .or. nu /= nu_prev ) then
          iMat=iMat+1
          matData(1,imat)=e
          matData(2,imat)=nu
          matData(3,imat)=rho
          E_prev = e
          nu_prev = nu
       end if
       elMat(iel)  = iMat
    end do


    !**************f2c****** Properties to CFEM *****************


    nProp = 0
    prop = 0.0_dp
    prop_prev = -1.0_dp

    do iel = 1, sam%nel
       nenod=sam%mpmnpc(iel+1)-sam%mpmnpc(iel)
       if( nenod == 2 ) then ! beam element
          call ffl_getBeamSection (prop,iel,ierr)
          if ( prop(4) /= prop_prev(4) .or. &
               prop(5) /= prop_prev(5) .or. &
               prop(6) /= prop_prev(6) .or. &
               prop(7) /= prop_prev(7) .or. &
               prop(8) /= prop_prev(8) ) then
             nProp =nProp + 1
             prop_prev = prop
          end if
       else if( nenod == 3 .or. nenod == 4) then !shell element
          call ffl_getthick (prop(1),iel,ierr)
          if (prop(1) /= prop_prev(1)) then
             nProp =nProp + 1
             prop_prev(1) = prop(1)
          end if
       end if
    end do


    allocate(elProp(sam%nel),propData(14,nProp),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('cfem_model 2')
       return
    end if

    elProp    = 0
    propData  = 0.0_dp
    prop_prev = -1.0_dp

    iProp = 0
    do iel = 1, sam%nel
       nenod=sam%mpmnpc(iel+1)-sam%mpmnpc(iel)
       if( nenod == 2 ) then
          call ffl_getBeamSection (prop,iel,ierr)
          if ( prop(4) /= prop_prev(4) .or.  &
               prop(5) /= prop_prev(5) .or.  &
               prop(6) /= prop_prev(6) .or.  &
               prop(7) /= prop_prev(7) .or.  &
               prop(8) /= prop_prev(8) ) then
             iProp=iProp+1
             propData(:,iProp) = prop
             prop_prev = prop
          end if
       else if( nenod == 3 .or. nenod == 4) then !shell element
          call ffl_getthick (prop(1),iel,ierr)
          if (prop(1) /= prop_prev(1)) then
             iProp=iProp+1
             propData(1,iProp)=prop(1)
             prop_prev(1) = prop(1)
          end if
       end if
       elProp(iel) = iProp
    end do

    num_stiff = sam%mpar(40)
    write(lterm,"('---> Allocating f2c arrays num_stiff=',i3)") num_stiff
    allocate(stiff_f2c(sam%nsky*num_stiff), &
         &   force_f2c(sam%ndof*num_stiff), &
         &   disp_f2c(sam%ndof*num_stiff), &
         &   force_disp_f2c(sam%ndof*num_stiff,2), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('cfem_model 3')
       return
    end if

    stiff_f2c = 0.0_dp
    force_f2c = 0.0_dp
    disp_f2c  = 0.0_dp

    if ( lterm /= 6 ) then
       bRedirectStdout = 1
    else
       bRedirectStdout = 0
    end if

    call Link2CFEM(tncoor(1,1), tncoor(1,2), tncoor(1,3), &
         sam%mpar(1),sam%minex(1), sam%mpmnpc(1), sam%mmnpc(1), &
         sam%msc(1), sam%madof(1), sam%meqn(1), &
         msky(1), elMat(1), elProp(1), matData(1,1), propData(1,1), &
         nMat, nProp, stiff_f2c(1), force_f2c(1), disp_f2c(1), &
         num_stiff, num_stiff, additonalCfemInput, bRedirectStdout)

    force_disp_f2c(:,1)=force_f2c
    force_disp_f2c(:,2)=disp_f2c

    deallocate(force_f2c, disp_f2c)

  end subroutine cfem_model

#endif

end module cfemFedemModule
