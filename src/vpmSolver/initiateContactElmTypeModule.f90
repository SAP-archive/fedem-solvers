!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateContactElmTypeModule.f90
!> @brief Initialization of contact element objects from the solver input file.

!!==============================================================================
!> @brief Initialization of contact element objects from the solver input file.

module initiateContactElmTypeModule

  implicit none

  private :: initiateContactElement


contains

  !!============================================================================
  !> @brief Initializes contact elements with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triads Array of all triads in the model
  !> @param[in] cSurfs Array of all contact surfaces in the model
  !> @param[in] frictionSets Array of all friction parameter sets in the model
  !> @param[in] baseSprings Array of all spring objects in the model
  !> @param[in] baseDampers Array of all damper objects in the model
  !> @param[out] cElems Array of all contact elements in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date January 2002
  !>
  !> @author Knut Morten Okstad
  !> @date 6 Oct 2005

  subroutine ReadContactElements (infp,triads,cSurfs,frictionSets, &
       &                          baseSprings,baseDampers,cElems,err)

    use TriadTypeModule         , only : TriadType, dp, GetPtrToId
    use ContactSurfaceModule    , only : GliderCurveType, CurvePointType
    use ContactSurfaceModule    , only : GetPtrToId
    use ContactElementTypeModule, only : ContactElementType
    use ContactElementTypeModule, only : DeallocateContactElements
    use ContactElementTypeModule, only : NullifyContactElement, GetPtrToId
    use ContactElementTypeModule, only : ignoreEccN, ignoreEccF
    use FrictionTypeModule      , only : FrictionParameterType
    use FrictionTypeModule      , only : InitializeFriction
    use IdTypeModule            , only : IdType, ldesc_p
    use IdTypeModule            , only : initId, getId, ReportInputError
    use SpringTypeModule        , only : SpringBaseType, GetPtrToId
    use DamperTypeModule        , only : DamperBaseType, GetPtrToId
    use inputUtilities          , only : iuGetNumberOfEntries
    use inputUtilities          , only : iuSetPosAtNextEntry
    use progressModule          , only : lterm
    use reportErrorModule       , only : debugFileOnly_p, warning_p, note_p
    use reportErrorModule       , only : reportError, AllocationError
    use FFaCmdLineArgInterface  , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface  , only : ffa_cmdlinearg_getint

    integer                    , intent(in)  :: infp
    type(TriadType)            , intent(in)  :: triads(:)
    type(GliderCurveType)      , intent(in)  :: cSurfs(:)
    type(FrictionParameterType), intent(in)  :: frictionSets(:)
    type(SpringBaseType)       , intent(in)  :: baseSprings(:)
    type(DamperBaseType)       , intent(in)  :: baseDampers(:)
    type(ContactElementType)   , pointer     :: cElems(:)
    integer                    , intent(out) :: err

    !! Local variables
    integer :: i, idIn, nCElem, nDOFs, iprint, stat, lerr
    logical :: allSec, allRest, allJoint, allLen, allVel, ignoreElm

    !! Define the CONTACT_ELEMENT namelist
    character(len=30)  :: type
    character(ldesc_p) :: extDescr
    integer            :: id, extId(10), springId(6), damperId(6)
    integer            :: frictionSpringId, frictionSetId
    integer            :: followerTriad, contactSurface, saveVar(4)
    real(dp)           :: thickness, width, radius

    namelist /CONTACT_ELEMENT/ type, extDescr, id, extId, springId, damperId, &
         &                     frictionSpringId, frictionSetId, &
         &                     followerTriad, contactSurface, &
         &                     thickness, width, radius, saveVar

    !! --- Logic section ---

    nCElem = iuGetNumberOfEntries(infp,'&CONTACT_ELEMENT',err)
    if (err /= 0) goto 900
    write(lterm,*) 'Number of &CONTACT_ELEMENT =',nCElem

    call DeallocateContactElements (cElems)
    allocate(cElems(nCElem),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadContactElements 10')
       return
    end if

    call ffa_cmdlinearg_getbool ('ignoreCamEvecFgrad',ignoreEccN)
    call ffa_cmdlinearg_getbool ('ignoreCamEvecF',ignoreEccF)
    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allRestartVars',allRest)
    call ffa_cmdlinearg_getbool ('allJointVars',allJoint)
    call ffa_cmdlinearg_getbool ('allLengthVars',allLen)
    call ffa_cmdlinearg_getbool ('allVelVars',allVel)
    call ffa_cmdlinearg_getint ('debug',iprint)

    do idIn = 1, nCElem

       call NullifyContactElement (cElems(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&CONTACT_ELEMENT')) then
          err = err - 1
          call ReportInputError ('CONTACT_ELEMENT',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; type=''
       springId=0; damperId=0; frictionSpringId=0; frictionSetId=0
       followerTriad=0; contactSurface=0; saveVar=0
       thickness=-1.0_dp; width=-1.0_dp; radius=-1.0_dp

       read(infp,nml=CONTACT_ELEMENT,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('CONTACT_ELEMENT',idIn)
          cycle
       end if

       call initId (cElems(idIn)%id,id,extId,extDescr,stat)
       cElems(idIn)%thickness = thickness
       cElems(idIn)%width     = width
       cElems(idIn)%radius    = radius
       ignoreElm = .false.

       !! Initialize springs and dampers
       lerr = err
       do i = 1, size(cElems(idIn)%springs)
          if (springId(i) > 0) then
             cElems(idIn)%springs(i)%p => GetPtrToId(baseSprings,springId(i))
             if (.not.associated(cElems(idIn)%springs(i)%p)) then
                err = err - 1
                call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
                     &                 msg='Invalid base spring Id')
                cycle
             else if (cElems(idIn)%springs(i)%p%dof >= 0) then
                err = err - 1
                call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
                     msg='Reference to a base spring that already is in use')
                cycle
             end if
             cElems(idIn)%springs(i)%p%dof = 0
             cElems(idIn)%springs(i)%p%length => cElems(idIn)%cVar(i,1)
          else if (i == 1) then
             ignoreElm = .true.
             call reportError (warning_p, &
                  'Contact element'//getId(cElems(idIn)%id), &
                  'No contact spring specified in the local X-direction.', &
                  'This element is ignored by the Dynamics Solver.')
             exit
          end if
       end do

       do i = 1, size(cElems(idIn)%dampers)
          if (damperId(i) > 0 .and. .not.ignoreElm) then
             cElems(idIn)%dampers(i)%p => GetPtrToId(baseDampers,damperId(i))
             if (.not.associated(cElems(idIn)%dampers(i)%p)) then
                err = err - 1
                call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
                     &                 msg='Invalid base damper Id')
                cycle
             else if (cElems(idIn)%dampers(i)%p%isActive) then
                err = err - 1
                call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
                     msg='Reference to a base damper that already is in use')
                cycle
             end if
             if (associated(cElems(idIn)%dampers(i)%p%spr)) then
                if (associated(cElems(idIn)%springs(i)%p)) then
                   !! Connect to associated spring for deformational damping
                   cElems(idIn)%dampers(i)%p%spr => cElems(idIn)%springs(i)%p
                else
                   nullify(cElems(idIn)%dampers(i)%p%spr)
                   call reportError (note_p,'Switching off deformational'// &
                        ' damping for Contact Damper'// &
                        trim(getId(cElems(idIn)%dampers(i)%p%id))// &
                        ' since no associated spring')
                end if
             end if
             cElems(idIn)%dampers(i)%p%isActive = .true.
             if (radius > 0.0_dp .and. i < 3) then ! Radial contact
                cElems(idIn)%nonRadialDamping = .true. ! Use rectangular damping
                !! Must allocate separate velocity and length variables when the
                !! associated contact element variables are in the radial system
                allocate(cElems(idIn)%dampers(i)%p%velocity, &
                     &   cElems(idIn)%dampers(i)%p%length, STAT=stat)
                if (stat /= 0) then
                   err = AllocationError('ReadContactElements 20')
                   return
                end if
                cElems(idIn)%dampers(i)%p%velocity = 0.0_dp
                cElems(idIn)%dampers(i)%p%length   = 0.0_dp
             else if (associated(cElems(idIn)%dampers(i)%p%spr)) then
                !! Allocate separate velocity variable for deformational damper
                allocate(cElems(idIn)%dampers(i)%p%velocity, STAT=stat)
                if (stat /= 0) then
                   err = AllocationError('ReadContactElements 19')
                   return
                end if
                cElems(idIn)%dampers(i)%p%velocity = 0.0_dp
                cElems(idIn)%dampers(i)%p%length   => cElems(idIn)%cVar(i,1)
             else
                cElems(idIn)%dampers(i)%p%velocity => cElems(idIn)%cVar(i,2)
                cElems(idIn)%dampers(i)%p%length   => cElems(idIn)%cVar(i,1)
             end if
          end if
       end do
       if (err < lerr) cycle

       if (cElems(idIn)%nonRadialDamping) then
          !! Need a separate position matrix for the rectangular dampers
          !! when the springs are in radial/angular directions
          allocate(cElems(idIn)%DPosInG(3,3),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadContactElements 30')
             return
          end if
       else
          !! Springs and dampers share the same position matrix
          cElems(idIn)%DPosInG => cElems(idIn)%CPosInG
       end if

       !! Set the possible friction
       if (frictionSetId > 0 .and. .not.ignoreElm) then
          call InitializeFriction (cElems(idIn)%friction, &
               &                   frictionSets,frictionSetId,saveVar(3:4),stat)
          if (stat /= 0) then
             err = err - 1
             call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
                  &                 msg='Invalid friction set Id')
             cycle
          end if
       end if

       if (frictionSpringId > 0 .and. associated(cElems(idIn)%friction)) then
          cElems(idIn)%friction%spr => GetPtrToId(baseSprings,frictionSpringId)
          if (.not.associated(cElems(idIn)%friction%spr)) then
             err = err - 1
             call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
                  &                 msg='Invalid friction base spring Id')
             cycle
          else if (cElems(idIn)%friction%spr%dof >= 0) then
             err = err - 1
             call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
                  msg='Reference to a base spring that already is in use')
             cycle
          end if
          cElems(idIn)%friction%spr%dof = 3
          !! Connect the spring length to the slider motion
          cElems(idIn)%friction%spr%length => cElems(idIn)%cVar(3,1)
       end if

       !! Follower triad
       cElems(idIn)%triad1 => GetPtrToId(triads,followerTriad)
       if (.not.associated(cElems(idIn)%triad1)) then
          err = err - 1
          call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
               &                msg='Invalid follower triad Id')
          cycle
       end if
       !! Set samNodNum temporarily to flag that the triad should receive
       !! a proper node number later, also when the triad is grounded
       cElems(idIn)%triad1%samNodNum = 1

       !! Contact surface
       cElems(idIn)%cSurf => GetPtrToId(cSurfs,contactSurface)
       if (.not.associated(cElems(idIn)%cSurf)) then
          err = err - 1
          call ReportInputError ('CONTACT_ELEMENT',id=cElems(idIn)%id, &
               &                 msg='Invalid contact surface Id')
          cycle
       end if

       if (ignoreElm) cycle ! Ignore this element

       !! Define dimension of element matrices
       nDOFs = cElems(idIn)%triad1%nDOFs
       do i = 1, size(cElems(idIn)%cSurf%CPoints)
          nDOFs = nDOFs + cElems(idIn)%cSurf%CPoints(i)%triad%nDOFs
       end do
       cElems(idIn)%nDOFs = nDOFs

       !! Determine which contact variables will be saved. Settings in the
       !! solver input file are overruled by command-line arguments, if any.
       do i = 1, 2
          cElems(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec .or. allJoint
       end do
       if (allRest) cElems(idIn)%saveVar(1) = .true.
       if (allLen) cElems(idIn)%saveVar(1) = .true.
       if (allVel) cElems(idIn)%saveVar(2) = .true.

       !! Compute some initial contact properties
#ifdef FT_DEBUG
       call initiateContactElement (cElems(idIn),cElems(idIn)%cSurf, &
            &                       nCElem < 4 .or. extId(1) == -iprint)
#else
       call initiateContactElement (cElems(idIn),cElems(idIn)%cSurf)
#endif

    end do

900 if (err < 0) call reportError (debugFileOnly_p,'ReadContactElements')

  end subroutine ReadContactElements


  !!============================================================================
  !> @brief Initializes a contact element after reading the solver input file.
  !>
  !> @param cElem The contact element to initializes
  !> @param cSurf The contact surface which is used by @a cElem
  !> @param[in] doDebug If .true., print extra output to separate debug file
  !>
  !> @details This subroutine determines which contact curve domain the follower
  !> initially is in contact with, or which domain the follower initially is
  !> closest to if not initially in contact.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Oct 2005

  subroutine initiateContactElement (cElem,cSurf,doDebug)

    use ContactElementTypeModule, only : ContactElementType, dp
    use ContactSurfaceModule    , only : GliderCurveType
    use CurveTypeModule         , only : GetPosOnCircSection
    use IdTypeModule            , only : getId
    use reportErrorModule       , only : reportError, warning_p
#ifdef FT_DEBUG
    use dbgUnitsModule          , only : dbgCurve
    use fileUtilitiesModule     , only : getDBGfile
#endif

    type(ContactElementType), intent(inout) :: cElem
    type(GliderCurveType)   , intent(inout) :: cSurf
    logical, optional       , intent(in)    :: doDebug

    !! Local variables
    logical  :: isInDomain, isInPreviousDomain
    integer  :: i, j, nCurvePoints, nCircSections
    real(dp) :: Xsi, XsiPrev, dX, dXmin, dXPrev(6), CPosPrev(3,4)

    character(len=2), parameter :: dofs(3:6) = (/'Tz','Rx','Ry','Rz'/)

    !! --- Logic section ---

    if (present(doDebug)) then
#ifdef FT_DEBUG
       if (doDebug) then
          dbgCurve = getDBGfile(7,'contactCurve.dbg')
          write(dbgCurve,"(/'===== Contact Element',A)") trim(getId(cElem%id))
          write(dbgCurve,"( '      Follower',1P3E13.5)") cElem%triad1%ur(:,4)
       else
          dbgCurve = 0
       end if
#endif
    end if

    nCurvePoints = size(cSurf%CPoints)
    if (cSurf%isLooping) then
       nCircSections = nCurvePoints
    else
       nCircSections = nCurvePoints - 1
    end if
    isInDomain = .false.
    isInPreviousDomain = .false.

    dxMin    = 0.0_dp
    XsiPrev  = 0.0_dp
    dXPrev   = 0.0_dp
    CPosPrev = 0.0_dp
    if (cSurf%isExtended) then
       call reportError (warning_p,'Contact element'//getId(cElem%id), &
            'The isExtended flag is not yet supported (ignored).')
    end if

    !! Check all curve section domains sequentially
    do i = 1, nCircSections
       j = 1 + mod(i,nCurvePoints)
#ifdef FT_DEBUG
       if (dbgCurve > 0) write(dbgCurve,"(' ---- Checking domain:',2I4)") i, j
#endif

       call GetPosOnCircSection (cElem%triad1%ur, cElem%cVarPrev, &
            &                    cSurf%CPoints(i), cSurf%CPoints(j), &
            &                    cElem%thickness, cElem%width, cElem%radius, &
            &                    isInDomain, cElem%CPosInG, &
            &                    Xsi, cElem%cVar(:,1))

       dX = abs(cElem%cVar(1,1))
       if (isInDomain) then

          !! We are in contact with current domain
          if (isInPreviousDomain) then

             !! We are in contact with the previous domain too,
             !! now choose the one having the smallest dX value
             if (dX < abs(dXprev(1))) then
                cElem%lastDomain = i
             else
                cElem%lastDomain = i-1
                cElem%CPosInG    = CPosPrev
                cElem%cVar(:,1)  = dXPrev
                Xsi              = XsiPrev
             end if
             exit ! Exit the search loop

          else if (i == nCircSections .and. .not.cSurf%isLooping) then

             !! We are in the last domain of a non-looping contact curve
             cElem%lastDomain = i

          else

             !! This is the first domain that we detect contact for, we need to
             !! check the subsequent domain too to find which one is the closest
             isInPreviousDomain = .true.
             CPosPrev = cElem%CPosInG
             dXPrev   = cElem%cVar(:,1)
             XsiPrev  = Xsi

          end if

       else if (isInPreviousDomain) then

          !! We are in contact with the previous domain only
          isInDomain       = .true.
          cElem%lastDomain = i-1
          cElem%CPosInG    = CPosPrev
          cElem%cVar(:,1)  = dXPrev
          Xsi              = XsiPrev
          exit ! Exit the search loop

       else if (i == 1 .or. dX < dXmin) then

          !! Not yet in contact, but this domain is the closest so far...
          dXmin = dX
          cElem%lastDomain = i

       end if

    end do

    if (isInDomain) then
       !! Set the initial slider variable
       i = cElem%lastDomain
       j = 1 + mod(i,nCurvePoints)
       cElem%cVar(3,1) = (1.0_dp-Xsi) * cSurf%CPoints(i)%slideVar &
            &          +         Xsi  * cSurf%CPoints(j)%slideVar
       cElem%cVarPrev = cElem%cVar(1:3,1)
       cElem%isActive = .true.

       !! Establish the actual position matrix of the follower with respect to
       !! the follower triad, such that the rotational contact variables
       !! initially become zero (in case rotational springs have been assigned)
       cElem%FPosInT(:,1:3) = matmul(transpose(celem%triad1%ur(:,1:3)), &
            &                        cElem%CPosInG(:,1:3))
       cElem%cVar(4:6,1) = 0.0_dp

       !! Initiate the stress free lengths of the slider and friction springs,
       !! if any, such that they are stress free in the modeling configuration
       if (associated(cElem%springs(3)%p)) then
          cElem%springs(3)%p%l0 = cElem%springs(3)%p%l0 + cElem%cVar(3,1)
       end if
       if (associated(cElem%friction)) then
          if (associated(cElem%friction%spr)) then
             cElem%friction%spr%l0 = cElem%cVar(3,1)
          end if
       end if

#ifdef FT_DEBUG
       if (dbgCurve > 0) then
          write(dbgCurve,"('Found domain i,j, Xsi, dX,dY,dZ :',2I3,1P4E12.3)") &
               &                         i,j, Xsi, cElem%cVar(1:3,1)
       end if
#endif
    else
       !! Permanently deactivate slider and rotational springs when the follower
       !! is initially outside the domain, since l0 then is unknown
       do i = 3, 6
          if (associated(cElem%springs(i)%p)) then
             cElem%springs(i)%p%dof = -i
             nullify(cElem%springs(i)%p)
             call reportError (warning_p,'Contact element'//getId(cElem%id), &
                  'The '//dofs(i)//'-spring is ignored because '// &
                  'the follower is initially outside the domain.')
          end if
       end do

#ifdef FT_DEBUG
       if (dbgCurve > 0) then
          write(dbgCurve,"('Not in domain, closest =',I4,1PE12.3)") &
               &                           cElem%lastDomain, dXmin
       end if
#endif
    end if

  end subroutine initiateContactElement

end module initiateContactElmTypeModule
