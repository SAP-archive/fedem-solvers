!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module initiateBushingElmTypeModule

  implicit none

contains

  subroutine ReadBushingElements (infp,triads,baseSprings,baseDampers, &
       &                          bElems,err)

    !!==========================================================================
    !! Initialize the bushing elements with data from the solver input file.
    !!
    !! Programmer : Knut Morten Okstad                Date/Rev : July 2002 / 1.0
    !!==========================================================================

    use BushingElementTypeModule, only : BushingElementType
    use BushingElementTypeModule, only : DeallocateBushingElements
    use BushingElementTypeModule, only : NullifyBushingElement
    use TriadTypeModule         , only : TriadType, dp, GetPtrToId
    use SpringTypeModule        , only : SpringBaseType, GetPtrToId
    use DamperTypeModule        , only : DamperBaseType, GetPtrToId
    use IdTypeModule            , only : ldesc_p, initId, ReportInputError
    use inputUtilities          , only : iuGetNumberOfEntries
    use inputUtilities          , only : iuSetPosAtNextEntry
    use progressModule          , only : lterm
    use reportErrorModule       , only : reportError, debugFileOnly_p
    use reportErrorModule       , only : AllocationError

    integer                 , intent(in)  :: infp
    type(TriadType)         , intent(in)  :: triads(:)
    type(SpringBaseType)    , intent(in)  :: baseSprings(:)
    type(DamperBaseType)    , intent(in)  :: baseDampers(:)
    type(BushingElementType), pointer     :: bElems(:)
    integer                 , intent(out) :: err

    !! Local variables
    integer :: i, j, idIn, nBElem, nSpr, nDmp, nDOF, nSpokes, stat, lerr
    type(SpringBaseType), pointer :: sprBase
    type(DamperBaseType), pointer :: dmpBase

    !! Define the BUSHING_ELEMENT namelist
    character(ldesc_p) :: extDescr
    integer, parameter :: MAX_MTRIAD = 1000
    integer :: id, extId(10), nTriads, triadId(0:MAX_MTRIAD)
    integer :: springId(6), damperId(6)

    namelist /BUSHING_ELEMENT/ extDescr, id, extId, nTriads, triadId, &
         &                     springId, damperId

    !! --- Logic section ---

    nBElem = iuGetNumberOfEntries(infp,'&BUSHING_ELEMENT',err)
    if (err /= 0) goto 900
    write(lterm,*) 'Number of &BUSHING_ELEMENT =',nBElem

    call DeallocateBushingElements (bElems)
    allocate(bElems(nBElem),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadBushingElements 10')
       return
    end if

    do idIn = 1, nBElem

       call NullifyBushingElement (bElems(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&BUSHING_ELEMENT')) then
          err = err - 1
          call ReportInputError ('BUSHING_ELEMENT',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       springId=0; damperId=0; nTriads=0; triadId=0

       read(infp,nml=BUSHING_ELEMENT,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('BUSHING_ELEMENT',idIn)
          cycle
       end if

       call initId (bElems(idIn)%id,id,extId,extDescr,stat)

       nSpokes = nTriads-1
       allocate(bElems(idIn)%spokes(nSpokes),STAT=stat)
       if (stat < 0) then
          err = AllocationError('ReadBushingElement 20')
          return
       end if

       bElems(idIn)%centre => GetPtrToId(triads,triadId(0))
       if (.not. associated(bElems(idIn)%centre)) then
          err = err - 1
          call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
               msg='Invalid centre triad Id')
          cycle
       else if (bElems(idIn)%centre%nDOFs < 6) then
          err = err - 1
          call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
               msg='The centre triad has to be a 6-DOF triad')
          cycle
       end if

       nSpr = 0
       nDmp = 0
       do i = 1, 6
          if (springId(i) > 0) nSpr = nSpr + 1
          if (damperId(i) > 0) nDmp = nDmp + 1
       end do
       nDOF = bElems(idIn)%centre%nDOFs

       lerr = err
       do j = 1, nSpokes
          allocate(bElems(idIn)%spokes(j)%springs(nSpr), &
               &   bElems(idIn)%spokes(j)%dampers(nDmp), STAT=stat)
          if (stat < 0) then
             err = AllocationError('ReadBushingElement 30')
             return
          end if
          bElems(idIn)%spokes(j)%elmDof = nDOF + 1
          bElems(idIn)%spokes(j)%length = 0.0_dp
          bElems(idIn)%spokes(j)%Tlg = 0.0_dp
          bElems(idIn)%spokes(j)%triad => GetPtrToId(triads,triadId(j))
          if (associated(bElems(idIn)%spokes(j)%triad)) then
             nDOF = nDOF + bElems(idIn)%spokes(j)%triad%nDOFs
          else
             err = err - 1
             call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
                  msg='Invalid triad Id')
          end if
       end do
       bElems(idIn)%nDOFs = nDOF
       if (err < lerr) cycle

       nSpr = 0
       nDmp = 0
       do i = 1, 6
          if (springId(i) > 0) then
             nSpr = nSpr + 1

             sprBase => GetPtrToId(baseSprings,springId(i))
             if (.not. associated(sprBase)) then
                err = err - 1
                call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
                     msg='Invalid spring Id')
             else if (sprBase%dof /= -1) then
                err = err - 1
                call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
                     msg='Spring used more than once')
             else
                sprBase%dof = i
                do j = 1, nSpokes
                   call initiateSpring (bElems(idIn)%spokes(j)%springs(nSpr), &
                        &               bElems(idIn)%centre, &
                        &               bElems(idIn)%spokes(j)%triad,stat)
                   if (stat < 0) then
                      err = stat
                      return
                   end if
                end do
                sprBase%dof = -10-i
             end if

          end if
          if (damperId(i) > 0) then
             nDmp = nDmp + 1

             dmpBase => GetPtrToId(baseDampers,damperId(i))
             if (.not. associated(dmpBase)) then
                err = err - 1
                call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
                     msg='Invalid damper Id')
             else if (dmpBase%dof /= -1) then
                err = err - 1
                call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
                     msg='Damper used more than once')
             else
                dmpBase%dof = i
                do j = 1, nSpokes
                   call initiateDamper (bElems(idIn)%spokes(j)%dampers(nDmp), &
                        &               stat)
                   if (stat < 0) then
                      err = stat
                      return
                   end if
                end do
                dmpBase%dof = -10-i
             end if

          end if
       end do

       if (nSpokes < 1) then
          err = err - 1
          call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
               msg='The element must be connected to at least two triads')
       else if (nSpr < 1 .and. nDmp < 1) then
          err = err - 1
          call ReportInputError ('BUSHING_ELEM',id=bElems(idIn)%id, &
               msg='At least one spring or damper must be specified')
       end if

    end do

900 if (err < 0) call reportError (debugFileOnly_p,'ReadBushingElements')

  contains

    subroutine InitiateSpring (spring,triad1,triad2,ierr)

      use SpringTypeModule, only : nullifySpring
      use rotationModule  , only : deltaRot

      type(SpringBaseType), intent(out) :: spring
      type(TriadType)     , intent(in)  :: triad1, triad2
      integer             , intent(out) :: ierr

      real(dp) :: delta(3)

      call nullifySpring (spring)

      spring%dof                 =  sprBase%dof
      spring%stiffnessFunction   => sprBase%stiffnessFunction
      spring%forceFunction       => sprBase%forceFunction
      spring%s0                  =  sprBase%s0
      spring%s1                  =  sprBase%s1
      spring%stiffScaleEnginePos => sprBase%stiffScaleEnginePos
      spring%stiffScaleEngineNeg => sprBase%stiffScaleEngineNeg
      spring%saveVar             =  sprBase%saveVar

      if (spring%dof < 4) then
         !! Find the initial spring length
         spring%l0 = triad2%ur(spring%dof,4) - triad1%ur(spring%dof,4)
      else
         !! Find the initial spring angle
         delta = deltaRot(triad1%ur(:,1:3),triad2%ur(:,1:3))
         spring%l0 = delta(spring%dof-3)
      end if

      allocate(spring%length,STAT=ierr)
      if (ierr == 0) then
         spring%allocatedLength = .true.
      else
         ierr = AllocationError('ReadBushingElement 40')
      end if

    end subroutine InitiateSpring

    subroutine InitiateDamper (damper,ierr)

      use DamperTypeModule, only : nullifyDamper

      type(DamperBaseType), intent(out) :: damper
      integer             , intent(out) :: ierr

      call nullifyDamper (damper)

      damper%dof              =  dmpBase%dof
      damper%coeffFunction    => dmpBase%coeffFunction
      damper%forceFunction    => dmpBase%forceFunction
      damper%dmp0             =  dmpBase%dmp0
      damper%dmp1             =  dmpBase%dmp1
      damper%coeffScaleEngine => dmpBase%coeffScaleEngine
      damper%saveVar          =  dmpBase%saveVar

      allocate(damper%length,damper%velocity,STAT=ierr)
      if (ierr /= 0) ierr = AllocationError('ReadBushingElement 50')

    end subroutine InitiateDamper

  end subroutine ReadBushingElements

end module initiateBushingElmTypeModule
