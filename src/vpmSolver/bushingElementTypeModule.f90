!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module BushingElementTypeModule

  use TriadTypeModule , only : TriadType, IdType, dp
  use SpringTypeModule, only : SpringBaseType
  use DamperTypeModule, only : DamperBaseType

  implicit none

  type SpokeType

     type(TriadType)     , pointer :: triad
     type(SpringBaseType), pointer :: springs(:)
     type(DamperBaseType), pointer :: dampers(:)

     integer  :: elmDof   !! The first element DOF for this spoke
     real(dp) :: length   !! Distance between the centre and spoke triad
     real(dp) :: Tlg(3,3) !! Local-to-global transformation matrix

  end type SpokeType


  type BushingElementType

     type(IdType) :: id  !! General identification data

     integer :: samElNum !! Element number for SAM reference (MPMNPC)
     integer :: nDOFs    !! = centre%nDOFs + sum(spokes(j)%triad)

     type(TriadType), pointer :: centre
     type(SpokeType), pointer :: spokes(:)

  end type BushingElementType

  interface WriteObject
     module procedure WriteBushingElementType
  end interface


contains

  subroutine WriteBushingElementType (bElem,io,complexity)

    !!==========================================================================
    !! Standard routine for writing an object to io.
    !!
    !! Programmer : Knut Morten okstad
    !! date/rev   : 11 Jul 2002/1.0
    !!==========================================================================

    use IdTypeModule    , only : writeId
    use SpringTypeModule, only : writeObject
    use DamperTypeModule, only : writeObject

    type(BushingElementType), intent(in) :: bElem
    integer                 , intent(in) :: io
    integer, optional       , intent(in) :: complexity

    !! Local variables
    integer :: i, j, n

    !! --- Logic section ---

    write(io,'(A)') 'BushingElement','{'
    call writeId (bElem%id,io)

    write(io,*) 'samElNum    =', bElem%samElnum
    write(io,*) 'nDOFs       =', bElem%nDOFs

    if (associated(bElem%centre)) then
       write(io,*) 'centre(id)  = ', bElem%centre%id%baseId
    end if

    if (associated(bElem%spokes)) then
       n = size(bElem%spokes)
       write(io,*) 'triads(id)  = ',(bElem%spokes(i)%triad%id%baseId,i=1,n)
       write(io,*) 'elmDofs     = ',(bElem%spokes(i)%elmDof,i=1,n)
    else
       n = 0
    end if

    if (present(complexity) .and. n > 0) then
       if (complexity >= 1) then
          do i = 1, n

             write(io,*)
             write(io,*) 'Spoke number',i
             if (complexity >= 2) then
                write(io,*) 'length      =', bElem%spokes(i)%length
                write(io,6) 'Tlg         =',(bElem%spokes(i)%Tlg(j,:),j=1,3)
             end if

             if (associated(bElem%spokes(i)%springs)) then
                do j = 1, size(bElem%spokes(i)%springs)
                   call writeObject (bElem%spokes(i)%springs(j),io,complexity)
                end do
             end if

             if (associated(bElem%spokes(i)%dampers)) then
                do j = 1, size(bElem%spokes(i)%dampers)
                   call writeObject (bElem%spokes(i)%dampers(j),io,complexity)
                end do
             end if

          end do
       end if
    end if

    write(io,'(A)') '}'

6   format(1X,A13,3F12.8/(14X,3F12.8))

  end subroutine WriteBushingElementType


  subroutine NullifyBushingElement (bElem,deallocating)

    !!==========================================================================
    !! Initialize the BushingElementType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 11 Jul 2002/1.0
    !!==========================================================================

    use IdTypeModule, only : nullifyId

    type(BushingElementType), intent(out) :: bElem
    logical, optional       , intent(in)  :: deallocating

    !! --- Logic section ---

    call nullifyId (bElem%id)

    bElem%samElNum = 0
    bElem%nDOFs = 0

    nullify(bElem%centre)
    if (present(deallocating)) then
       nullify(bElem%spokes)
    else
       !! Allocate initially to zero length, such that the size() operator works
       allocate(bElem%spokes(0))
    end if

  end subroutine NullifyBushingElement


  subroutine DeallocateBushingElement (bElem)

    !!==========================================================================
    !! Deallocate the BushingElementType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 23 Jan 2017/1.0
    !!==========================================================================

    use IdTypeModule, only : deallocateId

    type(BushingElementType), intent(inout) :: bElem

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    call deallocateId (bElem%id)

    do i = 1, size(bElem%spokes)
       do j = 1, size(belem%spokes(i)%springs)
          deallocate(belem%spokes(i)%springs(j)%length)
       end do
       do j = 1, size(belem%spokes(i)%dampers)
          deallocate(belem%spokes(i)%dampers(j)%length)
          deallocate(belem%spokes(i)%dampers(j)%velocity)
       end do
       deallocate(belem%spokes(i)%springs,belem%spokes(i)%dampers)
    end do
    deallocate(bElem%spokes)

    call nullifyBushingElement (bElem,.true.)

  end subroutine DeallocateBushingElement


  subroutine deallocateBushingElements (bElms)
    type(BushingElementType), pointer :: bElms(:)
    integer :: i
    do i = 1, size(bElms)
       call DeallocateBushingElement (bElms(i))
    end do
    deallocate(bElms)
    nullify(bElms)
  end subroutine deallocateBushingElements

end module BushingElementTypeModule
