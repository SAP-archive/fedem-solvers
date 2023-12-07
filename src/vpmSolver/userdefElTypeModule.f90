!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file userdefElTypeModule.f90
!> @brief User-defined element data container.

!!==============================================================================
!> @brief Module with data types representing user-defined element objects.
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the user-defined element data.

module UserdefElTypeModule

  use TriadTypeModule   , only : TriadPtrType, IdType, dp
  use FunctionTypeModule, only : EnginePtrType
  use SupElTypeModule   , only : HydroDynType

  implicit none

  !> @brief Data type representing a user-defined element.
  type UserdefElType

     type(IdType) :: id  !< General identification data

     integer :: type     !< Element type identifier
     integer :: samElNum !< Element number for SAM reference
     integer :: nTotDOFs !< Number of total DOFs (NDIM)

     type(TriadPtrType), pointer :: triads(:) !< All triads on the element
     real(dp), pointer :: TrUndeformed(:,:,:) !< Undeformed position matrices

     type(EnginePtrType), pointer :: engines(:) !< Time-dependent parameters

     real(dp) :: ePot0 !< Initial potential energy
     real(dp) :: ePot  !< Current potential energy (relative to initial)
     real(dp) :: eKin  !< Current kinetic energy
     real(dp) :: eStr  !< Current strain energy
     real(dp) :: eDmp  !< Energy loss from damping

     real(dp) :: Tlg(3,4) !< Local-to-global transformation matrix

     real(dp), pointer :: Nmat(:,:) !< Newton matrix
     real(dp), pointer :: Mmat(:,:) !< Structural mass matrix
     real(dp), pointer :: Cmat(:,:) !< Structural damping matrix
     real(dp), pointer :: KMat(:,:) !< Tangent stiffness matrix

     real(dp), pointer :: Q (:) !< External Forces
     real(dp), pointer :: FS(:) !< Forces related to the stiffness matrix
     real(dp), pointer :: FD(:) !< Forces related to the damping matrix
     real(dp), pointer :: FI(:) !< Forces related to the inertia matrix

     integer , pointer :: iwork(:) !< Integer work area
     real(dp), pointer :: rwork(:) !< Real work area

     type(HydroDynType), pointer :: hydyn !< Hydrodynamics data

  end type UserdefElType

  !> @brief Data type representing a user-defined element pointer.
  !> @details This data type is used to construct arrays of user-defined
  !> element objects where each array element is a pointer to a user-defined
  !> element object, and not the objects themselves.
  type UserdefElPtrType
     type(UserdefElType), pointer :: p !< Pointer to a user-defined element
  end type UserdefElPtrType

  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdUserdefEl
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteUserdefElType
  end interface

  private :: GetPtrToIdUserdefEl, WriteUserdefElType


contains

  !!============================================================================
  !> @brief Returns pointer to (first) user-defined element with specified ID.
  !>
  !> @param[in] array Array of userdefeltypemodule::userdefeltype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[in] userId If @e .true., search for a user ID instead
  !>
  !> @details If the user-defined element is not found, NULL is returned.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  function GetPtrToIdUserdefEl (array,id,userId) result(ptr)
    type(UserdefElType), pointer :: ptr

    integer            , intent(in)           :: id
    type(UserdefElType), intent(in), target   :: array(:)
    logical            , intent(in), optional :: userId

    !! Local variables
    integer :: i
    logical :: searchUserId

    !! --- Logic section ---

    if (present(userId)) then
       searchUserId = userId
    else
       searchUserId = .false.
    end if

    do i = 1, size(array)
       if (searchUserId) then
          if (array(i)%id%userId /= id) cycle
       else
          if (array(i)%id%baseId /= id) cycle
       end if
       ptr => array(i)
       return
    end do

    nullify(ptr)
    if (searchUserId) then
       write(*,*) '*** GetPtrToIdUserdefEl returned nullified, userId =',id
    else
       write(*,*) '*** GetPtrToIdUserdefEl returned nullified, baseId =',id
    end if

  end function GetPtrToIdUserdefEl


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] elm The userdefeltypemodule::userdefeltype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine WriteUserdefElType (elm,io,complexity)

    use IdTypeModule, only : writeId

    type(UserdefElType), intent(in) :: elm
    integer            , intent(in) :: io
    integer,   optional, intent(in) :: complexity

    !! Local variables
    integer :: i, j ,n

    !! --- Logic section ---

    write(io,'(A)') 'Userdefined element','{'
    call writeId (elm%id,io)

    write(io,*) 'type        =', elm%type
    write(io,*) 'samElNum    =', elm%samElNum
    write(io,*) 'nTotDOFs    =', elm%nTotDOFs

    if (associated(elm%hydyn)) then
       write(io,*) 'bodyIndex   =', elm%hydyn%bodyIndex
       write(io,3) 'Morison     =', elm%hydyn%Morison
    end if

    write(io,2) 'Tlg         =', (elm%Tlg(i,:),i=1,3)

    if (present(complexity)) then
       if (complexity >= 1) then

          n = size(elm%triads)
          write(io,9) 'triads(id)  =', (elm%triads(j)%p%id%baseId,j=1,n)
          write(io,9) 'dofStart    =', (elm%triads(j)%firstDOF,j=1,n)

          n = size(elm%TrUndeformed,3)
          do i = 1, n
             write(io,'(A,i3,A/4ES15.6E3/4ES15.6E3/4ES15.6E3)')&
                  & ' TrUndeformed, node',i,' =', elm%TrUndeformed(1,:,i), &
                  &         elm%TrUndeformed(2,:,i), elm%TrUndeformed(3,:,i)
          end do

          n = size(elm%engines)
          write(io,9) 'engines(id) =', (elm%engines(j)%p%id%baseId,j=1,n)

       end if
       if (complexity >= 2) then

          if (associated(elm%hydyn)) then
             write(io,3) 'wn            =', elm%hydyn%wn
             write(io,3) 'As            =', elm%hydyn%As
             write(io,3) 'Vb            =', elm%hydyn%Vb
             write(io,3) 'C0s           =', elm%hydyn%C0s
             write(io,3) 'C0b           =', elm%hydyn%C0b
             write(io,3) 'B             =', elm%hydyn%B
             write(io,3) 'M             =', elm%hydyn%M
             write(io,3) 'D             =', elm%hydyn%D
             write(io,3) 'S             =', elm%hydyn%S
          end if

          write(io,3) 'ePot0       =', elm%ePot0
          write(io,3) 'ePot        =', elm%ePot
          write(io,3) 'eKin        =', elm%eKin
          write(io,3) 'eStr        =', elm%eStr
          write(io,3) 'eDmp        =', elm%eDmp

       end if
       if (complexity >= 3 .or. complexity < 0) then

          write(io,*)
          write(io,*) 'size(triads)        =', size(elm%triads)
          write(io,*) 'size(TrUndeformed)  =', size(elm%TrUndeformed,1) &
               &                             , size(elm%TrUndeformed,2) &
               &                             , size(elm%TrUndeformed,3)
          if (associated(elm%Nmat)) then
             write(io,*) 'size(Nmat)          =', size(elm%Nmat,1) &
                  &                             , size(elm%Nmat,2)
          end if
          if (associated(elm%Mmat)) then
             write(io,*) 'size(Mmat)          =', size(elm%Mmat,1) &
                  &                             , size(elm%Mmat,2)
          end if
          if (associated(elm%Cmat)) then
             write(io,*) 'size(Cmat)          =', size(elm%Cmat,1) &
                  &                             , size(elm%Cmat,2)
          end if
          if (associated(elm%KMat)) then
             write(io,*) 'size(KMat)          =', size(elm%KMat,1) &
                  &                             , size(elm%KMat,2)
          end if
          if (associated(elm%Q)) then
             write(io,*) 'size(Q)             =', size(elm%Q)
          end if
          if (associated(elm%FS)) then
             write(io,*) 'size(FS)            =', size(elm%FS)
          end if
          if (associated(elm%FD)) then
             write(io,*) 'size(FD)            =', size(elm%FD)
          end if
          if (associated(elm%FI)) then
             write(io,*) 'size(FI)            =', size(elm%FI)
          end if

          if (associated(elm%iwork)) then
             write(io,*) 'size(iwork)         =', size(elm%iwork)
          end if
          if (associated(elm%rwork)) then
             write(io,*) 'size(rwork)         =', size(elm%rwork)
          end if

       end if
    end if

    write(io,'(A)') '}'

2   format(1X,A,3(/4ES15.6E3))
3   format(1X,A,8ES15.6E3)
9   format(1X,A,10I9/(15X,10I9))

  end subroutine WriteUserdefElType


  !!============================================================================
  !> @brief Initializes a user-defined element object.
  !>
  !> @param elm The userdefeltypemodule::userdefeltype object to initialize
  !> @param[in] deallocating If .true., the pointers are nullified
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Mar 2014

  subroutine NullifyUserdefEl (elm,deallocating)

    use IdTypeModule, only : nullifyId

    type(UserdefElType), intent(out) :: elm
    logical, optional  , intent(in)  :: deallocating

    !! --- Logic section ---

    call nullifyId (elm%id)

    elm%type     = 0
    elm%samElNum = 0
    elm%nTotDOFs = 0

    if (present(deallocating)) then
       nullify(elm%triads)
       nullify(elm%TrUndeformed)
       nullify(elm%engines)
    else
       allocate(elm%triads(0))
       allocate(elm%TrUndeformed(0,0,0))
       allocate(elm%engines(0))
    end if

    elm%ePot0 = 0.0_dp
    elm%ePot = 0.0_dp
    elm%eKin = 0.0_dp
    elm%eStr = 0.0_dp
    elm%eStr = 0.0_dp
    elm%eDmp = 0.0_dp
    elm%Tlg  = 0.0_dp

    nullify(elm%Nmat)
    nullify(elm%Mmat)
    nullify(elm%Cmat)
    nullify(elm%KMat)
    nullify(elm%Q)
    nullify(elm%FS)
    nullify(elm%FD)
    nullify(elm%FI)
    nullify(elm%iwork)
    nullify(elm%rwork)
    nullify(elm%hydyn)

  end subroutine NullifyUserdefEl


  !!============================================================================
  !> @brief Deallocates a user-defined element object.
  !>
  !> @param elm The userdefeltypemodule::userdefeltype object to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateUserdefEl (elm)

    use IdTypeModule, only : deallocateId

    type(UserdefElType), intent(inout) :: elm

    !! --- Logic section ---

    call deallocateId (elm%id)

    if (associated(elm%triads))       deallocate(elm%triads)
    if (associated(elm%TrUndeformed)) deallocate(elm%TrUndeformed)
    if (associated(elm%engines))      deallocate(elm%engines)
    if (associated(elm%Nmat))         deallocate(elm%Nmat)
    if (associated(elm%Mmat))         deallocate(elm%Mmat)
    if (associated(elm%Cmat))         deallocate(elm%Cmat)
    if (associated(elm%Kmat))         deallocate(elm%KMat)
    if (associated(elm%Q))            deallocate(elm%Q)
    if (associated(elm%FS))           deallocate(elm%FS)
    if (associated(elm%FD))           deallocate(elm%FD)
    if (associated(elm%FI))           deallocate(elm%FI)
    if (associated(elm%iwork))        deallocate(elm%iwork)
    if (associated(elm%rwork))        deallocate(elm%rwork)
    if (associated(elm%hydyn))        deallocate(elm%hydyn)

    call nullifyUserdefEl (elm,.true.)

  end subroutine DeallocateUserdefEl


  !!============================================================================
  !> @brief Deallocates all user-defined element object.
  !>
  !> @param elms Array of all userdefeltypemodule::userdefeltype objects
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateUdElms (elms)

    type(UserdefElType), pointer :: elms(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(elms)
       call DeallocateUserdefEl (elms(i))
    end do

    deallocate(elms)
    nullify(elms)

  end subroutine deallocateUdElms

end module UserdefElTypeModule
