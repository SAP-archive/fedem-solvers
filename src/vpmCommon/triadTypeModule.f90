!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file triadTypeModule.f90
!> @brief Triad object data container.

!!==============================================================================
!> @brief Module with data types representing triad objects.
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the triad data.

module TriadTypeModule

  use KindModule  , only : dp
  use IdTypeModule, only : IdType

  implicit none

  !> @brief Data type representing a triad object.
  type TriadType

     type(IdType) :: id !< General identification data

     integer :: samNodNum !< Node number for SAM reference (madof)
     integer :: nDOFs     !< 0 (grounded), 3 or 6 DOFs is allowed
     integer :: inFluid   !< Flag/index to fluid data for hydrodynamics
     logical :: dependent !< If @e .true., this is a dependent joint triad
     logical :: isFEnode  !< If @e .true., this triad is used by finite elements
     logical :: isSEnode  !< If @e .true., this triad is used by a superelement

     integer , pointer :: sysDOF      !< System DOF no. for the 1st DOF in triad
     integer , pointer :: eqNumber(:) !< Global equation numbers for all DOFs

     !> @brief Boundary condition status flag for each DOF.
     !> @details 0 is fixed, 1 is free, 2 is fixed during initial static
     !> equilibrium iterations and eigenvalue calculations, and free otherwise.
     integer , pointer :: BC(:)

     !> @brief Global-to-system direction transformation matrix.
     !> @details If associated, the system directions of the triad are defined
     !> by this 3x3 matrix. If nullified, the system directions coinside
     !> with the global axis directions.
     real(dp), pointer :: sysDirInG(:,:)

     real(dp), pointer :: dragPar(:,:) !< Local drag coefficients

     real(dp) :: ur(3,4)     !< Position matrix at current time step
     real(dp) :: urPrev(3,4) !< Position matrix at previous time step

     real(dp), pointer :: ur_def(:)    !< Deformational displacements
     real(dp), pointer :: urd(:)       !< Velocities in global directions
     real(dp), pointer :: urdd(:)      !< Accelerations in global directions
     real(dp), pointer :: scoord       !< Running coordinate along 1D beams
     real(dp), pointer :: kappa(:)     !< Curvature along 1D beams
     real(dp), pointer :: Moment(:)    !< Bending moments along 1D beams
     real(dp), pointer :: nodeForce(:) !< Nodal forces in system directions
     real(dp), pointer :: supForce(:)  !< Superelement forces, system directions
     real(dp), pointer :: aeroForce(:) !< Aerodynamic forces, local directions
     real(dp), pointer :: force(:,:)   !< Force contributions (s, d, i and ex)

     logical :: savePos    !< Flag indicating whether position should be saved
     logical :: saveVar(9) !< Flags indicating which variables should be saved

  end type TriadType


  !> @brief Data type representing a triad pointer.
  !> @details This data type is used to construct arrays of triads where each
  !> element is a pointer to a triad object, and not the objects themselves.
  !> Mostly used for representing element connectivities.
  type TriadPtrType
     type(TriadType), pointer :: p !< Pointer to a triad object
     integer :: firstDOF !< Index to first DOF in element using this triad
  end type TriadPtrType


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdTriad
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteTriadType
  end interface

  !> @brief Transforms a triad quantity from global to system directions.
  interface transGlobToSys
     module procedure transVecGlobToSys
     module procedure transMatGlobToSys
  end interface

  !> @brief Transforms a triad quantity from system to global directions.
  interface transSysToGlob
     module procedure transVecSysToGlob
  end interface

  !> @brief Transforms a triad quantity from triad to system directions.
  interface transTriadToSys
     module procedure transVecTriadToSys
     module procedure transMatTriadToSys
  end interface

  !> @brief Transforms a triad quantity from system to triad directions.
  interface transSysToTriad
     module procedure transVecSysToTriad
  end interface

  !> @brief Transforms a triad quantity from triad to global directions.
  interface transTriadToGlob
     module procedure transVecTriadToGlob
     module procedure transMatTriadToGlob
  end interface

  !> @brief Transforms a triad quantity from global to triad directions.
  interface transGlobToTriad
     module procedure transVecGlobToTriad
  end interface

  private :: GetPtrToIdTriad, WriteTriadType
  private :: transVecGlobToSys, transMatGlobToSys, transVecSysToGlob
  private :: transVecTriadToSys, transMatTriadToSys, transVecSysToTriad
  private :: transVecTriadToGlob, transMatTriadToGlob, transVecGlobToTriad


contains

  !!============================================================================
  !> @brief Returns pointer to (first) triad with specified ID.
  !>
  !> @param[in] array Array of triadtypemodule::triadtype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[out] index The array index of the found object
  !> @param[in] userId If @e .true., search for a user ID instead
  !>
  !> @details If the triad is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen / Knut Morten Okstad
  !>
  !> @date 1 Dec 2009

  function GetPtrToIdTriad (array,id,index,userId) result(ptr)

    integer        , intent(in)            :: id
    type(TriadType), intent(in) , target   :: array(:)
    integer        , intent(out), optional :: index
    logical        , intent(in) , optional :: userId
    type(TriadType), pointer               :: ptr

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
       if (present(index)) index = i
       return
    end do

    nullify(ptr)
    if (searchUserId) then
       write(*,*) '*** GetPtrToIdTriad returned nullified, userId =',id
    else
       write(*,*) '*** GetPtrToIdTriad returned nullified, baseId =',id
    end if

  end function GetPtrToIdTriad


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 27 Sep 1998

  subroutine WriteTriadType (triad,io,complexity)

    use IdTypeModule, only : writeId

    type(TriadType) , intent(in) :: triad
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'Triad','{'
    call writeId (triad%id,io)

    write(io,*) 'samNodNum   =', triad%samNodNum
    write(io,*) 'nDOFs       =', triad%nDOFs
    write(io,*) 'inFluid     =', triad%inFluid
    write(io,*) 'dependent?  =', triad%dependent
    write(io,*) 'isFEnode?   =', triad%isFEnode
    write(io,*) 'isSEnode?   =', triad%isSEnode

    if (associated(triad%sysDOF))   write(io,*) 'sysDOF      =', triad%sysDOF
    if (associated(triad%eqNumber)) write(io,*) 'eqNumber    =', triad%eqNumber
    if (associated(triad%BC))       write(io,*) 'BC          =', triad%BC

    if (associated(triad%sysDirInG)) then
       write(io,600) ' sysDir      =', (triad%sysDirInG(i,:),i=1,3)
    end if
    if (associated(triad%dragPar)) then
       write(io,600) ' dragPar     =', triad%dragPar
    end if

    write(io,*) 'savePos     =', triad%savePos
    if (triad%nDOFs > 0) write(io,*) 'saveVar     =', triad%saveVar

    if (associated(triad%scoord)) then
       write(io,600) ' scoord      =', triad%scoord
    end if
    if (associated(triad%kappa)) then
       write(io,602) ' kappa       =', triad%kappa
    end if
    if (associated(triad%Moment)) then
       write(io,602) ' Moment      =', triad%Moment
    end if
    write(io,601) ' ur          =', (triad%ur(i,:),i=1,3)

    if (present(complexity) .and. triad%nDOFs > 0) then
       if (complexity >= 2) then

          write(io,601) ' urPrev       =', (triad%urPrev(i,:),i=1,3)
          if (associated(triad%ur_def)) then
             write(io,602) ' ur_def      =', triad%ur_def
          end if
          if (associated(triad%urd)) then
             write(io,602) ' urd         =', triad%urd
          end if
          if (associated(triad%urdd)) then
             write(io,602) ' urdd        =', triad%urdd
          end if
          if (associated(triad%nodeForce)) then
             write(io,602) ' nodeForce   =', triad%nodeForce
          end if
          if (associated(triad%supForce)) then
             write(io,602) ' supForce    =', triad%supForce
          end if
          if (associated(triad%aeroForce)) then
             write(io,602) ' aeroForce   =', triad%aeroForce
          end if

       end if
       if (complexity >= 3 .or. complexity < 0) then

          write(io,*)
          if (associated(triad%eqNumber)) then
             write(io,*) 'size(eqNumber)  =', size(triad%eqNumber)
          end if
          if (associated(triad%BC)) then
             write(io,*) 'size(BC)        =', size(triad%BC)
          end if
          if (associated(triad%sysDirInG)) then
             write(io,*) 'size(sysDirInG) =', size(triad%sysDirInG,1) &
                  &                         , size(triad%sysDirInG,2)
          end if
          if (associated(triad%dragPar)) then
             write(io,*) 'size(dragPar)   =', size(triad%dragPar,1) &
                  &                         , size(triad%dragPar,2)
          end if
          write(io,*) 'size(ur)        =', size(triad%ur,1), size(triad%ur,2)
          write(io,*) 'size(urPrev)    =', size(triad%urPrev,1) &
               &                         , size(triad%urPrev,2)
          if (associated(triad%ur_def)) then
             write(io,*) 'size(ur_def)    =', size(triad%ur_def)
          end if
          if (associated(triad%urd)) then
             write(io,*) 'size(urd)       =', size(triad%urd)
          end if
          if (associated(triad%urdd)) then
             write(io,*) 'size(urdd)      =', size(triad%urdd)
          end if
          if (associated(triad%nodeForce)) then
             write(io,*) 'size(nodeForce) =', size(triad%nodeForce)
          end if
          if (associated(triad%supForce)) then
             write(io,*) 'size(supForce)  =', size(triad%supForce)
          end if
          if (associated(triad%aeroForce)) then
             write(io,*) 'size(aeroForce) =', size(triad%aeroForce)
          end if
          if (associated(triad%force)) then
             write(io,*) 'size(force)     =', size(triad%force,1) &
                  &                         , size(triad%force,2)
          end if

       end if
    end if

    write(io,'(A)') '}'

600 format(A14,3ES15.6E3/(14X,3ES15.6E3))
601 format(A14,4ES15.6E3/(14X,4ES15.6E3))
602 format(A14,9ES15.6E3)

  end subroutine WriteTriadType


  !!============================================================================
  !> @brief Initializes a triad object.
  !>
  !> @param triad The triadtypemodule::triadtype object to initialize
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 13 Oct 2000

  subroutine NullifyTriad (triad)

    use IdTypeModule, only : nullifyId

    type(TriadType), intent(out) :: triad

    !! --- Logic section ---

    call nullifyId (triad%id)

    triad%samNodNum = 0
    triad%nDOFs = 0
    triad%inFluid = 0
    triad%dependent = .false.
    triad%isFEnode = .false.
    triad%isSEnode = .false.
    triad%ur = 0.0_dp
    triad%urPrev = 0.0_dp
    triad%savePos = .false.
    triad%saveVar = .false.

    nullify(triad%sysDOF)
    nullify(triad%eqNumber)
    nullify(triad%BC)
    nullify(triad%sysDirInG)
    nullify(triad%dragPar)
    nullify(triad%ur_def)
    nullify(triad%urd)
    nullify(triad%urdd)
    nullify(triad%scoord)
    nullify(triad%kappa)
    nullify(triad%Moment)
    nullify(triad%nodeForce)
    nullify(triad%supForce)
    nullify(triad%aeroForce)
    nullify(triad%force)

  end subroutine NullifyTriad


  !!============================================================================
  !> @brief Deallocates a triad object.
  !>
  !> @param triad The triadtypemodule::triadtype object to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 23 Jan 2017

  subroutine DeallocateTriad (triad)

    use IdTypeModule, only : deallocateId

    type(TriadType), target, intent(inout) :: triad

    !! --- Logic section ---

    call deallocateId (triad%id)

    if (associated(triad%BC)) deallocate(triad%BC)
    if (associated(triad%sysDirInG,triad%urPrev(:,1:3))) then
       nullify(triad%sysDirInG)
    else if (associated(triad%sysDirInG,triad%ur(:,1:3))) then
       nullify(triad%sysDirInG)
    else if (associated(triad%sysDirInG)) then
       deallocate(triad%sysDirInG)
    end if
    if (associated(triad%dragPar)) deallocate(triad%dragPar)
    if (associated(triad%urd))     deallocate(triad%urd)
    if (associated(triad%urdd))    deallocate(triad%urdd)
    if (associated(triad%scoord))  deallocate(triad%scoord)
    if (associated(triad%kappa))   deallocate(triad%kappa)
    if (associated(triad%Moment))  deallocate(triad%Moment)
    if (associated(triad%supForce)) then
       if (associated(triad%nodeForce,triad%supForce)) then
          nullify(triad%nodeForce)
       end if
       deallocate(triad%supForce)
    end if
    if (associated(triad%nodeForce)) deallocate(triad%nodeForce)
    if (associated(triad%aeroForce)) then
       deallocate(triad%ur_def)
    else if (associated(triad%ur_def) .and. .not.triad%isSEnode) then
       deallocate(triad%ur_def)
    end if
    if (associated(triad%force)) deallocate(triad%force)

    call nullifyTriad (triad)

  end subroutine DeallocateTriad


  !!============================================================================
  !> @brief Deallocates an array of triad objects.
  !>
  !> @param triads The triadtypemodule::triadtype objects to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 23 Jan 2017

  subroutine deallocateTriads (triads)

    type(TriadType), pointer :: triads(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(triads)
       call DeallocateTriad (triads(i))
    end do

    deallocate(triads)
    nullify(triads)

  end subroutine deallocateTriads


  !!============================================================================
  !> @brief Initializes all triad forces to zero.
  !>
  !> @param[out] triads The triadtypemodule::triadtype objects to initialize
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Nov 2008

  subroutine clearTriadForces (triads)

    type(TriadType), intent(out) :: triads(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(triads)
       if (associated(triads(i)%nodeForce)) then
          triads(i)%nodeForce = 0.0_dp
       end if
       if (associated(triads(i)%supForce)) then
          triads(i)%supForce = 0.0_dp
       end if
       if (associated(triads(i)%force)) then
          triads(i)%force = 0.0_dp
       end if
    end do

  end subroutine clearTriadForces


  !!============================================================================
  !> @brief Allocates force vectors for frs-output of superelement forces.
  !>
  !> @param triad The triadtypemodule::triadtype object to allocate for
  !> @param ierr Error flag, decremented on allocation error(s)
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Oct 2009

  subroutine AllocateTriadForces (triad,ierr)

    use reportErrorModule, only : allocationError

    type(TriadType), intent(inout) :: triad
    integer        , intent(inout) :: ierr

    !! Local variables
    integer :: stat

    !! --- Logic section ---

    triad%isSEnode = .true.
    if (triad%nDOFs < 1) return

    if (triad%saveVar(3) .or. triad%saveVar(6)) then
       if (.not. associated(triad%supForce)) then
          allocate(triad%supForce(triad%nDOFs),STAT=stat)
          if (stat /= 0) then
             ierr = ierr + AllocationError('AllocateTriadForces 1')
             return
          end if
          triad%supForce = 0.0_dp
       end if
       if (.not. associated(triad%nodeForce)) then
          triad%nodeForce => triad%supForce
       end if
    end if

    if (triad%saveVar(8)) then
       if (.not. associated(triad%force)) then
          allocate(triad%force(triad%nDOFs,4),STAT=stat)
          if (stat /= 0) then
             ierr = ierr + AllocationError('AllocateTriadForces 2')
             return
          end if
          triad%force = 0.0_dp
       end if
    end if

  end subroutine AllocateTriadForces


  !!============================================================================
  !> @brief Allocates the total nodal force vector associated with a triad.
  !>
  !> @param triad The triadtypemodule::triadtype object to allocate for
  !> @param ierr Error flag, decremented on allocation error(s)
  !>
  !> @details Total nodal force vectors are allocated only if the triad is
  !> connected to other force-giving objects in addition to (or instead of)
  !> a superelement. If it is only connected to a superelement,
  !> the @a nodeForce is identical to the @a supForce vector.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Oct 2009

  subroutine AllocateNodeForce (triad,ierr)

    use reportErrorModule, only : allocationError

    type(TriadType), intent(inout) :: triad
    integer        , intent(inout) :: ierr

    !! Local variables
    integer :: stat

    !! --- Logic section ---

    !! Do nothing for grounded triads
    if (triad%nDOFs < 1) return

    !! Do nothing if nodeForce already points to something else
    if (associated(triad%nodeForce)) then
       if (.not. associated(triad%nodeForce,triad%supForce)) return
    end if

    !! This triad is connected to both a superelement, and another force-giving
    !! mechanism object. Allocate a separate array for the total forces in the
    !! triad to distinguish them from the superelement (or section) force.
    allocate(triad%nodeForce(triad%nDOFs),STAT=stat)
    if (stat /= 0) ierr = ierr + allocationError('AllocateNodeForce')

  end subroutine AllocateNodeForce


  !!============================================================================
  !> @brief Updates the total nodal force vector associated with a triad.
  !>
  !> @param triad The triadtypemodule::triadtype object to update
  !> @param[in] force Nodal force to add to the triad
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Oct 2009

  subroutine UpdateNodeForce (triad,force)

    type(TriadType), intent(inout) :: triad
    real(dp),target, intent(in)    :: force(:)

    !! Local variables
    integer :: n

    !! --- Logic section ---

    if (.not. associated(triad%nodeForce)) return
    if (associated(triad%nodeForce,force)) return

    n = min(triad%nDOFs,size(force))
    if (n < 1) return

    triad%nodeForce(1:n) = triad%nodeForce(1:n) + force(1:n)

  end subroutine UpdateNodeForce


  !!============================================================================
  !> @brief Transforms a system vector to global directions.
  !>
  !> @param[in] triads All triad objects in the model
  !> @param inc System displacement vector
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 16 Nov 2001

  subroutine TransSysVecToGlobal (triads,inc)

    type(TriadType), intent(in)    :: triads(:)
    real(dp)       , intent(inout) :: inc(:)

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    do i = 1, size(triads)
       if (triads(i)%nDOFs == 0) cycle ! fixed triad

       j = triads(i)%sysDOF
       k = j + triads(i)%nDOFs-1

       call transVecSysToGlob (triads(i),inc(j:k))

    end do

  end subroutine TransSysVecToGlobal


  !!============================================================================
  !> @brief Updates the position matrices for all triads.
  !>
  !> @param triads All triad objects in the model
  !> @param[in] inc Incremental global displacement vector
  !> @param[in] useTotalInc If @e .true., update w.r.t. previous configuration
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 16 Nov 2001

  subroutine IncTriadsPos (triads,inc,useTotalInc)

    use rotationModule, only : vec_to_mat

    type(TriadType)  , intent(inout) :: triads(:)
    real(dp)         , intent(in)    :: inc(:)
    logical, optional, intent(in)    :: useTotalInc

    !! Local variables
    logical  :: fromPrevious
    integer  :: i, j
    real(dp) :: R(3,3) ! Rotation matrix for small rotational increments

    !! --- Logic section ---

    if (present(useTotalInc)) then
       fromPrevious = useTotalInc
    else
       fromPrevious = .false.
    end if

    do i = 1, size(triads)
       if (triads(i)%nDOFs == 0) cycle ! fixed triad

       j = triads(i)%sysDOF

       !! Add translatory part of the increment to the triad position
       if (fromPrevious) then
          triads(i)%ur(:,4) = triads(i)%urPrev(:,4) + inc(j:j+2)
       else
          triads(i)%ur(:,4) = triads(i)%ur(:,4) + inc(j:j+2)
       end if

       !! Add rotational part of the increment to the triad position
       if (triads(i)%nDOFs == 6) then
          call vec_to_mat (inc(j+3:j+5),R)
          if (fromPrevious) then
             triads(i)%ur(:,1:3) = matmul(R,triads(i)%urPrev(:,1:3))
          else
             triads(i)%ur(:,1:3) = matmul(R,triads(i)%ur(:,1:3))
          end if
       end if

    end do

  end subroutine IncTriadsPos


  !!============================================================================
  !> @brief Returns the current incremental displacement in a triad.
  !>
  !> @param triad The triadtypemodule::triadtype object to get displacement for
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Apr 2008

  function GetTriadIncrement (triad) result(inc)

    use rotationModule, only : deltaRot

    type(TriadType), intent(in) :: triad
    real(dp)                    :: inc(triad%nDOFs)

    !! --- Logic section ---

    if (triad%nDOFs >= 3) then
       inc(1:3) = triad%ur(:,4) - triad%urPrev(:,4)
       if (triad%nDOFs >= 6) then
          inc(4:6) = deltaRot(triad%urPrev(:,1:3),triad%ur(:,1:3))
       end if
    end if

  end function GetTriadIncrement


  !!============================================================================
  !> @brief Sets global velocity and acceleration for all triads.
  !>
  !> @param triads All triadtypemodule:triadtype objects in the model
  !> @param[in] velGlobal System velocity vector
  !> @param[in] accGlobal System acceleration vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Jun 2002

  subroutine SetTriadsVelAcc (triads,velGlobal,accGlobal)

    type(TriadType), intent(inout) :: triads(:)
    real(dp)       , intent(in)    :: velGlobal(:), accGlobal(:)

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    do i = 1, size(triads)
       if (triads(i)%nDOFs == 0) cycle ! fixed triad

       j = triads(i)%sysDOF
       k = j + triads(i)%nDOFs-1

       triads(i)%urd  = velGlobal(j:k)
       triads(i)%urdd = accGlobal(j:k)

    end do

  end subroutine SetTriadsVelAcc


  !!============================================================================
  !> @brief Gets current velocity and acceleration from all triads.
  !>
  !> @param[in] triads All triadtypemodule:triadtype objects in the model
  !> @param[out] velGlobal System velocity vector
  !> @param[out] accGlobal System acceleration vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Oct 2002

  subroutine GetTriadsVelAcc (triads,velGlobal,accGlobal)

    type(TriadType), intent(in)  :: triads(:)
    real(dp)       , intent(out) :: velGlobal(:), accGlobal(:)

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    do i = 1, size(triads)
       if (triads(i)%nDOFs == 0) cycle ! fixed triad

       j = triads(i)%sysDOF
       k = j + triads(i)%nDOFs-1

       velGlobal(j:k) = triads(i)%urd
       accGlobal(j:k) = triads(i)%urdd

    end do

  end subroutine GetTriadsVelAcc


  !!============================================================================
  !> @brief Checks if a triad has a local coordinate system attached.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to check
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Oct 2000

  function hasLocalDirections (triad) result(answer)

    type(TriadType), intent(in) :: triad
    logical                     :: answer

    !! --- Logic section ---

    answer = associated(triad%sysDirInG)

  end function hasLocalDirections


  !!============================================================================
  !> @brief Transforms a nodal vector from global to system directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param v The nodal vector to transform
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Jul 2002

  subroutine transVecGlobToSys (triad,v)

    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: v(:)

    !! --- Logic section ---

    if (.not. associated(triad%sysDirInG)) return

    if (size(v) >= 3) v(1:3) = matmul(v(1:3),triad%sysDirInG)
    if (size(v) >= 6) v(4:6) = matmul(v(4:6),triad%sysDirInG)

  end subroutine transVecGlobToSys


  !!============================================================================
  !> @brief Returns a component of a vector transformed to system directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param[in] ldof Local index of the component to return
  !> @param[in] v The nodal vector to transform
  !>
  !> @details The @a ldof'th component of the vector @a v is transformed
  !> form global to system directions on given triad.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jun 2015

  function getSysValue (triad,ldof,v) result(lval)

    type(TriadType), intent(in) :: triad
    integer        , intent(in) :: ldof
    real(dp)       , intent(in) :: v(:)
    real(dp)                    :: lval

    !! --- Logic section ---

    if (ldof < 1 .or. ldof > size(v)) then
       lval = 0.0_dp
    else if (.not. associated(triad%sysDirInG)) then
       lval = v(ldof)
    else if (ldof <= 3) then
       lval = dot_product(v(1:3),triad%sysDirInG(:,ldof))
    else if (ldof <= 6) then
       lval = dot_product(v(4:6),triad%sysDirInG(:,ldof-3))
    else
       lval = v(ldof)
    end if

  end function getSysValue


  !!============================================================================
  !> @brief Transforms a triad vector from system to global directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param v The nodal vector to transform
  !>
  !> @author Knut Morten Okstad / Bjorn Haugen
  !>
  !> @date 15 Jul 2002

  subroutine transVecSysToGlob (triad,v)

    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: v(:)

    !! --- Logic section ---

    if (.not. associated(triad%sysDirInG)) return

    if (size(v) >= 3) v(1:3) = matmul(triad%sysDirInG,v(1:3))
    if (size(v) >= 6) v(4:6) = matmul(triad%sysDirInG,v(4:6))

  end subroutine transVecSysToGlob


  !!============================================================================
  !> @brief Transforms a triad vector from system to global directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param[in] v The nodal vector to transform
  !> @return The transformed nodal vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad / Bjorn Haugen
  !>
  !> @date 15 Jul 2002

  function transVSysToGlob (triad,v) result(u)

    type(TriadType), intent(in) :: triad
    real(dp)       , intent(in) :: v(:)
    real(dp)                    :: u(size(v))

    !! --- Logic section ---

    if (associated(triad%sysDirInG)) then
       if (size(v) >= 3) u(1:3) = matmul(triad%sysDirInG,v(1:3))
       if (size(v) >= 6) u(4:6) = matmul(triad%sysDirInG,v(4:6))
    else
       u = v
    end if

  end function transVSysToGlob


  !!============================================================================
  !> @brief Transforms a triad vector from triad to system directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param[in] v The nodal vector to transform
  !> @return The transformed nodal vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Jul 2002

  function transVTriadToSys (triad,v) result(u)

    type(TriadType), intent(in) :: triad
    real(dp)       , intent(in) :: v(:)
    real(dp)                    :: u(size(v))

    !! --- Logic section ---

    if (size(v) >= 3) u(1:3) = matmul(triad%ur(:,1:3),v(1:3))
    if (size(v) >= 6) u(4:6) = matmul(triad%ur(:,1:3),v(4:6))

    call transVecGlobToSys (triad,u)

  end function transVTriadToSys

  !> @brief Transforms a triad vector from triad to system directions.
  subroutine transVecTriadToSys (triad,v)
    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: v(:)
    v = transVTriadToSys(triad,v)
  end subroutine transVecTriadToSys


  !!============================================================================
  !> @brief Transforms a triad vector from system to triad directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param[in] v The nodal vector to transform
  !> @return The transformed nodal vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Aug 2007

  function transVSysToTriad (triad,v) result(u)

    type(TriadType), intent(in) :: triad
    real(dp)       , intent(in) :: v(:)
    real(dp)                    :: u(size(v))

    !! --- Logic section ---

    if (associated(triad%sysDirInG)) then
       if (size(v) >= 3) u(1:3) = matmul(triad%sysDirInG,v(1:3))
       if (size(v) >= 6) u(4:6) = matmul(triad%sysDirInG,v(4:6))
    else
       u = v
    end if

    if (size(v) >= 3) u(1:3) = matmul(u(1:3),triad%ur(:,1:3))
    if (size(v) >= 6) u(4:6) = matmul(u(4:6),triad%ur(:,1:3))

  end function transVSysToTriad

  !> @brief Transforms a triad vector from system to triad directions.
  subroutine transVecSysToTriad (triad,v)
    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: v(:)
    v = transVSysToTriad(triad,v)
  end subroutine transVecSysToTriad


  !!============================================================================
  !> @brief Transforms a triad vector from triad to global directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param[in] v The nodal vector to transform
  !> @return The transformed nodal vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Jul 2002

  function transVTriadToGlob (triad,v) result(u)

    type(TriadType), intent(in) :: triad
    real(dp)       , intent(in) :: v(:)
    real(dp)                    :: u(size(v))

    !! --- Logic section ---

    if (size(v) >= 3) u(1:3) = matmul(triad%ur(:,1:3),v(1:3))
    if (size(v) >= 6) u(4:6) = matmul(triad%ur(:,1:3),v(4:6))

  end function transVTriadToGlob

  !> @brief Transforms a triad vector from triad to global directions.
  subroutine transVecTriadToGlob (triad,v)
    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: v(:)
    v = transVTriadToGlob(triad,v)
  end subroutine transVecTriadToGlob


  !!============================================================================
  !> @brief Transforms a triad vector from global to triad directions.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param[in] v The nodal vector to transform
  !> @return The transformed nodal vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Dec 2000

  function transVGlobToTriad (triad,v) result(u)

    type(TriadType), intent(in) :: triad
    real(dp)       , intent(in) :: v(:)
    real(dp)                    :: u(size(v))

    !! --- Logic section ---

    if (size(v) >= 3) u(1:3) = matmul(v(1:3),triad%ur(:,1:3))
    if (size(v) >= 6) u(4:6) = matmul(v(4:6),triad%ur(:,1:3))

  end function transVGlobToTriad

  !> @brief Transforms a triad vector from global to triad directions.
  subroutine transVecGlobToTriad (triad,v)
    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: v(:)
    v = transVGlobToTriad(triad,v)
  end subroutine transVecGlobToTriad


  !!============================================================================
  !> @brief Transforms a matrix from global to system directions for a triad.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param A The matrix to transform
  !> @param[in] pos The position in the matrix @a A for the given triad
  !> @param[in] n Number of components (3 or 6) to transform
  !>
  !> @author Knut Morten Okstad
  !> @date 11 Nov 2005
  !>
  !> @author Knut Morten Okstad
  !> @date 3 Nov 2014

  subroutine transMatGlobToSys (triad,A,pos,n)

    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: A(:,:)
    integer        , intent(in)    :: pos, n

    !! Local variables
    integer  :: k, m, ierr
    real(dp) :: tmp(3)

    !! --- Logic section ---

    if (.not. associated(triad%sysDirInG)) return

    m = size(A,1)
    if (m >= 3 .and. pos >= 1 .and. pos <= m-2 .and. n >= 3) then
       do k = pos, min(pos+n-1,m)-2, 3
          call MATTRA (triad%sysDirInG,A,tmp,m,3,k,0,ierr)
       end do
    end if

  end subroutine transMatGlobToSys


  !!============================================================================
  !> @brief Transforms a matrix from triad to system directions for a triad.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param A The matrix to transform
  !> @param[in] pos The position in the matrix @a A for the given triad
  !> @param[in] n Number of components (3 or 6) to transform
  !>
  !> @author Knut Morten Okstad
  !> @date 11 Nov 2005
  !>
  !> @author Knut Morten Okstad
  !> @date 3 Nov 2014

  subroutine transMatTriadToSys (triad,A,pos,n)

    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: A(:,:)
    integer        , intent(in)    :: pos, n

    !! Local variables
    integer  :: k, m, ierr
    real(dp) :: tmp(3), Rsys(3,3)

    !! --- Logic section ---

    m = size(A,1)
    if (m >= 3 .and. pos >= 1 .and. pos <= m-2 .and. n >= 3) then

       Rsys = transpose(triad%ur(:,1:3))
       if (associated(triad%sysDirInG)) then
          Rsys = matmul(Rsys,triad%sysDirInG)
       end if

       do k = pos, min(pos+n-1,m)-2, 3
          call MATTRA (Rsys,A,tmp,m,3,k,0,ierr)
       end do

    end if

  end subroutine transMatTriadToSys


  !!============================================================================
  !> @brief Transforms a matrix from triad to global directions for a triad.
  !>
  !> @param[in] triad The triadtypemodule::triadtype object to transform for
  !> @param A The matrix to transform
  !> @param[in] pos The position in the matrix @a A for the given triad
  !> @param[in] n Number of components (3 or 6) to transform
  !>
  !> @author Knut Morten Okstad
  !> @date 11 Nov 2005
  !>
  !> @author Knut Morten Okstad
  !> @date 3 Nov 2014

  subroutine transMatTriadToGlob (triad,A,pos,n)

    type(TriadType), intent(in)    :: triad
    real(dp)       , intent(inout) :: A(:,:)
    integer        , intent(in)    :: pos, n

    !! Local variables
    integer  :: k, m, ierr
    real(dp) :: tmp(3), T(3,3)

    !! --- Logic section ---

    m = size(A,1)
    if (m >= 3 .and. pos >= 1 .and. pos <= m-2 .and. n >= 3) then
       T = transpose(triad%ur(:,1:3))
       do k = pos, min(pos+n-1,m)-2, 3
          call MATTRA (T,A,tmp,m,3,k,0,ierr)
       end do
    end if

  end subroutine transMatTriadToGlob

end module TriadTypeModule


!!==============================================================================
!> @brief Module with a namelist for reading triad data.

module TriadNamelistModule

  use KindModule  , only : dp
  use IdTypeModule, only : ldesc_p

  implicit none

  !! Define the TRIAD namelist
  integer :: id         !< Base ID of the triad
  integer :: extId(10)  !< User ID path of the triad
  integer :: nDOFs      !< Number of nodal DOFs (0, 3 or 6)
  integer :: sysDir     !< System directions flag (0-3)
  integer :: BC(6)      !< Boundary condition codes
  integer :: savePos    !< Flag indicating whether position should be saved
  integer :: saveVar(7) !< Flags indicating which variables should be saved

  character(ldesc_p) :: extDescr !< User description

  real(dp) :: ur(4,3)         !< Initial position matrix (undeformed state)
  real(dp) :: urd(6)          !< Initial velocities
  real(dp) :: urdd(6)         !< Initial accelerations
  real(dp) :: dragParams(3,3) !< Drag parameters

  namelist /TRIAD/ id, extId, extDescr, nDOFs, sysDir, ur, urd, urdd, &
       &           dragParams, BC, savePos, saveVar

contains

  !> @cond FULL_DOC
  !> @brief Reads the TRIAD namelist from the given file unit.
  subroutine read_TRIAD (infp,stat)

    integer, intent(in)  :: infp
    integer, intent(out) :: stat

    !! Default values
    id=0; extId=0; extDescr=''; nDOFs=0; sysDir=0
    ur=0.0_dp; urd=0.0_dp; urdd=0.0_dp; dragParams=0.0_dp; BC=1
    savePos=0; saveVar=0

    read(infp,nml=TRIAD,iostat=stat)

  end subroutine read_TRIAD
  !> @endcond

end module TriadNamelistModule
