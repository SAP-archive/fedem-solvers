!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module TireTypeModule

  use KindModule     , only : dp, lfnam_p
  use TriadTypeModule, only : TriadType, IdType
  use RoadTypeModule , only : RoadType

  implicit none

  character(len=10), parameter :: tireTypes_p(6) = (/ 'MF-TYRE   ', &
       &                                              'MF-MC-TYRE', &
       &                                              'FTIRE     ', &
       &                                              'SWIFT     ', &
       &                                              'JD-TIRE   ', &
       &                                              'FEDEM     '  /)

  character(len=3), parameter :: interfaces_p(2) = (/ 'STI', &
       &                                              'CTI'  /)

  integer, parameter :: MF_p      = 1, &
       &                MF_MC_p   = 2, &
       &                FTIRE_p   = 3, &
       &                SWIFT_p   = 4, &
       &                JD_TIRE_p = 5, &
       &                FEDEM_p   = 6

  integer, parameter :: STI_p = 1, &
       &                CTI_p = 2


  type STIapiType

     character(len=256) :: tyrmod ! Tire model description returned by DTYRE
     integer            :: iswtch ! Tire characteristics
     !                            ! =  4 : Steady state forces/torques in XYZ
     !                            ! = 14 : Dynamic forces/torques in XYZ
     !                            ! <  0 : Swapped tire characteristics

     !! Special STI arrays
     integer           :: ndeqvar   ! number of tire state variables
     real(dp), pointer :: deqini(:) ! initial value for the state variables
     real(dp), pointer :: deqvar(:), deqvar_prev(:) ! state variable values
     real(dp), pointer :: deqder(:), deqder_prev(:) ! state variable gradients

     integer           :: nvars     ! Dimension of varinf
     real(dp), pointer :: varinf(:) ! Tire model variables

     integer           :: ntyparr
     real(dp), pointer :: typarr(:)

     integer           :: nwork
     real(dp), pointer :: wrkarr(:)

     integer           :: niwork
     integer,  pointer :: iwrkarr(:)

  end type STIapiType


  type TireType

     type(IdType) :: id    ! General identification data

     integer  :: idIn      ! Index into the global tire array
     integer  :: samElNum  ! Element number for SAM reference (MPMNPC)
     integer  :: type, api ! Tire model type and solver interface

     type(TriadType) , pointer :: WCtriad  ! Wheel Carrier (independent triad)
     type(TriadType) , pointer :: RimTriad ! Spindel triad (dependent triad)
     real(dp)        , pointer :: pJvar(:) ! Joint variables for tire rotation
     !                                     ! relative to wheel carrier

     type(STIapiType), pointer :: sti      ! Data used only by the STI-interface
     type(RoadType)  , pointer :: road     ! Road data

     real(dp) :: radius, defl      ! Tire radius and radial deflection
     real(dp) :: roadContactInG(3) ! Position of tire contact with road
     logical  :: isInContact       ! On the road or not
     logical  :: isSliding         ! Exceeding friction angle or not

     !! Wheel Carrier definitions
     logical  :: hasEccVec  ! Somewhat redundant, but used to speed up
     !                      ! calculations when WC is in origin of triad system
     real(dp) :: Zoffset    ! offset along rotation axis
     real(dp) :: WCinT(3,4) ! Wheel Carrier relative to triad
     !                      ! set at initialization and constant during dynamics
     !                      ! (BH: possible z-offset relative to the rim triad)
     real(dp) :: WCinG(3,4) ! Wheel Carrier in global system, calculated from
     !                      ! triad position and WC position relative to triad
     logical  :: WCYalongZ  ! Y-axis of WC system along Z-axis of revolute joint

     real(dp) :: rStiff(2)  ! Radial stiffness = rStiff(1) + rStiff(2)*defl
     real(dp) :: rDamp      ! Radial damping (assumed constant)
     real(dp) :: yStiff     ! Transverse stiffness for JD tire
     real(dp) :: yDamp(2)   ! Transverse and rotational damping for JD tire
     real(dp) :: mass       ! Portion of tire mass that is fixed to the rim
     real(dp) :: Iyy, Izz   ! Portion of tire inertia that is fixed to the rim

     real(dp) :: forceInWC(6)    ! Force from tire to rim at wheel center
     real(dp) :: forceInG(6)     ! Force from tire to rim in global coordinates
     real(dp) :: forceInGPrev(6) ! Same as above, previous time step

     real(dp) :: ePot0      ! Initial potential energy
     real(dp) :: ePot       ! Potential energy associated with rim-fixed mass
     real(dp) :: eKin       ! Kinetic energy associated with rim-fixed mass
     real(dp) :: eInp       ! Energy contribution from this tire

     logical  :: saveVar(9) ! Flags indicating which variables should be saved

     integer            :: nCharTireDataFileName
     character(lfnam_p) :: tireDataFileName

  end type TireType


  type ScaleToType
     real(dp) :: kg, m, N, s ! Scaling factors from modeling units to SI
  end type ScaleToType


  interface WriteObject
     module procedure WriteTireType
  end interface


contains

  subroutine WriteTireType (tire,io,complexity)

    !!==========================================================================
    !! Standard routine for writing an object to io.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Aug 2002/1.0
    !!==========================================================================

    use IdTypeModule, only : writeId

    type(TireType)  , intent(in) :: tire
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'Tire','{'
    call writeId (tire%id,io)

    write(io,*) 'samElNum    = ', tire%samElNum
    write(io,*) 'type        = ', tire%type
    write(io,*) 'API         = ', tire%api
    write(io,*) 'radius      = ', tire%radius
    write(io,*) 'roadContact = ', tire%roadContactInG
    write(io,*) 'Zoffset     = ', tire%Zoffset
    write(io,*) 'rStiff      = ', tire%rStiff
    write(io,*) 'rDamp       = ', tire%rDamp
    if (tire%type == JD_TIRE_p) then
       write(io,*) 'yStiff      = ', tire%yStiff
       write(io,*) 'yDamp       = ', tire%yDamp
    end if
    write(io,*) 'mass        = ', tire%mass
    write(io,*) 'Iyy, Izz    = ', tire%Iyy, tire%Izz

    if (associated(tire%WCTriad)) then
       write(io,*) 'WCTriad(id) = ', tire%WCTriad%id%baseId
    end if
    if (associated(tire%RimTriad)) then
       write(io,*) 'RimTriad(id)= ', tire%RimTriad%id%baseId
    end if
    if (associated(tire%road)) then
       write(io,*) 'road(id)    = ', tire%road%id%baseId
    end if

    write(io,601) ' WCinT       =', (tire%WCinT(i,:),i=1,3)
    write(io,601) ' WCinG       =', (tire%WCinG(i,:),i=1,3)

    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,602) ' forceInWC   =', tire%forceInWC
          write(io,602) ' forceInG    =', tire%forceInG
          write(io,*) 'deflection  = ', tire%defl
          write(io,*) 'Epot0       = ', tire%ePot0
          write(io,*) 'Epot        = ', tire%ePot
          write(io,*) 'Ekin        = ', tire%eKin
          write(io,*) 'Einp        = ', tire%eInp
       end if
    end if

    write(io,'(A)') '}'

601 format(A14,4ES15.6E3/(14X,4ES15.6E3))
602 format(A14,9ES15.6E3)

  end subroutine WriteTireType


  subroutine NullifyTire (tire)

    !!==========================================================================
    !! Initialize the TireType object.
    !!
    !! Programmer : Bjorn Haugen                         date/rev : Aug 2002/1.0
    !!==========================================================================

    use IdTypeModule, only : nullifyId

    type(TireType), intent(out) :: tire

    !! --- Logic section ---

    call nullifyId (tire%id)

    tire%idIn = 0
    tire%samElNum = 0
    tire%type = 0
    tire%api = 0

    nullify(tire%WCTriad)
    nullify(tire%RimTriad)
    nullify(tire%pJvar)
    nullify(tire%sti)
    nullify(tire%road)

    tire%radius = 0.0_dp
    tire%defl = 0.0_dp
    tire%roadContactInG = 0.0_dp
    tire%isInContact = .false.
    tire%isSliding = .false.

    tire%hasEccVec = .false.
    tire%Zoffset = 0.0_dp
    tire%WCinT = 0.0_dp
    tire%WCinG = 0.0_dp
    tire%WCYalongZ = .true.
    tire%rStiff = 0.0_dp
    tire%rDamp = 0.0_dp
    tire%yStiff = 0.0_dp
    tire%yDamp = 0.0_dp
    tire%mass = 0.0_dp
    tire%Iyy = 0.0_dp
    tire%Izz = 0.0_dp
    tire%forceInWC = 0.0_dp
    tire%forceInG = 0.0_dp
    tire%forceInGPrev = 0.0_dp
    tire%ePot0 = 0.0_dp
    tire%ePot = 0.0_dp
    tire%eKin = 0.0_dp
    tire%eInp = 0.0_dp

    tire%saveVar = .false.
    tire%nCharTireDataFileName = 0
    tire%tireDataFileName = ''

  end subroutine NullifyTire


  subroutine NullifySTIapi (tireChar,sti)

    !!==========================================================================
    !! Initialize the STIapiType object.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Jan 2003/1.0
    !!==========================================================================

    integer         , intent(in)  :: tireChar
    type(STIapiType), intent(out) :: sti

    !! --- Logic section ---

    sti%tyrmod = ''
    sti%iswtch = tireChar
    sti%ndeqvar = 0
    sti%nvars = 0
    sti%ntyparr = 0
    sti%nwork = 0
    sti%niwork = 0
    nullify(sti%deqini)
    nullify(sti%deqvar)
    nullify(sti%deqder)
    nullify(sti%deqvar_prev)
    nullify(sti%deqder_prev)
    nullify(sti%varinf)
    nullify(sti%typarr)
    nullify(sti%wrkarr)
    nullify(sti%iwrkarr)

  end subroutine NullifySTIapi


  subroutine AllocateSTIapi (sti,ierr)

    !!==========================================================================
    !! Allocate arrays of the STIapiType object.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Jan 2003/1.0
    !!==========================================================================

    type(STIapiType), intent(inout) :: sti
    integer         , intent(out)   :: ierr

    !! --- Logic section ---

    allocate(sti%deqini(sti%ndeqvar), &
        &    sti%deqvar(sti%ndeqvar), &
        &    sti%deqder(sti%ndeqvar), &
        &    sti%deqvar_prev(sti%ndeqvar), &
        &    sti%deqder_prev(sti%ndeqvar), &
        &    sti%varinf(sti%nvars), &
        &    sti%typarr(sti%ntyparr), &
        &    sti%wrkarr(sti%nwork), &
        &    sti%iwrkarr(sti%niwork), &
        &    STAT=ierr)
    if (ierr /= 0) return

    sti%deqini = 0.0_dp
    sti%deqvar = 0.0_dp
    sti%deqder = 0.0_dp
    sti%deqvar_prev = 0.0_dp
    sti%deqder_prev = 0.0_dp
    sti%varinf = 0.0_dp
    sti%typarr = 0.0_dp

  end subroutine AllocateSTIapi


  subroutine DeallocateTire (tire)

    !!==========================================================================
    !! Deallocate the TireType object.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Jan 2017/1.0
    !!==========================================================================

    use IdTypeModule, only : deallocateId

    type(TireType), intent(inout) :: tire

    !! --- Logic section ---

    call deallocateId (tire%id)

    if (associated(tire%sti)) then
       if (associated(tire%sti%varinf)) then
          deallocate(tire%sti%deqini)
          deallocate(tire%sti%deqvar)
          deallocate(tire%sti%deqder)
          deallocate(tire%sti%deqvar_prev)
          deallocate(tire%sti%deqder_prev)
          deallocate(tire%sti%varinf)
          deallocate(tire%sti%typarr)
          deallocate(tire%sti%wrkarr)
          deallocate(tire%sti%iwrkarr)
       end if
       deallocate(tire%sti)
    end if

    call nullifyTire (tire)

  end subroutine DeallocateTire


  subroutine deallocateTires (tires)
    type(TireType), pointer :: tires(:)
    integer :: i
    do i = 1, size(tires)
       call DeallocateTire (tires(i))
    end do
    deallocate(tires)
    nullify(tires)
  end subroutine deallocateTires


  function hasTireElements (tires,tireModel)

    !!==========================================================================
    !! Checks the presence of tires in the mechanism.
    !!
    !! Programmer: Jens Lien                              date/rev: Oct 2002/1.0
    !!             Knut Morten Okstad                               May 2003/2.0
    !!==========================================================================

    type(TireType)  , intent(in) :: tires(:)
    integer,optional, intent(in) :: tireModel
    logical                      :: hasTireElements

    !! Local variables
    integer :: i

    !! --- Logic section ---

    hasTireElements = size(tires) > 0

    if (hasTireElements .and. present(tireModel)) then

       do i = 1, size(tires)
          if (tires(i)%type == tireModel) return
       end do

       hasTireElements = .false. ! No tires of the specified type in mechanism

    end if

  end function hasTireElements

end module TireTypeModule
