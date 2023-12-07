!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file windTurbineTypeModule.f90
!> @brief Wind turbine data container.

!!==============================================================================
!> @brief Module with wind turbine data containers.
!>
!> @details This module contains data containers for all components that
!> together define a wind turbine object. It also contains some subroutines
!> for accessing or manipulating the wind turbine data.

module windTurbineTypeModule

  use SupElTypeModule, only : SupElType, IdType, dp
  use TriadTypeModule, only : TriadType
#if FT_HAS_AERODYN > 12
  use AeroDyn, only : AeroConfig, AllAeroMarkers, AllAeroLoads, AeroLoadsOptions
#endif

  implicit none

  !> @brief Data type describing a turbine blade element.
  type BladeElem
     real(dp) :: R           !< Distance from blade root to aerodynamic center
     real(dp) :: DR          !< Length of each blade element (along pitch axis)
     real(dp) :: ACoffset(2) !< Offset from pitch axis to aerodynamic center
  end type BladeElem

  !> @brief Data type describing the root node of a turbine blade.
  type BladeRoot
     type(TriadType), pointer :: cone  !< Unpitched blade coordinate system
     type(TriadType), pointer :: triad !< Pitched blade coordinate system
     type(TriadType), pointer :: tip   !< Blade tip coordinate system
     real(dp)       , pointer :: pitch !< Current pitch angle
     real(dp) :: preCone !< Precone angle of the blade
  end type BladeRoot

  !> @brief Data type describing a turbine blade node.
  type BladeNode
     type(TriadType), pointer :: triad !< Aerodynamic center of this element
     real(dp) :: aeroForce(3)     !< Aerodynamic force on this element
     real(dp) :: aeroForcePrev(6) !< Global aerodynamic force in previous step
  end type BladeNode

  !> @brief Data type describing a wind turbine configuration.
  type TurbineConfig
     type(IdType) :: id  !< General identification data
     real(dp) :: PtfmRef !< Platform reference height
     real(dp) :: RefHH   !< Reference height of the Hub
     real(dp) :: HubRad  !< Preconed Hub radius
     real(dp) :: TipRad  !< Preconed blade-tip radius
     real(dp)       , pointer :: windPt(:,:) !< Userdefined wind sampling points
     real(dp)       , pointer :: genDOF(:)   !< Generator DOF (ur, urd, urdd)
     type(TriadType), pointer :: tower       !< Tower base triad
     type(TriadType), pointer :: nacelle     !< Nacelle triad
     type(TriadType), pointer :: shaft       !< Shaft triad
     type(TriadType), pointer :: azimuth     !< Azimuth triad
     type(TriadType), pointer :: hub         !< Hub triad
     type(SupElType), pointer :: hubEl       !< Hub superelement
     type(BladeRoot), pointer :: blade(:)    !< The root node of each blade
     type(BladeNode), pointer :: node(:,:)   !< Aerodynamic nodes of the blades
     type(BladeElem), pointer :: bladeElm(:) !< Blade geometry data for AeroDyn
#if FT_HAS_AERODYN > 12
     type(AeroConfig)       :: ADinterface !< AeroDyn bodies
     type(AllAeroMarkers)   :: ADmarkers   !< AeroDyn nodes
     type(AllAeroLoads)     :: ADloads     !< AeroDyn loads
     type(AeroLoadsOptions) :: ADoptions   !< AeroDyn load options
#endif
  end type TurbineConfig

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteTurbineConfig
  end interface


contains

  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] turb Wind turbine configuration object
  !> @param[in] io File unit number to write to
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2009

  subroutine WriteTurbineConfig (turb,io)

    use IdTypeModule, only : writeId

    type(TurbineConfig), intent(in) :: turb
    integer            , intent(in) :: io

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'TurbineConfiguration','{'
    call writeId (turb%id,io)

    write(io,*) 'PtfmRef     =', turb%PtfmRef
    write(io,*) 'RefHH       =', turb%RefHH
    write(io,*) 'HubRad      =', turb%HubRad
    write(io,*) 'TipRad      =', turb%TipRad
    if (associated(turb%tower)) then
       write(io,*) 'tower(id)   =', turb%tower%id%baseId
    end if
    if (associated(turb%nacelle)) then
       write(io,*) 'nacelle(id) =', turb%nacelle%id%baseId
    end if
    if (associated(turb%shaft)) then
       write(io,*) 'shaft(id)   =', turb%shaft%id%baseId
    end if
    if (associated(turb%azimuth)) then
       write(io,*) 'azimuth(id) =', turb%azimuth%id%baseId
    end if
    if (associated(turb%hub)) then
       write(io,*) 'hub(id)     =', turb%hub%id%baseId
    end if
    if (associated(turb%hubEl)) then
       write(io,*) 'hubEl(id)   =', turb%hubEl%id%baseId
    end if
    do i = 1, size(turb%blade)
       if ( associated(turb%blade(i)%triad) .and. &
            associated(turb%blade(i)%tip) ) then
          write(io,600) i, turb%blade(i)%triad%id%baseId, &
               &           turb%blade(i)%tip%id%baseId
       end if
       write(io,610) 'preCone',i, turb%blade(i)%preCone
       write(io,610) 'pitch',i, turb%blade(i)%preCone
    end do
    if (associated(turb%genDOF)) then
       write(io,*) 'genDOF      =', turb%genDOF(1:2)
    end if
    if (associated(turb%windPt)) then
       write(io,620) (i,turb%windPt(:,i),i=1,size(turb%windPt,2))
    end if

    write(io,'(A/)') '}'

600 format(' blade',I1,'(id)  =',I12)
610 format(1X,A,I1,T14,'=',1PE12.5)
620 format(' windPt',I2,4X,'=',1P3E13.5)

  end subroutine WriteTurbineConfig


  !!============================================================================
  !> @brief Initializes the TurbineConfig object.
  !>
  !> @param[out] turb Wind turbine configuration object
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2009

  subroutine NullifyTurbine (turb)

    use IdTypeModule, only : nullifyId

    type(TurbineConfig), intent(out) :: turb

    !! --- Logic section ---

    call nullifyId (turb%id)

    turb%PtfmRef = 0.0_dp
    turb%RefHH   = 0.0_dp
    turb%HubRad  = 0.0_dp
    turb%TipRad  = 0.0_dp

    nullify(turb%windPt)
    nullify(turb%genDOF)
    nullify(turb%tower)
    nullify(turb%nacelle)
    nullify(turb%shaft)
    nullify(turb%azimuth)
    nullify(turb%hub)
    nullify(turb%hubEl)
    nullify(turb%blade)
    nullify(turb%node)
    nullify(turb%bladeElm)

  end subroutine NullifyTurbine


  !!============================================================================
  !> @brief Deallocates the TurbineConfig object.
  !>
  !> @param turb Wind turbine configuration object
  !>
  !> @author Knut Morten Okstad
  !>
  !> @callergraph
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateTurbine (turb)

    use IdTypeModule, only : deallocateId

    type(TurbineConfig), intent(inout) :: turb

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call deallocateId (turb%id)

    do i = 1, size(turb%blade)
       if (associated(turb%blade(i)%tip)) then
          if (associated(turb%blade(i)%tip%ur_def)) then
             deallocate(turb%blade(i)%tip%ur_def)
          end if
       end if
    end do

    if (associated(turb%windPt))   deallocate(turb%windPt)
    if (associated(turb%blade))    deallocate(turb%blade)
    if (associated(turb%node))     deallocate(turb%node)
    if (associated(turb%bladeElm)) deallocate(turb%bladeElm)
#if FT_HAS_AERODYN > 12
    if (allocated(turb%ADoptions%SetMulTabLoc)) then
       deallocate(turb%ADoptions%SetMulTabLoc)
    end if
    if (allocated(turb%ADinterface%Blade)) then
       deallocate(turb%ADinterface%Blade)
    end if
#endif

    call nullifyTurbine (turb)

  end subroutine DeallocateTurbine

end module windTurbineTypeModule
