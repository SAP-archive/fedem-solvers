!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module JD_tireRoutinesModule

  use KindModule        , only : dp
  use TriadTypeModule   , only : TriadType
  use FunctionTypeModule, only : EngineType

  implicit none

  type MarkerType
     type(TriadType) , pointer :: pTriad
     real(dp)        , pointer :: pWCinT(:,:)
     real(dp)        , pointer :: pWCinG(:,:)
     type(EngineType), pointer :: xTrEng, yTrEng, zTrEng
  end type MarkerType

  type(MarkerType), allocatable :: markers(:)


contains

  subroutine jdt_allocate_arrays (nJDT,stat)

    !!=======================================================
    !! Allocate internal arrays for JD tires
    !!=======================================================

    use jdt_fdm, only : jdt_allocate_fdm
    use jdt_inp, only : jdt_allocate_inp
    use jdt_var, only : jdt_allocate_var

    integer, intent(in)  :: nJDT
    integer, intent(out) :: stat

    !! --- Logic section ---

    allocate(markers(nJDT), STAT=stat)
    if (stat == 0) call jdt_allocate_fdm(nJDT,stat)
    if (stat == 0) call jdt_allocate_inp(nJDT,stat)
    if (stat == 0) call jdt_allocate_var(nJDT,stat)

  end subroutine jdt_allocate_arrays


  subroutine jdt_InitializeMarkerData (idIn,tire,road,stat)

    !!=======================================================
    !! Initialize tire and road pointers in the markers array
    !!=======================================================

    use TireTypeModule, only : TireType, RoadType

    integer       , intent(in)         :: idIn
    type(TireType), intent(in), target :: tire
    type(RoadType), intent(in), target :: road
    integer       , intent(out)        :: stat

    !! --- Logic section ---

    if (idIn < 1 .or. idIn > size(markers)) then
       stat = min(-1,-idIn)
       return
    end if

    markers(idIn)%pTriad => tire%WCtriad
    markers(idIn)%pWCinT => tire%WCinT
    markers(idIn)%pWCinG => tire%WCinG

    if (associated(road%xTrEng)) then
       markers(idIn)%xTrEng => road%xTrEng
    else
       nullify(markers(idIn)%xTrEng)
    end if

    if (associated(road%yTrEng)) then
       markers(idIn)%yTrEng => road%yTrEng
    else
       nullify(markers(idIn)%yTrEng)
    end if

    if (associated(road%zTrEng)) then
       markers(idIn)%zTrEng => road%zTrEng
    else
       nullify(markers(idIn)%zTrEng)
    end if

  end subroutine jdt_InitializeMarkerData


  subroutine jdt_stiffness (lDyn,tire,s,scaleK,eM)

    !!=======================================================
    !! Calculate element stiffness matrix for a JD tire
    !!=======================================================

    use TireTypeModule        , only : TireType, dp
    use JDT_VAR               , only : stiffM
    use ManipMatrixModule     , only : dyadic_product
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgTire
    use ManipMatrixModule     , only : writeObject
#endif
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    logical       , intent(in)  :: lDyn
    type(TireType), intent(in)  :: tire
    real(dp)      , intent(in)  :: s(:), scaleK
    real(dp)      , intent(out) :: eM(:,:)

    !! Local variables
    integer :: JDTstiff

    !! --- Logic section ---

    eM = 0.0_dp
    call ffa_cmdlinearg_getint ('JDTstiffness',JDTstiff)
    if (JDTstiff <= 0 .and. lDyn) return
    if (JDTstiff == 2 .and. lDyn) then

#ifdef FT_DEBUG
       write(dbgTire,"('==== New Stiffness, lDynamics = ',L3,' ====')") lDyn
#endif
       eM = stiffM(:,:,tire%idIn)*scaleK
       eM(2,2) = eM(2,2) + tire%yStiff*scaleK

    else
#ifdef FT_DEBUG
       write(dbgTire,"('==== Old Stiffness, lDynamics = ',L3,' ====')") lDyn
#endif
       eM(1:3,1:3) = dyadic_product(tire%WCinG(:,1))*s(1) &
            &      + dyadic_product(tire%WCinG(:,2))*s(2) &
            &      + dyadic_product(tire%WCinG(:,3))*s(3)
    end if

#ifdef FT_DEBUG
    write(dbgTire,"('WCinG X :',1p,3e12.3)") tire%WCinG(:,1)
    write(dbgTire,"('WCinG Y :',1p,3e12.3)") tire%WCinG(:,2)
    write(dbgTire,"('WCinG Z :',1p,3e12.3)") tire%WCinG(:,3)
    write(dbgTire,"('s(1..3) =',1p,3e12.3)") s(1:3)
    call writeObject(eM,dbgTire,'JDTire stiffness matrix')
#endif

  end subroutine jdt_stiffness

end module JD_tireRoutinesModule
