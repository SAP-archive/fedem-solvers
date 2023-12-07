!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ErrorFlag

  !!============================================================================
  !! This module is a dummy implementation of the ErrorFlag module belonging to
  !! the FE linear solver from DNVS.
  !!============================================================================

  implicit none

  type Error_Flag
     integer :: dummy
  end type Error_Flag

contains

  function ZeroPivotError (ErrFlag)
    type(Error_Flag), intent(in) :: ErrFlag
    logical                      :: ZeroPivotError
    print *,' ** ZeroPivotError dummy: ',ErrFlag%dummy
    ZeroPivotError = .false.
  end function ZeroPivotError

  function NegativePivotError (ErrFlag)
    type(Error_Flag), intent(in) :: ErrFlag
    logical                      :: NegativePivotError
    print *,' ** NegativePivotError dummy: ',ErrFlag%dummy
    NegativePivotError = .false.
  end function NegativePivotError

  function Get_ErrStat (ErrFlag,Routine)
    type(Error_Flag), pointer              :: ErrFlag
    character(len=*), intent(in), optional :: Routine
    logical                                :: Get_ErrStat
    print *,' ** Get_ErrStat dummy: ',Routine,ErrFlag%dummy
    Get_ErrStat = .false.
  end function Get_ErrStat

  subroutine Set_MemoryAlloc_Error (ErrFlag,Routine,Errnum)
    type(Error_Flag), intent(inout) :: ErrFlag
    character(len=*), intent(in)    :: Routine
    integer         , intent(in)    :: Errnum
    print *,' ** Set_MemoryAlloc_Error dummy: ',Routine,ErrFlag%dummy,Errnum
  end subroutine Set_MemoryAlloc_Error

  subroutine Set_Internal_Error (ErrFlag,Routine,Errnum)
    type(Error_Flag), intent(inout) :: ErrFlag
    character(len=*), intent(in)    :: Routine
    integer         , intent(in)    :: Errnum
    print *,' ** Set_Internal_Error dummy: ',Routine,ErrFlag%dummy,Errnum
  end subroutine Set_Internal_Error

  subroutine Set_Traceback (ErrFlag,Routine)
    type(Error_Flag), intent(inout) :: ErrFlag
    character(len=*), intent(in)    :: Routine
    print *,' ** Set_Traceback dummy: ',Routine,ErrFlag%dummy
  end subroutine Set_Traceback

  subroutine Print_ErrStat (ErrFlag,LPU)
    type(Error_Flag), intent(in) :: ErrFlag
    integer         , intent(in) :: LPU
    print *,' ** Print_ErrStat dummy:',LPU,ErrFlag%dummy
  end subroutine Print_ErrStat

  subroutine Destroy_Error_Flag (ErrFlag)
    Type(Error_Flag), pointer :: ErrFlag
    if (associated(ErrFlag)) deallocate(ErrFlag)
  end subroutine Destroy_Error_Flag

end module ErrorFlag
