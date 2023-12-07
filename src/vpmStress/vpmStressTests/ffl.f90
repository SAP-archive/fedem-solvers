!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @brief Dummy FFlLib interface for the purpose of stress unit testing.

module FFlLinkHandlerInterface

  implicit none

  integer, private, parameter :: dp = kind(1.0D0)

contains

   function ffl_getelmid (iel)
     integer, intent(in) :: iel
     integer             :: ffl_getelmid
     ffl_getelmid = iel
   end function ffl_getelmid

   subroutine ffl_getcoor (x, y, z, iel, ierr)
     real(dp), intent(out) :: x(*), y(*), z(*)
     integer , intent(in)  :: iel
     integer , intent(out) :: ierr
     integer :: i, n
     select case (iel)
     case (1); n = 6
        x(1:n) = 0.0_dp
        y(1:n) = 0.0_dp
        z(1:n) = 0.0_dp
        x(2)   = 1.0_dp
        y(3)   = 1.0_dp
        x(4:5) = 0.5_dp
        y(5:6) = 0.5_dp
     case (2); n = 8
        x(1:n) = 0.0_dp
        y(1:n) = 0.0_dp
        z(1:n) = 0.0_dp
        x(2)   = 0.5_dp
        x(3:5) = 1.0_dp
        x(6)   = 0.5_dp
        y(4)   = 0.5_dp
        y(5:7) = 1.0_dp
        y(8)   = 0.5_dp
     case (11); n = 3
        x(1:n) = 0.0_dp
        y(1:n) = 0.0_dp
        z(1:n) = 0.0_dp
        x(2)   = 1.0_dp
        y(3)   = 1.0_dp
     case (12); n = 4
        x(1:n) = 0.0_dp
        y(1:n) = 0.0_dp
        z(1:n) = 0.0_dp
        x(2:3) = 1.0_dp
        y(3:4) = 1.0_dp
     case default
        n = 0
     end select
     ierr = 0
     print 100,iel,(x(i),y(i),z(i),i=1,n)
100  format(/'Element',I3,' :',3F6.2/(12X,3F6.2))
   end subroutine ffl_getcoor

   subroutine ffl_getmat (E, nu, rho, iel, ierr)
     real(dp), intent(out) :: E, nu, rho
     integer , intent(in)  :: iel
     integer , intent(out) :: ierr
     E = 2.1e11_dp
     nu = 0.3_dp
     rho = 7850.0_dp
     ierr = 0
     print 100,iel,E,nu,rho
100  format('Element',I3,' :',1P3E10.2)
   end subroutine ffl_getmat

   subroutine ffl_getthick (Th, iel, ierr)
     real(dp), intent(out) :: Th(*)
     integer , intent(in)  :: iel
     integer , intent(out) :: ierr
     integer :: n
     select case (iel)
     case (1); n = 6
     case (2); n = 8
     case(11); n = 3
     case(12); n = 4
     case default
        return
     end select
     Th(1:n) = 0.01_dp
     ierr = 0
     print 100,iel,Th(1)
100  format('Element',I3,' :',F6.2)
   end subroutine ffl_getthick

   subroutine ffl_getpinflags (pA, pB, iel, ierr)
     integer , intent(in)  :: iel
     integer , intent(out) :: pA, pB, ierr
     pA = 0
     pB = 0
     ierr = 0
     print 100,iel,pA,pB
100  format('Element',I3,' :',2I2)
   end subroutine ffl_getpinflags

   subroutine ffl_getbeamsection (sec, iel, ierr)
     real(dp), intent(out) :: sec(1)
     integer , intent(in)  :: iel
     integer , intent(out) :: ierr
     sec = 0.0_dp
     ierr = -9999
     print 100,iel,sec
100  format('Element',I3,' :',F4.1)
   end subroutine ffl_getbeamsection

end module FFlLinkHandlerInterface
