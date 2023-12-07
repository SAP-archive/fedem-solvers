!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ResStressModule

  implicit none

  private

#ifdef FT_HAS_VKI
  interface

     subroutine getShellStress (elmno, mnpc, Sigma, nstrp, IERR)
       use KindModule, only : dp
       integer , intent(in)    :: elmno, mnpc
       real(dp), intent(inout) :: Sigma
       integer , intent(inout) :: nstrp
       integer , intent(out)   :: IERR
     end subroutine getShellStress

     subroutine getSolidStress (elmno, mnpc, Sigma, nstrp, IERR)
       use KindModule, only : dp
       integer , intent(in)    :: elmno, mnpc
       real(dp), intent(inout) :: Sigma
       integer , intent(inout) :: nstrp
       integer , intent(out)   :: IERR
     end subroutine getSolidStress

     subroutine readResStress (int IERR)
       integer , intent(out)   :: IERR
     end subroutine readResStress

  end interface
#endif

  public :: addResStress, getSolidSurfaceStress, readResStress


contains

  subroutine addResStress (sam, Sigma, IEL, NENOD, NSTRP, IPRINT, LPU, IERR)

    !!==========================================================================
    !! Add residual stresses to the stress tensor at predefined results points
    !! (for shells) or element nodes (for solids), for the specified element.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 29 Aug 2005/1.0
    !!==========================================================================

    use SamModule              , only : SamType, dp
    use FFlLinkHandlerInterface, only : ffl_getElmId

    type(SamType), intent(in)    :: sam
    real(dp)     , intent(inout) :: Sigma(:)
    integer      , intent(in)    :: IEL, IPRINT, LPU
    integer      , intent(inout) :: NENOD, NSTRP
    integer      , intent(inout) :: IERR

    !! Local variables
    integer :: n, elmNo, ipnc1, ipnc2, lerr

    !! --- Logic section ---

    elmNo = ffl_getElmId(iel)
    if (elmno < 1) return

    ipnc1 = sam%mpmnpc(iel)
    ipnc2 = sam%mpmnpc(iel+1)
    select case (sam%melcon(iel))

    case (21:22) ! Thin shell elements
#ifdef FT_HAS_VKI
       call getShellStress (elmno,sam%mmnpc(ipnc1),Sigma(1),nstrp,lerr)
       if (nenod == 0) nenod = ipnc2-ipnc1
#else
       lerr = -99
#endif

    case (41:46) ! Solid elements
#ifdef FT_HAS_VKI
       call getSolidStress (elmno,sam%mmnpc(ipnc1),Sigma(1),nstrp,lerr)
       if (nenod == 0) nenod = nstrp
#else
       lerr = -99
#endif

    case default ! Silently ignore all the other element types
       return

    end select

    if (lerr < 0) then
       ierr = lerr
       call resStressError ('addResStress',elmNo)
       return
    else if (nenod /= ipnc2-ipnc1) then
       ierr = -1
       call resStressError ('addResStress',elmNo,ipnc2-ipnc1,nenod)
       return
    end if

    if (iprint == 0) return

    !! Print out stress results to the log-file

    select case (sam%melcon(iel))

    case (21:22) ! Thin shell elements
       write(LPU,620) elmNo,(n,Sigma(3*n-2:3*n),n=1,nstrp)

    case (41:46) ! Solid elements
       write(LPU,630) elmNo,(n,Sigma(6*n-5:6*n),n=1,nstrp)

    case default
       return

    end select

620 format(//5X,'RESIDUAL STRESSES FOR SHELL ELEMENT',I8,' :' &
         & //5X,'Point#',4X,'sigma_xx',5X,'sigma_yy',5X,'sigma_xy' &
         &  /1P,(I9,2X,3E13.5))

630 format(//5X,'RESIDUAL STRESSES FOR SOLID ELEMENT',I8,' :' &
         & //5X,'Node #',4X,'sigma_xx',5X,'sigma_yy',5X,'sigma_zz', &
         &               5X,'sigma_xy',5X,'sigma_xz',5X,'sigma_yz' &
         &  /1P,(I9,2X,6E13.5))

  end subroutine addResStress


  subroutine getSolidSurfaceStress (sam, iel, mnpc, Sigma, istat)

    !!==========================================================================
    !! Compute residual stresses on the middle of the surface defined by the
    !! nodes MNPC on the given element IEL.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 31 Aug 2005/1.0
    !!==========================================================================

    use SamModule              , only : SamType, dp
    use FFlLinkHandlerInterface, only : ffl_getElmId

    type(SamType), intent(in)  :: sam
    integer      , intent(in)  :: iel, mnpc(:)
    real(dp)     , intent(out) :: Sigma(6)
    integer      , intent(out) :: istat

    !! Local variables
    integer, parameter :: maxnod_p = 20
    integer            :: i, j, elmNo, ipnc1, ipnc2, nenod
    real(dp)           :: SigNo(6,maxnod_p)

    !! --- Logic section ---

    istat = 0
    Sigma = 0.0_dp

    if (iel < 1 .or. iel > size(sam%melcon)) return ! Element index out of range
    if (sam%melcon(iel)/10 /= 4)             return ! Only for solid elements

    elmNo = abs(ffl_getElmId(iel))
    if (elmno == 0)                          return ! Non-existing element

    ipnc1 = sam%mpmnpc(iel)
    ipnc2 = sam%mpmnpc(iel+1)
    SigNo = 0.0_dp
    nenod = 0
#ifdef FT_HAS_VKI
    call getSolidStress (elmNo,sam%mmnpc(ipnc1),SigNo(1,1),nenod,istat)
#else
    istat = -99
#endif
    if (istat < 0) then
       call resStressError ('getSolidSurfaceStress',elmNo)
       return
    else if (nenod /= ipnc2-ipnc1 .or. nenod > maxnod_p) then
       istat = -1
       call resStressError ('getSolidSurfaceStress',elmNo,ipnc2-ipnc1,nenod)
       return
    end if

    !! Calculate mid-point surface stress by averaging associated nodal values

    do i = 1, size(mnpc)
       do j = 1, nenod
          if (sam%mmnpc(ipnc1+j-1) == mnpc(i)) then
             istat = istat + 1
             Sigma = Sigma + SigNo(:,j)
             exit
          end if
       end do
    end do

    if (istat < size(mnpc)) then
       istat = -size(mnpc)
       call resStressError ('getSolidSurfaceStress',elmNo,istat)
    else
       Sigma = Sigma / real(istat,dp)
    end if

  end subroutine getSolidSurfaceStress


  subroutine resStressError (prgnam,elmNo,n1,n2)

    !!==========================================================================
    !! Print an error message on residual stress calculation failure.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Feb 2006/1.0
    !!==========================================================================

    use ReportErrorModule, only : reportError, error_p, debugFileOnly_p

    character(len=*) , intent(in) :: prgnam
    integer          , intent(in) :: elmNo
    integer, optional, intent(in) :: n1, n2

    !! Local variables
    character(len=128) :: msg1, msg2

    !! --- Logic section ---

    if (present(n1) .and. present(n2)) then
       write(msg1,691) elmNo
       write(msg2,692) n1, n2
       call reportError (error_p,msg1,msg2,'The specified residual '// &
            'stress file is not compatible with the current Fedem model')
    else if (present(n1)) then
       write(msg1,693) elmNo
       call reportError (error_p,msg1)
    else
       write(msg1,694) elmNo
       call reportError (error_p,msg1)
    end if
    call reportError (debugFileOnly_p,prgnam)

691 format('Node topology mismatch detected for element',i8)
692 format('There are',i4,' element nodes in Fedem, but', &
         & i4,' element nodes on the external file')
693 format('Strain coat topology error for element',i8)
694 format('Stress calculation failed for element',i8)

  end subroutine resStressError

#ifndef FT_HAS_VKI
  subroutine readResStress (IERR)
    integer , intent(out)   :: IERR
    IERR = -99
    print *,'*** readResStress: Compiled without VKI support.'
  end subroutine readResStress
#endif

end module ResStressModule
