!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module saveNastranModule2

  implicit none

contains

  subroutine findFENodeNumsForTriads (sam,sup,ierr)

    !!==========================================================================
    !! Find the corresponding FE node numbers for all triads on current part
    !! by comparing the global nodal coordinates. Store in triad%samNodNum.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 23 Oct 2001/1.0
    !!==========================================================================

    use SamModule              , only : SamType, dp
    use SupElTypeModule        , only : SupElType
    use IdTypeModule           , only : getId
    use manipMatrixModule      , only : matmul34
    use reportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getNodalCoor

    type(SamType)  , intent(in)    :: sam
    type(SupElType), intent(inout) :: sup
    integer        , intent(out)   :: ierr

    !! Local variables
    integer            :: i, j
    real(dp)           :: X(3), XG(3)
    character(len=128) :: errmsg

    !! --- Logic section ---

    do i = 1, sam%nnod
       if (sam%mnnn(i) > 1) then
          call ffl_getNodalCoor (X(1),X(2),X(3),i,ierr)
          if (ierr < 0) exit

          XG = matmul34(sup%supTr,X)
          do j = 1, size(sup%triads)
             if (isEqual(sup%triads(j)%p%ur(:,4),XG)) then
                sup%triads(j)%p%samNodNum = sam%minex(i)
                goto 100
             end if
          end do
          write(errMsg,"(' for node ',I6,', X = ',1P3E13.5)") sam%minex(i), XG
          call reportError (error_p,'Could not find a matching triad'//errmsg)
          ierr = ierr - 1
100       continue
       end if
    end do

    do j = 1, size(sup%triads)
       if (sup%triads(j)%p%samNodNum == 0) then
          write(errMsg,"(' for triad ',a,', X = ',1P3E13.5)") &
               &       trim(getId(sup%triads(j)%p%id)), sup%triads(j)%p%ur(:,4)
          call reportError (error_p,'Could not find a matching node'//errmsg)
          ierr = ierr - 1
       end if
    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'findFENodeNumsForTriads')

  contains

    logical function isEqual (X,Y)
      real(dp), intent(in) :: X(3), Y(3)
      real(dp), parameter  :: eps_p = 1.0e-4_dp
      real(dp)             :: tolX
      integer              :: i
      isEqual = .false.
      tolX = eps_p * sqrt(max(1.0e-24_dp,dot_product(X,X),dot_product(Y,Y)))
      do i = 1,3
         if (abs(X(i)-Y(i)) > tolX) return
      end do
      isEqual = .true.
    end function isEqual

  end subroutine findFENodeNumsForTriads


  subroutine writeSupelDeformation (istep,sup,ierr)

    !!==========================================================================
    !! Write deformatinal displacements to Nastran file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 16 Oct 2001/1.0
    !!==========================================================================

    use SupElTypeModule    , only : SupElType
    use IdTypeModule       , only : getId
    use versionModule      , only : getCurrentDate, progVer, buildDate
    use fileUtilitiesModule, only : findUnitNumber
    use reportErrorModule  , only : reportError, error_p

    integer        , intent(in)  :: istep
    type(SupElType), intent(in)  :: sup
    integer        , intent(out) :: ierr

    !! Local variables
    integer           :: c, i, j, lpu
    character(len=8)  :: chl, cht
    character(len=32) :: chname

    !! --- Logic section ---

    lpu = findUnitNumber(10)
    write(chl,'(i8)') sup%id%baseId
    write(cht,'(i8)') istep
    chname = 'adisp_L'//trim(adjustl(chl))//'_t'//trim(adjustl(cht))//'.bdf'
    open(lpu,FILE=chname,STATUS='UNKNOWN',IOSTAT=ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open output file '//chname, &
            &            addString='writeSupElDeformation')
       return
    end if

    call getCurrentDate (chname)
    write(lpu,600) trim(getId(sup%id))
    write(lpu,601) trim(chname)
    write(lpu,602) 'fedem_stress '//trim(progVer)//' '//trim(buildDate)

    do i = 1, size(sup%triads)
       j = sup%triads(i)%firstDOF
       do c = 1, sup%triads(i)%p%nDOFs
          write(lpu,610) iStep,sup%triads(i)%p%samNodNum,c,sup%finit(j)
          j = j + 1
       end do
    end do
    write(lpu,*)
    close(lpu)

600 format('$ Supernode deformations for Part ',a/'$')
601 format('$ Creation date: ',a)
602 format('$ Created by: ',a/'$')
610 format('SPCD*   ',3i16,1pe16.9)

  end subroutine writeSupElDeformation

end module saveNastranModule2
