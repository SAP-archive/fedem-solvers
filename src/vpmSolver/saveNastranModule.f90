!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module saveNastranModule

  implicit none

  private :: writeSupelDeformation


contains

  subroutine writeSolverDB3 (sys,mech,ierr)

    !!==========================================================================
    !! Administer writing of results from fedem_solver for external recovery.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Jun 2005/1.0
    !!==========================================================================

    use SystemTypeModule   , only : SystemType
    use MechanismTypeModule, only : MechanismType
    use ReportErrorModule  , only : debugFileOnly_p, reportError

    type(SystemType)   , intent(in)  :: sys
    type(MechanismType), intent(in)  :: mech
    integer            , intent(out) :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(mech%sups)
       if (associated(mech%sups(i)%nodeId)) then
          call writeSupelDeformation (sys%nStep,sys%time,mech%sups(i),ierr)
          if (ierr < 0) goto 500
       end if
    end do

    return

500 call reportError (debugFileOnly_p,'writeSolverDB3')

  end subroutine writeSolverDB3


  subroutine writeSupelDeformation (istep,time,sup,ierr)

    !!==========================================================================
    !! Write deformational displacements to Nastran bulk data file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Jun 2005/2.0
    !!==========================================================================

    use KindModule         , only : dp, i8
    use IdTypeModule       , only : getId, StrId
    use SupElTypeModule    , only : SupElType
    use VersionModule      , only : progVer, buildDate, getCurrentDate
    use FileUtilitiesModule, only : findUnitNumber
    use ReportErrorModule  , only : note_p, error_p, reportError

    integer(i8)    , intent(in)  :: iStep
    real(dp)       , intent(in)  :: time
    type(SupElType), intent(in)  :: sup
    integer        , intent(out) :: ierr

    !! Local variables
    integer           :: c, i, j, lpu
    character(len=64) :: chpart, chname

    !! --- Logic section ---

    write(chname,"(i20)") iStep
    if (sup%id%descr == '') then
       chpart = 'adisp_P'//adjustl(StrId(sup%id%baseId))
    else
       chpart = 'adisp_'//trim(adjustl(sup%id%descr))
    end if
    chname = trim(chpart)//'_t'//trim(adjustl(chname))//'.bdf'
    call reportError (note_p,'Saving supernode deformations for Part'// &
         &            trim(getId(sup%id))//' to bulk data file '//chname)

    lpu = findUnitNumber(10)
    open(lpu,FILE=chname,STATUS='UNKNOWN',IOSTAT=ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open output file '//chname)
       return
    end if

    call getCurrentDate (chname)
    if (iStep < 10000000_i8) then
       write(lpu,600) trim(getId(sup%id)),iStep,time
    else
       write(lpu,605) trim(getId(sup%id)),time
    end if
    write(lpu,601) trim(chname)
    write(lpu,602) 'fedem_solver '//trim(progVer)//' '//trim(buildDate)

    do i = 1, size(sup%nodeId)
       c = 0
       do j = 0, sup%triads(i)%p%nDOFs
          c = c + 1
          write(lpu,610) iStep,sup%nodeId(i),c, &
               &         sup%finit(sup%triads(i)%firstDOF+j)
       end do
    end do

    close(lpu)

600 format('$ Supernode deformations for Part ',a/'$ Time step:',i8,1pe12.5/'$')
605 format('$ Supernode deformations for Part ',a/'$ Time:',1pe12.5/'$')
601 format('$ Creation date: ',a)
602 format('$ Created by: ',a/'$')
610 format('SPCD*',i19,2i16,1pe16.9)

  end subroutine writeSupElDeformation

end module saveNastranModule
