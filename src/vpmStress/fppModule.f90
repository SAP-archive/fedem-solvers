!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module fppModule

  use PointerKindModule, only : ptr

  implicit none

  integer(ptr), save, private :: fppObjectHandle = 0
  integer(ptr), save, private :: msdHandle = 0

  private :: reportFPPerror


contains

  subroutine reportFPPerror (subr,iret,ierr)

    !!==========================================================================
    !! Output an error message from the fpp processor.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 Apr 2008/1.0
    !!==========================================================================

    use ReportErrorModule           , only : internalError, reportError, error_p
    use FPPprocessorHandlerInterface, only : fpp_error

    character(len=*), intent(in)    :: subr
    integer         , intent(in)    :: iret
    integer         , intent(inout) :: ierr

    !! Local variables
    character(len=128) :: errMsg

    !! --- Logic section ---

    if (iret < 0) then
       ierr = ierr + internalError('reported by '//subr)
    else if (iret == 0) then
       ierr = ierr - 1
       call fpp_error (errMsg)
       call reportError (error_p,errMsg,addString='reported by '//subr)
    end if

  end subroutine reportFPPerror


  subroutine fppOpen (fppFile,supId,stopTime,nMat,ierr)

    !!==========================================================================
    !! Open the fpp-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 Apr 2008/1.0
    !!==========================================================================

    use KindModule                  , only : dp
    use IdTypeModule                , only : IdType
    use ReportErrorModule           , only : reportError, error_p
    use FPPprocessorHandlerInterface, only : nf_fpplib_create, nf_fpplib_open
    use FPPprocessorHandlerInterface, only : nf_fpplib_init
    use FPPprocessorHandlerInterface, only : nf_fpplib_get_data_msd

    character(len=*), intent(in)  :: fppFile
    type(IdType)    , intent(in)  :: supId
    real(dp)        , intent(in)  :: stopTime
    integer         , intent(in)  :: nMat
    integer         , intent(out) :: ierr

    !! Local variables
    integer :: iret

    !! --- Logic section ---

    fppObjectHandle = nf_fpplib_create()
    if (fppObjectHandle == 0) then
       ierr = -1
       call reportError (error_p,'Failure creating fpp object', &
            &            addString='reported by fppOpen')
       return
    end if

    call nf_fpplib_open (fppObjectHandle,trim(fppFile),iret)
    if (iret /= 1) ierr = ierr - 1

    call nf_fpplib_init (fppObjectHandle,supId%baseId,supId%userId, &
         &               1,stopTime,nMat,iret)
    if (iret /= 1) ierr = ierr - 1

    call nf_fpplib_get_data_msd (fppObjectHandle,msdHandle,iret)
    if (iret /= 1) ierr = ierr - 1

    if (ierr /= 0) then
       call reportError (error_p,'Failure initializing fpp-file '//fppFile, &
            &            addString='reported by fppOpen')
    end if

  end subroutine fppOpen


  subroutine fppClose (ierr)

    !!==========================================================================
    !! Close the fpp-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 Apr 2008/1.0
    !!==========================================================================

    use FPPprocessorHandlerInterface, only : nf_fpplib_close, nf_fpplib_release

    integer, intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0
    call nf_fpplib_close (fppObjectHandle,0,ierr)
    call nf_fpplib_release (fppObjectHandle,ierr)

  end subroutine fppClose


  subroutine fppInit (strCoats,ierr)

    !!==========================================================================
    !! Initialize the fpp-file part of the strain coat elements.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 Apr 2008/1.0
    !!==========================================================================

    use IdTypeModule                , only : GetId, StrId
    use StrainCoatModule            , only : StrainCoatType
    use ReportErrorModule           , only : reportError, error_p
    use FPPprocessorHandlerInterface, only : fpp_delete_all, fpp_create
    use FPPprocessorHandlerInterface, only : fpp_set_fpp_handle
    use FPPprocessorHandlerInterface, only : fpp_set_data_msd
    use FPPprocessorHandlerInterface, only : fpp_set_element
    use FPPprocessorHandlerInterface, only : fpp_init
    use FFaProfilerInterface        , only : ffa_startTimer, ffa_stopTimer

    type(StrainCoatType), intent(inout) :: strCoats(:)
    integer             , intent(out)   :: ierr

    !! Local variables
    integer :: i, j, iret

    !! --- Logic section ---

    ierr = 0
    call fpp_delete_all

    do i = 1, size(strCoats)
       do j = 1, size(strCoats(i)%fppProcessorHandle)

          strCoats(i)%fppProcessorHandle(j) = fpp_create()
          call fpp_set_fpp_handle (strCoats(i)%fppProcessorHandle(j), &
               &                   fppObjectHandle, iret)
          if (iret /= 1) call reportFPPerror ('fpp_set_fpp_handle',iret,ierr)

          call fpp_set_data_msd (strCoats(i)%fppProcessorHandle(j), &
               &                 msdHandle, iret)
          if (iret /= 1) call reportFPPerror ('fpp_set_data_msd',iret,ierr)

          call fpp_set_element (strCoats(i)%fppProcessorHandle(j), &
               &                strCoats(i)%tensorRosette%id%userId, &
               &                strCoats(i)%results(j)%resultSet, &
               &                strCoats(i)%results(j)%matGroup, iret)
          if (iret /= 1) call reportFPPerror ('fpp_set_element',iret,ierr)

          call ffa_starttimer ('fpp_init')
          call fpp_init (strCoats(i)%fppProcessorHandle(j), iret)
          call ffa_stoptimer ('fpp_init')
          if (iret /= 1) call reportFPPerror ('fpp_init',iret,ierr)

          if (ierr /= 0) then
             call reportError (error_p,'Failure initializing fpp', &
                  'Strain coat element:'//GetId(strCoats(i)%tensorRosette%id), &
                  'fpp Handle ='//StrId(strCoats(i)%fppProcessorHandle(j)), &
                  addString='reported by fppInit')
             return
          end if

       end do
    end do

  end subroutine fppInit


  subroutine fppProcess (strCoats,iBuf,ierr)

    !!==========================================================================
    !! Calculate the fpp results.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 Apr 2008/1.0
    !!==========================================================================

    use IdTypeModule                , only : GetId, StrId
    use StrainCoatModule            , only : StrainCoatType, dp
    use ReportErrorModule           , only : reportError, error_p
    use FPPprocessorHandlerInterface, only : fpp_process
    use FFaProfilerInterface        , only : ffa_startTimer, ffa_stopTimer

    type(StrainCoatType), intent(inout) :: strCoats(:)
    integer             , intent(in)    :: iBuf
    integer             , intent(out)   :: ierr

    !! Local variables
    integer :: i, j, iret

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(strCoats)
       do j = 1, size(strCoats(i)%fppProcessorHandle)

          call ffa_starttimer ('fpp_process')
          call fpp_process (strCoats(i)%fppProcessorHandle(j), &
               &            strCoats(i)%results(j)%fppValue(1), &
               &            iBuf, 0.0_dp, 0.0_dp, 0, 0.0_dp, &
               &            0.0_dp, 0.0_dp, 0, iret)
          call ffa_stoptimer ('fpp_process')

          if (iret /= 1) then
             call reportFPPerror ('fpp_process',iret,ierr)
             call reportError (error_p,'Failure processing fpp', &
                  'Strain coat element:'//GetId(strCoats(i)%tensorRosette%id), &
                  'fpp Handle ='//StrId(strCoats(i)%fppProcessorHandle(j)), &
                  addString='reported by fppProcess')
             return
          end if

       end do
    end do

  end subroutine fppProcess


  subroutine fppProcessLast (strCoats,iBuf,biAxialGate,ierr)

    !!==========================================================================
    !! Calculate the fpp results and close the time loop
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 Apr 2008/1.0
    !!==========================================================================

    use IdTypeModule                , only : GetId, StrId
    use StrainCoatModule            , only : StrainCoatType, dp
    use ReportErrorModule           , only : reportError, error_p
    use FPPprocessorHandlerInterface, only : fpp_process
    use FFaProfilerInterface        , only : ffa_startTimer, ffa_stopTimer

    type(StrainCoatType), intent(inout) :: strCoats(:)
    integer             , intent(in)    :: iBuf
    real(dp)            , intent(in)    :: biAxialGate
    integer             , intent(out)   :: ierr

    !! Local variables
    integer :: i, j, iret, jBuf, kBuf

    !! --- Logic section ---

    ierr = 0
    jBuf = iBuf
    kBuf = 1
    do i = 1, size(strCoats)
       do j = 1, size(strCoats(i)%fppProcessorHandle)

          if (jBuf == 0) then
             !! Can't call fpp_process with nPnt = 0
             !! so process the last point twice
             jBuf = 1
             kBuf = size(strCoats(i)%results(j)%fppValue)
          end if

          call ffa_starttimer ('fpp_process')
          call fpp_process (strCoats(i)%fppProcessorHandle(j), &
               &            strCoats(i)%results(j)%fppValue(kBuf), jBuf, &
               &            strCoats(i)%results(j)%popAngle, &
               &            strCoats(i)%results(j)%angSpread, &
               &            strCoats(i)%results(j)%nBiaxial, &
               &            strCoats(i)%results(j)%biAxialSum, &
               &            strCoats(i)%results(j)%biAxialSqr, &
               &            biAxialGate, 1, iret)
          call ffa_stoptimer ('fpp_process')
          if (iret /= 1) then
             ierr = -1
             call reportFPPerror ('fpp_process',iret,ierr)
             call reportError (error_p,'Failure closing fpp', &
                  'Strain coat element:'//GetId(strCoats(i)%tensorRosette%id), &
                  'fpp Handle ='//StrId(strCoats(i)%fppProcessorHandle(j)), &
                  addString='reported by fppProcessLast')
             return
          end if

       end do
    end do

  end subroutine fppProcessLast

end module fppModule
