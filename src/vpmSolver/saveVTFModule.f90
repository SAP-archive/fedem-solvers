!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module saveVTFModule

  !!============================================================================
  !! This module contains an interface to the Ceetron VTF API for automatic
  !! export of superelement transformations from the Fedem solver to a VTF-file.
  !! The transformations are appended to an existing VTF-file which already
  !! should contain the model geometry.
  !!============================================================================

  implicit none

  public :: openVTF, closeVTF, writeSolverVTF

  private

  integer, save :: vtfFile = -999 ! VTF-file handle
  integer, save :: vtfInfo = -999 ! State info block handle
  integer, save :: nStep   = 0    ! Number of steps saved on the VTF-file

#ifdef FT_HAS_VTF
  interface

     !!=========================================================================
     !! All functions and subroutines in this interface are defined in the
     !! Ceetron VTF API library. However, since we are using F90, an explicit
     !! interface is preferred instead of the F77-style external declaractions
     !! found in the file VTFAPIFortran.f (see /devlibs/F3P_src/Ceetron/source
     !! and /devlibs/F3P_src/Ceetron/doc for the documentation).
     !!=========================================================================

     subroutine VTFFFileSetOutputDebugError (iFlag)
       integer, intent(in) :: iFlag
     end subroutine VTFFFileSetOutputDebugError

     function VTFFFileAppendFile (fName)
       character(len=*), intent(in) :: fName
       integer :: VTFFFileAppendFile
     end function VTFFFileAppendFile

     function VTFFFileCloseFile (iFile)
       integer, intent(in) :: iFile
       integer :: VTFFFileCloseFile
     end function VTFFFileCloseFile

     subroutine VTFFFileDelete (iFile)
       integer, intent(in) :: iFile
     end subroutine VTFFFileDelete

     function VTFFFileWriteBlock (iFile,iBlock)
       integer, intent(in) :: iFile, iBlock
       integer :: VTFFFileWriteBlock
     end function VTFFFileWriteBlock

     function VTFFTransBlockCreate (idBlock)
       integer, intent(in) :: idBlock
       integer :: VTFFTransBlockCreate
     end function VTFFTransBlockCreate

     function VTFFTransBlockSetResBlocks (iBlock,idBlocks,nBlocks,iStep)
       integer, intent(in) :: iBlock, idBlocks, nBlocks, iStep
       integer :: VTFFTransBlockSetResBlocks
     end function VTFFTransBlockSetResBlocks

     subroutine VTFFTransBlockSetName (iBlock,name)
       integer         , intent(in) :: iBlock
       character(len=*), intent(in) :: name
     end subroutine VTFFTransBlockSetName

     subroutine VTFFTransBlockDelete (iBlock)
       integer, intent(in) :: iBlock
     end subroutine VTFFTransBlockDelete

     function VTFFStateInfoBlockCreate()
       integer :: VTFFStateInfoBlockCreate
     end function VTFFStateInfoBlockCreate

     function VTFFStateInfoBlockSetStepData (iBlock,iStep,refVal,refType,name)
       use kindModule, only : sp
       integer         , intent(in) :: iBlock, iStep, refType
       real(sp)        , intent(in) :: refVal
       character(len=*), intent(in) :: name
       integer :: VTFFStateInfoBlockSetStepData
     end function VTFFStateInfoBlockSetStepData

     subroutine VTFFStateInfoBlockDelete (iBlock)
       integer, intent(in) :: iBlock
     end subroutine VTFFStateInfoBlockDelete

     function VTFFMatrixResBlockCreate (idBlock)
       integer, intent(in) :: idBlock
       integer :: VTFFMatrixResBlockCreate
     end function VTFFMatrixResBlockCreate

     function VTFFMatrixResBlockSetMatrix (iBlock,Tmat)
       use kindModule, only : sp
       integer , intent(in) :: iBlock
       real(sp), intent(in) :: Tmat
       integer :: VTFFMatrixResBlockSetMatrix
     end function VTFFMatrixResBlockSetMatrix

     subroutine VTFFMatrixResBlockSetMapToElemBlockID (iBlock,iBlockId)
       integer, intent(in) :: iBlock, iBlockId
     end subroutine VTFFMatrixResBlockSetMapToElemBlockID

     subroutine VTFFMatrixResBlockDelete (iBlock)
       integer, intent(in) :: iBlock
     end subroutine VTFFMatrixResBlockDelete

  end interface
#endif


contains

  subroutine openVTF (vtfName,ierr)

    !!==========================================================================
    !! Open the specified VTF-file for append.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jan 2006/1.0
    !!==========================================================================

#ifdef FT_HAS_VTF
    use reportErrorModule, only : warning_p, reportError
#endif

    character(len=*), intent(in)  :: vtfName
    integer         , intent(out) :: ierr

    !! --- Logic section ---

    ierr  = 0
    nStep = 0

#ifdef FT_HAS_VTF
    call VTFFFileSetOutputDebugError (1)
    vtfFile = VTFFFileAppendFile(trim(vtfName))
    if (vtfFile < 0) then
       ierr = -vtfFile
       call reportError (warning_p,'Could not open VTF-file '//vtfName, &
            'Output to VTF is therefore disabled in this simulation.')
       return
    end if
    vtfInfo = VTFFStateInfoBlockCreate()
#else
    print *,' ** openVTF dummy: ',trim(vtfName)
#endif

  end subroutine openVTF


  subroutine writeSolverVTF (sys,sups,ierr)

    !!==========================================================================
    !! Write response data for current time step to the VTF-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 9 Feg 2006/1.0
    !!==========================================================================

    use SystemTypeModule , only : SystemType
    use SupElTypeModule  , only : SupElType
    use reportErrorModule, only : debugFileOnly_p, reportError

    type(SystemType), intent(in)  :: sys
    type(SupElType) , intent(in)  :: sups(:)
    integer         , intent(out) :: ierr

    !! --- Logic section ---

    call writeStepVTF (sys,ierr)
    if (ierr < 0) goto 900

    call writeSupElTransformVTF (sups,ierr)
    if (ierr == 0) return

900 continue
    call reportError (debugFileOnly_p,'writeSolverVTF')

  end subroutine writeSolverVTF


  subroutine writeStepVTF (sys,ierr)

    !!==========================================================================
    !! Write time step info for current time step to the VTF-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 9 Feb 2006/1.0
    !!==========================================================================

    use SystemTypeModule , only : SystemType
#ifdef FT_HAS_VTF
    use kindModule       , only : sp
    use reportErrorModule, only : error_p, debugFileOnly_p, reportError
#endif

    type(SystemType), intent(in)  :: sys
    integer         , intent(out) :: ierr

    !! Local variables
#ifdef FT_HAS_VTF
    character(len=32) :: cn
#endif

    !! --- Logic section ---

    ierr = 0
    if (vtfFile < 0) return

    nStep = nStep + 1

#ifdef FT_HAS_VTF
    write(cn,"(I20)") sys%nStep
    cn = 'Step: '//adjustl(cn)
    ierr = VTFFStateInfoBlockSetStepData(vtfInfo,nStep,real(sys%time,sp),0,cn)-1
    if (ierr < 0) then
       call reportError (error_p,'Error defining state info block')
       call reportError (debugFileOnly_p,'writeStepVTF',ierr=ierr+1)
    end if
#else
    print *,' ** writeStepVTF dummy: ',sys%nStep
#endif

  end subroutine writeStepVTF


  subroutine writeSupElTransformVTF (sups,ierr)

    !!==========================================================================
    !! Write current superelement transformation matrices to the VTF-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jan 2006/1.0
    !!==========================================================================

    use SupElTypeModule  , only : SupElType
#ifdef FT_HAS_VTF
    use kindModule       , only : sp
    use reportErrorModule, only : error_p, debugFileOnly_p, reportError
#endif

    type(SupElType), intent(in)  :: sups(:)
    integer        , intent(out) :: ierr

    !! Local variables
#ifdef FT_HAS_VTF
    integer  :: i, iBlock, idBlock = 0
    real(sp) :: Tmat(3,4)
#endif

    !! --- Logic section ---

    ierr = 0
    if (vtfFile < 0) return

#ifdef FT_HAS_VTF
    do i = 1, size(sups)

       idBlock = idBlock + 1
       iBlock  = VTFFMatrixResBlockCreate(idBlock)
       if (iBlock < 0) then
          ierr = iBlock - 1
          call reportError (error_p,'Error creating matrix result block')
          goto 900
       end if

       Tmat = sups(i)%supTr ! cast to float
       ierr = VTFFMatrixResBlockSetMatrix(iBlock,Tmat(1,1)) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error defining matrix result block')
          goto 900
       end if

       call VTFFMatrixResBlockSetMapToElemBlockID (iBlock,sups(i)%id%baseId)

       ierr = VTFFFileWriteBlock(vtfFile,iBlock) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error writing matrix result block')
          goto 900
       end if

       call VTFFMatrixResBlockDelete (iBlock)

    end do
#else
    print *,' ** writeSupElTransformVTF dummy: ',size(sups)
#endif

    return

#ifdef FT_HAS_VTF
900 call reportError (debugFileOnly_p,'writeSupElTransformVTF',ierr=ierr+1)
#endif

  end subroutine writeSupElTransformVTF


  subroutine closeVTF (nParts,ierr)

    !!==========================================================================
    !! Write transformation block definition to the VTF-file and close the file.
    !! This subroutine is invoked after the time integration has terminated and
    !! the number of time steps are known. It is assumed that the matrix block
    !! numbers associated with step number "i" are given as (i-1)*nPart+1, ...
    !! i*nPart, where nPart equals the number of superelements in the model.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jan 2006/1.0
    !!==========================================================================

#ifdef FT_HAS_VTF
    use scratchArrayModule, only : getIntegerScratchArray
    use reportErrorModule , only : error_p, debugFileOnly_p, reportError
#endif

    integer, intent(in)  :: nParts
    integer, intent(out) :: ierr

    !! Local variables
#ifdef FT_HAS_VTF
    integer          :: i, iBlock
    integer, pointer :: idBlocks(:)
#endif

    !! --- Logic section ---

    ierr = 0
    if (vtfFile < 0) return

#ifdef FT_HAS_VTF
    idBlocks => getIntegerScratchArray(nParts,ierr)
    if (ierr < 0) goto 900

    !! Initialize the result block numbers for the first time step
    do i = 1, nParts
       idBlocks(i) = i
    end do

    !! Define and write the transformation block

    iBlock = VTFFTransBlockCreate(1)
    if (iBlock < 0) then
       ierr = iBlock - 1
       call reportError (error_p,'Error creating transformation block')
       goto 900
    end if

    call VTFFTransBlockSetName (iBlock,'Rigid body motion')

    do i = 1, nStep
       ierr = VTFFTransBlockSetResBlocks(iBlock,idBlocks(1),nParts,i) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error defining transformation block')
          goto 900
       end if
       idBlocks = idBlocks + nParts
    end do

    ierr = VTFFFileWriteBlock(vtfFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing transformation block')
       goto 900
    end if

    call VTFFTransBlockDelete (iBlock)

    !! Write the state info block

    ierr = VTFFFileWriteBlock(vtfFile,vtfInfo) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing state info block')
       goto 900
    end if

    call VTFFStateInfoBlockDelete (vtfInfo)

    !! Finally close the VTF-file

    ierr = VTFFFileCloseFile(vtfFile) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error closing VTF-file')
       goto 900
    end if

    call VTFFFileDelete (vtfFile)
#else
    print *,' ** closeVTF dummy: ',nParts
#endif

    vtfFile = -999
    vtfInfo = -999
    nStep = 0
    return

#ifdef FT_HAS_VTF
900 call reportError (debugFileOnly_p,'closeVTF',ierr=ierr+1)
#endif

  end subroutine closeVTF

end module saveVTFModule
