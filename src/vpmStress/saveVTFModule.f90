!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module saveVTFModule2

  !!============================================================================
  !! This module contains an interface to the Ceetron VTF API for automatic
  !! export of recovery results from Fedem to a VTF-file. The results are
  !! appended to an already existing VTF-file such that several superelements
  !! may be exported to the same file. Note that only results are covered by
  !! this module. The geometry of the superelement in question is assumed
  !! to have been exported to the specified VTF-file before the functions in
  !! this module are invoked.
  !!============================================================================

  use kindModule, only : sp, dp

  implicit none

#ifdef FT_HAS_VTF
  interface

     !!=========================================================================
     !! All functions and subroutines in this interface are defined in the
     !! Ceetron VTF API library. However, since we are using F90, an explicit
     !! interface is preferred instead of the F77-style external declarations
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

     function VTFFDispBlockCreate (idBlock)
       integer, intent(in) :: idBlock
       integer :: VTFFDispBlockCreate
     end function VTFFDispBlockCreate

     function VTFFDispBlockSetResBlocks (iBlock,idBlocks,nBlocks,iStep)
       integer, intent(in) :: iBlock, idBlocks, nBlocks, iStep
       integer :: VTFFDispBlockSetResBlocks
     end function VTFFDispBlockSetResBlocks

     subroutine VTFFDispBlockSetName (iBlock,name)
       integer         , intent(in) :: iBlock
       character(len=*), intent(in) :: name
     end subroutine VTFFDispBlockSetName

     subroutine VTFFDispBlockSetRelative (iBlock,relflag)
       integer, intent(in) :: iBlock, relflag
     end subroutine VTFFDispBlockSetRelative

     subroutine VTFFDispBlockDelete (iBlock)
       integer, intent(in) :: iBlock
     end subroutine VTFFDispBlockDelete

     function VTFFScalBlockCreate (idBlock)
       integer, intent(in) :: idBlock
       integer :: VTFFScalBlockCreate
     end function VTFFScalBlockCreate

     function VTFFScalBlockSetResBlocks (iBlock,idBlocks,nBlocks,iStep)
       integer, intent(in) :: iBlock, idBlocks, nBlocks, iStep
       integer :: VTFFScalBlockSetResBlocks
     end function VTFFScalBlockSetResBlocks

     subroutine VTFFScalBlockSetName (iBlock,name)
       integer         , intent(in) :: iBlock
       character(len=*), intent(in) :: name
     end subroutine VTFFScalBlockSetName

     subroutine VTFFScalBlockDelete (iBlock)
       integer, intent(in) :: iBlock
     end subroutine VTFFScalBlockDelete

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

     function VTFFResBlockCreate (idBlock,ndim,imaptyp,imapflag)
       integer, intent(in) :: idBlock, ndim, imaptyp, imapflag
       integer :: VTFFResBlockCreate
     end function VTFFResBlockCreate

     subroutine VTFFResBlockSetMapToBlockID (iBlock,isup)
       integer, intent(in) :: iBlock, isup
     end subroutine VTFFResBlockSetMapToBlockID

     function VTFFResBlockSetNumRes (iBlock,nRes)
       integer, intent(in) :: iBlock, nRes
       integer :: VTFFResBlockSetNumRes
     end function VTFFResBlockSetNumRes

     function VTFFResBlockAddRes (iBlock,fScalar)
       use kindModule, only : sp
       integer , intent(in) :: iBlock
       real(sp), intent(in) :: fScalar
       integer :: VTFFResBlockAddRes
     end function VTFFResBlockAddRes

     function VTFFResBlockSetRes3D (iBlock,res,nres)
       use kindModule, only : sp
       integer , intent(in) :: iBlock, nres
       real(sp), intent(in) :: res
       integer :: VTFFResBlockSetRes3D
     end function VTFFResBlockSetRes3D

     subroutine VTFFResBlockDelete (iBlock)
       integer, intent(in) :: iBlock
     end subroutine VTFFResBlockDelete

  end interface
#endif

  interface writeFooterVTF
     module procedure writeFooterStress
     module procedure writeFooterModes
     module procedure writeFooterMode
  end interface

  private :: writeFooterStress, writeFooterModes, writeFooterMode
#ifdef FT_HAS_VTF
  private :: StrMode
#endif


  !! These parameters are taken from VTFAPIFortran.f
  integer , parameter :: VTFA_RESMAP_NODE         = 0
  integer , parameter :: VTFA_RESMAP_ELEMENT      = 1
  integer , parameter :: VTFA_RESMAP_ELEMENT_NODE = 3
  real(sp), parameter :: VTFA_UNDEFINED           = 3.4e38_sp

  !! Result set options for multi-surface results
  integer , parameter :: VTF_RESSET_ABSMAX  = -1
  integer , parameter :: VTF_RESSET_ABSMIN  = -2
  integer , parameter :: VTF_RESSET_AVERAGE = -3
  integer , parameter :: VTF_RESSET_MAX     = -4
  integer , parameter :: VTF_RESSET_MAXDIFF = -5
  integer , parameter :: VTF_RESSET_MIN     = -6


contains

#ifdef FT_HAS_VTF

  subroutine writeTimeStepVTF (iBlock,iStep,time,ierr)

    !!==========================================================================
    !! Write time step information of current time step to the VTF-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 9 Feb 2006/1.0
    !!==========================================================================

    use reportErrorModule, only : error_p, debugFileOnly_p, reportError

    integer , intent(in)  :: iBlock, iStep
    real(dp), intent(in)  :: time
    integer , intent(out) :: ierr

    !! Local variables
    integer, save     :: nStep = 0
    character(len=16) :: cStep

    !! --- Logic section ---

    ierr = 0
    if (iBlock < 0) return

    nStep = nStep + 1
    if (iStep == 0) then
       cStep = 'Initial state'
    else
       write(cStep,"('Step:',I6)") iStep
    end if
    ierr = VTFFStateInfoBlockSetStepData(iBlock,nStep,real(time,sp),0,cStep) - 1
    if (ierr == 0) return

    call reportError (error_p,'Error defining state info block')
    call reportError (debugFileOnly_p,'writeTimeStepVTF',ierr=ierr+1)

  end subroutine writeTimeStepVTF


  subroutine writeSupElTransformVTF (iFile,idBlock,supEl,ierr)

    !!==========================================================================
    !! Write current superelement transformation matrix to the VTF-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jan 2006/1.0
    !!==========================================================================

    use SupElTypeModule  , only : SupElType
    use reportErrorModule, only : error_p, debugFileOnly_p, reportError

    integer        , intent(in)  :: iFile, idBlock
    type(SupElType), intent(in)  :: supEl
    integer        , intent(out) :: ierr

    !! Local variables
    integer  :: iBlock
    real(sp) :: Tmat(3,4)

    !! --- Logic section ---

    iBlock = VTFFMatrixResBlockCreate(idBlock)
    if (iBlock < 0) then
       ierr = iBlock - 1
       call reportError (error_p,'Error creating matrix result block')
       goto 900
    end if

    Tmat = supEl%supTr ! cast to float
    ierr = VTFFMatrixResBlockSetMatrix(iBlock,Tmat(1,1)) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error defining matrix result block')
       goto 900
    end if

    call VTFFMatrixResBlockSetMapToElemBlockID (iBlock,supEl%id%baseId)

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing matrix result block')
       goto 900
    end if

    call VTFFMatrixResBlockDelete (iBlock)

    return

900 continue
    call reportError (debugFileOnly_p,'writeSupElTransformVTF',ierr=ierr+1)

  end subroutine writeSupElTransformVTF


  subroutine writeDisplacementVTF (iFile,idBlock,supEl,sam,sv,dscale,ierr)

    !!==========================================================================
    !! Write expanded superelement deformations to the VTF-file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jan 2006/1.0
    !!==========================================================================

    use SamModule         , only : SamType
    use SupElTypeModule   , only : SupElType
    use scratchArrayModule, only : getIntegerScratchArray, getSingleScratchArray
    use reportErrorModule , only : error_p, debugFileOnly_p, reportError

    integer        , intent(in)  :: iFile, idBlock
    type(SupElType), intent(in)  :: supEl
    type(SamType)  , intent(in)  :: sam
    real(dp)       , intent(in)  :: sv(:), dscale
    integer        , intent(out) :: ierr

    !! Local variables
    integer           :: i, j, inod, iBlock
    real(sp), pointer :: disVec(:)

    !! --- Logic section ---

    disVec => getSingleScratchArray(3*sam%nnod,ierr)
    if (ierr < 0) goto 900

    iBlock = VTFFResBlockCreate(idBlock,3,VTFA_RESMAP_NODE,0)
    if (iBlock < 0) then
       ierr = iBlock - 1
       call reportError (error_p,'Error creating displacement block')
       goto 900
    end if

    inod = 0
    do i = 1, sam%nnod
       j = sam%madof(i)
       if (sam%minex(i) > 0) then
          inod = inod + 1
          disVec(3*inod-2:3*inod) = sv(j:j+2) * dscale ! cast to float
       end if
    end do
    if (inod > 0) then
       ierr = VTFFResBlockSetRes3D(iBlock,disVec(1),inod) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error defining displacement block')
          goto 900
       end if
    end if

    call VTFFResBlockSetMapToBlockID (iBlock,supEl%id%baseId)

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing displacement block')
       goto 900
    end if

    call VTFFResBlockDelete (iBlock)

    return

900 continue
    call reportError (debugFileOnly_p,'writeDisplacementVTF',ierr=ierr+1)

  end subroutine writeDisplacementVTF


  subroutine writeStressVTF (iBlock,nenod,nrset,irset,irmap,stress,ierr)

    !!==========================================================================
    !! Write stress results to the VTF-file for current element. For multi-
    !! surface results, the actual result value written is selected based on
    !! the value of the irset argument. The available choices correspond to
    !! the "Result Set" pull-down menus in the "Animation" property view.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jan 2006/1.0
    !!==========================================================================

    use reportErrorModule, only : error_p, debugFileOnly_p, reportError

    integer , intent(in)  :: iBlock, nenod, nrset, irset, irmap
    real(dp), intent(in)  :: stress(nenod,nrset)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, j
    real(sp) :: sval, smin, smax, savg
    !!                                  1  2  3  4  5  6  7  8  9 10
    integer, parameter :: T10m(10) = (/ 1, 3, 5,10, 2, 4, 6, 7, 8, 9 /)
    integer, parameter :: P15m(15) = (/ 1, 3, 5,10,12,14, 2, 4, 6,11, &
         &                             13,15, 7, 8, 9/)
    integer, parameter :: H20m(20) = (/ 1, 3, 5, 7,13,15,17,19, 2, 4, &
         &                              6, 8,14,16,18,20, 9,10,11,12 /)

    !! --- Logic section ---

    savg = 0.0_sp
    do i = 1, nenod
       if (nenod == 10) then
          j = T10m(i)
       else if (nenod == 15) then
          j = P15m(i)
       else if (nenod == 20) then
          j = H20m(i)
       else
          j = i
       end if
       if (nrset <= 1) then ! Basic (solid stress, or single-surface shell)

          sval = stress(j,1)

       else if (irset > 0 .and. irset <= nrset) then ! Top/bottom shell surface

          sval = stress(j,irset)

       else if (irset == VTF_RESSET_ABSMAX) then

          smin = minval(stress(j,:))
          smax = maxval(stress(j,:))
          if (-smin > smax) then
             sval = smin
          else
             sval = smax
          end if

       else if (irset == VTF_RESSET_ABSMIN) then

          smin = minval(stress(j,:))
          smax = maxval(stress(j,:))
          if (-smin < smax) then
             sval = smin
          else
             sval = smax
          end if

       else if (irset == VTF_RESSET_AVERAGE) then

          sval = sum(stress(j,:)) / real(nrset,sp)

       else if (irset == VTF_RESSET_MAX) then

          sval = maxval(stress(j,:))

       else if (irset == VTF_RESSET_MAXDIFF) then

          sval = maxval(stress(j,:)) - minval(stress(j,:))

       else if (irset == VTF_RESSET_MIN) then

          sval = minval(stress(j,:))

       else
          sval = VTFA_UNDEFINED
       end if

       if (irmap == VTFA_RESMAP_ELEMENT_NODE) then
          ierr = VTFFResBlockAddRes(iBlock,sval) - 1
          if (ierr < 0) goto 900
       else if (sval /= VTFA_UNDEFINED) then
          savg = savg + sval / real(nenod,sp)
       else
          savg = sval
          exit
       end if

    end do

    if (irmap == VTFA_RESMAP_ELEMENT) then
       ierr = VTFFResBlockAddRes(iBlock,savg) - 1
       if (ierr < 0) goto 900
    end if

    return

900 continue
    call reportError (error_p,'Error defining scalar result block')
    call reportError (debugFileOnly_p,'writeStressVTF',ierr=ierr+1)

  end subroutine writeStressVTF


  subroutine writeDummyStressVTF (iFile,idBlock,iSup,ierr)

    !!==========================================================================
    !! Write a dummy (empty) scalar block for the given superelement.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 07 Mar 2006/1.0
    !!==========================================================================

    use reportErrorModule, only : error_p, reportError

    integer, intent(in)  :: iFile, idBlock, iSup
    integer, intent(out) :: ierr

    !! Local variables
    integer :: iBlock

    !! --- Logic section ---

    iBlock = VTFFResBlockCreate(idBlock,1,VTFA_RESMAP_ELEMENT,0)
    if (iBlock < 0) then
       ierr = iBlock
       call reportError (error_p,'Error creating scalar result block', &
            &            addString='writeDummyStressVTF',ierr=ierr)
       return
    end if

    call VTFFResBlockSetMapToBlockID (iBlock,isup)
    ierr = VTFFResBlockSetNumRes(iBlock,1) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error allocating scalar result block', &
            &            addString='writeDummyStressVTF',ierr=ierr+1)
       return
    end if

    ierr = VTFFResBlockAddRes(iBlock,VTFA_UNDEFINED) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error defining scalar result block', &
            &            addString='writeDummyStressVTF',ierr=ierr+1)
       return
    end if

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing scalar result block', &
            &            addString='writeDummyStressVTF',ierr=ierr+1)
       return
    end if

    call VTFFResBlockDelete (iBlock)

  end subroutine writeDummyStressVTF


  subroutine writeFooterStress (iFile,iInfo,nParts,idOffset,nStep, &
       &                        lVTFAt0,lDeformation,sclName,ierr)

    !!==========================================================================
    !! Write result block definitions to the VTF-file. This subroutine is
    !! invoked only after the last superelement (part) has been processed.
    !! It is assumed that all the previous parts have the same amount of blocks
    !! and that the first block number is equal to 1 + idOffset*(i-1)/(nPart-1)
    !! for part number "i".
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Jan 2006/1.0
    !!==========================================================================

    use scratchArrayModule, only : getIntegerScratchArray
    use reportErrorModule , only : error_p, debugFileOnly_p, reportError

    integer         , intent(in)  :: iFile, iInfo, nParts, idOffset, nStep
    logical         , intent(in)  :: lVTFAt0, lDeformation
    character(len=*), intent(in)  :: sclName
    integer         , intent(out) :: ierr

    !! Local variables
    integer          :: i, iinc, iStep, iBlock
    integer, pointer :: idBlocks(:)

    !! --- Logic section ---

    idBlocks => getIntegerScratchArray(nParts,ierr)
    if (ierr < 0) goto 900

    !! Initialize the result block numbers for the first time step
    idBlocks(1) = 1
    if (nParts > 1) then
       iinc = idOffset/(nParts-1)
       do i = 2, nParts
          idBlocks(i) = idBlocks(i-1) + iinc
       end do
    end if

    !! Define and write the transformation block

    iBlock = VTFFTransBlockCreate(1)
    if (iBlock < 0) then
       ierr = iBlock
       call reportError (error_p,'Error creating transformation block')
       goto 900
    end if

    call VTFFTransBlockSetName (iBlock,'Rigid body motion')

    do iStep = 1, nStep
       ierr = VTFFTransBlockSetResBlocks(iBlock,idBlocks(1),nParts,iStep) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error defining transformation block')
          goto 900
       end if
       idBlocks = idBlocks + 1
    end do
    idBlocks = idBlocks - nStep

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing transformation block')
       goto 900
    end if

    call VTFFTransBlockDelete (iBlock)

    if (lDeformation) then

       !! Define and write the deformation block

       iBlock = VTFFDispBlockCreate(1)
       if (iBlock < 0) then
          ierr = iBlock - 1
          call reportError (error_p,'Error creating displacement block')
          goto 900
       end if

       call VTFFDispBlockSetName (iBlock,'Deformation')
       call VTFFDispBlockSetRelative (iBlock,1)

       if (lVTFAt0) then
          !! The initial (residual) state does not have any deformation
          iinc = 2
          idBlocks = idBlocks + 2
       else
          iinc = 1
       end if
       do iStep = iinc, nStep
          ierr = VTFFDispBlockSetResBlocks(iBlock,idBlocks(1),nParts,iStep) - 1
          if (ierr < 0) then
             call reportError (error_p,'Error defining displacement block')
             goto 900
          end if
          idBlocks = idBlocks + 2
       end do
       idBlocks = idBlocks - 2*nStep + 1

       ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error writing displacement block')
          goto 900
       end if

       call VTFFDispBlockDelete (iBlock)

    else
       idBlocks = idBlocks + 1
    end if

    !! Define and write the scalar result block

    iBlock = VTFFScalBlockCreate(1)
    if (iBlock < 0) then
       ierr = iBlock - 1
       call reportError (error_p,'Error creating scalar result block')
       goto 900
    end if

    call VTFFScalBlockSetName (iBlock,sclName)

    do iStep = 1, nStep
       ierr = VTFFScalBlockSetResBlocks(iBlock,idBlocks(1),nParts,iStep) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error defining scalar result block')
          goto 900
       end if
       idBlocks = idBlocks + 2
    end do

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing scalar result block')
       goto 900
    end if

    call VTFFScalBlockDelete (iBlock)

    if (iInfo < 0) return

    !! Write the state info block

    ierr = VTFFFileWriteBlock(iFile,iInfo) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing state info block')
       goto 900
    end if

    call VTFFStateInfoBlockDelete (iInfo)

    return

900 continue
    call reportError (debugFileOnly_p,'writeFooterVTF',ierr=ierr+1)

  end subroutine writeFooterStress


  subroutine writeFooterModes (iFile,iInfo,idOffset,nParts,nStep,nMode, &
       &                       modeNum,modeIds,freqs,ierr)

    !!==========================================================================
    !! Write result block definitions to the VTF-file. This subroutine is
    !! invoked only after the last superelement (part) has been processed.
    !! It is assumed that all the previous parts have the same amount of blocks
    !! and that the first block number is equal to 1 + idOffset*(i-1)/(nPart-1)
    !! for part number "i".
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 31 Jan 2006/1.0
    !!==========================================================================

    use scratchArrayModule, only : getIntegerScratchArray
    use reportErrorModule , only : error_p, debugFileOnly_p, reportError

    integer , intent(in)  :: iFile, iInfo, idOffset, nParts, nStep, nMode
    integer , intent(in)  :: modeNum(:), modeIds(:,:)
    real(dp), intent(in)  :: freqs(:)
    integer , intent(out) :: ierr

    !! Local variables
    character(len=8) :: chMode
    integer          :: i, iinc, iStep, iMode, iBlock
    integer, pointer :: idBlocks(:)

    !! --- Logic section ---

    idBlocks => getIntegerScratchArray(nParts,ierr)
    if (ierr < 0) goto 900

    !! Initialize the result block numbers for the first time step

    idBlocks(1) = 1
    if (nParts > 1) then
       iinc = idOffset/(nParts-1)
       do i = 2, nParts
          idBlocks(i) = idBlocks(i-1) + iinc
       end do
    end if

    !! Define and write the transformation block

    iBlock = VTFFTransBlockCreate(1)
    if (iBlock < 0) then
       ierr = iBlock
       call reportError (error_p,'Error creating transformation block')
       goto 900
    end if

    call VTFFTransBlockSetName (iBlock,'Rigid body motion')

    do iStep = 1, nStep
       ierr = VTFFTransBlockSetResBlocks(iBlock,idBlocks(1),nParts,iStep) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error defining transformation block')
          goto 900
       end if
       idBlocks = idBlocks + 1
    end do

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing transformation block')
       goto 900
    end if

    call VTFFTransBlockDelete (iBlock)

    !! Define and write a deformation block for the dynamic response

    iBlock = VTFFDispBlockCreate(1)
    if (iBlock < 0) then
       ierr = iBlock - 1
       call reportError (error_p,'Error creating displacement block')
       goto 900
    end if

    call VTFFDispBlockSetName (iBlock,'Dynamic response')
    call VTFFDispBlockSetRelative (iBlock,1)

    idBlocks = idBlocks - nStep
    do iStep = 1, nStep
       ierr = VTFFDispBlockSetResBlocks(iBlock,idBlocks(1),nParts,iStep) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error defining displacement block')
          goto 900
       end if
       idBlocks = idBlocks + 1 + count(modeIds(iStep,1:nMode) > 0)
    end do

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing displacement block')
       goto 900
    end if

    call VTFFDispBlockDelete (iBlock)

    !! Define and write a deformation block for each expanded mode shape

    do iMode = 1, nMode
       if (.not. any(modeIds(:,iMode) > 0)) cycle

       iBlock = VTFFDispBlockCreate(1+iMode)
       if (iBlock < 0) then
          ierr = iBlock - 1
          call reportError (error_p,'Error creating displacement block')
          goto 900
       end if

       call VTFFDispBlockSetName (iBlock,StrMode(modeNum(iMode),freqs(iMode)))
       call VTFFDispBlockSetRelative (iBlock,1)

       do iStep = 1, nStep
          if (modeIds(iStep,iMode) == 0) cycle

          idBlocks(1) = modeIds(iStep,iMode) - idOffset
          if (nParts > 1) then
             iinc = idOffset/(nParts-1)
             do i = 2, nParts
                idBlocks(i) = idBlocks(i-1) + iinc
             end do
          end if

          ierr = VTFFDispBlockSetResBlocks(iBlock,idBlocks(1),nParts,iStep) - 1
          if (ierr < 0) then
             call reportError (error_p,'Error defining displacement block')
             goto 900
          end if

       end do

       ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error writing displacement block')
          goto 900
       end if

       call VTFFDispBlockDelete (iBlock)

    end do

    if (iInfo < 0) return

    !! Write the state info block

    ierr = VTFFFileWriteBlock(iFile,iInfo) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing state info block')
       goto 900
    end if

    call VTFFStateInfoBlockDelete (iInfo)

    return

900 continue
    call reportError (debugFileOnly_p,'writeFooterVTF',ierr=ierr+1)

  end subroutine writeFooterModes


  subroutine writeFooterMode (iFile,idOffset,nParts,modeNum,modeId,freq,ierr)

    !!==========================================================================
    !! Write result block definitions to the VTF-file. This subroutine is
    !! invoked only after the last superelement (part) has been processed.
    !! It is assumed that all the previous parts have the same amount of blocks
    !! and that the first block number is equal to 1 + idOffset*(i-1)/(nPart-1)
    !! for part number "i".
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 31 Jan 2006/1.0
    !!==========================================================================

    use scratchArrayModule, only : getIntegerScratchArray
    use reportErrorModule , only : error_p, debugFileOnly_p, reportError

    integer , intent(in)  :: iFile, idOffset, nParts, modeNum, modeId
    real(dp), intent(in)  :: freq
    integer , intent(out) :: ierr

    !! Local variables
    integer          :: i, iinc, iBlock
    integer, pointer :: idBlocks(:)

    !! --- Logic section ---

    idBlocks => getIntegerScratchArray(nParts,ierr)
    if (ierr < 0) goto 900

    !! Initialize the result block numbers

    idBlocks(1) = 1
    if (nParts > 1) then
       iinc = idOffset/(nParts-1)
       do i = 2, nParts
          idBlocks(i) = idBlocks(i-1) + iinc
       end do
    end if

    !! Define and write the transformation block

    iBlock = VTFFTransBlockCreate(1)
    if (iBlock < 0) then
       ierr = iBlock
       call reportError (error_p,'Error creating transformation block')
       goto 900
    end if

    call VTFFTransBlockSetName (iBlock,'Rigid body motion')

    ierr = VTFFTransBlockSetResBlocks(iBlock,idBlocks(1),nParts,1) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error defining transformation block')
       goto 900
    end if

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing transformation block')
       goto 900
    end if

    call VTFFTransBlockDelete (iBlock)

    if (modeId <= 0) return

    !! Define and write a deformation block for the expanded mode shape

    iBlock = VTFFDispBlockCreate(1)
    if (iBlock < 0) then
       ierr = iBlock - 1
       call reportError (error_p,'Error creating displacement block')
       goto 900
    end if

    call VTFFDispBlockSetName (iBlock,StrMode(modeNum,freq))
    call VTFFDispBlockSetRelative (iBlock,1)

    ierr = VTFFDispBlockSetResBlocks(iBlock,idBlocks(1),nParts,1) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error defining displacement block')
       goto 900
    end if

    ierr = VTFFFileWriteBlock(iFile,iBlock) - 1
    if (ierr < 0) then
       call reportError (error_p,'Error writing displacement block')
       goto 900
    end if

    call VTFFDispBlockDelete (iBlock)

    return

900 continue
    call reportError (debugFileOnly_p,'writeFooterVTF',ierr=ierr+1)

  end subroutine writeFooterMode


  function StrMode (mode,freq)

    !!==========================================================================
    !! Returns a string identifying the given mode number and eigenfrequency.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 20 Sep 2006/1.0
    !!==========================================================================

    integer , intent(in) :: mode
    real(dp), intent(in) :: freq
    character(len=24)    :: StrMode

    !! --- Logic section ---

    if (freq <= 0.0_dp) then
       write(StrMode,"('Mode',i4)") mode
    else if (freq < 10.0_dp) then
       write(StrMode,"('Mode',i4,': ',F7.4,' Hz')") mode,freq
    else if (freq < 100.0_dp) then
       write(StrMode,"('Mode',i4,': ',F7.3,' Hz')") mode,freq
    else if (freq < 1000.0_dp) then
       write(StrMode,"('Mode',i4,': ',F7.2,' Hz')") mode,freq
    else if (freq < 10000.0_dp) then
       write(StrMode,"('Mode',i4,': ',F7.1,' Hz')") mode,freq
    else if (freq < 100000.0_dp) then
       write(StrMode,"('Mode',i4,': ',I6,' Hz')") mode,nint(freq)
    else
       write(StrMode,"('Mode',i4,': ',1PE8.1,' Hz')") mode,freq
    end if

  end function StrMode

#else

  !!============================================================================
  !! Dummy implementation of the VTF interface.
  !! This is used on platforms where result export to VTF is not supported.
  !!============================================================================

  subroutine VTFFFileSetOutputDebugError (iFlag)
    integer, intent(in) :: iFlag
    print *,'VTFFileSetOutputDebugError dummy: ',iFlag
  end subroutine VTFFFileSetOutputDebugError

  function VTFFFileAppendFile (fName)
    use reportErrorModule, only : note_p, reportError
    character(len=*), intent(in) :: fName
    integer :: VTFFFileAppendFile
    VTFFFileAppendFile = 0
    call reportError (note_p,'Export to VTF is not available in this version',&
         &            addString=trim(fName)//' not written')
  end function VTFFFileAppendFile

  function VTFFFileCloseFile (iFile)
    integer, intent(in) :: iFile
    integer :: VTFFFileCloseFile
    VTFFFileCloseFile = 1
    print *,'VTFFileCloseFile dummy: ',iFile
  end function VTFFFileCloseFile

  subroutine VTFFFileDelete (iFile)
    integer, intent(in) :: iFile
    print *,'VTFFileDelete dummy: ',iFile
  end subroutine VTFFFileDelete

  function VTFFFileWriteBlock (iFile,iBlock)
    integer, intent(in) :: iFile, iBlock
    integer :: VTFFFileWriteBlock
    VTFFFileWriteBlock = 1
    print *,'VTFFileWriteBlock dummy: ',iFile,iBlock
  end function VTFFFileWriteBlock

  function VTFFStateInfoBlockCreate()
    integer :: VTFFStateInfoBlockCreate
    VTFFStateInfoBlockCreate = 0
    print *,'VTFStateInfoBlockCreate dummy:'
  end function VTFFStateInfoBlockCreate

  function VTFFResBlockCreate (idBlock,ndim,imaptyp,imapflag)
    integer, intent(in) :: idBlock, ndim, imaptyp, imapflag
    integer :: VTFFResBlockCreate
    VTFFResBlockCreate = 1
    print *,'VTFResBlockCreate dummy: ',idBlock,ndim,imaptyp,imapflag
  end function VTFFResBlockCreate

  subroutine VTFFResBlockSetMapToBlockID (iBlock,isup)
    integer, intent(in) :: iBlock, isup
    print *,'VTFResBlockSetMapToBlockID dummy: ',iBlock,isup
  end subroutine VTFFResBlockSetMapToBlockID

  function VTFFResBlockSetNumRes (iBlock,nRes)
    integer, intent(in) :: iBlock, nRes
    integer :: VTFFResBlockSetNumRes
    VTFFResBlockSetNumRes = 1
    print *,'VTFResBlockSetNumRes dummy: ',iBlock,nRes
  end function VTFFResBlockSetNumRes

  subroutine VTFFResBlockDelete (iBlock)
    integer, intent(in) :: iBlock
    print *,'VTFResBlockDelete dummy: ',iBlock
  end subroutine VTFFResBlockDelete

  subroutine writeTimeStepVTF (iBlock,iStep,time,ierr)
    integer , intent(in)  :: iBlock, iStep
    real(dp), intent(in)  :: time
    integer , intent(out) :: ierr
    ierr = 0
    print *,'writeTimeStepVTF dummy: ',iBlock,iStep,time
  end subroutine writeTimeStepVTF

  subroutine writeSupElTransformVTF (iFile,idBlock,supEl,ierr)
    use SupElTypeModule, only : SupElType
    integer        , intent(in)  :: iFile, idBlock
    type(SupElType), intent(in)  :: supEl
    integer        , intent(out) :: ierr
    ierr = 0
    print *,'writeSupElTransformVTF dummy: ',iFile,idBlock,supEl%id%userId
  end subroutine writeSupElTransformVTF

  subroutine writeDisplacementVTF (iFile,idBlock,supEl,sam,sv,dscale,ierr)
    use SamModule      , only : SamType
    use SupElTypeModule, only : SupElType
    integer        , intent(in)  :: iFile, idBlock
    type(SupElType), intent(in)  :: supEl
    type(SamType)  , intent(in)  :: sam
    real(dp)       , intent(in)  :: sv(:), dscale
    integer        , intent(out) :: ierr
    ierr = 0
    print *,'writeDisplacementVTF dummy: ',iFile,idBlock,supEl%id%userId, &
         &                                 sam%nel,size(sv),dscale
  end subroutine writeDisplacementVTF

  subroutine writeStressVTF (iBlock,nenod,nrset,irset,irmap,sigma,ierr)
    integer , intent(in)  :: iBlock, nenod, nrset, irset, irmap
    real(dp), intent(in)  :: sigma(nenod,nrset)
    integer , intent(out) :: ierr
    ierr = 0
    print *,'writeStressVTF dummy: ',iBlock,nenod,nrset,irset,irmap,size(sigma)
  end subroutine writeStressVTF

  subroutine writeDummyStressVTF (iFile,idBlock,iSup,ierr)
    integer, intent(in)  :: iFile, idBlock, iSup
    integer, intent(out) :: ierr
    ierr = 0
    print *,'writeDummyStressVTF dummy: ',iFile,idBlock,iSup
  end subroutine writeDummyStressVTF

  subroutine writeFooterStress (iFile,iInfo,nParts,idOffset,nStep, &
       &                        lVTFAt0,lDeformation,sclName,ierr)
    integer         , intent(in)  :: iFile, iInfo, nParts, idOffset, nStep
    logical         , intent(in)  :: lVTFAt0, lDeformation
    character(len=*), intent(in)  :: sclName
    integer         , intent(out) :: ierr
    ierr = 0
    print *,'writeFooterStress dummy: ',iFile,iInfo,nParts,idOffset,nStep, &
         &                              lVTFAt0,lDeformation,sclName
  end subroutine writeFooterStress

  subroutine writeFooterModes (iFile,iInfo,idOffset,nParts,nStep,nMode, &
       &                       modeNum,modeIds,freqs,ierr)
    integer , intent(in)  :: iFile, iInfo, idOffset, nParts, nStep, nMode
    integer , intent(in)  :: modeNum(:), modeIds(:,:)
    real(dp), intent(in)  :: freqs(:)
    integer , intent(out) :: ierr
    ierr = 0
    print *,'writeFooterModes dummy: ',iFile,iInfo,idOffset, &
         &                             nParts,nStep,nMode, &
         &                             size(modeNum),size(modeIds),size(freqs)
  end subroutine writeFooterModes

  subroutine writeFooterMode (iFile,idOffset,nParts,modeNum,modeId,freq,ierr)
    integer , intent(in)  :: iFile, idOffset, nParts, modeNum, modeId
    real(dp), intent(in)  :: freq
    integer , intent(out) :: ierr
    ierr = 0
    print *,'writeFooterMode dummy: ',iFile,idOffset,nParts,modeNum,modeId,freq
  end subroutine writeFooterMode

#endif

end module saveVTFModule2
