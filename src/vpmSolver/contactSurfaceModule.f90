!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file contactSurfaceModule.f90
!>
!> @brief Representation of contact surfaces.

!!==============================================================================
!> @brief Module with data types and subroutines representing contact surfaces.

module ContactSurfaceModule

  use IdTypeModule        , only : IdType
  use CurvePointTypeModule, only : CurvePointType, dp

  implicit none

  !> @brief Data type representing the glider curve in point-to-path joints.
  type GliderCurveType

     type(IdType)                  :: id         !< General identification data
     type(CurvePointType), pointer :: CPoints(:) !< Glider curve definition

     logical  :: isExtended !< Extended beyond the start and stop triads?
     logical  :: isLooping  !< Closed glider curve?
     real(dp) :: loopLength !< Total curve length in case of closed curve

  end type GliderCurveType


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdGliderCurve
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteGliderCurveType
  end interface


contains

  !!============================================================================
  !> @brief Returns pointer to (first) curve with specified ID.
  !>
  !> @param[in] array Array of contactsurfacemodule::glidercurvetype
  !>                  objects to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[out] index The array index of the found object
  !>
  !> @details If the glider curve is not found, NULL is returned.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Oct 2005

  function GetPtrToIdGliderCurve (array,id,index) result(ptr)
    type(GliderCurveType), pointer :: ptr

    integer              , intent(in)            :: id
    type(GliderCurveType), intent(in) , target   :: array(:)
    integer              , intent(out), optional :: index

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          if (present(index)) index = i
          return
       end if
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToIdGliderCurve returned nullified, baseId =',id

  end function GetPtrToIdGliderCurve


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] curve The contactsurfacemodule::glidercurvetype object to write
  !> @param[in] io File unit number to write to
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Oct 2005

  subroutine WriteGliderCurveType (curve,io)

    use IdTypeModule, only : writeId

    type(GliderCurveType), intent(in) :: curve
    integer              , intent(in) :: io

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'GliderCurve','{'
    call writeId (curve%id,io)

    if (associated(curve%CPoints)) then
       write(io,*) 'triads(id) = ', (curve%CPoints(i)%triad%id%baseId, &
            &                        i=1,size(curve%CPoints))
    end if

    write(io,*) 'isExtended? = ', curve%isExtended
    if (curve%isLooping) then
       write(io,*) 'loopLength  = ', curve%loopLength
    end if

    write(io,'(A)') '}'

  end subroutine WriteGliderCurveType


  !!============================================================================
  !> @brief Initializes a glider curve object.
  !>
  !> @param curve The contactsurfacemodule::glidercurvetype object to initialize
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Oct 2005

  pure subroutine NullifyGliderCurve (curve)

    use IdTypeModule, only : nullifyId

    type(GliderCurveType), intent(out) :: curve

    !! --- Logic section ---

    call nullifyId (curve%id)

    nullify(curve%CPoints)

    curve%isExtended = .false.
    curve%isLooping  = .false.
    curve%loopLength = 0.0_dp

  end subroutine NullifyGliderCurve


  !!============================================================================
  !> @brief Deallocates a glider curve object.
  !>
  !> @param curve The contactsurfacemodule::glidercurvetype object to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  elemental subroutine DeallocateGliderCurve (curve)

    use IdTypeModule, only : deallocateId

    type(GliderCurveType), intent(inout) :: curve

    !! --- Logic section ---

    call deallocateId (curve%id)

    if (associated(curve%CPoints)) deallocate(curve%CPoints)

    call nullifyGliderCurve (curve)

  end subroutine DeallocateGliderCurve


  !!============================================================================
  !> @brief Deallocates an array of glider curve objects.
  !>
  !> @param crvs The contactsurfacemodule::glidercurvetype objects to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 23 Jan 2017

  subroutine deallocateGliderCurves (crvs)
    type(GliderCurveType), pointer :: crvs(:)

    !! --- Logic section ---

    call DeallocateGliderCurve (crvs)

    deallocate(crvs)
    nullify(crvs)

  end subroutine deallocateGliderCurves


  !!============================================================================
  !> @brief Initializes the glider curves with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triads Array of all triads in the model
  !> @param[out] curves Array of all glider curve objects in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Oct 2005

  subroutine ReadGliderCurves (infp,triads,curves,err)

    use TriadTypeModule       , only : TriadType, dp, GetPtrToId
    use IdTypeModule          , only : IdType, ldesc_p
    use IdTypeModule          , only : initId, getId, ReportInputError
    use manipMatrixModule     , only : matmul34, invert34
    use inputUtilities        , only : iuGetNumberOfEntries
    use inputUtilities        , only : iuSetPosAtNextEntry
    use rotationModule        , only : orthonorm3
    use progressModule        , only : lterm
    use reportErrorModule     , only : reportError, debugFileOnly_p, error_p
    use reportErrorModule     , only : AllocationError
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer              , intent(in)  :: infp
    type(TriadType)      , intent(in)  :: triads(:)
    type(GliderCurveType), pointer     :: curves(:)
    integer              , intent(out) :: err

    !! Local variables
    integer                       :: i, j, idIn, nCurves, stat, lerr
    real(dp)                      :: PCG(3,4), secVec(3), phiQuart
    type(CurvePointType), pointer :: curveP1, curveP2

    !! Define the MASTER_CURVE namelist
    character(ldesc_p) :: extDescr
    integer, parameter :: MAX_MTRIAD = 1000
    integer            :: id, extId(10), nTriads, triadIds(MAX_MTRIAD)
    integer            :: isExtended, isLooping
    real(dp)           :: loopLength

    namelist /MASTER_CURVE/ extDescr, id, extId, nTriads, triadIds, &
         &                  isExtended, isLooping, loopLength

    !! Define the MASTER_POS namelist
    integer  :: masterId, triadId
    real(dp) :: PosInGlobal(4,3), curvature, slideVarVal, upVec(3)

    namelist /MASTER_POS/ masterId, triadId, PosInGlobal, &
         &                curvature, slideVarVal, upVec

    !! --- Logic section ---

    nCurves = iuGetNumberOfEntries(infp,'&MASTER_CURVE',err)
    if (err /= 0) goto 900
    write(lterm,*) 'Number of &MASTER_CURVE =',nCurves

    call DeallocateGliderCurves (curves)
    allocate(curves(nCurves),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadGliderCurves 10')
       return
    end if

    do idIn = 1, nCurves

       call NullifyGliderCurve (curves(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&MASTER_CURVE')) then
          err = err - 1
          call ReportInputError ('MASTER_CURVE',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       nTriads=0; triadIds=0; isExtended=0; isLooping=0; loopLength=0.0_dp

       read(infp,nml=MASTER_CURVE,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('MASTER_CURVE',idIn)
          cycle
       end if

       call initId (curves(idIn)%id,id,extId,extDescr,stat)
       curves(idIn)%isExtended = isExtended > 0 ! beyond beginning and end
       curves(idIn)%isLooping  = isLooping  > 0 ! i.e. circular

       if (curves(idIn)%isLooping) curves(idIn)%loopLength = loopLength

       !! Allocate array of curve definition points
       allocate(curves(idIn)%CPoints(nTriads),STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadGliderCurve 20')
          return
       end if

       !! Set pointer to glider triad for each curve point
       lerr = err
       do i = 1, nTriads
          curves(idIn)%CPoints(i)%triad => GetPtrToId(triads,triadIds(i))
          if (.not. associated(curves(idIn)%CPoints(i)%triad)) then
             err = err - 1
             call ReportInputError ('MASTER_CURVE',id=curves(idIn)%id, &
                  &                 msg='Invalid triad Id')
          end if
          !! Set samNodNum temporarily to flag that the triad should receive
          !! a proper node number later, also when the triad is grounded
          curves(idIn)%CPoints(i)%triad%samNodNum = 1
       end do
       if (err < lerr) cycle

       !! Read the curve orientation at each triad

       do i = 1, nTriads
          if (.not. iuSetPosAtNextEntry(infp,'&MASTER_POS')) then
             err = err - 1
             call ReportCurveError (curves(idIn)%id,i)
             cycle
          end if

          !! Default values
          masterId=0; triadId=0; PosInGlobal=0.0_dp
          curvature=0.0_dp; slideVarVal=0.0_dp; upVec=0.0_dp
          read(infp,nml=MASTER_POS,iostat=stat)
          if (stat /= 0) then
             err = err - 1
             call ReportCurveError (curves(idIn)%id,i)
             cycle
          else if (masterId /= id) then
             err = err - 1
             call ReportCurveError (curves(idIn)%id,triadId,masterId)
             cycle
          end if

          !! The MASTER_POS records on the input file are not necessarily in the
          !! same order as the curve point triads. Therefore we must search.
          nullify(curveP2)
          do j = 1, nTriads
             curveP1 => curves(idIn)%CPoints(j)
             if (triadId == curveP1%triad%id%baseId) then
                if (j < nTriads) then
                   curveP2 => curves(idIn)%CPoints(j+1)
                else
                   curveP2 => curves(idIn)%CPoints(1)
                end if
                exit
             end if
          end do
          if (.not. associated(curveP2)) then
             err = err - 1
             call ReportCurveError (curves(idIn)%id,triadId,masterId)
             cycle
          end if

          !! Position of curve point in Global system
          PCG = transpose(PosInGlobal) ! due to reading row by row
          call orthonorm3 (PCG(:,1:3)) ! due to possible round off

          !! Position of curve point in Triad system
          curveP1%CPosInT   = matmul34(invert34(curveP1%triad%ur),PCG)
          curveP1%slideVar  = slideVarVal
          curveP1%curvature = curvature
          curveP1%h0        = 0.0_dp
          curveP1%upVecInT1 = matmul(upVec,curveP1%triad%ur(:,1:3)) ! = R'*upVec
          curveP2%upVecInT2 = matmul(upVec,curveP2%triad%ur(:,1:3)) ! = R'*upVec

       end do
       if (err < lerr) cycle

       if (ffa_cmdlinearg_isTrue('constantArcHeight')) then

          !! Compute initial arc heights from secant
          do i = 1, nTriads
             curveP1 => curves(idIn)%CPoints(i)
             if (i < nTriads) then
                curveP2    => curves(idIn)%CPoints(i+1)
                slideVarVal = curveP2%slideVar - curveP1%slideVar
             else
                curveP2    => curves(idIn)%CPoints(1)
                slideVarVal = loopLength - curveP1%slideVar
             end if
             phiQuart   = slideVarVal * curveP1%curvature * 0.25_dp
             secVec     = matmul34(curveP2%triad%ur,curveP2%CPosInT(:,4)) &
                  &     - matmul34(curveP1%triad%ur,curveP1%CPosInT(:,4))
             curveP1%h0 = tan(phiQuart)*sqrt(dot_product(secVec,secVec))*0.5_dp
          end do

       end if

    end do

900 if (err < 0) call reportError (debugFileOnly_p,'ReadGliderCurves')

  contains

    !> @brief Produces an error message related to reading a namelist.
    subroutine ReportCurveError (id,tid,jid)
      type(IdType)    , intent(in) :: id
      integer         , intent(in) :: tid
      integer,optional, intent(in) :: jid
      character(len=128)           :: errMsg
      if (present(jid)) then
         write(errMsg,*) 'Invalid MASTER_POS specification for glider Curve', &
              &          trim(getId(id)),': GliderCurveId =',jid, &
              &          ' triadId =',tid
      else
         write(errMsg,*) 'Could not read MASTER_POS namelist number',tid, &
              &          'for glider Curve',trim(getId(id))
      end if
      call reportError (error_p,errMsg)
    end subroutine ReportCurveError

  end subroutine ReadGliderCurves

end module ContactSurfaceModule
