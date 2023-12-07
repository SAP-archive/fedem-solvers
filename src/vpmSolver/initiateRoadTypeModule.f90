!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module InitiateRoadTypeModule

  implicit none

  private

  public :: ReadRoads


contains

  subroutine ReadRoads (infp,functions,engines,roads,err)

    !!==========================================================================
    !! Initiates the road type with data from the solver input file.
    !!
    !! Programmer : Trond Arne Svidal                    date/rev : Aug 2002/1.0
    !!==========================================================================

    use KindModule        , only : dp, pi_p, lfnam_p
    use IdTypeModule      , only : ldesc_p, initId, ReportInputError
    use RoadTypeModule    , only : RoadType, DeallocateRoads, NullifyRoad
    use FunctionTypeModule, only : FunctionType, EngineType, GetPtrToId
    use inputUtilities    , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule    , only : lterm
    use reportErrorModule , only : AllocationError
    use reportErrorModule , only : reportError, debugFileOnly_p

    integer           , intent(in)  :: infp
    type(FunctionType), intent(in)  :: functions(:)
    type(EngineType)  , intent(in)  :: engines(:)
    type(RoadType)    , pointer     :: roads(:)
    integer           , intent(out) :: err

    !! Local variables
    integer :: idIn, nRoads, stat

    !! Define the ROAD namelist
    character(ldesc_p) :: extDescr
    character(lfnam_p) :: roadDataFileName
    integer            :: id, extId(10), roadFuncId
    integer            :: roadXengId, roadYengId, roadZengId
    real(dp)           :: Xoffset, Zshift, Theta
    character(len=30)  :: type
    logical            :: ThetaInRad
    namelist /ROAD/ id, extId, extDescr, type, roadFuncId, &
         &          roadXengId, roadYengId, roadZengId, &
         &          Xoffset, Zshift, Theta, ThetaInRad, roadDataFileName

    !! --- Logic section ---

    nRoads = iuGetNumberOfEntries(infp,'&ROAD',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &ROAD =',nRoads

    call DeallocateRoads (roads)
    allocate(roads(nRoads),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadRoads 10')
       return
    end if

    do idIn = 1,nRoads

       call NullifyRoad (roads(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&ROAD')) then
          err = err - 1
          call ReportInputError ('ROAD',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; type=''
       roadFuncId=0; roadXengId=0; roadYengId=0; roadZengId=0
       Xoffset=0.0_dp; Zshift=0.0_dp; Theta=0.0_dp; ThetaInRad=.true.
       roadDataFileName=''

       read(infp,nml=ROAD,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('ROAD',idIn)
          cycle
       end if

       call initId (roads(idIn)%id,id,extId,extDescr,stat)
       roads(idIn)%idIn = idIn

       if (roadFuncId > 0) then
          roads(idIn)%nRopar = 5
          allocate(roads(idIn)%roparr(roads(idIn)%nRopar),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadRoads 20')
             return
          end if
          !! Connect the specified function to this road
          if (ThetaInRad) theta = theta*180.0_dp/pi_p
          call ConnectRoadFunction (roadFuncId,Xoffset,Zshift,Theta,functions, &
               &                    roads(idIn),stat)
       else
          !! Connect the specified surface file to this road
          call ConnectRoadSurface (roadDataFileName,roads(idIn),stat)
       end if

       !! Connect the specified engines (if any) to this road
       call ConnectRoadEngines (roadXengId,roadYengId,roadZengId,engines, &
            &                   roads(idIn),stat)

       if (stat < 0) call ReportInputError ('ROAD',idIn,roads(idIn)%id)
       err = err + stat

    end do

    if (err < 0) call reportError (debugFileOnly_p,'ReadRoads')

  end subroutine ReadRoads


  subroutine ConnectRoadFunction (roadFuncId,X0,Z0,theta,functions,road,err)

    !!==========================================================================
    !! Connect functions to the given RoadType object.
    !!
    !! Programmer : Trond Arne Svidal                    date/rev : Aug 2002/1.0
    !!              Knut Morten Okstad                              Mar 2002/2.0
    !!==========================================================================

    use KindModule        , only : dp, pi_p
    use RoadTypeModule    , only : RoadType, FUNCTION_ROAD_p
    use FunctionTypeModule, only : FunctionType, GetPtrToId

    integer           , intent(in)    :: roadFuncId
    real(dp)          , intent(in)    :: X0, Z0, Theta
    type(FunctionType), intent(in)    :: functions(:)
    type(RoadType)    , intent(inout) :: road
    integer           , intent(inout) :: err

    !! --- Logic section ---

    err = 0
    road%type = FUNCTION_ROAD_p
    !! Connect the road function
    road%roadFunc => GetPtrToId(functions,roadFuncId)
    if (.not. associated(road%roadFunc)) err = err - 1

    !! Initialize road data array
    road%roparr(1) = 0.0_dp
    road%roparr(2) = cos(theta*pi_p/180.0_dp)
    road%roparr(3) = sin(theta*pi_p/180.0_dp)
    road%roparr(4) = -X0
    road%roparr(5) = Z0

  end subroutine ConnectRoadFunction


  subroutine ConnectRoadSurface (roadSurfaceFile,road,err)

    !!==========================================================================
    !! Connect a surface to the given RoadType object.
    !!
    !! Programmer : Trond Arne Svidal                    date/rev : Aug 2002/1.0
    !!==========================================================================

    use RoadTypeModule      , only : RoadType, SURFACE_ROAD_p
    use FFaFilePathInterface, only : FFa_checkPath

    character(len=*), intent(inout) :: roadSurfaceFile
    type(RoadType)  , intent(inout) :: road
    integer         , intent(inout) :: err

    !! --- Logic section ---

    road%type = SURFACE_ROAD_p
    if (len(roadSurfaceFile) > 0) then
       !! Connect the surface file
       call FFa_checkPath (roadSurfaceFile)
       road%RoadDataFileName = roadSurfaceFile
       road%nCharRoadDataFileName = len_trim(roadSurfaceFile)
    else
       err = err - 1
    end if

  end subroutine ConnectRoadSurface


  subroutine ConnectRoadEngines (rdXengId,rdYengId,rdZengId,engines,road,err)

    !!==========================================================================
    !! Connect engines to the given RoadType object.
    !!
    !! Programmer : Trond Arne Svidal                    date/rev : Sep 2002/1.0
    !!==========================================================================

    use RoadTypeModule    , only : RoadType
    use FunctionTypeModule, only : EngineType, GetPtrToId

    integer         , intent(in)    :: rdXengId, rdYengId, rdZengId
    type(EngineType), intent(in)    :: engines(:)
    type(RoadType)  , intent(inout) :: road
    integer         , intent(inout) :: err

    !! --- Logic section ---

    if (rdXengId > 0) then
       !! Connect the x engine
       road%xTrEng => GetPtrToId(engines,rdXengId)
       if (.not. associated(road%xTrEng)) err = err - 1
    else
       nullify(road%xTrEng)
    end if

    if (rdYengId > 0) then
       !! Connect the y engine
       road%yTrEng => GetPtrToId(engines,rdYengId)
       if (.not. associated(road%yTrEng)) err = err - 1
    else
       nullify(road%yTrEng)
    end if

    if (rdZengId > 0) then
       !! Connect the z engine
       road%zTrEng => GetPtrToId(engines,rdZengId)
       if (.not. associated(road%zTrEng)) err = err - 1
    else
       nullify(road%zTrEng)
    end if

  end subroutine ConnectRoadEngines

end module InitiateRoadTypeModule
