!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module RoadRoutinesModule

  use RoadTypeModule, only : RoadType, dp

  implicit none

  !! Global array of all roads. To be used as a back door into the data
  !! structure for the ROAD subroutine, which has to have a preset argument
  !! list and be F77-callable.
  type(RoadType), save, pointer :: rtm_gRoads(:)

  real(dp), save :: scaleToSI ! Scaling factor to the SI lengt unit (meter)


contains

  subroutine functionRoad (road, dis, ropar, &
       &                   nroad, z, dz, ddz, dflag, prsurf, ierr)

    !!==========================================================================
    !! Defines a road stretching out along a line parallel to the XY-plane.
    !! A function defines the profile of the road. The road elevation and
    !! associated derivatives are computed as follows:
    !!
    !!    Z = z_0 + F(x_loc),  x_loc = x_0 + X*cos(theta) + Y*sin(theta)
    !!
    !!    dZ/dX = F'(x_loc) * cos(theta)
    !!    dZ/dY = F'(x_loc) * sin(theta)
    !!
    !!    d2Z/dX2  = F'(x_loc) * cos(theta)*cos(theta)
    !!    d2Z/dY2  = F'(x_loc) * sin(theta)*sin(theta)
    !!    d2Z/dXdY = F'(x_loc) * cos(theta)*sin(theta)
    !!
    !! where X and Y are the global coordinates of current road contact point,
    !! theta is the angle between the global X-axis and the road direction, and
    !! x_0 and z_0 are offsets along the local X-axis and the vertical axis,
    !! respectively. The array ropar is assumed to contain the following data:
    !!
    !!    ropar(1) = current tire radius
    !!    ropar(2) = cos(theta)
    !!    ropar(3) = sin(theta)
    !!    ropar(4) = x_0
    !!    ropar(5) = z_0
    !!
    !! Programmer : Trond Arne Svidal                    date/rev : Jul 2002/1.0
    !!              Trond Arne Svidal                               Oct 2002/2.0
    !!              Knut Morten Okstad                              Mar 2003/3.0
    !!==========================================================================

    use IdTypeModule        , only : getId
    use FunctionTypeModule  , only : FunctionValue, FunctionDerivative
    use EngineRoutinesModule, only : EngineValue
    use reportErrorModule   , only : reportError, error_p, warning_p

    type(RoadType), intent(inout) :: road
    real(dp)      , intent(in)    :: dis(2), ropar(:)
    integer       , intent(out)   :: nroad, dflag, ierr
    real(dp)      , intent(out)   :: z, dz(2), ddz(3), prsurf(:)

    !! Local variables
    real(dp), parameter :: eps2der_p = 1.0e-3_dp
    real(dp), parameter :: nroadThreshold_p = 1.0e-6_dp
    real(dp)            :: localRoadRad, X, Y, dZdX, d2ZdX2

    !! --- Logic section ---

    ierr  = 0 ! 0 : No error.
    !           1 : Warning.
    !           2 : Error, do not use result.
    !           3 : Fatal error, calling program should stop.

    nroad = 0 ! 0 : Not possible to tell whether the road is even or not.
    !           1 : Even road.
    !           2 : Uneven road.
    !           3 : Rough road.

    dflag = 0 ! 0 : No road derivatives are available.
    !           1 : Only first derivatives of the road are available.
    !           2 : Both 1st and 2nd derivatives of the road are available.

    z   = 0.0_dp
    dz  = 0.0_dp
    ddz = 0.0_dp

    if (size(ropar) < 5) then
       ierr = 3 ! Fatal error
       call reportError (error_p,'Road'//getId(road%id),'Invalid road data', &
            &            addString='functionRoad')
       return
    end if

    !! --- Find new position ---

    if (associated(road%xTrEng)) then
       X = dis(1) - EngineValue(road%xTrEng,ierr)
       if (ierr /= 0) then
          ierr = 3 ! Fatal error
          call reportError (error_p,'Road'//getId(road%id), &
               &            'Evaluating X-translation engine failed', &
               &            addString='functionRoad')
          return
       end if
    else
       X = dis(1)
    end if

    if (associated(road%yTrEng)) then
       Y = dis(2) - EngineValue(road%yTrEng,ierr)
       if (ierr /= 0) then
          ierr = 3 ! Fatal error
          call reportError (error_p,'Road'//getId(road%id), &
               &            'Evaluating Y-translation engine failed', &
               &            addString='functionRoad')
          return
       end if
    else
       Y = dis(2)
    end if


    !! --- Calculate the road elevation z ---

    X = ropar(4) + X*ropar(2) + Y*ropar(3)
    z = ropar(5) + FunctionValue(road%roadFunc,X,ierr)
    if (ierr /= 0) then
       ierr = 3 ! Fatal error
       call reportError (error_p,'Road'//getId(road%id), &
            &            'Evaluating road function failed', &
            &            addString='functionRoad')
       return
    end if

    !! --- Find new position ---

    if (associated(road%zTrEng)) then
       z = z + EngineValue(road%zTrEng,ierr)
       if (ierr /= 0) then
          ierr = 3 ! Fatal error
          call reportError (error_p,'Road'//getId(road%id), &
               &            'Evaluating Z-translation engine failed', &
               &            addString='functionRoad')
          return
       end if
    end if

    !! --- Calculate dz ---

    dZdX = FunctionDerivative(road%roadFunc,X,1,ierr)
    if (ierr /= 0) then
       ierr = 1 ! Warning
       call reportError (warning_p,'Road'//getId(road%id), &
            &            'Evaluating road dz failed')
    else
       dflag = 1
       dz(1) = dZdX*ropar(2)
       dz(2) = dZdX*ropar(3)
    end if

    !! --- Calculate ddz ---

    if (dflag == 1) then
       d2ZdX2 = FunctionDerivative(road%roadFunc,X,2,ierr,eps2der_p)
       if (ierr /= 0) then
          ierr = 1 ! Warning, only 1st derivative available
          call reportError (warning_p,'Road'//getId(road%id), &
               &            'Evaluating road ddz failed')
       else
          dflag  = 2
          ddz(1) = d2ZdX2*ropar(2)*ropar(2)
          ddz(2) = d2ZdX2*ropar(3)*ropar(3)
          ddz(3) = d2ZdX2*ropar(2)*ropar(3)
       end if
    else
       d2ZdX2 = 0.0_dp
    end if

    !! --- Calculate nroad ---

    if (dflag == 2) then ! both dz and ddz are available
       localRoadRad = max(abs(dZdX),abs(d2ZdX2))
       if (localRoadRad < nroadThreshold_p) then
          nroad = 1 ! Even road
       else
          !! The road is not flat and must be compared with the tire radius.
          !! Strictly speaking, nroad will depend on the direction between the
          !! tire and road. Consider the worst case first, along the function.
          localRoadRad = 1.0_dp/kappa(dZdX,d2ZdX2)
          if (localRoadRad > ropar(1)) then
             nroad = 2 ! Uneven road
          else
             nroad = 3 ! Rough road
          end if
       end if
    else
       ierr = 1 ! dz and ddz not available, cannot to say anything about nroad
       call reportError (warning_p,'Road'//getId(road%id), &
            &            'Evaluating road nroad failed')
    end if

    !! --- Return psurf ---
    prsurf = 0.0_dp
    !! TODO TAS: It is not known what psurf should contain at this moment 8/8-02

  contains

    !! Curvature for a 2D curve on the form y = f(x):
    !! kappa = (d2y/dx2)/(1+(dy/dx)^2)^(3/2)
    function kappa (firstDerVal,secDerVal)
      real(dp),intent(in) :: firstDerVal, secDerVal
      real(dp)            :: kappa
      kappa = secDerVal / (1.0_dp + firstDerVal*firstDerVal)**1.5_dp
    end function kappa

  end subroutine functionRoad


  subroutine surfaceRoad (dis, ropar, nroad, z, dz, ddz, dflag, prsurf, ierr)

    !!==========================================================================
    !! Defines a 3D road surface with the surface in XY plane. The surface is
    !! defined on a ftl-file (or equivalent). The tire inputs X and Y, and the
    !! routine returns Z, its 1st and 2nd derivatives in earth fixed CS.
    !!
    !! Programmer : Trond Arne Svidal                    date/rev : Jul 2002/1.0
    !!==========================================================================

    real(dp), intent(in)  :: dis(2), ropar(:)
    integer , intent(out) :: nroad, dflag, ierr
    real(dp), intent(out) :: z, dz(2), ddz(3), prsurf(:)

    !! --- Logic section ---

    ierr   = 0
    nroad  = 0
    dflag  = 0
    z      = 0.0_dp
    dz     = 0.0_dp
    ddz    = 0.0_dp
    prsurf = 0.0_dp
    print *,' ** surfaceRoad dummy: ',dis,size(ropar)

  end subroutine surfaceRoad

end module RoadRoutinesModule


!! Callback routine for road evaluation. Since it is used as argument to the
!! DTYRE call (which is in F77), it has to to be defined outside any modules.

subroutine FedemRoadMain (time, dis, iflag, jflag, idtyre, idroad, &
       &                  nropar, ropar, nchrds, chrdst, nprsur, &
       &                  nroad, z, dz, ddz, dflag, prsurf, ierr)

  !=============================================================================
  !
  !     This subroutine handles the road. It is used by the
  !     STI/TYDEX tire model, release 1.2.
  !
  !     Input arguments:
  !     ----------------
  !     time   - Current simulation time [s] (For error messages, etc.)
  !     dis    - x and y coordinates in earth fixed CS [m]
  !     iflag  - Specifies if derivative of the road surface (dz, ddz) are
  !              are needed by the tire model
  !     jflag  - Specifies initialization phase for each wheel (each tire is
  !              handled independently, there are no links between the different
  !              tires on a vehicle)
  !              =  0 : Normal mode.
  !              =  1 : First initialization before starting the simulation
  !                     (once for each wheel)
  !              = 99 : Final call of the road model
  !     idtyre - Tire ID ("Serial number" of tire at the vehicle)
  !     idroad - Specifies which road model method is used
  !     nropar - Specifies the dimension of array ropar [diff]
  !     ropar  - Contains parameters for the road model
  !     nchrds - Specifies number of characters in chrdst
  !     chrdst - File name of the road file
  !     nprsur - Dimension of psurf
  !
  !     Output arguments:
  !     -----------------
  !     nroad  - Short description of road, can be used to change complexity of
  !              tire model or save computational time
  !              = 0 : No road data available
  !              = 1 : Even road (z = constant, dz = ddz = 0)
  !              = 2 : Uneven road (wavelength of the road > diameter of tire)
  !              = 3 : Rough road  (wavelength of the road < diameter of tire)
  !     z      - z component of the road in earth fixed CS [m]
  !     dz     - 1st derivative of the z component of the road (if dflag >= 1)
  !              in the earth fixed CS
  !              dz(1) = dz/dx
  !              dz(2) = dz/dy
  !     ddz    - 2nd derivative of the z component of the road (if dflag = 2)
  !              in the earth fixed CS [1/m]
  !              ddz(1) = d2z/dx2
  !              ddz(2) = d2z/dy2
  !              ddz(3) = d2z/dxdy
  !     dflag  - Specifies number of derivatives available
  !              = 0 : No derivative available
  !              = 1 : 1st derivative available (dz)
  !              = 2 : 1st and 2nd derivative available (dz, ddz)
  !     prsurf - Properties of the road, characteristics at x, y position [diff]
  !     ierr   - Error flag
  !              = 0 : No error
  !              = 1 : Warning
  !              = 2 : Error, do not use result
  !              = 3 : Fatal error, calling program should stop
  !
  ! Comments : Units are specified by the STI interface.
  !            diff = units specified by the tire model.
  !
  ! Programmer : Trond Arne Svidal                       date/rev : Jul 2002/1.0
  !=============================================================================

  use RoadTypeModule     , only : dp, FUNCTION_ROAD_P, SURFACE_ROAD_P
  use RoadRoutinesModule , only : functionRoad, surfaceRoad
  use RoadRoutinesModule , only : rtm_gRoads, scaleToSI
  use reportErrorModule  , only : internalError, reportError, debugFileOnly_p
#ifdef FT_DEBUG
  use fileUtilitiesModule, only : getDBGfile
  use dbgUnitsModule     , only : dbgRoad
#endif

  implicit none

  integer       , intent(in)  :: nropar, nchrds, nprsur
  real(dp)      , intent(in)  :: time, dis(2), ropar(nropar)
  integer       , intent(in)  :: iflag, jflag, idtyre, idroad
  character*(*) , intent(in)  :: chrdst
  integer       , intent(out) :: nroad, dflag, ierr
  real(dp)      , intent(out) :: z, dz(2), ddz(3), prsurf(nprsur)

  !! --- Logic section ---

  if (idroad < 0 .or. idroad > size(rtm_gRoads)) then
     print *,'*** FedemRoadMain: ',time,iflag,jflag,idtyre,chrdst(1:nchrds)
     ierr = 4 + internalError('FedemRoadMain: Invalid road index')
     return
#ifdef FT_DEBUG
  else if (dbgRoad == 0) then
     dbgRoad = getDBGfile(15,"road.dbg")
#endif
  end if

  select case (rtm_gRoads(idroad)%type)

  case (FUNCTION_ROAD_P)

#ifdef FT_DEBUG
     write(dbgRoad,"('FedemRoadMain: functionRoad')")
#endif
     call functionRoad (rtm_gRoads(idroad), dis/scaleToSI, ropar, &
          &             nroad, z, dz, ddz, dflag, prsurf, ierr)

  case (SURFACE_ROAD_P)

#ifdef FT_DEBUG
     write(dbgRoad,"('FedemRoadMain: surfaceRoad')")
#endif
     call surfaceRoad (dis/scaleToSI, ropar, &
          &            nroad, z, dz, ddz, dflag, prsurf, ierr)

  case default

     ierr = 4 + internalError('FedemRoadMain: Invalid road type')
     return

  end select

  if (ierr == 3) then
     call reportError (debugFileOnly_p,'FedemRoadMain')
  else if (scaleToSI > 0.0_dp) then
     !! Scale to SI units, [m] and [1/m] respectively
     z   = z   * scaleToSI
     ddz = ddz / scaleToSI
  end if

#ifdef FT_DEBUG
  write(dbgRoad,"(1p,'idTire=',i4,' time=',e12.3,' dis=',2e12.3,' z=',e12.3)") &
       &        idtyre, time, dis, z
#endif

end subroutine FedemRoadMain
