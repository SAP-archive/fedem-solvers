!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module JD_namelistModule

  use InputUtilities, only : iuSetPosAtFirstEntry

  implicit none

contains

  logical function isNameListFile( file )
    integer, intent(in) :: file

    isNameListFile = iuSetPosAtFirstEntry(file,'&STIFF_VERT')

  end function isNameListFile

  !! === Vertical damping read and write

  subroutine read_vertDamping_nml( file, cdampz, cdefl )
    integer ,intent(in) :: file
    real*8  ,intent(out):: cdampz, cdefl

    namelist /DAMP_VERT/ cdampz, cdefl

    cdampz = 0.0d0
    cdefl  = 0.0d0
    if (iuSetPosAtFirstEntry(file,'&DAMP_VERT')) then
       read(file,nml=DAMP_VERT)
    end if

  end subroutine read_vertDamping_nml

  subroutine write_vertDamping_nml( file, cdampz, cdefl )
    integer, intent(in) :: file
    real*8,  intent(in) :: cdampz, cdefl

    namelist /DAMP_VERT/ cdampz, cdefl

    write(file,nml=DAMP_VERT)
    write(file,*)''

  end subroutine write_vertDamping_nml


  !! === Vertical stiffness read and write

  subroutine read_vertStiff_nml(file, ivert, akv, defl, force, actReg, nLdPt, anSeg_dp )

    use jdt_max, only : maxCPT

    integer, intent(in) :: file

    integer ivert,nLdPt,anSeg
    real*8 akv, defl(maxCPT), force(maxCPT), actReg, anseg_dp
    namelist/STIFF_VERT/ ivert,akv, nLdPt, defl, force, actReg, anSeg

    ivert = 1
    akv   = 0.0d0
    defl  = 0.0d0
    force = 0.0d0
    actReg = 0.0d0
    nLdPt  = 0
    anSeg = 0
    if (iuSetPosAtFirstEntry(file,'&STIFF_VERT')) then
       read(file,nml=STIFF_VERT)
    end if
    anseg_dp = dble(anseg)

  end subroutine read_vertStiff_nml

  subroutine write_vertStiff_nml(file,ivert,akv,defl,force,actReg,nLdPt,AnSeg_dp )

    use jdt_max, only : maxCPT

    integer, intent(in) :: file

    integer ivert,nLdPt,anseg
    real*8 akv, defl(maxCPT), force(maxCPT), actReg, anSeg_dp
    namelist/STIFF_VERT/ ivert, akv, nLdPt, defl, force, actReg, anSeg

    anseg = int(anseg_dp)
    write(file,nml=STIFF_VERT)
    write(file,*)''

  end subroutine write_vertStiff_nml


  !! === Tire Radius data

  subroutine read_tire_radius_data_nml(file, wrad, yarm, rr )

    integer, intent(in) :: file
    real*8 , intent(out):: wrad, yarm, rr
    namelist/TIRE_RADIUS_DATA/  wrad, yarm, rr

    wrad  = 0.0d0
    yarm  = 0.0d0
    rr    = 0.0d0

    if (iuSetPosAtFirstEntry(file,'&TIRE_RADIUS_DATA')) then
       read(file,nml=TIRE_RADIUS_DATA)
    end if

  end subroutine read_tire_radius_data_nml

  subroutine write_tire_radius_data_nml(file, wrad, yarm, rr )

    integer, intent(in) :: file
    real*8 , intent(in) :: wrad, yarm, rr
    namelist/TIRE_RADIUS_DATA/  wrad, yarm, rr

    write(file,nml=TIRE_RADIUS_DATA)
    write(file,*)''

  end subroutine write_tire_radius_data_nml


  !! data for kinematic differential model (should not be used in Fedem, really

  subroutine read_differential_data_nml(file, angvel, alat, trdw )

    integer, intent(in) :: file
    real*8 , intent(out):: angvel, alat, trdw
    namelist/DIFFERENTIAL_DATA/ angvel, alat, trdw

    angvel  = 0.0d0
    alat    = 0.0d0
    trdw    = 0.0d0

    if (iuSetPosAtFirstEntry(file,'&DIFFERENTIAL_DATA')) then
       read(file,nml=DIFFERENTIAL_DATA)
    end if

  end subroutine read_differential_data_nml

  subroutine write_differential_data_nml(file, angvel, alat, trdw )

    integer, intent(in) :: file
    real*8 , intent(in) :: angvel, alat, trdw
    namelist/DIFFERENTIAL_DATA/ angvel, alat, trdw

    write(file,nml=DIFFERENTIAL_DATA)
    write(file,*)''

  end subroutine write_differential_data_nml


  !! Road surface data
  subroutine read_surface_data_nml(file,IROTAT,IPOWER,ISURF,ICONST,CI,B,D,H,DELTA)

    integer, intent(in) :: file
    integer, intent(out):: irotat,ipower,isurf,iconst
    real*8,  intent(out):: CI,B,D,H,DELTA
    namelist/SURFACE_DATA/ IROTAT,IPOWER,ISURF,ICONST,CI,B,D,H,DELTA

    irotat = 0 !! Rotation not considered
    ipower = 1 !! consider rolling resistance, but not traction
    isurf  = 2 !! hard surface
    iconst = 1 !! radial tire
    CI = 0.0d0
    B = 0.0d0
    D = 0.0d0
    H = 0.0d0
    DELTA = 0.0d0

    if (iuSetPosAtFirstEntry(file,'&SURFACE_DATA')) then
       read(file,nml=SURFACE_DATA)
    end if

  end subroutine read_surface_data_nml

  subroutine write_surface_data_nml(file,IROTAT,IPOWER,ISURF,ICONST,CI,B,D,H,DELTA)

    integer, intent(in) :: file
    integer, intent(in) :: irotat,ipower,isurf,iconst
    real*8,  intent(in) :: CI,B,D,H,DELTA
    namelist/SURFACE_DATA/ IROTAT,IPOWER,ISURF,ICONST,CI,B,D,H,DELTA

    write(file,nml=SURFACE_DATA)
    write(file,*)''

  end subroutine write_surface_data_nml


  !! Lateral stiffness data

  subroutine read_StiffDampLat_nml(file, ilat, akl, clat, yse_from_fedem, yse, &
       &                       iyfric, amul, alexp, &
       &                       irollm, aklat)

    integer, intent(in) :: file
    logical, intent(out):: yse_from_fedem
    integer, intent(out):: ilat, iyfric, irollm
    real*8,  intent(out):: akl, clat, yse, amul, alexp, aklat
    namelist/STIFFDAMP_LAT/ ilat, akl, clat, yse_from_fedem, yse, iyfric, amul, alexp, &
         &                  irollm, aklat

    ilat = 0
    irollm = 0
    akl = 0.0d0
    clat = 0.0d0
    yse_from_fedem = .false.
    yse = 0.0d0
    iyfric = 0
    amul = 0.0d0
    alexp = 0.0d0
    irollm = 0
    aklat = 0.0d0

    if (iuSetPosAtFirstEntry(file,'&STIFFDAMP_LAT')) then
       read(file,nml=STIFFDAMP_LAT)
    end if

  end subroutine read_StiffDampLat_nml

  subroutine write_StiffDampLat_nml(file, ilat, akl, clat, yse_from_fedem, yse, &
       &                       iyfric, amul, alexp, &
       &                       irollm, aklat)

    integer, intent(in) :: file
    integer, intent(in) :: ilat, iyfric, irollm
    logical, intent(in) :: yse_from_fedem
    real*8,  intent(in) :: akl, clat, yse, amul, alexp, aklat
    namelist/STIFFDAMP_LAT/ ilat, akl, clat, yse_from_fedem, yse, iyfric, amul, alexp, &
         &                  irollm, aklat

    logical :: niceprint = .true.
    character :: c

    if ( .not. niceprint ) then
       write(file,nml=STIFFDAMP_LAT)
    else
       write(file,"('&STIFFDAMP_LAT')")
       write(file,"('   ILAT   = ',i2,',')") ilat
       if (ilat == 1) then
          c = ' '; else;  c = '!'
       endif
       write(file,"(A,'   AKL    = ',1p,e12.3,',')") c, akl
       write(file,"(A,'   CLAT   = ',1p,e12.3,',')") c,clat
       if (yse_from_fedem) then
          write(file,"(A,'   YSE_FROM_FEDEM = T,')") c
       else
          write(file,"(A,'   YSE_FROM_FEDEM = F,')") c
       end if
       if (ilat == 1) then; c = ' '; else;  c = '!'
       endif
       write(file,"(A,'   YSE    = ',1p,e12.3,',')") c, yse
       write(file,"(A,'   IYFRIC = ',i2,',')")       c, iyfric
       if ((ilat == 1 .and. iyfric == 1) .or. ilat == 2) then
          c = ' '; else;  c = '!'
       endif
       write(file,"(A,'   AMUL   = ',1p,e12.3,',')") c, amul
       if (ilat == 2) then
          c = ' '; else;  c = '!'
       endif
       write(file,"(A,'   ALEXP  = ',1p,e12.3,',')") c, alexp
       c = ' '
       write(file,"(A,'   IROLLM = ',i2,',')")       c, irollm
       if (irollm == 1) then
          c = ' '; else;  c = '!'
       endif
       write(file,"(A,'   AKLAT  = ',1p,e12.3,',')") c, aklat
       write(file,"('/')")
    endif
    write(file,*)''

  end subroutine write_StiffDampLat_nml


  !! Longitudinal stiffness and damping

  subroutine read_StiffDampLong_nml(file, &
       &     ilong,akfa,cfa,xse,xarm, &
       &     ixfric,amufa)

    integer, intent(in) :: file
    integer, intent(out):: ilong, ixfric
    real*8,  intent(out):: akfa,cfa,xse,xarm,amufa
    namelist/STIFFDAMP_LONG/ ilong,akfa,cfa,xse,xarm, ixfric,amufa

    ilong = 0
    akfa = 0.0d0
    cfa  = 0.0d0
    xse  = 0.0d0
    xarm = 0.0d0
    ixfric = 0
    amufa = 0.0d0

    if (iuSetPosAtFirstEntry(file,'&STIFFDAMP_LONG')) then
       read(file,nml=STIFFDAMP_LONG)
    end if

  end subroutine read_StiffDampLong_nml

  subroutine write_StiffDampLong_nml(file, &
       &     ilong,akfa,cfa,xse,xarm, &
       &     ixfric,amufa)

    integer, intent(in) :: file
    integer, intent(in) :: ilong, ixfric
    real*8,  intent(in) :: akfa,cfa,xse,xarm,amufa
    namelist/STIFFDAMP_LONG/ ilong,akfa,cfa,xse,xarm, ixfric,amufa

    write(file,nml=STIFFDAMP_LONG)
    write(file,*)''

  end subroutine write_StiffDampLong_nml

end module JD_namelistModule
