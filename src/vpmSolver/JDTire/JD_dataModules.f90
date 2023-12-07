!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module JDT_FDM

  use KindModule, only : dp

  INTEGER              :: outFile, echoFile
  INTEGER, allocatable :: TireFiles(:)
  LOGICAL, allocatable :: hasReadTire(:), YSE_from_fedem(:)
  !! Either the road is a JD road file or it is a Fedem defined road function
  integer, allocatable :: RoadFiles(:), FedemRoadIdin(:)

  REAL*8,  allocatable :: yFricCoeff(:), yFricExp(:)
  REAL*8,  allocatable :: yFricSlipAngle(:), rollRadius(:)

  real(dp), allocatable :: forcePrevL(:,:), torquePrevL(:,:)
  real(dp), allocatable :: forcePrevLOut(:,:), torquePrevLout(:,:)
  real(dp), allocatable :: posPrevG(:,:), velPrevG(:,:)

contains

  subroutine jdt_allocate_fdm (nJDT,stat)
    integer, intent(in)  :: nJDT
    integer, intent(out) :: stat

    stat = 0
    if (stat==0) allocate( TireFiles(nJDT), RoadFiles(nJDT), FedemRoadIdin(nJDT),stat=stat)
    if (stat==0) allocate( hasReadTire(nJDT), YSE_from_fedem(nJDT), stat=stat)
    if (stat==0) allocate( yFricCoeff(nJDT), yFricExp(nJDT), stat=stat)
    if (stat==0) allocate( yFricSlipAngle(nJDT), rollRadius(nJDT), stat=stat)

    if (stat==0) allocate( forcePrevL(3,nJDT), torquePrevL(3,nJDT), stat=stat)
    if (stat==0) allocate( forcePrevLOut(3,nJDT), torquePrevLout(3,nJDT), stat=stat)
    if (stat==0) allocate( posPrevG(3,nJDT), velPrevG(3,nJDT), stat=stat)

  end subroutine jdt_allocate_fdm

end module JDT_FDM


module JDT_INP

  INTEGER              ::  NTIRES,IRAMKR,IGMKR,IENG,NTSPTS, NSHOCK
  INTEGER, allocatable :: IMRKR(:)
  INTEGER, allocatable :: IROTAT(:),ISTEER(:),IPOWER(:)
  INTEGER, allocatable :: IVERT(:),ILAT(:),IROLLM(:)
  INTEGER, allocatable :: ILONG(:),ISURF(:),ISPLIN(:)
  INTEGER, allocatable :: IWMKR(:),IVMKR(:),IWCMKR(:)
  INTEGER, allocatable :: IYFRIC(:),IXFRIC(:),ICONST(:)
  INTEGER, allocatable :: IWROT(:),JWROT(:)
  INTEGER, allocatable :: NEPTS(:),NCPTS(:)
  INTEGER, allocatable :: ISHCK(:),JSHCK(:)

  REAL*8              :: ETA
  REAL*8, allocatable :: B(:),D(:),RR(:)
  REAL*8, allocatable :: ANGVEL(:),ALAT(:),TRDW(:)
  REAL*8, allocatable :: CI(:),H(:),DELTA(:)
  REAL*8, allocatable :: AKV(:),CDAMPZ(:),CDEFL(:),AKL(:)
  REAL*8, allocatable :: CLAT(:),AMUL(:),YSE(:),ALEXP(:)
  REAL*8, allocatable :: YARM(:),AKLAT(:),AKFA(:)
  REAL*8, allocatable :: CFA(:),AMUFA(:),XSE(:),XARM(:)
  REAL*8, allocatable :: ENGSPD(:),ENGTOR(:),VELEX(:,:)
  REAL*8, allocatable :: FOREX(:,:),VELCP(:,:),FORCP(:,:)

contains

  subroutine jdt_allocate_inp (nJDT,stat)
    integer, intent(in)  :: nJDT
    integer, intent(out) :: stat

    stat = 0
    if (stat == 0) allocate( IMRKR(nJDT), stat=stat)
    if (stat == 0) allocate( IROTAT(nJDT),ISTEER(nJDT),IPOWER(nJDT), stat=stat)
    if (stat == 0) allocate( IVERT(nJDT),ILAT(nJDT),IROLLM(nJDT), stat=stat)
    if (stat == 0) allocate( ILONG(nJDT),ISURF(nJDT),ISPLIN(nJDT), stat=stat)
    if (stat == 0) allocate( IWMKR(nJDT),IVMKR(nJDT),IWCMKR(nJDT), stat=stat)
    if (stat == 0) allocate( IYFRIC(nJDT),IXFRIC(nJDT),ICONST(nJDT), stat=stat)
    if (stat == 0) allocate( IWROT(nJDT),JWROT(nJDT), stat=stat)
    if (stat == 0) allocate( NEPTS(nJDT),NCPTS(nJDT), stat=stat)
    if (stat == 0) allocate( ISHCK(nJDT),JSHCK(nJDT), stat=stat)

    if (stat == 0) allocate( B(nJDT),D(nJDT),RR(nJDT), stat=stat)
    if (stat == 0) allocate( ANGVEL(nJDT),ALAT(nJDT),TRDW(nJDT), stat=stat)
    if (stat == 0) allocate( CI(nJDT),H(nJDT),DELTA(nJDT), stat=stat)
    if (stat == 0) allocate( AKV(nJDT),CDAMPZ(nJDT),CDEFL(nJDT),AKL(nJDT), stat=stat)
    if (stat == 0) allocate( CLAT(nJDT),AMUL(nJDT),YSE(nJDT),ALEXP(nJDT), stat=stat)
    if (stat == 0) allocate( YARM(nJDT),AKLAT(nJDT),AKFA(nJDT), stat=stat)
    if (stat == 0) allocate( CFA(nJDT),AMUFA(nJDT),XSE(nJDT),XARM(nJDT), stat=stat)
    if (stat == 0) allocate( ENGSPD(nJDT),ENGTOR(nJDT),VELEX(10,nJDT), stat=stat)
    if (stat == 0) allocate( FOREX(10,nJDT),VELCP(10,nJDT),FORCP(10,nJDT), stat=stat)

  end subroutine jdt_allocate_inp

end module JDT_INP


module JDT_MAX

  integer, parameter :: maxRDP = 70000 ! max number of road points in the model
  integer, parameter :: maxCPT =    10 ! max number of force-deflection points

end module JDT_MAX


module JDT_VAR

  REAL*8, allocatable :: WRAD(:),ANSEG(:),ACTREG(:),ASEGK(:,:)
  REAL*8, allocatable :: X1(:),Y1(:),SM(:,:),TMAG(:),tmag2(:),springStiff(:)
  REAL*8, allocatable :: SSEGA(:,:),CSEGA(:,:),STA(:,:)
  REAL*8, allocatable :: HGT(:,:),FX(:),FZ(:)
  REAL*8, allocatable :: DEFL(:,:),AREA(:,:)
  REAL*8, allocatable :: FORCE(:,:),DTHET(:)
  REAL*8, allocatable :: stiffM(:,:,:)

  INTEGER              :: LTIRE
  INTEGER, allocatable :: NLDPT(:),NTSTA(:),ITR(:)

contains

  subroutine jdt_allocate_var(nJDT, stat)

    use JDT_MAX

    integer, intent(in)  :: nJDT
    integer, intent(out) :: stat

    stat = 0
    if (stat==0) allocate( WRAD(nJDT),ANSEG(nJDT),ACTREG(nJDT),ASEGK(36,nJDT), stat=stat)
    if (stat==0) allocate( X1(nJDT),Y1(nJDT),SM(36,nJDT),TMAG(36),tmag2(36),springStiff(36),stat=stat)
    if (stat==0) allocate( SSEGA(36,nJDT),CSEGA(36,nJDT),STA(maxRDP,nJDT), stat=stat)
    if (stat==0) allocate( HGT(maxRDP,nJDT),FX(nJDT),FZ(nJDT), stat=stat)
    if (stat==0) allocate( DEFL(maxCPT,nJDT),AREA(maxCPT,nJDT), stat=stat)
    if (stat==0) allocate( FORCE(maxCPT,nJDT),DTHET(nJDT), stat=stat)
    if (stat==0) allocate( NLDPT(nJDT),NTSTA(nJDT),ITR(nJDT), stat=stat)
    if (stat==0) allocate( stiffM(6,6,nJDT), stat=stat)

    stiffM = 0.0d0

  end subroutine jdt_allocate_var

end module JDT_VAR
