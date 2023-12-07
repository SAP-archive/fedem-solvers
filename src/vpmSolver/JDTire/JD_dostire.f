C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE READINPUT_FDM_NML(ID,TIME,PAR,NPAR,ierr)

      use JD_namelistModule
      use jdt_fdm
      use jdt_inp
      use jdt_var
      use fileUtilitiesModule

      IMPLICIT NONE

      integer id, npar, ierr
      REAL*8 TIME, PAR(*)
C
      integer IFPASS
      real*8 RTOD
      real*8 PI
      real*8 POT
      real*8 DTR
      DATA IFPASS /0/
      DATA RTOD /57.29578D0/
      DATA PI/3.1415926535897932D0/
      DATA POT/1.5707963267948D0/
      DATA DTR/0.01745329/

      integer i, j, k, npts, nseg, nsegd2, nsegp1, imex, nmlFile
      real*8 rseg1, dseg1, rseg, segang, wchgt, r2dfl, adfl, bdfl
      real*8 xzero, yzero, ttad, ybar, al1, al2, a1, a2, ttar, x2, y2
      real*8 x3, y3, xk, yk, absxk, arg, x24, y24, deltax
C
C     === Logic section
      ierr = 0
C
C     READ THE TIRE DATA FILE ON THE FIRST PASS THROUGH THIS SUBROUTINE
C
      WRITE(outFile,2)
 2    FORMAT(//,3X,' ======== TIRE DATA FILE ========',/)
c
C     IENG = 0 - DON'T CONSIDER ENGINE TORQUE-SPEED RELATION
C     1 - DO CONSIDER ENGINE TORQUE-SPEED RELATION
C
      IENG = 0
      WRITE(outFile,"('IENG',I3)") IENG
C
C     GROUND ORIGIN MARKER
      IGMKR  = 0
C     IRAMKR = MARKER AT REAR AXLE CENTER
      IRAMKR = 0
C
C     READ THE NUMBER OF TIRES
C
      NTIRES = NTIRES + 1
      WRITE (outFile,"('Number of tires so far ',i3)") NTIRES
C
C     READ APPROPRIATE DATA FOR EACH TIRE
C
C     Dissever tire number from input identifier
      J = ID - (ID/100)*100
      rewind(tireFiles(j))
C
C     IMRKR(J) = MARKER AT WHICH FORCE ACTUALLY ACTS
C
      IMRKR(J) = J
      WRITE(outFile,20)'IMRKR',IMRKR(J)
 20   FORMAT(A,X,10I5)
C
C     IWMKR(J) = MARKER AT CENTER OF WHEEL
C     IVMKR(J) = REFERENCE MARKER ON PART TO WHICH WHEEL IS ATTACHED
C
      IWMKR(J) = J
      IVMKR(J) = J
      WRITE(outFile,20)'IWMKR IVMKR',IWMKR(J),IVMKR(J)
C
C     WRAD(J) = UNDEFLECTED TIRE RADIUS
C     YARM(J) = MOMENT ARM OF LATERAL FORCE - POSITIVE IF FORCE ACTS
C               BELOW WHEEL CENTER
      call read_tire_radius_data_nml(TireFiles(j),WRAD(J),YARM(J),RR(J))
      call write_tire_radius_data_nml(outFile,WRAD(J),YARM(J),RR(J))
      rollRadius(J) = wrad(j)
      rollRadius(J) = rr(j)
C
C     CDAMPZ(J) = VERTICAL DAMPING RATE OF TIRE
C     CDEFL(J) = TIRE DEFLECTION AT WHICH DAMPING REACHES MAXIMUM VALUE
C
      call read_vertDamping_nml(TireFiles(j),cdampz(j),cdefl(j))
      call write_vertDamping_nml(outFile,cdampz(j),cdefl(j))
C
C     IROTAT(J) = 0 - ROTATION NOT CONSIDERED
C     1 - ROTATION CONSIDERED USING KINEMATIC MODEL OF DIFFERENTIAL
C     2 - ROTATION DETERMINED FROM MODEL
C     ISTEER(J) = 0 - NO STEERING INPUT (RIDE ANALYSIS)
C     1 - USER INPUT - WHEEL DOESN'T PHYSICALLY STEER IN MODEL
C     2 - STEER ANGLE DETERMINED FROM MODEL
C     IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C     1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C     2 - POWERED WHEEL - CONSIDER TRACTION & ROLLING RESISTANCE
C     ILONG(J)  = 0 - NO FORE-AFT SPRING AND DAMPER
C     1 - USE FORE-AFT SPRING AND DAMPER

      isteer(j) = 0
      write(outFile,"('Setting ISTEER to 0')")

C     INPUT IF TRACTION AND/OR ROLLING RESISTANCE ARE CONSIDERED
C
C     ISURF(J)  = 1 - SOIL SURFACE    = 2 - HARD SURFACE
C     ICONST(J) = 1 - BIAS TIRE CONSTRUCTION  = 2 - RADIAL TIRE CONSTRUCTION
C     CI(J) = CONE INDEX (FOR SOIL) OR K (FOR HARD SURFACE)
C     B(J) = UNLOADED TIRE SECTION WIDTH
C     D(J) = UNLOADED TIRE DIAMETER
C     RR(J) = TIRE ROLLING RADIUS
C     H(J) = TIRE SECTION HEIGHT
C     DELTA(J) = LOADED TIRE DEFLECTION
C

      call read_surface_data_nml(TireFiles(j),IROTAT(J),IPOWER(J),
     *     ISURF(J),ICONST(J),CI(J),B(J),D(J),H(J),DELTA(J))
      call write_surface_data_nml(outFile,IROTAT(J),IPOWER(J),
     *     ISURF(J),ICONST(J),CI(J),B(J),D(J),H(J),DELTA(J))
C     INPUT DATA IF KINEMATIC DIFFERENTIAL MODEL IS USED
C
C     ANGVEL(J) = ANGULAR VELOCITY OF AXLE (RAD/SEC)
C     ALAT(J) = LATERAL POSITION OF WHEEL ON AXLE (IN VEHICLE SYSTEM)
C     TRDW(J) = EFFECTIVE TREAD WIDTH OF WHEELS ON AXLE
C     RR(J) = ROLLING RADIUS
C
      call read_differential_data_nml(TireFiles(j),
     *      ANGVEL(J),ALAT(J),TRDW(J))
Cbh      call write_differential_data_nml(outFile,
Cbh     *      ANGVEL(J),ALAT(J),TRDW(J))
C
C     INPUT DATA IF WHEEL ANGULAR VELOCITIES ARE FOUND FROM MODEL
C
C     IWROT(J) = I MARKER MAKING UP WHEEL-CHASSIS REVOLUTE JOINT
C     JWROT(J) = J MARKER MAKING UP WHEEL-CHASSIS REVOLUTE JOINT
C
      IF(IROTAT(J).EQ.2) then
         iwrot(j) = j
         jwrot(j) = 0
         WRITE(outFile,20)'IWROT,JWROT :', IWROT(J),JWROT(J)
      end if
C
C     IWCMKR(J) = MARKER USED TO CALCULATE STEER ANGLE
C
      IWCMKR(J) = 0
C
C     INPUT IF TRACTION AND/OR ROLLING RESISTANCE ARE CONSIDERED
C
C     ISURF(J)  = 1 - SOIL SURFACE    = 2 - HARD SURFACE
C     ICONST(J) = 1 - BIAS TIRE CONSTRUCTION  = 2 - RADIAL TIRE CONSTRUCTION
C     CI(J) = CONE INDEX (FOR SOIL) OR K (FOR HARD SURFACE)
C     B(J) = UNLOADED TIRE SECTION WIDTH
C     D(J) = UNLOADED TIRE DIAMETER
C     RR(J) = TIRE ROLLING RADIUS
C     H(J) = TIRE SECTION HEIGHT
C     DELTA(J) = LOADED TIRE DEFLECTION
C

      call read_surface_data_nml(TireFiles(j),IROTAT(J),IPOWER(J),
     *     ISURF(J),ICONST(J),CI(J),B(J),D(J),H(J),DELTA(J))
      call write_surface_data_nml(outFile,IROTAT(J),IPOWER(J),
     *     ISURF(J),ICONST(J),CI(J),B(J),D(J),H(J),DELTA(J))

C     INPUT DEPENDING ON VERTICAL TIRE FORCE MODEL

C     IVERT(J)  = 1 - LINEAR SPRING MODEL FOR VERTICAL FORCE
C     2 - FORCE CALCULATED USING LOAD-DEFLECTION RELATION
C     3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C     AKV(J) = VERTICAL SPRING RATE OF TIRE
C     NLDPT(J) = NUMBER OF LOAD-DEFLECTION POINTS
C     DEFL(I,J) = ARRAY OF TIRE DEFLECTIONS FOR TIRE J
C     FORCE(I,J) = ARRAY OF CORRESPONDING TIRE VERTICAL LOADS FOR TIRE J
C     ANSEG(J) = NUMBER OF SEGMENTS IN ACTIVE REGION
C     ACTREG(J) = SIZE OF ACTIVE REGION (IN DEGREES)
      call read_vertStiff_nml(TireFiles(j),ivert(j),akv(j),defl(1,j),
     * force(1,j),actReg(j),nLdPt(j),anseg(j))
      call write_vertStiff_nml(outFile,ivert(j),akv(j),defl(1,j),
     * force(1,j),actReg(j),nLdPt(j),anseg(j))
C
      NPTS=NLDPT(J)
C
C     CALCULATE SOME INITIAL VALUES FOR THE DISTRIBUTED CONTACT MODEL
C
      IF(IVERT(J).EQ.3) THEN
         ITR(J)=2
         NSEG=ANSEG(J)
         SEGANG=ACTREG(J)/ANSEG(J)
         DSEG1=-(ACTREG(J)-SEGANG)/2.0
         RSEG1=DTR*DSEG1
         SSEGA(1,J)=DSIN(RSEG1)
         CSEGA(1,J)=DCOS(RSEG1)
         SM(1,J)=-CSEGA(1,J)/SSEGA(1,J)
         DO 101 I=2,NSEG
            RSEG=DTR*(DSEG1+(I-1)*SEGANG)
            SSEGA(I,J)=DSIN(RSEG)
            CSEGA(I,J)=DCOS(RSEG)
            SM(I,J)=-CSEGA(I,J)/SSEGA(I,J)
 101     CONTINUE
C     TRANSLATE THE TIRE MODEL INDICATOR
C     NOTE: LTIRE=0 - INDIVIDUAL SPRING MODEL
C     LTIRE=1 - DEFLECTED AREA MODEL
C     NOTE THE TIREF SUBROUTINE IMPLEMENTS BOTH MODELS BUT HERE THE
C     DEFLECTED AREA MODEL IS THE ONLY ONE USED
         LTIRE=1
         IF(LTIRE.EQ.1) GO TO 6001
C     CALCULATIONS RELATED TO THE ORIGINAL TIRE MODEL
         NSEGD2=NSEG/2
         NSEGP1=NSEGD2+1
C     READ IN THE SEGMENT SPRING RATES FOR THE LEFT HALF OF THE ACTIVE
C     REGION
         READ(TireFiles(j),*) (ASEGK(I,J),I=1,NSEGD2)
C     SET THE SEGMENT SPRING RATES FOR THE RIGHT HALF OF THE REGION
         DO 6002 I=NSEGP1,NSEG
            ASEGK(I,J)=ASEGK(NSEG-I+1,J)
 6002    CONTINUE
         GO TO 6003
C     CALCULATIONS RELATED TO THE NEW TIRE MODEL
 6001    CONTINUE
C     DETERMINE THE DEFLECTED AREA CORRESPONDING TO EACH DEFLECTION
         DTHET(J)=SEGANG*DTR
         AREA(1,J)=0.0
         DO 6004 I=2,NPTS
            WCHGT=WRAD(J)-DEFL(I,J)
            R2DFL=WRAD(J)**2
            ADFL=DSQRT(R2DFL-WCHGT**2)
            BDFL=DATAN2(WCHGT,ADFL)
            AREA(I,J)=(PI*R2DFL)/2.0-(WCHGT*ADFL+R2DFL*BDFL)
 6004    CONTINUE
 6003    CONTINUE
      ENDIF
C
C     READ THE INDICATOR FOR THE MEXICAN BUMP FROM ROAD FILE
C     bh------------------------------------------begin read ROAD data
C
      IF(IVERT(J).EQ.3 .and. roadFiles(j) /= 0 ) THEN
         rewind(roadFiles(j))
         READ(roadFiles(j),*) IMEX
         WRITE(outFile,20) 'IMEX :',IMEX
         IF(IMEX.EQ.1) THEN
C
            READ(roadFiles(j),*)STA(1,J),STA(25,J)
            WRITE(outFile,19)'STA first and last :',
     1           STA(1,J),STA(25,J)
 19         FORMAT(A,X,1p,6E12.3,/,(10X,6F20.4))
            HGT(1,J)=0.0D0
            HGT(25,J)=0.0D0
C
            READ(roadFiles(j),*)XZERO,YZERO,TTAD,YBAR
            WRITE(outFile,19)'XZERO,YZERO,TTAD,YBAR :',
     1           XZERO,YZERO,TTAD,YBAR
            TTAR=TTAD*DTR
C
            READ(roadFiles(j),*)AL1,AL2,A1,A2
            WRITE(outFile,19)'AL1,AL2,A1,A2 :', AL1,AL2,A1,A2
C
            NTSTA(J)=25
            WRITE(outFile,20)'NTSTA :',NTSTA(J)
            NPTS=NTSTA(J)
C
            HGT(2,J)=0.0D0
            X2=-(AL1+AL2)
            Y2=(YBAR-YZERO-X2*DSIN(TTAR))/DCOS(TTAR)
            STA(2,J)=XZERO+X2*DCOS(TTAR)-Y2*DSIN(TTAR)
C
            HGT(3,J)=-A1
            X3=-AL2
            Y3=(YBAR-YZERO-X3*DSIN(TTAR))/DCOS(TTAR)
            STA(3,J)=XZERO+X3*DCOS(TTAR)-Y3*DSIN(TTAR)
C
            DELTAX=AL2/10.0D0
            DO 99 K=4,23
               XK=-AL2+(K-3)*DELTAX
               ABSXK=DABS(XK)
               ARG=(PI*ABSXK)/AL2
               HGT(K,J)=(A2-A1)-(A2/2.0D0)*(1.0D0-DCOS(ARG))
               YK=(YBAR-YZERO-XK*DSIN(TTAR))/DCOS(TTAR)
               STA(K,J)=XZERO+XK*DCOS(TTAR)-YK*DSIN(TTAR)
 99         CONTINUE
C
            HGT(24,J)=0.0D0
            X24=AL1+AL2
            Y24=(YBAR-YZERO-X24*DSIN(TTAR))/DCOS(TTAR)
            STA(24,J)=XZERO+X24*DCOS(TTAR)-Y24*DSIN(TTAR)
C
            write(outfile,"('------ STA ----- HGT ------')")
            do k = 1,npts
               write(outfile,"(1p,e12.3,'  ',e12.3)") sta(k,j),hgt(k,j)
            end do
            write(outfile,"('--------')")

         ENDIF
C
         IF(IMEX.EQ.0) THEN
C
C     NTSTA(J) = NUMBER OF TERRAIN STATIONS FOR WHEEL J
C
            READ(roadFiles(j),*)NTSTA(J)
            WRITE(outFile,20)'NTSTA :',NTSTA(J)
            NPTS=NTSTA(J)

            if ( npts > size(sta,1) ) then
               WRITE(outFile,
     +              "('Error :   Number of road points to large:'
     +              ,i10,' > max size ',i10 )") npts, size(sta,1)
               ierr = -1
               return
            end if
C
C     STA(K,J) = ARRAY OF TERRAIN STATIONS FOR TIRE J
C     HGT(K,J) = ARRAY OF CORRESPONDING ELEVATIONS FOR TIRE J
C
            READ(roadFiles(j),*)(STA(K,J),K=1,NPTS)
            READ(roadFiles(j),*)(HGT(K,J),K=1,NPTS)
            WRITE(outFile,19)'HGT :',(HGT(K,J),K=1,NPTS)

            write(outfile,"('------ STA ----- HGT ------')")
            do k = 1,npts
               write(outfile,"(1p,e12.3,'  ',e12.3)") sta(k,j),hgt(k,j)
            end do
            write(outfile,"('--------')")

CbhCbh  Remove temp
Cbh            do ii = j+1,size(tireFiles)
Cbh               ntsta(ii) = npts
Cbh#ifdef FT_DEBUG
Cbh               write(dbgTire,"('Setting NTSTA(',i2,') = ',I5)")
Cbh     *              ii, ntsta(ii)
Cbh#endif
Cbh               write(*      ,"('Setting NTSTA(',i2,') = ',I5)")
Cbh     *              ii, ntsta(ii)
Cbh               do k = 1,npts
Cbh                  sta(k,ii) = sta(k,j)
Cbh                  hgt(k,ii) = hgt(k,j)
Cbh               end do
Cbh            end do

         ENDIF
      ENDIF

C     bh------------------------------------------end read ROAD data
C
C     ILAT(J)   = 1 - LATERAL SPRING AND DAMPER (RIDE ANALYSIS)
C     2 - EQUATION FOR LATERAL FORCE-SLIP ANGLE RELATION
C     3 - CARPET PLOT DATA FOR LATERAL FORCE AS FUNCTION OF SLIP ANGLE AND VERTICAL FORCE
C     AKL(J) = LATERAL SPRING RATE OF TIRE
C     CLAT(J) = LATERAL DAMPING RATE OF TIRE
C     YSE(J) = LATERAL POSITION OF WHEEL CENTER FOR STATIC EQUILIBRIUM
C
C     IYFRIC(J) = 0 - DON'T CONSIDER FRICTION LIMIT ON LATERAL FORCE
C     = 1 - DO CONSIDER FRICTION LIMIT ON LATERAL FORCE
C     AMUL(J) = LIMITING FRICTION COEFFICIENT FOR LATERAL FORCE
C     AMUL(J) = MAXIMUM LATERAL FORCE COEFFICIENT
C     ALEXP(J) = EXPONENT IN LATERAL FORCE COEFF VS SLIP ANGLE RELATION
C     IROLLM(J) = 0 - DON'T CONSIDER ROLL MOMENT RESULTING FROM VERTICAL
C     FORCE AND LATERAL DEFLECTION OF TIRES IN HANDLING SITUATION
C     = 1 - DO CONSIDER THE ABOVE
C     AKLAT(J) = LATERAL STIFFNESS USED TO FIND LATERAL DEFLECTION
      call read_StiffDampLat_nml(TireFiles(j),
     *     ILAT(J),AKL(J),CLAT(J),YSE_from_fedem(J), YSE(J),
     *     IYFRIC(J),AMUL(J),ALEXP(J),
     *     IROLLM(J),AKLAT(J))
C
      if ( YSE_from_fedem(J) ) then
        write(outfile,"('YSE will be initialized from fedem data')")
      else
        write(outfile,"('Keeping YSE from tire data file')")
      end if

      yFricCoeff(j) = amul(j)
      yFricExp(j)   = alexp(j)
C
C     Spline number from ADAMS is not possible
      IF(ILAT(J).EQ.3) then
         ierr = -1
         return
      end if
C
C     AKFA(J) = FORE-AFT SPRING RATE OF TIRE
C     CFA(J) = FORE-AFT DAMPING RATE OF TIRE
C     XSE(J) = FORE-AFT POSITION OF WHEEL CENTER FOR STATIC EQUILIBRIUM
C     WRT REAR AXLE
C     XARM(J) = MOMENT ARM OF FORCE - POSITIVE IF FORCE ACTS BELOW WHEEL CENTER
C     IXFRIC(J) = 0 - DON'T CONSIDER FRICTION LIMIT ON FORE-AFT FORCE
C     = 1 - DO CONSIDER FRICTION LIMIT ON FORE-AFT FORCE
C     AMUFA(J) = LIMITING FRICTION COEFFICIENT FOR FORE-AFT FORCE/
C
      call read_StiffDampLong_nml(TireFiles(j),
     *     ILONG(J),AKFA(J),CFA(J),XSE(J),XARM(J),
     *     IXFRIC(J),AMUFA(J))

C
C     READ THE NUMBER OF SHOCK ABSORBERS
C
      NSHOCK = 0

      if ( j == 1 ) then

         nmlfile = findUnitNumber(60)
         open( UNIT=nmlFile, FILE="JDT_input_example.nml")

Cbh         call write_differential_data_nml(nmlFile,angvel(j),
Cbh     *        alat(j),trdw(j))

         call write_tire_radius_data_nml(nmlFile,wrad(j),yarm(j),rr(j))

         call write_surface_data_nml(nmlFile,IROTAT(J),IPOWER(J),
     *        ISURF(j),ICONST(j),CI(J),B(J),D(J),H(J),DELTA(J))

         call write_vertDamping_nml(nmlFile,CDAMPZ(J),CDEFL(J))

         call write_vertStiff_nml(nmlFile,ivert(j),akv(j),defl(1,j),
     *        force(1,j),actReg(j),nLdPt(j),anseg(j))

         call write_StiffDampLat_nml(nmlFile,
     *        ILAT(J),AKL(J),CLAT(J),YSE_from_fedem(J),YSE(J),
     *        IYFRIC(J),AMUL(J),ALEXP(J),
     *        IROLLM(J),AKLAT(J))

         call write_StiffDampLong_nml(nmlFile,
     *        ILONG(J),AKFA(J),CFA(J),XSE(J),XARM(J),
     *        IXFRIC(J),AMUFA(J))

         close( UNIT=nmlFile )

      end if

C
      END
C================END SUBROUTINE READINPUT_FDM_NML


      SUBROUTINE READINPUT_FDM(ID,TIME,PAR,NPAR,ierr)

      use JD_namelistModule
      use jdt_inp
      use jdt_var
      use jdt_fdm
      use fileUtilitiesModule

      implicit none

      integer id, npar, ierr
      real*8 TIME, PAR(*)
C
      integer IFPASS
      real*8 RTOD
      real*8 PI
      real*8 POT
      real*8 DTR
      DATA IFPASS /0/
      DATA RTOD /57.29578D0/
      DATA PI/3.1415926535897932D0/
      DATA POT/1.5707963267948D0/
      DATA DTR/0.01745329/

      integer i, j, k, npts, nseg, nsegd2, nsegp1, imex, nmlFile
      real*8 rseg1, dseg1, rseg, segang, wchgt, r2dfl, adfl, bdfl
      real*8 xzero, yzero, ttad, ybar, al1, al2, a1, a2, ttar, x2, y2
      real*8 x3, y3, xk, yk, absxk, arg, x24, y24, deltax
C
      logical isOpen
C
C     === Logic section
      ierr = 0

C     Dissever tire number from input identifier
      J = ID - (ID/100)*100
      rewind(tireFiles(j))
C
      inquire(FILE="TireEcho.nml", OPENED=isOpen)
      if ( .not. isOpen ) then
         open( UNIT=echoFile, FILE="TireEcho.nml")
      endif
C
C     READ THE TIRE DATA FILE ON THE FIRST PASS THROUGH THIS SUBROUTINE
C
      WRITE(outFile,2)
 2    FORMAT(//,3X,' ======== TIRE DATA FILE ========',/)
      write(echoFile,"('======== Tire data for tire number',i3)") j
c
C     IENG = 0 - DON'T CONSIDER ENGINE TORQUE-SPEED RELATION
C     1 - DO CONSIDER ENGINE TORQUE-SPEED RELATION
C
      IENG = 0
      WRITE(outFile,16)IENG
C
C     IRAMKR = MARKER AT REAR AXLE CENTER
C
      IGMKR  = 0
      IRAMKR = 0
C     bh Set IRAMKR, IGMKR
C     bh      READ(TireFiles(j),*) IRAMKR
C     bh      WRITE(outFile,16) IRAMKR
 16   FORMAT(10X,6I5)
C
C     GROUND ORIGIN MARKER
C
C     bh      READ(TireFiles(j),*) IGMKR
C     bh      WRITE(outFile,16) IGMKR
C
C     READ THE NUMBER OF TIRES
C
      NTIRES = NTIRES + 1
      WRITE (outFile,"('Number of tires so far ',i3)") NTIRES
C
C     READ APPROPRIATE DATA FOR EACH TIRE
C
C     IMRKR(J) = MARKER AT WHICH FORCE ACTUALLY ACTS
C
      IMRKR(J) = J
      WRITE(outFile,20)'IMRKR',IMRKR(J)
 20   FORMAT(A,X,10I5)
C
C     IWMKR(J) = MARKER AT CENTER OF WHEEL
C     IVMKR(J) = REFERENCE MARKER ON PART TO WHICH WHEEL IS ATTACHED
C
      IWMKR(J) = J
      IVMKR(J) = J
      WRITE(outFile,20)'IWMKR IVMKR',IWMKR(J),IVMKR(J)
C
C     WRAD(J) = UNDEFLECTED TIRE RADIUS
C
      write(*,"('tire number ',I2,' input file number ',I3)")
     *      J, TireFiles(j)
      READ(TireFiles(j),*)WRAD(J)
      WRITE(outFile,19)'WRAD :',WRAD(J)
 19   FORMAT(A,X,1p,6E12.3,/,(10X,6F20.4))
      rollRadius(J) = wrad(j)
C
C     CDAMPZ(J) = VERTICAL DAMPING RATE OF TIRE
C     CDEFL(J) = TIRE DEFLECTION AT WHICH DAMPING REACHES MAXIMUM VALUE
C
      READ(TireFiles(j),*)CDAMPZ(J),CDEFL(J)
      WRITE(outFile,19)'CDAMP,CDEFL :',CDAMPZ(J),CDEFL(J)
C
C     YARM(J) = MOMENT ARM OF LATERAL FORCE - POSITIVE IF FORCE ACTS
C     BELOW WHEEL CENTER
C
      READ(TireFiles(j),*)YARM(J)
      WRITE(outFile,19)'YARM :',YARM(J)
C
C     IROTAT(J) = 0 - ROTATION NOT CONSIDERED
C     1 - ROTATION CONSIDERED USING KINEMATIC MODEL OF DIFFERENTIAL
C     2 - ROTATION DETERMINED FROM MODEL
C     ISTEER(J) = 0 - NO STEERING INPUT (RIDE ANALYSIS)
C     1 - USER INPUT - WHEEL DOESN'T PHYSICALLY STEER IN MODEL
C     2 - STEER ANGLE DETERMINED FROM MODEL
C     IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C     1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C     2 - POWERED WHEEL - CONSIDER TRACTION & ROLLING RESISTANCE
C     IVERT(J)  = 1 - LINEAR SPRING MODEL FOR VERTICAL FORCE
C     2 - FORCE CALCULATED USING LOAD-DEFLECTION RELATION
C     3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C     ILAT(J)   = 1 - LATERAL SPRING AND DAMPER (RIDE ANALYSIS)
C     2 - EQUATION FOR LATERAL FORCE-SLIP ANGLE RELATION
C     3 - CARPET PLOT DATA FOR LATERAL FORCE AS FUNCTION OF
C     SLIP ANGLE AND VERTICAL FORCE
C     IROLLM(J) = 0 - DON'T CONSIDER ROLL MOMENT RESULTING FROM VERTICAL
C     FORCE AND LATERAL DEFLECTION OF TIRES IN HANDLING
C     SITUATION
C     = 1 - DO CONSIDER THE ABOVE
C     ILONG(J)  = 0 - NO FORE-AFT SPRING AND DAMPER
C     1 - USE FORE-AFT SPRING AND DAMPER

      READ(TireFiles(j),*)IROTAT(J),ISTEER(J),IPOWER(J),
     1     IVERT(J),ILAT(J),IROLLM(J),ILONG(J)
      WRITE(outFile,20)
     1     'IROTAT,ISTEER,IPOWER,IVERT,ILAT,IROLLM,ILONG :',
     2     IROTAT(J),ISTEER(J),IPOWER(J),IVERT(J),ILAT(J),
     3     IROLLM(J),ILONG(J)
C
C     INPUT DATA IF KINEMATIC DIFFERENTIAL MODEL IS USED
C
C     ANGVEL(J) = ANGULAR VELOCITY OF AXLE (RAD/SEC)
C     ALAT(J) = LATERAL POSITION OF WHEEL ON AXLE (IN VEHICLE SYSTEM)
C     TRDW(J) = EFFECTIVE TREAD WIDTH OF WHEELS ON AXLE
C     RR(J) = ROLLING RADIUS
C
      IF(IROTAT(J).EQ.1) then
         READ(TireFiles(j),*)ANGVEL(J),ALAT(J),TRDW(J),RR(J)
         rollRadius(J) = rr(j)
         WRITE(outFile,19)'ANGVEL,ALAT,TRDW,RR :',
     *        ANGVEL(J),ALAT(J),TRDW(J),RR(J)
         call write_differential_data_nml(echoFile,angvel(j),
     *        alat(j),trdw(j))
      end if
C
C     INPUT DATA IF WHEEL ANGULAR VELOCITIES ARE FOUND FROM MODEL
C
C     IWROT(J) = I MARKER MAKING UP WHEEL-CHASSIS REVOLUTE JOINT
C     JWROT(J) = J MARKER MAKING UP WHEEL-CHASSIS REVOLUTE JOINT
C
      IF(IROTAT(J).EQ.2) then
         iwrot(j) = j
         jwrot(j) = 0
         IF(IROTAT(J).EQ.2)WRITE(outFile,20)'IWROT,JWROT :',
     1        IWROT(J),JWROT(J)
      end if
C
C     IWCMKR(J) = MARKER USED TO CALCULATE STEER ANGLE
C
      IF(ISTEER(J).EQ.2)READ(TireFiles(j),*)IWCMKR(J)
      IF(ISTEER(J).EQ.2)WRITE(outFile,20)'IWCMKR :',IWCMKR(J)
C
C     INPUT IF TRACTION AND/OR ROLLING RESISTANCE ARE CONSIDERED
C
C     ISURF(J) = 1 - SOIL SURFACE
C     = 2 - HARD SURFACE
C
      IF(IPOWER(J).NE.0)READ(TireFiles(j),*)ISURF(J)
      IF(IPOWER(J).NE.0)WRITE(outFile,20)'ISURF :',ISURF(J)
C
C     ICONST(J) = 1 - BIAS TIRE CONSTRUCTION
C     = 2 - RADIAL TIRE CONSTRUCTION
C
      IF(IPOWER(J).NE.0)READ(TireFiles(j),*)ICONST(J)
      IF(IPOWER(J).NE.0)WRITE(outFile,20)'ICONST :',ICONST(J)
C
C     CI(J) = CONE INDEX (FOR SOIL) OR K (FOR HARD SURFACE)
C     B(J) = UNLOADED TIRE SECTION WIDTH
C     D(J) = UNLOADED TIRE DIAMETER
C     RR(J) = TIRE ROLLING RADIUS
C     H(J) = TIRE SECTION HEIGHT
C     DELTA(J) = LOADED TIRE DEFLECTION
C
      IF(IPOWER(J).EQ.1.AND.ISURF(J).EQ.1) then
         READ(TireFiles(j),*)CI(J),B(J),D(J),H(J),DELTA(J)
         WRITE(outFile,19)'CI,B,D,H,DELTA :',
     1        CI(J),B(J),D(J),H(J),DELTA(J)
      end if
C
      IF(IPOWER(J).EQ.2.AND.ISURF(J).EQ.1) then
         READ(TireFiles(j),*)CI(J),B(J),D(J),H(J),DELTA(J),RR(J)
         WRITE(outFile,19)'CI,B,D,H,DELTA,RR :',
     1        CI(J),B(J),D(J),H(J),DELTA(J),RR(J)
         rollRadius(J) = rr(j)
      end if
C
      IF(IPOWER(J).EQ.2.AND.ISURF(J).EQ.2) then
         READ(TireFiles(j),*)CI(J),B(J),D(J),RR(J)
         WRITE(outFile,19)'CI,B,D,RR :',
     1        CI(J),B(J),D(J),RR(J)
         rollRadius(J) = rr(j)
      end if

      call write_tire_radius_data_nml(echoFile,wrad(j),yarm(j),rr(j))
      call write_surface_data_nml(echoFile,IROTAT(J),IPOWER(J),
     *     ISURF(j),ICONST(j),CI(J),B(J),D(J),H(J),DELTA(J))
C
C     INPUT DEPENDING ON VERTICAL TIRE FORCE MODEL
C
C
C     AKV(J) = VERTICAL SPRING RATE OF TIRE
C
      IF(IVERT(J).EQ.1) then
         READ(TireFiles(j),*)AKV(J)
         WRITE(outFile,19)'AKV :',AKV(J)
      end if
C
C     NLDPT(J) = NUMBER OF LOAD-DEFLECTION POINTS
C
      IF(IVERT(J).NE.1)then
         READ(TireFiles(j),*)NLDPT(J)
         WRITE(outFile,20)'NLDPT :',NLDPT(J)
      endif
      NPTS=NLDPT(J)
C
C     DEFL(I,J) = ARRAY OF TIRE DEFLECTIONS FOR TIRE J
C     FORCE(I,J) = ARRAY OF CORRESPONDING TIRE VERTICAL LOADS FOR TIRE J
C
      IF(IVERT(J).NE.1)READ(TireFiles(j),*)(DEFL(I,J),I=1,NPTS)
      IF(IVERT(J).NE.1)WRITE(outFile,19)'DEFL :',
     1     (DEFL(I,J),I=1,NPTS)
C
      IF(IVERT(J).NE.1)READ(TireFiles(j),*)(FORCE(I,J),I=1,NPTS)
      IF(IVERT(J).NE.1)WRITE(outFile,19)'FORCE :',
     1     (FORCE(I,J),I=1,NPTS)
C
C     ANSEG(J) = NUMBER OF SEGMENTS IN ACTIVE REGION
C     ACTREG(J) = SIZE OF ACTIVE REGION (IN DEGREES)
C
      IF(IVERT(J).EQ.3)READ(TireFiles(j),*)ANSEG(J),ACTREG(J)
      IF(IVERT(J).EQ.3)WRITE(outFile,19)'ANSEG,ACTREG',
     1     ANSEG(J),ACTREG(J)

C     Write echo of input to namelist format
      call write_vertDamping_nml(echoFile,CDAMPZ(J),CDEFL(J))
      call write_vertStiff_nml(echoFile,ivert(j),akv(j),defl(1,j),
     *     force(1,j),actReg(j),nLdPt(j),anseg(j))
C
C     CALCULATE SOME INITIAL VALUES FOR THE DISTRIBUTED CONTACT MODEL
C
      IF(IVERT(J).EQ.3) THEN
         ITR(J)=2
         NSEG=ANSEG(J)
         SEGANG=ACTREG(J)/ANSEG(J)
         DSEG1=-(ACTREG(J)-SEGANG)/2.0
         RSEG1=DTR*DSEG1
         SSEGA(1,J)=DSIN(RSEG1)
         CSEGA(1,J)=DCOS(RSEG1)
         SM(1,J)=-CSEGA(1,J)/SSEGA(1,J)
         DO 101 I=2,NSEG
            RSEG=DTR*(DSEG1+(I-1)*SEGANG)
            SSEGA(I,J)=DSIN(RSEG)
            CSEGA(I,J)=DCOS(RSEG)
            SM(I,J)=-CSEGA(I,J)/SSEGA(I,J)
 101     CONTINUE
C     TRANSLATE THE TIRE MODEL INDICATOR
C     NOTE: LTIRE=0 - INDIVIDUAL SPRING MODEL
C     LTIRE=1 - DEFLECTED AREA MODEL
C     NOTE THE TIREF SUBROUTINE IMPLEMENTS BOTH MODELS BUT HERE THE
C     DEFLECTED AREA MODEL IS THE ONLY ONE USED
         LTIRE=1
         IF(LTIRE.EQ.1) GO TO 6001
C     CALCULATIONS RELATED TO THE ORIGINAL TIRE MODEL
         NSEGD2=NSEG/2
         NSEGP1=NSEGD2+1
C     READ IN THE SEGMENT SPRING RATES FOR THE LEFT HALF OF THE ACTIVE
C     REGION
         READ(TireFiles(j),*) (ASEGK(I,J),I=1,NSEGD2)
C     SET THE SEGMENT SPRING RATES FOR THE RIGHT HALF OF THE REGION
         DO 6002 I=NSEGP1,NSEG
            ASEGK(I,J)=ASEGK(NSEG-I+1,J)
 6002    CONTINUE
         GO TO 6003
C     CALCULATIONS RELATED TO THE NEW TIRE MODEL
 6001    CONTINUE
C     DETERMINE THE DEFLECTED AREA CORRESPONDING TO EACH DEFLECTION
         DTHET(J)=SEGANG*DTR
         AREA(1,J)=0.0
         DO 6004 I=2,NPTS
            WCHGT=WRAD(J)-DEFL(I,J)
            R2DFL=WRAD(J)**2
            ADFL=DSQRT(R2DFL-WCHGT**2)
            BDFL=DATAN2(WCHGT,ADFL)
            AREA(I,J)=(PI*R2DFL)/2.0-(WCHGT*ADFL+R2DFL*BDFL)
 6004    CONTINUE
 6003    CONTINUE
      ENDIF
C
C     READ THE INDICATOR FOR THE MEXICAN BUMP FROM ROAD FILE
C     bh  This should be read from a separate file
C     bh------------------------------------------begin read ROAD data
C
      IF(IVERT(J).EQ.3 .and. roadFiles(j) /= 0 ) THEN
         rewind(roadFiles(j))
         READ(roadFiles(j),*) IMEX
         WRITE(outFile,20) 'IMEX :',IMEX
         IF(IMEX.EQ.1) THEN
C
            READ(roadFiles(j),*)STA(1,J),STA(25,J)
            WRITE(outFile,19)'STA first and last :',
     1           STA(1,J),STA(25,J)
            HGT(1,J)=0.0D0
            HGT(25,J)=0.0D0
C
            READ(roadFiles(j),*)XZERO,YZERO,TTAD,YBAR
            WRITE(outFile,19)'XZERO,YZERO,TTAD,YBAR :',
     1           XZERO,YZERO,TTAD,YBAR
            TTAR=TTAD*DTR
C
            READ(roadFiles(j),*)AL1,AL2,A1,A2
            WRITE(outFile,19)'AL1,AL2,A1,A2 :', AL1,AL2,A1,A2
C
            NTSTA(J)=25
            WRITE(outFile,20)'NTSTA :',NTSTA(J)
            NPTS=NTSTA(J)
C
            HGT(2,J)=0.0D0
            X2=-(AL1+AL2)
            Y2=(YBAR-YZERO-X2*DSIN(TTAR))/DCOS(TTAR)
            STA(2,J)=XZERO+X2*DCOS(TTAR)-Y2*DSIN(TTAR)
C
            HGT(3,J)=-A1
            X3=-AL2
            Y3=(YBAR-YZERO-X3*DSIN(TTAR))/DCOS(TTAR)
            STA(3,J)=XZERO+X3*DCOS(TTAR)-Y3*DSIN(TTAR)
C
            DELTAX=AL2/10.0D0
            DO 99 K=4,23
               XK=-AL2+(K-3)*DELTAX
               ABSXK=DABS(XK)
               ARG=(PI*ABSXK)/AL2
               HGT(K,J)=(A2-A1)-(A2/2.0D0)*(1.0D0-DCOS(ARG))
               YK=(YBAR-YZERO-XK*DSIN(TTAR))/DCOS(TTAR)
               STA(K,J)=XZERO+XK*DCOS(TTAR)-YK*DSIN(TTAR)
 99         CONTINUE
C
            HGT(24,J)=0.0D0
            X24=AL1+AL2
            Y24=(YBAR-YZERO-X24*DSIN(TTAR))/DCOS(TTAR)
            STA(24,J)=XZERO+X24*DCOS(TTAR)-Y24*DSIN(TTAR)
C
            write(outfile,"('------ STA ----- HGT ------')")
            do k = 1,npts
               write(outfile,"(1p,e12.3,'  ',e12.3)") sta(k,j),hgt(k,j)
            end do
            write(outfile,"('--------')")

         ENDIF
C
         IF(IMEX.EQ.0) THEN
C
C     NTSTA(J) = NUMBER OF TERRAIN STATIONS FOR WHEEL J
C
            READ(roadFiles(j),*)NTSTA(J)
            WRITE(outFile,20)'NTSTA :',NTSTA(J)
            NPTS=NTSTA(J)

            if ( npts > size(sta,1) ) then
               WRITE(outFile,
     +              "('Error :   Number of road points to large:'
     +              ,i10,' > max size ',i10 )") npts, size(sta,1)
               ierr = -1
               return
            end if
C
C     STA(K,J) = ARRAY OF TERRAIN STATIONS FOR TIRE J
C     HGT(K,J) = ARRAY OF CORRESPONDING ELEVATIONS FOR TIRE J
C
            READ(roadFiles(j),*)(STA(K,J),K=1,NPTS)
            READ(roadFiles(j),*)(HGT(K,J),K=1,NPTS)

            write(outfile,"('------ STA ----- HGT ------')")
            do k = 1,npts
               write(outfile,"(1p,e12.3,'  ',e12.3)") sta(k,j),hgt(k,j)
            end do
            write(outfile,"('--------')")

CbhCbh  Remove temp
Cbh            do ii = j+1,size(tireFiles)
Cbh               ntsta(ii) = npts
Cbh#ifdef FT_DEBUG
Cbh               write(dbgTire,"('Setting NTSTA(',i2,') = ',I5)")
Cbh     *                         ii, ntsta(ii)
Cbh#endif
Cbh               write(*      ,"('Setting NTSTA(',i2,') = ',I5)")
Cbh     *                         ii, ntsta(ii)
Cbh               do k = 1,npts
Cbh                  sta(k,ii) = sta(k,j)
Cbh                  hgt(k,ii) = hgt(k,j)
Cbh               end do
Cbh            end do

         ENDIF
      ENDIF

C     bh------------------------------------------end read ROAD data
C
C     AKL(J) = LATERAL SPRING RATE OF TIRE
C     CLAT(J) = LATERAL DAMPING RATE OF TIRE
C     YSE(J) = LATERAL POSITION OF WHEEL CENTER FOR STATIC EQUILIBRIUM
C
      IF(ILAT(J).EQ.1)READ(TireFiles(j),*)AKL(J),CLAT(J),YSE(J)
      IF(ILAT(J).EQ.1)WRITE(outFile,19)'AKL,CLAT,YSE :',
     1     AKL(J),CLAT(J),YSE(J)
C
C     IYFRIC(J) = 0 - DON'T CONSIDER FRICTION LIMIT ON LATERAL FORCE
C     = 1 - DO CONSIDER FRICTION LIMIT ON LATERAL FORCE
C
      IF(ILAT(J).EQ.1)READ(TireFiles(j),*)IYFRIC(J)
      IF(ILAT(J).EQ.1)WRITE(outFile,20)'IYFRIC :',IYFRIC(J)
C
C     AMUL(J) = LIMITING FRICTION COEFFICIENT FOR LATERAL FORCE
C
      IF(ILAT(J).EQ.1.AND.IYFRIC(J).EQ.1)READ(TireFiles(j),*)AMUL(J)
      IF(ILAT(J).EQ.1.AND.IYFRIC(J).EQ.1)
     1     WRITE(outFile,19)'AMUL :',AMUL(J)
C
C     AMUL(J) = MAXIMUM LATERAL FORCE COEFFICIENT
C     ALEXP(J) = EXPONENT IN LATERAL FORCE COEFF VS SLIP ANGLE RELATION
C
      IF(ILAT(J).EQ.2)READ(TireFiles(j),*)AMUL(J),ALEXP(J)
      IF(ILAT(J).EQ.2)WRITE(outFile,19)
     1     'AMUL,ALEXP :',AMUL(J),ALEXP(J)

      yFricCoeff(j) = amul(j)
      yFricExp(j)   = alexp(j)
C
C     ISPLIN(J) = SPLINE NUMBER TO USE FOR CALCULATING LATERAL FORCE
C
      IF(ILAT(J).EQ.3)READ(TireFiles(j),*)ISPLIN(J)
      IF(ILAT(J).EQ.3)WRITE(outFile,20)'ISPLIN :',ISPLIN(J)
C
C     AKLAT(J) = LATERAL STIFFNESS USED TO FIND LATERAL DEFLECTION
C
      IF(IROLLM(J).EQ.1)READ(TireFiles(j),*)AKLAT(J)
      IF(IROLLM(J).EQ.1)WRITE(outFile,19)'AKLAT :',AKLAT(J)
C
      YSE_from_fedem(J) = .false.
      call write_StiffDampLat_nml(echoFile,
     *     ILAT(J),AKL(J),CLAT(J),YSE_from_fedem(J),YSE(J),
     *     IYFRIC(J),AMUL(J),ALEXP(J),
     *     IROLLM(J),AKLAT(J))

C     AKFA(J) = FORE-AFT SPRING RATE OF TIRE
C     CFA(J) = FORE-AFT DAMPING RATE OF TIRE
C     XSE(J) = FORE-AFT POSITION OF WHEEL CENTER FOR STATIC EQUILIBRIUM
C     WRT REAR AXLE
C     XARM(J) = MOMENT ARM OF FORCE - POSITIVE IF FORCE ACTS BELOW
C     WHEEL CENTER
C
      IF(ILONG(J).EQ.1)
     1     READ(TireFiles(j),*)AKFA(J),CFA(J),XSE(J),XARM(J)
      IF(ILONG(J).EQ.1)
     1     WRITE(outFile,19) 'AKFA,CFA,XSE,XARM :',
     2     AKFA(J),CFA(J),XSE(J),XARM(J)
C
C     IXFRIC(J) = 0 - DON'T CONSIDER FRICTION LIMIT ON FORE-AFT FORCE
C     = 1 - DO CONSIDER FRICTION LIMIT ON FORE-AFT FORCE
C
      IF(ILONG(J).EQ.1)READ(TireFiles(j),*)IXFRIC(J)
      IF(ILONG(J).EQ.1)WRITE(outFile,20)'IXFRIC :',IXFRIC(J)
C
C     AMUFA(J) = LIMITING FRICTION COEFFICIENT FOR FORE-AFT FORCE/
C
      IF(ILONG(J).EQ.1.AND.IXFRIC(J).EQ.1)
     1     READ(TireFiles(j),*)AMUFA(J)
      IF(ILONG(J).EQ.1.AND.IXFRIC(J).EQ.1)
     1     WRITE(outFile,19)'AMUFA',AMUFA(J)

      call write_StiffDampLong_nml(echoFile,
     *     ILONG(J),AKFA(J),CFA(J),XSE(J),XARM(J),
     *     IXFRIC(J),AMUFA(J))
C
C
C     READ THE NUMBER OF SHOCK ABSORBERS
C
      NSHOCK = 0

      if ( j == 1 ) then

         nmlfile = findUnitNumber(60)
         open( UNIT=nmlFile, FILE="JDT_input_example.nml")

Cbh         call write_differential_data_nml(nmlFile,angvel(j),
Cbh     *        alat(j),trdw(j))

         call write_tire_radius_data_nml(nmlFile,wrad(j),yarm(j),rr(j))

         call write_surface_data_nml(nmlFile,IROTAT(J),IPOWER(J),
     *        ISURF(j),ICONST(j),CI(J),B(J),D(J),H(J),DELTA(J))

         call write_vertDamping_nml(nmlFile,CDAMPZ(J),CDEFL(J))

         call write_vertStiff_nml(nmlFile,ivert(j),akv(j),defl(1,j),
     *        force(1,j),actReg(j),nLdPt(j),anseg(j))

         call write_StiffDampLat_nml(nmlFile,
     *        ILAT(J),AKL(J),CLAT(J),YSE_from_fedem(J),YSE(J),
     *        IYFRIC(J),AMUL(J),ALEXP(J),
     *        IROLLM(J),AKLAT(J))

         call write_StiffDampLong_nml(nmlFile,
     *        ILONG(J),AKFA(J),CFA(J),XSE(J),XARM(J),
     *        IXFRIC(J),AMUFA(J))

         close( UNIT=nmlFile )

      end if

C
      END
C================END SUBROUTINE READINPUT_FDM

C
C     ADAMS USER FORCE SUBROUTINE IMPLEMENTING OFF ROAD TIRE MODEL
C     DAVID SMITH - DEERE & COMPANY TECHNICAL CENTER
C     SEPTEMBER 1991
C
C*********************************************************************
C
      SUBROUTINE SFOSUB(ID,TIME,PAR,NPAR,DFLAG,IFLAG,FMAG,ierr)

      use JD_namelistModule
      use jdt_inp
      use jdt_fdm
      use jdt_var

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION PAR(*), DATA(6)
      DIMENSION VECTOR(10), SPLVAL(3)
      LOGICAL DFLAG, ERRFLG, IFLAG
C
      DATA IFPASS /0/
      DATA RTOD /57.29578D0/
      DATA PI/3.1415926535897932D0/
      DATA POT/1.5707963267948D0/
      DATA DTR/0.01745329/

      real*8 zeroDelta(3)

      zeroDelta(1) = 0.0d0
      zeroDelta(2) = 0.0d0
      zeroDelta(3) = 0.0d0
C
C     Dissever tire number from input identifier

      j = id - (id/100)*100
      if (.not.hasReadTire(j) ) then
         if ( isNameListFile(tireFiles(j)) ) then
            call readinput_fdm_nml(ID,TIME,PAR,NPAR,ierr)
         else
            call readInput_fdm(ID,TIME,PAR,NPAR,ierr)
         endif
         hasReadTire(j) = .true.
         if (ierr < 0 ) return
      end if
C
C
C INITIALIZE THE FORCE MAGNITUDE TO ZERO
C
      FMAG = 0.0D0
C
      IF ( ID .GT. 200 )       GO TO 140
C
C            L O N G I T U D I N A L      F O R C E S
C
C    IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C                1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C                2 - POWERED WHEEL - CONSIDER TRACTION & ROLLING RESISTANCE
C    IVERT(J)  = 3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C    ILONG(J)  = 0 - NO FORE-AFT SPRING AND DAMPER
C                1 - USE FORE-AFT SPRING AND DAMPER
C
C THREE TYPES OF FORE-AFT OR LONGITUDINAL FORCES ARE CONSIDERED IN THIS
C  SECTION.  THESE ARE:
C    1. ROLLING RESISTANCE AND TRACTIVE FORCES
C    2. THE HORIZONTAL FORCE COMPONENT PREDICTED BY THE DISTRIBUTED
C        CONTACT MODEL
C    3. A LINEAR SPRING AND DAMPER IN THE FORE-AFT DIRECTION
C
C  THIS SECTION DEFINES THE ADAMS MARKERS USED TO DETERMINE THE FORCE
C
      IF(IFLAG) THEN
       J=ID-100
       VECTOR(1)=IWMKR(J)
       VECTOR(2)=IGMKR
       VECTOR(3)=IMRKR(J)
       NV=3
        IF(IPOWER(J).EQ.2.AND.IROTAT(J).EQ.1) THEN
         VECTOR(4)=IVMKR(J)
         NV=4
        ENDIF
        IF(IPOWER(J).EQ.2.AND.IROTAT(J).EQ.2) THEN
         VECTOR(4)=IVMKR(J)
         VECTOR(5)=IWROT(J)
         VECTOR(6)=JWROT(J)
         NV=6
        ENDIF
        IF(ILONG(J).EQ.1) THEN
         VECTOR(4)=IRAMKR
         NV=4
        ENDIF
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 70
C
C CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 60
      J=ID-100
C
C ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C  POSITION) WHERE KX=PAR(2) AND DESIRED X POSITION IS PAR(3)
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3))
   60 RETURN
C
C ***************
C
   70 CONTINUE
C
C     DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = ID - 100
C
C IF THE OPTIONS BEING USED DON'T INCLUDE A FORE-AFT FORCE, RETURN
C
      IF(IPOWER(J).EQ.0.AND.IVERT(J).NE.3.AND.ILONG(J).EQ.0) THEN
       FMAG=0.0
       RETURN
      ENDIF
C
C CALCULATE THE VERTICAL FORCE
C
C FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
       X1(J)=DATA(1)
       Y1(J)=DATA(3)
      ENDIF
C
C FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
       VF=0.0D0
       FMAG=0.0D0
       RETURN
      ENDIF
C
C FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
       CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
       VFORC=FZ(J)
C IF THE DISTRIBUTED CONTACT MODEL IS BEING USED, INCLUDE THE FORE-AFT
C  COMPONENT
       FMAG=-FX(J)
       CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C  IF SO, SET FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C IF THEY ARE TO BE CONSIDERED, FIND THE TRACTIVE AND ROLLING RESISTANCE
C   FORCES
C
      IF(IPOWER(J).NE.0) THEN
C
C FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
      CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
      UW = DATA (1)
C
C FOR POWERED WHEELS, FIND THE SLIPPAGE
C
      IF(IPOWER(J).EQ.2) THEN
C
C FIRST FIND THE ANGULAR VELOCITY OF THE WHEEL
C
      IF(IROTAT(J).EQ.1) THEN
C
C ESTIMATE THE ANGULAR VELOCITY OF THE DRIVE WHEEL USING A SIMPLE
C MODEL OF THE DIFFERENTIAL.
C
C FIND THE PITCH AND YAW VELOCITY OF THE CHASSIS
C
       CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
       WY = DATA(5)
       DPHI = DATA (6)
C
C ESTIMATE THE ANGULAR VELOCITY OF THE WHEEL USING A KINEMATIC MODEL
C  OF THE DIFFERENTIAL
C
      OMEGA =ANGVEL(J)-(DSIGN(1.0D0,ALAT(J))*TRDW(J)*DPHI)/(2.0D0*RR(J))
      ENDIF
C
      IF(IROTAT(J).EQ.2) THEN
C
C FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
       CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
       OMEGA=DATA(6)

#ifdef FT_DEBUG
       write(88,"('0 rot vel (3 is wheel speed):',1p,3e12.3)")
     +            data(4),data(5),data(6)
#endif
C
C FIND THE CHASSIS PITCH VELOCITY FROM THE MODEL
C
       CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
       WY = DATA(5)
      ENDIF
C
C FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = RR ( J ) * (OMEGA + WY)
      ABSVWC = DABS ( VELWC )
      ABSUW = DABS ( UW )
C
C FIND THE LONGITUDINAL SLIPPAGE
C
      SLIP = ( ABSVWC - UW ) / DMAX1 ( ABSVWC, ABSUW )
      ASLIP = DABS ( SLIP )
      ENDIF
C
C CALCULATE THE ROLLING RESISTANCE COEFFICIENT, RHO
C
       IF(ISURF(J).EQ.1) THEN
        CN = (CI(J)*B(J)*D(J))/VF
        BN1 = 1.0 + 5.0 * VDEFL/H(J)
        BN2 = 1.0 + 3.0 * B(J)/D(J)
        BN = CN * BN1/BN2
         IF(IPOWER(J).EQ.1.AND.ICONST(J).EQ.1) THEN
          RHO=1.0/BN + 0.04
         ENDIF
         IF(IPOWER(J).EQ.1.AND.ICONST(J).EQ.2) THEN
          RHO=0.9/BN + 0.0325
         ENDIF
         IF(IPOWER(J).EQ.2.AND.ICONST(J).EQ.1) THEN
          RHO=1.0/BN + 0.04 + 0.5*ASLIP/DSQRT(BN)
         ENDIF
         IF(IPOWER(J).EQ.2.AND.ICONST(J).EQ.2) THEN
          RHO=0.9/BN + 0.0325 + 0.5*ASLIP/DSQRT(BN)
         ENDIF
       ENDIF
       IF(ISURF(J).EQ.2.AND.ICONST(J).EQ.1) RHO=0.02
       IF(ISURF(J).EQ.2.AND.ICONST(J).EQ.2) RHO=0.01625
C
C FINALLY CALCULATE THE ROLLING RESISTANCE FORCE
C
      TF = RHO * VF
C
C THE ROLLING RESISTANCE IS ASSUMED TO ACT IN THE DIRECTION OPPOSITE TO
C  THE DIRECTION OF THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL.
C
      FMAG = FMAG + DSIGN ( TF, UW )
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED
C
      IF(IPOWER(J).EQ.2) THEN
C
C FIND THE TRACTIVE COEFFICIENT
C
      IF(ISURF(J).EQ.1) THEN
        IF(ICONST(J).EQ.1) POWER1 = -0.1*BN
        IF(ICONST(J).EQ.2) POWER1 = -0.1*BN
       Q1 = 1.0 - DEXP(POWER1)
        IF(ICONST(J).EQ.1) POWER2 = -7.5*ASLIP
        IF(ICONST(J).EQ.2) POWER2 = -9.5*ASLIP
       Q2 = 1.0 - DEXP(POWER2)
        IF(ICONST(J).EQ.1) AMU = 0.88 * Q1 * Q2 + 0.04
        IF(ICONST(J).EQ.2) AMU = 0.88 * Q1 * Q2 + 0.0325
      ENDIF
      IF(ISURF(J).EQ.2) THEN
       CN=(CI(J)*B(J)*D(J))/VF
       POWER3 = -CN*ASLIP
C      IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.02
       IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3))
       IF(ICONST(J).EQ.1) AMU = DSIGN(AMU,SLIP) + .02
C      IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.01625
       IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3))
       IF(ICONST(J).EQ.2) AMU = DSIGN(AMU,SLIP) + 0.01625
      ENDIF
C
C FINALLY THE TRACTIVE FORCE
C
      TRAC = AMU * VF
C
C COMBINE THE ROLLING RESISTANCE WITH THE TRACTIVE FORCE
C
C     FMAG = FMAG - DSIGN ( TRAC, SLIP )
      FMAG = FMAG - TRAC
      ENDIF
C
      ENDIF
C
C ALLOW THE OPTION OF INCLUDING A FORE-AFT SPRING AND DAMPER
C
      IF(ILONG(J).EQ.1) THEN
C
C FIND THE FORE-AFT POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      XSPACE=DATA(1)
C
C FIND THE FORE-AFT POSITION OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
      CALL INFO('DISP',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
      XRA=DATA(1)
C
C FIND THE FORE-AFT VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      XDOT=DATA(1)
C
C FIND THE FORE-AFT VELOCITY OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
      CALL INFO('VEL',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
      XRADOT=DATA(1)
C
C FIND THE FORE-AFT FORCE
C
      XF=AKFA(J)*(XSPACE-(XRA+XSE(J)))+CFA(J)*(XDOT-XRADOT)
C
C CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
      IF(IXFRIC(J).EQ.1) THEN
       XFS=XF
       REQMU=XF/VF
       REQMU=DABS(REQMU)
       IF(REQMU.GE.AMUFA(J)) THEN
        XF=AMUFA(J)*VF
        XF=DSIGN(XF,XFS)
       ENDIF
      ENDIF
C
C ADD THE FORE-AFT SPRING-DAMPER FORCE TO FORM THE TOTAL FORE-AFT FORCE
C
      FMAG=FMAG+XF
      ENDIF
      RETURN
C
C **************
C
  140 IF ( ID .GT. 300 )  GO TO 230
C
C
C          L A T E R A L       F O R C E S
C
C    ISTEER(J) = 0 - NO STEERING INPUT (RIDE ANALYSIS)
C                1 - USER INPUT - WHEEL DOESN'T PHYSICALLY STEER IN MODEL
C                2 - STEER ANGLE DETERMINED FROM MODEL
C    ILAT(J)   = 1 - LATERAL SPRING AND DAMPER (RIDE ANALYSIS)
C                2 - EQUATION FOR LATERAL FORCE-SLIP ANGLE RELATION
C                3 - CARPET PLOT DATA FOR LATERAL FORCE AS FUNCTION OF
C                    SLIP ANGLE AND LATERAL FORCE
C
C DETERMINE THE LATERAL FORCE ON THE WHEEL
C
C TELL ADAMS THE MARKERS USED TO DETERMINE THE LATERAL FORCE
C
      IF(IFLAG)THEN
       J=ID-200
       VECTOR(1)=IWMKR(J)
       VECTOR(2)=IGMKR
       VECTOR(3)=IMRKR(J)
       NV=3
        IF(ILAT(J).NE.1.AND.ISTEER(J).NE.2) THEN
         VECTOR(4)=IVMKR(J)
         NV=4
        ENDIF
        IF(ILAT(J).NE.1.AND.ISTEER(J).EQ.2) THEN
         VECTOR(4)=IVMKR(J)
         VECTOR(5)=IWCMKR(J)
         NV=5
        ENDIF
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  THIS SECTION CALCULATES A LATERAL FORCE ON THE TIRE TO
C   PUT THE TRACTOR IN STATIC EQUILIBRIUM
C
      IF ( TIME .GT. 0.0D0 ) GO TO 160
      IF(PAR(1).NE.1) GO TO 150
       J=ID-200
C
C ASSUME LINEAR SPRING IN LATERAL DIRECTION
C  (FMAG=KY*(DATA(2)-DESIRED Y POSITION) WHERE KY=PAR(2) AND DESIRED Y
C   POSITION IS PAR(3))
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG =  - PAR(2) * ( DATA ( 2 ) - PAR ( 3 ) )
  150 RETURN
C
C *************
  160 CONTINUE
C
C     DETERMINE THE LATERAL FORCE IN THE DYNAMIC SITUATION
C
C
C  FIND WHEEL INDEX J
C
      J=ID-200
C
      IF(ILAT(J).NE.1.OR.(ILAT(J).EQ.1.AND.IYFRIC(J).EQ.1)) THEN
C
C CALCULATE THE VERTICAL FORCE
C
C FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
       X1(J)=DATA(1)
       Y1(J)=DATA(3)
      ENDIF
C
C FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
       VF=0.0D0
       FMAG=0.0D0
       RETURN
      ENDIF
C
C FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
       CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
       VFORC=FZ(J)
       CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C  IF SO, SET LATERAL FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
      ENDIF
C
C FOR RIDE ANALYSIS, A LATERAL SPRING AND DAMPER IS APPRPORIATE
C
      IF(ILAT(J).EQ.1) THEN
C
C FIND THE LATERAL POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      YSPACE=DATA(2)
C
C FIND THE LATERAL VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      YDOT=DATA(2)
C
C FIND THE LATERAL FORCE
C
      YF=-AKL(J)*(YSPACE-YSE(J))-CLAT(J)*YDOT
C
C CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
      IF(IYFRIC(J).EQ.1) THEN
       YFS=YF
       REQMU=YF/VF
       REQMU=DABS(REQMU)
       IF(REQMU.GE.AMUL(J)) THEN
        YF=AMUL(J)*VF
        YF=DSIGN(YF,YFS)
       ENDIF
      ENDIF
C
      FMAG=YF
      RETURN
      ENDIF
C
C FOR HANDLING STUDIES, THE LATERAL FORCE IS BASED ON LATERAL FORCE-
C  SLIP ANGLE RELATIONS
C
      IF(ILAT(J).NE.1) THEN
C
C FIND THE LONGITUDINAL AND LATERAL VELOCITY COMPONENTS
C
      CALL INFO ( 'VEL', IWMKR ( J ), IGMKR, IVMKR ( J ), DATA, ERRFLG )
      UV = DATA ( 1 )
      VV = DATA ( 2 )
#ifdef FT_DEBUG
      write(88,
     + "('Time:',1p,e12.3,'   Local velocities:',3e12.3)")
     +     time, data(1), data(2), data(3)
#endif
C
C FIND THE STEERING ANGLE OF THE STEERED WHEELS
C
      STRANG = 0.0D0
      IF (ISTEER(J).EQ.2) THEN
      CALL INFO ( 'DISP', IWCMKR (J), IVMKR (J),IVMKR(J), DATA, ERRFLG)
      STRANG = DATAN2 ( DATA (2), DATA(1) )
      ENDIF
C
C FIND THE SLIP ANGLE OF THE WHEEL
C
      drvAng = DATAN2 ( VV, UV )
      SLPANG = ( DATAN2 ( VV, UV ) - STRANG ) * RTOD
      ABSSA = DABS ( SLPANG )
C
C FIND THE LATERAL FORCE
C
      IF(ILAT(J).EQ.2) SIDEF=AMUL(J)*(1.0-DEXP(-ALEXP(J)*ABSSA))*VF
      IF(ILAT(J).EQ.3) THEN
      CALL AKISPL(ABSSA,VF,ISPLIN(J),0,SPLVAL,ERRFLG)
      SIDEF=SPLVAL(1)
      ENDIF
C
C USE THE SIGN OF THE SLIP ANGLE TO FIND THE SIGN OF THE FORCE
C
      FMAG = -DSIGN ( SIDEF, SLPANG )

#ifdef FT_DEBUG
      write(88,
     + "('Time, Steer ang, drive ang, slip ang, fmag:'1p,5e12.3,//)")
     +     time, strang*RTOD, drvang*RTOD, slpang, fmag
#endif

      ENDIF
C
      RETURN
C
C *************
C
 230  IF ( ID .GT. 400 )    GO TO 250
C
C
C         V E R T I C A L      F O R C E S
C
C
C     DETERMINE THE VERTICAL FORCE ON THE WHEEL
C
C    IVERT(J)  = 1 - LINEAR SPRING MODEL FOR VERTICAL FORCE
C                2 - FORCE CALCULATED USING LOAD-DEFLECTION RELATION
C                3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C
      J = ID - 300
      IF(IFLAG)THEN
       VECTOR(1)=IWMKR(J)
       VECTOR(2)=IGMKR
       NV=2
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C CALCULATE THE VERTICAL FORCE
C
C FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
       X1(J)=DATA(1)
       Y1(J)=DATA(3)
      ENDIF
C
C FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
       VF=0.0D0
       FMAG=0.0D0
       RETURN
      ENDIF
C
C FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
       CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
       VFORC=FZ(J)
       CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
      IF ( VF .LT. 0.0D0 )   VF = 0.0D0
      FMAG = VF
      RETURN
C
C *************
C
  250 IF(ID.GT.500)GO TO 260
C
C      R O L L   M O M E N T
C
C          ROLL MOMENT OF LATERAL TIRE FORCES AND (OPTIONALLY)
C          OF VERTICAL FORCE DUE TO LATERAL TIRE DEFLECTION
C
C     INFORM ADAMS OF DEPENDENCIES ON MARKERS
C
      IF(IFLAG)THEN
       J=ID-400
       VECTOR(1)=IWMKR(J)
       VECTOR(2)=IGMKR
       VECTOR(3)=IMRKR(J)
       NV=3
        IF(ILAT(J).NE.1.AND.ISTEER(J).NE.2) THEN
         VECTOR(4)=IVMKR(J)
         NV=4
        ENDIF
        IF(ILAT(J).NE.1.AND.ISTEER(J).EQ.2) THEN
         VECTOR(4)=IVMKR(J)
         VECTOR(5)=IWCMKR(J)
         NV=5
        ENDIF
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C     CALCULATIONS FOR STATIC EQUILIBRIUM
C
C  THIS SECTION CALCULATES A LATERAL FORCE ON THE TIRE TO
C   PUT THE TRACTOR IN STATIC EQUILIBRIUM
C
      IF ( TIME .GT. 0.0D0 ) GO TO 760
      IF(PAR(1).NE.1) GO TO 750
       J=ID-400
C
C ASSUME LINEAR SPRING IN LATERAL DIRECTION
C  (FMAG=KY*(DATA(2)-DESIRED Y POSITION) WHERE KY=PAR(2) AND DESIRED Y
C   POSITION IS PAR(3))
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG =  - PAR(2) * ( DATA ( 2 ) - PAR ( 3 ) )
      FMAG=FMAG*YARM(J)
  750 RETURN
  760 CONTINUE
C
C  FIND WHEEL INDEX J
C
      J=ID-400
C
      IF(ILAT(J).NE.1.OR.(ILAT(J).EQ.1.AND.IYFRIC(J).EQ.1)) THEN
C
C CALCULATE THE VERTICAL FORCE
C
C FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
       X1(J)=DATA(1)
       Y1(J)=DATA(3)
      ENDIF
C
C FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
       VF=0.0D0
       FMAG=0.0D0
       RETURN
      ENDIF
C
C FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
       CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
       VFORC=FZ(J)
       CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C  IF SO, SET MOMENT TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
      ENDIF
C
C FOR RIDE ANALYSIS, A LATERAL SPRING AND DAMPER IS APPROPRIATE
C
      IF(ILAT(J).EQ.1) THEN
C
C FIND THE LATERAL POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      YSPACE=DATA(2)
C
C FIND THE LATERAL VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      YDOT=DATA(2)
C
C FIND THE LATERAL FORCE
C
      YF=-AKL(J)*(YSPACE-YSE(J))-CLAT(J)*YDOT
C
C CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
      IF(IYFRIC(J).EQ.1) THEN
       YFS=YF
       REQMU=YF/VF
       REQMU=DABS(REQMU)
       IF(REQMU.GE.AMUL(J)) THEN
        YF=AMUL(J)*VF
        YF=DSIGN(YF,YFS)
       ENDIF
      ENDIF
      ENDIF
C
C FOR HANDLING STUDIES, THE LATERAL FORCE IS BASED ON LATERAL FORCE-
C  SLIP ANGLE RELATIONS
C
      IF(ILAT(J).NE.1) THEN
C
C FIND THE LONGITUDINAL AND LATERAL VELOCITY COMPONENTS
C
      CALL INFO ( 'VEL', IWMKR ( J ), IGMKR, IVMKR ( J ), DATA, ERRFLG )
      UV = DATA ( 1 )
      VV = DATA ( 2 )
C
C FIND THE STEERING ANGLE OF THE STEERED WHEELS
C
      STRANG = 0.0D0
      IF (ISTEER(J).EQ.2) THEN
      CALL INFO ( 'DISP', IWCMKR (J), IVMKR (J),IVMKR(J), DATA, ERRFLG)
      STRANG = DATAN2 ( DATA (2), DATA(1) )
      ENDIF
C
C FIND THE SLIP ANGLE OF THE WHEEL
C
      SLPANG = ( DATAN2 ( VV, dabs(UV) ) - STRANG ) * RTOD
      ABSSA = DABS ( SLPANG )
      yFricSlipAngle(J) = slpAng
C
C FIND THE LATERAL FORCE
C
      IF(ILAT(J).EQ.2) SIDEF=AMUL(J)*(1.0-DEXP(-ALEXP(J)*ABSSA))*VF
      IF(ILAT(J).EQ.3) THEN
      CALL AKISPL(ABSSA,VF,ISPLIN(J),0,SPLVAL,ERRFLG)
      SIDEF=SPLVAL(1)
      ENDIF
C
C USE THE SIGN OF THE SLIP ANGLE TO FIND THE SIGN OF THE FORCE
C
      YF = -DSIGN ( SIDEF, SLPANG )
      ENDIF
C
C  DETERMINE THE ROLL MOMENT OF THE LATERAL FORCE
C
      FMAG=YF*YARM(J)
C
C  IF THE TIRE LATERAL DEFLECTION EFFECT IS TO BE CONSIDERED,
C   PERFORM THE ADDITIONAL COMPUTATIONS
C
      IF(IROLLM(J).EQ.1) THEN
C
C  DETERMINE THE LATERAL DEFLECTION
C
      FDEFL=YF/AKLAT(J)
C
C  CALCULATE THE MOMENT AND ADD IT TO THE MOMENT OF THE LATERAL FORCE
C
      FMAG=FMAG+VF*FDEFL
      ENDIF
      RETURN
C
  260 IF(ID.GT.600)GO TO 270
C
C            P I T C H   M O M E N T
C
C    IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C                1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C                2 - POWERED WHEEL - CONSIDER TRACTION & ROLLING RESISTANCE
C    IVERT(J)  = 3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C    ILONG(J)  = 0 - NO FORE-AFT SPRING AND DAMPER
C                1 - USE FORE-AFT SPRING AND DAMPER
C
C  THIS SECTION DEFINES THE ADAMS MARKERS USED TO DETERMINE THE FORCE
C
      IF(IFLAG) THEN
       J=ID-500
       VECTOR(1)=IWMKR(J)
       VECTOR(2)=IGMKR
       VECTOR(3)=IMRKR(J)
       NV=3
        IF(IPOWER(J).EQ.2.AND.IROTAT(J).EQ.1) THEN
         VECTOR(4)=IVMKR(J)
         NV=4
        ENDIF
        IF(IPOWER(J).EQ.2.AND.IROTAT(J).EQ.2) THEN
         VECTOR(4)=IVMKR(J)
         VECTOR(5)=IWROT(J)
         VECTOR(6)=JWROT(J)
         NV=6
        ENDIF
        IF(ILONG(J).EQ.1) THEN
         VECTOR(4)=IRAMKR
         NV=4
        ENDIF
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 870
C
C CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 860
      J=ID-500
C
C ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C  POSITION) WHERE KX=PAR(2), DESIRED X POSITION IS PAR(3) AND THE
C  MOMENT ARM OF THE FORCE IS PAR(4)
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3))*PAR(4)
  860 RETURN
C
C ***************
C
  870 CONTINUE
C
C      DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = ID - 500
C
C IF THE OPTIONS BEING USED DON'T REQUIRE CONSIDERATION OF THE MOMENT
C  OF THE FORE-AFT FORCE, SET THE MOMENT TO ZERO
C
      IF(IPOWER(J).NE.2.AND.ILONG(J).EQ.0) THEN
       FMAG=0.0
       RETURN
      ENDIF
C
C CALCULATE THE VERTICAL FORCE
C
C FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
       X1(J)=DATA(1)
       Y1(J)=DATA(3)
      ENDIF
C
C FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
       VF=0.0D0
       FMAG=0.0D0
       RETURN
      ENDIF
C
C FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
       CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
       VFORC=FZ(J)
       CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C  IF SO, SET MOMENT OF FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C IF IT IS TO BE CONSIDERED, FIND THE TRACTIVE FORCE
C
      IF(IPOWER(J).EQ.2) THEN
C
C FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
      CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
      UW = DATA (1)
C
C FOR POWERED WHEELS, FIND THE SLIPPAGE
C
C FIRST FIND THE ANGULAR VELOCITY OF THE WHEEL
C
      IF(IROTAT(J).EQ.1) THEN
C
C ESTIMATE THE ANGULAR VELOCITY OF THE DRIVE WHEEL USING A SIMPLE
C MODEL OF THE DIFFERENTIAL.
C
C FIND THE YAW VELOCITY OF THE CHASSIS
C
       CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
       WY = DATA(5)
       DPHI = DATA (6)
C
C ESTIMATE THE ANGULAR VELOCITY OF THE WHEEL USING A KINEMATIC MODEL
C  OF THE DIFFERENTIAL
C
      OMEGA =ANGVEL(J)-(DSIGN(1.0D0,ALAT(J))*TRDW(J)*DPHI)/(2.0D0*RR(J))
      ENDIF
C
      IF(IROTAT(J).EQ.2) THEN
C
C FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
       CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
       OMEGA=DATA(6)
#ifdef FT_DEBUG
       write(88,"('2 rot vel (3 is wheel speed):',1p,3e12.3)")
     +            data(4),data(5),data(6)
#endif
C
C FIND THE CHASSIS PITCH VELOCITY FROM THE MODEL
C
       CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
       WY = DATA(5)
      ENDIF
C
C FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = RR ( J ) * (OMEGA + WY)
C     VELWC = RR ( J ) * OMEGA
      ABSVWC = DABS ( VELWC )
      ABSUW = DABS ( UW )
C
C FIND THE LONGITUDINAL SLIPPAGE
C
      SLIP = ( ABSVWC - UW ) / DMAX1 ( ABSVWC, ABSUW )
      ASLIP = DABS ( SLIP )
C
C CALCULATE THE MOBILITY NUMBER
C
       IF(ISURF(J).EQ.1) THEN
        CN = (CI(J)*B(J)*D(J))/VF
        BN1 = 1.0 + 5.0 * VDEFL/H(J)
        BN2 = 1.0 + 3.0 * B(J)/D(J)
        BN = CN * BN1/BN2
       ENDIF
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED
C
C FIND THE TRACTIVE COEFFICIENT
C
      IF(ISURF(J).EQ.1) THEN
        IF(ICONST(J).EQ.1) POWER1 = -0.1*BN
        IF(ICONST(J).EQ.2) POWER1 = -0.1*BN
       Q1 = 1.0 - DEXP(POWER1)
        IF(ICONST(J).EQ.1) POWER2 = -7.5*ASLIP
        IF(ICONST(J).EQ.2) POWER2 = -9.5*ASLIP
       Q2 = 1.0 - DEXP(POWER2)
        IF(ICONST(J).EQ.1) AMU = 0.88 * Q1 * Q2 + 0.04
        IF(ICONST(J).EQ.2) AMU = 0.88 * Q1 * Q2 + 0.0325
      ENDIF
      IF(ISURF(J).EQ.2) THEN
       CN=(CI(J)*B(J)*D(J))/VF
       POWER3 = -CN*ASLIP
C      IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.02
       IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3))
       IF(ICONST(J).EQ.1) AMU = DSIGN(AMU,SLIP) + 0.02
C      IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.01625
       IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3))
       IF(ICONST(J).EQ.2) AMU = DSIGN(AMU,SLIP) + 0.01625
      ENDIF
C
C FINALLY THE TRACTIVE FORCE
C
      TRAC = AMU * VF
C
C THEN THE TORQUE IS THE PRODUCT OF THE TRACTIVE FORCE
C  AND THE ROLLING RADIUS
C
      FMAG = - TRAC * RR(J)
C
      ENDIF
C
C ALLOW THE OPTION OF INCLUDING A FORE-AFT SPRING AND DAMPER
C
      IF(ILONG(J).EQ.1) THEN
C
C FIND THE FORE-AFT POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      XSPACE=DATA(1)
C
C FIND THE FORE-AFT POSITION OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
      CALL INFO('DISP',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
      XRA=DATA(1)
C
C FIND THE FORE-AFT VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
      CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      XDOT=DATA(1)
C
C FIND THE FORE-AFT VELOCITY OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
      CALL INFO('VEL',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
      XRADOT=DATA(1)
C
C FIND THE FORE-AFT FORCE
C
      XF=AKFA(J)*(XSPACE-(XRA+XSE(J)))+CFA(J)*(XDOT-XRADOT)
C
C CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
      IF(IXFRIC(J).EQ.1) THEN
       XFS=XF
       REQMU=XF/VF
       REQMU=DABS(REQMU)
       IF(REQMU.GE.AMUFA(J)) THEN
        XF=AMUFA(J)*VF
        XF=DSIGN(XF,XFS)
       ENDIF
      ENDIF
C
C ADD THE MOMENT OF THE FORE-AFT SPRING FORCE TO THE TRACTIVE FORCE
C  MOMENT (IF THE LATTER IS NONZERO)
C
      FMAG=FMAG+XF*XARM(J)
      ENDIF
      RETURN
C
C    THIS SECTION COULD BE USED TO CALCULATE A MOMENT ABOUT A VERTICAL
C     AXIS.  AN EXAMPLE OF SUCH A TORQUE WOULD BE AN ALIGNING TORQUE.
C     BECAUSE OF LIMITED DATA CONCERNING THE ALIGNING TORQUES ACTING
C     ON OFF-ROAD TIRES, SUCH CALCULATIONS HAVE NOT BEEN INCLUDED HERE.
C
C     Y A W   M O M E N T
C
  270 IF(ID.GT.700)GO TO 280
C
C    SET MOMENT TO ZERO IF ID IS GREATER THAN 600 BUT LESS THAN 700
C
      FMAG=0.0D0
      RETURN
C
C      E N G I N E   T O R Q U E
C
C
  280 IF(ID.GT.800)GO TO 290
C
      IF(IFLAG) THEN
       VECTOR(1)=PAR(1)
       VECTOR(2)=PAR(2)
       NV=2
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  THIS SECTION IS USED FOR STATIC EQUILIBRIUM
C
      IF ( TIME .GT. 0.0D0 )   GO TO 281
C
C ENGINE TORQUE IS ZERO FOR STATIC EQUILIBRIUM
C
      FMAG=0.0
C
  281 CONTINUE
C
C  FIND THE ENGINE SPEED
C
      IMARKR=PAR(1)
      JMARKR=PAR(2)
      CALL INFO('VEL',IMARKR,JMARKR,JMARKR,DATA,ERRFLG)
      ESPEED=DATA(6)*(60.0/(2.0*PI))
C
C  INTERPOLATE TO FIND THE ENGINE TORQUE
C
      CALL VBM(ESPEED,ETORQ,ENGSPD,ENGTOR,NTSPTS)
      ETORQ=ETA*ETORQ
      FMAG=ETORQ
      RETURN
C
C      STATIC EQUILIBRIUM TORQUES ON ENGINE AND DRIVE WHEELS
C
  290 IF(ID.GT.900) GO TO 300
C
      IF(IFLAG) THEN
       VECTOR(1)=PAR(1)
       VECTOR(2)=PAR(2)
       NV=2
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  THIS SECTION IS USED FOR STATIC EQUILIBRIUM ONLY
C
      IF ( TIME .GT. 0.0D0 )   GO TO 291
C
C ASSUME A TORSIONAL SPRING "CENTERS" THE INERTIA FOR STATIC EQUILIBRIUM
C
      IMARKR=PAR(1)
      JMARKR=PAR(2)
      CALL INFO('DISP',IMARKR,JMARKR,JMARKR,DATA,ERRFLG)
      FMAG = - PAR(3)*DATA(4)
C
C ***************
C
  291 CONTINUE
C
C THE CENTERING TORQUE IS ZERO OTHERWISE
C
      FMAG = 0.0
      RETURN
C
  300 IF ( ID .GT. 1000) GO TO 310
C
C            L O N G I T U D I N A L   F O R C E S-B R A K I N G
C
C    IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C                1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C                2 - POWERED WHEEL - CONSIDER TRACTION/BRAKING FORCE
C
C TWO TYPES OF FORE-AFT OR LONGITUDINAL FORCES ARE CONSIDERED IN THIS
C  SECTION.  THESE ARE:
C    1. ROLLING RESISTANCE FORCES ON UNPOWERED WHEELS
C    2. TRACTIVE/BRAKING FORCES ON POWERED/BRAKED WHEELS
C THESE MODELS ARE DERIVED FROM AN EARLIER DRAM MODEL FOR COMBINE BRAKING 036700
C  THIS SECTION DEFINES THE ADAMS MARKERS USED TO DETERMINE THE FORCE
C
      IF(IFLAG) THEN
       J=ID-900
       VECTOR(1)=IWMKR(J)
       VECTOR(2)=IGMKR
       VECTOR(3)=IMRKR(J)
       VECTOR(4)=IWROT(J)
       VECTOR(5)=JWROT(J)
       NV=5
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 970
C
C CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 960
      J=ID-900
C
C ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C  POSITION) WHERE KX=PAR(2) AND DESIRED X POSITION IS PAR(3)
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3))
  960 RETURN
C
C ***************
C
  970 CONTINUE
C
C      DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = ID - 900
C
C CALCULATE THE VERTICAL FORCE
C
C FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
       X1(J)=DATA(1)
       Y1(J)=DATA(3)
      ENDIF
C
C FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
       VF=0.0D0
       FMAG=0.0D0
       RETURN
      ENDIF
C
C FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
       CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
       VFORC=FZ(J)
C IF THE DISTRIBUTED CONTACT MODEL IS BEING USED, INCLUDE THE FORE-AFT
C  COMPONENT
       FMAG=-FX(J)
       CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C  IF SO, SET FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C  FIND THE ROLLING RESISTANCE OR TRACTIVE/BRAKING FORCE
C
C FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
      CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
      UW = DATA (1)
      ABSUW = DABS ( UW )
C
C FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
      CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
      OMEGA=DATA(6)
#ifdef FT_DEBUG
       write(88,"('3 rot vel (3 is wheel speed):',1p,3e12.3)")
     +            data(4),data(5),data(6)
#endif
C
C CALCULATE THE ROLLING RESISTANCE FORCE ON AN UNPOWERED WHEEL
C
      IF(IPOWER(J).EQ.1) THEN
C
C FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = PAR(6) * OMEGA
      ABSVWC = DABS ( VELWC )
      VT = ABSUW - ABSVWC
      VTMAG = DABS(VT)
C
C CALCULATE THE ROLLING RESISTANCE COEFFICIENT, RHO
C
      IF(VTMAG.LE.PAR(5)) THEN
       ARG=POT*(VT/PAR(5))
       RHO=PAR(4)*DSIN(ARG)
      ENDIF
      IF(VTMAG.GT.PAR(5)) RHO=DSIGN(PAR(4),VT)
      RHO=RHO*DSIGN(1.0D0,OMEGA)
C
C FINALLY CALCULATE THE ROLLING RESISTANCE FORCE
C
      TF = RHO * VF
      FMAG = TF
      RETURN
      ENDIF
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED OR BRAKED
C
      IF(IPOWER(J).EQ.2) THEN
C
C FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = RR ( J ) * OMEGA
      ABSVWC = DABS ( VELWC )
C
C FIND THE LONGITUDINAL SLIPPAGE
C
      RELVEL = ABSVWC - ABSUW
      SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW, PAR(2))
C     SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW )
C     SLIP = RELVEL / ABSUW
      ASLIP = DABS ( SLIP )
C
C FIND THE TRACTIVE COEFFICIENT
C
      CN=(CI(J)*B(J)*D(J))/VF
      POWER3 = -CN*ASLIP
      AMU = 1.02 * (1.0 - DEXP(POWER3))
      AMU = DSIGN(AMU,SLIP)
C
C FINALLY THE TRACTIVE FORCE
C
      TRAC = AMU * VF * DSIGN(1.0D0,OMEGA)
      FMAG = - TRAC
      RETURN
      ENDIF
C
  310 IF ( ID .GT. 1100) GO TO 320
C
C            P I T C H  M O M E N T- A L T E R N A T E   M O D E L
C
C    IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C                1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C                2 - POWERED WHEEL - CONSIDER TRACTION/BRAKING FORCE
C
C TWO TYPES OF FORE-AFT OR LONGITUDINAL FORCES ARE CONSIDERED IN THIS
C  SECTION.  THESE ARE:
C    1. ROLLING RESISTANCE FORCES ON UNPOWERED WHEELS
C    2. TRACTIVE/BRAKING FORCES ON POWERED/BRAKED WHEELS
C THESE MODELS ARE DERIVED FROM AN EARLIER DRAM MODEL
C  THIS SECTION DEFINES THE ADAMS MARKERS USED TO DETERMINE THE FORCE
C
      IF(IFLAG) THEN
       J=ID-1000
       VECTOR(1)=IWMKR(J)
       VECTOR(2)=IGMKR
       VECTOR(3)=IMRKR(J)
       VECTOR(4)=IWROT(J)
       VECTOR(5)=JWROT(J)
       NV=5
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 1070
C
C CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 1060
      J=ID-1000
C
C ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C  POSITION) WHERE KX=PAR(2) AND DESIRED X POSITION IS PAR(3).
C PAR(4) = MOMENT ARM OF FORCE
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3)) * PAR(4)
 1060 RETURN
C
C ***************
C
 1070 CONTINUE
C
C      DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = ID - 1000
C
C CALCULATE THE VERTICAL FORCE
C
C FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
       X1(J)=DATA(1)
       Y1(J)=DATA(3)
      ENDIF
C
C FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
       VF=0.0D0
       FMAG=0.0D0
       RETURN
      ENDIF
C
C FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
       CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
       VFORC=FZ(J)
C IF THE DISTRIBUTED CONTACT MODEL IS BEING USED, INCLUDE THE FORE-AFT
C  COMPONENT
       FMAG=-FX(J)
       CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C  IF SO, SET FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C  FIND THE ROLLING RESISTANCE OR TRACTIVE/BRAKING FORCE
C
C FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
      CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
      UW = DATA (1)
      ABSUW = DABS ( UW )
C
C FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
      CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
      OMEGA=DATA(6)
#ifdef FT_DEBUG
       write(88,"('4 rot vel (3 is wheel speed):',1p,3e12.3)")
     +            data(4),data(5),data(6)
#endif
C
C CALCULATE THE ROLLING RESISTANCE FORCE ON AN UNPOWERED WHEEL
C
      IF(IPOWER(J).EQ.1) THEN
C
C FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = PAR(4) * OMEGA
      ABSVWC = DABS ( VELWC )
      VT = ABSUW - ABSVWC
      VTMAG = DABS(VT)
C
C CALCULATE THE ROLLING RESISTANCE COEFFICIENT, RHO
C
      IF(VTMAG.LE.PAR(6)) THEN
       ARG=POT*(VT/PAR(6))
       RHO=PAR(5)*DSIN(ARG)
      ENDIF
      IF(VTMAG.GT.PAR(6)) RHO=DSIGN(PAR(5),VT)
      RHO=RHO*DSIGN(1.0D0,OMEGA)
C
C FINALLY CALCULATE THE ROLLING RESISTANCE FORCE AND ITS MOMENT
C
      TF = RHO * VF
      FMAG = TF * PAR(4)
      RETURN
      ENDIF
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED OR BRAKED
C
      IF(IPOWER(J).EQ.2) THEN
C
C FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = RR ( J ) * OMEGA
      ABSVWC = DABS ( VELWC )
C
C FIND THE LONGITUDINAL SLIPPAGE
C
      RELVEL = ABSVWC - ABSUW
      SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW, PAR(2))
C     SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW )
C     SLIP = RELVEL / ABSUW
      ASLIP = DABS ( SLIP )
C
C FIND THE TRACTIVE COEFFICIENT
C
      CN=(CI(J)*B(J)*D(J))/VF
      POWER3 = -CN*ASLIP
      AMU = 1.02 * (1.0 - DEXP(POWER3))
      AMU = DSIGN(AMU,SLIP)
C
C FINALLY THE TRACTIVE FORCE AND ITS MOMENT
C
      TRAC = AMU * VF * DSIGN(1.0D0,OMEGA)
      FMAG = - TRAC * RR(J)
      RETURN
      ENDIF
C
C
  320 IF ( ID .GT. 1200) GO TO 330
C
C            B R A K I N G   M O M E N T
C
C  THIS SECTION DEFINES THE ADAMS MARKERS USED TO DETERMINE THE MOMENT
C
      IF(IFLAG) THEN
       J=ID-1100
       VECTOR(1)=IWROT(J)
       VECTOR(2)=JWROT(J)
       NV=2
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C      DETERMINE THE BRAKING MOMENT FOR THE DYNAMIC CONDITION
C
C FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = ID - 1100
C
C FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
      CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
      OMEGA=DATA(6)
      ABSAV=DABS(OMEGA)
#ifdef FT_DEBUG
       write(88,"('5 rot vel (3 is wheel speed):',1p,3e12.3)")
     +            data(4),data(5),data(6)
#endif
C
C FIND THE MAXIMUM TORQUE AS A FUNCTION OF TIME
C
      OMEGA1=PI/PAR(3)
      TBEGIN=PAR(4)
      TEND=PAR(4)+PAR(3)
      IF(TIME.LE.TBEGIN) TMAX=0.0D0
      IF(TIME.GT.TBEGIN.AND.TIME.LE.TEND) THEN
        TSTAR=TIME-TBEGIN
        TMAX=(PAR(1)/2.0D0)*(1.0D0-DCOS(OMEGA1*TSTAR))
      ENDIF
      IF(TIME.GT.TEND) TMAX=PAR(1)
C
C FIND THE BRAKING TORQUE AS A FUNCTION OF THE RELATIVE ANGULAR VELOCITY
C
      IF(ABSAV.LE.PAR(2)) THEN
       OMEGA2=POT*(OMEGA/PAR(2))
       TBRAKE=-TMAX*DSIN(OMEGA2)
      ENDIF
      IF(ABSAV.GT.PAR(2)) TBRAKE=-DSIGN(TMAX,OMEGA)
C
C FINALLY CALCULATE THE MOMENT
C
      FMAG = TBRAKE
      RETURN
C
C   NO 1200 SFORCES CURRENTLY EXIST
C
  330 IF(ID.GT.1300) GO TO 340
C
C    SET SFORCE TO ZERO IF ID IS GREATER THAN 1200 BUT LESS THAN 1300
C
      FMAG=0.0D0
      RETURN
C
C   S H O C K   A B S O R B E R S
C
  340 IF(ID.GT.1400) GO TO 240
C
      J = ID - 1300
      IF(IFLAG)THEN
       VECTOR(1)=ISHCK(J)
       VECTOR(2)=JSHCK(J)
       VECTOR(3)=IGMKR
       NV=3
       CALL FNCDEP('SFOSUB','MARKER',ID,VECTOR,NV,ERRFLG)
      ENDIF
C
C  FIRST FIND THE RELATIVE VELOCITY IN GLOBAL COORDINATES
C
      CALL INFO('VEL',ISHCK(J),JSHCK(J),IGMKR,DATA,ERRFLG)
      XVEL=DATA(1)
      YVEL=DATA(2)
      ZVEL=DATA(3)
C
C  NEXT FIND THE RELATIVE DISPLACEMENT IN GLOBAL COORDINATES
C
      CALL INFO('DISP',ISHCK(J),JSHCK(J),IGMKR,DATA,ERRFLG)
      XDISP=DATA(1)
      YDISP=DATA(2)
      ZDISP=DATA(3)
      DIST=DSQRT(XDISP**2+YDISP**2+ZDISP**2)
C
C  TAKE THE DOT PRODUCT TO CALCULATE THE RELATIVE VELOCITY
C
      RVEL=(XVEL*XDISP+YVEL*YDISP+ZVEL*ZDISP)/DIST
      ARVEL=DABS(RVEL)
C
      IF(RVEL.LE.0.0D0)GO TO 6107
C
C  SHOCK ABSORBER IS IN EXTENSION
C
      CALL VBM2(J,ARVEL,DFORCE,VELEX,FOREX,NEPTS)
      GO TO 6108
C
C  SHOCK ABSORBER IS IN COMPRESSION
C
 6107 CALL VBM2(J,ARVEL,DFORCE,VELCP,FORCP,NCPTS)
C
 6108 IF(ARVEL.EQ.0.0D0)GO TO 6109
      FMAG=-(RVEL/ARVEL)*DFORCE
      RETURN
 6109 FMAG=0.0D0
      RETURN
C **************
C
C
  240 CONTINUE
C
C  LONGITUDINAL OR LATERAL FORCE IS ZERO IF VERTICAL FORCE IS ZERO
C
      FMAG = 0.0D0
      RETURN
      END
C
C****************************************************************
C
      SUBROUTINE REQSUB (IDR, TIME, PAR, NPAR, IFLAG, RESULT )
C
C     USER REQUEST SUBROUTINE
C
      use jdt_inp
      use jdt_var

      IMPLICIT REAL*8 (A-H,O-Z)

      CHARACTER*30 AFILE
      DIMENSION PAR (*), RESULT(8), DATA(6), SPLVAL(3)
      DIMENSION XDIST(10),YDIST(10),ZDIST(10),XFORC(10),YFORC(10)
      DIMENSION ZFORC(10),XMOM(10),YMOM(10),ZMOM(10)
      DIMENSION TIMEX(1000)
C     DIMENSION XDBAR(1000),YDBAR(1000)
      DIMENSION HITCHX(1000),HITCHY(1000),HITCHZ(1000)
      DIMENSION ROX(1000),ROY(1000),ROZ(1000)
      DIMENSION RIX(1000),RIY(1000),RIZ(1000)
      DIMENSION ALIX(1000),ALIY(1000),ALIZ(1000)
      DIMENSION ALOX(1000),ALOY(1000),ALOZ(1000)
      DIMENSION ACC(2200)
      LOGICAL ERRFLG, IFLAG
C
      DATA IFPSS/0/
      DATA IFPASS/0/
      DATA JFPASS/0/
      DATA JCOUNT/1/
      DATA ICOUNT/1/
      DATA PI/3.1415926535897932D0/
      DATA POT/1.5707963267948D0/
      DATA RTOD /57.29578D0/
      DATA DTR/0.01745329/

      real*8 zeroDelta(3)

      zeroDelta(1) = 0.0d0
      zeroDelta(2) = 0.0d0
      zeroDelta(3) = 0.0d0
C
C
C     INITIALIZE THE FORCE MAGNITUDE TO ZERO
C
      FMAG = 0.0D0
C
      IF ( IDR.GT. 200 )       GO TO 140
C
C     L O N G I T U D I N A L      F O R C E S
C
C     IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C     1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C     2 - POWERED WHEEL - CONSIDER TRACTION & ROLLING RESISTANCE
C     IVERT(J)  = 3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C     ILONG(J)  = 0 - NO FORE-AFT SPRING AND DAMPER
C     1 - USE FORE-AFT SPRING AND DAMPER
C
C     THREE TYPES OF FORE-AFT OR LONGITUDINAL FORCES ARE CONSIDERED IN THIS
C     SECTION.  THESE ARE:
C     1. ROLLING RESISTANCE AND TRACTIVE FORCES
C     2. THE HORIZONTAL FORCE COMPONENT PREDICTED BY THE DISTRIBUTED
C     CONTACT MODEL
C     3. A LINEAR SPRING AND DAMPER IN THE FORE-AFT DIRECTION
C
C     IF THE OPTIONS BEING USED DON'T INCLUDE A FORE-AFT FORCE, RETURN
C
C     FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = IDR- 100
      IF(IPOWER(J).EQ.0.AND.IVERT(J).NE.3.AND.ILONG(J).EQ.0) THEN
         FMAG=0.0
         RETURN
      ENDIF
C
C     THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 70
C
C     CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 60
      J=IDR-100
C
C     ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C     POSITION) WHERE KX=PAR(2) AND DESIRED X POSITION IS PAR(3)
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3))
      RESULT(2)=DATA(1)
      RESULT(4)=DATA(1)-PAR(3)
      RESULT(8)=FMAG
 60   RETURN
C
C     ***************
C
 70   CONTINUE
C
C     DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C     FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = IDR- 100
C
C     CALCULATE THE VERTICAL FORCE
C
C     FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
         X1(J)=DATA(1)
         Y1(J)=DATA(3)
      ENDIF
C
C     FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C     FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
         VF=0.0D0
         FMAG=0.0D0
         RETURN
      ENDIF
C
C     FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
         CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
         VFORC=FZ(J)
C     IF THE DISTRIBUTED CONTACT MODEL IS BEING USED, INCLUDE THE FORE-AFT
C     COMPONENT
         FMAG=-FX(J)
         CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
         RESULT(4)=FMAG
      ENDIF
C
C     FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C     FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C     FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C     CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C     IF SO, SET FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C     IF THEY ARE TO BE CONSIDERED, FIND THE TRACTIVE AND ROLLING RESISTANCE
C     FORCES
C
      IF(IPOWER(J).NE.0) THEN
C
C     FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
         CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
         UW = DATA (1)
C
C     FOR POWERED WHEELS, FIND THE SLIPPAGE
C
         IF(IPOWER(J).EQ.2) THEN
C
C     FIRST FIND THE ANGULAR VELOCITY OF THE WHEEL
C
            IF(IROTAT(J).EQ.1) THEN
C
C     ESTIMATE THE ANGULAR VELOCITY OF THE DRIVE WHEEL USING A SIMPLE
C     MODEL OF THE DIFFERENTIAL.
C
C     FIND THE YAW VELOCITY OF THE CHASSIS
C
               CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
               WY = DATA(5)
               DPHI = DATA (6)
C
C     ESTIMATE THE ANGULAR VELOCITY OF THE WHEEL USING A KINEMATIC MODEL
C     OF THE DIFFERENTIAL
C
               OMEGA =ANGVEL(J)-(DSIGN(1.0D0,ALAT(J))*TRDW(J)*DPHI)
     +                /(2.0D0*RR(J))
            ENDIF
C
            IF(IROTAT(J).EQ.2) THEN
C
C     FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
               CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
               OMEGA=DATA(6)

#ifdef FT_DEBUG
               write(88,"('6 rot vel (3 is wheel speed):',1p,3e12.3)")
     +              data(4),data(5),data(6)
#endif
C
C     FIND THE CHASSIS PITCH VELOCITY FROM THE MODEL
C
               CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
               WY = DATA(5)
            ENDIF
C
C     FIND THE NO SLIP FORWARD VELOCITY
C
            VELWC = RR ( J ) * (OMEGA + WY)
C     VELWC = RR ( J ) * OMEGA
            ABSVWC = DABS ( VELWC )
            ABSUW = DABS ( UW )
C
C     FIND THE LONGITUDINAL SLIPPAGE
C
            SLIP = ( ABSVWC - UW ) / DMAX1 ( ABSVWC, ABSUW )
            RESULT(2)=SLIP
            ASLIP = DABS ( SLIP )
         ENDIF
C
C     CALCULATE THE ROLLING RESISTANCE COEFFICIENT, RHO
C
         IF(ISURF(J).EQ.1) THEN
            CN = (CI(J)*B(J)*D(J))/VF
            BN1 = 1.0 + 5.0 * VDEFL/H(J)
            BN2 = 1.0 + 3.0 * B(J)/D(J)
            BN = CN * BN1/BN2
            IF(IPOWER(J).EQ.1.AND.ICONST(J).EQ.1) THEN
               RHO=1.0/BN + 0.04
            ENDIF
            IF(IPOWER(J).EQ.1.AND.ICONST(J).EQ.2) THEN
               RHO=0.9/BN + 0.0325
            ENDIF
            IF(IPOWER(J).EQ.2.AND.ICONST(J).EQ.1) THEN
               RHO=1.0/BN + 0.04 + 0.5*ASLIP/DSQRT(BN)
            ENDIF
            IF(IPOWER(J).EQ.2.AND.ICONST(J).EQ.2) THEN
               RHO=0.9/BN + 0.0325 + 0.5*ASLIP/DSQRT(BN)
            ENDIF
         ENDIF
         IF(ISURF(J).EQ.2.AND.ICONST(J).EQ.1) RHO=0.02
         IF(ISURF(J).EQ.2.AND.ICONST(J).EQ.2) RHO=0.01625
         IF(IPOWER(J).EQ.1)RESULT(3)=RHO
C
C     FINALLY CALCULATE THE ROLLING RESISTANCE FORCE
C
         TF = RHO * VF
         RESULT(5)=VF
         RESULT(6)=DSIGN(TF,UW)
C
C     THE ROLLING RESISTANCE IS ASSUMED TO ACT IN THE DIRECTION OPPOSITE TO
C     THE DIRECTION OF THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL.
C
         FMAG = FMAG + DSIGN ( TF, UW )
         RESULT(7)=FMAG
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED
C
         IF(IPOWER(J).EQ.2) THEN
C
C     FIND THE TRACTIVE COEFFICIENT
C
         IF(ISURF(J).EQ.1) THEN
            IF(ICONST(J).EQ.1) POWER1 = -0.1*BN
            IF(ICONST(J).EQ.2) POWER1 = -0.1*BN
            Q1 = 1.0 - DEXP(POWER1)
            IF(ICONST(J).EQ.1) POWER2 = -7.5*ASLIP
            IF(ICONST(J).EQ.2) POWER2 = -9.5*ASLIP
            Q2 = 1.0 - DEXP(POWER2)
            IF(ICONST(J).EQ.1) AMU = 0.88 * Q1 * Q2 + 0.04
            IF(ICONST(J).EQ.2) AMU = 0.88 * Q1 * Q2 + 0.0325
         ENDIF
         IF(ISURF(J).EQ.2) THEN
            CN=(CI(J)*B(J)*D(J))/VF
            RESULT(1)=CN
            POWER3 = -CN*ASLIP
C     IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.02
            IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3))
            IF(ICONST(J).EQ.1) AMU = DSIGN(AMU,SLIP) + 0.02
C     IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.01625
            IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3))
            IF(ICONST(J).EQ.2) AMU = DSIGN(AMU,SLIP) + 0.01625
         ENDIF
C
C     FINALLY THE TRACTIVE FORCE
C
         TRAC = AMU * VF
         IF(ISURF(J).EQ.1)RESULT(1)=BN
         IF(ISURF(J).EQ.2)RESULT(1)=CN
         RESULT(3)=AMU
         RESULT(5)=VF
         RESULT(7)=TRAC
C
C     COMBINE THE ROLLING RESISTANCE WITH THE TRACTIVE FORCE
C
         FMAG = FMAG - TRAC
         RESULT(8)=FMAG
      ENDIF
C
      ENDIF
C
C     ALLOW THE OPTION OF INCLUDING A FORE-AFT SPRING AND DAMPER
C
      IF(ILONG(J).EQ.1) THEN
C
C     FIND THE FORE-AFT POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         XSPACE=DATA(1)
C
C     FIND THE FORE-AFT POSITION OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
         CALL INFO('DISP',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
         XRA=DATA(1)
C
C     FIND THE FORE-AFT VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         XDOT=DATA(1)
C
C     FIND THE FORE-AFT VELOCITY OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
         CALL INFO('VEL',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
         XRADOT=DATA(1)
C
C     FIND THE FORE-AFT FORCE
C
         XF=AKFA(J)*(XSPACE-(XRA+XSE(J)))+CFA(J)*(XDOT-XRADOT)
C
C     CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
         IF(IXFRIC(J).EQ.1) THEN
            XFS=XF
            REQMU=XF/VF
            REQMU=DABS(REQMU)
            IF(REQMU.GE.AMUFA(J)) THEN
               XF=AMUFA(J)*VF
               XF=DSIGN(XF,XFS)
            ENDIF
         ENDIF
C
C     FORM THE RESULTS ARRAY
C
         RESULT(1)=AMUFA(J)
         RESULT(2)=XSPACE
         RESULT(3)=XDOT
         RESULT(4)=XSPACE-(XRA+XSE(J))
         RESULT(5)=REQMU
         RESULT(6)=XDOT-XRADOT
         RESULT(7)=VF
         RESULT(8)=XF
C
C     ADD THE FORE-AFT SPRING-DAMPER FORCE TO FORM THE TOTAL FORE-AFT FORCE
C
         FMAG=FMAG+XF
      ENDIF
      RETURN
C
C     **************
C
 140  IF ( IDR.GT. 300 )  GO TO 230
C
C
C     L A T E R A L       F O R C E S
C
C     ISTEER(J) = 0 - NO STEERING INPUT (RIDE ANALYSIS)
C     1 - USER INPUT - WHEEL DOESN'T PHYSICALLY STEER IN MODEL
C     2 - STEER ANGLE DETERMINED FROM MODEL
C     ILAT(J)   = 1 - LATERAL SPRING AND DAMPER (RIDE ANALYSIS)
C     2 - EQUATION FOR LATERAL FORCE-SLIP ANGLE RELATION
C     3 - CARPET PLOT DATA FOR LATERAL FORCE AS FUNCTION OF
C     SLIP ANGLE AND LATERAL FORCE
C
C     DETERMINE THE LATERAL FORCE ON THE WHEEL
C
C     THIS SECTION CALCULATES A LATERAL FORCE ON THE TIRE TO
C     PUT THE TRACTOR IN STATIC EQUILIBRIUM
C
      IF ( TIME .GT. 0.0D0 ) GO TO 160
      IF(PAR(1).NE.1) GO TO 150
      J=IDR-200
C
C     ASSUME LINEAR SPRING IN LATERAL DIRECTION
C     (FMAG=KY*(DATA(2)-DESIRED Y POSITION) WHERE KY=PAR(2) AND DESIRED Y
C     POSITION IS PAR(3))
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG =  - PAR(2) * ( DATA ( 2 ) - PAR ( 3 ) )
      RESULT(2)=DATA(2)
      RESULT(4)=DATA(2)-PAR(3)
      RESULT(7)=FMAG
 150  RETURN
C
C     *************
 160  CONTINUE
C
C     DETERMINE THE LATERAL FORCE IN THE DYNAMIC SITUATION
C
C
C     FIND WHEEL INDEX J
C
      J=IDR-200
C
C     CALCULATE THE VERTICAL FORCE
C
      IF(ILAT(J).NE.1.OR.(ILAT(J).EQ.1.AND.IYFRIC(J).EQ.1)) THEN
C
C     FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
         CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         AXLHGT=DATA(3)
         IF(IVERT(J).EQ.3) THEN
            X1(J)=DATA(1)
            Y1(J)=DATA(3)
         ENDIF
C
C     FIND THE TIRE DEFLECTION
C
         VDISP=AXLHGT-WRAD(J)
C
C     FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
         IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
            VF=0.0D0
            FMAG=0.0D0
            RETURN
         ENDIF
C
C     FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
         IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
         IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
         IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
         IF(IVERT(J).EQ.3) THEN
            CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
            VFORC=FZ(J)
            CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
         ENDIF
C
C     FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
         CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         VVEL=DATA(3)
C
C     FIND THE VERTICAL DAMPING COEFFICIENT
C
         CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C     FIND THE VERTICAL FORCE
C
         VF = VFORC - CDAMP * VVEL
C
C     CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C     IF SO, SET LATERAL FORCE TO ZERO
C
         IF(VF.LE.0.0D0) GO TO 240
      ENDIF
C
C     FOR RIDE ANALYSIS, A LATERAL SPRING AND DAMPER IS APPRPORIATE
C
      IF(ILAT(J).EQ.1) THEN
C
C     FIND THE LATERAL POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         YSPACE=DATA(2)
         RESULT(2)=YSPACE
C
C     FIND THE LATERAL VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         YDOT=DATA(2)
         RESULT(3)=YDOT
C
C     FIND THE LATERAL FORCE
C
         YF=-AKL(J)*(YSPACE-YSE(J))-CLAT(J)*YDOT
         RESULT(4)=(YSPACE-YSE(J))
C
C     CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
         IF(IYFRIC(J).EQ.1) THEN
            RESULT(1)=AMUL(J)
            RESULT(6)=VF
            YFS=YF
            REQMU=YF/VF
            REQMU=DABS(REQMU)
            RESULT(8)=REQMU
            IF(REQMU.GE.AMUL(J)) THEN
               YF=AMUL(J)*VF
               YF=DSIGN(YF,YFS)
            ENDIF
         ENDIF
C
         FMAG=YF
         RESULT(7)=YF
         RETURN
      ENDIF
C
C     FOR HANDLING STUDIES, THE LATERAL FORCE IS BASED ON LATERAL FORCE-
C     SLIP ANGLE RELATIONS
C
      IF(ILAT(J).NE.1) THEN
C
C     FIND THE LONGITUDINAL AND LATERAL VELOCITY COMPONENTS
C
         CALL INFO ( 'VEL', IWMKR(J), IGMKR, IVMKR(J), DATA, ERRFLG )
         UV = DATA ( 1 )
         VV = DATA ( 2 )
         RESULT(1)=UV
         RESULT(2)=VV
C
C     FIND THE STEERING ANGLE OF THE STEERED WHEELS
C
         STRANG = 0.0D0
         IF (ISTEER(J).EQ.2) THEN
         CALL INFO ('DISP', IWCMKR(J), IVMKR(J),IVMKR(J), DATA, ERRFLG)
         STRANG = DATAN2 ( DATA (2), DATA(1) )
         RESULT(3)=STRANG*RTOD
      ENDIF
C
C     FIND THE SLIP ANGLE OF THE WHEEL
C
      SLPANG = ( DATAN2 ( VV, UV ) - STRANG ) * RTOD
      RESULT(4)=SLPANG
      ABSSA = DABS ( SLPANG )
C
C     FIND THE LATERAL FORCE
C
      IF(ILAT(J).EQ.2) THEN
         SIDEF=AMUL(J)*(1.0-DEXP(-ALEXP(J)*ABSSA))*VF
         RESULT(8)=AMUL(J)*(1.0-DEXP(-ALEXP(J)*ABSSA))
      ENDIF
      IF(ILAT(J).EQ.3) THEN
         CALL AKISPL(ABSSA,VF,ISPLIN(J),0,SPLVAL,ERRFLG)
         SIDEF=SPLVAL(1)
         RESULT(8)=SIDEF/VF
      ENDIF
C
C     USE THE SIGN OF THE SLIP ANGLE TO FIND THE SIGN OF THE FORCE
C
      FMAG = -DSIGN ( SIDEF, SLPANG )
      RESULT(6)=VF
      RESULT(7)=FMAG
      ENDIF
C
      RETURN
C
C     *************
C
 230  IF ( IDR.GT. 400 )    GO TO 250
C
C
C     V E R T I C A L      F O R C E S
C
C
C     DETERMINE THE VERTICAL FORCE ON THE WHEEL
C
C     IVERT(J)  = 1 - LINEAR SPRING MODEL FOR VERTICAL FORCE
C     2 - FORCE CALCULATED USING LOAD-DEFLECTION RELATION
C     3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C
      J = IDR- 300
C
C     CALCULATE THE VERTICAL FORCE
C
C     FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      XWHL=DATA(1)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
         X1(J)=DATA(1)
         Y1(J)=DATA(3)
      ENDIF
C
C     FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C     FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
         VF=0.0D0
         FMAG=0.0D0
         RETURN
      ENDIF
C
C     FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
         CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
         VFORC=FZ(J)
         RESULT(8)=-FX(J)
         CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C     FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C     FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C     FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
      IF ( VF .LT. 0.0D0 )   VF = 0.0D0
      FMAG = VF
      RESULT(1)=VFORC
      RESULT(2)=XWHL
      RESULT(3)=AXLHGT
      RESULT(4)=VVEL
      RESULT(5)=CDAMP
      RESULT(6)=VDEFL
      RESULT(7)=VF
      RETURN
C
C     *************
C
 250  IF(IDR.GT.500)GO TO 260
C
C     R O L L   M O M E N T
C
C     ROLL MOMENT OF LATERAL TIRE FORCES AND (OPTIONALLY)
C     OF VERTICAL FORCE DUE TO LATERAL TIRE DEFLECTION
C
C     CALCULATIONS FOR STATIC EQUILIBRIUM
C
C     THIS SECTION CALCULATES A LATERAL FORCE ON THE TIRE TO
C     PUT THE TRACTOR IN STATIC EQUILIBRIUM
C
      IF ( TIME .GT. 0.0D0 ) GO TO 760
      IF(PAR(1).NE.1) GO TO 750
C
C     ASSUME LINEAR SPRING IN LATERAL DIRECTION
C     (FMAG=KY*(DATA(2)-DESIRED Y POSITION) WHERE KY=PAR(2) AND DESIRED Y
C     POSITION IS PAR(3))
C
      J=IDR-400
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG =  - PAR(2) * ( DATA ( 2 ) - PAR ( 3 ) )
      RESULT(3)=FMAG
      FMAG=FMAG*YARM(J)
      RESULT(4)=FMAG
 750  RETURN
C
C     FIND WHEEL INDEX J
C
 760  J=IDR-400
C
C     CALCULATE THE VERTICAL FORCE
C
      IF(ILAT(J).NE.1.OR.(ILAT(J).EQ.1.AND.IYFRIC(J).EQ.1)) THEN
C
C     FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
         CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         AXLHGT=DATA(3)
         IF(IVERT(J).EQ.3) THEN
            X1(J)=DATA(1)
            Y1(J)=DATA(3)
         ENDIF
C
C     FIND THE TIRE DEFLECTION
C
         VDISP=AXLHGT-WRAD(J)
C
C     FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
         IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
            VF=0.0D0
            FMAG=0.0D0
            RETURN
         ENDIF
C
C     FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
         IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
         IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
         IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
         IF(IVERT(J).EQ.3) THEN
            CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
            VFORC=FZ(J)
            CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
         ENDIF
C
C     FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
         CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         VVEL=DATA(3)
C
C     FIND THE VERTICAL DAMPING COEFFICIENT
C
         CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C     FIND THE VERTICAL FORCE
C
         VF = VFORC - CDAMP * VVEL
C
C     CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C     IF SO, SET MOMENT TO ZERO
C
         IF(VF.LE.0.0D0) GO TO 240
      ENDIF
C
C     FOR RIDE ANALYSIS, A LATERAL SPRING AND DAMPER IS APPROPRIATE
C
      IF(ILAT(J).EQ.1) THEN
C
C     FIND THE LATERAL POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         YSPACE=DATA(2)
C
C     FIND THE LATERAL VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         YDOT=DATA(2)
C
C     FIND THE LATERAL FORCE
C
         YF=-AKL(J)*(YSPACE-YSE(J))-CLAT(J)*YDOT
C
C     CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
         IF(IYFRIC(J).EQ.1) THEN
            YFS=YF
            REQMU=YF/VF
            REQMU=DABS(REQMU)
            IF(REQMU.GE.AMUL(J)) THEN
               YF=AMUL(J)*VF
               YF=DSIGN(YF,YFS)
            ENDIF
         ENDIF
      ENDIF
C
C     FOR HANDLING STUDIES, THE LATERAL FORCE IS BASED ON LATERAL FORCE-
C     SLIP ANGLE RELATIONS
C
      IF(ILAT(J).NE.1) THEN
C
C     FIND THE LONGITUDINAL AND LATERAL VELOCITY COMPONENTS
C
         CALL INFO ('VEL', IWMKR(J), IGMKR, IVMKR(J), DATA, ERRFLG )
         UV = DATA ( 1 )
         VV = DATA ( 2 )
C
C     FIND THE STEERING ANGLE OF THE STEERED WHEELS
C
         STRANG = 0.0D0
         IF (ISTEER(J).EQ.2) THEN
         CALL INFO ( 'DISP', IWCMKR(J), IVMKR(J),IVMKR(J), DATA, ERRFLG)
         STRANG = DATAN2 ( DATA (2), DATA(1) )
      ENDIF
C
C     FIND THE SLIP ANGLE OF THE WHEEL
C
      SLPANG = ( DATAN2 ( VV, UV ) - STRANG ) * RTOD
      ABSSA = DABS ( SLPANG )
C
C     FIND THE LATERAL FORCE
C
      IF(ILAT(J).EQ.2) SIDEF=AMUL(J)*(1.0-DEXP(-ALEXP(J)*ABSSA))*VF
      IF(ILAT(J).EQ.3) THEN
         CALL AKISPL(ABSSA,VF,ISPLIN(J),0,SPLVAL,ERRFLG)
         SIDEF=SPLVAL(1)
      ENDIF
C
C     USE THE SIGN OF THE SLIP ANGLE TO FIND THE SIGN OF THE FORCE
C
      YF = -DSIGN ( SIDEF, SLPANG )
      ENDIF
C
C     DETERMINE THE ROLL MOMENT OF THE LATERAL FORCE
C
      FMAG=YF*YARM(J)
      RESULT(2)=VF
      RESULT(3)=YF
      RESULT(4)=FMAG
C
C     IF THE TIRE LATERAL DEFLECTION EFFECT IS TO BE CONSIDERED,
C     PERFORM THE ADDITIONAL COMPUTATIONS
C
      IF(IROLLM(J).EQ.1) THEN
C
C     DETERMINE THE LATERAL DEFLECTION
C
         FDEFL=YF/AKLAT(J)
         RESULT(6)=FDEFL
C
C     CALCULATE THE MOMENT AND ADD IT TO THE MOMENT OF THE LATERAL FORCE
C
         FMAG=FMAG+VF*FDEFL
         RESULT(7)=VF*FDEFL
      ENDIF
      RETURN
C
 260  IF(IDR.GT.600)GO TO 270
C
C     P I T C H   M O M E N T
C
C     IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C     1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C     2 - POWERED WHEEL - CONSIDER TRACTION & ROLLING RESISTANCE
C     IVERT(J)  = 3 - DISTRIBUTED CONTACT MODEL FOR VERTICAL FORCE
C     ILONG(J)  = 0 - NO FORE-AFT SPRING AND DAMPER
C     1 - USE FORE-AFT SPRING AND DAMPER
C
C     THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 870
C
C     CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 860
      J=IDR-500
C
C     ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C     POSITION) WHERE KX=PAR(2), DESIRED X POSITION IS PAR(3) AND THE
C     MOMENT ARM OF THE FORCE IS PAR(4)
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3))
      RESULT(6)=FMAG
      FMAG=FMAG*PAR(4)
      RESULT(7)=FMAG
 860  RETURN
C
C     ***************
C
 870  CONTINUE
C
C     DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C     FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = IDR- 500
C
C     IF THE OPTIONS BEING USED DON'T REQUIRE CONSIDERATION OF THE MOMENT
C     OF THE FORE-AFT FORCE, SET THE MOMENT TO ZERO
C
      IF(IPOWER(J).NE.2.AND.ILONG(J).EQ.0) THEN
         FMAG=0.0
         RETURN
      ENDIF
C
C     CALCULATE THE VERTICAL FORCE
C
C     FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
         X1(J)=DATA(1)
         Y1(J)=DATA(3)
      ENDIF
C
C     FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C     FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
         VF=0.0D0
         FMAG=0.0D0
         RETURN
      ENDIF
C
C     FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
         CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
         VFORC=FZ(J)
         CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C     FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C     FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C     FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C     CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C     IF SO, SET MOMENT OF FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C     IF IT IS TO BE CONSIDERED, FIND THE TRACTIVE FORCE
C
      IF(IPOWER(J).EQ.2) THEN
C
C     FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
         CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
         UW = DATA (1)
C
C     FOR POWERED WHEELS, FIND THE SLIPPAGE
C
C     FIRST FIND THE ANGULAR VELOCITY OF THE WHEEL
C
         IF(IROTAT(J).EQ.1) THEN
C
C     ESTIMATE THE ANGULAR VELOCITY OF THE DRIVE WHEEL USING A SIMPLE
C     MODEL OF THE DIFFERENTIAL.
C
C     FIND THE YAW VELOCITY OF THE CHASSIS
C
            CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
            WY = DATA(5)
            DPHI = DATA (6)
C
C     ESTIMATE THE ANGULAR VELOCITY OF THE WHEEL USING A KINEMATIC MODEL
C     OF THE DIFFERENTIAL
C
            OMEGA =ANGVEL(J)-(DSIGN(1.0D0,ALAT(J))*TRDW(J)*DPHI)
     +             /(2.0D0*RR(J))
         ENDIF
C
         IF(IROTAT(J).EQ.2) THEN
C
C     FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
            CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
            OMEGA=DATA(6)
#ifdef FT_DEBUG
            write(88,"('7 rot vel (3 is wheel speed):',1p,3e12.3)")
     +           data(4),data(5),data(6)
#endif
C
C     FIND THE CHASSIS PITCH VELOCITY FROM THE MODEL
C
            CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
            WY = DATA(5)
         ENDIF
C
C     FIND THE NO SLIP FORWARD VELOCITY
C
         VELWC = RR ( J ) * (OMEGA + WY)
         ABSVWC = DABS ( VELWC )
         ABSUW = DABS ( UW )
C
C     FIND THE LONGITUDINAL SLIPPAGE
C
         SLIP = ( ABSVWC - UW ) / DMAX1 ( ABSVWC, ABSUW )
         ASLIP = DABS ( SLIP )
C
C     CALCULATE THE MOBILITY NUMBER
C
         IF(ISURF(J).EQ.1) THEN
            CN = (CI(J)*B(J)*D(J))/VF
            BN1 = 1.0 + 5.0 * VDEFL/H(J)
            BN2 = 1.0 + 3.0 * B(J)/D(J)
            BN = CN * BN1/BN2
         ENDIF
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED
C
C     FIND THE TRACTIVE COEFFICIENT
C
         IF(ISURF(J).EQ.1) THEN
            IF(ICONST(J).EQ.1) POWER1 = -0.1*BN
            IF(ICONST(J).EQ.2) POWER1 = -0.1*BN
            Q1 = 1.0 - DEXP(POWER1)
            IF(ICONST(J).EQ.1) POWER2 = -7.5*ASLIP
            IF(ICONST(J).EQ.2) POWER2 = -9.5*ASLIP
            Q2 = 1.0 - DEXP(POWER2)
            IF(ICONST(J).EQ.1) AMU = 0.88 * Q1 * Q2 + 0.04
            IF(ICONST(J).EQ.2) AMU = 0.88 * Q1 * Q2 + 0.0325
         ENDIF
         IF(ISURF(J).EQ.2) THEN
            CN=(CI(J)*B(J)*D(J))/VF
            POWER3 = -CN*ASLIP
C     IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.02
            IF(ICONST(J).EQ.1) AMU = 1.02 * (1.0 - DEXP(POWER3))
            IF(ICONST(J).EQ.1) AMU = DSIGN(AMU,SLIP) + 0.02
C     IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3)) + 0.01625
            IF(ICONST(J).EQ.2) AMU = 1.02 * (1.0 - DEXP(POWER3))
            IF(ICONST(J).EQ.2) AMU = DSIGN(AMU,SLIP) + 0.01625
         ENDIF
C
C     FINALLY THE TRACTIVE FORCE
C
         TRAC = AMU * VF
C
C     THEN THE TORQUE IS THE PRODUCT OF THE TRACTIVE FORCE
C     AND THE ROLLING RADIUS
C
         FMAG = - TRAC * RR(J)
         RESULT(1)=VF
         RESULT(2)=UW
         RESULT(3)=OMEGA
         RESULT(4)=WY
         RESULT(5)=SLIP
         IF(ISURF(J).EQ.1)RESULT(6)=BN
         IF(ISURF(J).EQ.2)RESULT(6)=CN
         RESULT(7)=-TRAC
         RESULT(8)=FMAG
C
      ENDIF
C
C     ALLOW THE OPTION OF INCLUDING A FORE-AFT SPRING AND DAMPER
C
      IF(ILONG(J).EQ.1) THEN
C
C     FIND THE FORE-AFT POSITION OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         XSPACE=DATA(1)
C
C     FIND THE FORE-AFT POSITION OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
         CALL INFO('DISP',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
         XRA=DATA(1)
C
C     FIND THE FORE-AFT VELOCITY OF THE WHEEL IN THE GLOBAL SYSTEM
C
         CALL INFO('VEL',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
         XDOT=DATA(1)
C
C     FIND THE FORE-AFT VELOCITY OF THE REAR AXLE IN THE GLOBAL SYSTEM
C
         CALL INFO('VEL',IRAMKR,IGMKR,IGMKR,DATA,ERRFLG)
         XRADOT=DATA(1)
C
C     FIND THE FORE-AFT FORCE
C
         XF=AKFA(J)*(XSPACE-(XRA+XSE(J)))+CFA(J)*(XDOT-XRADOT)
C
C     CHECK TO SEE IF FRICTION LIMITS THE FORCE
C
         IF(IXFRIC(J).EQ.1) THEN
            XFS=XF
            REQMU=XF/VF
            REQMU=DABS(REQMU)
            IF(REQMU.GE.AMUFA(J)) THEN
               XF=AMUFA(J)*VF
               XF=DSIGN(XF,XFS)
            ENDIF
         ENDIF
C
C     ADD THE MOMENT OF THE FORE-AFT SPRING FORCE TO THE TRACTIVE FORCE
C     MOMENT (IF THE LATTER IS NONZERO)
C
         FMAG=FMAG+XF*XARM(J)
         RESULT(2)=VF
         RESULT(6)=XF
         RESULT(7)=XF*XARM(J)
      ENDIF
      RETURN
C
C     THIS SECTION COULD BE USED TO CALCULATE A MOMENT ABOUT A VERTICAL
C     AXIS.  AN EXAMPLE OF SUCH A TORQUE WOULD BE AN ALIGNING TORQUE.
C     BECAUSE OF LIMITED DATA CONCERNING THE ALIGNING TORQUES ACTING
C     ON OFF-ROAD TIRES, SUCH CALCULATIONS HAVE NOT BEEN INCLUDED HERE.
C
C     Y A W   M O M E N T
C
 270  IF(IDR.GT.700)GO TO 280
C
C     SET MOMENT TO ZERO IF ID IS GREATER THAN 600 BUT LESS THAN 700
C
      FMAG=0.0D0
      RETURN
C
C     E N G I N E   T O R Q U E
C
C
 280  IF(IDR.GT.800)GO TO 290
C
C     THIS SECTION IS USED FOR STATIC EQUILIBRIUM
C
      IF ( TIME .GT. 0.0D0 )   GO TO 281
C
C     ENGINE TORQUE IS ZERO FOR STATIC EQUILIBRIUM
C
      RESULT(2)=0.0
      RESULT(3)=0.0
      RESULT(4)=0.0
C
 281  CONTINUE
C
C     FIND THE ENGINE SPEED
C
      IMARKR=PAR(1)
      JMARKR=PAR(2)
      CALL INFO('VEL',IMARKR,JMARKR,JMARKR,DATA,ERRFLG)
      ESPEED=DATA(6)*(60.0/(2.0*PI))
      RESULT(2)=ESPEED
C
C     INTERPOLATE TO FIND THE ENGINE TORQUE
C
      CALL VBM(ESPEED,ETORQ,ENGSPD,ENGTOR,NTSPTS)
      RESULT(3)=ETORQ
      ETORQ=ETA*ETORQ
      RESULT(4)=ETORQ
      RETURN
C
C     T R A C T I V E      E F F I C I E N C Y
C
 290  IF(IDR.GT.900) GO TO 300
C
C     FIND THE FORWARD VELOCITY
C
      IMARKR=IRAMKR
      JMARKR=IGMKR
      IRMKR=IRAMKR
      CALL INFO('VEL',IMARKR,JMARKR,IRMKR,DATA,ERRFLG)
      V=DABS(DATA(1))
      RESULT(3)=V
C
C     FIND THE NET TRACTIVE FORCE AND INPUT POWER BY SUMMING OVER
C     THE NUMBER OF TIRES
C
      TPULL=0.0
      TOTALI=0.0
C
      DO 291 J=1,NTIRES
C
C     FIND THE NET TRACTIVE FORCE
C
         IMARKR=PAR(J)
         JMARKR=0
         IRMKR=IRAMKR
         CALL INFO('FORC',IMARKR,JMARKR,IRMKR,DATA,ERRFLG)
         PULL=DABS(DATA(1))
C
C     FIND THE TOTAL NET TRACTIVE FORCE
C
         TPULL=TPULL+PULL
C
C     FIND THE LOADING TORQUE
C
         IPAR=J+8
         IMARKR=PAR(IPAR)
         JMARKR=0
         IRMKR=IWROT(J)
         CALL INFO('FORC',IMARKR,JMARKR,IRMKR,DATA,ERRFLG)
         TORQUE=DABS(DATA(6))
C
C     FIND THE WHEEL ANGULAR VELOCITY
C
         IMARKR=IWROT(J)
         JMARKR=JWROT(J)
         IRMKR=JWROT(J)
         CALL INFO('VEL',IMARKR,JMARKR,IRMKR,DATA,ERRFLG)
         OMEGA=DABS(DATA(6))

#ifdef FT_DEBUG
         write(88,"('8 rot vel (3 is wheel speed):',1p,3e12.3)")
     +        data(4),data(5),data(6)
#endif
C
C     FIND THE AXLE POWER
C
         TOTALI=TOTALI+TORQUE*OMEGA
C
C     FIND THE WHEEL SLIPPAGE
C
C     FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
         CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
         UW = DATA (1)
C
C     FIND THE CHASSIS PITCH VELOCITY FROM THE MODEL
C
         CALL INFO( 'VEL', IVMKR(J), IGMKR, IGMKR, DATA, ERRFLG )
         WY = DATA(5)
C
C     FIND THE NO SLIP FORWARD VELOCITY
C
         VELWC = RR ( J ) * (OMEGA + WY)
         ABSVWC = DABS ( VELWC )
         ABSUW = DABS ( UW )
C
C     FIND THE LONGITUDINAL SLIPPAGE
C
         SLIP = ( ABSVWC - UW ) / DMAX1 ( ABSVWC, ABSUW )
         IF(J.EQ.1)RESULT(4)=SLIP
         IF(J.EQ.5)RESULT(6)=SLIP
C
 291  CONTINUE
C
C     FINALLY CALCULATE THE TRACTIVE EFFICIENCY
C
      IF(TOTALI.NE.0.0)TE=(TPULL*V)/TOTALI
      IF(TOTALI.EQ.0.0)TE=0.0
      RESULT(2)=TPULL
      RESULT(7)=TOTALI
      RESULT(8)=TE
      RETURN
C
C
 300  IF ( IDR.GT. 1000) GO TO 310
C
C     L O N G I T U D I N A L   F O R C E S-B R A K I N G
C
C     IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C     1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C     2 - POWERED WHEEL - CONSIDER TRACTION/BRAKING FORCE
C
C     TWO TYPES OF FORE-AFT OR LONGITUDINAL FORCES ARE CONSIDERED IN THIS
C     SECTION.  THESE ARE:
C     1. ROLLING RESISTANCE FORCES ON UNPOWERED WHEELS
C     2. TRACTIVE/BRAKING FORCES ON POWERED/BRAKED WHEELS
C     THESE MODELS ARE DERIVED FROM AN EARLIER DRAM MODEL FOR COMBINE BRAKING 036700
C
C     THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 970
C
C     CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 960
      J=IDR-900
C
C     ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C     POSITION) WHERE KX=PAR(2) AND DESIRED X POSITION IS PAR(3)
C
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3))
 960  RETURN
C
C     ***************
C
 970  CONTINUE
C
C     DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C     FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = IDR- 900
C
C     CALCULATE THE VERTICAL FORCE
C
C     FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
         X1(J)=DATA(1)
         Y1(J)=DATA(3)
      ENDIF
C
C     FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C     FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
         VF=0.0D0
         FMAG=0.0D0
         RETURN
      ENDIF
C
C     FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
         CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
         VFORC=FZ(J)
C     IF THE DISTRIBUTED CONTACT MODEL IS BEING USED, INCLUDE THE FORE-AFT
C     COMPONENT
         FMAG=-FX(J)
         CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C     FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C     FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C     FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C     CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C     IF SO, SET FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C     FIND THE ROLLING RESISTANCE OR TRACTIVE/BRAKING FORCE
C
C     FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
      CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
      UW = DATA (1)
      ABSUW = DABS ( UW )
C
C     FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
      CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
      OMEGA=DATA(6)

#ifdef FT_DEBUG
      write(88,"('9 rot vel (3 is wheel speed):',1p,3e12.3)")
     +     data(4),data(5),data(6)
#endif
C
C     CALCULATE THE ROLLING RESISTANCE FORCE ON AN UNPOWERED WHEEL
C
      IF(IPOWER(J).EQ.1) THEN
C
C     FIND THE NO SLIP FORWARD VELOCITY
C
         VELWC = PAR(6) * OMEGA
         ABSVWC = DABS ( VELWC )
         VT = ABSUW - ABSVWC
         VTMAG = DABS(VT)
C
C     CALCULATE THE ROLLING RESISTANCE COEFFICIENT, RHO
C
         IF(VTMAG.LE.PAR(5)) THEN
            ARG=POT*(VT/PAR(5))
            RHO=PAR(4)*DSIN(ARG)
         ENDIF
         IF(VTMAG.GT.PAR(5)) RHO=DSIGN(PAR(4),VT)
         RHO=RHO*DSIGN(1.0D0,OMEGA)
C
C     FINALLY CALCULATE THE ROLLING RESISTANCE FORCE
C
         TF = RHO * VF
         FMAG = TF
         RESULT(2)=OMEGA
         RESULT(3)=UW
         RESULT(4)=VT
         RESULT(6)=RHO
         RESULT(7)=VF
         RESULT(8)=FMAG
         RETURN
      ENDIF
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED OR BRAKED
C
      IF(IPOWER(J).EQ.2) THEN
C
C     FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = RR ( J ) * OMEGA
      ABSVWC = DABS ( VELWC )
C
C     FIND THE LONGITUDINAL SLIPPAGE
C
      RELVEL = ABSVWC - ABSUW
      SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW, PAR(2))
C     SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW )
C     SLIP = RELVEL / ABSUW
      ASLIP = DABS ( SLIP )
C
C     FIND THE TRACTIVE COEFFICIENT
C
      CN=(CI(J)*B(J)*D(J))/VF
      POWER3 = -CN*ASLIP
      AMU = 1.02 * (1.0 - DEXP(POWER3))
      AMU = DSIGN(AMU,SLIP)
C
C     FINALLY THE TRACTIVE FORCE
C
      TRAC = AMU * VF * DSIGN(1.0D0,OMEGA)
      FMAG = - TRAC
      RESULT(1)=CN
      RESULT(2)=VELWC
      RESULT(3)=UW
      RESULT(4)=SLIP
      RESULT(6)=AMU
      RESULT(7)=VF
      RESULT(8)=FMAG
      RETURN
      ENDIF
C
 310  IF ( IDR.GT. 1100) GO TO 320
C
C     P I T C H  M O M E N T- A L T E R N A T E   M O D E L
C
C     IPOWER(J) = 0 - DON'T CONSIDER ROLLING RESISTANCE AND TRACTION FORCES
C     1 - UNPOWERED WHEEL - CONSIDER ROLLING RESISTANCE
C     2 - POWERED WHEEL - CONSIDER TRACTION/BRAKING FORCE
C
C     TWO TYPES OF FORE-AFT OR LONGITUDINAL FORCES ARE CONSIDERED IN THIS
C     SECTION.  THESE ARE:
C     1. ROLLING RESISTANCE FORCES ON UNPOWERED WHEELS
C     2. TRACTIVE/BRAKING FORCES ON POWERED/BRAKED WHEELS
C     THESE MODELS ARE DERIVED FROM AN EARLIER DRAM MODEL
C
C     THIS SECTION IS USED FOR STATIC EQUILIBRIUM IN THE F-A DIRECTION
C
      IF ( TIME .GT. 0.0D0 )   GO TO 1070
C
C     CALCULATE THE FORE-AFT FORCES TO PLACE TRACTOR IN STATIC EQUILIBRIUM
C
      IF (PAR(1).NE.1) GO TO 1060
      J=IDR-1000
C
C     ASSUME SPRING IN FORE-AFT DIRECTION (FMAG = KX*(DATA(1)-DESIRED X
C     POSITION) WHERE KX=PAR(2) AND DESIRED X POSITION IS PAR(3).
C     PAR(4) = MOMENT ARM OF FORCE
      CALL INFO('DISP',IMRKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      FMAG = PAR(2) * (DATA(1)-PAR(3)) * PAR(4)
 1060 RETURN
C
C     ***************
C
 1070 CONTINUE
C
C     DETERMINE THE LONGITUDINAL FORCE FOR THE DYNAMIC CONDITION
C
C     FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = IDR- 1000
C
C     CALCULATE THE VERTICAL FORCE
C
C     FIRST FIND THE VERTICAL DISPLACEMENT OF THE WHEEL CENTER WRT GROUND
C
      CALL INFO('DISP',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      AXLHGT=DATA(3)
      IF(IVERT(J).EQ.3) THEN
         X1(J)=DATA(1)
         Y1(J)=DATA(3)
      ENDIF
C
C     FIND THE TIRE DEFLECTION
C
      VDISP=AXLHGT-WRAD(J)
C
C     FOR IVERT(J) = 1 OR 2, THE VERTICAL FORCE IS 0 IF VDISP .GE. 0
C
      IF(VDISP.GE.0.0.AND.IVERT(J).NE.3) THEN
         VF=0.0D0
         FMAG=0.0D0
         RETURN
      ENDIF
C
C     FIND THE SPRING FORCE PORTION OF THE VERTICAL FORCE
C
      IF(VDISP.LT.0.0D0.AND.IVERT(J).NE.3)VDEFL=DABS(VDISP)
      IF(IVERT(J).EQ.1)VFORC=AKV(J)*VDEFL
      IF(IVERT(J).EQ.2)CALL VBM2(J,VDEFL,VFORC,DEFL,FORCE,NLDPT)
      IF(IVERT(J).EQ.3) THEN
         CALL TIREF_fdm(J,time,zeroDelta,fx(j),fz(j))
         VFORC=FZ(J)
C     IF THE DISTRIBUTED CONTACT MODEL IS BEING USED, INCLUDE THE FORE-AFT
C     COMPONENT
         FMAG=-FX(J)
         CALL VBM2(J,FZ(J),VDEFL,FORCE,DEFL,NLDPT)
      ENDIF
C
C     FIND THE VERTICAL VELOCITY OF THE WHEEL CENTER
C
      CALL INFO('VEL',IWMKR(J),IGMKR,IGMKR,DATA,ERRFLG)
      VVEL=DATA(3)
C
C     FIND THE VERTICAL DAMPING COEFFICIENT
C
      CDAMP=CDAMPZ(J)*(1.0D0-DEXP(-3.0D0*VDEFL/CDEFL(J)))
C
C     FIND THE VERTICAL FORCE
C
      VF = VFORC - CDAMP * VVEL
C
C     CHECK TO SEE IF WHEEL HAS LEFT GROUND.
C     IF SO, SET FORE-AFT FORCE TO ZERO
C
      IF(VF.LE.0.0D0) GO TO 240
C
C     FIND THE ROLLING RESISTANCE OR TRACTIVE/BRAKING FORCE
C
C     FIRST DETERMINE THE VELOCITY COMPONENT IN THE PLANE OF THE WHEEL
C
      CALL INFO ( 'VEL', IWMKR(J), IGMKR, IWMKR(J), DATA, ERRFLG )
      UW = DATA (1)
      ABSUW = DABS ( UW )
C
C     FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
      CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
      OMEGA=DATA(6)
#ifdef FT_DEBUG
      write(88,"('10 rot vel (3 is wheel speed):',1p,3e12.3)")
     +     data(4),data(5),data(6)
#endif
C
C     CALCULATE THE ROLLING RESISTANCE FORCE ON AN UNPOWERED WHEEL
C
      IF(IPOWER(J).EQ.1) THEN
C
C     FIND THE NO SLIP FORWARD VELOCITY
C
         VELWC = PAR(4) * OMEGA
         ABSVWC = DABS ( VELWC )
         VT = ABSUW - ABSVWC
         VTMAG = DABS(VT)
C
C     CALCULATE THE ROLLING RESISTANCE COEFFICIENT, RHO
C
         IF(VTMAG.LE.PAR(6)) THEN
            ARG=POT*(VT/PAR(6))
            RHO=PAR(5)*DSIN(ARG)
         ENDIF
         IF(VTMAG.GT.PAR(6)) RHO=DSIGN(PAR(5),VT)
         RHO=RHO*DSIGN(1.0D0,OMEGA)
C
C     FINALLY CALCULATE THE ROLLING RESISTANCE FORCE AND ITS MOMENT
C
         TF = RHO * VF
         FMAG = TF * PAR(4)
         RESULT(2)=OMEGA
         RESULT(3)=UW
         RESULT(4)=VT
         RESULT(6)=RHO
         RESULT(7)=VF
         RESULT(8)=FMAG
         RETURN
      ENDIF
C
C     DETERMINE THE TRACTIVE FORCE IF THE WHEEL IS POWERED OR BRAKED
C
      IF(IPOWER(J).EQ.2) THEN
C
C     FIND THE NO SLIP FORWARD VELOCITY
C
      VELWC = RR ( J ) * OMEGA
      ABSVWC = DABS ( VELWC )
C
C     FIND THE LONGITUDINAL SLIPPAGE
C
      RELVEL = ABSVWC - ABSUW
      SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW, PAR(2))
C     SLIP = RELVEL / DMAX1 ( ABSVWC, ABSUW )
C     SLIP = RELVEL / ABSUW
      ASLIP = DABS ( SLIP )
C
C     FIND THE TRACTIVE COEFFICIENT
C
      CN=(CI(J)*B(J)*D(J))/VF
      POWER3 = -CN*ASLIP
      AMU = 1.02 * (1.0 - DEXP(POWER3))
      AMU = DSIGN(AMU,SLIP)
C
C     FINALLY THE TRACTIVE FORCE AND ITS MOMENT
C
      TRAC = AMU * VF * DSIGN(1.0D0,OMEGA)
      FMAG = - TRAC * RR(J)
      RESULT(1)=CN
      RESULT(2)=VELWC
      RESULT(3)=UW
      RESULT(4)=SLIP
      RESULT(6)=AMU
      RESULT(7)=VF
      RESULT(8)=FMAG
      RETURN
      ENDIF
C
C
 320  IF ( IDR.GE. 1200) GO TO 330
C
C     B R A K I N G   M O M E N T
C
C     DETERMINE THE BRAKING MOMENT FOR THE DYNAMIC CONDITION
C
C     FIND THE WHEEL INDEX J FROM THE FORCE ID
C
      J = IDR- 1100
C
C     FIND THE WHEEL ANGULAR VELOCITY FROM THE MODEL
C
      CALL INFO('VEL',IWROT(J),JWROT(J),JWROT(J),DATA,ERRFLG)
      OMEGA=DATA(6)
      ABSAV=DABS(OMEGA)

#ifdef FT_DEBUG
      write(88,"('11 rot vel (3 is wheel speed):',1p,3e12.3)")
     +     data(4),data(5),data(6)
#endif
C
C     FIND THE MAXIMUM TORQUE AS A FUNCTION OF TIME
C
      OMEGA1=PI/PAR(3)
      TBEGIN=PAR(4)
      TEND=PAR(4)+PAR(3)
      IF(TIME.LE.TBEGIN) TMAX=0.0D0
      IF(TIME.GT.TBEGIN.AND.TIME.LE.TEND) THEN
      TSTAR=TIME-TBEGIN
      TMAX=(PAR(1)/2.0D0)*(1.0D0-DCOS(OMEGA1*TSTAR))
      ENDIF
      IF(TIME.GT.TEND) TMAX=PAR(1)
C
C     FIND THE BRAKING TORQUE AS A FUNCTION OF THE RELATIVE ANGULAR VELOCITY00039100
C
      IF(ABSAV.LE.PAR(2)) THEN
         OMEGA2=POT*(OMEGA/PAR(2))
         TBRAKE=-TMAX*DSIN(OMEGA2)
      ENDIF
      IF(ABSAV.GT.PAR(2)) TBRAKE=-DSIGN(TMAX,OMEGA)
C
C     FINALLY CALCULATE THE MOMENT
C
      FMAG = TBRAKE
      RESULT(2)=TMAX
      RESULT(3)=OMEGA
      RESULT(4)=FMAG
      RETURN
C
C     *************
C
 330  IF ( IDR.GT. 1300) GO TO 6106
C
C     O U T P U T     F O R C E S   F O R     F E A
C
C     PAR(1)=REFERENCE MARKER AT CG OF PART
C     PAR(2)=ACCELERATION OF GRAVITY
C     PAR(3)=MARKER INDICATING DIRECTION OF GRAVITY
C     PAR(4)=MASS OF PART
C     PAR(5)=IXX OF PART
C     PAR(6)=IYY OF PART
C     PAR(7)=IZZ OF PART
C     PAR(8)=IXY OF PART
C     PAR(9)=IXZ OF PART
C     PAR(10)=IYZ OF PART
C     PAR(11)=NUMBER OF APPLIED FORCES=N
C     PAR(12,14,16,ETC)=I MARKER FOR FORCES
C     PAR(13,15,17,ETC)=J MARKER FOR FORCES
C
C
C     FIND THE MOMENT ARMS OF THE FORCES
C
      IF(IDR.EQ.1200.AND.IFPASS.EQ.0) THEN
c
c     open the CAEDS finite element output file
c
C     CFILE='                              '
C     write(*,'(A)') 'Enter file name for CAEDS input file'
C     read(*,'(A)') CFILE
C     OPEN(34,FILE=CFILE)
         OPEN(33,FILE='baler.fea',STATUS='OLD')
         OPEN(34,FILE='baler.fea2',STATUS='OLD')
         NFORCS=PAR(11)
         IRM=PAR(1)
         DO 338 I=1,NFORCS
            N1=12+2*(I-1)
            IMARK=PAR(N1)
            CALL INFO('DISP',IMARK,IRM,IRM,DATA,ERRFLG)
            XDIST(I)=DATA(1)
            YDIST(I)=DATA(2)
            ZDIST(I)=DATA(3)
            WRITE(34,335) XDIST(I),YDIST(I),ZDIST(I)
 338     CONTINUE
         IFPASS=1
         RETURN
      ENDIF
C
      IF(IDR.EQ.1201) WRITE(33,334)TIME
 334  FORMAT(5X,F7.3)
C
C     FIND THE COMPONENTS OF THE GRAVITY VECTOR
C
      IRM=PAR(1)
      IGRAVM=PAR(3)
      CALL INFO('DISP',IGRAVM,IGMKR,IRM,DATA,ERRFLG)
      GX=PAR(2)*DATA(1)*PAR(4)/1000.0D0
      GY=PAR(2)*DATA(2)*PAR(4)/1000.0D0
      GZ=PAR(2)*DATA(3)*PAR(4)/1000.0D0
      GX1=PAR(2)*DATA(1)
      GY1=PAR(2)*DATA(2)
      GZ1=PAR(2)*DATA(3)
      IF(IDR.EQ.1201) THEN
         RESULT(2)=GX
         RESULT(3)=GY
         RESULT(4)=GZ
         WRITE(33,335) GX,GY,GZ
 335     FORMAT(8X,6F16.2)
         RETURN
      ENDIF
C
C     FIND THE ANGULAR VELOCITY OF THE PART
C
      CALL INFO('VEL',IRM,IGMKR,IRM,DATA,ERRFLG)
      VX=DATA(1)
      VY=DATA(2)
      VZ=DATA(3)
      WX=DATA(4)
      WY=DATA(5)
      WZ=DATA(6)
      IF(IDR.EQ.1202) THEN
         RESULT(2)=VX
         RESULT(3)=VY
         RESULT(4)=VZ
         RESULT(6)=WX
         RESULT(7)=WY
         RESULT(8)=WZ
         WXRPS=WX/(2.0D0*PI)
         WYRPS=WY/(2.0D0*PI)
         WZRPS=WZ/(2.0D0*PI)
         WRITE(33,339) VX,VY,VZ,WX,WY,WZ
 339     FORMAT(8X,3F16.3,3F16.6)
         RETURN
      ENDIF
C
C     FIND THE ACCELERATION OF THE CG
C
      CALL INFO('ACC',IRM,IGMKR,IRM,DATA,ERRFLG)
      AX=DATA(1)
      AY=DATA(2)
      AZ=DATA(3)
      WXDOT=DATA(4)
      WYDOT=DATA(5)
      WZDOT=DATA(6)
      IF(IDR.EQ.1203) THEN
         RESULT(2)=AX
         RESULT(3)=AY
         RESULT(4)=AZ
         RESULT(6)=WXDOT
         RESULT(7)=WYDOT
         RESULT(8)=WZDOT
         WXDRPS=WXDOT/(2.0D0*PI)
         WYDRPS=WYDOT/(2.0D0*PI)
         WZDRPS=WZDOT/(2.0D0*PI)
         WRITE(33,339) AX,AY,AZ,WXDOT,WYDOT,WZDOT
         RETURN
      ENDIF
C
C     FIND THE MOMENT ARMS OF THE FORCES
C
      NFORCS=PAR(11)
      DO 331 I=1,NFORCS
         N1=12+2*(I-1)
         IMARK=PAR(N1)
         CALL INFO('DISP',IMARK,IRM,IRM,DATA,ERRFLG)
         XDIST(I)=DATA(1)
         YDIST(I)=DATA(2)
         ZDIST(I)=DATA(3)
 331  CONTINUE
C
      IF(IDR.EQ.1204) THEN
         DO 336 I=1,NFORCS
            WRITE(33,335) XDIST(I),YDIST(I),ZDIST(I)
 336     CONTINUE
         RETURN
      ENDIF
C
C     FIND THE FORCES AND MOMENTS
C
      DO 332 I=1,NFORCS
         N1=12+2*(I-1)
         N2=13+2*(I-1)
         IMARK=PAR(N1)
         JMARK=PAR(N2)
         CALL INFO('FORC',IMARK,JMARK,IRM,DATA,ERRFLG)
         XFORC(I)=DATA(1)
         YFORC(I)=DATA(2)
         ZFORC(I)=DATA(3)
         XMOM(I)=DATA(4)
         YMOM(I)=DATA(5)
         ZMOM(I)=DATA(6)
 332  CONTINUE
C
      IF(IDR.EQ.1205) THEN
         RESULT(2)=XFORC(1)
         RESULT(3)=YFORC(1)
         RESULT(4)=ZFORC(1)
         RESULT(6)=XMOM(1)
         RESULT(7)=YMOM(1)
         RESULT(8)=ZMOM(1)
         WRITE(33,335)XFORC(1),YFORC(1),ZFORC(1),XMOM(1),YMOM(1),ZMOM(1)
         RETURN
      ENDIF
C
      IF(IDR.EQ.1206) THEN
         RESULT(2)=XFORC(2)
         RESULT(3)=YFORC(2)
         RESULT(4)=ZFORC(2)
         RESULT(6)=XMOM(2)
         RESULT(7)=YMOM(2)
         RESULT(8)=ZMOM(2)
         WRITE(33,335)XFORC(2),YFORC(2),ZFORC(2),XMOM(2),YMOM(2),ZMOM(2)
         RETURN
      ENDIF
C
      IF(IDR.EQ.1207) THEN
         RESULT(2)=XFORC(3)
         RESULT(3)=YFORC(3)
         RESULT(4)=ZFORC(3)
         RESULT(6)=XMOM(3)
         RESULT(7)=YMOM(3)
         RESULT(8)=ZMOM(3)
         WRITE(33,335)XFORC(3),YFORC(3),ZFORC(3),XMOM(3),YMOM(3),ZMOM(3)
         RETURN
      ENDIF
C
      IF(IDR.EQ.1208) THEN
         RESULT(2)=XFORC(4)
         RESULT(3)=YFORC(4)
         RESULT(4)=ZFORC(4)
         RESULT(6)=XMOM(4)
         RESULT(7)=YMOM(4)
         RESULT(8)=ZMOM(4)
         WRITE(33,335)XFORC(4),YFORC(4),ZFORC(4),XMOM(4),YMOM(4),ZMOM(4)
         RETURN
      ENDIF
C
      IF(IDR.EQ.1209) THEN
         RESULT(2)=XFORC(5)
         RESULT(3)=YFORC(5)
         RESULT(4)=ZFORC(5)
         RESULT(6)=XMOM(5)
         RESULT(7)=YMOM(5)
         RESULT(8)=ZMOM(5)
         WRITE(33,335)XFORC(5),YFORC(5),ZFORC(5),XMOM(5),YMOM(5),ZMOM(5)
         RETURN
      ENDIF
C
      IF(IDR.EQ.1210) THEN
         RESULT(2)=XFORC(6)
         RESULT(3)=YFORC(6)
         RESULT(4)=ZFORC(6)
         RESULT(6)=XMOM(6)
         RESULT(7)=YMOM(6)
         RESULT(8)=ZMOM(6)
         WRITE(33,335)XFORC(6),YFORC(6),ZFORC(6),XMOM(6),YMOM(6),ZMOM(6)
         RETURN
      ENDIF
C
C     SUM THE FORCES AND MOMENTS
C
      SUMFX=0.0D0
      SUMFY=0.0D0
      SUMFZ=0.0D0
      SUMMX=0.0D0
      SUMMY=0.0D0
      SUMMZ=0.0D0
      DO 333 I=1,NFORCS
         SUMFX=SUMFX+XFORC(I)
         SUMFY=SUMFY+YFORC(I)
         SUMFZ=SUMFZ+ZFORC(I)
         SUMMX=SUMMX+XMOM(I)+ZFORC(I)*YDIST(I)-YFORC(I)*ZDIST(I)
         SUMMY=SUMMY+YMOM(I)+XFORC(I)*ZDIST(I)-ZFORC(I)*XDIST(I)
         SUMMZ=SUMMZ+ZMOM(I)+YFORC(I)*XDIST(I)-XFORC(I)*YDIST(I)
 333  CONTINUE
C
      IF(IDR.EQ.1211) THEN
         RESULT(2)=SUMFX
         RESULT(3)=SUMFY
         RESULT(4)=SUMFZ
         RESULT(6)=SUMMX
         RESULT(7)=SUMMY
         RESULT(8)=SUMMZ
         WRITE(33,335) SUMFX,SUMFY,SUMFZ,SUMMX,SUMMY,SUMMZ
         RETURN
      ENDIF
C
C     ADD THE GRAVITATIONAL FORCE
C
      SUMFX=SUMFX+GX
      SUMFY=SUMFY+GY
      SUMFZ=SUMFZ+GZ
C
      IF(IDR.EQ.1212) THEN
         RESULT(2)=SUMFX
         RESULT(3)=SUMFY
         RESULT(4)=SUMFZ
         RESULT(6)=SUMMX
         RESULT(7)=SUMMY
         RESULT(8)=SUMMZ
         WRITE(33,335) SUMFX,SUMFY,SUMFZ,SUMMX,SUMMY,SUMMZ
         RETURN
      ENDIF
C
C     FIND THE "INERTIAL" FORCES AND MOMENTS
C
      AMASS=PAR(4)
      AIXX=PAR(5)
      AIYY=PAR(6)
      AIZZ=PAR(7)
      AIXY=PAR(8)
      AIXZ=PAR(9)
      AIYZ=PAR(10)
      XINERF=AMASS*AX/1000.0D0
      YINERF=AMASS*AY/1000.0D0
      ZINERF=AMASS*AZ/1000.0D0
      XINERM=(AIXX*WXDOT+(AIZZ-AIYY)*WY*WZ+AIXY*(WZ*WX-WYDOT)
     1     -AIXZ*(WZDOT+WY*WX)-AIYZ*(WY**2-WZ**2))/1000.0D0
      YINERM=(AIYY*WYDOT+(AIXX-AIZZ)*WZ*WX+AIYZ*(WX*WY-WZDOT)
     1     -AIXY*(WXDOT+WZ*WY)-AIXZ*(WZ**2-WX**2))/1000.0D0
      ZINERM=(AIZZ*WZDOT+(AIYY-AIXX)*WX*WY+AIXZ*(WY*WZ-WXDOT)
     1     -AIYZ*(WYDOT+WX*WZ)-AIXY*(WX**2-WY**2))/1000.0D0
C
      IF(IDR.EQ.1213) THEN
         RESULT(2)=XINERF
         RESULT(3)=YINERF
         RESULT(4)=ZINERF
         RESULT(6)=XINERM
         RESULT(7)=YINERM
         RESULT(8)=ZINERM
         WRITE(33,335) XINERF,YINERF,ZINERF,XINERM,YINERM,ZINERM
C
C     WRITE THE OUTPUT FILE FOR THE DURABILITY ANALYSIS
C
         WRITE(34,345)TIME
     1        ,XFORC(1),YFORC(1),ZFORC(1)
     2        ,XFORC(2),YFORC(2),ZFORC(2),XMOM(2),YMOM(2),ZMOM(2)
     3        ,XFORC(3),YFORC(3),ZFORC(3),XMOM(3),YMOM(3),ZMOM(3)
C    4,GX,GY,GZ,AX,AY,AZ,WXDOT,WYDOT,WZDOT,WX,WY,WZ
 345     FORMAT(1P16E13.5)
         IFPASS=IFPASS+1
         RETURN
      ENDIF
C
      RETURN
C
 6106 IF(IDR.GT.1400) GO TO 6200
C
C     S H O C K   A B S O R B E R S
C
      J=IDR-1300
C
C     FIRST FIND THE RELATIVE VELOCITY IN GLOBAL COORDINATES
C
      CALL INFO('VEL',ISHCK(J),JSHCK(J),IGMKR,DATA,ERRFLG)
      XVEL=DATA(1)
      YVEL=DATA(2)
      ZVEL=DATA(3)
C
C     NEXT FIND THE RELATIVE DISPLACEMENT IN GLOBAL COORDINATES
C
      CALL INFO('DISP',ISHCK(J),JSHCK(J),IGMKR,DATA,ERRFLG)
      XDISP=DATA(1)
      YDISP=DATA(2)
      ZDISP=DATA(3)
      DIST=DSQRT(XDISP**2+YDISP**2+ZDISP**2)
C
C     TAKE THE DOT PRODUCT TO CALCULATE THE RELATIVE VELOCITY
C
      RVEL=(XVEL*XDISP+YVEL*YDISP+ZVEL*ZDISP)/DIST
      ARVEL=DABS(RVEL)
C
      IF(RVEL.LE.0.0D0)GO TO 6107
C
C     SHOCK ABSORBER IS IN EXTENSION
C
      CALL VBM2(J,ARVEL,DFORCE,VELEX,FOREX,NEPTS)
      GO TO 6108
C
C     SHOCK ABSORBER IS IN COMPRESSION
C
 6107 CALL VBM2(J,ARVEL,DFORCE,VELCP,FORCP,NCPTS)
C
 6108 IF(ARVEL.EQ.0.0D0)GO TO 6109
      FMAG=-(RVEL/ARVEL)*DFORCE
      RESULT(2)=RVEL
      RESULT(3)=DFORCE
      RESULT(4)=FMAG
      RETURN
 6109 FMAG=0.0D0
      RETURN

 6200 IF(IDR.GT.1500) GO TO 6300
C
C     O U T P U T   A C C E L E R A T I O N   T I M E    H I S T O R Y
C
      IF(IFPSS.EQ.0) THEN
         ICOUNT=0
         IFPSS=1
      ENDIF
C
C     FIRST FIND THE ACCELERATION OF THE MARKER IN QUESTION
C
      IMKR=PAR(1)
      CALL INFO('ACC',IMKR,IGMKR,IGMKR,DATA,ERRFLG)
      RESULT(2)=DATA(1)
      RESULT(3)=DATA(2)
      RESULT(4)=DATA(3)
      RESULT(6)=DATA(4)
      RESULT(7)=DATA(5)
      RESULT(8)=DATA(6)
C
C     STORE THE DESIRED ACCELERATION IN AN ARRAY
C
      ICOUNT=ICOUNT+1
      RESULT(1)=ICOUNT
      ACC(ICOUNT)=DATA(3)
C
C     AT THE END OF THE RUN, WRITE OUT THE STORED ACCELERATION
C
C
      IEND=PAR(2)
      IF(ICOUNT.EQ.IEND) THEN
         AFILE='                              '
         WRITE(*,'(A)')'Enter the file name of the spectrum file'
         read(*,'(A)') AFILE
         OPEN(32,FILE=AFILE)
         WRITE(32,6205) (ACC(I),I=1,IEND)
      ENDIF
 6205 FORMAT(6F12.1)
      RETURN
C
 6300 IF(IDR.GT.1600) GO TO 6400
C
C     FIND TURNING RADIUS
C
      IMR=PAR(1)
      JMR=PAR(2)
      I1=PAR(3)
      I2=PAR(4)
      I3=PAR(5)
      CALL INFO('DISP',IMR,JMR,JMR,DATA,ERRFLG)
      RESULT(2)=DATA(1)
      RESULT(3)=DATA(2)
      IF(JCOUNT.EQ.I1)THEN
         X1C=DATA(1)
         Y1C=DATA(2)
      ENDIF
      IF(JCOUNT.EQ.I2)THEN
         X2C=DATA(1)
         Y2C=DATA(2)
      ENDIF
      IF(JCOUNT.EQ.I3)THEN
         X3C=DATA(1)
         Y3C=DATA(2)
         CALL CIRCL(X1C,Y1C,X2C,Y2C,X3C,Y3C,XC,YC,RADIUS)
         RESULT(6)=XC
         RESULT(7)=YC
         RESULT(8)=RADIUS
      ENDIF
      JCOUNT=JCOUNT+1
      RETURN
 6400 IF(IDR.GT.1700) GO TO 240
C
C     SKIP ON THE FIRST CALL TO REQSUB
C
      IF(IFPSS.EQ.0) THEN
         IFPSS=1
         RETURN
      ENDIF
C
C     FIRST FIND THE DISPLACEMENTS
C
      TIMEX(ICOUNT)=TIME
      IMKR=PAR(1)
      JMKR=PAR(2)
      CALL INFO('DISP',IMKR,JMKR,JMKR,DATA,ERRFLG)
      HITCHX(ICOUNT)=DATA(1)
      HITCHY(ICOUNT)=DATA(2)
      HITCHZ(ICOUNT)=DATA(3)
      IMKR=PAR(3)
      JMKR=PAR(4)
      CALL INFO('DISP',IMKR,JMKR,JMKR,DATA,ERRFLG)
      ROX(ICOUNT)=DATA(1)
      ROY(ICOUNT)=DATA(2)
      ROZ(ICOUNT)=DATA(3)
      IMKR=PAR(5)
      JMKR=PAR(6)
      CALL INFO('DISP',IMKR,JMKR,JMKR,DATA,ERRFLG)
      RIX(ICOUNT)=DATA(1)
      RIY(ICOUNT)=DATA(2)
      RIZ(ICOUNT)=DATA(3)
      IMKR=PAR(7)
      JMKR=PAR(8)
      CALL INFO('DISP',IMKR,JMKR,JMKR,DATA,ERRFLG)
      ALIX(ICOUNT)=DATA(1)
      ALIY(ICOUNT)=DATA(2)
      ALIZ(ICOUNT)=DATA(3)
      IMKR=PAR(9)
      JMKR=PAR(10)
      CALL INFO('DISP',IMKR,JMKR,JMKR,DATA,ERRFLG)
      ALOX(ICOUNT)=DATA(1)
      ALOY(ICOUNT)=DATA(2)
      ALOZ(ICOUNT)=DATA(3)
C
C     WRITE THESE OUT
C
      RESULT(2)=HITCHY(ICOUNT)
      RESULT(3)=ROY(ICOUNT)
      RESULT(4)=RIY(ICOUNT)
      RESULT(6)=ALIY(ICOUNT)
      RESULT(7)=ALOY(ICOUNT)
C
C
C
      IEND=PAR(11)
 1086 FORMAT(F6.2,F10.5)
      IF(ICOUNT.EQ.IEND) THEN
         AFILE='hitchx.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),HITCHX(K),K=1,IEND)
         CLOSE(35)
         AFILE='hitchy.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),HITCHY(K),K=1,IEND)
         CLOSE(35)
         AFILE='hitchz.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),HITCHZ(K),K=1,IEND)
         CLOSE(35)
         AFILE='rox.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ROX(K),K=1,IEND)
         CLOSE(35)
         AFILE='roy.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ROY(K),K=1,IEND)
         CLOSE(35)
         AFILE='roz.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ROZ(K),K=1,IEND)
         CLOSE(35)
         AFILE='rix.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),RIX(K),K=1,IEND)
         CLOSE(35)
         AFILE='riy.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),RIY(K),K=1,IEND)
         CLOSE(35)
         AFILE='riz.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),RIZ(K),K=1,IEND)
         CLOSE(35)
         AFILE='lix.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ALIX(K),K=1,IEND)
         CLOSE(35)
         AFILE='liy.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ALIY(K),K=1,IEND)
         CLOSE(35)
         AFILE='liz.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ALIZ(K),K=1,IEND)
         CLOSE(35)
         AFILE='lox.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ALOX(K),K=1,IEND)
         CLOSE(35)
         AFILE='loy.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ALOY(K),K=1,IEND)
         CLOSE(35)
         AFILE='loz.dat'
         OPEN(35,FILE=AFILE)
         WRITE(35,1086)(TIMEX(K),ALOZ(K),K=1,IEND)
         CLOSE(35)
      ENDIF
      ICOUNT=ICOUNT+1
      RETURN
C
C     LONGITUDINAL OR LATERAL FORCE IS ZERO IF VERTICAL FORCE IS ZERO
C
 240  CONTINUE
      FMAG = 0.0D0
      RETURN
C
C
C     WRITE ERROR MESSAGE
C
C     1000 WRITE(1,990) TIME
C     990 FORMAT(1X,'WARNING--ERROR IN RSUB AT TIME=',1PE13.5)
C
      END
C
C****************************************************************
C
C
C******************************************************************
C
C
C ********************************************************************
C
C
C****************************************************************
C
      SUBROUTINE CALCTH ( XI, OMEGAT, AMPD2, OMEGA, X, DX, D2X )
C
C   IT IS CALLED BY UGEN.FORT
C
      IMPLICIT REAL*8 ( A-H, O-Z )
      SINWT = DSIN ( OMEGAT )
      COSWT = DCOS ( OMEGAT )
      X = XI + AMPD2 * ( 1.0D0 - COSWT )
      DX = AMPD2 * OMEGA * SINWT
      D2X = AMPD2 * ( OMEGA**2 ) * COSWT
      RETURN
      END
C
C *****************************************************************
C
      SUBROUTINE VBM(P,BM,PRES,BULKM,NPTS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PRES(*),BULKM(*)
      NPM1=NPTS-1

      DO 2 I=1,NPM1
      IF(P.GE.PRES(I).AND.P.LE.PRES(I+1))ICOUNT=I+1
      IF(P.GE.PRES(I).AND.P.LE.PRES(I+1))GO TO 3
    2 CONTINUE

      IF(P.LE.PRES(1))BM=BULKM(1)
      IF(P.GE.PRES(NPTS))BM=BULKM(NPTS)
      RETURN
C
    3 SLOPE=(P-PRES(ICOUNT-1))/(PRES(ICOUNT)-PRES(ICOUNT-1))
      BM=BULKM(ICOUNT-1)+SLOPE*(BULKM(ICOUNT)-BULKM(ICOUNT-1))
      RETURN
      END
C
      SUBROUTINE VBM2(J,P,BM,PRES,BULKM,NPTS)
      use jdt_max
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PRES(maxCPT,*),BULKM(maxCPT,*),NPTS(*)
      NPM1=NPTS(J)-1

      DO 2 I=1,NPM1
      IF(P.GE.PRES(I,J).AND.P.LE.PRES(I+1,J))ICOUNT=I+1
      IF(P.GE.PRES(I,J).AND.P.LE.PRES(I+1,J))GO TO 3
    2 CONTINUE

      IF(P.LE.PRES(1,J))BM=BULKM(1,J)
      IF(P.GE.PRES(NPTS(J),J))BM=BULKM(NPTS(J),J)
      RETURN
C
    3 SLOPE=(P-PRES(ICOUNT-1,J))/(PRES(ICOUNT,J)-PRES(ICOUNT-1,J))
      BM=BULKM(ICOUNT-1,J)+SLOPE*(BULKM(ICOUNT,J)-BULKM(ICOUNT-1,J))
      RETURN
      END
C

      SUBROUTINE TIREF(J)
C     TIREF USES A DISTRIBUTED CONTACT TIRE MODEL TO PREDICT THE FORCE
C     ACTING ON A TIRE AS IT TRAVELS OVER ROUGH GROUND
      use jdt_var
      IMPLICIT REAL*8 (A-H,O-Z)
C     FIND THE NUMBER, ITR, OF THE TERRAIN STATION IMMEDIATELY AHEAD OF THE
C     CENTER OF WHEEL J
      ITRM1=ITR(J)-3
      IF(ITR(J).LE.3) ITRM1=1
      NSTAM1=NTSTA(J)-1
      DO 5410 I=ITRM1,NSTAM1
         IF(X1(J).GE.STA(I,J).AND.X1(J).LE.STA(I+1,J)) ITR(J)=I+1
         IF(X1(J).GE.STA(I,J).AND.X1(J).LE.STA(I+1,J)) GO TO 5411
 5410 CONTINUE
C     CALCULATE THE LONGITUDINAL AND VERTICAL TIRE SPRING FORCE COMPONENTS
C     THE SEGMENTED TIRE MODEL
 5411 FX(J)=0.0
      FZ(J)=0.0
      NSEG=ANSEG(J)
      NSEGD2=NSEG/2
C     CALCULATE THE LENGTH, TMAG(ISEG), OF EACH SEGMENT SPRING. IF THE SEGMENT
C     SPRING IS NOT COMPRESSED BY THE TERRAIN, THIS LENGTH IS SET EQUAL TO THE
C     UNDEFLECTED WHEEL RADIUS. IF THE SEGMENT SPRING IS COMPRESSED, THE
C     LENGTH IS THE DISTANCE BETWEEN THE WHEEL CENTER AND THE POINT OF INTER-
C     SECTION OF THE SPRING WITH THE TERRAIN.
C     FIRST PERFORM THIS DETERMINATION FOR THE SEGMENT SPRINGS LYING BEHIND THE
C     WHEEL CENTER
      KLAST=0
      DO 4001 I=1,NSEGD2
         ISEG=NSEGD2+1-I
         DO 4000 K=1,500
            XTEST1=X1(J)-STA(ITR(J)-KLAST-K,J)
            IF(DABS(XTEST1).LE.0.01) GO TO 4000
            SP=(Y1(J)-HGT(ITR(J)-KLAST-K,J))/XTEST1
            IF(SP.GT.SM(ISEG,J)) GO TO 4000
            X2=STA(ITR(J)-KLAST-K,J)
            Y2=HGT(ITR(J)-KLAST-K,J)
            X3=STA(ITR(J)-KLAST-K+1,J)
            Y3=HGT(ITR(J)-KLAST-K+1,J)
            ST=(Y3-Y2)/(X3-X2)
            XX=(Y2-Y1(J)+SM(ISEG,J)*X1(J)-ST*X2)/(SM(ISEG,J)-ST)
            YY=ST*(XX-X2)+Y2
            DIST=DSQRT((X1(J)-XX)**2+(Y1(J)-YY)**2)
            IF(DIST.GE.WRAD(J)) TMAG(ISEG)=WRAD(J)
            IF(DIST.LT.WRAD(J)) TMAG(ISEG)=DIST
            KLAST=KLAST+K-1
            GO TO 4001
 4000    CONTINUE
 4001 CONTINUE
C     NEXT PERFORM THIS DETERMINATION FOR THE SEGMENT SPRINGS LYING AHEAD OF THE
C     WHEEL CENTER
      NSTART=NSEGD2+1
      KLAST=-1
      DO 4003 I=NSTART,NSEG
         ISEG=I
         DO 4002 K=1,500
            XTEST2=STA(ITR(J)+KLAST+K,J)-X1(J)
            IF(DABS(XTEST2).LE.0.01) GO TO 4002
            SP=(Y1(J)-HGT(ITR(J)+KLAST+K,J))/XTEST2
            IF(SP.GT.DABS(SM(ISEG,J))) GO TO 4002
            X2=STA(ITR(J)+KLAST+K-1,J)
            Y2=HGT(ITR(J)+KLAST+K-1,J)
            X3=STA(ITR(J)+KLAST+K,J)
            Y3=HGT(ITR(J)+KLAST+K,J)
            ST=(Y3-Y2)/(X3-X2)
            XX=(Y2-Y1(J)+SM(ISEG,J)*X1(J)-ST*X2)/(SM(ISEG,J)-ST)
            YY=ST*(XX-X2)+Y2
            DIST=DSQRT((X1(J)-XX)**2+(Y1(J)-YY)**2)
            IF(DIST.GE.WRAD(J)) TMAG(ISEG)=WRAD(J)
            IF(DIST.LT.WRAD(J)) TMAG(ISEG)=DIST
            KLAST=KLAST+K-1
            GO TO 4003
 4002    CONTINUE
 4003 CONTINUE
C     CALCULATE THE SPRING FORCE ACTING IN EACH SEGMENT
      DAREA=0.0
      DO 5358 I=1,NSEG
         TMAG(I)=WRAD(J)-TMAG(I)
         IF(LTIRE.EQ.0) TMAG(I)=TMAG(I)*ASEGK(I,J)
         IF(LTIRE.EQ.1) DAREA=DAREA+(WRAD(J)*TMAG(I)-0.5*TMAG(I)**2)
     1                        *DTHET(J)
         IF(LTIRE.EQ.1) TMAG(I)=
     1                     TABLE1(DEFL,FORCE,TMAG(I),2,NLDPT(J),1,J)
 5358 CONTINUE
C     SUM THE LONGITUDINAL AND VERTICAL COMPONENTS OF EACH SEGMENT SPRING FORCE
C     TO FIND THE NET LONGITUDINAL AND VERTICAL TIRE SPRING FORCE COMPONENTS
      DO 5359 I=1,NSEG
         FX(J)=FX(J)-TMAG(I)*SSEGA(I,J)
         FZ(J)=FZ(J)+TMAG(I)*CSEGA(I,J)
 5359 CONTINUE
C     IF THE INTERPOLATION METHOD IS TO BE USED, FIND THE NET LONGITUDINAL
C     AND VERTICAL TIRE SPRING FORCE COMPONENTS
      IF(LTIRE.EQ.0) GO TO 5360
      IF(FZ(J).LE.0.0) GO TO 5360
      ABSFZ=DABS(FZ(J))
      ALPH=DATAN2(FX(J),ABSFZ)
      FMAG=TABLE1(AREA,FORCE,DAREA,2,NLDPT(J),1,J)
      FX(J)=FMAG*DSIN(ALPH)
      FZ(J)=FMAG*DCOS(ALPH)
 5360 CONTINUE
C     CHECK TO SEE IF THE VERTICAL TIRE SPRING FORCE IS NEGATIVE INDICATING
C     THE TIRE HAS LOST CONTACT WITH THE GROUND. IF SO, SET ALL THE TIRE FORCE
C     COMPONENTS TO ZERO.
      IF(FZ(J).LE.0.0) THEN
         FX(J)=0.0
         FZ(J)=0.0
      ENDIF
      RETURN
      END


      FUNCTION TABLE1(X,Y,XARG,IDEG,NDIM,JMIN,K)
C
C  THIS FUNCTION USES LAGRANGIAN INTERPOLATION OF DEGREE IDEG TO
C         DETERMINE A VALUE OF Y CORRESPONDING TO XARG, BASED UPON
C         THE NDIM PAIRED VALUES OF X AND Y.
C         INTERPOLATION USES VALUES SUBSCRIPTED.GE.JMIN.
C  EXAMPLE USE:
C         TWO ARRAYS, TEMP(I) AND HUMID(I), IN SOME PROGRAM CONTAIN
C             CONSECUTIVE VALUES OF TWO VARIABLES. THERE ARE NDIM
C             OF EACH VARIABLE, WITH THOSE HAVING EQUAL SUBSCRIPTS
C             CORRESPONDING TO ONE ANOTHER. TO OBTAIN THE INTERPOLATED
C             VALUE, HVALUE, CORRESPONDING TO A SPECIFIC VALUE WITHIN
C             THE TEMP(I) RANGE, SAY TVALUE, THE FOLLOWING ASSIGNMENT
C             STATEMENT WOULD BE USED IN THE PROGRAM WHICH HAS THE TWO
C             ARRAYS STORED:
C                  HVALUE = TABLE(TEMP,HUMID,TVALUE,IDEG,NDIM,JMIN)
C
C  *** CAUTION ***   X AND Y MUST BE DIMENSIONED NDIM IN MAIN PROGRAM
C
      use jdt_max
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(maxCPT,*),Y(maxCPT,*)
C  SEARCH FOR THE INTERVAL CONTAINING XARG.
#ifdef FT_DEBUG
      IF(XARG.GT.X(NDIM,K)) WRITE(17,40)
   40 FORMAT(5X,'*** DEFLECTION (OR AREA) IS BEYOND THE RANGE FOR WHICH
     1DATA IS AVAILABLE')
#endif
      IF(IDEG.LE.NDIM-JMIN) GO TO 9
      WRITE(17,1)
    1 FORMAT(5X,'*** INSUFFICIENT NUMBER OF TABULATED POINTS FOR DESIRED
     1 DEGREE OF INTERPOLATION')
    9 J=JMIN-1
   10 IF(J.EQ.NDIM) GO TO 300
      J=J+1
      IF(X(J,K).EQ.XARG) GO TO 20
      IF(J.GE.NDIM-1) GO TO 300
      IF(X(J,K).LT.XARG.AND.X(J+1,K).GE.XARG) GO TO 30
      IF(X(J,K).GT.XARG.AND.X(J+1,K).LE.XARG) GO TO 30
      GO TO 10
   20 TABLE1=Y(J,K)
      RETURN
   30 JHALF=(IDEG+1)/2
      IF(JHALF.GE.J) GO TO 31
      J1=J-JHALF
      IF(J1+IDEG.LE.NDIM) GO TO 32
  300 J1=NDIM-IDEG
      GO TO 32
   31 J1=JMIN
   32 JMAX=J1+IDEG
C  START THE LAGRANGIAN INTERPOLATION
      FACTOR=1.0
      DO 33 J=J1,JMAX
      IF(X(J,K).EQ.XARG) GO TO 20
   33 FACTOR=FACTOR*(XARG-X(J,K))
      YEST=0.0
      DO 35 I=J1,JMAX
      TERM=Y(I,K)*FACTOR/(XARG-X(I,K))
      DO 34 J=J1,JMAX
   34 IF(I.NE.J) TERM=TERM/(X(I,K)-X(J,K))
   35 YEST=YEST+TERM
      TABLE1=YEST
      RETURN
      END
C
C
C
      SUBROUTINE CIRCL(X1,Y1,X2,Y2,X3,Y3,XC,YC,RADIUS)
C     THIS SUBROUTINE COMPUTES THE CENTER OF A CIRCLE XC,YC
C     GIVEN THREE POINTS ON THE CIRCLE X1,Y1  X2,Y2  X3,Y3
      IMPLICIT REAL*8(A-H,O-Z)
      UX1=X2-X1
      UY1=Y2-Y1
      DEN=SQRT(UX1*UX1+UY1*UY1)
      UX1=UX1/DEN
      UY1=UY1/DEN
      D1=(-UX1*(X2+X1)-UY1*(Y2+Y1))/2.0
      UX2=X3-X2
      UY2=Y3-Y2
      DEN=SQRT(UX2*UX2+UY2*UY2)
      UX2=UX2/DEN
      UY2=UY2/DEN
      D2=(-UX2*(X3+X2)-UY2*(Y3+Y2))/2.0
      DEN=UX1*UY2-UX2*UY1
      XC=(-D1*UY2+D2*UY1)/DEN
      YC=(-D2*UX1+D1*UX2)/DEN
      RADIUS=SQRT((X1-XC)**2+(Y1-YC)**2)
      RETURN
      END


      SUBROUTINE FNCDEP(SUB,TYPE,ID,VECTOR,NV,ERRFLG)
      DOUBLE PRECISION VECTOR(NV)
      DOUBLE PRECISION DATA(6)
      INTEGER          ID
      INTEGER          NV
      INTEGER          IPAR
      INTEGER          I
      INTEGER          NDUMMY
      LOGICAL          ERRFLG
      CHARACTER*(*) SUB
      CHARACTER*(*) TYPE
C-----------------------------------------------------------------------
      DO 10 I=1,NV
         IPAR = VECTOR(I)
         CALL SYSARY('DISP',IPAR,1,DATA, NDUMMY, ERRFLG)
         CALL SYSARY('VEL ',IPAR,1,DATA, NDUMMY, ERRFLG)
  10  CONTINUE
      RETURN
      END
