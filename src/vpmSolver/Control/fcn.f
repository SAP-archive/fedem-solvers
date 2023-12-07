C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
C **********************************************************************

C         *****    ****   *   *
C         *       *       **  *
C         ****    *       * * *
C         *       *       *  **
C         *        ****   *   *

C **********************************************************************

C  Name         : FCN

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Compute functional values for a specified system

C  Method       : Each module in the control system library is given a number
C                 which refers to a specific subroutine. When information on
C                 module types and connections is provided, the actual modules
C                 will be called.

C  Call         : CALL FCN(IOP,T,TNM,NREG,NMOD,NTD,VREG,VDREG,VPREG,DELAY,
C                >         TOL,MSTAT,MMOD,MPMTOP,MMTOP,MPRPAR,RPAR,IERR,
C                >         LIM,NFCN,ZCTRL0)

C  Arguments    : IOP    - integer, operation code
C                        = 1 : initialization, determination of MSTAT
C                       ^= 1 : function evaluation
C                 T      - real, simulation time
C                 TNM    - real, maximun end point for next integration step
C                 NREG   - integer, no. of variables
C                 NMOD   - integer, no. of modules
C                 NTD    - integer, no. of time delays
C                 VREG   - real array(NREG), the most recent values of the
C                          control variables
C                 VDREG  - real array(NREG), function values
C                 VPREG  - real array(NREG), previous values of variables
C                 DELAY  - real array(*), time history of delayed variables
C                 TOL    - real array(NREG), local error tolerance
C                 MSTAT  - integer array(NREG), status flags for variables
C                 MMOD   - integer array(NREG), module number where each
C                          variable is calculated
C                 MPMTOP - integer array(NMOD+1), pointers to the content of
C                          MMTOP
C                 MMTOP  - integer array(MPMTOP(NMOD+1)), topology vector
C                 MPRPAR - integer array(NMOD+1), pointers to the content of
C                          RPAR
C                 RPAR   - real array(MPRPAR(NMOD+1)), module parameters
C                 IERR   - integer array(2), error flag
C                 IERR(1)= -2 : Could not solve function
C                 IERR(1)= -1 : non-existent module type
C                        =  0 : no error
C                 IERR(2)= module type number
C                 LIM    - integer array(2), limit crossing indicator
C                 SIGN(LIM(1)) - indicates crossing direction
C                 ABS(LIM(1)) - parameter number for limit value
C                 LIM(2) - crossing variable number
C                 NFCN   - integer, no. of function evaluations

C  Input        : IOP,T,TNM,NREG,NMOD,NTD,VREG,VPREG,DELAY,
C                 MSTAT (when IOP ^= 1), MPMTOP,MMTOP,MPRPAR,RPAR,NFCN

C  Output       : TNM,MSTAT & MMOD (when IOP = 1), VDREG, IERR, LIM, NFCN

C  Work area    : None

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 MODxx - for function evaluation of module no. xx

C  Fault react. : Return with a non-zero error flag

C  Note         : Nothing

C **********************************************************************

      SUBROUTINE FCN(IOP, T, TNM, NREG, NMOD, NTD, VREG, VDREG, VPREG,
     >               DELAY, TOL, MSTAT, MMOD, MPMTOP, MMTOP, MPRPAR,
     >               RPAR, IERR, LIM, NFCN, ZCTRL0)

      IMPLICIT NONE

      INTEGER          IOP, NREG, NMOD, NTD
      DOUBLE PRECISION T, TNM, VREG(NREG), VDREG(NREG), VPREG(NREG),
     >                 DELAY(*), TOL(NREG), RPAR(*)
      INTEGER          MSTAT(NREG), MMOD(NREG), MPMTOP(NMOD+1),
     >                 MPRPAR(NMOD+1), MMTOP(*), IERR(2), NFCN, LIM(2)
      LOGICAL          ZCTRL0

      INTEGER          ITYPE, MF, MODLIM, MODNO, MPAR, MPLUS, MTYP,
     >                 N, NALG, NDIM, NDLAY, NINPUT, NN, NOUTPT, NSTATE
      DOUBLE PRECISION TAU, TAUMOD, TMF

C ............................. Set tolerance multiplication factor for
C ..............................piecewise linear modules
      PARAMETER ( TMF = 0.001D0 )
C ............................. Initialize time delay counter
      NDLAY = 0
C ---> 02.02.10/TI ............ Initialise time of earliest discontinuity
      TAU = 1.0d0
C ............................. For each module ...
      DO 999 MODNO = 1, NMOD
C ............................. find the pointer to the first node ...
         MF   = MPMTOP(MODNO)
C ............................. find the pointer to the first parameter ...
         MPAR = MPRPAR(MODNO) + 1
C ............................. which follows type of module ...
         MTYP = IDINT(RPAR(MPAR-1))
C ............................. and move to that type of module
         IF (MTYP .LT. 1 .OR. MTYP .GT. 110) GOTO 211
         GOTO (101,102,103,104,105,106,107,107,107,107,
     >         111,112,113,113,113,113,113,113,113,113,
     >         121,122,123,124,125,125,125,125,125,125,
     >         131,132,133,134,135,136,137,138,138,138,
     >         141,142,143,144,145,146,147,147,147,147,
     >         147,147,147,147,147,147,147,147,147,147,
     >         147,147,147,147,147,147,147,147,147,147,
     >         147,147,147,147,147,147,147,147,147,147,
     >         147,147,147,147,147,147,147,147,147,147,
     >         147,147,147,147,147,147,147,147,147,147,
     >         2000,2000,2000,2000,2000,2000,2000,2000,2000,2000) MTYP

C ............................. Module type 1 - Comparator
  101    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD1(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),  VREG(MMTOP(MF+2)),
     >        VDREG(MMTOP(MF)), VDREG(MMTOP(MF+1)), VDREG(MMTOP(MF+2)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)), MSTAT(MMTOP(MF+2)),
     >        MMOD(MMTOP(MF)),  MMOD(MMTOP(MF+1)),  MMOD(MMTOP(MF+2)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 2 - Adder
  102    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD2(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),  VREG(MMTOP(MF+2)),
     >        VDREG(MMTOP(MF)), VDREG(MMTOP(MF+1)), VDREG(MMTOP(MF+2)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)), MSTAT(MMTOP(MF+2)),
     >        MMOD(MMTOP(MF)),  MMOD(MMTOP(MF+1)),  MMOD(MMTOP(MF+2)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 3 - Amplifier
  103    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD3(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 4 - Integrator
  104    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD4(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >                          MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR))
         GOTO 999

C ............................. Module type 5 - Derivator
  105    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD5(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),  VREG(MMTOP(MF+2)),
     >                          VDREG(MMTOP(MF+1)), VDREG(MMTOP(MF+2)),
     >                          MSTAT(MMTOP(MF+1)), MSTAT(MMTOP(MF+2)),
     >                          MMOD(MMTOP(MF+1)),  MMOD(MMTOP(MF+2)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 6 - Multiplier
  106    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD6(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),  VREG(MMTOP(MF+2)),
     >        VDREG(MMTOP(MF)), VDREG(MMTOP(MF+1)), VDREG(MMTOP(MF+2)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)), MSTAT(MMTOP(MF+2)),
     >        MMOD(MMTOP(MF)),  MMOD(MMTOP(MF+1)),  MMOD(MMTOP(MF+2)),
     >        ZCTRL0)
         GOTO 999

C ............................. Module type 7 - Power function
  107    CONTINUE

C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD7(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 8 - 10
         IERR(1) = -1
         IERR(2) = MTYP
         RETURN

C ............................. Module type 11 - Time delay
  111    CONTINUE
         NDLAY = NDLAY + 1
         CALL MOD11(IOP, MODNO, T, NDLAY, NTD, DELAY,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 12 - Sample and hold
  112    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD12(IOP, MODNO, T, TNM, VPREG(MMTOP(MF)),
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 13 - 20
  113    CONTINUE
         IERR(1) = -1
         IERR(2) = MTYP
         RETURN

C ............................. Module type 21 - Logical switch
  121    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         TAUMOD = 2.0d0*TAU
         CALL MOD21(IOP, MODNO, MODLIM, TAUMOD,
     >        TMF*TOL(MMTOP(MF)), VPREG(MMTOP(MF)),
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         IF ( MODLIM .NE. 0 ) THEN
            IF ( TAUMOD .LT. TAU ) THEN
               LIM(1) = SIGN(MPAR + IABS(MODLIM) - 1, MODLIM)
               LIM(2) = MMTOP(MF)
               TAU = TAUMOD
            END IF
         END IF
         GOTO 999

C ............................. Module type 22 - Limitation
  122    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         TAUMOD = 2.0d0*TAU
         CALL MOD22(IOP, MODNO, MODLIM, TAUMOD,
     >        TMF*TOL(MMTOP(MF)), VPREG(MMTOP(MF)),
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         IF ( MODLIM .NE. 0 ) THEN
            IF ( TAUMOD .LT. TAU ) THEN
               LIM(1) = SIGN(MPAR + IABS(MODLIM) - 1, MODLIM)
               LIM(2) = MMTOP(MF)
               TAU = TAUMOD
            END IF
         END IF
         GOTO 999

C ............................. Module type 23 - Dead zone
  123    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         TAUMOD = 2.0d0*TAU
         CALL MOD23(IOP, MODNO, MODLIM, TAUMOD,
     >        TMF*TOL(MMTOP(MF)), VPREG(MMTOP(MF)),
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         IF ( MODLIM .NE. 0 ) THEN
            IF ( TAUMOD .LT. TAU ) THEN
               LIM(1) = SIGN(MPAR + IABS(MODLIM) - 1, MODLIM)
               LIM(2) = MMTOP(MF)
               TAU = TAUMOD
            END IF
         END IF
         GOTO 999

C ............................. Module type 24 - Hysteresis
  124    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD24(IOP, MODNO, T, TNM, TMF*TOL(MMTOP(MF)),
     >        VPREG(MMTOP(MF)), VPREG(MMTOP(MF+1)),
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >        MSTAT(MMTOP(MF)), MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 25 - 30
  125    CONTINUE
         IERR(1) = -1
         IERR(2) = MTYP
         RETURN

C ............................. Module type 31 - PI-controller
  131    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD31(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),  VREG(MMTOP(MF+2)),
     >                          VDREG(MMTOP(MF+1)), VDREG(MMTOP(MF+2)),
     >                          MSTAT(MMTOP(MF+1)), MSTAT(MMTOP(MF+2)),
     >                          MMOD(MMTOP(MF+1)),  MMOD(MMTOP(MF+2)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 32 - Propotional + limited
C ............................. Integral controller
  132    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD32(IOP, MODNO,  VREG(MMTOP(MF)),
     >        VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)), VREG(MMTOP(MF+3)),
     >        VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),VDREG(MMTOP(MF+3)),
     >        MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),MSTAT(MMTOP(MF+3)),
     >        MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)), MMOD(MMTOP(MF+3)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 33 - PD-controller
  133    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD33(IOP, MODNO, VREG(MMTOP(MF)),
     >        VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)), VREG(MMTOP(MF+3)),
     >        VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),VDREG(MMTOP(MF+3)),
     >        MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),MSTAT(MMTOP(MF+3)),
     >        MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)), MMOD(MMTOP(MF+3)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 34 - Proportional + limited
C ............................. Derivative controller
  134    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD34(IOP, MODNO, VREG(MMTOP(MF)),
     >        VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)), VREG(MMTOP(MF+3)),
     >        VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),VDREG(MMTOP(MF+3)),
     >        MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),MSTAT(MMTOP(MF+3)),
     >        MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)), MMOD(MMTOP(MF+3)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 35 - PID-controller
  135    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD35(IOP, MODNO,
     >        VREG(MMTOP(MF)),   VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)),
     >        VREG(MMTOP(MF+3)), VREG(MMTOP(MF+4)),
     >                           VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),
     >        VDREG(MMTOP(MF+3)),VDREG(MMTOP(MF+4)),
     >                           MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),
     >        MSTAT(MMTOP(MF+3)),MSTAT(MMTOP(MF+4)),
     >                           MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)),
     >        MMOD(MMTOP(MF+3)), MMOD(MMTOP(MF+4)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 36 - Proportional + Integral
C ............................. + limited Derivative controller
  136    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD36(IOP, MODNO,
     >        VREG(MMTOP(MF)),   VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)),
     >        VREG(MMTOP(MF+3)), VREG(MMTOP(MF+4)),
     >                           VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),
     >        VDREG(MMTOP(MF+3)),VDREG(MMTOP(MF+4)),
     >                           MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),
     >        MSTAT(MMTOP(MF+3)),MSTAT(MMTOP(MF+4)),
     >                           MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)),
     >        MMOD(MMTOP(MF+3)), MMOD(MMTOP(MF+4)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 37 - Proportional + limited
C ............................. Integral + limited Derivative controller
  137    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD37(IOP, MODNO,
     >        VREG(MMTOP(MF)),   VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)),
     >        VREG(MMTOP(MF+3)), VREG(MMTOP(MF+4)),
     >                           VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),
     >        VDREG(MMTOP(MF+3)),VDREG(MMTOP(MF+4)),
     >                           MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),
     >        MSTAT(MMTOP(MF+3)),MSTAT(MMTOP(MF+4)),
     >                           MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)),
     >        MMOD(MMTOP(MF+3)), MMOD(MMTOP(MF+4)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 38 - 40
  138    CONTINUE
         IERR(1) = -1
         IERR(2) = MTYP
         RETURN

C ............................. Module type 41 - Real pole
  141    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD41(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),
     >                          VDREG(MMTOP(MF+1)),
     >                          MSTAT(MMTOP(MF+1)),
     >                          MMOD(MMTOP(MF+1)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 42 - Complex conjugate poles
  142    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD42(IOP, MODNO,
     >        VREG(MMTOP(MF)),  VREG(MMTOP(MF+1)),  VREG(MMTOP(MF+2)),
     >                          VDREG(MMTOP(MF+1)), VDREG(MMTOP(MF+2)),
     >                          MSTAT(MMTOP(MF+1)), MSTAT(MMTOP(MF+2)),
     >                          MMOD(MMTOP(MF+1)),  MMOD(MMTOP(MF+2)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 43 - First order element
  143    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD43(IOP, MODNO,  VREG(MMTOP(MF)),
     >        VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)), VREG(MMTOP(MF+3)),
     >        VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),VDREG(MMTOP(MF+3)),
     >        MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),MSTAT(MMTOP(MF+3)),
     >        MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)), MMOD(MMTOP(MF+3)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 44 - Second order element
  144    CONTINUE
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD44(IOP, MODNO,
     >        VREG(MMTOP(MF)),   VREG(MMTOP(MF+1)), VREG(MMTOP(MF+2)),
     >        VREG(MMTOP(MF+3)), VREG(MMTOP(MF+4)), VREG(MMTOP(MF+5)),
     >                           VDREG(MMTOP(MF+1)),VDREG(MMTOP(MF+2)),
     >        VDREG(MMTOP(MF+3)),VDREG(MMTOP(MF+4)),VDREG(MMTOP(MF+5)),
     >                           MSTAT(MMTOP(MF+1)),MSTAT(MMTOP(MF+2)),
     >        MSTAT(MMTOP(MF+3)),MSTAT(MMTOP(MF+4)),MSTAT(MMTOP(MF+5)),
     >                           MMOD(MMTOP(MF+1)), MMOD(MMTOP(MF+2)),
     >        MMOD(MMTOP(MF+3)), MMOD(MMTOP(MF+4)), MMOD(MMTOP(MF+5)),
     >        RPAR(MPAR), ZCTRL0)
         GOTO 999

C ............................. Module type 45 - Linear control system
  145    CONTINUE
C ............................. Determine the number of inputs, state vari-
C ............................. ables and outputs
         ITYPE = IDINT(RPAR(MPAR))
         IF (ITYPE .EQ. 0) THEN
            NINPUT = IDINT(RPAR(MPAR+1))
            NSTATE = IDINT(RPAR(MPAR+2))
            NOUTPT = IDINT(RPAR(MPAR+3))
            MPLUS  = 4
         ELSE
            NINPUT = 1
            NSTATE = IDINT(RPAR(MPAR+1))
            NOUTPT = 1
            MPLUS  = 2
         END IF

C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD45(IOP, MODNO, VREG, VDREG, MSTAT, MMOD, MMTOP, MF,
     >        NINPUT, NSTATE, NOUTPT, RPAR(MPAR+MPLUS), ITYPE)
         GOTO 999

C ............................. Module type 46 - General transfer function
  146    CONTINUE
C ............................. The first parameter is the order. The leading
C ............................. coefficient in the denominator is 1.
         NDIM = IDINT(RPAR(MPAR))
C ............................. Invoke the module for initialization ...
C ............................. or computation of functional values
         CALL MOD46(IOP, MODNO, VREG, VDREG, MSTAT, MMOD, MMTOP, MF,
     >        NDIM, RPAR(MPAR+1))
         GOTO 999

C ............................. Module type 47 - 100
  147    CONTINUE
         IERR(1) = -1
         IERR(2) = MTYP
         RETURN

C ............................. Module type 101 - 110 - User defined models
 2000    CONTINUE
C ............................. In case of initialization, ...
         IF (IOP .EQ. 1) THEN
C ............................. Read the number of inputs, state variables and
C ............................. algebraic variables
            NINPUT = IDINT(RPAR(MPAR))
            NSTATE = IDINT(RPAR(MPAR+1))
            NALG   = IDINT(RPAR(MPAR+2))
C ............................. Initialize the status mstat and mmod
            NN = MF + NINPUT - 1
            DO 2001 N = 1,NSTATE
               MSTAT(MMTOP(NN+N)) = 2
               MMOD(MMTOP(NN+N)) = MODNO
 2001       CONTINUE
            NN = NN + NSTATE
            DO 2002 N = 1,NALG
               MSTAT(MMTOP(NN+N)) = 3
               MMOD(MMTOP(NN+N)) = MODNO
 2002       CONTINUE
         ELSE
C ............................. When the mode is set to execution, call
C ............................. the routine containing the user model
Ckmo         GOTO (201,202,203,204,205,206,207,208,209,210) MTYP-100
            goto 211
         END IF
         GOTO 999
Ckmo User-defined models currently removed
C 201  CALL MOD101(VREG, VDREG, MMTOP, MF, RPAR(MPAR))
C      GOTO 999
C      ...
C 210  CALL MOD110(VREG, VDREG, MMTOP, MF, RPAR(MPAR))
C      GOTO 999
C ............................. Module type > 110
  211 CONTINUE
      IERR(1) = -1
      IERR(2) = MTYP
      RETURN

  999 CONTINUE

C ............................. Update the number of function calls
      NFCN = NFCN + 1

      RETURN
      END
