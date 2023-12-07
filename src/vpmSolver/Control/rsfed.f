C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C

C **********************************************************************
C
C  Files included:
C
C  rsfed.f
C  rsinit.f
C  sstate.f
C  oneda.f
C  newton.f
C  equil.f
C  calfcn.f
C  beuler.f
C  lobatt.f
C  dlay.f
C  outmes.f

C **********************************************************************

C               ****     ****   *****   *****   ***
C               *   *   *       *       *       *   *
C               ****     ***    ****    ****    *   *
C               * *         *   *       *       *   *
C               *   *   ****    *       *****   ***

C **********************************************************************

C  Name         : RSFED

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Compute the state of a control system in loop with a mecha-
C                 nical system, FEDEM. The two systems are integrated numeri-
C                 cally using different methods. Steady state computation is
C                 also included.

C  Method       : RSFED is called from FEDEM for initialization/ steady state
C                 computation or once per iteration in FEDEM to receive the
C                 equilibrium state of the control system. The state of mecha-
C                 nical variables and external inputs are provided from the
C                 calling routine. Algebraic and state variables are integrat-
C                 ed using an imbedded pair of numerical methods, which allows
C                 for estimation of the truncation error.

C  Call         : CALL RSFED(IOP,T,DT,NREG,NMOD,VREG,MPIREG,IREG,MPRREG,RREG,
C                >           DELAY,NDDIM,IERR,ISTAT,
C                >           ZCTRL0,DCSPLIT,tolFactCtrl,iprint)
C
C  Arguments    : IOP    - integer, operation code
C                        = 1 : initialize the control system
C                        = 2 : find steady state
C                        = 3 : initialize for a new time step and integrate
C                              one step forward
C                        = 4 : integrate one step forward after time step
C                              reduction
C                        = 5 : integrate one step forward
C                 T      - real, end point of the integration step
C                 DT     - real, actual time step
C                 NREG   - integer, no. of variables
C                 NMOD   - integer, no. of modules
C                 VREG   - real array(NREG), the most recent accepted values
C                          of the control state variables
C                 MPIREG - integer array(8), pointers to the content of IREG
C                 IREG   - integer array(12 +5*NREG +2*NMOD +MPMTOP(NMOD+1)),
C                          working area for integer variables,
C                        = ICNT(10) - integer parameters
C                             IOP - operation code
C                             NMOD - no. of modules
C                             LIM(1) - limit crossing indicator
C                             LIM(2) - parameter number showing limit value
C                             NTD - number of time delays
C                             NEWDT - changed time step indicator
C ---> 92.04.24/TI            VARDT - variable step integration flag
C ---> 02.02.09/TI            DISCON - last parameter no. in limit crossing
C                          MSTAT(NREG) - status flags for variables
C                          MMOD(NREG) - module where each variable is computed
C                          IPVT1(NREG) - pivot indices for BNEW
C                          IPVT2(2*NREG) - pivot indices for LNEW
C                          MPMTOP(NMOD+1) - pointers to the content of MMTOP
C                          MPRPAR(NMOD+1) - pointers to the content of RPAR
C                          MMTOP(*) - topology vector
C                 MPRREG - integer array(12), pointers to the content of RREG
C                 RREG   - real array(10 + NREG*(12+6*NREG) + MPRPAR(NMOD+1)),
C                          working area for real variables,
C                        = RCNT(10) - real parameters
C                             T     - end point of the integration step
C                             TNM   - max. end point of next integration step
C                             DT    - previous time step
C                             AEPS  - absolute error tolerance
C                             REPS  - relative error tolerance
C                             ACCU  - accuracy parameter for iteration
C                             PERT  - relative perturbation for computation of
C                                     numerical Jacobian
C                             TRUNC - truncation error
C                          VPREG(NREG) - previous values of control variables
C                          VBREG(NREG) - numerical solution using method 1
C                          VLREG(2*NREG) - numerical solution using method 2
C                          VDREG(NREG) - function values
C                          VHREG(4*NREG) - help vector (increments, old values)
C                          TOL(2*NREG) - local error tolerance
C                          ERR(NREG) - estimated truncation error
C                          AJAC(NREGxNREG) - numerical Jacobian
C                          ABE(NREGxNREG) - iteration matrix for method 1
C                          ALO(2*NREGx2*NREG) - iteration matrix for method 2
C                          RPAR(*) - module parameters
C                 DELAY  - real array(*), time history of delay variables
C                 NDDIM  - integer, delay array dimension
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning
C                 ISTAT  - integer array(10), statistics vector
C                        1 : function evaluations for initialization
C                        2 : function eval. for steady state and integration
C                        3 : Jacobian evaluations
C                        4 : LU-factorizations for steady state
C                        5 : LU-factorizations for Backward Euler
C                        6 : LU-factorizations for Lobatto IIIC
C                        7 : iterations for steady state
C                        8 : iterations in Backward Euler
C                        9 : iterations in Lobatto IIIC
C                       10 : times solving a linear system

C  Input        : IOP,T,DT,NREG,NMOD,VREG,DELAY

C  Output       : VREG,IERR,ISTAT,DELAY

C  Work area    : MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 RSINIT - for initalization
C                 SSTATE - for steady state computation
C                 ONEDA  - for integration one step forward
C                 DLAY   - for updating of the delay buffer

C  Fault react. : Return with non-zero error flag

C  Note         : Nothing

C **********************************************************************

      SUBROUTINE RSFED(JOP, T, DT, NREG, NMOD, VREG, MPIREG, IREG,
     >     MPRREG, RREG, DELAY, NDDIM, IERR, ISTAT,
     >     ZCTRL0, DCSPLIT, tolFactCtrl, iprint)
CDEC$ ATTRIBUTES DLLEXPORT :: RSFED

      use ioModule, only : initIO

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VREG(NREG), RREG(*), DELAY(*)
      INTEGER MPIREG(*), IREG(*), MPRREG(*), IERR(2), ISTAT(10)
      LOGICAL ZCTRL0, DCSPLIT

      IERR(1) = 0
C ---> 02.02.18/TI ............ Initialise time splitting flag and
C ............................. new Jacobian/iteration flag
      IREDUC = 0
      IMODE  = 0
C ............................. Get the error tolerances from the work area
      AEPS = RREG(4)
      REPS = RREG(5)
C ---> 08.03.29/KMO ........... Handling of DELAY buffer reallocation
      iop = abs(jop)
      lev = min(jop,0)
      if (lev.lt.0 .and. iop.le.2) goto 130
      if (lev.lt.0 .and. iop.eq.3) goto 315
C ............................. Save the operation mode in the work area
      IREG(1) = IOP
C ---> 02.02.18/TI ............ Initialise discontinuity parameter and
C ............................. time splitting tolerance
      IREG(8) = 0
      RREG(9) = 1000.0*tiny(1.0d0)
C ............................. Initialize number of iterations per step
      IF (IOP .LT. 5) THEN
         IREG(9) = 0
      ELSE
         IREG(9) = IREG(9) + 1
      END IF
C ............................. Initialize max. end point of next step
C     RREG(2) = T + DT

      GOTO (100,100,300,400,500) IOP

  100 CONTINUE
      call initIO
C ............................. If the operation mode is set to initalization
C ............................. save the number of modules in the work area ..
      IREG(2) = NMOD
      IREG(6) = 0
C ............................. invoke a subroutine for initialization
      RREG(2) = T
      CALL RSINIT(T, RREG(2), NREG, NMOD, IREG(5), VREG,
     >            RREG(MPRREG(5)), IREG(MPIREG(2)), IREG(MPIREG(3)),
     >            IREG(MPIREG(6)), IREG(MPIREG(8)), IREG(MPIREG(7)),
     >            RREG(MPRREG(12)), DELAY, IERR, ISTAT, ZCTRL0)
C ---> 02.02.18/TI ............ Initialise last step's previous values
      DO 120 N = 1,NREG
         RREG(MPRREG(6)+2*NREG+N-1) = VREG(N)
         RREG(MPRREG(6)+3*NREG+N-1) = VREG(N)
  120 CONTINUE
C ............................. Initialise delay buffer
  130 CONTINUE
      IF (IREG(5) .GT. 0) THEN
         CALL DLAY(IREG(5), T, DT, RREG(2), IREG(2), VREG,
     >             IREG(MPIREG(7)), RREG(MPRREG(12)), IREG(MPIREG(6)),
     >             IREG(MPIREG(8)), 2, lev, DELAY, NDDIM, IERR)
         IF (IERR(1) .NE. 0) RETURN
      END IF
      GOTO 540

C ---> 02.09.12/TI............. TT #1151: SSTATE does not work in some cases
C  200 CONTINUE
C ............................. Determine local error tolerance
C      DO 210 N = 1,NREG
C         RREG(MPRREG(7)+N-1) = 0.5d0 * (AEPS + REPS * DABS(VREG(N)))
C  210 CONTINUE
C
C ............................. If the operation mode is set to steady state
C      IF (IOP .GT. 2) GOTO 300
C ............................. invoke a subroutine for steady state
C      CALL SSTATE(NREG,   VREG,  RREG(MPRREG(7)),
C     >            IREG(MPIREG(2)),MPIREG,IREG,
C     >            MPRREG, RREG,  DELAY,
C     >            IERR,  ISTAT, ZCTRL0, tolFactCtrl, iprint)
C
C      RETURN

  300 CONTINUE
C ............................. When the operation mode indicates a reduced
C ............................. time step ...
C ............................. Update the array of previous values

      DO 310 N = 1,NREG
         IF (IREG(MPIREG(2)+N-1) .EQ. 1) THEN
            RREG(MPRREG(2)+N-1) = RREG(MPRREG(6)+2*NREG+N-1)
         ELSE
            IF (IOP < 5 .AND. IREG(6) > 0) THEN
               RREG(MPRREG(6)+2*NREG+N-1) = 2.0*VREG(N)
     >                                    - RREG(MPRREG(2)+N-1)
            END IF
            RREG(MPRREG(2)+N-1) = VREG(N)
         END IF
         RREG(MPRREG(6)+3*NREG+N-1) = RREG(MPRREG(2)+N-1)
  310 CONTINUE
C ............................. Put initial values on top of delay buffer and
C ............................. eventually move the other variables downwards
 315  CONTINUE
      IF (IREG(5) .GT. 0) THEN
         CALL DLAY(IREG(5), T, DT, RREG(2), IREG(2), VREG,
     >             IREG(MPIREG(7)), RREG(MPRREG(12)), IREG(MPIREG(6)),
     >             IREG(MPIREG(8)), 2, lev, DELAY, NDDIM, IERR)
         IF (IERR(1) .NE. 0) RETURN
      END IF
C ............................. Determine local error tolerance
      DO 320 N = 1,NREG
         RREG(MPRREG(7)+N-1) = 0.5d0 * (AEPS + REPS * DABS(VREG(N)))
         RREG(MPRREG(7)+N-1+NREG) = RREG(MPRREG(7)+N-1)
  320 CONTINUE

  400 CONTINUE
C ............................. When the operation mode indicates a new or
C ............................. reduced time step
C ............................. Determine a starting vector for the iterations
C ............................. and indicate new Jacobian or iteration matrix

      DTP = RREG(3)
      IF (ABS(DTP-DT) .GT. 1.0D-12) IREG(6) = 1
      IF (DTP .GT. 0.0D0) THEN
         DO 410 N = 1,NREG
            IF (IREG(MPIREG(2)+N-1) .NE. 1) THEN
               VREGP = RREG(MPRREG(6)+2*NREG+N-1)
               IF (IOP < 4 .AND. IREG(6) > 0) THEN
                  VREG(N) = VREG(N) + (VREGP - VREG(N)) * DT/DTP
               ELSE
                  VREG(N) = VREGP
               END IF
            END IF
  410    CONTINUE
      END IF
      IF (IREG(6) > 0) THEN
         IF (IOP < 5 .AND. DTP .GT. 0.0) THEN
            IMODE = 2
         ELSE
            IMODE = 1
         END IF
      END IF
      IREG(6) = 0

C ............................. Save the current time and time step
      IF (IOP < 5) THEN
         RREG(1) = T
         RREG(3) = DT
      END IF

  500 CONTINUE
      IF (IREG(6) .NE. 0) GOTO 400
C ............................. Integrate the variables one time step forward
  510 CALL ONEDA(IMODE,           NREG,            RREG(MPRREG(2)),
     >           VREG,            RREG(MPRREG(3)), RREG(MPRREG(4)),
     >           RREG(MPRREG(6)), RREG(MPRREG(7)),
     >           IREG(MPIREG(2)), MPIREG,          IREG,
     >           MPRREG,          RREG,            DELAY,
     >           IERR,            ISTAT,           ZCTRL0,
     >           DCSPLIT,         tolFactCtrl,     iprint)
C ---> 10.08.10/KMO ........... Avoid infinite loop due to the GOTO 510 below
      if (ierr(1) .lt. 0) return

C ---> 02.02.10/TI ............ If a discontinuity is discovered, split
C ---> ........................ the time step. Find discontinuity by
C ---> ........................ interpolation and integrate from there.
      IF ( IREG(7) .EQ. 0 .AND. IREG(8) .NE. 0 ) THEN
         TAU = (RREG(2) - T + DT)/DT
         RREG(3) = (1.0 - TAU)*DT
         DO 520 N = 1,NREG
            VAR  = (1.0-TAU)*(1.0-TAU)*RREG(MPRREG(2)+ N-1)
     >           + TAU*(1.0-TAU)*RREG(MPRREG(4)+ N-1) + TAU*VREG(N)
            IF ( IREG(MPIREG(2)+ N-1) .NE. 1 ) THEN
               IF (IREDUC .EQ. 0) RREG(MPRREG(6)+2*NREG+N-1) = VREG(N)
               VREG(N) = VAR + (VAR-RREG(MPRREG(2)+N-1))*(1.0-TAU)/TAU
            END IF
            RREG(MPRREG(2)+ N-1) = VAR
  520    CONTINUE
         IREDUC = 1
         IMODE = 1
         GOTO 510
      END IF
      IF ( IREDUC .NE. 0 ) THEN
         IREDUC = 0
         IREG(6) = 1
         RREG(3) = DT
         DO 530 N = 1,NREG
            RREG(MPRREG(2)+ N-1) = RREG(MPRREG(6)+3*NREG+N-1)
  530    CONTINUE
      END IF
C <---

  540 CONTINUE

C ............................. Unless a discontinuity is discovered,
C ............................. save the inputs from this time step
      DO 550 N = 1,NREG
         IF ( IREG(MPIREG(2)+ N-1) .EQ. 1 ) THEN
            RREG(MPRREG(6)+2*NREG+N-1) = VREG(N)
         ELSE IF ( IREG(6) .EQ. 0 ) THEN
            RREG(MPRREG(6)+2*NREG+N-1) = 2.0*VREG(N) -
     >           RREG(MPRREG(6)+3*NREG+N-1)
         END IF
  550 CONTINUE

      RETURN
      END


C **********************************************************************

C               ****     ****   *   *   *   *   *****
C               *   *   *       *   **  *   *     *
C               ****     ***    *   * * *   *     *
C               * *         *   *   *  **   *     *
C               *   *   ****    *   *   *   *     *

C **********************************************************************

C  Name         : RSINIT

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Find values of constants which are configuration dependent
C                 and unchanged during a simulation. The legality of the con-
C                 figuration is checked.

C  Method       : The values are found by activating the modules through call
C                 to the central function subroutine

C  Call         : CALL RSINIT(T,TNM,NREG,NMOD,NTD,VREG,VDREG,MSTAT,MMOD,
C                >            MPMTOP,MMTOP,MPRPAR,RPAR,DELAY,IERR,ISTAT,
C                >            ZCTRL0)

C  Arguments    : T      - real, initial simulation time
C                 TNM    - real, maximum end point for next integration step
C                 NREG   - integer, no. of variables
C                 NMOD   - integer, no. of modules
C                 NTD    - integer, no. of time delays in configuration
C                 VREG   - real array(NREG), the vector of variables
C                 VDREG  - real array(NREG), function values
C                 MSTAT  - integer array(NREG), status flags for variables
C                 MMOD   - integer array(NREG), module number where each
C                          variable is calculated
C                 MPMTOP - integer array(NMOD+1), pointers to the content of
C                          MMTOP
C                 MMTOP  - integer array(MPMTOP(NMOD+1)), topology vector
C                 MPRPAR - integer array(NMOD+1), pointers to the content of
C                          RPAR
C                 RPAR   - real array(MPRPAR(NMOD+1)), module parameters
C                 DELAY  - real array(*), time history of delayed variables
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning
C                 ISTAT  - integer array(10), statistics vector
C                          ISTAT(1) = no. of function evaluation

C  Input        : T,TNM,NREG,NMOD,VREG,VDREG,MSTAT,MMOD,MPMTOP,MMTOP,
C                 > MPRPAR,RPAR,ISTAT

C  Output       : TNM,MSTAT,NTD,MMOD,DELAY,IERR,ISTAT

C  Work area    : None

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 FCN - for initalization in the modules

C  Fault react. : Return with non-zero error flag

C  Note         : Nothing

C **********************************************************************

      SUBROUTINE RSINIT(T, TNM, NREG, NMOD, NTD, VREG, VDREG, MSTAT,
     >     MMOD, MPMTOP, MMTOP, MPRPAR, RPAR, DELAY, IERR, ISTAT,
     >     ZCTRL0)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VREG(NREG), VDREG(NREG), RPAR(*), DELAY(*)
      INTEGER MSTAT(NREG), MMOD(NREG), MPMTOP(*), MMTOP(*), MPRPAR(*),
     >     IERR(2), ISTAT(10)
      LOGICAL ZCTRL0

C ............................. Reset the status vector in case it is set
      DO 10 N = 1,NREG
         IF (IABS(MSTAT(N)) .NE. 1) MSTAT(N) = 0
   10 CONTINUE
C ............................. Initialize the number of zeros in MSTAT
      NZB = NMOD
C ............................. Reset number of time delays in configuration
  100 NTD = 0
C ............................. Call function routine to determine the status
C ............................. vector
      CALL FCN(1, T, TNM, NREG, NMOD, NTD, VREG, VDREG, VREG, DELAY,
     >         VREG, MSTAT, MMOD, MPMTOP, MMTOP, MPRPAR, RPAR, IERR,
     >         IERR, ISTAT, ZCTRL0)

C ............................. Check the error flag
      IF (IERR(1) .LT. 0) RETURN
C ............................. Insert initial time in time delay buffer
      IF (NTD .GT. 0) DELAY(1) = T
C ............................. Count up the number of unlocated outputs
      NZA = 0
      DO 110 N = 1,NREG
         IF (MSTAT(N) .EQ. 0) NZA = NZA + 1
  110 CONTINUE
C ............................. If the number of unlocated outputs ...
C ............................. are zero, all status flags are set ...
      IF (NZA .EQ. 0) RETURN
C ............................. has decreased, repeat the function call ...
      IF (NZA .LT. NZB) THEN
         NZB = NZA
         GOTO 100
      END IF
C ............................. is unchanged, an algebraic loop has occurred
C ............................. and must be resolved
C ............................. Find the first unlocated output node ...
      DO 120 N = 1,NREG
         IF (MSTAT(N) .EQ. 0) GOTO 122
  120 CONTINUE
  122 NODE = N
C ............................. and set the corresponding status flag
      MSTAT(NODE) = 3
C ............................. Find a module where this variable can be an
C ............................. output
      MSTOP = MPMTOP(NMOD+1) - 1
      DO 130 M = 1,MSTOP
         IF (MMTOP(M) .EQ. NODE) GOTO 132
  130 CONTINUE
  132 MFIRST = M
      DO 134 M = 1,NMOD
         IF (MPMTOP(M+1) .GT. MFIRST) GOTO 136
  134 CONTINUE
  136 MODNO = MPMTOP(M)
C ............................. and allocate the variable to this module
      MMOD(NODE) = MODNO
      NZB = NZB - 1
      GOTO 100

      END


C **********************************************************************

C                ****    ****   *****    ***    *****   *****
C               *       *         *     *   *     *     *
C                ***     ***      *     *****     *     ****
C                   *       *     *     *   *     *     *
C               ****    ****      *     *   *     *     *****

C **********************************************************************

C  Name         : SSTATE

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Find the initial steady state of the configuration

C  Method       : The steady state equations are obtained by setting the time
C                 derivative to zero. The resulting equations are solved by
C                 Newton iteration, for which the LU-decomposed Jacobian is
C                 needed.

C  Call         : CALL SSTATE(NREG,VREG,TOL,MSTATE,MPIREG,IREG,MPRREG,RREG,
C                >            DELAY,IERR,ISTAT,ZCTRL0,tolFactCtrl,iprint)

C  Arguments    : NREG   - integer, no. of variables
C                 VREG   - real array(NREG), the vector of variables
C                 TOL    - real array(NREG), local error tolerance
C                 MSTAT  - integer array(NREG), status flags for variables
C                 MPIREG - integer array(*), pointers to the content of IREG,
C                          see calling routine
C                 IREG   - integer array(*), working area for integer variab-
C                          les, see calling routine
C                 MPRREG - integer array(*), pointers to the content of RREG,
C                          see calling routine
C                 RREG   - real array(*), working area for real variables,
C                          see calling routine
C                 DELAY  - real array(*), time history of delayed variables
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning
C                 ISTAT  - integer array(10), statistics vector
C                          ISTAT(2) = no. of functional evaluations for
C                                     steady state and integration

C  Input        : NREG,VREG,TOL,ISTAT

C  Output       : VREG,DELAY,IERR

C  Work area    : MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:

C  Fault react. : Return with non-zero error flag

C  Note         : Very preliminary documentation

C **********************************************************************

      SUBROUTINE SSTATE(NREG,   VREG,   TOL,   MSTAT,  MPIREG, IREG,
     >     MPRREG, RREG,   DELAY, IERR,   ISTAT,
     >     ZCTRL0, tolFactCtrl, iprint)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VREG(NREG), TOL(NREG), RREG(*), DELAY(*)
      INTEGER MSTAT(NREG), MPIREG(*), IREG(*), MPRREG(*), IERR(2),
     >     ISTAT(10)
      LOGICAL ZCTRL0
      EXTERNAL CALFCN

C ............................. Derive the LU-decomposed
C ............................. numerically computed Jacobian
      CALL NEWTON(0, NREG, RREG(3), RREG(7), VREG, VREG,
     >            RREG(MPRREG(5)), RREG(MPRREG(6)), IREG(MPIREG(2)),
     >            IREG(MPIREG(4)), IREG(MPIREG(5)), RREG(MPRREG(9)),
     >            RREG(MPRREG(10)), RREG(MPRREG(11)),
     >            MPIREG, IREG, MPRREG, RREG, DELAY, IERR, ISTAT,
     >            ZCTRL0, tolFactCtrl, iprint)

C ............................. Check the error flag
      IF (IERR(1) .LT. 0) RETURN

C ............................. Determine the steady state solution
      CALL EQUIL(NREG, RREG(6), VREG, RREG(MPRREG(10)), CALFCN,
     >           IREG(MPIREG(4)), RREG(MPRREG(6)), TOL, MSTAT ,
     >           MPIREG, IREG, MPRREG, RREG, DELAY, IERR, ISTAT,
     >           ZCTRL0)

C ............................. Check the error flag
      IERR(1) = IERR(1)*3
      IF (IERR(1) .LT. 0) RETURN

      RETURN
      END


C **********************************************************************

C                ***    *   *   *****   ***      ***
C               *   *   **  *   *       *   *   *   *
C               *   *   * * *   ****    *   *   *****
C               *   *   *  **   *       *   *   *   *
C                ***    *   *   *****   ***     *   *

C **********************************************************************

C  Name         : ONEDA

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Integrate the differential/algebraic system one step for-
C                 ward.

C  Method       : In order to estimate the truncation error two methods of
C                 different order are implemented. These are Backward Euler
C                 (1'st order) and Lobatto method IIIC (2'nd order), which
C                 are both implicit. The composite method proceeds with local
C                 extrapolation. The numerical equations are solved using
C                 Newton iteration, for which the Jacobian and the LU-decom-
C                 posed Newton matrices are needed. The maximum norm is used
C                 to express the truncation error.

C  Call         : CALL ONEDA(IMODE,NREG,VPREG,VREG,VBREG,VLREG,VHREG,TOL,
C                >      MSTAT,MPIREG,IREG,MPRREG,RREG,DELAY,IERR,ISTAT,
C                >      ZCTRL0,DCSPLIT,tolFactCtrl,iprint)

C  Arguments    : IMODE  - integer, operation mode
C                        = 0 : normal mode, using old iteration matrices
C                        = 1 : a new Jacobian is needed
C                        = 2 : a new step length requires new iteration
C                              matrices
C                 NREG   - integer, no. of variables
C                 VPREG  - real array(NREG), the previous values of the con-
C                          trol variables
C                 VREG   - real array(NREG), the most recent values of the
C                          control variables
C                 VBREG  - real array(NREG), the numerical solution using
C                          Backward Euler
C                 VLREG  - real array(2*NREG), the numerical solution using
C                          Lobatto IIIC
C                 VHREG  - real array(2*NREG), work area for subroutines
C                 TOL    - real array(2*NREG), local error tolerance
C                 MSTAT  - integer array(NREG), status flags for variables
C                 MPIREG - integer array(*), pointers to the content of IREG,
C                          see calling routine
C                 IREG   - integer array(*), working area for integer variab-
C                          les, see calling routine
C                 MPRREG - integer array(*), pointers to the content of RREG,
C                          see calling routine
C                 RREG   - real array(*), working area for real variables,
C                          see calling routine
C                 DELAY  - real array(*), time history of delayed variables
C                 IERR   - integer array(2), error flag
C                          IERR(1) = -6 : truncation error limit exceeded
C                          IERR(1) = -5 : divergence in Lobatto IIIC
C                          IERR(1) = -4 : divergence in Backward Euler
C                          IERR(1) =  0 : no error
C                          IERR(1) =  4 : slow convergence in Backward Euler
C                          IERR(1) =  5 : slow convergence in Lobatto IIIC
C                 ISTAT  - integer array(10), statistics vector
C                          ISTAT(2) = no. of function evaluations for steady
C                                     state and integration

C  Input        : IMODE,NREG,TOL,VPREG,VREG,MSTAT,ISTAT

C  Output       : VREG,IERR,ISTAT,DELAY

C  Work area    : VBREG,VLREG,VHREG,ERR,MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 NEWTON - for numerical calculation of the Jacobian and
C                          generation of iteration matrices for Backward
C                          Euler and Lobatto IIIC, including LU-factorization
C                 EQUIL  - for iteration to equilibrium
C                 DLAY   - for updating of the delay buffer

C  Fault react. : Return with non-zero error flag

C  Note         : Nothing

C **********************************************************************

      SUBROUTINE ONEDA(IMODE,  NREG,   VPREG,
     >                 VREG,   VBREG,  VLREG,  VHREG,  TOL,
     >                 MSTAT,  MPIREG, IREG,
     >                 MPRREG, RREG,   DELAY,
     >                 IERR,   ISTAT,  ZCTRL0, DCSPLIT,
     >                 tolFactCtrl, iprint)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VPREG(NREG), VREG(NREG), VBREG(NREG),
     >     VLREG(2*NREG), VHREG(2*NREG), TOL(2*NREG), RREG(*), DELAY(*)
      INTEGER MSTAT(NREG), MPIREG(*), IREG(*), MPRREG(*), IERR(2),
     >     ISTAT(10)
      LOGICAL ZCTRL0, DCSPLIT, NEWJAC, DISCON
      EXTERNAL BEULER, LOBATT

C ............................. Initialize new Jacobian and discontinuity flags
      NEWJAC = .FALSE.
      DISCON = .FALSE.
C ............................. Save maximum end point of next time step
      TNM = RREG(2)

C ............................. In normal mode, ommit generation of new
C ............................. iteration matrices
      IF (IMODE .EQ. 0) GOTO 110

C ............................. Form new iteration matrices
  100 CALL NEWTON(IMODE, NREG, RREG(3), RREG(7), VREG, VPREG,
     >            RREG(MPRREG(5)), VHREG, IREG(MPIREG(2)),
     >            IREG(MPIREG(4)), IREG(MPIREG(5)), RREG(MPRREG(9)),
     >            RREG(MPRREG(10)), RREG(MPRREG(11)),
     >            MPIREG, IREG, MPRREG, RREG, DELAY, IERR, ISTAT,
     >            ZCTRL0, tolFactCtrl, iprint)

C ............................. Restore maximum end point of next time step
      RREG(2) = TNM
C ............................. Check error flag
      IF (IERR(1) .LT. 0) RETURN
C ............................. In case of new Jacobian, set the flag
      IF (IMODE .EQ. 1) NEWJAC = .TRUE.

  110 CONTINUE
C ............................. Determine starting value for iteration in
C ............................. Backward Euler
      DO 120 N = 1,NREG
         VBREG(N) = VREG(N)
  120 CONTINUE

C ............................. Find the Backward Euler solution
      CALL EQUIL(NREG, RREG(6), VBREG, RREG(MPRREG(10)), BEULER,
     >           IREG(MPIREG(4)), VHREG, TOL, MSTAT,
     >           MPIREG, IREG, MPRREG, RREG, DELAY, IERR, ISTAT,
     >           ZCTRL0)
      IF (IERR(1) .LT. -1) RETURN
      IERR(1) = IERR(1)*4

C ............................. If a new Jacobian is requested and it is
C ............................. not fresh, renew the Jacobian
      IF (.NOT.NEWJAC .AND. IERR(1) .NE. 0) THEN
         IERR(1) = 0
         IMODE = 1
         GOTO 100
      END IF
C ............................. In case of new Jacobian and divergence, return
      IF (NEWJAC .AND. IERR(1) .LT. 0) RETURN

C ............................. Reset limit crossing indicator
      IREG(3) = 0

C ............................. Determine starting values for iteration in
C ............................. Lobatto IIIC
      DO 130 N = 1,NREG
         VLREG(N) = VPREG(N)
         VLREG(NREG+N) = VBREG(N)
  130 CONTINUE

C ............................. Find the Lobatto IIIC solution
      CALL EQUIL(2*NREG, RREG(6), VLREG, RREG(MPRREG(11)), LOBATT,
     >           IREG(MPIREG(5)), VHREG, TOL, MSTAT,
     >           MPIREG, IREG, MPRREG, RREG, DELAY, IERR, ISTAT,
     >           ZCTRL0)
      IF (IERR(1) .LT. -1) RETURN
      IERR(1) = IERR(1)*5

C ............................. If a new Jacobian is requested and it is
C ............................. not fresh, renew the Jacobian
      IF (.NOT.NEWJAC .AND. IERR(1) .NE. 0) THEN
         IERR(1) = 0
         IMODE = 1
         GOTO 100
      END IF
C ............................. In case of new Jacobian and divergence, return
      IF (NEWJAC .AND. IERR(1) .LT. 0) RETURN

C ............................. Transfer the result to the solution vector
      DO 140 N = 1,NREG
         VREG(N) = VLREG(NREG+N)
  140 CONTINUE

C ............................. Estimate the local truncation error
      TRUNC = 0.0d0
      DO 150 N = 1,NREG
C         TRUNC = MAX(TRUNC, DABS(VREG(N)-VBREG(N))/TOL(N))
         ABSDIF = DABS(VREG(N)-VBREG(N))/TOL(N)
         IF (ABSDIF .GT. TRUNC) THEN
            TRUNC = ABSDIF
            IERR(2) = N
         END IF
 150  CONTINUE

C ............................. If the estimated truncation error is above
C ............................. the tolerance, set the error flag before
C ............................. return
C ---> 92.04.23/TI ............ Changed sign on IERR
      IF (TRUNC .GT. 1) IERR(1) = 6
      RREG(8) = TRUNC

      IF (IREG(5) .GT. 0) THEN
C ............................. Put delayed variables on top of delay buffer
         CALL DLAY(IREG(5), RREG(1), RREG(3), RREG(2), IREG(2), VLREG,
     >             IREG(MPIREG(7)), RREG(MPRREG(12)), IREG(MPIREG(6)),
     >             IREG(MPIREG(8)), 2, 1, DELAY, 1, IERR)
         CALL DLAY(IREG(5), RREG(1), RREG(3), RREG(2), IREG(2), VREG,
     >             IREG(MPIREG(7)), RREG(MPRREG(12)), IREG(MPIREG(6)),
     >             IREG(MPIREG(8)), 2, 2, DELAY, 1, IERR)
      END IF

C ............................. Set end of next time step according to
C ............................. discontinuities or truncation error
      IF (DCSPLIT .AND. IREG(3) .NE. 0 .AND.
     >    IABS(IREG(3)) .NE. IREG(8)) THEN
         Y2P = RREG(MPRREG(6) + 2*NREG + IREG(4) - 1)
         Y2 = VREG(IREG(4))
         YLIM = RREG(MPRREG(12) + IABS(IREG(3)) - 1)
         YTOL = RREG(MPRREG(7) + IREG(4) - 1)
         IF ( IREG(9) > 10 .AND. DABS(Y2P - YLIM) < YTOL .OR.
     >        IREG(3) < 0 .AND. Y2 > YLIM .OR.
     >        IREG(3) > 0 .AND. Y2 < YLIM ) THEN
            IREG(3) = 0
         ELSE
            DISCON = .TRUE.
            YP = VPREG(IREG(4))
            Y1 = VLREG(IREG(4))
            A1 = YP - Y1
            A2 = YP - Y2
            A3 = YP - YLIM

            IF (IABS(MSTAT(IREG(4))) .EQ. 1) THEN
               RREG(2) = RREG(1) + (A3 / A2 - 1.0d0) * RREG(3)
            ELSE
               TAU = 1.0d0
               if ( abs(A1) < RREG(9)) then
                  if ( abs(A2) > RREG(9)) TAU = A3/A2
               else
                  B = 0.5d0 * (A1 + A2)
                  ARG = 1.0d0 - A1*A3 /(B*B)
                  IF (ARG .LT. 0.0d0) THEN
                     IERR(1) = -7
                     IERR(2) = IREG(4)
                     IREG(8) = 0
                     return
                  ELSE
                     TAU = B*(1.0d0 - SQRT(ARG))/A1
                  END IF
               end if
               RREG(2) = RREG(1) - (1.0d0 - TAU)*RREG(3)
               IF ( TAU > 1.0d0 - RREG(9)) THEN
                  IREG(8) = 0
               ELSE
                  IREG(8) = IABS(IREG(3))
                  IF ( TAU < RREG(9) ) RREG(2) = RREG(1)
               END IF
            END IF
         END IF
      END IF
      IF (.NOT.DISCON) THEN
         IREG(8) = 0
         IF (TRUNC .GT. 1.0) THEN
            RREG(2) = RREG(1) - RREG(3) + 0.8d0*RREG(3)/SQRT(TRUNC)
         ELSE IF (TRUNC .GT. 0.0) THEN
            RREG(2) = RREG(1) + 0.8d0*RREG(3)/SQRT(TRUNC)
         END IF
      END IF
      RETURN
      END


C **********************************************************************

C               *   *   *****   *   *   *****    ***    *   *
C               **  *   *       *   *     *     *   *   **  *
C               * * *   ****    * * *     *     *   *   * * *
C               *  **   *       ** **     *     *   *   *  **
C               *   *   *****   *   *     *      ***    *   *

C  *********************************************************************

C  Name         : NEWTON

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Find the iteration matrices for steady state or numerical
C                 integration using Backward Euler and Lobatto IIIC.

C  Method       : The status flags are used to separate between inputs, state
C                 variables and algebraic variables. The steady state itera-
C                 tion matrix is the Jacobian, except that the column for
C                 inputs are cleared. The iteration matrices for Backward
C                 Euler and Lobatto IIIC are found from the Jacobian. The li-
C                 near dimension of the matrix for Lobatto IIIC is twice that
C                 for Backward Euler. In terms of block matrices the iteration
C                 matrix for Lobatto IIIC is antisymmetric and the two blocks
C                 on the diagonal are identical and, except for the time step,
C                 equal to the iteration matrix for Backward Euler.

C  Call         : CALL NEWTON(IMODE,NREG,H,DELTA,VREG,VDREG,VHREG,MSTAT,IPVBE,
C                >     IPVLO,AJAC,ABE,ALO,MPIREG,IREG,MPRREG,RREG,DELAY,IERR,
C                >     ISTAT,ZCTRL0,tolFactCtrl,iprint)

C  Arguments    : IMODE  - integer, operation mode
C                        = 0 : Steady state
C                        = 1 : A new Jacobian is needed
C                        > 1 : Iteration matrices based on the last Jacobian
C                 NREG   - integer, dimension of variable vector
C                 H      - real, time step
C                 DELTA  - real, relative perturbation
C                 VREG   - real array(NREG), the vector of variables
C                 VPREG  - real array(NREG), the vector of previous variables
C                 VDREG  - real array(NREG), function values
C                 VHREG  - real array(NREG), pertubed function values
C                 MSTAT  - integer array(NREG), status flags for variables
C                 IPVBE  - integer array(NREG), pivot indices for ABE
C                 IPVLO  - integer array(2*NREG), pivot indices for ALO
C                 AJAC   - real array(NREGxNREG), negative numerical Jacobian
C                          or iteration matrix for steady state
C                 ABE    - real array(NREGxNREG), iteration matrix for Back-
C                          ward Euler
C                 ALO    - real array(2*NREGx2*NREG), iteration matrix for
C                          Lobatto IIIC
C                 MPIREG - integer array(*), pointers to the content of IREG,
C                          see calling routine
C                 IREG   - integer array(*), working area for integer variab-
C                          les, see calling routine
C                 MPRREG - integer array(*), pointers to the content of RREG,
C                          see calling routine
C                 RREG   - real array(*), working area for real variables,
C                          see calling routine
C                 DELAY  - real array(*), time history of delayed variables
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning
C                 ISTAT  - integer array(10), statistcs vector
C                          ISTAT(3) = no. of Jacobi evaluations
C                          ISTAT(4) = no. of LU-factorizat. for steady state
C                          ISTAT(5) = no. of LU-factorizat. for Backward Euler
C                          ISTAT(6) = no. of LU-factorizat. for Lobatto IIIC

C  Input        : IMODE,NREG,H,DELTA,VREG,VDREG,VHREG,MSTAT,ISTAT,AJAC

C  Output       : IPVBE,IPVLO,AJAC,ABE,ALO,DELAY,IERR,ISTAT

C  Work area    : MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 FCN    - for functional evaluation
C                 DGETOL - for LU-factorization (LAPACK wrapper)

C  Fault react. : Return with a non-zero error flag

C  Note         : The routine is dependent on the actual numerical methods

C **********************************************************************

      SUBROUTINE NEWTON(IMODE, NREG, H, DELTA, VREG, VPREG, VDREG,
     >                  VHREG, MSTAT, IPVBE, IPVLO, AJAC, ABE, ALO,
     >                  MPIREG, IREG, MPRREG, RREG, DELAY,
     >                  IERR, ISTAT, ZCTRL0, tolFactCtrl, iprint)

      use ioModule, only : IO

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VREG(NREG), VPREG(NREG), VDREG(NREG),
     >     VHREG(NREG), AJAC(NREG,NREG), ABE(NREG,NREG),
     >     ALO(2*NREG,2*NREG), RREG(*), DELAY(*)
      INTEGER MSTAT(NREG), IPVBE(NREG), IPVLO(2*NREG), MPIREG(*),
     >     IREG(*), MPRREG(*), IERR(2), ISTAT(10)
      LOGICAL ZCTRL0, LPRINT

      LPRINT = iprint .eq. -98
C ............................. If a new Jacobian is needed, ...
      IF (IMODE .GT. 1) GOTO 200
C ............................. Find the unpertubed function values
      CALL FCN(IREG(1), RREG(1), RREG(2), NREG, IREG(2), IREG(5),
     >         VREG, VDREG, RREG(MPRREG(2)), DELAY, RREG(MPRREG(7)),
     >         MSTAT, IREG(MPIREG(3)), IREG(MPIREG(6)), IREG(MPIREG(8)),
     >         IREG(MPIREG(7)), RREG(MPRREG(12)),
     >         IERR, IREG(3), ISTAT(2), ZCTRL0)

C ............................. Check error flag
      IF (IERR(1) .LT. 0) RETURN
C ............................. For each variable, ...
      DO 110 N = 1,NREG
C ............................. perturbate the value, ...
         IF (VREG(N) .GT. 1.0E-4) THEN
C            DY = DELTA * VREG(N)
            DY = DELTA * DABS(VREG(N))
         ELSE
            DY = DELTA
         END IF
         IF (VREG(N) .LT. VPREG(N)) DY = -DY
         VREGN = VREG(N)
         VREG(N) = VREG(N) + DY
C ............................. invoke the function subroutine ...
         CALL FCN(IREG(1), RREG(1), RREG(2), NREG, IREG(2), IREG(5),
     >        VREG, VHREG, RREG(MPRREG(2)), DELAY, RREG(MPRREG(7)),
     >        MSTAT, IREG(MPIREG(3)), IREG(MPIREG(6)), IREG(MPIREG(8)),
     >        IREG(MPIREG(7)), RREG(MPRREG(12)),
     >        IERR, IREG(3), ISTAT(2), ZCTRL0)

C ............................. and compute the corresponding column in the
C ............................. negative Jacobian which is also the iteration
C ............................. matrix for steady state:

C ............................. l   I       0     0  l
C ............................. l                    l
C ............................. l   df     df     df l
C ............................. l - --   - --   - -- l
C ............................. l   du     dx     dz l
C ............................. l                    l
C ............................. l   dg     dg     dg l
C ............................. l - --   - --   - -- l
C ............................. l   du     dx     dz l

         DO 100 NN = 1,NREG
            AJAC(NN,N) = (VDREG(NN) - VHREG(NN))/DY
  100    CONTINUE
C ............................. Reset the value for this variable
C         VREG(N) = VREG(N) - DY
         VREG(N) = VREGN

  110 CONTINUE

C ............................. Insert identity row for inputs
      DO 114 N = 1,NREG
         IF (IABS(MSTAT(N)) .EQ. 1) THEN
            DO 112 NN = 1,NREG
               AJAC(N,NN) = 0.0d0
  112       CONTINUE
            AJAC(N,N) = 1.0d0
         END IF
  114 CONTINUE

C ............................. Debug print of Jacobian
      if (LPRINT) call WriMat(ajac,nreg,nreg,IO,'Jacobian')

C ............................. Update the number of Jacobi evaluations
      ISTAT(3) = ISTAT(3) + 1
C ---> 02.02.08/TI ............ Zero eventual discontinuity discovery
      IREG(3) = 0
C ............................. When the operation mode is steady state ..
      IF (IMODE .EQ. 0) THEN
C ............................. Transfer the Jacobian to ABE
         DO 121 N = 1,NREG
            DO 120 NN = 1,NREG
               ABE(N,NN) = AJAC(N,NN)
  120       CONTINUE
  121    CONTINUE

C ............................. LU factorize the Jacobian for steady state
         call dgetol(abe,nreg,nreg,ipvbe,tolFactCtrl,
     >               'Jacobian for steady state',LPRINT,IO,IERR)

C ............................. Update the no. of LU-fact. for steady state
         ISTAT(4) = ISTAT(4) + 1
         RETURN
      END IF
C ............................. Find the iteration matrix for Backward Euler:

C ............................. l    I           0         0 l
C ............................. l                            l
C ............................. l     df        df        df l
C ............................. l - h*--  I - h*--    - h*-- l
C ............................. l     du        dx        dz l
C ............................. l                            l
C ............................. l    dg         dg        dg l
C ............................. l  - --       - --      - -- l
C ............................. l    du         dx        dz l

  200 CONTINUE
C ............................. Clear iteration matrix for Backward Euler
      DO 202 NN = 1,NREG
         DO 201 N = 1,NREG
            ABE(NN,N) = 0.0d0
  201    CONTINUE
  202 CONTINUE

C ............................. For each variable ...
      DO 220 NN = 1,NREG
C ............................. find the elements in the corresponding
C ............................. row depending on whether the variable is
C ............................. an input
         IF (IABS(MSTAT(NN)) .EQ. 1) THEN
            ABE(NN,NN) = 1.0d0
C ............................. a state variable
         ELSE IF (MSTAT(NN) .EQ. 2) THEN
            DO 210 N = 1,NREG
               ABE(NN,N) = H * AJAC(NN,N)
  210       CONTINUE
            ABE(NN,NN) = ABE(NN,NN) + 1.0d0
C ............................. or an algebraic variable
         ELSE IF (MSTAT(NN) .EQ. 3) THEN
            DO 212 N = 1,NREG
               ABE(NN,N) = AJAC(NN,N)
  212       CONTINUE
         END IF
  220 CONTINUE

C ............................. Find the iteration matrix for Lobatto IIIC:

C ..... l     I          0         0  .     0          0         0  l
C ..... l                             .                             l
C ..... l   h df       h df      h df .   h df       h df      h df l
C ..... l - -*--   I - -*--    - -*-- .   -*--      -*--       -*-- l
C ..... l   2 du       2 dx      2 dz .   2 du       2 dx      2 dz l
C ..... l                             .                             l
C ..... l                dg        dg .                             l
C ..... l     0        - --      - -- .     0          0         0  l
C ..... l                dx        dz .                             l
C ..... l . . . . . . . . . . . . . . . . . . . . . . . . . . . . . l
C ..... l     0          0         0  .     I          0         0  l
C ..... l                             .                             l
C ..... l   h df       h df      h df .   h df       h df      h df l
C ..... l - -*--     - -*--    - -*-- . - -*--   I - -*--    - -*-- l
C ..... l   2 du       2 dx      2 dz .   2 du       2 dx      2 dz l
C ..... l                             .                             l
C ..... l                             .     dg         dg        dg l
C ..... l     0          0         0  .   - --       - --      - -- l
C ..... l                             .     du         dx        dz l


C ............................. Clear iteration matrix for Lobatto IIIC
      DO 302 NN = 1,NREG
         DO 301 N = 1,NREG
            ALO(NN,N) = 0.0d0
            ALO(NN,N+NREG) = 0.0d0
            ALO(NN+NREG,N) = 0.0d0
            ALO(NN+NREG,N+NREG) = 0.0d0
  301    CONTINUE
  302 CONTINUE

C ............................. For each variable ...
      DO 320 NN = 1,NREG
C ............................. find the elements in the corresponding
C ............................. row depending on whether the variable is
C ............................. an input,
         IF (IABS(MSTAT(NN)) .EQ. 1) THEN
            ALO(NN,NN) = 1.0d0
            ALO(NN+NREG,NN+NREG) = 1.0d0
C ............................. a state variable
         ELSE IF (MSTAT(NN) .EQ. 2) THEN
            DO 310 N = 1,NREG
               ALO(NN,N) = 0.5d0 * H * AJAC(NN,N)
               ALO(NN,N+NREG) = - 0.5d0 * H * AJAC(NN,N)
               ALO(NN+NREG,N) = 0.5d0 * H * AJAC(NN,N)
               ALO(NN+NREG,N+NREG) = 0.5d0 * H * AJAC(NN,N)
  310      CONTINUE
           ALO(NN,NN) = ALO(NN,NN) + 1.0d0
           ALO(NN+NREG,NN+NREG) = ALO(NN+NREG,NN+NREG) + 1.0d0
C ............................. or an algebraic variable
        ELSE IF (MSTAT(NN) .EQ. 3) THEN
           DO 312 N = 1,NREG
              ALO(NN,N) = AJAC(NN,N)
              ALO(NN+NREG,N+NREG) = AJAC(NN,N)
  312      CONTINUE
        END IF
  320 CONTINUE

C ..................................... LU-factorize the Newton matrix for
C ..................................... Backward Euler ...
      call dgetol(abe,nreg,nreg,ipvbe,tolFactCtrl,
     >            'Newton matrix for Backward Euler',LPRINT,IO,IERR)
      if (ierr(1) .lt. 0) return

C ..................................... Update the no. of LU-factorizations
C ..................................... for Backward Euler
      ISTAT(5) = ISTAT(5) + 1

C ..................................... LU-factorize the Newton matrix for
C ..................................... Lobatto IIIC
      call dgetol(alo,2*nreg,2*nreg,ipvlo,tolFactCtrl,
     >            'Newton matrix for Lobatto IIIC',LPRINT,IO,IERR)

C ..................................... Update the no. of LU-factorizations
C ..................................... for Lobatto IIIC
      ISTAT(6) = ISTAT(6) + 1

      RETURN
      END


C **********************************************************************

C               *****    ***    *   *   *   *
C               *       *   *   *   *   *   *
C               ****    *   *   *   *   *   *
C               *       * * *   *   *   *   *
C               *****    ** *    ***    *   *****

C **********************************************************************

C  Name         : EQUIL

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Find the solution of the general equation F(y) = 0

C  Method       : The solution is found through Newton iteration using a LU-
C                 decomposed iteration matrix. The accuracy and rate of con-
C                 vergence are controlled.

C  Call         : CALL EQUIL(NDIM,ACCU,Y,A,FRHS,IPVT,DY,TOL,MSTAT,MPIREG,
C                >           IREG,MPRREG,RREG,DELAY,IERR,ISTAT,ZCTRL0)

C  Arguments    : NDIM   - integer, dimension of the vector y
C                 ACCU   - real, accuracy parameter
C                 Y      - real array(NDIM), solution vector
C                 A      - real array(NDIMxNDIM), iteration matrix
C                 FRHS   - external subroutine, right hand side of numerical
C                          equation
C                 IPVT   - integer array(NDIM), pivot indices for LU-factor-
C                          ized iteration matrix
C                 DY     - real array(NDIM), solution increment per iteration
C                 TOL    - real array(NDIM), local error tolerance
C                 MSTAT  - integer array(NDIM), variable type indicators
C                 MPIREG - integer array(*), pointers to the content of IREG,
C                          see calling routine
C                 IREG   - integer array(*), working area for integer variab-
C                          les, see calling routine
C                 MPRREG - integer array(*), pointers to the content of RREG,
C                          see calling routine
C                 RREG   - real array(*), working area for real variables,
C                          see calling routine
C                 DELAY  - real array(*), time history of delayed variables
C                 IERR   - integer array(2), error flag
C                          IERR(1) = -1 : divergence
C                          IERR(1) =  0 : no error
C                          IERR(1) = +1 : slow convergence
C                          IERR(2) - variable no. in maximum norm
C                 ISTAT  - integer array(10), statistics vector
C                          ISTAT(10) - no. of times solving a linear system

C  Input        : NDIM,ACCU,Y,A,IPVT,TOL,FRHS,ISTAT

C  Output       : Y,DELAY,IERR,ISTAT

C  Work area    : DY,MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 FCN         - for function evaluation
C                 BEULER - for right hand side using Backward Euler
C                 LOBATT - for right hand side using Lobatto IIIC
C                 DGETRS - for solution of the Newton iteration equation

C  Fault react. : Return with non-zero error flag

C  Note         : Nothing

C **********************************************************************

      SUBROUTINE EQUIL(NDIM, ACCU, Y, A, FRHS, IPVT, DY, TOL, MSTAT,
     >     MPIREG, IREG, MPRREG, RREG, DELAY, IERR, ISTAT, ZCTRL0)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION Y(NDIM), DY(NDIM), TOL(NDIM), A(NDIM,NDIM),
     >     RREG(*), DELAY(*)

      INTEGER IPVT(NDIM), MSTAT(NDIM), MPIREG(*), IREG(*), MPRREG(*),
     >     IERR(2),ISTAT(10)

      LOGICAL ZCTRL0

C ............................. Initialize rate of convergence and iteration
C ............................. number
      DYMAXP = 0.0d0
      ROC = 0.5d0
      NIT = 0
C     nWarn = 0
C ............................. Reset the number of slow convergences
  100 CONTINUE
      NSLOW = 0

C ............................. Start iteration loop
  110 CONTINUE
      NIT = NIT + 1
C ............................. Transfer variables to increment vector
      DO 120 N = 1,NDIM
         DY(N) = Y(N)
  120 CONTINUE

C ............................. Invoke subroutine to get rhs = F(y(k))
      CALL FRHS(NDIM, RREG(3), RREG(MPRREG(2)), DY, RREG(MPRREG(5)),
     >          IREG(MPIREG(2)), MPIREG, IREG, MPRREG, RREG, DELAY,
     >          IERR, ISTAT, ZCTRL0)
C ............................. Check the error flag
      IF (IERR(1) .LT. 0) RETURN

C ............................. Solve the system
C .............................   A * (y(k+1) - y(k)) = F(y(k))
      CALL DGETRS('N', NDIM, 1, A, NDIM, IPVT, DY, NDIM, IERR)
C ............................. Check the error flag
      IF (IERR(1) .LT. 0) RETURN
C ............................. Update the no. of times solving a system
      ISTAT(10) = ISTAT(10) + 1
C ............................. Remove increments for input variables
      DO 124 N = 1,NDIM
         IF (IABS(MSTAT(N)) .EQ. 1) DY(N) = 0.0d0
  124 CONTINUE

C ............................. Update the variables
      DO 130 N = 1,NDIM
         Y(N) = Y(N) + DY(N)
  130 CONTINUE

C ............................. Find the maximum norm of the increments
      DYMAX = 0.0d0
      DO 140 N = 1,NDIM
         ABSDY = DABS(DY(N))/TOL(N)
         IF (ABSDY .GT. DYMAX) THEN
            DYMAX = ABSDY
            IERR(2) = N
         END IF
  140 CONTINUE
      IF (DYMAX .LE. 0.0d0) RETURN

C ............................. From second iteration compute rate of conver-
C ............................. gence
      IF (NIT .EQ. 1) GOTO 110
      ROC = DYMAX/DYMAXP
C ............................. If the solution is accurate enough, return
      IF (DYMAX .LT. 0.5d0*ACCU/ROC) GOTO 200

C ............................. In case of convergence, continue iteration
      DYMAXP = DYMAX
      IF (ROC .LE. 0.5) GOTO 100
C ............................. In case of divergence, set error flag to
C ............................. "divergence" and return
C ---> 02.02.10/TI ............ This is not abnormal, warnings are skipped
C      if (ROC.gt.10.0D0) then
C         nWarn = nWarn+1
C         if ( nWarn < 10 ) then
C            write(IO,*) '  ** Non-convergence in control module,',
C     >             ' Time =',RREG(1),' Iteration no.',NIT,' ROC =',ROC
C         else if (nWarn == 10) then
C            write(IO,*) '  ** Non-convergence in control module,',
C     >             ' Time =',RREG(1),' Iteration no.',NIT,' ROC =',ROC,
C     >             ' (ignoring further warnings)'
C         end if
C      end if
C <---

      IF (NIT.GT.200) THEN
         IERR(1) = -1
         RETURN
      END IF
C ............................. If the iteration shows slow convergence, but
C ............................. no more than three times, try another
C ............................. iteration
      IF (NSLOW .LT. 3) THEN
         NSLOW = NSLOW + 1
         GOTO 110
      END IF
C ............................. Otherwise, set error flag to "slow conver-
C ............................. gence" and return
      IERR(1) = 1
  200 CONTINUE
C ---> 02.02.10/TI ............ This warning is skipped along with those above
C      if ( nWarn > 0 ) then
C         write(IO,*) '  ** Convergence in control module after ', NIT,
C     >               ' iterations'
C      end if
C <---

      RETURN
      END


C **********************************************************************

C                ****    ***    *       *****    ****   *   *
C               *       *   *   *       *       *       **  *
C               *       *****   *       ****    *       * * *
C               *       *   *   *       *       *       *  *
C                ****   *   *   *****   *        ****   *   *

C **********************************************************************

C  Name         : CALFCN

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Provide the right hand side of the numerical equation for
C                 steady state

C  Method       : Inputs, state variables and algebraic variables are included
C                 in the vector of variables. This is essentially a calling
C                 routine for the function evaluation routine FCN. Input
C                 components are not affected.

C  Call         : CALL CALFCN(NREG,DT,VPREG,VREG,VDREG,MSTAT,MPIREG,IREG,
C                >            MPRREG,RREG,DELAY,IERR,ISTAT,ZCTRL0)

C  Arguments    : NREG   - integer, dimension of the variable vector
C                 DT     - real, time step
C                 VPREG  - real array(NREG), the previous values of the vari-
C                          ables
C                 VREG   - real array(NREG), on input: the most recent itera-
C                          tion value using this integration method, on out-
C                          put: the right hand side of the numerical equation
C                 VDREG  - real array(NREG), function values
C                 MSTAT  - integer array(NREG), status flags for variables
C                 MPIREG - integer array(*), pointers to the content of IREG,
C                          see calling routine
C                 IREG   - integer array(*), working area for integer variab-
C                          les, see calling routine
C                 MPRREG - integer array(*), pointers to the content of RREG,
C                          see calling routine
C                 RREG   - real array(*), working area for real variables,
C                          see calling routine
C                 DELAY  - real array(*), time history of delayed varables
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning
C                 ISTAT  - integer array(10), statistics vector
C                          ISTAT(7) = no. of iterations for steady state

C  Input        : NREG,DT,VPREG,VREG,VDREG,MSTAT,ISTAT

C  Output       : VREG,IERR,ISTAT

C  Work area    : MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 FCN - for function evaluation

C  Fault react. : None

C  Note         : Nothing

C **********************************************************************

      SUBROUTINE CALFCN(NREG, DT, VPREG, VREG, VDREG, MSTAT, MPIREG,
     >     IREG, MPRREG, RREG, DELAY, IERR, ISTAT, ZCTRL0)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VPREG(NREG), VREG(NREG), VDREG(NREG), RREG(*),
     >     DELAY(*)
      INTEGER MSTAT(NREG), MPIREG(*), IREG(*), MPRREG(*),
     >     IERR(2), ISTAT(10)
      LOGICAL ZCTRL0

C     The right hand side of the numerical equation for steady state is:

C     0 = -u(n+1) + s(t(n+1))
C     0 =  f(u(n+1),x(n+1),z(n+1))
C     0 =  g(u(n+1),x(n+1),z(n+1))

C ............................. Invoke the function routine to obtain the
C ............................. right hand side
      CALL FCN(IREG(1), RREG(1), RREG(2), NREG, IREG(2), IREG(5),
     >     VREG, VDREG, RREG(MPRREG(2)), DELAY, RREG(MPRREG(7)),
     >     MSTAT, IREG(MPIREG(3)), IREG(MPIREG(6)), IREG(MPIREG(8)),
     >     IREG(MPIREG(7)), RREG(MPRREG(12)), IERR, IREG(3), ISTAT(2),
     >     ZCTRL0)

C ............................. Check the error vector
      IF (IERR(1) .LT. 0) RETURN

C ............................. For each variable ...
      DO 100 N = 1,NREG
C ............................. identify variable status
         ITYP = IABS(MSTAT(N))
C ............................. When the variable is an input,
         IF (ITYP .EQ. 1) THEN
C ............................. set rhs to zero ...
            VREG(N) = 0.0d0
C ............................. and when the variable is a state variable
C ............................. or an algebraic variable,
         ELSE IF (ITYP .GT. 1) THEN
C ............................. assign the functional value to the rhs
            VREG(N) = VDREG(N)
         ELSE
C ............................. will never get here, to suppress gfortran warning
            VREG(N) = DT*VPREG(N)
         END IF
  100 CONTINUE

C ............................. Update the no. of ierations for steady state
      ISTAT(7) = ISTAT(7) + 1

      RETURN
      END


C **********************************************************************

C               ****    *****   *   *   *       *****   ****
C               *   *   *       *   *   *       *       *   *
C               ****    ****    *   *   *       ****    ****
C               *   *   *       *   *   *       *       * *
C               ****    *****    ***    *****   *****   *   *

C **********************************************************************

C  Name         : BEULER

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Provide the right hand side of the numerical equation using
C                 Backward Euler

C  Method       : Inputs, state variables and algebraic variables are included
C                 in the vector of variables. The right hand side components
C                 are a combination of function values, time step and previous
C                 values. Input components are not affected.

C  Call         : CALL BEULER(NREG,DT,VPREG,VBREG,VDREG,MSTAT,MPIREG,IREG,
C                >            MPRREG,RREG,DELAY,IERR,ISTAT,ZCTRL0)

C  Arguments    : NREG   - integer, dimension of the variable vector
C                 DT     - real, time step
C                 VPREG  - real array(NREG), the previous values of the vari-
C                          ables
C                 VBREG  - real array(NREG), on input: the most recent itera-
C                          tion value using this integration method, on out-
C                          put: the right hand side of the numerical equation
C                 VDREG  - real array(NREG), function values
C                 MSTAT  - integer array(NREG), status flags for variables
C                 MPIREG - integer array(*), pointers to the content of IREG,
C                          see calling routine
C                 IREG   - integer array(*), working area for integer variab-
C                          les, see calling routine
C                 MPRREG - integer array(*), pointers to the content of RREG,
C                          see calling routine
C                 RREG   - real array(*), working area for real variables,
C                          see calling routine
C                 DELAY  - real array(*), time history of delayed variables
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning
C                 ISTAT  - integer array(10), statistics vector
C                          ISTAT(8) = no. of iterations with Backward Euler

C  Input        : NREG,DT,VPREG,VBREG,VDREG,MSTAT,ISTAT

C  Output       : VBREG,IERR,ISTAT

C  Work area    : MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                 FCN - for function evaluation

C  Note         : Nothing

C  Fault react. : None

C **********************************************************************

      SUBROUTINE BEULER(NREG, DT, VPREG, VBREG, VDREG, MSTAT, MPIREG,
     >                  IREG, MPRREG, RREG, DELAY, IERR, ISTAT,
     >                  ZCTRL0)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VPREG(NREG), VBREG(NREG), VDREG(NREG), RREG(*),
     >                 DELAY(*)
      INTEGER MSTAT(NREG), MPIREG(*), IREG(*), MPRREG(*),
     >        IERR(2), ISTAT(10)
      LOGICAL ZCTRL0

C     The right hand side of the numerical equation for Backward Euler is:

C     0 = -U(n+1) + s(t(n+1))
C     0 = -X(n+1) + x(n) + h*f(U(n+1),X(n+1),Z(n+1))
C     0 =  g(U(n+1),X(n+1),Z(n+1))

C ............................. Save maximum end point of next time step
      TNM = RREG(2)
C ............................. Invoke the function routine to obtain the
C ............................. function values
      CALL FCN(IREG(1), RREG(1), RREG(2), NREG, IREG(2), IREG(5),
     >     VBREG, VDREG, RREG(MPRREG(2)), DELAY, RREG(MPRREG(7)),
     >     MSTAT, IREG(MPIREG(3)), IREG(MPIREG(6)), IREG(MPIREG(8)),
     >     IREG(MPIREG(7)), RREG(MPRREG(12)), IERR, IREG(3), ISTAT(2),
     >     ZCTRL0)

C ............................. Restore maximum end point of next time step
      RREG(2) = TNM
C ............................. Check the error flag
      IF (IERR(1) .LT. 0) RETURN

C ............................. For each variable ...
      DO 100 N = 1,NREG
C ............................. identify variable status
         ITYP = IABS(MSTAT(N))
C ............................. When the variable is an input,
         IF (ITYP .EQ. 1) THEN
C ............................. set rhs to zero
            VBREG(N) = 0.0d0
C ............................. When the variable is a state variable,
         ELSE IF (ITYP .EQ. 2) THEN
C ............................. compute the the rhs expression
            VBREG(N) = VPREG(N) - VBREG(N) + DT*VDREG(N)
C ............................. When the variable is an algebraic variable,
         ELSE IF (ITYP .EQ. 3) THEN
C ............................. assign the functional value to the rhs
            VBREG(N) = VDREG(N)
         END IF
  100 CONTINUE

C ............................. Update the no. of iterations with Backward
C ............................. Euler
      ISTAT(8) = ISTAT(8) + 1

      RETURN
      END


C **********************************************************************

C               *        ***    ****     ***    *****   *****
C               *       *   *   *   *   *   *     *       *
C               *       *   *   ****    *****     *       *
C               *       *   *   *   *   *   *     *       *
C               *****    ***    ****    *   *     *       *

C **********************************************************************

C  Name         : LOBATT

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Provide the right hand side of the numerical equation using
C                 Lobatto method IIIC

C  Method       : Inputs, state variables and algebraic variables are included
C                 in the vector of variables. The right hand side components
C                 are a combination of function values, time step and previous
C                 values. Input components are not affected.

C  Call         : CALL LOBATT(NREG,DT,VPREG,VLREG,VDREG,MSTAT,MPIREG,IREG,
C                >            MPRREG,RREG,DELAY,IERR,ISTAT,ZCTRL0)

C  Arguments    : NDIM   - integer, dimension of the variable vector
C                 DT     - real, time step
C                 VPREG  - real array(NDIM/2), the previous values of the
C                          variables
C                 VLREG  - real array(3*NDIM/2), on input: the most recent
C                          iteration value using this integration method, on
C                          output: the right hand side of the numerical
C                          equation. The last third is a help array.
C                 VDREG  - real array(NDIM/2), function values
C                 MSTAT  - integer array(NDIM/2), status flags for variables
C                 MPIREG - integer array(*), pointers to the content of IREG,
C                          see calling routine
C                 IREG   - integer array(*), working area for integer variab-
C                          les, see calling routine
C                 MPRREG - integer array(*), pointers to the content of RREG,
C                          see calling routine
C                 RREG   - real array(*), working area for real variables,
C                          see calling routine
C                 DELAY  - real array(*), time history of delayed variables
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning
C                 ISTAT  - integer array(10), statistics vector
C                          ISTAT(9) = no. of iterations with Lobatto IIIC

C  Input        : NDIM,DT,VPREG,VLREG,VDREG,MSTAT,ISTAT

C  Output       : VLREG,IERR,ISTAT

C  Work area    : MPIREG,IREG,MPRREG,RREG

C  Common       : No common

C  Restriction  : Functions and subroutines called:
C                  FCN - for function evaluation

C  Fault react. : None

C  Note         : Nothing

C **********************************************************************

      SUBROUTINE LOBATT(NDIM, DT, VPREG, VLREG, VDREG, MSTAT, MPIREG,
     >                  IREG, MPRREG, RREG, DELAY, IERR, ISTAT,
     >                  ZCTRL0)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VPREG(NDIM/2), VLREG(3*NDIM/2), VDREG(NDIM/2),
     >                 RREG(*), DELAY(*)
      INTEGER MSTAT(NDIM/2), MPIREG(*), IREG(*), MPRREG(*),
     >        IERR(2), ISTAT(10)
      LOGICAL ZCTRL0

C The right hand side of the numerical equation for Lobatto IIIC is:

C     0 = -U(1) + s(t(n))
C     0 = -X(1) + x(n) + 0.5*h*(f(U(1),X(1),Z(1)) - f(U(2),X(2),Z(2)))
C     0 =  g(U(1),X(1),Z(1))
C     0 = -U(2) + s(t(n+1))
C     0 = -X(2) + x(n) + 0.5*h*(f(U(1),X(1),Z(1)) + f(U(2),X(2),Z(2)))
C     0 =  g(U(2),X(2),Z(2))

C ............................. The dimension of the numerical equation is
C ............................. twice the dimension of the system equation
      NREG = NDIM/2

C ............................. Save maximum end point of next time step
      TNM = RREG(2)
C ............................. Invoke the function routine to obtain the
C ............................. first part of the function values
      CALL FCN(IREG(1), RREG(1), RREG(2), NREG, IREG(2), IREG(5),
     >     VLREG, VDREG, RREG(MPRREG(2)), DELAY, RREG(MPRREG(7)),
     >     MSTAT, IREG(MPIREG(3)), IREG(MPIREG(6)), IREG(MPIREG(8)),
     >     IREG(MPIREG(7)), RREG(MPRREG(12)), IERR, IREG(3), ISTAT(2),
     >     ZCTRL0)

C ............................. Restore maximum end point of next time step
      RREG(2) = TNM
C ............................. Check the error flag.
      IF (IERR(1) .LT. 0) RETURN

C ............................. For each variable ...
      DO 100 N = 1,NREG
C ............................. identify variable status
         ITYP = IABS(MSTAT(N))
C ............................. When the variable is an input,
         IF (ITYP .EQ. 1) THEN
C ............................. set rhs to zero
            VLREG(N) = 0.0d0
C ............................. When the variable is a state variable,
         ELSE IF (ITYP .EQ. 2) THEN
C ............................. compute the first part of the rhs expression
            VLREG(N) = VPREG(N) - VLREG(N) + 0.5d0*DT*VDREG(N)
            VLREG(NDIM+N) = VPREG(N) -VLREG(NREG+N) +0.5d0*DT*VDREG(N)
C ............................. When the variable is an algebraic variable,
         ELSE IF (ITYP .EQ. 3) THEN
C ............................. assign the functional value to the rhs
            VLREG(N) = VDREG(N)
         END IF
  100 CONTINUE

C ............................. Invoke the function routine to obtain the
C ............................. second part of the function values
      CALL FCN(IREG(1), RREG(1), RREG(2), NREG, IREG(2), IREG(5),
     >     VLREG(NREG+1), VDREG, RREG(MPRREG(2)), DELAY,
     >     RREG(MPRREG(7)), MSTAT, IREG(MPIREG(3)), IREG(MPIREG(6)),
     >     IREG(MPIREG(8)), IREG(MPIREG(7)), RREG(MPRREG(12)), IERR,
     >     IREG(3), ISTAT(2), ZCTRL0)

C ............................. Check the error flag
      IF (IERR(1) .LT. 0) RETURN

C ............................. For each variable ...
      DO 200 N = 1,NREG
C ............................. identify variable status
         ITYP = IABS(MSTAT(N))
C ............................. When the variable is an input,
         IF (ITYP .EQ. 1) THEN
C ............................. set rhs to zero
            VLREG(NREG+N) = 0.0d0
C ............................. When the variable is a state variable,
         ELSE IF (ITYP .EQ. 2) THEN
C ............................. compute the first part of the rhs expression
            VLREG(N) = VLREG(N) - 0.5d0*DT*VDREG(N)
            VLREG(NREG+N) = VLREG(NDIM+N) + 0.5d0*DT*VDREG(N)
C ............................. When the variable is an algebraic variable,
         ELSE IF (ITYP .EQ. 3) THEN
C ............................. assign the functional value to the rhs
            VLREG(NREG+N) = VDREG(N)
         END IF
  200 CONTINUE

C ............................. Update the no. of iterations with Loatto IIIC
      ISTAT(9) = ISTAT(9) + 1

      RETURN
      END


C **********************************************************************

C               ***     *        ***    *   *
C               *   *   *       *   *    * *
C               *   *   *       *****     *
C               *   *   *       *   *     *
C               ***     *****   *   *     *

C **********************************************************************

C  Name         : DLAY

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : Mars 1992

C  Last revision: No revision

C  Purpose      : Update delay buffer.

C  Method       : For each time step the value of time and delay variables
C                 are stored in one cluster and put on top of the delay buffer,
C                 shifting the old clusters backwards. Each variable require
C                 one buffer location per level of the RK method.
C                 Alternatively (LEV > 0) only the last values for level LEV
C                 are inserted on top of the buffer.

C  Call         : CALL DLAY(NTD,T,DT,TNM,NMOD,VREG,MPRPAR,RPAR,MPMTOP,MMTOP
C                >          M,LEV,DELAY,NDDIM,IERR)

C  Arguments    : NTD    - number of time delays
C                 T      - real, end point of the integration step
C                 DT     - real, actual time step
C                 TNM    - real, max. end point of the integration step
C                 NMOD   - integer, no. of modules
C                 VREG   - real array(NREG), the most recent accepted values
C                          of the control state variables
C                 MPRPAR(NMOD+1) - pointers to the content of RPAR
C                 RPAR(*) - module parameters
C                 MPMTOP(NMOD+1) - pointers to the content of MMTOP
C                 MMTOP(*) - topology vector
C                 M      - number of levels of the RK method
C                 LEV    - integer, update mode
C                        < 0 : re-entry after DELAY buffer reallocation
C                        = 0 : insert initial values and shift downwards
C                        > 0 : insert last values of level LEV on top of buffer
C                 DELAY(NDDIM) - real array, time history of delay variables
C                 NDDIM  - integer, delay array dimension
C                 IERR   - integer array(2), error flag
C                        < 0 : error has occurred
C                        = 0 : no error or warning
C                        > 0 : warning

C  Input        : NTD,T,DT,TNM,NMOD,VREG,MPRPAR,RPAR,MPMTOP,MMTOP,M,LEV
C                 DELAY,NDDIM

C  Output       : DELAY,TNM,IERR

C  Fault react. : Return with non-zero error flag

C **********************************************************************

      SUBROUTINE DLAY(NTD, T, DT, TNM, NMOD, VREG, MPRPAR, RPAR,
     >                MPMTOP, MMTOP, M, LEV, DELAY, NDDIM, IERR)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VREG(*), RPAR(*), DELAY(*), DBUF(9)
      INTEGER MPRPAR(*), MPMTOP(*), MMTOP(*), IERR(2)
      SAVE NDBUF, NTIME, DLYMIN, DLYMAX

C ---> 08.03.29/KMO ........... Check for retry after DELAY buffer reallocation
      if (lev .lt. 0) goto 145
C ............................. Initialize parameters
      NDBUF = M*NTD + 1
      NDLAY = 2 - M
      DLYMIN = DT
      DLYMAX = DLYMIN
      NTIME = 0
      DO 140 N = 1,NMOD
C ............................. Invoke each module and check the module type
         MTYP = IDINT(RPAR(MPRPAR(N)))
         IF (MTYP .NE. 11) GOTO 140

C ............................. Indices for delay buffer and variables vector
         NDLAY = NDLAY + M
         NVAR = MMTOP(MPMTOP(N))
C ............................. Eventually update only the last values
         IF (LEV .NE. 0) THEN
            DELAY(NDLAY + LEV - 1) = VREG(NVAR)
            GOTO 140
         END IF
C ............................. Put the delayed value in help buffers
C ............................. One value for each level in the RKM is required
         DO 100 MM = 1,9
            DBUF(MM) = VREG(NVAR)
  100    CONTINUE
C ............................. Save the shortest and longest time delays
         TDLAY = RPAR(MPRPAR(N)+1)
         DLYMIN = DMIN1(DLYMIN, TDLAY)
         DLYMAX = DMAX1(DLYMAX, TDLAY)

         NN = NDLAY - NDBUF
C ............................. Put the latest value on top of the buffer and
C ............................. shift the corresponding values backwards as far
C ............................. as necessary considering the delay magnitude.
  110    NN = NN + NDBUF
         DO 120 MM = 1,M
            BUF = DELAY(NN+MM-1)
            DELAY(NN+MM-1) = DBUF(MM)
            DBUF(MM) = BUF
  120    CONTINUE
         NOW = NN - NDLAY + 1 - NDBUF
         if (now .lt. 1) goto 110
         IF (DELAY(NOW) .GT. MAX(T-DLYMAX,DELAY(NOW+NDBUF))) GOTO 110
C ............................. Finish the shift of values
         DO 130 MM = 1,M
C ---> 02.01.16/TI ............ -1 included in indexing
            DELAY(NN + NDBUF + MM - 1) = DBUF(MM)
  130    CONTINUE
         NTIME = NOW + NDBUF
  140 CONTINUE
C
      IF (LEV .NE. 0) RETURN
C ............................. Check required space against the buffer length
  145 CONTINUE
      IF (NTIME + 3*NDBUF .GT. NDDIM) THEN
C ---> 08.03.29/KMO ........... Delay buffer too short, set IERR(2) = min size
C ............................. and return to calling routine for reallocation.
C ............................. Note: Had to add the 3* factor to avoid crash
C ............................. due to overwriting in the loops 110-130 above.
         IERR(1) = 999
         IERR(2) = NTIME + 3*NDBUF
         RETURN
      END IF
C ............................. Zero the next time value in order to indicate
C ............................. the working space.
      IF (T-DELAY(NTIME) .GT. DLYMAX) THEN
         DELAY(NTIME + NDBUF) = 0.0d0
      END IF
  150 CONTINUE
C ............................. Shift the time values correspondingly
      NTIME = NTIME - NDBUF
      IF (DELAY(NTIME) .LT. T) DELAY(NTIME + NDBUF) = DELAY(NTIME)
      IF (NTIME .GT. NDBUF) GOTO 150
      DELAY(1) = T
C ............................. Check max. end point of next integration step
      IF (DLYMIN .LT. DT) TNM = T - DT + DLYMIN

      RETURN
      END


C **********************************************************************
C
C                ***    *   *   *****   *   *   *****    ****
C               *   *   *   *     *     ** **   *       *
C               *   *   *   *     *     * * *   ****     ***
C               *   *   *   *     *     *   *   *           *
C                ***     ***      *     *   *   *****   ****
C
C **********************************************************************
C
C  Name         : OUTMES
C
C  Type         : Subroutine
C
C  Language     : FORTRAN-77
C
C  Computer     : VAX/VMS
C
C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control
C
C  Date         : November 1988
C
C  Last revision: No revision
C
C  Purpose      : Print out error messages and statistics
C
C  Method       : Company security
C
C  Call         : CALL OUTMES(T,IERR,ISTAT)
C
C  Arguments    : T     - real, integration time
C                 IERR  - integer array(2), error flag
C                       IERR(1) = -7 : discrete event not found
C                       IERR(1) = -6 : truncation error limit exceeded
C                       IERR(1) = -5 : divergence in Lobatto IIIC
C                       IERR(1) = -4 : divergence in Backward Euler
C                       IERR(1) = -3 : divergence in steady state
C                       IERR(1) = -2 : singular iteration matrix
C                       IERR(1) = -1 : non-existent module number
C                       IERR(1) =  0 : no error
C                       IERR(1) =  1 : not used
C                       IERR(1) =  2 : not used
C                       IERR(1) =  3 : slow convergence in steady state
C                       IERR(1) =  4 : slow convergence in Backward Euler
C                       IERR(1) =  5 : slow convergence in Lobatto IIIC
C ---> 92.04.24/TI      IERR(1) =  6 : truncation error limit exceeded
C                 ISTAT - integer array(10), statistics vector
C                       1 : function evaluations for initialization
C                       2 : function eval. for steady state and integration
C                       3 : Jacobian evaluations
C                       4 : LU-factorizations for steady state
C                       5 : LU-factorizations for Backward Euler
C                       6 : LU-factorizations for Lobatto IIIC
C                       7 : iterations for steady state
C                       8 : iterations in Backward Euler
C                       9 : iterations in Lobatto IIIC
C                      10 : times solving a linear system
C
C  Input        : IERR,ISTAT
C
C  Output       : None
C
C  Work area    : None
C
C  Common       : No common
C
C  Restriction  : Functions and subroutines called: None
C
C  Note         : Nothing
C
C  Fault react. : Error messages printout
C
C **********************************************************************

      SUBROUTINE OUTMES(T, IERR, ISTAT)
CDEC$ ATTRIBUTES DLLEXPORT :: OUTMES

      use ioModule, only : IO

      implicit none

      DOUBLE PRECISION T
      INTEGER IERR(2), ISTAT(10)

C     ..................................... Print error messages
      GOTO (110,109,107,105,103,102,101,200,200,200,104,106,108,109)
     >     8 + IERR(1)

  101 WRITE(IO,600) T
      WRITE(IO,610) 'NON-EXISTENT MODULE TYPE NO.',
     >     IERR(2)
      GOTO 200

  102 WRITE(IO,600) T
      WRITE(IO,610) 'ILL-CONDITIONED NEWTON MATRIX, DIAG. ELEMENT NO.',
     >     IERR(2)
      GOTO 200

  103 WRITE(IO,600) T
      WRITE(IO,610) 'DIVERGENCE FOR STEADY STATE, MAX IN VARIABLE NO.',
     >     IERR(2)
      GOTO 200

  104 WRITE(IO,600) T
      WRITE(IO,610) 'SLOW CONVERGENCE FOR STEADY STATE, MAX IN VAR. NO.'
     >     ,IERR(2)
      GOTO 200

  105 WRITE(IO,600) T
      WRITE(IO,610) 'DIVERGENCE IN BACKWARD EULER, MAX IN VARIABLE NO.',
     >     IERR(2)
      GOTO 200

  106 WRITE(IO,600) T
      WRITE(IO,610) 'SLOW CONVERGENCE IN BACKWARD EULER, MAX IN VAR. NO.
     >     ',IERR(2)
      GOTO 200

  107 WRITE(IO,600) T
      WRITE(IO,610) 'DIVERGENCE IN LOBATTO IIIC, MAX IN VARIABLE NO.',
     >     IERR(2)
      GOTO 200

  108 WRITE(IO,600) T
      WRITE(IO,610) 'SLOW CONVERGENCE IN LOBATTO IIIC, MAX IN VAR. NO.',
     >     IERR(2)
      GOTO 200

  109 WRITE(IO,600) T
      WRITE(IO,610) 'Truncation error limit exceeded by variable no.',
     >     IERR(2)
      GOTO 200

  110 WRITE(IO,601) T
      WRITE(IO,610) 'Failed to find discrete event time for var. no.',
     >     IERR(2)

  200 CONTINUE

C      ..................................... Print statistics
C      WRITE(IO,660) (ISTAT(I),I=1,10)
      WRITE(IO,690) ISTAT(2),ISTAT(3),ISTAT(5),ISTAT(10)

      RETURN

  600 FORMAT(' *** Error return from RSFED: Time =',1PE12.5)
  601 FORMAT('  ** Warning from RSFED: Time =',1PE12.5)
  610 FORMAT(5X,A,I6)
C  660 FORMAT('   ITIN  =',I4,'    NFCN  =',I4,'    JAC   =',I4
C     >      /'   LUSS  =',I4,'    LUBE  =',I4,'    LULO  =',I4
C     >      /'   ITSS  =',I4,'    ITBE  =',I4,'    ITLO  =',I4
C     >      /'   NSOLV =',I4)
  690 FORMAT(/2X,39('*')
     >       /I6,' FUNCTION EVALUATIONS'
     >       /I6,' JACOBI EVALUATIONS'
     >       /I6,' LU-FACTORIZATIONS'
     >       /I6,' ITERATION'
     >       /2X,39('*'))

      END
