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
C  mod1.f  - Comparator
C  mod2.f  - Adder
C  mod3.f  - Amplifier
C  mod4.f  - Integrator
C  mod5.f  - Limited derivator
C  mod6.f  - Multiplier
C  mod11.f - Time delay
C  mod12.f - Sample and hold
C  mod21.f - Logical switch
C  mod22.f - Limitation
C  mod23.f - Dead zone
C  mod24.f - Hysteresis
C  mod31.f - PI-controller
C  mod32.f - Proportional + limited Integral controller
C  mod33.f - PD-controller
C  mod34.f - Proportional + limited Derivative controller
C  mod35.f - PID-controller
C  mod36.f - Proportional + Integral + limited Derivative controller
C  mod37.f - Proportional + limited Integral + limited Derivative controller
C  mod41.f - Real pole
C  mod42.f - Complex conjugate poles
C  mod43.f - First order element
C  mod44.f - Second order element
C  mod45.f - Linear control system
C  mod46.f - General transfer function

C **********************************************************************

C               *   *    ***    ****
C               ** **   *   *   *   *
C               * * *   *   *   *   *   *   *   *   *
C               *   *   *   *   *   *     *       *
C               *   *    ***    ****    *   *   *   *

C **********************************************************************

C  Name         : MODxx

C  Type         : Subroutine

C  Language     : FORTRAN-77

C  Computer     : VAX/VMS

C  Implementor  : Torleif Iversen, SINTEF/Div. of Aut. Control

C  Date         : October 1988

C  Last revision: No revision

C  Purpose      : Compute functional values for module no. xx, or initialize
C                 the status flags MST1,...,MSTm for the variables, and the
C                 module assignment vector elements MM1,...,MMm.

C  Method       : The equations are written on an explicit form which means
C                 that a functional value is an expression for the derivative
C                 in case of state equations or zero in case of algebraic
C                 equations.

C  Call         : CALL MODxx(IOP,MODNO,Y1,..,Ym,YD1,..,YDm,MST1,..,MSTm,
C                >           MM1,...,MMm,PAR)

C  Arguments    : IOP   - integer, operation mode
C                       = 1 : initalization in order to set the status flags
C                      ^= 1 : functional evaluations will be carried out
C                 MODNO - integer, module number in configuration
C                 Yi    - real, i=<1,m>, the most recent value for variable i
C                 YDi   - real, i=<1,m>, the functional value for variable i
C                 MSTi  - integer, i=<1,m>, status flag for each variable.
C                         Dependent on the actual value, the following changes
C                         will take place during initialization:
C                       = 0 : change to 2 or 3 to indicate state- or algebraic
C                             variable respectively,
C                       = 1 : input, not changed
C                       = 2 : state variable, not changed
C                       = 3 : algebraic variable, may eventually be changed to
C                             a state variable
C                         During a functional evaluation the status flags are
C                         not used.
C                 MMi   - integer, i=<1,m>, module assignment elements.
C                         During initialization the module number will be ass-
C                         igned to MMi if the value of MSTi is changed.
C                         During execution MMi indicates where the functional
C                         values should be placed.
C                 PAR   - real array(*), parameter values

C  Input        : IOP,Y1,...,Ym,MST1,...,MSTm,PAR

C  Output       : YD1,...,YDm,MST1,...,MSTm (when IOP = 1)

C  Work area    : m

C  Common       : No common

C  Restriction  : Functions and subroutines called: None

C  Fault react. : None

C  Note         : This is a template documentation


C **********************************************************************
C ******        GROUP 0 : BASIC ELEMENTS                          ******
C **********************************************************************

      SUBROUTINE MOD1(IOP, MODNO, Y1, Y2, Y3, YD1, YD2, YD3, MST1,
     >                MST2, MST3, MM1, MM2, MM3, PAR, ZCTRL0)

C ............................. Comparator
C ............................. y3 = k*(y1 - y2)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Count variables where status is not set
         NOTSET = 0
         IF (MST1 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
         IF (MST2 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
         IF (MST3 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
C ............................. Unless more than one variable status remains,
C ............................. set the missing variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         IF (NOTSET .GT. 1) RETURN
         IF (MST1 .EQ. 0) THEN
            MST1 = 3
            MM1  = MODNO
            if (ZCTRL0) return
            Y1 = Y2 + Y3/PAR
         ELSE IF (MST2 .EQ. 0) THEN
            MST2 = 3
            MM2  = MODNO
            if (ZCTRL0) return
            Y2 = Y1 - Y3/PAR
         ELSE IF (MST3 .EQ. 0) THEN
            MST3 = 3
            MM3  = MODNO
            if (ZCTRL0) return
            Y3 = PAR*(Y1 - Y2)
         END IF
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output
         G = Y3 - PAR*(Y1 - Y2)
         IF (MM1 .EQ. MODNO) THEN
            YD1 = G/PAR
         ELSE IF (MM2 .EQ. MODNO) THEN
            YD2 = - G/PAR
         ELSE IF (MM3 .EQ. MODNO) THEN
            YD3 = - G
         END IF
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD2(IOP, MODNO, Y1, Y2, Y3, YD1, YD2, YD3, MST1,
     >                MST2, MST3, MM1, MM2, MM3, PAR, ZCTRL0)

C ............................. Adder
C ............................. y3 = k*(y1 + y2)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Count variables where status is not set
         NOTSET = 0
         IF (MST1 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
         IF (MST2 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
         IF (MST3 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
C ............................. Unless more than one variable status remains,
C ............................. set the missing variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         IF (NOTSET .GT. 1) RETURN
         IF (MST1 .EQ. 0) THEN
            MST1 = 3
            MM1  = MODNO
            if (ZCTRL0) return
            Y1 = Y3/PAR - Y2
         ELSE IF (MST2 .EQ. 0) THEN
            MST2 = 3
            MM2  = MODNO
            if (ZCTRL0) return
            Y2 = Y3/PAR - Y1
         ELSE IF (MST3 .EQ. 0) THEN
            MST3 = 3
            MM3  = MODNO
            if (ZCTRL0) return
            Y3 = PAR*(Y1 + Y2)
         END IF
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output
         G = Y3 - PAR*(Y1 + Y2)
         IF (MM1 .EQ. MODNO) THEN
            YD1 = G/PAR
         ELSE IF (MM2 .EQ. MODNO) THEN
            YD2 = G/PAR
         ELSE IF (MM3 .EQ. MODNO) THEN
            YD3 = - G
         END IF
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD3(IOP, MODNO, Y1, Y2, YD2, MST1, MST2, MM2,
     >                PAR, ZCTRL0)

C ............................. Amplifier
C ............................. 0 = y2 - k*y1

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         IF (MST1 .EQ. 0) THEN
            RETURN
         ELSE
            MST2 = 3
            MM2  = MODNO
            if (ZCTRL0) return
            Y2 = PAR*Y1
         END IF
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output
         YD2 = PAR*Y1 - Y2
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD4(IOP, MODNO, Y1, Y2, YD2, MST2, MM2, PAR)

C ............................. Integrator
C ............................. y2' = y1

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the output variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         MST2 = 2
         MM2  = MODNO
         Y2 = PAR
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output
         YD2 = Y1
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD5(IOP, MODNO, Y1, Y2, Y3, YD2, YD3,
     >                MST2, MST3, MM2, MM3, PAR, ZCTRL0)

C ............................. Limited derivator
C ............................. y2' = y3
C .............................  0  = Tu*y3 - y1 + y2

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set missing variable status and module number
C ---> 02.01.02/TI ............ (Zero initial state values assumed)
         MST2 = 2
         MST3 = 3
         MM2  = MODNO
         MM3  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = Y1/PAR
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output
         YD2 = Y3
         YD3 = (Y1-Y2)/PAR - Y3
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD6(IOP, MODNO, Y1, Y2, Y3, YD1, YD2, YD3,
     >                MST1, MST2, MST3, MM1, MM2, MM3, ZCTRL0)

C ............................. Multiplier
C ............................. 0 = y3 - y1*y2
C ---> 02.01.02/TI ............ Module rewritten to allow for division

      use ioModule, only : IO

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL ZCTRL0

      eps = 1000.0D0*tiny(0.0D0)
C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ---> 02.01.02/TI ............ Count variables where status is not set
         NOTSET = 0
         IF (MST1 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
         IF (MST2 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
         IF (MST3 .EQ. 0) THEN
            NOTSET = NOTSET + 1
         END IF
C ............................. Unless more than one variable status remains,
C ............................. set the missing variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         IF (NOTSET .GT. 1) RETURN
         IF (MST1 .EQ. 0) THEN
            MST1 = 3
            MM1  = MODNO
            if (ZCTRL0) return
            if ( Y2 < - eps .or. Y2 > eps ) then
               Y1 = Y3/Y2
            else
               write(IO,600) Y2
               Y1 = 0.0D0
            end if
         ELSE IF (MST2 .EQ. 0) THEN
            MST2 = 3
            MM2  = MODNO
            if (ZCTRL0) return
            if ( Y1 < - eps .or. Y1 > eps ) then
               Y2 = Y3/Y1
            else
               write(IO,600) Y1
               Y2 = 0.0D0
            end if
         ELSE IF (MST3 .EQ. 0) THEN
            MST3 = 3
            MM3  = MODNO
            if (ZCTRL0) return
            Y3 = Y1*Y2
         END IF
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output
         G = Y3 - Y1*Y2
         IF (MM1 .EQ. MODNO) THEN
            if ( Y2 < - eps .or. Y2 > eps ) then
               YD1 = G/Y2
            else
               write(IO,600) Y2
               YD1 = 0.0D0 - Y1
            end if
         ELSE IF (MM2 .EQ. MODNO) THEN
            if ( Y1 < - eps .or. Y1 > eps ) then
               YD2 = G/Y1
            else
               write(IO,600) Y1
               YD2 = 0.0D0 - Y2
            end if
         ELSE IF (MM3 .EQ. MODNO) THEN
            YD3 = - G
         END IF
      END IF

      RETURN
 600  FORMAT('  ** WARNING: Multiplier/division mode:',
     >       ' Divisor close zero',1PE12.5)
      END


C **********************************************************************

      SUBROUTINE MOD7(IOP, MODNO, Y1, Y2, YD2, MST1, MST2, MM2,
     >                PAR, ZCTRL0)

C ............................. Power function
C ............................. 0 = y2 - y1^par
C ---> 02.01.02/TI ............ Completely rewritten module

      use ioModule, only : IO

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL ZCTRL0

      eps = 10.0D0*tiny(0.0D0)

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         IF (MST1 .EQ. 0) RETURN
         MST2 = 3
         MM2  = MODNO
         if (ZCTRL0) return
         if ( Y1 >= 0.0D0 ) then
            if ( Y1 > eps .or. PAR > 0.0D0 ) then
               Y2 = Y1**PAR
            else
               write(IO,602) Y1
               Y2 = 0.0D0
            end if
         else
            if ( dabs(PAR-anint(PAR)) > 1.0D-6 .and. IOP == 2 ) then
               write(IO,601) PAR,anint(PAR)
            end if
            if ( Y1 < -eps .or. PAR > 0.0D0 ) then
               Y2 = Y1**anint(PAR)
            else
               write(IO,602) Y1
               Y2 = 0.0D0
            end if
         end if
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output

C ............................. Check sign of Y1 and PAR, and closeness
C ............................. to zero of Y1 when PAR is negative
         if ( Y1 >= 0.0D0 ) then
            if ( Y1 > eps .or. PAR > 0.0D0 ) then
               YD2 = Y1**PAR - Y2
            else
               write(IO,602) Y1
C              YD2 = eps**PAR - Y2
               YD2 = 0.0D0 - Y2
            end if
         else
            if ( dabs(PAR-anint(PAR)) > 1.0D-6 .and. IOP == 2 ) then
               write(IO,601) PAR,anint(PAR)
            end if
            if ( Y1 < -eps .or. PAR > 0.0D0 ) then
               YD2 = Y1**anint(PAR) - Y2
            else
               write(IO,602) Y1
C              YD2 = (-eps)**anint(PAR) - Y2
               YD2 = 0.0D0 - Y2
            end if
         end if
      END IF

      RETURN
 601  FORMAT('  ** WARNING: Power function: Negative argument.',
     >       ' Decimal exponent',1PE12.5,' made integer',I6)
 602  FORMAT('  ** WARNING: Power function: Inversion of value',
     >       ' close to zero',1PE12.5)
      END


C **********************************************************************
C ******        GROUP 1 : TIME DEPENDENT ELEMENTS                  ******
C **********************************************************************

      SUBROUTINE MOD11(IOP, MODNO, T, NDLAY, NTD, DELAY, Y1, Y2, YD2,
     >                 MST1, MST2, MM2, PAR, ZCTRL0)

C ............................. Time delay
C ............................. 0 = y2(t) - y1(t-tdlay)

C Parameters :
C     t     - simulation time
C     ndlay - delay element number in configuration
C     ntd   - number of delay elements in configuration
C     delay - time history of delayed variables
C     PAR   =  TDLAY- time delay value

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR, DELAY(*)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         IF (MST1 .EQ. 0) RETURN
         MST2 = 3
         MM2  = MODNO
         if (.not.ZCTRL0) Y2 = Y1
C ............................. Increase the number of delays and insert into
C ............................. delay buffer delayed time and initial state
         NTD = NDLAY
C        IF (DELAY(1) .LT. PAR) DELAY(1) = PAR
         DELAY(2*NDLAY) = Y1
         DELAY(2*NDLAY+1) = Y1
C ............................. Execution mode
      ELSE
C ............................. Set length of delay buffer per time step
         NDBUF = 2*NTD + 1
C ............................. Search down the delay buffer to find the
C ............................. time point for the delayed variable
         N = 1
  100    N = N + NDBUF
         IF (DELAY(N) .GT. T - PAR .AND. DELAY(N+NDBUF) .LT. DELAY(N))
     >        GOTO 100
C ............................. Save the upper end point of the time interval
         IF (DELAY(N) .LE. T - PAR) N = N - NDBUF
C ............................. Compute the relative distance from this point
         IF (DELAY(N) .LE. DELAY(N+NDBUF)) THEN
            TAU = 0.0D0
         ELSE
            TAU = (DELAY(N) - T + PAR) / (DELAY(N) - DELAY(N+NDBUF))
         END IF
C ............................. Get variable values for interpolation
         YLEV1 = DELAY(N + 2*NDLAY - 1)
         YLEV2 = DELAY(N + 2*NDLAY)
         YLEVP = DELAY(N + 2*NDLAY + NDBUF)
C ............................. Update output using interpolated value
C ............................. NB! METHOD DEPENDENT INTERPOLATION !!!
         YD2 = (1-TAU)*YLEV2 + TAU*(1-TAU)*YLEV1 + TAU*TAU*YLEVP - Y2
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD12(IOP, MODNO, T, TNM, Y1P, Y1, Y2, YD2,
     >                 MST1, MST2, MM2, PAR, ZCTRL0)

C ............................. Sample and hold
C ............................. 0 = y1(t') - y2(t),  t' = int(t/T)*T

C Parameters :
C     t    - simulation time
C     tnm  - maximum end point for next integration step
C     y1p  - previous input value
C     PAR(1) =  T    - sample period
C     PAR(2) =  m    - multiplicity of basic sample periods
C     PAR(3) =  phd  - phase delay in basic sample periods
C     PAR(4) =  tp   - previous simulation time
C     PAR(5) =  ts   - sample time
C     PAR(6) =  ys   - sampled value

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(6)
      LOGICAL ZCTRL0

      eps = 1000.0D0*tiny(0.0D0)
      tn = 0.0D0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number,
C ---> 02.01.02/TI ............ initial output value and initialise parameters
         IF (MST1 .EQ. 0) RETURN
         MST2 = 3
         MM2  = MODNO
         if (.not.ZCTRL0) Y2 = Y1
         IF (PAR(1) .LE. 0.0d0) PAR(1) = TNM - T
C TI 02.10.01 ................. The next two parameters are not really needed
         PAR(2) = 1.0D0
         PAR(3) = 0.0D0
         PAR(4) = T
         PAR(5) = T
         IF (PAR(3) .GE. 1.0)
     >        PAR(5) = T + (IDINT(PAR(3)) - IDINT(PAR(2)))*PAR(1)
         PAR(6) = Y1
         IF (PAR(2) .LT. 1.0) PAR(2) = 1.0D0
         TNM = DMIN1(TNM, T + PAR(1))

C ............................. Execution mode
      else if (par(1) .lt. eps) then
C ............................. No sample-and-hold at all
         yd2 = y1 - y2
      ELSE
C ............................. Check for new time step, reduced time
C ............................. step or iteration
         RLIM = PAR(1) * 1.0D-8
C ............................. Iteration
         IF (DABS(T - PAR(4)) .LT. RLIM) GOTO 100
C ............................. Reduced time step
         IF (T .LT.  PAR(4)) THEN
            IF (T .LT. PAR(5)) PAR(5) = PAR(5) - IDINT(PAR(2))*PAR(1)
         ELSE
C ............................. When the simulation time is advanced ...
C ............................. when the time step start from the sample time
C ............................. update sampled value
            IF (DABS(PAR(4) - PAR(5)) .LT. RLIM) PAR(6) = Y1P
         END IF
C ............................. Update previous simulation time
         PAR(4) = T

C ............................. When the simulation time is in the the next
C ............................. sample time, update the sample time and the
C ............................. most recent input value
         TN = PAR(5) + IDINT(PAR(2))*PAR(1)
C TI 02.10.01 ................. Next statement in order to handle sample
C ............................. intervals shorter than the time step
         IF (TN .LT. T) TN = T
         IF (DABS(T - TN) .LT. RLIM) THEN
C ............................. use the previous output value
            PAR(5) = T
            TN = T + PAR(1)
         END IF

C ............................. Update maximum end point for the next time
C ............................. step and update the output
  100    CONTINUE
         IF (TNM .GT. TN) TNM = TN
         YD2 = PAR(6) - Y2
      END IF

      RETURN
      END


C **********************************************************************
C ******        GROUP 2 : PIECEWISE CONTINUOUS FUNCTIONS          ******
C **********************************************************************

      SUBROUTINE MOD21(IOP, MODNO, LIM, TAU, TOL, Y1P, Y1, Y2, YD2,
     >                 MST1, MST2, MM2, PAR, ZCTRL0)

C ............................. Logical switch
C ............................. 0 = ylow - y2  for y1 < yon
C ............................. 0 = yup  - y2  for y1 > yon

C Parameters:
C     LIM - limit cross indicator
C           = -1 : decreasing crossing
C           =  0 : no crossing
C           = +1 : increasing crossing
C     TOL - limit cross tolerance
C     Y1p - previous input value
C     PAR(1) = Yon  - limit cross value for input
C     PAR(2) = Ylow - low value for output
C     PAR(3) = Yup  - high value for output

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(3)
      LOGICAL ZCTRL0

      LIM = 0
C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number,
C ---> 02.01.02/TI ............ initial output value and initialise parameters
         IF (MST1 .EQ. 0) RETURN
         MST2 = 3
         MM2  = MODNO
         if (ZCTRL0) return
         IF (Y1 .LE. PAR(1)) THEN
            Y2 = PAR(2)
         ELSE
            Y2 = PAR(3)
         END IF
C ............................. Execution mode
      ELSE
C ............................. Compute functional value ...
C ............................. and assign it to the proper output

C ............................. If previous value is below the cross band
         IF (Y1P .LT. PAR(1) - TOL) THEN
            YD2 = PAR(2) - Y2
C ............................. Check for positive crossing
            IF (Y1 .GT. PAR(1) + TOL) LIM = 1
C ............................. If previous value is above the cross band
         ELSE IF (Y1P .GT. PAR(1) + TOL) THEN
            YD2 = PAR(3) - Y2
C ............................. Check for negative crossing
            IF (Y1 .LT. PAR(1) - TOL) LIM = -1
C ............................. When the previous value is in the cross band
         ELSE
C ............................. Determine output without regarding tolerance
            IF (Y1 .LT. PAR(1)) THEN
               YD2 = PAR(2) - Y2
            ELSE
               YD2 = PAR(3) - Y2
            END IF
         END IF
C ---> 02.02.10/TI ............ Find time point of eventual discontinuity
         IF ( LIM .NE. 0 ) TAU = (PAR(1)-Y1P)/(Y1-Y1P)
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD22(IOP, MODNO, LIM, TAU, TOL, Y1P, Y1, Y2, YD2,
     >                 MST1, MST2, MM2, PAR, ZCTRL0)

C ............................. Limitation
C ............................. 0 = ylow - y2  for y1 < ylow
C ............................. 0 = y1 - y2  for ylow < y1 < yup
C ............................. 0 = yup  - y2  for y1 > yup

C Parameters:
C     LIM - limit cross indicator
C           = -2 : decreasing crossing through upper limit
C           = -1 : decreasing crossing through lower limit
C           =  0 : no crossing
C           = +1 : increasing crossing through lower limit
C           = +2 : increasing crossing through upper limit
C     TOL - limit cross tolerance
C     Y1p - previous input value
C     PAR(1) = Ylow - lower limit
C     PAR(2) = Yup  - upper limit

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(2)
      LOGICAL ZCTRL0

      LIM = 0
C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number
C ---> 02.01.02/TI ............ and initial output value
         IF (MST1 .EQ. 0) RETURN
         MST2 = 3
         MM2  = MODNO
         if (ZCTRL0) return
         Y2 = Y1
         IF (Y1 .LT. PAR(1)) THEN
            Y2 = PAR(1)
         ELSE IF (Y1 .GT. PAR(2)) THEN
            Y2 = PAR(2)
         END IF
C ............................. Execution mode
      ELSE
C ............................. Assign functional value to the proper output
C ............................. In the area of no effect ...
         YD2 = Y1 - Y2
C ............................. For values below lower limit
         IF ( Y1P .LT. PAR(1) + TOL .AND. Y1 .LT. PAR(1) .OR.
     >        Y1P .LT. PAR(1) - TOL) YD2 = PAR(1) - Y2
C ............................. For values above upper limit
         IF ( Y1P .GT. PAR(2) - TOL .AND. Y1 .GT. PAR(2) .OR.
     >        Y1P .GT. PAR(2) + TOL) YD2 = PAR(2) - Y2
C ............................. Check for positive crossing of lower limit
         IF (Y1 .GT. PAR(1) + TOL .AND. Y1P .LT. PAR(1) - TOL) LIM = 1
C ............................. Check for positive crossing of upper limit
         IF (Y1 .GT. PAR(2) + TOL .AND. Y1P .LT. PAR(2) - TOL) LIM = 2
C ............................. Check for negative crossing of upper limit
         IF (Y1 .LT. PAR(2) - TOL .AND. Y1P .GT. PAR(2) + TOL) LIM = -2
C ............................. Check for negative crossing of lower limit
         IF (Y1 .LT. PAR(1) - TOL .AND. Y1P .GT. PAR(1) + TOL) LIM = -1
C ---> 02.02.10/TI ............ Find time point of eventual discontinuity
         IF ( LIM .NE. 0 ) TAU = (PAR(IABS(LIM))-Y1P)/(Y1-Y1P)
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD23(IOP, MODNO, LIM, TAU, TOL, Y1P, Y1, Y2, YD2,
     >                 MST1, MST2, MM2, PAR, ZCTRL0)

C ............................. Dead zone
C ............................. 0 = y1 - ylow - y2  for y1 > ylow
C ............................. 0 = - y2  for ylow < y1 < yup
C ............................. 0 = y1 - yup - y2   for y1 < yup

C Parameters:
C     LIM - limit cross indicator
C           = -2 : decreasing crossing through upper limit
C           = -1 : decreasing crossing through lower limit
C           =  0 : no crossing
C           = +1 : increasing crossing through lower limit
C           = +2 : increasing crossing through upper limit
C     TOL - limit cross tolerance
C     Y1p - previous input value
C     PAR(1) = Ylow - lower limit
C     PAR(2) = Yup  - upper limit

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(2)
      LOGICAL ZCTRL0

      LIM = 0
C ............................. If the mode is set to initialization ...
      IF (IOP .EQ. 1) THEN
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number,
C ---> 02.01.02/TI ............ initial output value and initialise parameters
         IF (MST1 .EQ. 0) RETURN
         MST2 = 3
         MM2  = MODNO
         if (ZCTRL0) return
         IF (Y1 .LT. PAR(1)) THEN
            Y2 = Y1 - PAR(1)
         ELSE IF (Y1 .GT. PAR(2)) THEN
            Y2 = Y1 - PAR(2)
         ELSE
            Y2 = 0.0d0
         END IF
C ............................. Execution mode
      ELSE
C ............................. Compute functional value ...
C ............................. and assign it to the proper output
C ............................. In the area of no effect ...
         YD2 = - Y2
C ............................. For values below lower limit
         IF ( Y1P .LT. PAR(1) + TOL .AND. Y1 .LT. PAR(1) .OR.
     >        Y1P .LT. PAR(1) - TOL) YD2 = Y1 - PAR(1) - Y2
C ............................. For values above upper limit
         IF ( Y1P .GT. PAR(2) - TOL .AND. Y1 .GT. PAR(2) .OR.
     >        Y1P .GT. PAR(2) + TOL) YD2 = Y1 - PAR(2) - Y2
C ............................. Check for positive crossing of lower limit
         IF (Y1 .GT. PAR(1) + TOL .AND. Y1P .LT. PAR(1) - TOL) LIM = 1
C ............................. Check for positive crossing of upper limit
         IF (Y1 .GT. PAR(2) + TOL .AND. Y1P .LT. PAR(2) - TOL) LIM = 2
C ............................. Check for negative crossing of upper limit
         IF (Y1 .LT. PAR(2) - TOL .AND. Y1P .GT. PAR(2) + TOL) LIM = -2
C ............................. Check for negative crossing of lower limit
         IF (Y1 .LT. PAR(1) - TOL .AND. Y1P .GT. PAR(1) + TOL) LIM = -1
C ---> 02.02.10/TI ............ Find time point of eventual discontinuity
         IF ( LIM .NE. 0 ) TAU = (PAR(IABS(LIM))-Y1P)/(Y1-Y1P)
      END IF

      RETURN
      END


C **********************************************************************

      subroutine mod24(iop, modno, t, tnm, tol, y1p, y2p, y1, y2, yd2,
     >                 mst1, mst2, mm2, par, ZCTRL0)

c     ............................. Hysteresis
c     ...................... 0 = a*(y1-r) - y2  for y1 > y1p & y2p < a*(y1-r)
c     ...................... 0 = a*(y1-l) - y2  for y1 < y1p & y2p > a*(y1-l)
c     ...................... 0 = y2p - y2  otherwise

c     parameters:
c     t - simulation time
c     tnm - maximum end point for next integration step
c     tol - limit cross tolerance
c     Y1p - previous input value
c     Y2p - previous output value
c     par(1) = l    - left crossing with input axis
c     par(2) = r    - right crossing with input axis
c     par(3) = a    - hysteresis slope
c     par(4) = tpp  - pre-previous time
c     par(5) = y1pp - pre-previous input value
c     par(6) = tp   - previous time
c     par(7) = y1p  - previous input value
c     par(8) = tlst - last time
c     par(9) = ylst - last input value

      implicit none

      integer iop, modno, mst1, mst2, mm2
      double precision t, tnm, tol, y1p, y2p, y1, y2, yd2, par(9)
      logical ZCTRL0

      double precision zero, eps, half
      parameter ( zero = 0.0D0, eps = 1.0D-16, half = 0.5D0 )

      double precision arg, b, c, tlim, tp, tpp, y1pp, y2l, y2r, ylim

c     ............................. If the mode is set to initialization ...
      if (iop .eq. 1) then
C ............................. Unless the input variable status remains,
C ............................. set the missing variable status, module number,
C ---> 02.01.02/TI ............ initial output value and initialise history
C ---> ........................ of time and input
         IF (MST1 .EQ. 0) RETURN
         MST2 = 3
         MM2  = MODNO
         if (.not.ZCTRL0) Y2 = 0.0
         par(4) = t
         par(5) = y1
         par(6) = t
         par(7) = y1
         par(8) = t
         par(9) = y1
c     ............................. Execution mode
      else
c     ............................. Compute functional value ...
c     ............................. and assign it to the proper output
c     ............................. In the area of no effect ...
         if (t .gt. par(8)) then
            par(4) = par(6)
            par(5) = par(7)
            par(6) = par(8)
            par(7) = par(9)
         end if
         if (t .lt. par(8)) then
            tpp  = par(8)
            y1pp = par(9)
         else
            tpp  = par(4)
            y1pp = par(5)
         end if
         tp = par(6)
         b = (y1 - y1p)/(t - tp)
         tlim = t+t - tp
         if (dabs(tp-tpp) .gt. eps) then
            c = (b - (y1p - y1pp)/(tp-tpp))/(t-tpp)
            b = b - c*(t-tp)
            if (dabs(c) .gt. eps) tlim = tp - half*b/c
         else
            c = zero
         end if
         yd2 = y2 - y2p
         y2l = par(3)*(y1 - par(1))
         y2r = par(3)*(y1 - par(2))
c     ............................. If the direction is leftwards,
         if ( y1 .lt. y1p - tol .or.
     >        y1 .lt. y1p .and. y1pp .lt. y1p-tol) then
c     ............................. and horizontally,
            if (y2p .lt. y2l) then
c     ............................. check passage of upper right corner
               if ( dabs(y2p - par(3)*(y1p - par(2))) .lt. tol .and.
     >              tlim .gt. tp .and. tlim .lt. t) then
                  yd2 = y2 - par(3)*(y1 - par(2))
                  if (tlim .lt. tnm-tol) tnm = tlim
               end if
c     ............................. or if the direction is left and downwards,
            else
c     ............................. check passage of upper left corner
               if (y2p .lt. par(3)*(y1p - par(1) - tol)) then
                  ylim = y2p/par(3) + par(1)
                  if (dabs(c) .gt. eps) then
                     arg = b*b - 4*c*(y1p - ylim)
                     if (arg .ge. zero) tlim = tp - half*(sqrt(arg)+b)/c
                  else
                     if (dabs(b) .gt. eps) tlim = tp + (ylim - y1p)/b
                  end if
                  if (tlim .lt. t-tol) then
                     tnm = tlim
                  else
                     yd2 = y2 - par(3)*(y1 - par(1))
                  end if
c     ............................. else, check passage of lower left corner
               else
                  yd2 = y2 - par(3)*(y1 - par(1))
                  if (tlim .gt. tp .and. tlim .lt. t) then
                     if (tlim .lt. tnm-tol) tnm = tlim
                  end if
               end if
            end if
c     ............................. If the direction is rightwards,
         else if ( y1 .gt. y1p + tol .or.
     >             y1 .gt. y1p .and. y1pp .gt. y1p+tol ) then
c     ............................. and horizontally,
            if (y2p .gt. y2r) then
c     ............................. check passage of lower left corner
               if ( dabs(y2p - par(3)*(y1p - par(1))) .lt. tol  .and.
     >              tlim .gt. tp .and. tlim .lt. t ) then
                  yd2 = y2 - par(3)*(y1 - par(1))
                  if (tlim .lt. tnm-tol) tnm = tlim
               end if
c     ............................. or if the direction is right and upwards,
            else
c     ............................. check passage of lower right corner
               if (y2p .gt. par(3)*(y1p - par(2) + tol)) then
                  ylim = y2p/par(3) + par(2)
                  if (dabs(c) .gt. eps) then
                     arg = b*b - 4*c*(y1p - ylim)
                     if (arg .ge. zero) tlim = tp + half*(sqrt(arg)-b)/c
                  else
                     if (dabs(b) .gt. eps) tlim = tp + (ylim - y1p)/b
                  end if
                  if (tlim .lt. t-tol) then
                     tnm = tlim
                  else
                     yd2 = y2 - par(3)*(y1 - par(2))
                  end if
c     ............................. else, check passage of upper right corner
               else
                  yd2 = y2 - par(3)*(y1 - par(2))
                  if (tlim .gt. tp .and. tlim .lt. t) then
                     if (tlim .lt. tnm-tol) tnm = tlim
c     ............................. else, the direction is strictly upwards
                  end if
               end if
            end if
         end if
c     ............................. Save present time and input
         par(8) = t
         par(9) = y1
      end if

      return
      end


C **********************************************************************
C ******        GROUP 3 : COMPENSATOR ELEMENTS                          ******
C **********************************************************************

      SUBROUTINE MOD31(IOP, MODNO, Y1, Y2, Y3, YD2, YD3,
     >                 MST2, MST3, MM2, MM3, PAR, ZCTRL0)

C ............................. PI-controller element computing
C ............................. y2'= y1
C .............................  0 = Kp*(y1 + y2/Ti) - y3

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(2)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuess, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 3
         MM2  = MODNO
         MM3  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = PAR(1)*(Y1 + Y2/PAR(2)) - Y3
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD32(IOP, MODNO, Y1, Y2, Y3, Y4, YD2, YD3, YD4,
     >                 MST2, MST3, MST4, MM2, MM3, MM4, PAR, ZCTRL0)

C ............................. Proportional + limited Integral controller
C ............................. y2'= y1
C ............................. y3'= y4
C ............................. y4 = Kp*(y1 + y2/T2) - y3/T1

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(3)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = 0.0d0
         Y4 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = Y4
         YD4 = PAR(1)*(Y1 + Y2/PAR(2)) - Y3/PAR(3) - Y4
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD33(IOP, MODNO, Y1, Y2, Y3, Y4, YD2, YD3, YD4,
     >                 MST2, MST3, MST4, MM2, MM3, MM4, PAR, ZCTRL0)

C ............................. PD - controller
C ............................. y2'= y1
C ............................. y3'= y4
C .............................  0 = Kp*(Td*y1 + y2) - y3

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(2)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = 0.0d0
         Y4 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = Y4
         YD4 = PAR(1)*(PAR(2)*Y1 + Y2) - Y3
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD34(IOP, MODNO, Y1, Y2, Y3, Y4, YD2, YD3, YD4,
     >                 MST2, MST3, MST4,MM2, MM3, MM4, PAR, ZCTRL0)

C ............................. Proportional + limited Derivative controller
C ............................. y2'= y1
C ............................. y3'= y4
C ............................. y4 = (Kp*(T2*y1 + y2) - y3)/T1

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(3)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         if (ZCTRL0) return
         Y2 = (PAR(2) - PAR(3))*Y1
         Y3 = 0.0d0
         Y4 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = Y4
         YD4 = (PAR(1)*(PAR(3)*Y1 + Y2) - Y3)/PAR(2) - Y4
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD35(IOP, MODNO, Y1, Y2, Y3, Y4, Y5,
     >                 YD2, YD3, YD4, YD5, MST2, MST3, MST4, MST5,
     >                 MM2, MM3, MM4, MM5, PAR, ZCTRL0)

C ............................. PID - controller
C ............................. y2'= y1
C ............................. y3'= y2
C ............................. y4'= y5
C .............................  0 = Kp*(Td*y1 + y2 + y3/Ti) - y4

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(3)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 2
         MST5 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         MM5  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = 0.0d0
         Y4 = 0.0d0
         Y5 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = Y2
         YD4 = Y5
         YD5 = PAR(1)*(PAR(3)*Y1 + Y2 + Y3/PAR(2)) - Y4
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD36(IOP, MODNO, Y1, Y2, Y3, Y4, Y5,
     >                 YD2, YD3, YD4, YD5, MST2, MST3, MST4, MST5,
     >                 MM2, MM3, MM4, MM5, PAR, ZCTRL0)

C ............................. Proportional + Integral +
C ............................. limited Derivative controller
C ............................. y2'= y1
C ............................. y3'= Kp*(y1 + y2/Ti)
C ............................. y4'= y5
C .............................  0 = (Td*y3' + y3 - y4)/Tu - y5

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(4)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 2
         MST5 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         MM5  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = 0.0d0
         Y4 = 0.0d0
         Y5 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = PAR(1)*(Y1 + Y2/PAR(2))
         YD4 = Y5
         YD5 = (PAR(3)*YD3 + Y3 - Y4)/PAR(4) - Y5
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD37(IOP, MODNO, Y1, Y2, Y3, Y4, Y5,
     >                 YD2, YD3, YD4, YD5, MST2, MST3, MST4, MST5,
     >                 MM2, MM3, MM4, MM5, PAR, ZCTRL0)

C ............................. Proportional + limited Integral +
C ............................. limited Derivative controller
C ............................. y2'= y1
C ............................. y3'= Kp*(y1 + y2/Ti) - y3/Tl
C ............................. y4'= y5
C .............................  0 = (Td*y3' + y3 - y4)/Tu - y5

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(5)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 2
         MST5 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         MM5  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = 0.0d0
         Y4 = 0.0d0
         Y5 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = PAR(1)*(Y1 + Y2/PAR(2)) - Y3/PAR(4)
         YD4 = Y5
         YD5 = (PAR(3)*YD3 + Y3 - Y4)/PAR(5) - Y5
      END IF

      RETURN
      END


C **********************************************************************
C ******        GROUP 4 : GENERAL TRANSFER FUNCTIONS                  ******
C **********************************************************************

      SUBROUTINE MOD41(IOP, MODNO, Y1, Y2, YD2, MST2, MM2, PAR, ZCTRL0)

C ............................. Real pole (time constant)
C ............................. y2'= (k*y1 - y2)/T

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(2)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the output variable status, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MM2  = MODNO
         if (ZCTRL0) return
         Y2 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = (PAR(1)*Y1 - Y2)/PAR(2)
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD42(IOP, MODNO, Y1, Y2, Y3, YD2, YD3,
     >                 MST2, MST3, MM2, MM3, PAR, ZCTRL0)

C ............................. Complex conjugate poles (oscillator)
C ............................. y2'= k*y1 - 2*e*w*y2 - w*w*y3
C ............................. y3'= y2
C ............................. k - amplification
C ............................. w - undamped angular resonnance frequency
C ............................. e - relative damping coefficient

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(3)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MM2  = MODNO
         MM3  = MODNO
         if (ZCTRL0) return
         Y2 = 0.0d0
         Y3 = PAR(1)*Y1/(PAR(2)*PAR(2))
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = PAR(1)*Y1 - 2*PAR(2)*PAR(3)*Y2 - PAR(2)*PAR(2)*Y3
         YD3 = Y2
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD43(IOP, MODNO, Y1, Y2, Y3, Y4, YD2, YD3, YD4,
     >                 MST2, MST3, MST4, MM2, MM3, MM4, PAR, ZCTRL0)

C ............................. First order element
C ............................. y2'= y1
C ............................. y3'= y4
C .............................  0 = (k*(T2*y1 + y2) - y3)/T1 - y4

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(3)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         if (ZCTRL0) return
         Y2 = (PAR(2) - PAR(3))*Y1
         Y3 = 0.0d0
         Y4 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. Compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = Y4
         YD4 = (PAR(1)*(PAR(3)*Y1 + Y2) - Y3)/PAR(2) - Y4
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD44(IOP, MODNO, Y1, Y2, Y3, Y4, Y5, Y6,
     >                 YD2, YD3, YD4, YD5, YD6, MST2, MST3, MST4,
     >                 MST5, MST6, MM2, MM3, MM4, MM5, MM6, PAR, ZCTRL0)

C ............................. Second order element
C ............................. y2'= y1
C ............................. y3'= y2
C ............................. y4'= y6
C ............................. y5'= y4
C .............................  0 = (k*(T4*y1 + T3*y2 + y3) - T1*y4 - y5)/T1
C .............................      -y6

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PAR(5)
      LOGICAL ZCTRL0

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the state variable statuses, module number
C ---> 02.01.02/TI ............ and initial state values
         MST2 = 2
         MST3 = 2
         MST4 = 2
         MST5 = 2
         MST6 = 3
         MM2  = MODNO
         MM3  = MODNO
         MM4  = MODNO
         MM5  = MODNO
         MM6  = MODNO
         if (ZCTRL0) return
         Y2 = (PAR(2) - PAR(4))*Y1
         Y3 = (PAR(3) - PAR(5) - PAR(4)*(PAR(2) - PAR(4)))*Y1
         Y4 = 0.0d0
         Y5 = 0.0d0
         Y6 = PAR(1)*Y1
C ............................. Execution mode
      ELSE
C ............................. compute the functional values ...
C ............................. and assign them to the proper output(s)
         YD2 = Y1
         YD3 = Y2
         YD4 = Y6
         YD5 = Y4
         YD6 = (PAR(1)*(PAR(5)*Y1 + PAR(4)*Y2 + Y3) - PAR(2)*Y4
     >        -Y5)/PAR(3) - Y6
      END IF

      RETURN
      END


C **********************************************************************

      SUBROUTINE MOD45(IOP, MODNO, VREG, VDREG, MSTAT, MMOD, MMTOP,
     >                 MF, NU, NX, NY, PAR, ITYPE)

C ............................. Linear control system
C ............................. x' = Ax + Bu    y = Cx + Du

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION VREG(*), VDREG(*), PAR(*)
      INTEGER MSTAT(*), MMOD(*), MMTOP(*)

      MVU = MF - 1
      MVX = MVU + NU
      MVY = MVX + NX
      MPB = NX * NX
      MPC = MPB + NX * NU
      MPD = MPC + NY * NX

C ............................. Initialization mode
      IF (IOP .EQ. 1) THEN
C ............................. Set the missing variable status and insert
C ............................. the module number for ...
C ............................. state variables, ...
         M = MF + NU
         DO 10 N = 1,NX
            MSTAT(MMTOP(M)) = 2
            MMOD(MMTOP(M)) = MODNO
            M = M + 1
   10    CONTINUE
C ............................. and outputs (algebraic variables)
         DO 12 N = 1,NY
            MSTAT(MMTOP(M)) = 3
            MMOD(MMTOP(M)) = MODNO
            M = M + 1
   12    CONTINUE

C ............................. If the mode is set to execution ...
C ............................. and the matrices are full ...
      ELSE IF (ITYPE .EQ. 0) THEN
C ............................. compute the functional values for ...
C ............................. the state variables
         DO 24 N = 1,NX
            XD = 0.0D0
            MP = (N-1)*NX
            DO 20 NN = 1,NX
               XD = XD + PAR(MP+NN) * VREG(MMTOP(MVX+NN))
   20       CONTINUE
            MP = MPB + (N-1)*NU
            DO 22 NN = 1,NU
               XD = XD + PAR(MP+NN) * VREG(MMTOP(MVU+NN))
   22       CONTINUE
            VDREG(MMTOP(MVX+N)) = XD
   24    CONTINUE
C ............................. and the outputs (algebraic variables)

         DO 34 N = 1,NY
            YD = 0.0D0
            MP = MPC + (N-1)*NX
            DO 30 NN = 1,NX
               YD = YD + PAR(MP+NN) * VREG(MMTOP(MVX+NN))
   30       CONTINUE
            MP = MPD + (N-1)*NU
            DO 32 NN = 1,NU
               YD = YD + PAR(MP+NN) * VREG(MMTOP(MVU+NN))
   32       CONTINUE
            VDREG(MMTOP(MVY+N)) = YD - VREG(MMTOP(MVY+N))
   34    CONTINUE

C ............................. If the mode is set to execution ...
C ............................. and the system transfer function is:
C ............................. y(s)               b
C ............................. ---- = --------------------------
C ............................. u(s)   a(0) + a(1)s + ....+ s**NX
      ELSE IF (ITYPE .EQ. 1) THEN
C ............................. compute the functional values for the state
C ............................. variables and the output
         NXM = NX - 1
         XD = 0.0D0
         DO 40 N = 1,NX
            XD = XD - PAR(N) * VREG(MMTOP(MVX+N))
   40    CONTINUE
         XD = XD + PAR(NX+1) * VREG(MMTOP(MF))
         VDREG(MMTOP(MVX+1)) = XD
         DO 42 N = 1,NXM
            MV = MVX + N
            VDREG(MMTOP(MV+1)) = VREG(MMTOP(MV))
   42    CONTINUE
         VDREG(MMTOP(MVY+1)) = VREG(MMTOP(MVY)) - VREG(MMTOP(MVY+1))
      END IF

      RETURN
      END


c **********************************************************************

      subroutine mod46(iop, modno, vreg, vdreg, mstat, mmod, mmtop,
     >                 mf, nx, par)

c ............................. General transfer function
c .............................
c ............................. y(s)   b(n)s**n + ... + b(1)s + b(0)
c ............................. ---- = -----------------------------
c ............................. u(s)   a(n)s**n + ... + a(1)s + a(0)
c .............................
c ............................. where a(n) = 1

      implicit double precision (a-h, o-z)
      double precision vreg(*), vdreg(*), par(*)
      integer mstat(*), mmod(*), mmtop(*)

      mfp = mf

c ............................. Initialization mode
      if (iop .eq. 1) then
c ............................. Set the missing variable status and insert
c ............................. the module number for ...
c ............................. state variables, ...
         m = mf + 1
         do 10 n = 1,nx
            mstat(mmtop(m)) = 2
            mmod(mmtop(m)) = modno
            m = m + 1
   10    continue
c ............................. and outputs (algebraic variables)
         mstat(mmtop(m)) = 3
         mmod(mmtop(m)) = modno

c ............................. Execution mode
      else

c ............................. Compute the functional values for the state
c ............................. variables and the output
         nxm = nx - 1
         xd = vreg(mmtop(mf))
         do 20 n = 1,nx
            xd = xd - par(n) * vreg(mmtop(mfp+n))
   20    continue
         vdreg(mmtop(mfp+1)) = xd
         do 22 n = 1,nxm
            mv = mfp + n
            vdreg(mmtop(mv+1)) = vreg(mmtop(mv))
   22    continue
         yd = par(nx+1) * vdreg(mmtop(mfp+1))
         do 24 n = 1,nx
            yd = yd + par(nx+1+n) * vreg(mmtop(mfp+n))
   24    continue
         vdreg(mmtop(mfp+nx+1)) = yd - vreg(mmtop(mfp+nx+1))
      end if

      return
      end
