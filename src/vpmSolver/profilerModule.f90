!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file profilerModule.f90
!>
!> @brief Subroutines for profiling of various program parts.

!!==============================================================================
!> @brief Module with subroutines for profiling of various program parts.
!>
!> @details This module contains parameters used to identify the various program
!> parts to be profiled, and a subroutine for printing the profiling results.
!> Two subroutines for starting and stopping the timing of a task is imported
!> from the timermodule.

module ProfilerModule

  use TimerModule, only : startTimer, stopTimer

  implicit none

  integer, parameter :: tot_p  =  1 !< Total time
  integer, parameter :: ini_p  =  2 !< Data input and initializations
  integer, parameter :: asm_p  =  3 !< Finite element assembly
  integer, parameter :: sol_p  =  4 !< Linear equation solver
  integer, parameter :: eig_p  =  5 !< Eigenvalue solver
  integer, parameter :: upd_p  =  6 !< Configuration update
  integer, parameter :: sav_p  =  7 !< Results saving
  integer, parameter :: exp_p  =  8 !< Curve export
  integer, parameter :: ctrl_p =  9 !< Control systems
  integer, parameter :: tir_p  = 10 !< Tire models
  integer, parameter :: hyd_p  = 11 !< Hydrodynamic models
  integer, parameter :: wav_p  = 12 !< Wave kinematics evaluation
  integer, parameter :: aed_p  = 13 !< Aerodynamic models
  integer, parameter :: sup1_p = 14 !< Superelement matrix assembly
  integer, parameter :: sup2_p = 15 !< Superelement vector assembly
  integer, parameter :: sup3_p = 16 !< csAddEM (element matrix assembly)
  integer, parameter :: rec_p  = 17 !< Stress & Strain Gage recovery
  integer, parameter :: rec1_p = 18 !< Initialization of recovery data
  integer, parameter :: rec2_p = 19 !< Deformation recovery
  integer, parameter :: fra_p  = 20 !< Frequency response analysis
  integer, parameter :: oth_p  = 21 !< Other, i.e., total time minus all others

  !> @brief Total number of parts to be profiled
  integer, parameter :: nProfMod_p = oth_p-1

  private :: tot_p, oth_p


contains

  !!============================================================================
  !> @brief Prints out solver profiling information in a nicely formatted table.
  !>
  !> @param[in] lpu File unit number for res-file output
  !> @param[in] nSol Number of iterations (with one linear equation solve each)
  !> @param[in] nEig Number of eigenvalue solutions
  !> @param[in] nStep Number of time/load increments
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Jan 2003

  subroutine reportTiming (lpu,nSol,nEig,nStep)

    use KindModule , only : sp, i8
    use TimerModule, only : getTotalTime

    integer    , intent(in) :: lpu, nEig
    integer(i8), intent(in) :: nSol, nStep

    !! Local variables
    integer  :: i
    logical  :: lCtrl, lTire, lWave, lHdyn, lAdyn, lCexp, lSrec, lFrRA
    real(sp) :: dSol, dStep, totalTime(2,nProfMod_p+1)

    !! --- Logic section ---

    do i = 1, nProfMod_p
       totalTime(:,i) = getTotalTime(i)
    end do
    do i = 1, 2
       totalTime(i,oth_p) = totalTime(i,tot_p) - sum(totalTime(i,2:ctrl_p)) &
            &             - totalTime(i,aed_p) - sum(totalTime(i,rec_p:rec1_p))
    end do
    dSol  = real(max(1_i8,nSol),sp)
    dStep = real(max(1_i8,nStep),sp)
    lCtrl = totalTime(2,ctrl_p) >= 0.05 .and. nSol > 0_i8
    lTire = totalTime(2,tir_p)  >= 0.05 .and. nSol > 0_i8
    lWave = totalTime(2,wav_p)  >= 0.05 .and. nStep > 0_i8
    lHdyn = totalTime(2,hyd_p)  >= 0.05 .and. nStep > 0_i8
    lAdyn = totalTime(2,aed_p)  >= 0.05 .and. nStep > 0_i8
    lSrec = totalTime(2,rec_p)  >= 0.05 .and. nStep > 0_i8
    lFrRA = totalTime(2,fra_p)  >= 0.05
    lCexp = totalTime(2,exp_p)  >= 0.05

    !! The tire model times are also included in the configuration update times
    if (lTire) totalTime(:,upd_p) = totalTime(:,upd_p) - totalTime(:,tir_p)
    if (lHdyn) totalTime(:,upd_p) = totalTime(:,upd_p) - totalTime(:,hyd_p)
    !! The wave kinematics times are also included in the hydrodynamics times
    if (lWave) totalTime(:,hyd_p) = totalTime(:,hyd_p) - totalTime(:,wav_p)

    write(lpu,600) totalTime(:,ini_p), &
         &         totalTime(:,asm_p), totalTime(:,asm_p)/dSol, &
         &         totalTime(:,sup1_p),totalTime(:,sup1_p)/dSol, &
         &         totalTime(:,sup3_p),totalTime(:,sup3_p)/dSol, &
         &         totalTime(:,sup2_p),totalTime(:,sup2_p)/dSol, &
         &         totalTime(:,sol_p), totalTime(:,sol_p)/dSol
    if (nEig > 0)  write(lpu,610) totalTime(:,eig_p),  totalTime(:,eig_p)/nEig
    if (lCtrl)     write(lpu,620) totalTime(:,ctrl_p), totalTime(:,ctrl_p)/dSol
    if (lTire)     write(lpu,624) totalTime(:,tir_p),  totalTime(:,tir_p)/dSol
    if (lWave)     write(lpu,625) totalTime(:,wav_p),  totalTime(:,wav_p)/dStep
    if (lHdyn)     write(lpu,626) totalTime(:,hyd_p),  totalTime(:,hyd_p)/dStep
    if (lAdyn)     write(lpu,627) totalTime(:,aed_p),  totalTime(:,aed_p)/dStep
    write(lpu,630) totalTime(:,upd_p), totalTime(:,upd_p)/dSol, &
         &         totalTime(:,sav_p), totalTime(:,sav_p)/dStep
    if (lSrec) then
       if (totalTime(2,rec1_p) >= 0.05) write(lpu,631) totalTime(:,rec1_p)
       write(lpu,632) totalTime(:,rec_p),  totalTime(:,rec_p)/dStep, &
            &         totalTime(:,rec2_p), totalTime(:,rec2_p)/dStep
    end if
    if (lFrRA)     write(lpu,633) totalTime(:,fra_p)
    if (lCexp)     write(lpu,635) totalTime(:,exp_p)
    write(lpu,640) totalTime(:,oth_p), totalTime(:,oth_p)/dStep, &
         &         totalTime(:,tot_p), totalTime(:,tot_p)/dStep

600 format(// 41X,'>>>>>>   P R O F I L I N G   D A T A   <<<<<<' &
         &  / 40X,47('=') &
         & // 40X,'|     Total (sec)      |  Per iteration/step  |' &
         &  / 40X,'| Wall time  CPU time  | Wall time  CPU time  |' &
         &  /  4X,'+',35('-'),'+',2(22('-'),'+') &
         &  /  4X,'| Data input and initialisations ...|',2F10.1,'  |', &
         &                                                20X   ,'  |'  &
         &  /  4X,'| Finite element assembly ..........|',2F10.1,'  |', &
         &                                                2F10.3,'  |'  &
         &  /  4X,'|   Superelement matrices ..........|',2F10.1,'  |', &
         &                                                2F10.3,'  |'  &
         &  /  4X,'|     csAddEM ......................|',2F10.1,'  |', &
         &                                                2F10.3,'  |'  &
         &  /  4X,'|   Superelement force vectors .....|',2F10.1,'  |', &
         &                                                2F10.3,'  |'  &
         &  /  4X,'| Linear equation solution .........|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
610 format(    4X,'| Eigenvalue analysis ..............|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
620 format(    4X,'| Control systems ..................|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
624 format(    4X,'| Tire models ......................|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
625 format(    4X,'| Wave kinematics evaluations ......|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
626 format(    4X,'| Hydrodynamic calculations ........|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
627 format(    4X,'| Aerodynamic calculations (AeroDyn)|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
630 format(    4X,'| Configuration update .............|',2F10.1,'  |', &
         &                                                2F10.3,'  |'  &
         &  /  4X,'| Results saving ...................|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
631 format(    4X,'| Recovery matrix initialization ...|',2F10.1,'  |', &
         &                                                20X   ,'  |' )
632 format(    4X,'| Stress/strain recovery ...........|',2F10.1,'  |', &
         &                                                2F10.3,'  |', &
         &  /  4X,'|   Displacement recovery ..........|',2F10.1,'  |', &
         &                                                2F10.3,'  |' )
633 format(    4X,'| Frequency response analysis ......|',2F10.1,'  |', &
         &                                                20X   ,'  |' )
635 format(    4X,'| Curve export .....................|',2F10.1,'  |', &
         &                                                20X   ,'  |' )
640 format(    4X,'| Other ............................|',2F10.1,'  |', &
         &                                                2F10.3,'  |'  &
         &  /  4X,'+',35('-'),'+',2(22('-'),'+') &
         &  /  4X,'| Total ............................|',2F10.1,'  |', &
         &                                                2F10.3,'  |'  &
         &  /  4X,'+',35('-'),'+',2(22('-'),'+') / )

  end subroutine reportTiming

end module ProfilerModule
