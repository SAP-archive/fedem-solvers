!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

#ifdef FT_HAS_MKL
include "mkl_dfti.f90"
#endif

module FFTmodule

  !!============================================================================
  !! This module contains three subroutines for performing Fast Fourier
  !! transformations, using the Math Kernel Library from Intel.
  !!============================================================================

#ifdef FT_HAS_MKL
  use MKL_DFTI, only : DFTI_DESCRIPTOR
#endif

  implicit none

#ifdef FT_HAS_MKL
  type(DFTI_DESCRIPTOR), save, private, pointer :: hand => null()
#endif


contains

  subroutine initFFT (nSamp,ierr)

#ifdef FT_HAS_MKL
    use MKL_DFTI, only : DFTI_DOUBLE_R, DFTI_REAL
    use MKL_DFTI, only : DFTI_PLACEMENT, DFTI_NOT_INPLACE, DFTI_BACKWARD_SCALE
    use MKL_DFTI, only : DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX
    use MKL_DFTI, only : DftiCreateDescriptor, DftiFreeDescriptor
    use MKL_DFTI, only : DftiCommitDescriptor, DftiSetValue
#endif

    use reportErrorModule, only : reportError, error_p

    integer, intent(in)  :: nSamp
    integer, intent(out) :: ierr

    !! --- Logic section ---

#ifdef FT_HAS_MKL
    if (associated(hand)) then
       ierr = DftiFreeDescriptor(hand)
       if (ierr /= 0) goto 990
       nullify(hand)
    end if

    ierr = DftiCreateDescriptor(hand,DFTI_DOUBLE_R,DFTI_REAL,1,nSamp)
    if (ierr /= 0) goto 990

    ierr = DftiSetValue(hand,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
    if (ierr /= 0) goto 990

    ierr = DftiSetValue(hand,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)
    if (ierr /= 0) goto 990

    ierr = DftiSetValue(hand,DFTI_BACKWARD_SCALE,dble(1)/dble(nSamp))
    if (ierr /= 0) goto 990

    ierr = DftiCommitDescriptor(hand)
    if (ierr == 0) return
#else
    ierr = -nSamp
    goto 990
#endif

990 call reportError (error_p,'Failed to initialize FFT module')

  end subroutine initFFT


  subroutine closeFFT (ierr)

#ifdef FT_HAS_MKL
    use MKL_DFTI, only : DftiFreeDescriptor
#endif

    use reportErrorModule, only : reportError, error_p

    integer, intent(out) :: ierr

    !! --- Logic section ---

#ifdef FT_HAS_MKL
    if (associated(hand)) then
       ierr = DftiFreeDescriptor(hand)
       nullify(hand)
    else
       ierr = 0
    end if
#else
    ierr = -999
#endif
    if (ierr /= 0) call reportError (error_p,'Failed to close FFT module')

  end subroutine closeFFT


  subroutine doFFT (invers,dataR,dataC,ierr)

#ifdef FT_HAS_MKL
    use MKL_DFTI, only : DftiComputeForward, DftiComputeBackward
#endif

    use kindModule       , only : dp
    use reportErrorModule, only : reportError, error_p

    logical    , intent(in)    :: invers
    real(dp)   , intent(inout) :: dataR(:)
    complex(dp), intent(inout) :: dataC(:)
    integer    , intent(out)   :: ierr

    !! --- Logic section ---

#ifdef FT_HAS_MKL
    if (.not. associated(hand)) then
       ierr = 99
       call reportError (error_p,'FFT module is not initialized')
    else if (invers) then
       ierr = DftiComputeBackward(hand,dataC,dataR)
    else
       ierr = DftiComputeForward(hand,dataR,dataC)
    end if
#else
    if (invers) then
       ierr = -size(dataC)
    else
       ierr = -size(dataR)
    end if
#endif
    if (ierr /= 0) call reportError (error_p,'FFT calculation failed')

  end subroutine doFFT

end module FFTmodule


module FNVmodule

  !!============================================================================
  !! This module contains the core of the higher-order wave force computations
  !! on vertical cylinders, based on the FNV-method (Faltinsen, Newman, Vinje).
  !!
  !! Most of this material is a reimplementation of Matlab routines provided by
  !! Jørgen R. Krogstad, Statkraft.
  !!============================================================================

  use KindModule, only : dp, pi_p

  implicit none

  type MatrixPtr
     real(dp), pointer :: p(:,:)
  end type MatrixPtr

  private :: filter_time


contains

  subroutine filter_time (data,samp,LPcut,LPorder,LP,ierr)

    !!==========================================================================
    !! Performes Butterworth LP - filtering in the frequency domain.
    !! A transfer function of the filter is calculated.
    !! Then multiply and inverse Fourier transformation.
    !!
    !! Programmer : Knut Morten Okstad (based on filter_time.m from J.Krogstad)
    !! date/rev   : 21 Oct 2011/1.0
    !!==========================================================================

    use FFTmodule        , only : doFFT
    use reportErrorModule, only : allocationError, reportError, debugFileOnly_p

    real(dp), intent(inout) :: data(:)
    real(dp), intent(in)    :: samp, LPcut
    integer , intent(in)    :: LPorder
    logical , intent(in)    :: LP
    integer , intent(out)   :: ierr

    !! Local variables
    integer  :: i, nSamp, nFreq
    real(dp) :: delta, wfreq, ffunc

    complex(dp), allocatable :: respi(:), F(:)

    !! --- Logic section ---

    nSamp = size(data)
    nFreq = size(data)/2 + 1
    delta = (pi_p+pi_p) / samp

    allocate(F(nSamp),respi(nSamp),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('filter_time')
       return
    end if

    !! Fourier transform
    call doFFT (.false.,data,F,ierr)
    if (ierr /= 0) goto 990

    !! Calculate the Butterworth n-th order filter transfer function
    do i = 1, nFreq
       wfreq = real(i-1,dp)*delta
       ffunc = 1.0_dp + (wfreq/LPcut)**real(2*LPorder,dp)
       !! Mirror the filter function if in high-pass mode
       if (LP) then
          respi(i) = F(i) / ffunc
       else
          respi(i) = F(i) * (1.0_dp-1.0_dp/ffunc)
       end if
    end do

    !! Adjust mean values
    respi(1) = real(respi(1),dp)
    respi(nFreq) = real(respi(nFreq),dp)

    !! Mirror negative frequencies
    do i = 1, nSamp-nFreq
       respi(nFreq+i) = conjg(respi((nSamp+1)/2+1-i))
    end do

    !! Inverse fourier transform the signal
    call doFFT (.true.,data,respi,ierr)
    if (ierr /= 0) goto 990

    deallocate(F,respi)
    return

990 call reportError (debugFileOnly_p,'filter_time')
    deallocate(F,respi)

  end subroutine filter_time


  subroutine ring_kin_lin (X,Y,Z,waveAng,g,D,samp,dT,LPcut,iorder, &
       &                   Wavcomp,Wavekin,ierr)

    !!==========================================================================
    !! Calculate wave kinematics components at specified position in the
    !! frequency domain and inverse Fourier the kinematics to the time domain.
    !! Use linear theory up to the mean surface elevation.
    !!
    !! Wavcomp(1,:) = wave amplitudes
    !! Wavcomp(2,:) = angular frequencies
    !! Wavcomp(3,:) = phase angles
    !! Wavcomp(4,:) = wave numbers
    !!
    !! Wavekin(:,1) = uu     - horizontal velocity
    !! Wavekin(:,2) = dudt   - horizontal acceleration
    !! Wavekin(:,3) = dudtdz - horizontal acceleration gradient
    !! Wavekin(:,4) = ww     - vertical velocity
    !! Wavekin(:,5) = dudz   - horizontal velocity gradient
    !! Wavekin(:,6) = dudx   - horizontal velocity gradient
    !! Wavekin(:,7) = dwdt   - vertical acceleration
    !! Wavekin(:,8) = eta1   - first order wave elevation
    !! Wavekin(:,9) = eta2   - second order wave elevation
    !!
    !! Programmer : Knut Morten Okstad (based on ring_kin_lin.m from J.Krogstad)
    !! date/rev   : 21 Oct 2011/1.0
    !!==========================================================================

    use FFTmodule        , only : initFFT, closeFFT, doFFT
    use reportErrorModule, only : allocationError, reportError, debugFileOnly_p

    real(dp), intent(in)  :: X, Y, Z, waveAng, g, D, samp, dT, LPcut
    integer , intent(in)  :: iorder
    real(dp), intent(in)  :: Wavcomp(:,:)
    real(dp), pointer     :: Wavekin(:,:)
    integer , intent(out) :: ierr

    !! Local variables
    integer                  :: i, Nfreq, Nsamp
    real(dp)                 :: kx, ky, Phase, zz
    real(dp)   , allocatable :: tanD(:)
    complex(dp), allocatable :: F(:), FF(:)

    !! --- Logic section ---

    Nfreq = size(Wavcomp,2)
    Nsamp = int(samp/dT)
    if (Z > 0.0_dp) then ! evaluate wave kinematics at the surface
       allocate(F(Nfreq),FF(Nsamp),Wavekin(Nsamp,9),tanD(Nfreq),STAT=ierr)
    else ! only dudt, dudz and ww are needed (i.e., for 1st+2nd-order terms)
       allocate(F(Nfreq),FF(Nsamp),Wavekin(Nsamp,4),STAT=ierr)
    end if
    if (ierr /= 0) then
       ierr = allocationError('ring_kin_lin')
       return
    end if

    call initFFT (Nsamp,ierr)
    if (ierr /= 0) goto 990

    do i = 1, Nfreq
       kx = Wavcomp(4,i)*cos(waveAng) ! wave number directional x-component
       ky = Wavcomp(4,i)*sin(waveAng) ! wave number directional y-component
       Phase = Wavcomp(3,i) - (kx*X+ky*Y) ! sum of phase components
       F(i) = real(Nfreq,dp)*Wavcomp(1,i)*exp(cmplx(0,Phase,dp)) ! elevation
    end do

    if (Z > 0.0_dp) then

       !! -----------------------------------------
       !! --- Evaluate the 1/tanh(k*d) term
       !! -----------------------------------------

       do i = 1, Nfreq
          if (Wavcomp(4,i)*D > 0.0_dp .and. Wavcomp(4,i)*D < 500.0_dp) then
             tanD(i) = 1.0_dp/tanh(Wavcomp(4,i)*D)
          else
             tanD(i) = 0.0_dp
          end if
       end do

       !! -----------------------------------------
       !! --- Surface elevation, eta1
       !! -----------------------------------------

       FF(1:Nfreq) = F

       call iFFTandFilter (Wavekin(:,8))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Velocities uu
       !! -----------------------------------------

       FF(1:Nfreq) = Wavcomp(2,:)*F*tanD

       call iFFTandFilter (Wavekin(:,1))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Acceleration dudt
       !! -----------------------------------------

       FF(1:Nfreq) = cmplx(0,Wavcomp(2,:)**2.0_dp,dp)*F*tanD

       call iFFTandFilter (Wavekin(:,2))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Acceleration dudtdz
       !! -----------------------------------------

       FF(1:Nfreq) = cmplx(0,Wavcomp(2,:)**2.0_dp*Wavcomp(4,:),dp)*F

       call iFFTandFilter (Wavekin(:,3))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Velocities ww
       !! -----------------------------------------

       FF(1:Nfreq) = cmplx(0,Wavcomp(2,:),dp)*F

       call iFFTandFilter (Wavekin(:,4))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Kinematic dudz
       !! -----------------------------------------

       FF(1:Nfreq) = (Wavcomp(2,:)*Wavcomp(4,:))*F*tanD

       call iFFTandFilter (Wavekin(:,5))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Kinematic dudx
       !! -----------------------------------------

       FF(1:Nfreq) = cmplx(0,-Wavcomp(2,:)*Wavcomp(4,:),dp)*F*tanD

       call iFFTandFilter (Wavekin(:,6))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Kinematic dwdt
       !! -----------------------------------------

       FF(1:Nfreq) = -(Wavcomp(2,:)**2.0_dp)*F

       call iFFTandFilter (Wavekin(:,7))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- eta2 = -(1/g)*(0.5*(uu^2 + ww^2) + dwdt*eta1)
       !! -----------------------------------------

       Wavekin(:,9) = (-1.0_dp/g) * ( 0.5_dp*(Wavekin(:,1)**2.0_dp + &
            &                                 Wavekin(:,4)**2.0_dp) + &
            &                         Wavekin(:,7)*Wavekin(:,8) )

    else ! Z <= 0.0 i.e. calculation of 1st (and 2nd) order terms only
       zz = Z + D

       !! -----------------------------------------
       !! --- Velocities uu
       !! -----------------------------------------

       FF(1:Nfreq) = Wavcomp(2,:)*F
       do i = 1, Nfreq
          if (Wavcomp(4,i)*D > 0.0_dp .and. Wavcomp(4,i)*D < 500.0_dp) then
             FF(i) = FF(i) * cosh(Wavcomp(4,i)*zz)/sinh(Wavcomp(4,i)*D)
          else
             FF(i) = FF(i) * exp(Wavcomp(4,i)*Z)
          end if
       end do

       call iFFTandFilter (Wavekin(:,1))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Velocities ww
       !! -----------------------------------------

       FF(1:Nfreq) = cmplx(0,Wavcomp(2,:),dp)*F
       do i = 1, Nfreq
          if (Wavcomp(4,i)*D > 0.0_dp .and. Wavcomp(4,i)*D < 500.0_dp) then
             FF(i) = FF(i) * sinh(Wavcomp(4,i)*zz)/sinh(Wavcomp(4,i)*D)
          else
             FF(i) = FF(i) * exp(Wavcomp(4,i)*Z)
          end if
       end do

       call iFFTandFilter (Wavekin(:,2))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Acceleration dudt
       !! -----------------------------------------

       FF(1:Nfreq) = cmplx(0,Wavcomp(2,:)**2.0_dp,dp)*F
       do i = 1, Nfreq
          if (Wavcomp(4,i)*D > 0.0_dp .and. Wavcomp(4,i)*D < 500.0_dp) then
             FF(i) = FF(i) * cosh(Wavcomp(4,i)*zz)/sinh(Wavcomp(4,i)*D)
          else
             FF(i) = FF(i) * exp(Wavcomp(4,i)*Z)
          end if
       end do

       call iFFTandFilter (Wavekin(:,3))
       if (ierr /= 0) goto 990

       !! -----------------------------------------
       !! --- Kinematic dudz
       !! -----------------------------------------

       FF(1:Nfreq) = (Wavcomp(2,:)*Wavcomp(4,:))*F
       do i = 1, Nfreq
          if (Wavcomp(4,i)*D > 0.0_dp .and. Wavcomp(4,i)*D < 500.0_dp) then
             FF(i) = FF(i) * sinh(Wavcomp(4,i)*zz)/sinh(Wavcomp(4,i)*D)
          else
             FF(i) = FF(i) * exp(Wavcomp(4,i)*Z)
          end if
       end do

       call iFFTandFilter (Wavekin(:,4))
       if (ierr /= 0) goto 990

    end if

    call closeFFT (ierr)
    if (ierr /= 0) goto 990

    deallocate(F,FF)
    if (allocated(tanD)) deallocate(tanD)

    return

990 call reportError (debugFileOnly_p,'ring_kin_lin')
    deallocate(F,FF)
    if (allocated(tanD)) deallocate(tanD)

  contains

    subroutine iFFTandFilter (data)
      real(dp), intent(out) :: data(:)

      !! Adjust mean values
      FF(1)     = real(FF(1),dp)
      FF(Nfreq) = real(FF(Nfreq),dp)

      !! Mirror negative frequencies
      do i = 1, Nsamp-Nfreq
         FF(Nfreq+i) = conjg(FF((Nsamp+1)/2+1-i))
      end do

      !! In-phase component
      call doFFT (.true.,data,FF,ierr)
      if (ierr /= 0) return

      !! Apply low-pass filter
      call filter_time (data,samp,LPcut,iorder,.true.,ierr)

    end subroutine iFFTandFilter

  end subroutine ring_kin_lin


  subroutine reg_kin_lin (X,Y,Z,waveAng,g,D,samp,dT,Wavcomp,Wavekin,ierr)

    !!==========================================================================
    !! Calculate regular wave kinematics at specified position in time domain.
    !! Use linear theory up to the mean surface elevation.
    !!
    !! Wavcomp(1,1) = constant offset (zero frequency amplitude)
    !! Wavcomp(1,2) = wave amplitude
    !! Wavcomp(2,2) = angular frequency
    !! Wavcomp(3,2) = phase angle
    !! Wavcomp(4,2) = wave number
    !!
    !! Wavekin(:,1) = uu     - horizontal velocity
    !! Wavekin(:,2) = dudt   - horizontal acceleration
    !! Wavekin(:,3) = dudtdz - horizontal acceleration gradient
    !! Wavekin(:,4) = ww     - vertical velocity
    !! Wavekin(:,5) = dudz   - horizontal velocity gradient
    !! Wavekin(:,6) = dudx   - horizontal velocity gradient
    !! Wavekin(:,7) = dwdt   - vertical acceleration
    !! Wavekin(:,8) = eta1   - first order wave elevation
    !! Wavekin(:,9) = eta2   - second order wave elevation
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 24 Oct 2011/1.0
    !!==========================================================================

    use reportErrorModule, only : allocationError, internalError

    real(dp), intent(in)  :: X, Y, Z, waveAng, g, D, samp, dT
    real(dp), intent(in)  :: Wavcomp(:,:)
    real(dp), pointer     :: Wavekin(:,:)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, Nsamp
    real(dp) :: t, c, s, k, A, omega, Phase, Offset, tanD, chzD, shzD

    !! --- Logic section ---

    Nsamp = 1 + int(samp/dT)
    if (Z > 0.0_dp) then
       allocate(Wavekin(Nsamp,9),STAT=ierr)
    else ! only dudt, dudz and ww are needed (i.e., for 1st+2nd-order terms)
       allocate(Wavekin(Nsamp,4),STAT=ierr)
    end if
    if (ierr /= 0) then
       ierr = allocationError('reg_kin_lin')
       return
    end if

    i = size(Wavcomp,2)
    if (i == 1) then
       Offset = 0.0_dp
    else if (i == 2) then
       Offset = Wavcomp(1,1)
    else
       ierr = internalError('reg_kin_lin')
       return
    end if

    A     = Wavcomp(1,i)
    omega = Wavcomp(2,i)
    k     = Wavcomp(4,i)
    Phase = Wavcomp(3,i) - k*(X*cos(waveAng) + Y*sin(waveAng))

    if (Z > 0.0_dp) then
       tanD = 1.0/tanh(k*D)
       shzD = 0.0_dp
       chzD = 0.0_dp
    else
       tanD = 0.0_dp
       shzD = sinh(k*(Z+D))/sinh(k*D)
       chzD = cosh(k*(Z+D))/sinh(k*D)
    end if

    !! Generate the time series
    do i = 1, Nsamp
       t = real(i-1,dp)*dT
       s = sin(omega*t+Phase)
       c = cos(omega*t+Phase)
       if (Z > 0.0_dp) then
          Wavekin(i,1) =         omega*A*s*tanD ! u (horizontal velocity)
          Wavekin(i,2) =   omega*omega*A*c*tanD ! du/dt
          Wavekin(i,3) = k*omega*omega*A*c      ! d2u/(dtdz)
          Wavekin(i,4) =         omega*A*c      ! w (vertical velocity)
          Wavekin(i,5) =       k*omega*A*s      ! du/dz
          Wavekin(i,6) =      -k*omega*A*c*tanD ! du/dx
          Wavekin(i,7) =  -omega*omega*A*s      ! dw/dt
          Wavekin(i,8) =      Offset + A*s      ! eta1
          Wavekin(i,9) = (-1.0_dp/g) * ( 0.5_dp*(Wavekin(i,1)**2.0_dp + &
               &                                 Wavekin(i,4)**2.0_dp) + &
               &                         Wavekin(i,7)*Wavekin(i,8) ) ! eta2
       else
          Wavekin(i,1) =         omega*A*s*chzD ! u (horizontal velocity)
          Wavekin(i,2) =         omega*A*c*shzD ! w (vertical velocity)
          Wavekin(i,3) =   omega*omega*A*c*chzD ! du/dt
          Wavekin(i,4) =       k*omega*A*s*shzD ! du/dz
       end if

    end do

  end subroutine reg_kin_lin


  subroutine FNV12order (z,D,rho,Wavekin,F12,F12Tot,useFNV,ierr)

    !!==========================================================================
    !! Integrate the first- and second-order terms of the FNV wave force model.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 8 Nov 2011/1.0
    !!==========================================================================

    use reportErrorModule, only : allocationError

    real(dp)       , intent(in)  :: z(:), D, rho
    type(MatrixPtr), intent(in)  :: Wavekin(:)
    real(dp)       , pointer     :: F12(:,:), F12Tot(:,:)
    integer        , intent(in)  :: useFNV
    integer        , intent(out) :: ierr

    !! Local variables
    integer  :: i, j, nNod, nWC
    real(dp) :: C, dz, f0(2), f1(2), f2(2), fi(2)
    logical  :: Force2

    !! --- Logic section ---

    Force2 = mod(useFNV,1000)/100 > 0 ! Third last digit
    nNod = size(z)
    nWC  = size(Wavekin(1)%p,1)
    allocate(F12(nWC,nnod),F12Tot(nWC,2),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('FNV12order')
       return
    end if

    !! Integrate over z = z_1, z_1.5 for node 1
    !! Integrate over z = z_(i-0.5), z_(i+0.5) for node i, i = 2, nNod-1
    !! Integrate over z = z_(nNod-0.5), 0.0 for the last node nNod

    C = pi_p*rho * D*D * 0.25_dp

    !! uu   = Wavekin(i,1) ! horizontal velocity
    !! ww   = Wavekin(i,2) ! vertical velocity
    !! dudt = Wavekin(i,3) ! horizontal acceleration
    !! dudz = Wavekin(i,4) ! horizontal velocity gradient

    f2(2) = 0.0_dp
    do i = 1, nWC
       f1(1) = (C+C)*Wavekin(1)%p(i,3)
       f2(1) = (C+C)*Wavekin(2)%p(i,3)
       if (Force2) then
          f1(2) = C*Wavekin(1)%p(i,2)*Wavekin(1)%p(i,4)
          f2(2) = C*Wavekin(2)%p(i,2)*Wavekin(2)%p(i,4)
       end if
       dz = 0.5_dp*(z(2)-z(1))
       fi = (f1 + 0.5_dp*f2)*dz
       F12(i,1) = fi(1) + fi(2)
       F12Tot(i,:) = fi

       do j = 2, nNod-1
          f0 = f1
          f1 = f2
          f2(1) = (C+C)*Wavekin(j+1)%p(i,3)
          if (Force2) then
             f2(2) = C*Wavekin(j+1)%p(i,2)*Wavekin(j+1)%p(i,4)
          end if
          dz = 0.25_dp*(z(j+1)-z(j-1))
          fi = (f1 + 0.5_dp*(f0+f2))*dz
          F12(i,j) = fi(1) + fi(2)
          F12Tot(i,:) = F12Tot(i,:) + fi
       end do

       dz = 0.5_dp*(z(nNod)-z(nNod-1))
       fi = (f2 + 0.5_dp*f1)*dz
       F12(i,nNod) = fi(1) + fi(2)
       F12Tot(i,:) = F12Tot(i,:) + fi

    end do

  end subroutine FNV12order


  subroutine FNV23order (D,rho,g,LPcut,samp,Wavcomp,Wavekin, &
       &                 F3Tot,M4Tot,useFNV,ierr)

    !!==========================================================================
    !! Calculate time history of high frequency force, 2nd+3rd order or higher,
    !! in one specified nodal point on a vertical column.
    !!
    !! Programmer : Knut Morten Okstad (based on FNV23Order.m from J.Krogstad)
    !! date/rev   : 21 Oct 2011/1.0
    !!==========================================================================

    use FFTmodule        , only : initFFT, closeFFT, doFFT
    use reportErrorModule, only : allocationError
    use reportErrorModule, only : reportError, debugFileOnly_p, note_p, error_p

    real(dp), intent(in)  :: D, rho, g, LPcut, samp
    real(dp), intent(in)  :: Wavcomp(:,:), Wavekin(:,:)
    real(dp), pointer     :: F3Tot(:), M4Tot(:)
    integer , intent(in)  :: useFNV
    integer , intent(out) :: ierr

    !! Local variables
    integer                  :: i, iorder, nWC, Nfreq
    real(dp), parameter      :: Xsi_p = 4.0_dp
    real(dp)                 :: Rc1, Rc2, uu, ww, dudt, dwdt, dudx, dudz, dudtdz
    real(dp)                 :: eta1, eta2, F2Free, F3Tayl, F3Conv, F3Sec, F3Xsi
    complex(dp), allocatable :: FF(:)
    logical                  :: F2, F3, M4
    character(len=64)        :: msg

    !! --- Logic section ---

    !! Assuming useFNV is a 4- or 5 digit integer.
    !! 1&2: iorder, 3: Bool for 2nd order Force, 4: Bool for 3rd order Force,
    !! 5: Bool for 4th order Moment
    if (useFNV > 99999 .or. useFNV < 1000) then
       write(msg,"('FNV option out range:',I6)") useFNV
       call reportError(error_p,msg)
       ierr = -1
       return
    end if

    M4 = mod(useFNV,10) > 0      ! Last digit
    F3 = mod(useFNV,100)/10 > 0  ! Second last digit
    F2 = mod(useFNV,1000)/100 >0 ! Third last digit
    iorder = useFNV/1000         ! First digit(s) (one or two)
    write(msg,"('FNV options: iorder=',I2,' F2=',L1,' F3=',L1,' M4=',L1,' Cutoff=',F8.3,'rad/s')") iorder,F2,F3,M4,LPcut
    call reportError(note_p,msg)

    nWC   = size(Wavekin,1)
    Nfreq = size(Wavcomp,2)
    allocate(F3Tot(nWC),M4Tot(nWC),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('FNV23order')
       return
    end if

    Rc1 = pi_p*rho * D*D*0.25_dp
    Rc2 = Rc1 * D*0.5_dp

    do i = 1, nWC

       !! Fetch wave kinematics
       uu     = Wavekin(i,1) ! horizontal velocity
       dudt   = Wavekin(i,2) ! horizontal acceleration
       dudtdz = Wavekin(i,3) ! horizontal acceleration gradient
       ww     = Wavekin(i,4) ! vertical velocity
       dudz   = Wavekin(i,5) ! horizontal velocity gradient
       dudx   = Wavekin(i,6) ! horizontal velocity gradient
       dwdt   = Wavekin(i,7) ! vertical acceleration
       eta1   = Wavekin(i,8) ! first-order wave elevation
       eta2   = Wavekin(i,9) ! second-order wave elevation

       !! 2nd-order free surface contribution
       F2Free = (Rc1+Rc1)*dudt*eta1

       !! Free surface Taylor expansion, FNV-term
       F3Tayl = Rc1*dudtdz*eta1*eta1

       !! Free surface convection term, FNV-term
       F3Conv = Rc1*(2.0_dp*ww*dudz + uu*dudx)*eta1

       !! Extended 2nd-order free surface term to 2nd-order wave elevation
       F3Sec  = (Rc1+Rc1)*dudt*eta2

       !! Xsi surface contribution, infinite depth cylinder
       F3Xsi  = Xsi_p*Rc1/g * uu*uu*dudt

       !! Total 3rd-order force
       if (F2) then
          F3Tot(i) = F2Free
       else
          F3Tot(i) = 0.0_dp
       end if
       if (F3) then
          F3Tot(i) = F3Tot(i) + F3Tayl + F3Conv + F3Sec + F3Xsi
       end if

       !! Total 4th-order moment (see equation (4) in the Marintek report)
       if (M4) then
          M4Tot(i) = (Rc1*eta1*eta1 + (Rc2/g)*4.95037_dp*uu*uu)*dudt &
               &   + (Rc1+Rc1)*(-dudt*eta1*dwdt/g + dudtdz*eta1*eta1/3.0_dp + &
               &                ww*dudz*eta1/4.0_dp + uu*uu*dudt*Xsi_p/(g+g))*eta1
       else
          M4Tot(i) = 0.0_dp
       end if
    end do

    !! AR: According to Jørgen, the FNV force should not filtered since the kinematics are filtered
    return

    !! -----------------------------------------
    !! --- Fourier Transform
    !! -----------------------------------------

    if (Nfreq > 2) then ! Skip filtering for regular waves

       allocate(FF(nWC),STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('FNV23order')
          return
       end if

       !! Initialisation
       call initFFT (nWC,ierr)
       if (ierr /= 0) goto 990

       !! 3rd-order force Fourier components
       call doFFT (.false.,F3tot,FF,ierr)
       if (ierr /= 0) goto 990

       !! Perform low-pass filter
       call filter_time (F3tot,samp,2.0_dp*LPcut,iorder,.true.,ierr)
       if (ierr /= 0) goto 990
    end if

    if (Nfreq > 2 .and. M4) then

       !! 4th-order moment Fourier components
       call doFFT (.false.,M4tot,FF,ierr)
       if (ierr /= 0) goto 990

       !! Perform low-pass filter
       call filter_time (M4tot,samp,2.0_dp*LPcut,iorder,.true.,ierr)
       if (ierr /= 0) goto 990
    end if

    if (Nfreq > 2) then
       !! Cleaning up
       call closeFFT (ierr)
       if (ierr /= 0) goto 990

       deallocate(FF)
    end if

    return

990 call reportError (debugFileOnly_p,'FNV23order')
    if (allocated(FF)) deallocate(FF)

  end subroutine FNV23order

end module FNVmodule
