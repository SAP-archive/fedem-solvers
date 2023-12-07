!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module GenericPartModule

  implicit none

  private

  public :: ReadGenericPart, InitiateGenericParts


contains

  subroutine ReadGenericPart (infp,sup,ierr)

    !!==========================================================================
    !! Read stiffness and mass properties for a generic part from the solver
    !! input file and generate the equivalent superelement matrices.
    !! Optionally condense out the center of gravity (CoG) triad.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 16 Jan 2004/1.0
    !!==========================================================================

    use SupElTypeModule  , only : SupElType, dp
    use IdTypeModule     , only : getId, StrId, ReportInputError
    use DenseMatrixModule, only : solveAxB
    use inputUtilities   , only : iuSetPosAtNextEntry
    use reportErrorModule, only : reportError, error_p, debugFileOnly_p
    use reportErrorModule, only : internalError, allocationError

    integer        , intent(in)    :: infp
    type(SupElType), intent(inout) :: sup
    integer        , intent(out)   :: ierr

    !! Local variables
    integer           :: nDOFs
    real(dp), pointer :: Kgp(:,:), Mgp(:,:)
    real(dp)          :: Kii(6,6)

    !! Define the GENERIC_PART namelist
    integer, parameter :: maxTriads_p = 10000
    integer            :: supElId, isRigid
    real(dp)           :: mass, inertia(6), kt, kr, stiffness(6*maxTriads_p)
    namelist /GENERIC_PART/ supElId, isRigid, mass, inertia, kt, kr, stiffness

    !! --- Logic section ---

    if (sup%nExtNods > maxTriads_p) then
       ierr = internalError('ReadGenericPart: Too many triads in generic part')
       return
    else if (.not. iuSetPosAtNextEntry(infp,'&GENERIC_PART')) then
       ierr = 1
       call ReportInputError ('GENERIC_PART')
       call reportError (debugFileOnly_p,'ReadGenericPart')
       return
    end if

    !! Read mass and stiffness properties
    supElId=0; isRigid=0
    mass=0.0_dp; inertia=0.0_dp
    kt=0.0_dp; kr=0.0_dp; stiffness=0.0_dp
    read(infp,nml=GENERIC_PART,iostat=ierr)
    if (ierr /= 0) then
       ierr = 1
       call ReportInputError ('GENERIC_PART',id=sup%id)
       call reportError (debugFileOnly_p,'ReadGenericPart')
       return
    else if (supElId /= sup%id%baseId) then
       ierr = 2
       call ReportInputError ('GENERIC_PART',id=sup%id, &
            &                 msg='Invalid superelement baseId'//StrId(supElId))
       call reportError (debugFileOnly_p,'ReadGenericPart')
       return
    else if (isRigid == 1) then
       !! Estimate stiffness coefficients that makes the superelement rigid
       call EstimateRigidStiffness (sup, mass, inertia, kt, kr)
       sup%stressStiffFlag = 0 ! No stress stiffening
       sup%rigidFlag = 2
    else
       sup%rigidFlag = 1
    end if

    if (sup%triads(1)%p%nDOFs == 0) then
       !! The CoG triad should be condensed out, must allocate the full matrices
       nDOFs = sup%nTotDOFs + 6
       allocate(sup%Bgp(6,sup%nTotDOFs),Kgp(nDOFs,nDOFs),Mgp(6,6),STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('ReadGenericPart')
          return
       end if
       !! Set samNodNum temporarily to flag that the CoG triad should
       !! receive a proper node number later, also when it is condensed out
       sup%triads(1)%p%samNodNum = 1
    else
       !! The CoG triad should be retained in the dynamics
       Kgp => sup%KmMat
       Mgp => sup%Mmat
    end if

    !! Assemble the equivalent superelement stiffness matrix
    call AssembleStiffness (sup, Kgp, kt, kr, &
         &                  reshape(stiffness,(/6,sup%nExtNods-1/)))

    !! Assemble the equivalent superelement mass matrix
    call AssembleMass (Mgp, mass, inertia)

    !! Equivalent gravity force vectors associated with unit accelerations
    sup%fg = 0.0_dp
    sup%fg(1,1) = mass
    sup%fg(2,2) = mass
    sup%fg(3,3) = mass

    if (sup%triads(1)%p%nDOFs > 0) return ! Retain the CoG triad

    !! Eliminate the CoG triad DOFs through static condensation
    Kii =  Kgp(1:6,1:6)
    sup%Bgp = -Kgp(1:6,7:)
    call solveAxB (Kii,sup%Bgp,ierr)
    if (ierr > 0) then
       call reportError (error_p,'Singular stiffness matrix, zero pivot'// &
            &            ' detected in local DOF'//trim(StrId(ierr))// &
            &            ' for Generic Part'//getId(sup%id), &
            &            addString='ReadGenericPart')
    else if (ierr < 0) then
       call reportError (debugFileOnly_p,'ReadGenericPart'//getId(sup%id))
       return
    end if

    !! Reduced stiffness, mass and gravitation vectors
    sup%KmMat = Kgp(7:,7:) + matmul(Kgp(7:,1:6),sup%Bgp)
    sup%Mmat  = matmul(transpose(sup%Bgp),matmul(Mgp,sup%Bgp))
    sup%fg    = matmul(transpose(sup%Bgp),sup%fg(1:6,:))

    !! Release the temporary arrays
    deallocate(Kgp,Mgp)

  end subroutine ReadGenericPart


  subroutine EstimateRigidStiffness (sup,mass,inertia,kt,kr)

    !!==========================================================================
    !! Estimate stiffness coefficients to be assigned at each triad of the
    !! superelement, such that it behaves as completely rigid. The estimate
    !! is based on the given mass and a user-specified target eigenfrequency.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 28 Oct 2004/1.0
    !!==========================================================================

    use kindModule            , only : dp, epsDiv0_p, pi_p
    use SupElTypeModule       , only : SupElType
    use IdTypeModule          , only : getId
    use reportErrorModule     , only : getErrorFile, reportError, warning_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_doubleValue

    type(SupElType), intent(in)  :: sup
    real(dp)       , intent(in)  :: mass, inertia(:)
    real(dp)       , intent(out) :: kt, kr

    !! Local variables
    real(dp) :: mRot, omega

    !! --- Logic section ---

    if (mass < epsDiv0_p) then
       call reportError (warning_p,'Generic Part'//trim(getId(sup%id))// &
            &            ' has no mass and no stiffness')
       kt = 0.0_dp
       kr = 0.0_dp
       return
    end if

    !! Target angular eigenfrequency
    omega = ffa_cmdlinearg_doubleValue('targetFrequencyRigid') * 2.0_dp*pi_p

    !! Characteristic inertia (lumped and averaged)
    mRot = (sum(inertia(1:3)) + sum(inertia(4:6))*2.0_dp) / 3.0_dp

    !! If no inertia, estimate one based on the mass and bounding box size
    if (mRot < epsDiv0_p) mRot = mass * 0.25_dp*getBoundingBoxSize(sup)

    !! If also the bounding box is zero, estimate a size based on the mass
    if (mRot < epsDiv0_p) mRot = mass * 0.25_dp*estimateBallSize(mass)

    !! Compute the stiffness coefficients
    kt = omega*omega * mass
    kr = omega*omega * mRot
    write(getErrorFile(),600) trim(getId(sup%id)), kt, kr

600 format(/5X,'Estimated rigid stiffness coefficient for Generic Part',A &
         & /5X,'kt = ',1PE12.5,' kr = ',E12.5 )

  contains

    !! Computes the size of the bounding box of the given superelement.
    !! The square of the diameter is returned.
    function getBoundingBoxSize (sup)
      type(SupElType), intent(in) :: sup
      real(dp) :: getBoundingBoxSize, BB(3)
      integer  :: i
      do i = 1, 3
         BB(i) = maxval(sup%TrUndeformed(i,4,:)) &
              -  minval(sup%TrUndeformed(i,4,:))
      end do
      getBoundingBoxSize = BB(1)*BB(1) + BB(2)*BB(2) + BB(3)*BB(3)
    end function getBoundingBoxSize

    !! Computes the square of the diameter of a ball with the given mass
    !! and an assumed mass density half that of steel.
    function estimateBallSize (totalMass)
      use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble
      real(dp), intent(in) :: totalMass
      real(dp), parameter  :: density_p = 7850.0_dp ! Assumed steel [kg/m^3]
      real(dp) :: estimateBallSize, scaleToKG, scaleToM, ballVolume
      call ffa_cmdlinearg_getdouble('scaleToKG',scaleToKG)
      call ffa_cmdlinearg_getdouble('scaleToM',scaleToM)
      ballVolume = (totalMass*scaleToKG*2.0_dp/density_p) / scaleToM**3.0_dp
      estimateBallSize = (6.0_dp*ballVolume/pi_p) ** (2.0_dp/3.0_dp)
    end function estimateBallSize

  end subroutine EstimateRigidStiffness


  subroutine AssembleStiffness (sup,K,kt,kr,stiffness)

    !!==========================================================================
    !! Assemble the stiffness matrix for a generic part.
    !!
    !! The part of the stiffness matrix associated with the CoG and one of the
    !! other triads of the generic part is assumed to be:
    !!
    !! [K] = |  kt           -kt*spin(e)               -kt           0  |
    !!       |         kr + spin(e)'*kt*spin(e)   spin(e)'*kt      -kr  |
    !!       |                                          kt           0  |
    !!       |  symm.                                               kr  |
    !!
    !! where kt and kr are diagonal 3x3 stiffness matrices for translational and
    !! rotational DOFs, respectively, and e is the eccentricity vector giving
    !! the relative position of the hard-point with respect to the CoG.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Jan 2004/1.0
    !!==========================================================================

    use SupElTypeModule, only : SupElType, dp
    use rotationModule , only : EccExpand

    type(SupElType), intent(in)  :: sup
    real(dp)       , intent(out) :: K(:,:)
    real(dp)       , intent(in)  :: kt, kr, stiffness(:,:)

    !! Local variables
    integer  :: i, j
    real(dp) :: CG(3), eK(6,6), dK(6)

    !! --- Logic section ---

    K = 0.0_dp
    CG = sup%TrUndeformed(:,4,1) ! Triad 1 is at the CoG

    dK(1:3) = kt; dK(4:6) = kr ! Assume a common stiffness for all triads

    do i = 2, sup%nExtNods

       !! Individual stiffness coefficients?
       if (abs(kt)+abs(kr) <= 1.0e-16_dp) dK = stiffness(:,i-1)

       !! Compute eccentric element stiffness matrix for the nodal spring
       call EccExpand (sup%TrUndeformed(:,4,i)-CG,dK(1:3),eK)
       do j = 4, 6
          eK(j,j) = eK(j,j) + dK(j)
       end do

       !! Add into the superelement stiffness matrix
       K(1:6,1:6) = K(1:6,1:6) + eK
       K(6*i-5:6*i-3,1:6) = -eK(1:3,:)
       K(1:6,6*i-5:6*i-3) = -eK(:,1:3)
       do j = 4, 6
          K(6*i-6+j,j) = -dK(j)
          K(j,6*i-6+j) = -dK(j)
       end do
       do j = 1, 6
          K(6*i-6+j,6*i-6+j) = dK(j)
       end do

    end do

  end subroutine AssembleStiffness


  subroutine AssembleMass (M,mass,inertia)

    !!==========================================================================
    !! Assemble the mass matrix for a generic part.
    !! Assume that all mass is coupled to the first triad (CoG).
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 21 Jan 2004/1.0
    !!==========================================================================

    use KindModule, only : dp

    real(dp), intent(out) :: M(:,:)
    real(dp), intent(in)  :: mass, inertia(:)

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    M = 0.0_dp
    k = 3
    do i = 1, 3
       M(i,i) = mass
       M(3+i,3+i) = inertia(i)
       do j = i+1, 3
          k = k + 1
          M(3+i,3+j) = inertia(k)
          M(3+j,3+i) = inertia(k)
       end do
    end do

  end subroutine AssembleMass


  subroutine InitiateGenericParts (sups)

    !!==========================================================================
    !! Initiates some more generic part data after system initialization.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 28 Aug 2013/1.0
    !!==========================================================================

    use SupElTypeModule  , only : SupElType
    use IdTypeModule     , only : getId
    use reportErrorModule, only : getErrorFile

    type(SupElType), intent(inout) :: sups(:)

    !! Local variables
    integer :: i, j, lpu

    !! --- Logic section ---

    lpu = getErrorFile()
    do i = 1, size(sups)
       if ( sups(i)%rigidFlag > 0 .and. size(sups(i)%triads) > 1 .and. &
            sups(i)%triads(1)%p%nDOFs > 0 .and. &
            sups(i)%triads(1)%p%id%baseid < 0 ) then

          !! Calculate initial velocity and acceleration at the CoG triad
          !! by averaging the values at the other triads
          call findCoGVelAcc(sups(i)%triads(1)%p%urd, &
               &             sups(i)%triads(1)%p%urdd, &
               &             sups(i)%triads)

          write(lpu,600) 'Generic Part', trim(getId(sups(i)%id))
          write(lpu,610) 'Velocity     :', sups(i)%triads(1)%p%urd
          write(lpu,610) 'Acceleration :', sups(i)%triads(1)%p%urdd

       end if
    end do

600 format(/5X,'>>> Computed initial conditions at centre of gravity for ',2A )
610 format( 9X,A,1P10E13.5 )

  contains

    subroutine findCoGVelAcc (urd,urdd,triads)
      use TriadTypeModule, only : TriadPtrType, dp
      real(dp)          , intent(out) :: urd(:), urdd(:)
      type(TriadPtrType), intent(in)  :: triads(:)
      urd  = 0.0_dp
      urdd = 0.0_dp
      do j = 2, size(triads)
         if (triads(j)%p%nDOFs >= 3) then
            urd(1:3)  = urd(1:3)  + triads(j)%p%urd(1:3)
            urdd(1:3) = urdd(1:3) + triads(j)%p%urdd(1:3)
         end if
      end do
      urd  = urd  / real(size(triads)-1,dp)
      urdd = urdd / real(size(triads)-1,dp)
    end subroutine findCoGVelAcc

  end subroutine InitiateGenericParts

end module GenericPartModule
