!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module CorotUtilModule

  implicit none

  private

  public :: addKgrToKm, formShadowPosGrad


contains

  subroutine addKgrToKm (supEl,ktan,km,fint,ierr)

    !---------------------------------------------------------------------------
    ! Purpose:
    !     Calculate rot. geometric stiffness for the superelement based on the
    !     variation of its 3 reference nodes and add it to material stiffness.
    !
    ! Input:
    !     supEl : superelement data structure
    !     km    : material stiffness
    !     fint  : superelement force vector
    !
    ! Output:
    !     ktan  : tangential stiffness (material + rot. geometric stiffness)
    !     ierr  : error flag
    !---------------------------------------------------------------------------
    ! Programmer: Dag R Chr                                       date: 23/10/97
    !             Bjorn Haugen (significant upgrade)                    29/10/03
    !---------------------------------------------------------------------------

    use SupElTypeModule   , only : SupElType, dp
    use scratchArrayModule, only : getRealScratchMatrix
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SupElType), intent(in)  :: supEl
    real(dp)       , intent(in)  :: km(:,:), fint(:)
    real(dp)       , intent(out) :: ktan(:,:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer           :: i, j, ndof
    real(dp)          :: kgr
    real(dp), pointer :: fnm(:,:)

    !! --- Logic section ---

    ktan = 0.0_dp
    ndof = supel%nTotDOFs ! Number of triad dofs
    if (associated(supel%genDOFs)) then
       ndof = ndof - supel%genDOFs%nDOFs
    end if
    fnm => getRealScratchMatrix(3,ndof,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'addKgrToKm')
       return
    end if

    !! Calculate Fmn
    call formfnm (supEl%triads,fint,fnm)

    !! Multiply Kgr = Fnm*Gmat (stored in ktan)
    do i = 1, ndof
       do j = 1, ndof
          ktan(i,j) = dot_product(fnm(:,i),supEl%shadowPosGrad(4:6,j))
       end do
    end do

    !! Add kgr symmetrically to km, results in ktan
    do j = 1, ndof
       ktan(j,j) = km(j,j) - ktan(j,j)
       do i = 1, j-1
          kgr = 0.5_dp*(ktan(i,j)+ktan(j,i))
          ktan(i,j) = km(i,j) - kgr
          ktan(j,i) = km(j,i) - kgr
       end do
    end do

    if (supEl%nTotDofs > ndof) then

       !! Add the generalized DOFs (only diagonal terms are nonzero)
       ktan(:,ndof+1:) = 0.0_dp
       ktan(ndof+1:,:) = 0.0_dp
       do j = ndof+1, supEl%nTotDofs
          ktan(j,j) = km(j,j)
       end do

    end if

  end subroutine addKgrToKm


  subroutine gmat3nod (gmat,x,y,z,ierr)

    !---------------------------------------------------------------------------
    ! Purpose:
    !     Compute the rotation gradient matrix for a 3 node element.
    !
    ! Input:
    !     x,y,z : coordinates of the 3 nodes in global coordinate system
    !
    ! Output:
    !     gmat  : rotation gradients of rot_x,rot_y,rot_z in local system
    !             with respect to the local coordinate nodal degrees of freedom
    !             The nodal dof ordering for each node is: tx ty tz
    !     ierr  : error flag
    !---------------------------------------------------------------------------
    ! Programmer: Bjorn Haugen                                    date: 22/10/97
    !---------------------------------------------------------------------------

    use kindModule       , only : dp
    use manipMatrixModule, only : invert34, trans3P
    use reportErrorModule, only : reportError, error_p

    real(dp), intent(in)  :: x(3), y(3), z(3)
    real(dp), intent(out) :: gmat(3,3,3)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: tr(3,4), xloc(3), yloc(3), scr(3,3)

    !! --- Logic section ---

    do i = 1, 3
       scr(1,i) = x(i)
       scr(2,i) = y(i)
       scr(3,i) = z(i)
    end do

    !! Establish local coordinate system for the triangle
    tr = invert34(trans3P(scr(:,1),scr(:,2),scr(:,3),ierr=ierr))
    if (ierr < 0) then
       call reportError (error_p,'Cannot compute rotation gradient matrix', &
            &            'which is used in geometric stiffness calculation.', &
            &            'Check definition of superelement coordinate system', &
            &            addString='gmat3nod')
       return
    end if

    !! Compute coordinates in local system
    do i = 1, 3
       xloc(i) = tr(1,1)*x(i) + tr(1,2)*y(i) + tr(1,3)*z(i) + tr(1,4)
       yloc(i) = tr(2,1)*x(i) + tr(2,2)*y(i) + tr(2,3)*z(i) + tr(2,4)
    end do

    !! Form rotation gradient matrix in local system
    call gmat3nodLocal (gmat,xloc,yloc)
    do i = 1, 3
       scr = matmul(transpose(tr(:,1:3)),gmat(:,:,i))
       gmat(:,:,i) = matmul(scr,tr(:,1:3))
    end do

  end subroutine gmat3nod


  subroutine gmat3nodLocal (gmat,x,y)

    !---------------------------------------------------------------------------
    ! Purpose:
    !     Compute the rotation gradient matrix for a 3 node element.
    !     Element assumed in xy-plane, x-axis not neccesarily along side 1-2.
    !
    ! Input:
    !     x,y   : coordinates for the 3 nodes in local element coordinate system
    !
    ! Output:
    !     gmat  : rotation gradients of rot_x,rot_y,rot_z in local system
    !             with respect to the local coordinate nodal degrees of freedom
    !             The nodal dof ordering for each node is: tx ty tz
    !---------------------------------------------------------------------------
    ! Programmer: Dag R Chr                                       date: 22/10/97
    !---------------------------------------------------------------------------

    use kindModule, only : dp

    real(dp), intent(in)  :: x(3), y(3)
    real(dp), intent(out) :: gmat(3,9)

    !! Local variables
    real(dp) :: area2

    !! --- Logic section ---

    gmat = 0.0_dp

    !! Compute 2 x Area of the triangle
    area2 = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))

    !! Node 1
    gmat(1,3) = (x(3)-x(2))/area2
    gmat(2,3) = (y(3)-y(2))/area2
    gmat(3,2) = -1.0_dp/(x(2)-x(1))

    !! Node 2
    gmat(1,6) = (x(1)-x(3))/area2
    gmat(2,6) = (y(1)-y(3))/area2
    gmat(3,5) = 1.0_dp/(x(2)-x(1))

    !! Node 3
    gmat(1,9) = (x(2)-x(1))/area2
    gmat(2,9) = (y(2)-y(1))/area2

  end subroutine gmat3nodLocal


  subroutine mmat3nod (mmat,emat)

    !---------------------------------------------------------------------------
    ! Purpose:
    !     Calcutate the m-matrix in vref = m * vsup in local coordinates
    !     vref(9,1) - trans dofs only
    !
    ! Input:
    !     emat  : eccentricity matrix in local coordinates
    !
    ! Output:
    !     mmat  : m-matrix
    !---------------------------------------------------------------------------
    ! Programmer: Dag R Chr                                       date: 23/10/97
    !---------------------------------------------------------------------------

    use kindModule, only : dp

    real(dp), intent(out) :: mmat(3,18)
    real(dp), intent(in)  :: emat(3,3)

    !! Local variables
    integer :: jm, refnod

    !! --- Logic section ---

    mmat = 0.0_dp

    do refnod = 1, 3
       jm = (refnod-1)*6

       mmat(1,jm+1) = 1.0_dp
       mmat(2,jm+2) = 1.0_dp
       mmat(3,jm+3) = 1.0_dp

       mmat(1,jm+5) =  emat(3,refnod)
       mmat(1,jm+6) = -emat(2,refnod)
       mmat(2,jm+6) =  emat(1,refnod)

       mmat(2,jm+4) = -mmat(1,jm+5)
       mmat(3,jm+4) = -mmat(1,jm+6)
       mmat(3,jm+5) = -mmat(2,jm+6)

    end do

  end subroutine mmat3nod


  subroutine formfnm (triads,fint,fnm)

    !***************************************************************************
    !*
    !*  Purpose:
    !*     Form the Fn-matrix based on the internal force vector,
    !*     i.e., Fspin with both axial and moment contributions.
    !*
    !*  Input:
    !*     triads   : external nodes in superelement
    !*     fint     : internal force vector
    !*
    !*  Output:
    !*     fnm      : matrix of | Spin(n) | for each node ordered as columns
    !*                          | Spin(m) |
    !*
    !*  Coded by: Bjorn Haugen
    !***************************************************************************

    use TriadTypeModule, only : TriadPtrType, dp

    type(TriadPtrType), intent(in)  :: triads(:)
    real(dp)          , intent(in)  :: Fint(:)
    real(dp)          , intent(out) :: fnm(:,:)

    !! Local variables
    integer :: i, node, dofS, dofE

    !! --- Logic section ---

    fnm = 0.0_dp

    do node = 1, size(triads)
       dofS = triads(node)%firstDOF - 1
       dofE = dofS + triads(node)%p%nDOFs - 1
       do i = dofS, dofE, 3

          fnm(3,i+3) =  0.0_dp
          fnm(3,i+2) = -fint(i+1)
          fnm(3,i+1) =  fint(i+2)
          fnm(2,i+1) = -fint(i+3)
          fnm(2,i+2) =  0.0_dp
          fnm(2,i+3) =  fint(i+1)
          fnm(1,i+3) = -fint(i+2)
          fnm(1,i+2) =  fint(i+3)
          fnm(1,i+1) =  0.0_dp

       end do
    end do

  end subroutine formfnm


  subroutine formShadowPosGrad (supEl,dbgShadowPos,ierr)

    !!==========================================================================
    !! Forms a 6xNDOF gradient matrix used in the shadow position calculations.
    !!
    !! Programmer : Bjorn Haugen                           date/rev Nov 2003/1.0
    !!==========================================================================

    use kindModule                , only : dp, epsDiv0_p
    use SupElTypeModule           , only : SupElType
    use IdTypeModule              , only : getId
    use MassMatrixCorrectionModule, only : mmcRigAccelVectors
    use DenseMatrixModule         , only : solveAxB
    use rotationModule            , only : Spin
    use scratchArrayModule        , only : getRealScratchMatrix
    use manipMatrixModule         , only : invert34, matmul34, WriteObject
    use reportErrorModule         , only : error_p, empty_p, debugFileOnly_p
    use reportErrorModule         , only : reportError, AllocationError

    type(SupElType), intent(inout) :: supEl
    integer        , intent(in)    :: dbgShadowPos
    integer        , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, jg, jm, iNod, iRef, dofStart, dofEnd
    real(dp) :: wTran, wRot, Kmat(6,6), x(3), y(3), z(3)
    real(dp) :: eVec(3), eMat(3,3), gmat(3,9), mmat(3,18), invSuptr(3,4)
    real(dp), pointer :: rigDisp(:,:)
    character(len=2) :: cdof

    !! --- Logic section ---

    allocate(supEl%shadowPosGrad(6,supEl%nTotDOFs),stat=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('formShadowPosGrad')
       return
    end if

    if (dbgShadowPos > 0 .and. supEl%shadowPosAlg > 1) then
       write(dbgShadowPos,*)
       call WriteObject (supEl%Mmat,dbgShadowPos, &
            &            'Mass matrix for Part'//getId(supEl%id))
    end if

    supEl%shadowPosGrad = 0.0_dp
    select case (supEl%shadowPosAlg)

    case (1) ! ========== Shadow positioning algorithm based on 3-node triangle

       !! Need invSupTr for both emat and x,y,z
       invSupTr = invert34(supEl%supTr)

       !! Find emat and x, y and z
       do iRef = 1, 3
          inod = supEl%refTriad(iRef)
          !! e-vector in global system
          eVec = matmul(supEl%triads(inod)%p%ur(:,1:3),supEl%offset(:,iRef))
          !! e-vector in superelement system
          eMat(:,iRef) = matmul(invSupTr(:,1:3),eVec)
          !! reference point in global system
          eVec = eVec + supEl%triads(inod)%p%ur(:,4)
          !! reference point in superelement system
          eVec = matmul34(invSupTr,eVec)
          x(iRef) = eVec(1)
          y(iRef) = eVec(2)
          z(iRef) = eVec(3)
       end do

       !! Calculate gsupmat
       call gmat3nod (gmat,x,y,z,ierr)
       if (ierr < 0) goto 900

       call mmat3nod (mmat,emat)

       if (dbgShadowPos > 0) then
          call WriteObject (gmat,dbgShadowPos,'gmat')
          call WriteObject (mmat,dbgShadowPos,'mmat')
       end if

       do iRef = 1, 3
          inod = supEl%refTriad(iRef)
          dofStart = supEl%triads(inod)%firstDOF
          dofEnd   = dofStart + supEl%triads(inod)%p%nDOFs - 1
          do j = dofStart, dofEnd
             do i = 1, 3
                jg = (iRef-1)*3 + 1
                jm = (iRef-1)*6 + 1 + (j-dofStart)
                wRot = dot_product(gmat(i,jg:jg+2),mmat(:,jm))
                supEl%shadowPosGrad(i+3,j) = supEl%shadowPosGrad(i+3,j) + wRot
             end do
          end do
       end do

    case (2) ! ========== Shadow positioning algorithm based on nodal masses

       if (supEl%mass < epsDiv0_p) goto 910

       !! TODO,bh: add scaling of rotation contribution

       do inod = 1, supEl%nExtNods
          dofStart = supEl%triads(inod)%firstDOF
          dofEnd   = dofStart + supEl%triads(inod)%p%nDOFs - 1

          !! Get mass for this node
          !! Weight for translation
          wTran = ( supEl%Mmat(dofStart  ,dofStart  ) &
               &  + supEl%Mmat(dofStart+1,dofStart+1) &
               &  + supEl%Mmat(dofStart+2,dofStart+2) )/3.0_dp

          !! Translation component of a nodal translation
          do i = 1, 3
             supEl%shadowPosGrad(i,dofStart+i-1) = wTran
          end do

          !! Rotation component of a nodal translation
          eVec = supEl%TrUndeformed(:,4,inod) - supEl%posMassCenter
          call Spin (eVec,eMat)
          supEl%shadowPosGrad(4:6,dofStart:dofStart+2) = wTran*eMat

          if (dofEnd - dofStart >= 5) then
             !! Weight for rotation
             wRot = ( supEl%Mmat(dofStart+3,dofStart+3) &
                  & + supEl%Mmat(dofStart+4,dofStart+4) &
                  & + supEl%Mmat(dofStart+5,dofStart+5) )/3.0_dp

             do j = dofStart+3, dofEnd
                supEl%shadowPosGrad(j-dofStart+1,j) = wRot
             end do
          end if

       end do

       rigDisp => getRealScratchMatrix(supEl%nTotDofs,6,ierr)
       if (ierr < 0) goto 905

       call mmcRigAccelVectors (supEl,supEl%posMassCenter,rigDisp)
       Kmat = matmul(supEl%shadowPosGrad,rigDisp)

       if (dbgShadowPos > 0) then
          call WriteObject (rigdisp,dbgShadowPos,'Rigdisp')
          call WriteObject (supEl%shadowPosGrad,dbgShadowPos,'shadowPosGrad')
          call WriteObject (Kmat,dbgShadowPos,'Kmat for shadowPos, alg = 2')
       end if

       call solveAxB (Kmat,supEl%shadowPosGrad,ierr)
       if (ierr /= 0) goto 900

    case (3) ! ========== Shadow positioning algorithm based on zero rigid body
             !            component of inertia forces

       if (supEl%mass < epsDiv0_p) goto 910

       rigDisp => getRealScratchMatrix(supEl%nTotDofs,6,ierr)
       if (ierr < 0) goto 905

       call mmcRigAccelVectors (supEl,supEl%posMassCenter,rigDisp)
       supEl%shadowPosGrad = matmul(supEl%Mmat,rigdisp)
       Kmat = matmul(supEl%shadowPosGrad,rigDisp)

       if (dbgShadowPos > 0) then
          call WriteObject (rigdisp,dbgShadowPos,'Rigdisp')
          call WriteObject (supEl%shadowPosGrad,dbgShadowPos,'shadowPosGrad')
          call WriteObject (Kmat,dbgShadowPos,'Kmat for shadowPos, alg = 3')
       end if

       call solveAxB (Kmat,supEl%shadowPosGrad,ierr)
       if (ierr /= 0) goto 900

    end select

    if (dbgShadowPos > 0) then
       call WriteObject (supEl%shadowPosGrad,dbgShadowPos, &
            &            'shadowPosGrad for Part'//getId(supEl%id))
    end if

    return

900 if (ierr > 3) then
       cdof = 'R'//char(ichar('T')+ierr)
    else if (ierr > 0) then
       cdof = 'T'//char(ichar('W')+ierr)
    else
       call reportError (empty_p,'Cannot update reference coordinate system'// &
            &            ' for Part'//getId(supEl%id))
       goto 905
    end if
    call reportError (error_p,'Singular stiffness matrix detected when comp'// &
         &            'uting updated','reference coordinate system for Part'// &
         &            getId(supEl%id),'Detected zero pivot at the '// &
         &            cdof //' centre of gravity DOF for this part.', &
         &            'Please check mass matrix contributions to this DOF.')
905 call reportError (debugFileOnly_p,'formShadowPosGrad')
    return

910 ierr = 1
    call reportError (error_p,'The mass-based shadow positioning algorithm '// &
         'can not be used','for the mass-less Part'//getId(supEl%id), &
         'Switch to the triangle-based algorithm instead.', &
         addString='formShadowPosGrad')

  end subroutine formShadowPosGrad

end module CorotUtilModule
