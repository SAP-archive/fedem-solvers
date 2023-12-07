!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module BushingElementRoutinesModule

  implicit none

  private

  public :: addInBushingElementMat, addInBushingElementForces
  public :: updateBushingElements


contains

  subroutine addInBushingElementMat (scaleK,scaleC,Nmat,bElem,sam,ierr,Rhs)

    !!==========================================================================
    !! Calculates the stiffness- and/or damper matrix contributions for the
    !! bushing element and adds them into the system Newton matrix, multiplied
    !! by a factor, scaleK and scaleC, respectively.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 12 July 2002 / 1.0
    !!==========================================================================

    use SamModule               , only : SamType, dp
    use SysMatrixTypeModule     , only : SysMatrixType
    use BushingElementTypeModule, only : BushingElementType
    use TriadTypeModule         , only : transTriadToSys, transTriadToGlob
    use TriadTypeModule         , only : transGlobToSys
    use AsmExtensionModule      , only : csAddEM
    use scratchArrayModule      , only : getRealScratchMatrix
    use reportErrorModule       , only : reportError, debugFileOnly_p

    real(dp)                , intent(in)    :: scaleK, scaleC
    type(SysMatrixType)     , intent(inout) :: Nmat
    type(BushingElementType), intent(in)    :: bElem
    type(SamType)           , intent(in)    :: sam
    integer                 , intent(out)   :: ierr
    real(dp), optional      , intent(inout) :: Rhs(:)

    !! Local variables
    integer           :: i, j
    real(dp), pointer :: eMat(:,:)

    !! --- Logic section ---

    eMat => getRealScratchMatrix(bElem%nDOFs,bElem%nDOFs,ierr)
    if (ierr < 0) goto 900

    eMat = 0.0_dp
    do i = 1, size(bElem%spokes)
       if (scaleK > 0.0_dp) then

          !! Calculate stiffness matrix terms for this spoke
          do j = 1, size(bElem%spokes(i)%springs)
             call AddInSpokeMat (bElem%spokes(i), &
                  &              bElem%spokes(i)%springs(j)%dof, &
                  &              bElem%spokes(i)%springs(j)%stiffness*scaleK, &
                  &              eMat)
          end do

       end if
       if (scaleC > 0.0_dp) then

          !! Calculate damper matrix terms for this spoke
          do j = 1, size(bElem%spokes(i)%dampers)
             call AddInSpokeMat (bElem%spokes(i), &
                  &              bElem%spokes(i)%dampers(j)%dof, &
                  &              bElem%spokes(i)%dampers(j)%coeff*scaleC, &
                  &              eMat)
          end do

       end if
    end do

    !! Transform centre contributions to system directions
    call transTriadToSys (bElem%centre,eMat,1,6)

    !! Transform spoke contributions to system directions (translation only)
    do i = 1, size(bElem%spokes)
       j = bElem%spokes(i)%elmDof
       call transTriadToGlob (bElem%centre,eMat,j,3)
       call transGlobToSys (bElem%spokes(i)%triad,eMat,j,3)
    end do

    !! Add into system matrix
    call csAddEM (sam,bElem%samElnum,bElem%nDOFs,eMat,Nmat,ierr,Rhs)
    if (ierr == 0) return

900 call reportError (debugFileOnly_p,'addInBushingElementMat')

  end subroutine addInBushingElementMat


  subroutine addInSpokeMat (spoke,dof,coeff,elMat)

    !!==========================================================================
    !! Calculates the contributions to the bushing element coefficient matrix
    !! for the given local dof of the given spoke and adds them into elMat.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 12 July 2002 / 1.0
    !!==========================================================================

    use BushingElementTypeModule, only : SpokeType, dp

    type(SpokeType), intent(in)    :: spoke
    integer        , intent(in)    :: dof
    real(dp)       , intent(in)    :: coeff
    real(dp)       , intent(inout) :: elMat(:,:)

    !! Local variables
    integer  :: i, j, k
    real(dp) :: F(3)

    !! --- Logic section ---

    if (dof <= 3) then

       !! Translatory DOF, just add directly in
       i = dof
       j = spoke%elmDof + dof-1
       elMat(i,i) = elMat(i,i) + coeff
       elMat(i,j) = elMat(i,j) - coeff
       elMat(j,i) = elMat(j,i) - coeff
       elMat(j,j) = elMat(j,j) + coeff

    else

       !! Rotational DOF, transform the spoke triad contributions
       !! into equivalent translational terms
       f(1) =  0.0_dp
       f(2) = -coeff*spoke%Tlg(dof-3,2)/spoke%length
       f(3) =  coeff*spoke%Tlg(dof-3,3)/spoke%length
       F    = matmul(spoke%Tlg,f)

       i = dof
       j = spoke%elmDof
       k = j + 2
       elMat(i,i)   = elMat(i,i) + coeff
       elMat(i,j:k) = elMat(i,j:k) + F
       elMat(j:k,i) = elMat(j:k,i) + F
       elMat(i,1:3) = elMat(i,1:3) - F
       elMat(1:3,i) = elMat(1:3,i) - F

    end if

  end subroutine addInSpokeMat


  subroutine addInBushingElementForces (F,RF,bElem,sam,ierr, &
       &                                addSprings,addDampers)

    !!==========================================================================
    !! Calculates the contributions to the system force vector from a bushing
    !! element and adds them into the system vectors F and RF.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 12 July 2002 / 1.0
    !!==========================================================================

    use SamModule               , only : SamType, dp
    use BushingElementTypeModule, only : BushingElementType
    use TriadTypeModule         , only : transTriadToSys, transTriadToGlob
    use TriadTypeModule         , only : transGlobToSys
    use AsmExtensionModule      , only : csAddEV
    use scratchArrayModule      , only : getRealScratchArray
    use reportErrorModule       , only : reportError, debugFileOnly_p

    real(dp)                , intent(inout) :: F(:), RF(:)
    type(BushingElementType), intent(in)    :: bElem
    type(SamType)           , intent(in)    :: sam
    integer                 , intent(inout) :: ierr
    logical, optional       , intent(in)    :: addSprings, addDampers

    !! Local variables
    integer           :: i, j
    logical           :: doAddSprings, doAddDampers
    real(dp), pointer :: eVec(:)

    !! --- Logic section ---

    eVec => getRealScratchArray(bElem%nDOFs,ierr)
    if (ierr < 0) goto 900

    if (present(addSprings)) then
       doAddSprings = addSprings
    else
       doAddSprings = .false.
    end if
    if (present(addDampers)) then
       doAddDampers = addDampers
    else
       doAddDampers = .false.
    end if

    eVec = 0.0_dp
    do i = 1, size(bElem%spokes)
       if (doAddSprings) then

          !! Calculate stiffness forces for this spoke
          do j = 1, size(bElem%spokes(i)%springs)
             call AddInSpokeVec (bElem%spokes(i), &
                  &              bElem%spokes(i)%springs(j)%dof, &
                  &              bElem%spokes(i)%springs(j)%Force, &
                  &              eVec)
          end do

       end if
       if (doAddDampers) then

          !! Calculate damper forces for this spoke
          do j = 1, size(bElem%spokes(i)%dampers)
             call AddInSpokeVec (bElem%spokes(i), &
                  &              bElem%spokes(i)%dampers(j)%dof, &
                  &              bElem%spokes(i)%dampers(j)%Force, &
                  &              eVec)
          end do

       end if
    end do

    !! Transform centre contributions to system directions
    call transTriadToSys (bElem%centre,eVec(1:6))

    !! Transform spoke contributions to system directions (translations only)
    do i = 1, size(bElem%spokes)
       j = bElem%spokes(i)%elmDof
       call transTriadToGlob (bElem%centre,eVec(j:j+2))
       call transGlobToSys (bElem%spokes(i)%triad,eVec(j:j+2))
    end do

    !! Add into system vector
    call csAddEV (sam,bElem%samElnum,bElem%nDOFs,eVec,F,RF,ierr)
    if (ierr == 0) return

900 call reportError (debugFileOnly_p,'addInBushingElementForces')

  end subroutine addInBushingElementForces


  subroutine addInSpokeVec (spoke,dof,force,elVec)

    !!==========================================================================
    !! Calculates the contributions to the bushing element force vector
    !! for the given local dof of the given spoke and adds them into elVec.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 12 July 2002 / 1.0
    !!==========================================================================

    use BushingElementTypeModule, only : SpokeType, dp

    type(SpokeType), intent(in)    :: spoke
    integer        , intent(in)    :: dof
    real(dp)       , intent(in)    :: force
    real(dp)       , intent(inout) :: elVec(:)

    !! Local variables
    integer  :: i, j, k
    real(dp) :: F(3)

    !! --- Logic section ---

    if (dof <= 3) then

       !! Translatory DOF, just add directly in
       i = dof
       j = spoke%elmDof + dof-1
       elVec(i) = elVec(i) + force
       elVec(j) = elVec(j) - force

    else

       !! Rotational DOF, transform the spoke triad contributions
       !! into equivalent translational terms
       f(1) =  0.0_dp
       f(2) = -force*spoke%Tlg(dof-3,2)/spoke%length
       f(3) =  force*spoke%Tlg(dof-3,3)/spoke%length
       F    = matmul(spoke%Tlg,f)

       i = dof
       j = spoke%elmDof
       k = j + 2
       elVec(i)   = elVec(i) + force
       elVec(j:k) = elVec(j:k) + F
       elVec(1:3) = elVec(1:3) - F

    end if

  end subroutine addInSpokeVec


  subroutine updateBushingElements (bElems,ierr,springsOnly,dampersOnly)

    !!==========================================================================
    !! Updates all bushing elements.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 11 July 2002 / 1.0
    !!==========================================================================

    use BushingElementTypeModule, only : BushingElementType

    type(BushingElementType), intent(inout) :: bElems(:)
    integer                 , intent(inout) :: ierr
    logical, optional       , intent(in)    :: springsOnly, dampersOnly

    !! Local variables
    integer :: i
    logical :: doSprings, doDampers

    !! --- Logic section ---

    doSprings = .true.
    doDampers = .true.
    if (present(springsOnly)) then
       if (springsOnly) doDampers = .false.
    end if
    if (present(dampersOnly)) then
       if (dampersOnly) then
          doDampers = .true.
          doSprings = .false.
       end if
    end if

    do i = 1, size(bElems)
       call updateBushingElement (bElems(i),doSprings,doDampers,ierr)
    end do

  end subroutine updateBushingElements


  subroutine updateBushingElement (bElem,doSprings,doDampers,ierr)

    !!==========================================================================
    !! Updates the given bushing element.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 11 July 2002 / 1.0
    !!==========================================================================

    use BushingElementTypeModule, only : BushingElementType, dp
    use IdTypeModule            , only : getId
    use EngineRoutinesModule    , only : updateSpringBase, updateDamperBase
    use reportErrorModule       , only : reportError, error_p

    type(BushingElementType), intent(inout) :: bElem
    logical                 , intent(in)    :: doSprings, doDampers
    integer                 , intent(inout) :: ierr

    !! Local variables
    integer :: i, j, nSpokes, lerr

    !! --- Logic section ---

    nSpokes = size(bElem%spokes)

    !! Update the bushing geometry
    lerr = ierr
    do i = 1, nSpokes
       call updateSpoke (bElem%spokes(i),bElem%centre,doDampers,ierr)
    end do
    if (ierr < lerr) then
       call reportError (error_p,'Failed to update Bushing'//getId(bElem%id))
       return
    end if

    !! Compute scaling factors (TODO: add some intelligence here...)
    do i = 1, nSpokes
       if (doSprings) then
          do j = 1, size(bElem%spokes(i)%springs)
             bElem%spokes(i)%springs(j)%scale1 = 1.0_dp / dble(nSpokes)
          end do
       end if
       if (doDampers) then
          do j = 1, size(bElem%spokes(i)%dampers)
             bElem%spokes(i)%dampers(j)%scale1 = 1.0_dp / dble(nSpokes)
          end do
       end if
    end do

    !! Update the spring and damper variables
    do i = 1, nSpokes
       if (doSprings) then
          do j = 1, size(bElem%spokes(i)%springs)
             call updateSpringBase (bElem%spokes(i)%springs(j),ierr)
          end do
       end if
       if (doDampers) then
          do j = 1, size(bElem%spokes(i)%dampers)
             call updateDamperBase (bElem%spokes(i)%dampers(j),ierr)
          end do
       end if
    end do

    if (ierr < lerr) then
       call reportError (error_p,'Failed to update Bushing'//getId(bElem%id))
    end if

  end subroutine updateBushingElement


  subroutine updateSpoke (spoke,centre,lDynamics,ierr)

    !!==========================================================================
    !! Updates the length and orientation of the given spoke, and computes the
    !! length (and velocity) for all its component springs and dampers.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 11 July 2002 / 1.0
    !!==========================================================================

    use kindModule              , only : dp, epsDiv0_p
    use BushingElementTypeModule, only : SpokeType
    use TriadTypeModule         , only : TriadType, transVGlobToTriad
    use IdTypeModule            , only : getId
    use manipMatrixModule       , only : trans1V
    use rotationModule          , only : deltaRot
    use reportErrorModule       , only : internalError, reportError, error_p

    type(SpokeType), intent(inout) :: spoke
    type(TriadType), intent(in)    :: centre
    logical        , intent(in)    :: lDynamics
    integer        , intent(inout) :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: sDir(6), sVel(6)

    !! --- Logic section ---

    !! Position of the spoke triad relative to the center triad
    sDir(1:3) = spoke%triad%ur(:,4) - centre%ur(:,4)
    sDir(4:6) = deltaRot(centre%ur,spoke%triad%ur)

    if (lDynamics) then

       !! Relative velocity between the two triads
       if (centre%nDOFs == 6 .and. spoke%triad%nDOFs == 6) then
          sVel = spoke%triad%urd - centre%urd
       else if (centre%nDOFs == 6 .and. spoke%triad%nDOFs == 3) then
          sVel(1:3) = spoke%triad%urd(1:3) - centre%urd(1:3)
          sVel(4:6) = -centre%urd(1:3)
       else
          ierr = ierr + internalError('updateSpoke: Invalid triad nDOFs')
          return
       end if

    end if

    !! Length of the spoke
    spoke%length = sqrt(sDir(1)*sDir(1)+sDir(2)*sDir(2)+sDir(3)*sDir(3))
    if (spoke%length > epsDiv0_p) then
       !! Local-to-global transformation matrix for the spoke
       spoke%Tlg = matmul(transpose(centre%ur(:,1:3)),trans1V(sDir))
    else
       call reportError (error_p,'Triad'//trim(getId(centre%id))// &
            ' and Triad'//trim(getId(spoke%triad%id)), &
            'which define a spoke in current Bushing Element coincide.')
       ierr = ierr - 1
       return
    end if

    !! Transform the relative position to the local coordinate system
    !! of the bushing element (equal to that of the centre triad)
    sDir = transVGlobToTriad(centre,sDir)
    !! Assign length to the component springs
    do i = 1, size(spoke%springs)
       spoke%springs(i)%length = sDir(spoke%springs(i)%dof)
    end do

    if (lDynamics) then

       !! Transform the relative velocity to the local coordinate system
       !! of the bushing element (equal to that of the centre triad)
       sVel = transVGlobToTriad(centre,sVel)
       !! Assign length and velocity to the component dampers
       do i = 1, size(spoke%dampers)
          spoke%dampers(i)%length   = sDir(spoke%dampers(i)%dof)
          spoke%dampers(i)%velocity = sVel(spoke%dampers(i)%dof)
       end do

    end if

  end subroutine updateSpoke

end module BushingElementRoutinesModule
