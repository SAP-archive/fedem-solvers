!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file wavgmModule.f90
!>
!> @brief Weighted Average Motion constraint handling.

!!==============================================================================
!> @brief Module with subroutine for setting up linear multi-point constraints.
!>
!> @details This module contains a subroutine for processing the Weighted
!> Average Motion (WAVGM) elements (also known as RBE3 elements in Nastran).

module WAvgMotionModule

  implicit none

  private

  public :: wavgmConstrEqn


contains

  !!============================================================================
  !> @brief Computes constraint equation coefficients for a WAVGM element.
  !>
  !> @param[in] iel Element index
  !> @param[in] lDof Local index of dependent DOF to compute coefficients for
  !> @param[in] nM Number of independent nodes in current element
  !> @param[in] indC Nodal component indices (common for all nodes)
  !> @param[in] txc Table of x-coordinates for the nodal points
  !> @param[in] tyc Table of Y-coordinates for the nodal points
  !> @param[in] tzc Table of Z-coordinates for the nodal points
  !> @param[in] mnpc Matrix of nodal point correspondances for current element
  !> @param[in] epsX Relative geometric tolerance for WAVGM elements
  !> @param tolX Absolute geometric tolerance for each coordinate direction
  !> @param dX Element nodal coordinates relative to the centre of gravity
  !> @param work Work array
  !> @param[in] weight Independent DOF weights for current element
  !> @param[out] omega Resulting constraint equation coefficients
  !> @param[in] ipsw Print switch
  !> @param[in] lpu File unit number for res-file output
  !>
  !> @details The coefficients that couples a dependent DOF at the reference
  !> node of a Weighted AVerage Motion element to the independent nodal DOFs
  !> of the same element are computed.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 24 Sep 2002

  subroutine wavgmConstrEqn (iel,lDof,nM,indC,txc,tyc,tzc,mnpc, &
       &                     epsX,tolX,dX,work,weight,omega,ipsw,lpu)

    use kindModule            , only : dp, epsDiv0_p
    use manipMatrixModule     , only : writeObject
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer , intent(in)    :: iel, lDof, nM, mnpc(:), indC(:), ipsw, lpu
    real(dp), intent(in)    :: txc(:), tyc(:), tzc(:), weight(:), epsX
    real(dp), intent(inout) :: tolX(3), dX(:,:), work(:)
    real(dp), intent(out)   :: omega(:)

    !! Local variables
    logical           :: reComputeCG
    integer           :: i, j, k, iFrst, iLast, nndof
    integer, save     :: lastIC, lastEl = 0
    real(dp)          :: sumWM, tolWM
    character(len=32) :: label

    !! --- Logic section ---

    if (ffa_cmdlinearg_isTrue('useOldWAVGM')) then
       !! Use the old formulation, retained for testing purpose only
       call wavgmConstrEqn_old (lDof,nM,indC,tolX,dX,weight,omega,ipsw,lpu)
       return
    end if

    omega = 0.0_dp
    nndof = size(indC)
    if (lDof < 1 .or. lDof > nndof) return

    if (iel /= lastEl) then
       lastEl = iel
       lastIC = 0
    end if

    !!                        k
    !! Cyclic permutation:   / \
    !!                      i - j
    i = 1 + mod(lDof-1,3)
    j = 1 + mod(i,3)
    k = 1 + mod(j,3)

    if (indC(lDof) > 0) then

       !! Direct coupling of the lDof'th dependent DOF to
       !! corresponding independent DOFs

       iFrst = indC(lDof)
       iLast = iFrst + nM-1
       sumWM = sum(weight(iFrst:iLast))
       if (sumWM > epsDiv0_p) then
          call DAXPY (nM,1.0_dp/sumWM,weight(iFrst),1,omega(lDof),nndof)
       else if (ipsw > 0) then
          write(lpu,600) iel,lDof,lDof,sumWM,epsDiv0_p
       end if

    else if (indC(lDof) < 0) then

       !! Direct coupling with uniform weights
       call DCOPY (nM,1.0_dp/real(nM,dp),0,omega(lDof),nndof)

    end if

    reComputeCG = .not. ffa_cmdlinearg_isTrue('useUnWeightedCG')
    if (lDof > 3 .and. ffa_cmdlinearg_isTrue('useRotCpl')) then

       !! Coupling of rotational dependent DOF to translational DOFs j,k

       if (indC(j) /= 0) then
          iFrst = indC(j)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             if (iFrst > 0) then
                call computePointCoords (weight(iFrst:iLast))
             else ! uniform weights
                call computePointCoords ()
             end if
          end if
          if (iFrst > 0) then
             sumWM = sumSqrW(dX(j,2:1+nM),dX(k,2:1+nM),weight(iFrst:iLast))
             work(1:nM) = dX(k,2:1+nM)*weight(iFrst:iLast)
          else ! uniform weights
             sumWM = sumSqr(dX(j,2:1+nM),dX(k,2:1+nM)) * real(nM,dp)
             call DCOPY (nM,dX(k,2),size(dX,1),work(1),1)
          end if
          tolWM = tolX(j)*tolX(j) + tolX(k)*tolX(k)
          if (sumWM > tolWM) then
             call DAXPY (nM,-1.0_dp/sumWM,work(1),1,omega(j),nndof)
          else if (ipsw > 0) then
             write(lpu,600) iel,lDof,j,sumWM,tolWM
             write(lpu,620) work(1:nM) / sumWM
          end if
       end if
       if (indC(k) /= 0) then
          iFrst = indC(k)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             if (iFrst > 0) then
                call computePointCoords (weight(iFrst:iLast))
             else ! uniform weights
                call computePointCoords ()
             end if
          end if
          if (iFrst > 0) then
             sumWM = sumSqrW(dX(j,2:1+nM),dX(k,2:1+nM),weight(iFrst:iLast))
             work(1:nM) = dX(j,2:1+nM)*weight(iFrst:iLast)
          else ! uniform weights
             sumWM = sumSqr(dX(j,2:1+nM),dX(k,2:1+nM)) * real(nM,dp)
             call DCOPY (nM,dX(j,2),size(dX,1),work(1),1)
          end if
          tolWM = tolX(j)*tolX(j) + tolX(k)*tolX(k)
          if (sumWM > tolWM) then
             call DAXPY (nM,1.0_dp/sumWM,work(1),1,omega(k),nndof)
          else if (ipsw > 0) then
             write(lpu,600) iel,lDof,k,sumWM,tolWM
             write(lpu,620) work(1:nM) / sumWM
          end if
       end if

    else if (lDof < 4 .and. ffa_cmdlinearg_isTrue('useEccCpl')) then

       !! Coupling of translational dependent DOF to translational DOFs i,j
       !! and rotational DOF 3+k through the eccentricity e_j
       !! of the dependent node w.r.t. the independent node centroid

       if (indC(i) /= 0) then
          iFrst = indC(i)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             if (iFrst > 0) then
                call computePointCoords (weight(iFrst:iLast))
             else ! uniform weights
                call computePointCoords ()
             end if
          end if
          if (abs(dX(j,1)) > tolX(j)) then
             if (iFrst > 0) then
                sumWM = sumSqrW(dX(i,2:1+nM),dX(j,2:1+nM),weight(iFrst:iLast))
                work(1:nM) = dX(j,2:1+nM)*weight(iFrst:iLast)
             else
                sumWM = sumSqr(dX(i,2:1+nM),dX(j,2:1+nM)) * real(nM,dp)
                call DCOPY (nM,dX(j,2),size(dX,1),work(1),1)
             end if
             tolWM = tolX(i)*tolX(i) + tolX(j)*tolX(j)
             if (sumWM > tolWM) then
                call DAXPY (nM,dX(j,1)/sumWM,work(1),1,omega(i),nndof)
             else if (ipsw > 0) then
                write(lpu,600) iel,lDof,i,sumWM,tolWM
                write(lpu,620) work(1:nM) * dX(j,1)/sumWM
             end if
          else if (ipsw > 0) then
             write(lpu,610) iel,lDof,i,dX(j,1),tolX(j)
          end if
       end if
       if (indC(j) /= 0) then
          iFrst = indC(j)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             if (iFrst > 0) then
                call computePointCoords (weight(iFrst:iLast))
             else ! uniform weights
                call computePointCoords ()
             end if
          end if
          if (abs(dX(j,1)) > tolX(j)) then
             if (iFrst > 0) then
                sumWM = sumSqrW(dX(i,2:1+nM),dX(j,2:1+nM),weight(iFrst:iLast))
                work(1:nM) = dX(i,2:1+nM)*weight(iFrst:iLast)
             else
                sumWM = sumSqr(dX(i,2:1+nM),dX(j,2:1+nM)) * real(nM,dp)
                call DCOPY (nM,dX(i,2),size(dX,1),work(1),1)
             end if
             if (sumWM > tolWM) then
                call DAXPY (nM,-dX(j,1)/sumWM,work(1),1,omega(j),nndof)
             else if (ipsw > 0) then
                write(lpu,600) iel,lDof,j,sumWM,tolWM
                write(lpu,620) work(1:nM) * dX(j,1)/sumWM
             end if
          else if (ipsw > 0) then
             write(lpu,610) iel,lDof,j,dX(j,1),tolX(j)
          end if
       end if
       if (3+k <= nndof .and. indC(3+k) > 0) then
          iFrst = indC(3+k)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             call computePointCoords (weight(iFrst:iLast))
          end if
          if (abs(dX(j,1)) > tolX(j)) then
             sumWM = sum(weight(iFrst:iLast))
             if (sumWM > epsDiv0_p) then
                call DAXPY (nM,-dX(j,1)/sumWM,weight(iFrst),1,omega(3+k),nndof)
             else if (ipsw > 0) then
                write(lpu,600) iel,lDof,3+k,sumWM,epsDiv0_p
             end if
          else if (ipsw > 0) then
             write(lpu,610) iel,lDof,3+k,dX(j,1),tolX(j)
          end if
       end if

       !! Coupling of translational dependent DOF to translational DOFs i,k
       !! and rotational DOF 3+j through the eccentricity e_k
       !! of the dependent node w.r.t. the independent node centroid

       if (indC(i) /= 0) then
          iFrst = indC(i)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             if (iFrst > 0) then
                call computePointCoords (weight(iFrst:iLast))
             else ! uniform weights
                call computePointCoords ()
             end if
          end if
          if (abs(dX(k,1)) > tolX(k)) then
             if (iFrst > 0) then
                sumWM = sumSqrW(dX(k,2:1+nM),dX(i,2:1+nM),weight(iFrst:iLast))
                work(1:nM) = dX(k,2:1+nM)*weight(iFrst:iLast)
             else
                sumWM = sumSqr(dX(k,2:1+nM),dX(i,2:1+nM)) * real(nM,dp)
                call DCOPY (nM,dX(k,2),size(dX,1),work(1),1)
             end if
             tolWM = tolX(k)*tolX(k) + tolX(i)*tolX(i)
             if (sumWM > tolWM) then
                call DAXPY (nM,dX(k,1)/sumWM,work(1),1,omega(i),nndof)
             else if (ipsw > 0) then
                write(lpu,600) iel,lDof,i,sumWM,tolWM
                write(lpu,620) work(1:nM) * dX(k,1)/sumWM
             end if
          else if (ipsw > 0) then
             write(lpu,610) iel,lDof,i,dX(k,1),tolX(k)
          end if
       end if
       if (indC(k) /= 0) then
          iFrst = indC(k)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             if (iFrst > 0) then
                call computePointCoords (weight(iFrst:iLast))
             else ! uniform weights
                call computePointCoords ()
             end if
          end if
          if (abs(dX(k,1)) > tolX(k)) then
             if (iFrst > 0) then
                sumWM = sumSqrW(dX(k,2:1+nM),dX(i,2:1+nM),weight(iFrst:iLast))
                work(1:nM) = dX(i,2:1+nM)*weight(iFrst:iLast)
             else
                sumWM = sumSqr(dX(k,2:1+nM),dX(i,2:1+nM)) * real(nM,dp)
                call DCOPY (nM,dX(i,2),size(dX,1),work(1),1)
             end if
             tolWM = tolX(k)*tolX(k) + tolX(i)*tolX(i)
             if (sumWM > tolWM) then
                call DAXPY (nM,-dX(k,1)/sumWM,work(1),1,omega(k),nndof)
             else if (ipsw > 0) then
                write(lpu,600) iel,lDof,k,sumWM,tolWM
                write(lpu,620) work(1:nM) * dX(k,1)/sumWM
             end if
          else if (ipsw > 0) then
             write(lpu,610) iel,lDof,k,dX(k,1),tolX(k)
          end if
       end if
       if (3+j <= nndof .and. indC(3+j) > 0) then
          iFrst = indC(3+j)
          iLast = iFrst + nM-1
          if (reComputeCG .and. iFrst /= lastIC) then
             lastIC = iFrst
             call computePointCoords (weight(iFrst:iLast))
          end if
          if (abs(dX(k,1)) > tolX(k)) then
             sumWM = sum(weight(iFrst:iLast))
             if (sumWM > epsDiv0_p) then
                call DAXPY (nM,dX(k,1)/sumWM,weight(iFrst),1,omega(3+j),nndof)
             else if (ipsw > 0) then
                write(lpu,600) iel,lDof,3+j,sumWM,epsDiv0_p
             end if
          else if (ipsw > 0) then
             write(lpu,610) iel,lDof,3+j,dX(k,1),tolX(k)
          end if
       end if

    end if

    if (ipsw > 4) then
       write(lpu,*)
       write(label,"('     Omega for dependent DOF',I2,' :')") lDof
       call writeObject(omega,lpu,label,nndof)
    end if

600 format('  ** Warning: WAVGM element',I10, &
         & ': Ignored coupling by dependent DOF',I2,' to independent DOFs',I2 &
         / 14X, 'due to small denominator',1PE13.5,'  tolerance =',E12.5 )
610 format('  ** Warning: WAVGM element',I10, &
         & ': Ignored coupling by dependent DOF',I2,' to independent DOFs',I2 &
         / 14X, 'due to small eccentricity',1PE13.5,'  tolerance =',E12.5 )
620 format(14X,1P10E13.5/(14X,1P10E13.5))

  contains

    !> @brief Computes the nodal coordinates relative to the centre of gravity.
    !> @details If the weights (W) are present, a weighted CoG is used.
    !> This subroutine also recomputes the geometric tolerances (tolX).
    subroutine computePointCoords (W)
      real(dp), intent(in), optional :: W(:)
      integer  :: n
      real(dp) :: X0(4)
      X0 = 0.0_dp
      if (present(W)) then
         do n = 2, size(mnpc)
            X0(1) = X0(1) + txc(mnpc(n))*W(n-1)
            X0(2) = X0(2) + tyc(mnpc(n))*W(n-1)
            X0(3) = X0(3) + tzc(mnpc(n))*W(n-1)
         end do
         X0(4) = sum(W)
      else
         do n = 2, size(mnpc)
            X0(1) = X0(1) + txc(mnpc(n))
            X0(2) = X0(2) + tyc(mnpc(n))
            X0(3) = X0(3) + tzc(mnpc(n))
         end do
         X0(4) = real(nM,dp)
      end if
      X0(1:3) = X0(1:3) / X0(4)
      do n = 1, size(mnpc)
         dX(1,n) = txc(mnpc(n)) - X0(1)
         dX(2,n) = tyc(mnpc(n)) - X0(2)
         dX(3,n) = tzc(mnpc(n)) - X0(3)
      end do
      do n = 1, 3
         tolX(n) = max(maxval(-dX(n,1:nM))*epsX, &
              &        maxval( dX(n,1:nM))*epsX, epsDiv0_p)
      end do
      if (ipsw > 3) then
         if (present(W)) then
            write(lpu,699) tolX, X0, dX(:,1), (n,dX(:,1+n),W(n),n=1,nM)
         else
            write(lpu,699) tolX, X0, dX(:,1), (n,dX(:,1+n),1.0_dp,n=1,nM)
         end if
699      format(/25X,'X-coor',7X,'Y-coor',7X,'Z-coor',7X,'Weight', &
              & / 5X,'Tolerance       ',1P3E13.5, &
              & / 5X,'Centroid        ',1P4E13.5, &
              & / 5X,'Reference node  ',1P3E13.5, &
              & /(5X,'Element node', I4,1P4E13.5))
      end if
    end subroutine computePointCoords

    !> @brief Returns the weighted square sum of two arrays.
    !> @details The following sum is returned: @code
    !> Sum_i=1,nM {W(i)*(dX(i)^2+dY(i)^2)}
    !> @endcode
    function sumSqrW (dX,dY,W)
      real(dp),intent(in) :: dX(:), dY(:), W(:)
      real(dp)            :: sumSqrW
      sumSqrW = sum(W*(dX*dX+dY*dY))
    end function sumSqrW

    !> @brief Returns the square sum of two arrays.
    !> @details The following sum is returned: @code
    !> Sum_i=1,nM {dX(i)^2+dY(i)^2}
    !> @endcode
    function sumSqr (dX,dY)
      real(dp),intent(in) :: dX(:), dY(:)
      real(dp)            :: sumSqr
      sumSqr = sum(dX*dX+dY*dY)
    end function sumSqr

  end subroutine wavgmConstrEqn


  !!============================================================================
  !! TODO,kmo: This subroutine is obsolete, retained only for historical reason
  !! and for comparison with the new subroutine wavgmConstrEqn above.
  !> @cond NO_DOCUMENTATION
  !> @brief Computes constraint equation coefficients for a WAVGM element.
  !> @author Knut Morten Okstad
  !> @date 3 Sep 2002

  subroutine wavgmConstrEqn_old (lDof,nM,indC,tolX,dX,weight,omega,ipsw,lpu)

    use kindModule            , only : dp, epsDiv0_p
    use manipMatrixModule     , only : writeObject
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer , intent(in)  :: lDof, nM, indC(:), ipsw, lpu
    real(dp), intent(in)  :: tolX(3), dX(:,:), weight(:)
    real(dp), intent(out) :: omega(:)

    !! Local variables
    integer           :: i, j, k, nndof
    real(dp)          :: sumW, coeffs(nM)
    character(len=32) :: label

    !! --- Logic section ---

    omega = 0.0_dp
    nndof = size(indC)
    if (lDof < 1 .or. lDof > nndof) return

    !!                        k
    !! Cyclic permutation:   / \
    !!                      i - j
    i = 1 + mod(lDof-1,3)
    j = 1 + mod(i,3)
    k = 1 + mod(j,3)

    if (indC(lDof) > 0) then

       !! Direct coupling of the lDof'th dependent DOF to corresponding
       !! independent DOFs
       call DAXPY (nM,1.0_dp,weight(indC(lDof)),1,omega(lDof),nndof)

    end if
    if (lDof > 3 .and. ffa_cmdlinearg_isTrue('useRotCpl')) then

       !! Coupling of rotational dependent DOF to translational DOFs j,k

       if (indC(lDof) > 0) then
          sumW = sum(weight(indC(lDof):indC(lDof)+nM-1))
       else
          sumW = 0.0_dp
       end if

       if (indC(j) > 0) then
          call getRotationCoeff (tolX(k),dX(k,2:), &
               &                 weight(indC(j):indC(j)+nM-1),coeffs,sumW)
          if (ipsw > 5) write(lpu,600) lDof,j,coeffs
          call DAXPY (nM,-1.0_dp,coeffs(1),1,omega(j),nndof)
       end if
       if (indC(k) > 0) then
          call getRotationCoeff (tolX(j),dX(j,2:), &
               &                 weight(indC(k):indC(k)+nM-1),coeffs,sumW)
          if (ipsw > 5) write(lpu,600) lDof,k,coeffs
          call DAXPY (nM,1.0_dp,coeffs(1),1,omega(k),nndof)
       end if
       if (ipsw > 5) write(lpu,610) lDof,sumW
       if (sumW > epsDiv0_p) call DSCAL (nndof*nM,1.0_dp/sumW,omega(1),1)

    else if (lDof < 4 .and. ffa_cmdlinearg_isTrue('useEccCpl')) then
       if (abs(dX(j,1)) > tolX(j)) then

          !! Coupling of translational dependent DOF to translational DOFs i,j
          !! and rotational DOF 3+k through the eccentricity e_j
          !! of the dependent node w.r.t. the independent node centroid

          if (3+k <= nndof .and. indC(3+k) > 0) then
             sumW = sum(weight(indC(3+k):indC(3+k)+nM-1))
          else
             sumW = 0.0_dp
          end if
          if (indC(i) > 0) sumW = sumW + getRotationWeight(tolX(j),dX(j,2:), &
               &                         weight(indC(i):indC(i)+nM-1))
          if (indC(j) > 0) sumW = sumW + getRotationWeight(tolX(i),dX(i,2:), &
               &                         weight(indC(j):indC(j)+nM-1))
          if (ipsw > 5) write(lpu,610) lDof,sumW
          if (sumW > epsDiv0_p) then
             sumW = dX(j,1)/sumW
          else
             sumW = 0.0_dp
          end if

          if (3+k <= nndof .and. indC(3+k) > 0) then
             call DAXPY (nM,-sumW,weight(indC(3+k)),1,omega(3+k),nndof)
          end if
          if (indC(i) > 0) then
             call getRotationCoeff (tolX(j),dX(j,2:), &
                  &                 weight(indC(i):indC(i)+nM-1),coeffs)
             if (ipsw > 5) write(lpu,600) lDof,i,coeffs*dX(j,1)
             call DAXPY (nM,sumW,coeffs(1),1,omega(i),nndof)
          end if
          if (indC(j) > 0) then
             call getRotationCoeff (tolX(i),dX(i,2:), &
                  &                 weight(indC(j):indC(j)+nM-1),coeffs)
             if (ipsw > 5) write(lpu,600) lDof,j,coeffs*dX(j,1)
             call DAXPY (nM,-sumW,coeffs(1),1,omega(j),nndof)
          end if

       end if
       if (abs(dX(k,1)) > tolX(k)) then

          !! Coupling of translational dependent DOF to translational DOFs i,k
          !! and rotational DOF 3+j through the eccentricity e_k
          !! of the dependent node w.r.t. the independent node centroid

          if (3+j <= nndof .and. indC(3+j) > 0) then
             sumW = sum(weight(indC(3+j):indC(3+j)+nM-1))
          else
             sumW = 0.0_dp
          end if
          if (indC(i) > 0) sumW = sumW + getRotationWeight(tolX(k),dX(k,2:), &
               &                         weight(indC(i):indC(i)+nM-1))
          if (indC(k) > 0) sumW = sumW + getRotationWeight(tolX(i),dX(i,2:), &
               &                         weight(indC(k):indC(k)+nM-1))
          if (ipsw > 5) write(lpu,610) lDof,sumW
          if (sumW > epsDiv0_p) then
             sumW = dX(k,1)/sumW
          else
             sumW = 0.0_dp
          end if

          if (3+j <= nndof .and. indC(3+j) > 0) then
             call DAXPY (nM,sumW,weight(indC(3+j)),1,omega(3+j),nndof)
          end if
          if (indC(i) > 0) then
             call getRotationCoeff (tolX(k),dX(k,2:), &
                  &                 weight(indC(i):indC(i)+nM-1),coeffs)
             if (ipsw > 5) write(lpu,600) lDof,i,coeffs*sumW
             call DAXPY (nM,sumW,coeffs(1),1,omega(i),nndof)
          end if
          if (indC(k) > 0) then
             call getRotationCoeff (tolX(i),dX(i,2:), &
                  &                 weight(indC(k):indC(k)+nM-1),coeffs)
             if (ipsw > 5) write(lpu,600) lDof,k,coeffs*sumW
             call DAXPY (nM,-sumW,coeffs(1),1,omega(k),nndof)
          end if

       end if
    end if

    if (ipsw > 4) then
       write(lpu,*)
       write(label,"('     Omega for dependent DOF',I3,' :')") lDof
       call writeObject(omega,lpu,label,nndof)
    end if

600 format(/5X,'Coefficients coupling local dependent DOF',I3, &
         &     ' to local independent DOFs',I3,' :'/ (4X,1P,6E13.5) )
610 format(/5X,'Sum weights for local dependent DOF',I3,' :',1PE12.5 )

  contains

    !> @brief Computes weighting coefficients for a rotational dependent DOF.
    subroutine getRotationCoeff (tolX,dX,weight,coeff,sumW)
      real(dp), intent(in)  :: tolX, dX(:), weight(:)
      real(dp), intent(out) :: coeff(:)
      real(dp), intent(inout), optional :: sumW
      integer :: i, j
      do i = 1, nM
         if (abs(dX(i)) > tolX) then
            coeff(i) = 1.0_dp/dX(i)
            if (present(sumW)) sumW = sumW + weight(i)
         else
            coeff(i) = 0.0_dp
         end if
         do j = 1, nM
            if (abs(dX(j)) > tolX) coeff(i) = coeff(i) - weight(j)/dX(j)
         end do
         coeff(i) = coeff(i)*weight(i)
      end do
    end subroutine getRotationCoeff

    !> @brief Returns the sum of weights for a rotational dependent DOF.
    function getRotationWeight (tolX,dX,weight) result(sumW)
      real(dp), intent(in) :: tolX, dX(:), weight(:)
      real(dp) :: sumW
      integer  :: i
      sumW = 0.0_dp
      do i = 1, size(weight)
         if (abs(dX(i)) > tolX) sumW = sumW + weight(i)
      end do
    end function getRotationWeight

  end subroutine wavgmConstrEqn_old
  !> @endcond

end module WAvgMotionModule
