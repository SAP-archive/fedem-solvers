!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ModesRoutinesModule2

  implicit none

contains

  subroutine readFrequencies (nModes,modes,processThisMode,freq,ierr)

    !!==========================================================================
    !! Retrieve result database file pointers for system modal data.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 15 Nov 2000/1.0
    !!==========================================================================

    use KindModule           , only : dp
    use ReportErrorModule    , only : reportError, error_p, debugFileOnly_p
    use FFrExtractorInterface, only : ffr_findData

    integer , intent(in)  :: nModes, modes(:)
    logical , intent(in)  :: processThisMode(:)
    real(dp), intent(out) :: freq(:)
    integer , intent(out) :: ierr

    !! Local variables
    integer, parameter :: mechId_p = 2
    integer            :: j, lerr
    character(len=64)  :: chPath

    !! --- Logic section ---

    ierr = 0
    do j = 1, nModes
       if (.not. processThisMode(j)) cycle

       write(chPath,"('Eigenvalues|Mode',i3,'|Eigenfrequency')") modes(j)
       call ffr_findData (freq(j),1,trim(chPath),'Mechanism',mechId_p,lerr)
       if (lerr < 0) then
          ierr = ierr + 1
          call reportError (error_p,'Can not find eigenfrequency for '// &
               &            chPath(13:29))
          cycle
       end if
    end do

    if (ierr > 0) call reportError (debugFileOnly_p,'readFrequencies')

  end subroutine readFrequencies


  subroutine readModesPointers (sup,nModes,modes,processThisMode,pointers,ierr)

    !!==========================================================================
    !! Retrieve result database file pointers for system modal data.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 15 Nov 2000/1.0
    !!==========================================================================

    use SupElTypeModule      , only : SupElType
    use IdTypeModule         , only : getId
    use PointerKindModule    , only : ptr
    use ReportErrorModule    , only : reportError, error_p, debugFileOnly_p
    use FFrExtractorInterface, only : ffr_findPtr

    type(SupElType), intent(in)    :: sup
    integer        , intent(in)    :: nModes, modes(:)
    logical        , intent(inout) :: processThisMode(:,:)
    integer(ptr)   , intent(out)   :: pointers(:,:)
    integer        , intent(out)   :: ierr

    !! Local variables
    integer           :: i, j
    character(len=64) :: chPath

    !! --- Logic section ---

    ierr = 0
    do j = 1, nModes
       write(chPath,"('Eigenvectors|Mode',i3)") modes(j)
       do i = 1, sup%nExtNods
          call ffr_findPtr (trim(chPath),'Triad',sup%triads(i)%p%id%baseId, &
               &            pointers(i,j))
          if (pointers(i,j) == 0_ptr) then
             ierr = ierr + 1
             call reportError (error_p,'Can not find eigenvector components'// &
                  &            ' for Triad'//trim(getId(sup%triads(i)%p%id))// &
                  &            ' and '//chPath(14:20))
             pointers(:,j) = 0_ptr
             processThisMode(:,j) = .false. ! Skip this mode and continue
             exit
          end if
       end do
       i = sup%nExtNods+1
       if (sup%genDOFs%nDOFs > 0) then
          call ffr_findPtr (trim(chPath),'Part',sup%id%baseId,pointers(i,j))
          if (pointers(i,j) == 0_ptr .and. any(processThisMode(:,j))) then
             ierr = ierr + 1
             call reportError (error_p,'Can not find eigenvector components'// &
                  &            ' for Part'//trim(getId(sup%id))// &
                  &            ' (component modes) and '//chPath(14:20))
             pointers(:,j) = 0_ptr
             processThisMode(:,j) = .false. ! Skip this mode and continue
          end if
       else
          pointers(i,j) = 0_ptr
       end if
    end do

    if (ierr > 0) call reportError (debugFileOnly_p,'readModesPointers')

  end subroutine readModesPointers


  subroutine readSupElModes (iComp,nModes,processThisMode,pointers,sup, &
       &                     eigVec,eigFinit,ierr)

    !!==========================================================================
    !! Read the system eigenvectors of the modal analysis from the results
    !! database and calculate the associated local deformational displacements
    !! for the specified superelement.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 15 Nov 2000/1.0
    !!==========================================================================

    use SupElTypeModule      , only : SupElType, dp
    use IdTypeModule         , only : getId
    use PointerKindModule    , only : ptr
    use ManipMatrixModule    , only : invert34
    use ReportErrorModule    , only : reportError, error_p, debugFileOnly_p
    use FFrExtractorInterface, only : ffr_getData

    integer        , intent(in)    :: iComp, nModes
    logical        , intent(in)    :: processThisMode(:)
    integer(ptr)   , intent(in)    :: pointers(:,:)
    type(SupElType), intent(in)    :: sup
    real(dp)       , intent(inout) :: eigVec(:)
    real(dp)       , intent(out)   :: eigFinit(sup%nTotDofs,iComp,nModes)
    integer        , intent(out)   :: ierr

    !! Local variables
    integer            :: i, j, k, l, m, n, lerr
    real(dp)           :: tInv(3,4)
    character(len=128) :: errMsg

    !! --- Logic section ---

    ierr = 0
    tInv = invert34(sup%supTr)
    do j = 1, nModes
       if (.not. processThisMode(j)) cycle

       do i = 1, sup%nExtNods

          !! Read triad eigenvector components of this mode from database file
          n = sup%triads(i)%p%nDOFs
          call ffr_getData (eigVec(1),n*iComp,pointers(i,j),lerr)
          if (lerr /= 0) then
             ierr = ierr + 1
             write(errMsg,600) n*iComp,lerr+n*iComp
             call reportError (error_p,'Can not read eigenvector components'// &
                  &            ' for Triad'//getId(sup%triads(i)%p%id), &
                  &            addString=errMsg)
             cycle
          end if

          !! Find the corresponding superelement deformation vector
          k = sup%triads(i)%firstDOF
          m = 1
          do l = 1, iComp
             m = n*(l-1) + 1
             if (n >= 3) then
                eigFinit(k:k+2,l,j)  = matmul(tInv(:,1:3),eigVec(m:m+2))
             end if
             if (n >= 6) then
                eigFinit(k+3:k+5,l,j)= matmul(tInv(:,1:3),eigVec(m+3:m+5))
             end if
             m = m + n
          end do

       end do
       if (sup%genDOFs%nDOFs > 0) then

          !! Read and add generalized dofs
          i = sup%nExtNods+1
          n = sup%genDOFs%nDOFs
          call ffr_getData (eigVec(1),n*iComp,pointers(i,j),lerr)
          if (lerr /= 0) then
             ierr = ierr + 1
             write(errMsg,600) n*iComp,lerr+n*iComp
             call reportError (error_p,'Can not read eigenvector components '//&
                  &            'for Part'//getId(sup%id)//' (component modes)',&
                  &            addString=errMsg)
             cycle
          end if

          k = sup%genDOFs%firstDOF
          do l = 1, iComp
             eigFinit(k:k+n-1,l,j) = eigVec(n*(l-1)+1:n*l)
          end do

       end if
    end do

    if (ierr > 0) call reportError (debugFileOnly_p,'readSupElModes')

600 format('Mismatch between length of wanted array',i3, &
         & ' and actual variable size',i3,' on the results file')

  end subroutine readSupElModes


  subroutine calcStrainEnergyDensity (rdb,samSup,sv,iprint,lpu,err)

    !!==========================================================================
    !! Calculate and write scaled strain energy density to the result database.
    !!
    !! Programmer : Tommy Jorstad
    !! date/rev   : 14 May 2002/1.0
    !!==========================================================================

    use KindModule            , only : dp, hugeVal_p
    use RDBModule             , only : RDBType, writeRDB
    use SamModule             , only : SamType
    use ElStressModule        , only : extractEV, elStress
    use ReportErrorModule     , only : reportError, debugFileOnly_p
    use ReportErrorModule     , only : internalError
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(RDBType), intent(inout) :: rdb
    type(SamType), intent(in)    :: samSup
    real(dp)     , intent(in)    :: sv(:)
    integer      , intent(in)    :: lpu, iprint
    integer      , intent(out)   :: err

    !! Local variables
    integer, parameter :: maxstp = 100, maxnod = 8

    integer  :: iel, lerr, nenod, nstrp, n
    real(dp) :: SR(6*maxnod),SS(6*maxnod),Stress(6*maxstp),Strain(6*maxstp)
    real(dp) :: v(6*maxstp),resVec(maxstp)
    logical  :: lDouble

    !! --- Logic section ---

    call ffa_cmdlinearg_getbool ('double',lDouble)

    do iel = 1, samSup%nel

       !! Extract element nodal displacements
       call extractEV (iel,samSup,sv,v,err)
       if (err < 0) return

       !! Compute element stress/strain results
       call elStress (iel,samSup%melcon(iel),v,SR,SS,Stress,Strain, &
            &         nenod=nenod,nstrp=nstrp,ipsw=iprint,lpu=lpu,ierr=lerr)
       if (lerr < 0) then
          call reportError (debugFileOnly_p,'calcStrainEnergyDensity')
          err = lerr
          return
       else if (nstrp > maxstp) then
          err = internalError('calcStrainEnergyDensity: Too many stress points')
          return
       else if (nenod > maxnod) then
          err = internalError('calcStrainEnergyDensity: Too many element nodes')
          return
       end if

       if (lerr > 0) then ! Stress calculation failed, write hugeVal
          resVec = hugeVal_p
       else
          !! Compute scaled strain energy density
          v = Stress*Strain
          do n = 1, nstrp
             select case (samSup%melcon(iel))
             case (21:24)
                resVec(n) = sum(v(3*n-2:3*n)) + v(3*n)
             case (31:32,41:46)
                resVec(n) = sum(v(6*n-5:6*n)) + sum(v(6*n-2:6*n))
             end select
          end do
          if (iprint > 0) write(LPU,600) resVec(1:nstrp)
       end if

       !! Write the scaled strain energy density to results database
       call writeRDB (rdb,resVec(1:nstrp),err,lDouble)
       if (err < 0) then
          call reportError (debugFileOnly_p,'calcStrainEnergyDensity')
          return
       end if

    end do

600 format(/5X,'Scaled strain energy density at stress points:', &
         & /1P,(11X,6E13.5))

  end subroutine calcStrainEnergyDensity

end module ModesRoutinesModule2
