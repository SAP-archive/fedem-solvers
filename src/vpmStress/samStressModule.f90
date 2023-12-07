!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file samStressModule.f90
!> @brief Initialisation of the SAM data structure for recovery processes.

!!==============================================================================
!> @brief Initialisation of the SAM data structure for recovery processes.

module SamStressModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes the SAM data structure for the recovery processes.
  !>
  !> @param[in] filnam Name of binary file to read the SAM data from
  !> @param[out] sam Assembly management data for the superelement
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine basically reads the contents of the arrays in
  !> sammodule::samtype from the fsm-file that was produced by fedem_reducer.
  !> It also modifies a few of the arrays after reading from file,
  !> to fit special usage of the recovery processes.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Sep 2000

  subroutine initiateSAM (filnam,sam,lpu,ierr)

    use SamModule              , only : SamType, nullifySam, writeObject
    use SamModule              , only : invertEqPartition, check
    use allocationModule       , only : reAllocate
    use reportErrorModule      , only : reportError
    use reportErrorModule      , only : note_p, error_p, debugFileOnly_p
    use binaryDBInterface      , only : openBinaryDB, closeBinaryDB, read_p
    use binaryDBInterface      , only : readTagDB, readIntDB, readDoubleDB
    use FFlLinkHandlerInterface, only : ffl_getSize, ffl_getElmId

    character(len=*), intent(in)  :: filnam
    type(SamType)   , intent(out) :: sam
    integer         , intent(in)  :: lpu
    integer         , intent(out) :: ierr

    !! Local variables
    character(len=32) :: ctag
    integer :: iel, inod, ip, ifile, npar, nael, nanod, &
         &     nnod, nel, ndof, nmmnpc, nxnod, ncons

    !! --- Logic section ---

    call nullifySam (sam)
    call openBinaryDB (filnam,read_p,ifile,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Error opening SAM-file '//filnam)
       goto 900
    end if

    call readTagDB (ifile,ctag,npar,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Error reading file tag from '//filnam)
       goto 900
    else if (ctag /= '#SAM data') then
       call reportError (error_p,'File "'//trim(filnam)// &
            &            '" is not a SAM data file, tag='//trim(ctag))
       ierr = 1
       goto 900
    else
       write(lpu,"()")
       call reportError (note_p,'Reading FE data from: '//filnam)
    end if

    !! Initialize the Matrix of Parameters
    call readIntDB (ifile,npar,1,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Error reading MPAR size from file '//filnam)
       goto 900
    end if

    ierr = 0
    call reAllocate ('initiateSAM',sam%mpar,npar,ierr)
    if (ierr < 0) return

    call readIntDB (ifile,sam%mpar(1),npar,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Error reading MPAR from file '//filnam)
       goto 900
    end if

    sam%nnod   => sam%mpar(1)  !! N. of active nodes (nanod)
    sam%nel    => sam%mpar(2)  !! N. of elements
    sam%ndof   => sam%mpar(3)  !! N. of degrees of freedom
    sam%ndof1  => sam%mpar(4)  !! N. of 1-status dofs
    sam%ndof2  => sam%mpar(5)  !! N. of 2-status dofs
    sam%nspdof => sam%mpar(6)  !! N. of specified dofs
    sam%nceq   => sam%mpar(7)  !! N. of constraint equations
    sam%nsdof  => sam%mpar(8)  !! N. of suppressed dofs
    sam%npdof  => sam%mpar(9)  !! N. of prescribed dofs
    sam%nddof  => sam%mpar(10) !! N. of dependent dofs
    sam%neq    => sam%mpar(11) !! N. of equations
    sam%nsky   => sam%mpar(12) !! N. of elements in skyline array
    sam%nmmnpc => sam%mpar(15) !! N. of elements in mmnpc
    sam%nmmceq => sam%mpar(16) !! N. of elements in mmceq

    !! Get some FE model size parameters
    call ffl_getSize (nnod, nel, ndof, nmmnpc, ip, nxnod, &
         &            iel, inod, nael, nanod, & ! dummy arguments
         &            npar, ncons, ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failed to obtain FE model size parameters')
       goto 900
    end if

    nael = ierr
    ierr = 0 ! Check consistency with the SAM data
    if (nnod   /= sam%nnod)     ierr = ierr - 1
    if (nel    /= sam%nel)      ierr = ierr - 1
    if (ndof   /= sam%ndof)     ierr = ierr - 1
    if (nmmnpc  < sam%nmmnpc)   ierr = ierr - 1 ! Allow smaller size in SAM
    if (nxnod  /= sam%mpar(23)) ierr = ierr - 1
    if (ncons  /= sam%mpar(28)) ierr = ierr - 1
    if (ierr < 0) then
       call reportError (error_p,'The specified FE data file is not '// &
            'compatible with the internal SAM data.','Check that '//  &
            'the reduced FE data is consistent with the current model.')
       write(lpu,600) nnod, sam%nnod, nel, sam%nel, ndof, sam%ndof, &
            &         nmmnpc, sam%nmmnpc, nxnod, sam%mpar(23), &
            &         ncons, sam%mpar(28)
600    format(/12X,'FE part   SAM data' &
            & / 2X,'nnod  ', 2I11 &
            & / 2X,'nel   ', 2I11 &
            & / 2X,'ndof  ', 2I11 &
            & / 2X,'nmmnpc', 2I11 &
            & / 2X,'nxnod ', 2I11 &
            & / 2X,'ncons ', 2I11 /)
       goto 900
    end if

    !! Allocate and read in all SAM-arrays from file

    call reAllocate ('initiateSAM',sam%madof,sam%nnod+1,ierr)
    call reAllocate ('initiateSAM',sam%minex,sam%nnod,ierr)
    call reAllocate ('initiateSAM',sam%mnnn,sam%nnod,ierr)
    call reAllocate ('initiateSAM',sam%msc,sam%ndof,ierr)
    call reAllocate ('initiateSAM',sam%mpmnpc,sam%nel+1,ierr)
    call reAllocate ('initiateSAM',sam%mmnpc,sam%nmmnpc,ierr)
    call reAllocate ('initiateSAM',sam%melcon,sam%nel,ierr)
    call reAllocate ('initiateSAM',sam%mpmceq,sam%nceq+1,ierr)
    call reAllocate ('initiateSAM',sam%mmceq,sam%nmmceq,ierr)
    call reAllocate ('initiateSAM',sam%ttcc,sam%nmmceq,ierr)
    call reAllocate ('initiateSAM',sam%meqn,sam%ndof,ierr)
    call reAllocate ('initiateSAM',sam%meqn1,sam%ndof1,ierr)
    call reAllocate ('initiateSAM',sam%meqn2,sam%ndof2,ierr)
    if (ierr < 0) then
       call writeObject (sam,lpu)
       goto 900
    end if

    call readSAMarrays (ifile,sam,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Error reading SAM-arrays from file '//filnam)
       goto 900
    end if

    call closeBinaryDB (ifile,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Error closing SAM-file '//filnam)
       goto 900
    end if

    if (nael < nel) then

       !! Flag all active internal nodes by a negative value in sam%mnnn.
       !! This can then be used to restrict the deformation recovery to
       !! this group of nodes only, for efficiency.
       do iel = 1, nel
          if (ffl_getElmId(iel) > 0) then
             do ip = sam%mpmnpc(iel), sam%mpmnpc(iel+1)-1
                inod = sam%mmnpc(ip)
                if (sam%mnnn(inod) == 1) then
                   sam%mnnn(inod) = -1
                else if (sam%mnnn(inod) == 2) then
                   sam%mnnn(inod) = 0 ! Active external node
                   write(lpu,601) sam%minex(inod),ffl_getElmId(iel), &
                       &          ip-sam%mpmnpc(iel)+1
601                format(10X,'Group node',I8,' in element',I8, &
                       &  ' (local node',I3,') is marked as external')
                end if
             end do
          end if
       end do

       !! Second pass, for contraint elements with one of their dependent nodes
       !! marked, the associated independent nodes need to be marked as well
       do iel = 1, nel
          if (sam%melcon(iel) >= 61 .and. sam%melcon(iel) <= 63) then
             do ip = sam%mpmnpc(iel), sam%mpmnpc(iel+1)-1
                inod = sam%mmnpc(ip)
                if (sam%mnnn(inod) < 1) goto 100
             end do
             cycle
100          do ip = sam%mpmnpc(iel), sam%mpmnpc(iel+1)-1
                inod = sam%mmnpc(ip)
                if (sam%mnnn(inod) > 0) sam%mnnn(inod) = sam%mnnn(inod) - 2
             end do
          end if
       end do

       !! Now calculate running indices to the active internal group nodes
       nanod = 0
       do inod = 1, nnod
          if (sam%mnnn(inod) < 0) then
             nanod = nanod + 1
             sam%mnnn(inod) = -nanod
          end if
       end do

       sam%mpar(32) = nael
       sam%mpar(33) = nanod
       write(lpu,"()")
       write(lpu,610) nael,'elements',nel
       write(lpu,610) nanod,'internal nodes',nnod
610    format('Note :',I8,' active ',A,' in the specified group(s) ', &
           &  '(out of total',I8,').')

    end if
    if (sam%ndof2 > 0) then

       !! Form the control arrays dofPosIn1 and dofPosIn2
       call invertEqPartition (nael < nel, .true., sam, ierr)
       if (ierr < 0) goto 900

    end if

    !! Modify sam%msc such that it contains a running index of the external DOFs
    !! whereas the free internal DOFs will have status code zero and all fixed
    !! DOFs will have a negative value
    ip = 0
    do inod = 1, ndof
       if (sam%msc(inod) == 2) then
          ip = ip + 1
          sam%msc(inod) = ip
       else
          sam%msc(inod) = sam%msc(inod) - 1
       end if
    end do

    !! The system matrices are not used in fedem_stress
    call check (sam,lpu,ierr,.false.)
    if (ierr == 0) return

    call writeObject (sam,lpu,3)
900 call reportError (debugFileOnly_p,'initiateSAM')

  contains

    !> @brief Reads all SAM arrays from binary file.
    subroutine readSAMarrays (ifile,sam,ierr)
      integer      , intent(in)    :: ifile
      type(SamType), intent(inout) :: sam
      integer      , intent(out)   :: ierr

      call readIntDB (ifile,sam%madof(1),size(sam%madof),ierr)
      if (ierr < 0) return

      call readIntDB (ifile,sam%minex(1),size(sam%minex),ierr)
      if (ierr < 0) return

      call readIntDB (ifile,sam%mnnn(1),size(sam%mnnn),ierr)
      if (ierr < 0) return

      call readIntDB (ifile,sam%msc(1),size(sam%msc),ierr)
      if (ierr < 0) return

      call readIntDB (ifile,sam%mpmnpc(1),size(sam%mpmnpc),ierr)
      if (ierr < 0) return

      call readIntDB (ifile,sam%mmnpc(1),size(sam%mmnpc),ierr)
      if (ierr < 0) return

      call readIntDB (ifile,sam%melcon(1),size(sam%melcon),ierr)
      if (ierr < 0) return

      if (sam%nceq > 0) then
         call readIntDB (ifile,sam%mpmceq(1),size(sam%mpmceq),ierr)
         if (ierr < 0) return

         call readIntDB (ifile,sam%mmceq(1),size(sam%mmceq),ierr)
         if (ierr < 0) return

         call readDoubleDB (ifile,sam%ttcc(1),size(sam%ttcc),ierr)
         if (ierr < 0) return
      else
         sam%mpmceq(1) = 1
      end if

      call readIntDB (ifile,sam%meqn(1),size(sam%meqn),ierr)
      if (ierr < 0) return

      if (sam%ndof1 > 0) then
         call readIntDB (ifile,sam%meqn1(1),size(sam%meqn1),ierr)
         if (ierr < 0) return
      end if
      if (sam%ndof2 > 0) then
         call readIntDB (ifile,sam%meqn2(1),size(sam%meqn2),ierr)
      end if

    end subroutine readSAMarrays

  end subroutine initiateSAM

end module SamStressModule
