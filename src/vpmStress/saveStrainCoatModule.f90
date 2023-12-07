!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module saveStrainCoatModule

  use StrainCoatModule, only : maxCoats_p

  implicit none

  private

  integer, save :: id(15) = 0
  integer, save :: idSTRC(4,0:maxCoats_p) = 0

  public :: initiateStrainCoatHeader, finalizeStrainCoatHeader
  public :: writeElementsHeader, writeElementsDB
  public :: writeHistoryHeader, writeHistoryDB


contains

  subroutine writeHistoryHeader (resFileName,modelFileName,rdb,sup, &
       &                         strainCoats,ierr)

    !!==========================================================================
    !! Administer writing of history results database header for fedem_fpp.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 10 Dec 2003/1.0
    !!==========================================================================

    use RDBModule        , only : RDBType
    use SupElTypeModule  , only : SupElType
    use StrainCoatModule , only : StrainCoatType
    use reportErrorModule, only : reportError, debugFileOnly_p

    character(len=*)    , intent(in)    :: resFileName, modelFileName
    type(RDBType)       , intent(inout) :: rdb
    type(SupElType)     , intent(in)    :: sup
    type(StrainCoatType), intent(in)    :: strainCoats(:)
    integer             , intent(out)   :: ierr

    !! --- Logic section ---

    call initiateStrainCoatHeader (modelFileName,rdb,sup,ierr)
    if (ierr < 0) goto 500

    call writeElementsHeader (rdb,strainCoats,.true.,ierr)
    if (ierr < 0) goto 500

    call finalizeStrainCoatHeader (resFileName,rdb,ierr)
    if (ierr == 0) return

500 call reportError (debugFileOnly_p,'writeHistoryHeader')

  end subroutine writeHistoryHeader


  subroutine initiateStrainCoatHeader (modelFileName,rdb,sup,ierr)

    !!==========================================================================
    !! Administer writing of results database header for fedem_fpp.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use kindModule            , only : lfnam_p
    use RDBModule             , only : RDBType, nullifyRDB, idatd
    use RDBModule             , only : openHeaderFiles, writeTimeStepHeader
    use idTypeModule          , only : writeIdHeader
    use SupElTypeModule       , only : SupElType
    use StrainCoatModule      , only : StrainCoatType
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring

    character(len=*), intent(in)  :: modelFileName
    type(RDBType)   , intent(out) :: rdb
    type(SupElType) , intent(in)  :: sup
    integer         , intent(out) :: ierr

    !! Local variables
    character(lfnam_p) :: linkFileName

    !! --- Logic section ---

    call nullifyRDB (rdb)

    id = 0
    idSTRC = 0

    call ffa_cmdlinearg_getstring ('linkfile',linkFileName)
    call openHeaderFiles ('fedem_fpp',ierr, &
         &                'strain coat data base file', &
         &                modelFileName,linkFileName)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'initiateStrainCoatHeader')
       return
    end if

    call writeTimeStepHeader (rdb)
    call writeIdHeader ('Part',sup%id,idatd,.true.)
    write(idatd,600)

    rdb%nBstep = int(rdb%nBytes)
    rdb%nBytes = 0

600 format('  [;"Elements";')

  end subroutine initiateStrainCoatHeader


  subroutine finalizeStrainCoatHeader (chname,rdb,ierr)

    !!==========================================================================
    !! Terminate writing of results database header and move it to binary file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use RDBModule        , only : RDBType, openRDBfile, idatd
    use reportErrorModule, only : reportError, debugFileOnly_p

    character(len=*), intent(in)    :: chname
    type(RDBType)   , intent(inout) :: rdb
    integer         , intent(out)   :: ierr

    !! --- Logic section ---

    write(idatd,"('  ]'/'}')")
    call openRDBfile (rdb,ierr,chname,'#FEDEM strain coat data')
    if (ierr < 0) call reportError (debugFileOnly_p,'finalizeStrainCoatHeader')

  end subroutine finalizeStrainCoatHeader


  subroutine writeElementsHeader (rdb,strainCoats,writeHist,ierr)

    !!==========================================================================
    !! Write element result definitions to the temporary header files.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use RDBModule             , only : RDBType, openRDBfile
    use StrainCoatModule      , only : StrainCoatType
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(RDBType)       , intent(inout) :: rdb
    type(StrainCoatType), intent(in)    :: strainCoats(:)
    logical             , intent(in)    :: writeHist
    integer             , intent(out)   :: ierr

    !! Local variables
    logical :: lDouble
    integer :: i, nBits, nenod

    !! --- Logic section ---

    call ffa_cmdlinearg_getbool ('double',lDouble)
    if (lDouble) then
       nBits = 64 ! Write all results in double precision
    else
       nbits = 32 ! Write all results in single precision
    end if

    do i = 1, size(strainCoats)

       nenod = size(strainCoats(i)%tensorRosette%globalNodes)
       call writeStrainCoatHeader (rdb, strainCoats(i)%results, &
            &                      writeHist, strainCoats(i)%fppType, &
            &                      strainCoats(i)%tensorRosette%id%userId, &
            &                      nenod, nBits, ierr)
       if (ierr < 0) then
          call reportError (debugFileOnly_p,'writeElementsHeader')
          return
       end if

    end do

  end subroutine writeElementsHeader


  subroutine writeStrainCoatHeader (rdb,coats,writeHist,fppType, &
       &                            elmno,nenod,nbit,ierr)

    !!==========================================================================
    !! Write strain coat item group definition to the temporary header files.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use RDBModule        , only : RDBType, idatd, ivard, iitem
    use StrainCoatModule , only : CoatResultType
    use StrainCoatModule , only : maxNodes_p, elmtype_p, resultSet_p
    use reportErrorModule, only : internalError

    type(RDBType)       , intent(inout) :: rdb
    type(CoatResultType), intent(in)    :: coats(:)
    logical             , intent(in)    :: writeHist
    integer             , intent(in)    :: fppType, elmno, nenod, nbit
    integer             , intent(out)   :: ierr

    !! Local variables
    integer           :: i, iRset, itemG, jVar, nCoats
    character(len=16) :: chVars

    !! --- Logic section ---

    ierr = 0
    nCoats = size(coats)
    if (nCoats > maxCoats_p) then
       ierr = internalError('writeStrainCoatHeader: Too many strain coats')
       return
    else if (nCoats < 1) then
       return ! Empty strain coat (switched off with -surface option)
    end if
    if (nenod < 3 .or. nenod > maxNodes_p) then
       ierr = internalError('writeStrainCoatHeader: Invalid element type')
       return
    end if

    jVar = 1
    do i = 1, nCoats

       itemG = 1
       iRset = coats(i)%resultSet
       if (iRset < 0 .or. iRset > maxCoats_p) then
          ierr = internalError('writeStrainCoatHeader: invalid result set')
          return
       end if

       if (writeHist) then
          call writeVarDef (id(1),'Max principal stress','FORCE/AREA')
          call writeVarDef (id(2),'Min principal stress','FORCE/AREA')
          call writeVarDef (id(3),'Signed abs max stress','FORCE/AREA')
          call writeVarDef (id(4),'Max shear stress','FORCE/AREA')
          call writeVarDef (id(5),'Von Mises stress','FORCE/AREA')
          call writeVarDef (id(6),'Max principal strain','NONE')
          call writeVarDef (id(7),'Min principal strain','NONE')
          call writeVarDef (id(8),'Signed abs max strain','NONE')
          call writeVarDef (id(9),'Max shear strain','NONE')
          call writeVarDef (id(10),'Von Mises strain','NONE')
          call writeItemGroup (idSTRC(1,iRset),10)
       else
          call writeVarDef (id(1),'Max principal stress','FORCE/AREA')
          call writeVarDef (id(2),'Max shear stress','FORCE/AREA')
          call writeVarDef (id(3),'Max stress range','FORCE/AREA')
          call writeVarDef (id(4),'Max von Mises stress','FORCE/AREA')
          call writeVarDef (id(5),'Max principal strain','NONE')
          call writeVarDef (id(6),'Max shear strain','NONE')
          call writeVarDef (id(7),'Max strain range','NONE')
          call writeVarDef (id(8),'Max von Mises strain','NONE')
          if (fppType > 0) then
             itemG = 2
             call writeVarDef (id(9),'Damage','NONE')
             call writeVarDef (id(10),'Life (equnits)','TIME')
             call writeVarDef (id(11),'Life (repeats)','NONE')
          else
             id(9)  = -abs(id(9))
             id(10) = -abs(id(10))
             id(11) = -abs(id(11))
             rdb%nBstep = rdb%nBstep - 3*nbit/8
          end if
          call writeVarDef (id(12),'Angle spread','ANGLE')
          call writeVarDef (id(13),'Most popular angle','ANGLE')
          if (coats(i)%nBiAxial > 0) then
             itemG = itemG + 2
             call writeVarDef (id(14),'Mean bi-axiality','NONE')
             call writeVarDef (id(15),'Biaxiality standard deviation','NONE')
             call writeItemGroup (idSTRC(itemG,iRset),15)
          else
             call writeItemGroup (idSTRC(itemG,iRset),13)
          end if
       end if
       if (ierr < 0) return
       write(chVars(jVar:),600) idSTRC(itemG,iRset)
       jVar = jVar + 3
    end do

    if (fppType >= 0) then
       write(idatd,604) elmno,elmtype_p(nenod),trim(chVars)
    else
       !! Currently write STRCT3 always to match nCode's frs-files
       write(idatd,604) elmno,'STRCT3',trim(chVars)
    end if

600 format('[',i1,']')
604 format('    [;',i8,';[;"',a,'";[;"Element";',a,']]]')

  contains

    subroutine writeVarDef (idVar,name,unit)
      integer         , intent(inout) :: idVar
      character(len=*), intent(in)    :: name, unit
      if (idVar > 0 .or. ierr < 0) then
         return
      else if (idVar < 0) then
         idVar = -idVar
      else if (rdb%nVar > 98 .or. nbit > 99) then
         ierr = ierr + internalError('writeStrainCoatHeader: format overflow')
      else
         rdb%nVar = rdb%nVar + 1
         idVar = rdb%nVar
         if (idVar > 9) then
            write(ivard,602) idVar,name,unit,nbit
         else
            write(ivard,601) idVar,name,unit,nbit
         end if
      end if
601   format('<',i1,';"',a,'";',a,';FLOAT;',i2,';SCALAR>')
602   format('<',i2,';"',a,'";',a,';FLOAT;',i2,';SCALAR>')
    end subroutine writeVarDef

    subroutine writeItemGroup (idGr,nVar)
      integer, intent(inout) :: idGr
      integer, intent(in)    :: nVar
      integer                :: iVar
      rdb%nBstep = rdb%nBstep + nVar*nbit/8
      if (idGr > 0 .or. ierr < 0) then
         return
      else if (rdb%nIG > 8 .or. nVar > 99) then
         ierr = ierr + internalError('writeStrainCoatHeader: format overflow')
      else
         rdb%nIG = rdb%nIG + 1
         idGr = rdb%nIG
         write(iitem,600,advance='no') idGr,trim(resultSet_p(iRset))
         do iVar = 1, nVar
            if (id(iVar) > 9) then
               write(iitem,602,advance='no') id(iVar)
            else if (id(iVar) > 0) then
               write(iitem,601,advance='no') id(iVar)
            end if
         end do
         write(iitem,603)
      end if
600   format('[',i1,';"',a,'";')
601   format('<',i1,'>')
602   format('<',i2,'>')
603   format(']')
    end subroutine writeItemGroup

  end subroutine writeStrainCoatHeader


  subroutine writeHistoryDB (rdb,strainCoats,ierr)

    !!==========================================================================
    !! Write history results to the results database file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 10 Dec 2003/1.0
    !!==========================================================================

    use RDBModule             , only : RDBType
    use StrainCoatModule      , only : StrainCoatType
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(RDBType)       , intent(inout) :: rdb
    type(StrainCoatType), intent(in)    :: strainCoats(:)
    integer             , intent(out)   :: ierr

    !! Local variables
    logical :: lDouble
    integer :: i

    !! --- Logic section ---

    call ffa_cmdlinearg_getbool ('double',lDouble)

    do i = 1, size(strainCoats)

       call writeRosetteDB (rdb,strainCoats(i)%tensorRosette%data,lDouble,ierr)
       if (ierr < 0) then
          call reportError (debugFileOnly_p,'writeHistoryDB')
          return
       end if

    end do

  end subroutine writeHistoryDB


  subroutine writeRosetteDB (rdb,ros,lDouble,ierr)

    !!==========================================================================
    !! Write strain rosette history results to the results database file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use RDBModule          , only : RDBType, writeRDB
    use StrainRosetteModule, only : StrainRosetteType
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(RDBType)          , intent(inout) :: rdb
    type(StrainRosetteType), intent(in)    :: ros(:)
    logical                , intent(in)    :: lDouble
    integer                , intent(out)   :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(ros)

       call writeRDB (rdb,ros(i)%sigmaP,ierr,lDouble)
       call writeRDB (rdb,ros(i)%tauMax,ierr,lDouble)
       call writeRDB (rdb,ros(i)%sigmaVM,ierr,lDouble)

       call writeRDB (rdb,ros(i)%epsP,ierr,lDouble)
       call writeRDB (rdb,ros(i)%gammaMax,ierr,lDouble)
       call writeRDB (rdb,ros(i)%epsVM,ierr,lDouble)

    end do
    if (ierr < 0) call reportError (debugFileOnly_p,'writeStrainRosetteDB')

  end subroutine writeRosetteDB


  subroutine writeElementsDB (rdb,strainCoats,time,ierr)

    !!==========================================================================
    !! Write element results to the results database file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use kindModule            , only : dp
    use RDBModule             , only : RDBType
    use StrainCoatModule      , only : StrainCoatType
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(RDBType)       , intent(inout) :: rdb
    type(StrainCoatType), intent(in)    :: strainCoats(:)
    real(dp)            , intent(in)    :: time
    integer             , intent(out)   :: ierr

    !! Local variables
    logical :: lDouble
    integer :: i

    !! --- Logic section ---

    call ffa_cmdlinearg_getbool ('double',lDouble)

    do i = 1, size(strainCoats)

       call writeStrainCoatDB (rdb, strainCoats(i)%results, time, &
            &                  strainCoats(i)%fppType > 0, lDouble, ierr)
       if (ierr < 0) then
          call reportError (debugFileOnly_p,'writeElementsDB')
          return
       end if

    end do

  end subroutine writeElementsDB


  subroutine writeStrainCoatDB (rdb,coats,time,lDamage,lDouble,ierr)

    !!==========================================================================
    !! Write strain coat results to the results database file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use kindModule       , only : dp
    use RDBModule        , only : RDBType, writeRDB
    use StrainCoatModule , only : CoatResultType, BiAxMean, BiAxStdDev
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType)       , intent(inout) :: rdb
    type(CoatResultType), intent(in)    :: coats(:)
    real(dp)            , intent(in)    :: time
    logical             , intent(in)    :: lDamage, lDouble
    integer             , intent(out)   :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: damag, lifeR, lifeT

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(coats)

       if (abs(coats(i)%sigMax) > abs(coats(i)%sigMin)) then
          call writeRDB (rdb,coats(i)%sigMax,ierr,lDouble)
       else
          call writeRDB (rdb,coats(i)%sigMin,ierr,lDouble)
       end if
       call writeRDB (rdb,coats(i)%tauMax,ierr,lDouble)
       call writeRDB (rdb,coats(i)%strRange(1),ierr,lDouble)
       call writeRDB (rdb,coats(i)%vmsMax,ierr,lDouble)

       if (abs(coats(i)%epsMax) > abs(coats(i)%epsMin)) then
          call writeRDB (rdb,coats(i)%epsMax,ierr,lDouble)
       else
          call writeRDB (rdb,coats(i)%epsMin,ierr,lDouble)
       end if
       call writeRDB (rdb,coats(i)%gammaMax,ierr,lDouble)
       call writeRDB (rdb,coats(i)%strRange(2),ierr,lDouble)
       call writeRDB (rdb,coats(i)%vmeMax,ierr,lDouble)

       if (lDamage) then

          if (coats(i)%damage > 0.0_dp) then
             damag = coats(i)%damage
             lifeR = 1.0_dp / damag
             lifeT = time / damag
          else ! no damage = infinite life, insert "special" value instead
             damag = 1.0e20_dp
             lifeR = 1.0e20_dp
             lifeT = 1.0e20_dp
          end if
          call writeRDB (rdb,damag,ierr,lDouble)
          call writeRDB (rdb,lifeT,ierr,lDouble)
          call writeRDB (rdb,lifeR,ierr,lDouble)

       end if

       call writeRDB (rdb,coats(i)%angSpread,ierr,lDouble)
       call writeRDB (rdb,coats(i)%popAngle,ierr,lDouble)

       if (coats(i)%nBiAxial > 0) then

          call writeRDB (rdb,BiAxMean(coats(i)),ierr,lDouble)
          call writeRDB (rdb,BiAxStdDev(coats(i)),ierr,lDouble)

       end if

       if (ierr < 0) then
          call reportError (debugFileOnly_p,'writeStrainCoatDB')
          return
       end if

    end do

  end subroutine writeStrainCoatDB

end module saveStrainCoatModule
