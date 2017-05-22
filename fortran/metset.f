!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METSET           METeorological data structure SET
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL STRUCTURE SET READS AND UNPACKS THE METEOROLOGICAL
!   INDEX RECORD AND INITIALIZES THE FILE STRUCTURE DEFINITIONS.
!
! PROGRAM HISTORY LOG:
! LAST REVISION: 23 Dec 1998 (RRD)
!                20 Jan 1999 (RRD) - added status= on open statement
!                09 Feb 1999 (RRD) - test for max dir string length
!                29 Oct 1999 (RRD) - moved structure variables out of read
!                                  - fixed problem with single time files
!                26 Jul 2000 (RRD) - open meteo as read only
!                04 Sep 2000 (RRD) - fortran90 upgrade
!                29 Aug 2001 (RRD) - simultaneous multiple meteorology
!                23 Oct 2001 (RRD) - extended grid domains
!                09 Sep 2002 (RRD) - fortran coding standards
!                10 Apr 2003 (RRD) - added forecast hour decode 
!                09 Jun 2005 (RRD) - check for constant delta-t over file
!                20 Jul 2005 (RRD) - message file open kf=21 test
!                03 Jun 2008 (RRD) - embedded blanks dir/file
!                10 Dec 2008 (RRD) - common code for unregistered version
!                15 Jan 2009 (RRD) - replaced label with dname
!                13 Mar 2010 (TK) -  increase max header length from 5100 to 10000
!
! USAGE:  CALL METSET(KG,KT,KFOR)
!
!   INPUT ARGUMENT LIST:      see below
!   OUTPUT ARGUMENT LIST:     see below
!   INPUT FILES:              see funits module
!   OUTPUT FILES:             none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE METSET(KG,KT,KFOR)

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,       INTENT(IN)    :: kg           ! grid identfication number
  INTEGER,       INTENT(IN)    :: kt           ! data time period (1 or 2)
  INTEGER,       INTENT(IN)    :: kfor         ! 0=No_Forecasts 1=Use_all_Data    

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER,  PARAMETER   :: mlenh  = 10000  ! maximum length of header
  CHARACTER(mlenh)      :: header         ! extended index record header

  CHARACTER(2)   :: cgrid
  CHARACTER(4)   :: kvar
  CHARACTER(50)  :: label
  CHARACTER(80)  :: dname,fname

  LOGICAL        :: ftest

  INTEGER        :: n,nhl2,krec,nndx,nhl1,kol,nvar,k,l,nrec,nlvl
  INTEGER        :: kret,knx,kny,lenh,ldat,lrec,icx,klen,kunit,igrid

  INTEGER        :: DATE_TIME(8)
  CHARACTER(12)  :: REAL_CLOCK(3)
!dwen(20090822) ***************
  integer        :: ierr

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------


  INTERFACE
!-------------------------------------------------------------------------------
  SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
  INTEGER,  INTENT(OUT)  :: macc                 ! minutes since 1 Jan 1970
  END SUBROUTINE tm2min
!-------------------------------------------------------------------------------
  SUBROUTINE METOLD(KG,KT,KLEN)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)  :: kg       ! sequential grid identfication number
  INTEGER,  INTENT(IN)  :: kt       ! data time period identification (1 or 2)
  INTEGER,  INTENT(IN)  :: klen     ! length of directory string
  END SUBROUTINE metold
!-------------------------------------------------------------------------------
  END INTERFACE

!-------------------------------------------------------------------------------
! test for meteo file existence
!-------------------------------------------------------------------------------

  KLEN=LEN_TRIM(FILE(KG,KT)%DIR)
  IF(KLEN.LE.0)THEN
     WRITE(*,*)'Directory: ',FILE(KG,KT)%DIR
     WRITE(*,*)'String length exceeds character maximum'
     STOP 900
  END IF
  DNAME=FILE(KG,KT)%DIR
  FNAME=DNAME(1:KLEN)//FILE(KG,KT)%METEO
  INQUIRE(FILE=FNAME,EXIST=FTEST)
  IF(.NOT.FTEST)THEN
     WRITE(*,*)'Unable to find file: ',FILE(KG,KT)%METEO
     WRITE(*,*)'On local directory : ',FILE(KG,KT)%DIR(1:KLEN)
     WRITE(*,*)'Check input CONTROL file for correct values'
     STOP 900
  END IF

!-------------------------------------------------------------------------------
! set IO unit number
!-------------------------------------------------------------------------------

  KUNIT=FILE(KG,KT)%KUNIT

! open file to decode the standard label (50) plus the
! fixed portion (108) of the extended header
  OPEN(KUNIT, FILE=FNAME,STATUS='OLD',ACTION='READ',                           &
       RECL=158,ACCESS='DIRECT',FORM='UNFORMATTED')

! decode the standard portion of the index record
  READ(KUNIT,REC=1)LABEL,HEADER(1:108)
  READ(LABEL,'(5I2,2X,A2,A4)') FILE(KG,KT)%FIRST%YR, FILE(KG,KT)%FIRST%MO,     &
         FILE(KG,KT)%FIRST%DA, FILE(KG,KT)%FIRST%HR, FILE(KG,KT)%FIRST%IC,     &
         CGRID, KVAR

  IF(KVAR.NE.'INDX')THEN
     WRITE(*,*)'WARNING metset: Old format meteo data grid ',KG,KT
     CLOSE(KUNIT)
     CALL METOLD(KG,KT,KLEN)
     DREC(KG,KT)%TYPE=1
     RETURN
  END IF
  DREC(KG,KT)%TYPE=2

!-------------------------------------------------------------------------------
! decode extended portion of the header
!-------------------------------------------------------------------------------

  READ(HEADER(1:108),'(A4,I3,I2,12F7.0,3I3,I2,I4)')                          &
       GRID(KG,KT)%MODEL_ID, ICX,                   FILE(KG,KT)%FIRST%MN,    &
       GRID(KG,KT)%POLE_LAT, GRID(KG,KT)%POLE_LON,  GRID(KG,KT)%REF_LAT,     &
       GRID(KG,KT)%REF_LON,  GRID(KG,KT)%SIZE,      GRID(KG,KT)%ORIENT,      &
       GRID(KG,KT)%TANG_LAT, GRID(KG,KT)%SYNC_XP,   GRID(KG,KT)%SYNC_YP,     &
       GRID(KG,KT)%SYNC_LAT, GRID(KG,KT)%SYNC_LON,  GRID(KG,KT)%DUMMY,       &
       GRID(KG,KT)%NX,       GRID(KG,KT)%NY,        GRID(KG,KT)%NZ,          &
       DREC(KG,KT)%Z_FLAG,   LENH

! check number of levels against structure dimension
  IF((GRID(KG,KT)%NZ).GT.MLVL)THEN
     WRITE(*,*)'*ERROR* metset: exceeding DEFGRID dimension MLVL - ',mlvl
     WRITE(*,*)'        data file dimension - ',GRID(KG,KT)%NZ
     STOP 900
  END IF

! check header length against compiled dimension     
  IF(LENH.GT.MLENH)THEN
     WRITE(*,*)'*ERROR* metset: exceeding compiled header length: ',mlenh
     WRITE(*,*)'                                 Input data file: ',lenh
     STOP 900
  END IF

! grid id variable needed for old data sets or if either dimension of
! the grid domain is extended beyond 3 digits
  KNX=ICHAR(CGRID(1:1))
  KNY=ICHAR(CGRID(2:2))
  IF(KNX.GE.64.OR.KNY.GE.64)THEN
     GRID(KG,KT)%NX=(KNX-64)*1000+GRID(KG,KT)%NX
     GRID(KG,KT)%NY=(KNY-64)*1000+GRID(KG,KT)%NY
     GRID(KG,KT)%NUMBER=KNX*10+KNY
  ELSE
     READ(CGRID,'(I2)')IGRID
     GRID(KG,KT)%NUMBER=IGRID
  END IF

! close file and reopen with proper length
  CLOSE (KUNIT)
!************************************
!dwen(20090817)  LDAT = GRID(KG,KT)%NX*GRID(KG,KT)%NY
      !  change to 16 bit for RAMS and ECMWF (2 bytes rqd)
      IF (GRID(KG,kt)%MODEL_ID  == 'RAMS'  .OR.  GRID(KG,kt)%MODEL_ID(1:3)  == 'ECX' .or. &
           & GRID(KG,kt)%MODEL_ID  == 'DWRF') THEN
         LDAT=GRID(KG,kt)%NX*GRID(KG,kt)%NY*2
      ELSE
         LDAT=GRID(KG,kt)%NX*GRID(KG,kt)%NY
      END IF
!********************************* 

  FILE(KG,KT)%REC_LEN = LDAT+50
  LREC = FILE(KG,KT)%REC_LEN
  OPEN(KUNIT, FILE=FNAME,STATUS='OLD',ACTION='READ',RECL=LREC,  &
       ACCESS='DIRECT',FORM='UNFORMATTED')

! determine number of index records
  NNDX=LENH/LDAT+1

! read extended character string over multiple index records
  NHL1=1
  KREC=1
  DO N=1,NNDX
     NHL2=NHL1+LDAT-1
     IF(N.EQ.NNDX)NHL2=NHL1+(LENH-(NNDX-1)*LDAT)-1
     READ(KUNIT,REC=KREC)LABEL,HEADER(NHL1:NHL2)
     KREC=KREC+1
     NHL1=NHL2+1
  END DO

! loop through and decode the remainder of the index string
  KOL=109
  NREC=NNDX
  NLVL=GRID(KG,KT)%NZ
  DO L=1,NLVL
     READ(HEADER(KOL:KOL+7),'(F6.2,I2)')              &
          DREC(KG,KT)%HEIGHT(L),DREC(KG,KT)%NUM_VARB(L)
     KOL=KOL+8

     NVAR=DREC(KG,KT)%NUM_VARB(L)
!    check against compiled dimension     
     IF(NVAR.GT.MVAR)THEN
        WRITE(*,*)'*ERROR* metset: exceeding DEFGRID dimension MVAR - ',MVAR  
        STOP 900
     END IF

     DO K=1,NVAR
!****************************
! CHG (10/01/03) change from 8 to 16 bit
            IF(GRID(KG,kt)%MODEL_ID /= 'RAMS' .AND. GRID(KG,kt)%MODEL_ID(1:3) /= 'ECX' .AND. &
                 & GRID(KG,kt)%MODEL_ID  /= 'DWRF')THEN
               READ(HEADER(KOL:KOL+7),'(A4,I3)',iostat=ierr) DREC(KG,kt)%VARB_ID(K,L), DREC(KG,kt)%CHK_SUM(K,L)
               if (ierr .ne. 0) then
                  DREC(KG,kt)%CHK_SUM(K,L) = 0
                  READ(HEADER(KOL:KOL+7),'(A4,I3)',iostat=ierr) DREC(KG,kt)%VARB_ID(K,L)
                  if (ierr .ne. 0) stop 'variable ID read error in metset'
               endif
               KOL=KOL+8
            ELSE
               READ(HEADER(KOL:KOL+9),'(A4,I5)',iostat=ierr) DREC(KG,kt)%VARB_ID(K,L), DREC(KG,kt)%CHK_SUM(K,L)
               if (ierr .ne. 0) then
                  DREC(KG,kt)%CHK_SUM(K,L) = 0
                  READ(HEADER(KOL:KOL+9),'(A4,I5)',iostat=ierr) DREC(KG,kt)%VARB_ID(K,L)
                  if (ierr .ne. 0) stop 'variable ID read error in metset'
               endif
               KOL=KOL+10
            END IF
!********************************* 
!dwen(20090817)        READ(HEADER(KOL:KOL+7),'(A4,I3)')                    &
!dwen(20090817)             DREC(KG,KT)%VARB_ID(K,L),DREC(KG,KT)%CHK_SUM(K,L)
!dwen(20090817)        KOL=KOL+8
        NREC=NREC+1
     END DO
  END DO
  DREC(KG,KT)%REC_PER=NREC
  DREC(KG,KT)%OFFSET=0

! decode the forecast hour and minutes from the extended header
  READ(HEADER(1:108),'(4X,I3,I2)') FILE(KG,KT)%FIRST%IC, FILE(KG,KT)%FIRST%MN

! set the year's accumulated clock time in minutes
! for the first time period on the file
  CALL TM2MIN(FILE(KG,KT)%FIRST%YR, FILE(KG,KT)%FIRST%MO,                  &
              FILE(KG,KT)%FIRST%DA, FILE(KG,KT)%FIRST%HR,                  &
              FILE(KG,KT)%FIRST%MN, FILE(KG,KT)%FIRST%MACC)

!-------------------------------------------------------------------------------
! skip to the next time period index record to find the time
! interval between data periods (minutes)
!-------------------------------------------------------------------------------

  NREC=NREC+1

  READ(KUNIT,REC=NREC,ERR=900)LABEL,HEADER(1:108)
  READ(LABEL,'(5I2,4X,A4)') FILE(KG,KT)%LAST%YR,   FILE(KG,KT)%LAST%MO,        &
       FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,   FILE(KG,KT)%LAST%IC, KVAR
  IF(KVAR.NE.'INDX')THEN 
     WRITE(*,*)'*ERROR* metset: 2nd time period INDX record missing'
     RETURN
  END IF

! decode the forecast hour and minutes from the extended header
  READ(HEADER(1:108),'(4X,I3,I2)') FILE(KG,KT)%LAST%IC, FILE(KG,KT)%LAST%MN

! set the year's accumulated clock time in minutes
! for the next time period on the file
  CALL TM2MIN(FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO, FILE(KG,KT)%LAST%DA,   &
              FILE(KG,KT)%LAST%HR, FILE(KG,KT)%LAST%MN, FILE(KG,KT)%LAST%MACC)

! compute the minute difference between time periods
  DREC(KG,KT)%DELTA=FILE(KG,KT)%LAST%MACC-FILE(KG,KT)%FIRST%MACC
  K=FILE(KG,KT)%LAST%MACC

!-------------------------------------------------------------------------------
! continue on to find the last record
!-------------------------------------------------------------------------------

    NREC=NREC+DREC(KG,KT)%REC_PER
    eloop : DO WHILE (NREC.GE.1)

       READ(KUNIT,REC=NREC,IOSTAT=KRET) LABEL, HEADER(1:108)
       IF(KRET.NE.0)EXIT eloop
       READ(LABEL,'(5I2,4X,A4)',IOSTAT=KRET)                                   &
            FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO, FILE(KG,KT)%LAST%DA,     &
            FILE(KG,KT)%LAST%HR, FILE(KG,KT)%LAST%IC, KVAR
       IF(KRET.NE.0)EXIT eloop
       IF(KVAR.NE.'INDX') EXIT eloop
       READ(HEADER(1:108),'(7X,I2)')FILE(KG,KT)%LAST%MN
       NREC=NREC+DREC(KG,KT)%REC_PER

!      insure that time differences are the same through the file
       CALL TM2MIN(FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO,                   &
                   FILE(KG,KT)%LAST%DA, FILE(KG,KT)%LAST%HR,                   &
                   FILE(KG,KT)%LAST%MN, N)
       IF((N-K).NE.DREC(KG,KT)%DELTA)THEN
          WRITE(*,*)'*ERROR* metset: meteorological data time interval varies' 
          WRITE(*,*)' Changed from ',DREC(KG,KT)%DELTA,' min to ',(N-K),' min'
          WRITE(*,*)' At day/hr ',FILE(KG,KT)%LAST%DA,FILE(KG,KT)%LAST%HR
          STOP 900          
       ELSE
          K=N
       END IF
       
    END DO eloop
    FILE(KG,KT)%ENDREC=NREC-1

!   set the year's accumulated clock time in minutes
!   for the last time period on the file
    CALL TM2MIN(FILE(KG,KT)%LAST%YR, FILE(KG,KT)%LAST%MO, FILE(KG,KT)%LAST%DA, &
                FILE(KG,KT)%LAST%HR, FILE(KG,KT)%LAST%MN, FILE(KG,KT)%LAST%MACC)

!   current processor clock time
    real_clock(:)=' '
    CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),DATE_TIME)
    INQUIRE(UNIT=KF21,OPENED=FTEST)
    IF(FTEST)THEN
       WRITE(KF21,*)'Simulation Date (CCYYMMDD): ',REAL_CLOCK(1)
       WRITE(KF21,*)'Simulation Time (HHMMSS.S): ',REAL_CLOCK(2)
    END IF

!   trial version restriction, determine if forecast data are used
!   when kfor=0 cannot use forecast data; kfor=1 no restrictions
    IF(KFOR.EQ.0)THEN
       IF(FILE(KG,KT)%LAST%IC.GT.0)THEN
          CALL TM2MIN(DATE_TIME(1),DATE_TIME(2),DATE_TIME(3),DATE_TIME(5),   &
                      DATE_TIME(6),N)
!         last time in meteorology file greater than the current time is
!         not permitted in the trial version for unrestricted distribution
          IF(N.LT.FILE(KG,KT)%LAST%MACC)THEN           
             WRITE(*,*)'*ERROR* metset: the use of forecast meteorological data'
             WRITE(*,*)'  is not supported with this version. User registration'
             WRITE(*,*)'  is required for computations with forecast data files.'
             STOP 900      
          END IF
       END IF
    END IF
    RETURN

!   abnormal terminations
900 WRITE(*,*)'WARNING metset: Only one time period of meteo data'
    FILE(KG,KT)%LAST=FILE(KG,KT)%FIRST
    FILE(KG,KT)%ENDREC=NREC-1

!   set to dummy value to avoid division by zero
    DREC(KG,KT)%DELTA=1
    RETURN

END SUBROUTINE metset
