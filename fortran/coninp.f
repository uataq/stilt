!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONINP           READS THE GRIDDED HYSPLIT CONCENTRATION FILE
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ROUTINE TO READ GRIDDED CONCENTRATIONS AS USED BY MOST PLOT PROGRAMS
!   SELECTS ONE POLLUTANT FROM A MULTI-POLLUTANT FILE AND TRANSFERS DATA
!   AT ALL LEVELS FOR PROCESSING BY OTHER SUBROUTINES
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 09 Dec 1998 (RRD)
!                 20 Nov 2000 (RRD) - fortran90 upgrade
!                 01 Dec 2000 (RRD) - concentration packing
!                 21 Oct 2002 (RRD) - multiple species deposition summation
!                 22 Dec 2004 (RRD) - cpack=2 test
!
! USAGE: CALL CONINP(CPACK,CONC,TEMP,DEPT,HEIGHT,PTYPE,NTYP,NLVL,NUMP,NLAT,NLON,
!                    IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             unit 10 - gridded concentration data
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CONINP(CPACK,CONC,TEMP,DEPT,HEIGHT,PTYPE,NTYP,NLVL,NUMP,NLAT,NLON,  &
                  IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH)

  IMPLICIT NONE 

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,      INTENT(IN)    :: cpack                   ! packing flag   
  INTEGER,      INTENT(IN)    :: height (:)              ! levels of data
  INTEGER,      INTENT(IN)    :: ntyp                    ! number of types
  INTEGER,      INTENT(IN)    :: nlvl                    ! number of levels
  INTEGER,      INTENT(IN)    :: nump                    ! pollutant selected 
  INTEGER,      INTENT(IN)    :: nlat,nlon               ! input grid size 

  CHARACTER(4), INTENT(OUT)   :: ptype                   ! pollutant display
  INTEGER,      INTENT(OUT)   :: iyr,imo,ida,ihr,imn     ! start date/time
  INTEGER,      INTENT(OUT)   :: jyr,jmo,jda,jhr,jmn,jfh ! stop data/time
  REAL,         INTENT(OUT)   :: conc (:,:,:)            ! concentration array
  REAL,         INTENT(OUT)   :: dept (:,:)              ! deposition array
  REAL,         INTENT(OUT)   :: temp (:,:)              ! dummy array

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  CHARACTER(4)  :: PCHAR
  INTEGER(2)    :: IP,JP
  INTEGER       :: kp,kl,np,ii,jj,level,ifh,kret,nxyp

!-------------------------------------------------------------------------------

! sample start and stop times
  READ(10,END=500,ERR=500)IYR,IMO,IDA,IHR,IMN,IFH
  READ(10,END=500,ERR=500)JYR,JMO,JDA,JHR,JMN,JFH

  DO KP=1,NTYP
  DO KL=1,NLVL

     IF(CPACK.EQ.1)THEN
        TEMP = 0.0
        READ(10,IOSTAT=kret)PCHAR,LEVEL,NXYP,(IP,JP,TEMP(IP,JP),NP=1,NXYP)
     ELSEIF(CPACK.EQ.0)THEN
        READ(10,END=500,ERR=500)PCHAR,LEVEL,((TEMP(II,JJ),II=1,NLON),JJ=1,NLAT)
     ELSE
        WRITE(*,*)'*ERROR* coninp: unsupported packing option -',CPACK
        STOP 900 
     END IF

     IF(KP.EQ.NUMP)THEN
        PTYPE=PCHAR
     ELSEIF(NUMP.EQ.0)THEN
        PTYPE='SUM'
     END IF

!    if selected pollutant save values
     DO II=1,NLON
     DO JJ=1,NLAT

        IF( (KP.EQ.NUMP) .OR. (NUMP.EQ.0.AND.KP.EQ.1) )THEN
!          selected pollutant or initialize summation below
           CONC(II,JJ,KL)=TEMP(II,JJ)
        ELSEIF(NUMP.EQ.0)THEN
!          selected option to sum all pollutants
           CONC(II,JJ,KL)=CONC(II,JJ,KL)+TEMP(II,JJ)
        END IF

        IF(HEIGHT(KL).EQ.0)THEN
!          under all conditions sum deposition if available
           IF(KP.EQ.NUMP.OR.NUMP.EQ.0)DEPT(II,JJ)=DEPT(II,JJ)+TEMP(II,JJ)
        END IF

     END DO
     END DO

  END DO
  END DO

  RETURN

500 WRITE(*,*)'Premature concentration file termination'

END SUBROUTINE coninp
