!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONSET           CONcentration grid SET data entry
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION GRID SET IS THE DATA ENTRY FOR SAMPLING GRID LOCATION
!   LEVELS, AND TIME INFORMATION FOR START, STOP, AND DURATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Feb 1998 (RRD)
!                 07 May 1998 (RRD) - y2k mod
!                 20 Oct 1999 (RRD) - concentration grid dateline correction
!                 08 Nov 1999 (RRD) - concentration grid span correction
!                 17 Nov 2000 (RRD) - fortran90 upgrade
!                 21 Feb 2001 (RRD) - global longitude grid
!                 12 Oct 2001 (RRD) - fixed relative termimation time
!                 24 Oct 2001 (RRD) - simplified backward dispersion
!                 02 Nov 2001 (RRD) - shortened averaging time option
!                 18 Dec 2001 (RRD) - insure proper global grid
!                 11 Feb 2001 (RRD) - formatted STARTUP file
!                 12 Jul 2002 (RRD) - input data decoder
!                 30 Aug 2002 (RRD) - global grid span patch
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 03 Dec 2002 (RRD) - global grid span patch#2
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 13 Oct 2004 (RRD) - test for duplicate file names
!                 22 Dec 2004 (RRD) - non regular concentration grid
!                 06 Apr 2005 (BJS) - off grid source warning msg fix
!                 01 Jun 2005 (RRD) - removed minutes warning message
!                 04 Oct 2005 (RRD) - averaging minutes test if hr >1
!                 17 Jan 2007 (RRD) - leave <0 snap values <0
!                 19 Mar 2007 (RRD) - maximum for snap or average
!                 06 Jul 2007 (RRD) - heights real-integer
!                 30 Nov 2007 (RRD) - global definition 181->180
!                 03 Jun 2008 (RRD) - embedded blanks dir/file
!
! USAGE:  CALL CONSET(CONC,ZMDL,NUMGRD,IUNIT,CGSIZE,CPACK,OLAT,OLON,   &
!                     IBYR,IBMO,IBDA,IBHR,BACK)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            unit 5 or unit KF21 input CONTROL data file
!   OUTPUT FILES:           unit KF22 when input unit=5 then write STARTUP file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CONSET(CONC,ZMDL,NUMGRD,IUNIT,CGSIZE,CPACK,OLAT,OLON,     &
                  IBYR,IBMO,IBDA,IBHR,BACK)

  USE funits

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'         ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------
 
  TYPE(cset), INTENT(OUT)   :: conc(:)              ! each concentration grid 
  REAL,       INTENT(IN)    :: zmdl                 ! model top meters AGL
  INTEGER,    INTENT(IN)    :: numgrd               ! number concen grids
  INTEGER,    INTENT(IN)    :: iunit                ! unit number for input
  REAL,       INTENT(INOUT) :: cgsize               ! minimum grid spacing (km) 
  INTEGER,    INTENT(IN)    :: cpack                ! concentration packing
  REAL,       INTENT(IN)    :: olat, olon           ! starting point
  INTEGER,    INTENT(IN)    :: ibyr,ibmo,ibda,ibhr  ! starting time
  LOGICAL,    INTENT(IN)    :: back                 ! integration direction

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER,  ALLOCATABLE :: temp(:)
  REAL,    PARAMETER :: ympd = 111.198323        ! km/deg-lat 
  REAL,    PARAMETER :: dgpr = 57.295828         ! deg/radian
  INTEGER, PARAMETER :: maxxp = 360
  INTEGER, PARAMETER :: maxyp = 181

  INTEGER           :: ii,jj,imo,iyr,ida,imn,ihr,kl,kk,nlvl,knum,kdir
  REAL              :: clat2,clon2,xmpd,cnlat,cnlon,splonm,splatm,splat,splon

!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
  INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time
  END SUBROUTINE tm2day
!-------------------------------------------------------------------------------
  SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
  INTEGER,  INTENT(OUT)  :: macc                 ! minutes since 1 Jan 1970
  END SUBROUTINE tm2min
!-------------------------------------------------------------------------------
  SUBROUTINE DECODI(IUNIT,VAR1,VAR2,VAR3,VAR4,VAR5)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IUNIT   ! unit number
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR1
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR2
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR3
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR4
  INTEGER, OPTIONAL, INTENT(INOUT) :: VAR5
  END SUBROUTINE decodi
  END INTERFACE

!-------------------------------------------------------------------------------
! primary loop for the number of concentration grids
!-------------------------------------------------------------------------------

  DO KK=1,NUMGRD

     IF(IUNIT.EQ.5.AND.KK.GT.1)THEN
!       copy initialization from previous input
        CONC(KK)=CONC(KK-1)

     ELSE
!       generic defaults
        CNLAT=OLAT
        CNLON=OLON
        CONC(KK)%DELT_LAT=1.0
        CONC(KK)%DELT_LON=1.0
        CONC(KK)%LEVELS=1
        CONC(KK)%HEIGHT(1)=50

        CONC(KK)%START%YR=IBYR
        CONC(KK)%START%MO=IBMO
        CONC(KK)%START%DA=IBDA
        CONC(KK)%START%HR=IBHR
        CONC(KK)%START%MN=0
        CONC(KK)%START%IC=0

        IF(BACK)THEN
           CONC(KK)%STOP%YR=MOD(IBYR-1,100)
        ELSE
           CONC(KK)%STOP%YR=MOD(IBYR+1,100)
        END IF
        CONC(KK)%STOP%MO=IBMO
        CONC(KK)%STOP%DA=IBDA
        CONC(KK)%STOP%HR=IBHR
        CONC(KK)%STOP%MN=0
        CONC(KK)%STOP%IC=0 

        CONC(KK)%SNAP=0
        CONC(KK)%DELTA%HR=24
        CONC(KK)%DELTA%MN=0

        CONC(KK)%DIR='|main|sub|output|'
        CONC(KK)%FILE='file_name'
     END IF

!-------------------------------------------------------------------------------
! horizontal grid definition
!-------------------------------------------------------------------------------

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Enter the follow for grid #:',KK
     WRITE(*,*)'Center Latitude, Longitude'
     WRITE(*,*)CNLAT,CNLON
     READ(IUNIT,*)CNLAT,CNLON
  ELSE
     READ(IUNIT,*)CNLAT,CNLON
     IF(CNLAT.EQ.0.0.AND.CNLON.EQ.0.0)THEN
        CNLAT=OLAT
        CNLON=OLON
     END IF
  END IF
  IF(IUNIT.EQ.5)WRITE(KF22,'(2F10.3)')CNLAT,CNLON

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Grid spacing (deg) Latitude, Longitude'
     WRITE(*,*)CONC(KK)%DELT_LAT, CONC(KK)%DELT_LON
  END IF
  READ(IUNIT,*)CONC(KK)%DELT_LAT, CONC(KK)%DELT_LON
  IF(IUNIT.EQ.5)WRITE(KF22,'(2F10.3)')CONC(KK)%DELT_LAT, CONC(KK)%DELT_LON

! maximum concentration grid span
  SPLATM=MAXYP*CONC(KK)%DELT_LAT
  SPLONM=MAXXP*CONC(KK)%DELT_LON
  IF(IUNIT.EQ.5)THEN
     SPLAT=SPLATM
     SPLON=SPLONM
     WRITE(*,*)'Grid span (deg) Latitude, Longitide'
     WRITE(*,*)SPLAT,SPLON
     READ(IUNIT,*)SPLAT,SPLON
  ELSE
     READ(IUNIT,*)SPLAT,SPLON
  END IF
  IF(SPLAT.EQ.0.0) SPLAT=SPLATM
  IF(SPLON.EQ.0.0) SPLON=SPLONM
  IF(IUNIT.EQ.5)WRITE(KF22,'(2F10.3)')SPLAT,SPLON

! determine minimum grid spacing for all grids
! km/degree-longitude
  XMPD=YMPD*COS(CNLAT/DGPR)
  CONC(KK)%SIZE=MIN(XMPD*CONC(KK)%DELT_LON,YMPD*CONC(KK)%DELT_LAT)
! return minimum spacing over all grids
  CGSIZE=MIN(CONC(KK)%SIZE,CGSIZE)

! compute lower left corner
  CONC(KK)%X1Y1_LAT=MAX(-90.0,CNLAT-SPLAT/2.0)
  CONC(KK)%X1Y1_LON=CNLON-SPLON/2.0

! compute upper right corner
  CLAT2=MIN(90.0,CNLAT+SPLAT/2.0)
  CLON2=CNLON+SPLON/2.0

! special case for global grid
  IF(SPLON.GE.360.0)THEN
     CONC(KK)%X1Y1_LON=-180.0
     CLON2=180.0-CONC(KK)%DELT_LON
     CNLON=0.0
!    recompute span
     SPLON=CLON2-CONC(KK)%X1Y1_LON+1.0
  END IF

  IF(SPLAT.GE.180.0)THEN
     CONC(KK)%X1Y1_LAT=-90.0
     CLAT2=90.0
     CNLAT=0.0
!    recompute span
     SPLAT=CLAT2-CONC(KK)%X1Y1_LAT+1.0
  END IF

! determine number of grid points (not to exceed dimension)
  CONC(KK)%NUMB_LAT=1+NINT(SPLAT/CONC(KK)%DELT_LAT)
  CONC(KK)%NUMB_LON=1+NINT(SPLON/CONC(KK)%DELT_LON)

! dateline correction (RRD - 20/10/99)
  IF(CONC(KK)%X1Y1_LON.LT.-180.0) CONC(KK)%X1Y1_LON=CONC(KK)%X1Y1_LON+360.0

!-------------------------------------------------------------------------------
! Redefine concentration grid for non-regular grids (cpack=2). The span will
! be set such that the number of grid points equals one and the grid center
! will be moved to corresond with the corner point.
!-------------------------------------------------------------------------------

  IF(CPACK.EQ.2)THEN
     CONC(KK)%X1Y1_LAT=CNLAT
     CONC(KK)%X1Y1_LON=CNLON
     CONC(KK)%NUMB_LAT=1
     CONC(KK)%NUMB_LON=1
  END IF

!-------------------------------------------------------------------------------
! output file names and directory
!-------------------------------------------------------------------------------

  IF(IUNIT.EQ.5)THEN
     WRITE(*,'(A,I2,A)')' Enter grid #',KK,' directory (|...|)'
     KDIR=LEN_TRIM(ADJUSTL(CONC(KK)%DIR))
     WRITE(*,'(A)')CONC(KK)%DIR(1:KDIR)
  END IF
  READ(IUNIT,'(A)')CONC(KK)%DIR
  CONC(KK)%DIR=ADJUSTL(CONC(KK)%DIR)
  KDIR=LEN_TRIM(CONC(KK)%DIR)
  IF(IUNIT.EQ.5)WRITE(KF22,'(A)')CONC(KK)%DIR(1:KDIR)

  IF(IUNIT.EQ.5)THEN
     WRITE(*,'(A,I2,A)')' Enter grid #',KK,' file name (?????)'
     KNUM=LEN_TRIM(ADJUSTL(CONC(KK)%FILE))
     WRITE(*,'(A)')CONC(KK)%FILE(1:KNUM)
  END IF
  READ(IUNIT,'(A)')CONC(KK)%FILE
  CONC(KK)%FILE=ADJUSTL(CONC(KK)%FILE)
  KNUM=LEN_TRIM(CONC(KK)%FILE)
  IF(IUNIT.EQ.5)WRITE(KF22,'(A)')CONC(KK)%FILE(1:KNUM)

! test for duplicate directory and file names
  IF(KK.GE.2.AND.CPACK.NE.2)THEN
     IF(CONC(KK)%DIR(1:KDIR).EQ.CONC(KK-1)%DIR(1:KDIR).AND.   &
        CONC(KK)%FILE(1:KNUM).EQ.CONC(KK-1)%FILE(1:KNUM))THEN
        WRITE(*,*)'*ERROR* conset: Concentration dir/file names identical'
        WRITE(*,*)(KK  ),CONC(KK)%DIR(1:KDIR),  CONC(KK)%FILE(1:KNUM)
        WRITE(*,*)(KK-1),CONC(KK-1)%DIR(1:KDIR),CONC(KK-1)%FILE(1:KNUM)
        STOP 900
     END IF
  END IF

!-------------------------------------------------------------------------------
! vertical grid levels or spacing
!-------------------------------------------------------------------------------

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Number of vertical concentration levels'
     WRITE(*,*)CONC(KK)%LEVELS
  END IF
  READ(IUNIT,*)CONC(KK)%LEVELS
  IF(IUNIT.EQ.5)WRITE(KF22,'(I2.2)')CONC(KK)%LEVELS

! test limits
  IF((CONC(KK)%LEVELS).GT.99)THEN
      WRITE(*,*)'*ERROR* conset: Number of levels exceed limit'
      STOP 900
  END IF

  NLVL=CONC(KK)%LEVELS
  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Height of each level m AGL'
     WRITE(*,*)(CONC(KK)%HEIGHT(KL),KL=1,NLVL)
  END IF

! more generic section to avoid crash (06 Jul 2007)
! READ(IUNIT,*)(CONC(KK)%HEIGHT(KL),KL=1,NLVL)
  ALLOCATE (temp(nlvl))
  READ(IUNIT,*)temp

! IF(IUNIT.EQ.5) WRITE(KF22,'(20I6)')(CONC(KK)%HEIGHT(KL),KL=1,NLVL)
  IF(IUNIT.EQ.5) WRITE(KF22,'(20I6)') temp                           

  DO KL=1,NLVL
     CONC(KK)%HEIGHT(KL)=TEMP(KL) 
!    note level height "0" assumed to be for deposition
     IF(CONC(KK)%HEIGHT(KL).EQ.0.AND.KL.NE.1)THEN
        WRITE(*,*)'*ERROR* conset: height=0 should be level=1'
        STOP 900
     END IF
     IF(CONC(KK)%HEIGHT(KL).GT.INT(ZMDL))THEN
        WRITE(*,*)'*ERROR* conset: level above model top'
        STOP 900
     END IF
  END DO
  DEALLOCATE (temp)

!-------------------------------------------------------------------------------
! set sampling time intervals
!-------------------------------------------------------------------------------

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Sampling start time: year month day hour minute'
     WRITE(*,*)CONC(KK)%START%YR, CONC(KK)%START%MO,                    &
               CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN
     READ(*,*) CONC(KK)%START%YR, CONC(KK)%START%MO,                    &
               CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN
  ELSE

!    READ(IUNIT,*)IYR,IMO,IDA,IHR,IMN
     CALL DECODI(IUNIT,IYR,IMO,IDA,IHR,IMN)

     IF(IMO.EQ.0)THEN
        IF(BACK)THEN
           IDA=-ABS(IDA)
           IHR=-ABS(IHR)
        END IF
!       relative from file start
        CONC(KK)%START%DA=IBDA+IDA
        CONC(KK)%START%HR=IBHR+IHR
     ELSE
        CONC(KK)%START%YR=IYR
        CONC(KK)%START%MO=IMO
        CONC(KK)%START%DA=IDA
        CONC(KK)%START%HR=IHR
        CONC(KK)%START%MN=IMN
     END IF
  END IF

  IF(IUNIT.EQ.5)WRITE(KF22,'(5I3)')CONC(KK)%START%YR, CONC(KK)%START%MO,  &
              CONC(KK)%START%DA, CONC(KK)%START%HR, CONC(KK)%START%MN

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Sampling stop time: year month day hour minute'
     WRITE(*,*)CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,                      &
               CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN
     READ(*,*) CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,                      &
               CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN
  ELSE

!    READ(IUNIT,*)IYR,IMO,IDA,IHR,IMN
     CALL DECODI(IUNIT,IYR,IMO,IDA,IHR,IMN)

     IF(IYR+IMO+IDA+IHR+IMN.EQ.0)THEN
!       accept defaults
        CONTINUE 
     ELSEIF(IMO.EQ.0)THEN
        CONC(KK)%STOP%YR=IBYR
        IF(BACK)THEN
           IDA=-ABS(IDA)
           IHR=-ABS(IHR)
        END IF
!       relative from file start
        CONC(KK)%STOP%DA=IBDA+IDA
        CONC(KK)%STOP%HR=IBHR+IHR
     ELSE
!       direct from input - all values specified
        CONC(KK)%STOP%YR=IYR
        CONC(KK)%STOP%MO=IMO
        CONC(KK)%STOP%DA=IDA
        CONC(KK)%STOP%HR=IHR
        CONC(KK)%STOP%MN=IMN
     END IF
  END IF

  IF(IUNIT.EQ.5)WRITE(KF22,'(5I3)') CONC(KK)%STOP%YR, CONC(KK)%STOP%MO,  &
                CONC(KK)%STOP%DA, CONC(KK)%STOP%HR, CONC(KK)%STOP%MN

  IF(IUNIT.EQ.5)THEN
      WRITE(*,*)'Sampling interval: type hour minute'
      WRITE(*,*)CONC(KK)%SNAP, CONC(KK)%DELTA%HR, CONC(KK)%DELTA%MN
  END IF
  READ(IUNIT,*)CONC(KK)%SNAP, CONC(KK)%DELTA%HR, CONC(KK)%DELTA%MN

  IF(IUNIT.EQ.5)WRITE(KF22,'(3I3)')                                      &
     CONC(KK)%SNAP,CONC(KK)%DELTA%HR,CONC(KK)%DELTA%MN

!-------------------------------------------------------------------------------
! convert dates to accumulated time
!-------------------------------------------------------------------------------

  CALL TM2MIN(CONC(KK)%START%YR, CONC(KK)%START%MO, CONC(KK)%START%DA,  &
              CONC(KK)%START%HR, CONC(KK)%START%MN, CONC(KK)%START%MACC)

  CALL TM2MIN(CONC(KK)%STOP%YR, CONC(KK)%STOP%MO, CONC(KK)%STOP%DA,     &
              CONC(KK)%STOP%HR, CONC(KK)%STOP%MN, CONC(KK)%STOP%MACC)

! convert back to date to fix date+delta errors
  CALL TM2DAY(CONC(KK)%START%MACC, CONC(KK)%START%YR, CONC(KK)%START%MO, &
              CONC(KK)%START%DA,   CONC(KK)%START%HR, CONC(KK)%START%MN)

  CALL TM2DAY(CONC(KK)%STOP%MACC,  CONC(KK)%STOP%YR,  CONC(KK)%STOP%MO,  &
              CONC(KK)%STOP%DA,    CONC(KK)%STOP%HR,  CONC(KK)%STOP%MN)

!-------------------------------------------------------------------------------
! set remaining variables
!-------------------------------------------------------------------------------

! sampling intervals in years, months, days not used
  CONC(KK)%DELTA%YR=0
  CONC(KK)%DELTA%MO=0
  CONC(KK)%DELTA%DA=0

! convert sampling interval hours-minutes to minutes
  IF(CONC(KK)%DELTA%HR.GE.1.AND.CONC(KK)%DELTA%MN.NE.0)THEN
     WRITE(*,*)'WARNING conset: averaging >1 hr requires zero minutes field'
     WRITE(*,*)'Resetting averaging minutes to zero on grid #',KK
     CONC(KK)%DELTA%MN=0
  END IF
  CONC(KK)%DELTA%MACC=60*CONC(KK)%DELTA%HR+CONC(KK)%DELTA%MN

! save sampling start time in conc output file marker variable
  CONC(KK)%NOW=CONC(KK)%START

!-------------------------------------------------------------------------------
! check for consistency on snapshot and maximum value grids
!-------------------------------------------------------------------------------

! special case when SNAP<0 it specifies the averaging time in hours
! and DELTA becomes the interval at which that average is output
! modified 1/10/2007 to remain negative
  IF(CONC(KK)%SNAP.LT.0) CONC(KK)%SNAP=60*CONC(KK)%SNAP

! defining a maximum concentration grid requires a snapshot or average grid
  IF(CONC(KK)%SNAP.EQ.2)THEN 
     KNUM=0
     DO KL=1,NUMGRD 
        IF(CONC(KL)%SNAP.EQ.1)KNUM=KNUM+1
        IF(CONC(KL)%SNAP.EQ.2)KNUM=KNUM+1
     END DO
     IF(KNUM.EQ.0)THEN
        WRITE(*,*)'*ERROR* conset: maximum concentration grid defined for #',KK
        WRITE(*,*)'Requires definition of an identical snapshot or average grid!'
        STOP
     END IF
     IF(KNUM.GT.2)THEN
        WRITE(*,*)'*ERROR* conset: maximum concentration grid defined for #',KK
        WRITE(*,*)'Requires one snapshot or average grid ... not multiple!'
        STOP
     END IF
  END IF

!-------------------------------------------------------------------------------
! check if concentration grid covers source location (not required)
!-------------------------------------------------------------------------------

  II=NINT(1.0+(OLON-CONC(KK)%X1Y1_LON)/CONC(KK)%DELT_LON)
     IF(OLON.LT.0.0.AND.CONC(KK)%X1Y1_LON.GT.0.0)              &
     II=NINT(1.0+(OLON-CONC(KK)%X1Y1_LON+360)/CONC(KK)%DELT_LON)
  JJ=NINT(1.0+(OLAT-CONC(KK)%X1Y1_LAT)/CONC(KK)%DELT_LAT)
  IF(JJ.GT.CONC(KK)%NUMB_LAT.OR.JJ.LT.1.OR.      &
     II.GT.CONC(KK)%NUMB_LON.OR.II.LT.1) THEN
     WRITE(*,*)'WARNING conset: source point not on concentration grid #',KK
  END IF

END DO
!dwen(20090316):add
! JCL:
      WRITE(KFJCLMSG,*)'CONC(1)%START%YR:  ',CONC(1)%START%YR
      WRITE(KFJCLMSG,*)'CONC(1)%START%MO:  ',CONC(1)%START%MO
      WRITE(KFJCLMSG,*)'CONC(1)%START%DA:  ',CONC(1)%START%DA
      WRITE(KFJCLMSG,*)'CONC(1)%START%HR:  ',CONC(1)%START%HR
      WRITE(KFJCLMSG,*)'CONC(1)%START%MN:  ',CONC(1)%START%MN
      WRITE(KFJCLMSG,*)'CONC(1)%STOP%YR:  ',CONC(1)%STOP%YR
      WRITE(KFJCLMSG,*)'CONC(1)%STOP%MO:  ',CONC(1)%STOP%MO
      WRITE(KFJCLMSG,*)'CONC(1)%STOP%DA:  ',CONC(1)%STOP%DA
      WRITE(KFJCLMSG,*)'CONC(1)%STOP%HR:  ',CONC(1)%STOP%HR
      WRITE(KFJCLMSG,*)'CONC(1)%STOP%MN:  ',CONC(1)%STOP%MN

      WRITE(KFJCLMSG,*)'CONC(1)%START%MACC:  ',CONC(1)%START%MACC
      WRITE(KFJCLMSG,*)'CONC(1)%STOP%MACC:  ',CONC(1)%STOP%MACC

END SUBROUTINE conset
