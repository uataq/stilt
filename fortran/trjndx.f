!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TRJNDX           READ THE TRAJECTORY ENDPOINT FILE
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            READS THE HEADER INFORMATION OF THE TRAJECTORY ENDPOINTS
!            FILE AND SETS PROGRAM VARIABLES ACCORDINGLY
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 16 Dec 1998 (RRD)
!                 05 Jul 2000 (RRD) - return NDIAG variable
!                 12 Dec 2000 (RRD) - fortran90 upgrade
!                 23 Oct 2002 (RRD) - vertical hieghts for color
!                 19 Nov 2003 (RRD) - return diagnostic label
!                 31 May 2005 (RRD) - updated trajectory file format
!                 06 Oct 2005 (RRD) - four digit year test
!                 19 Oct 2005 (BS)  - labeling IDLBL
!
! USAGE: CALL TRJNDX(NTRAJ,TFMT,OLAT,OLON,OLVL,IBYR,IBMO,IBDA,IBHR,KV,TKOL,
!                    JFYR,JFMO,JFDA,JFHR,MODEL,DIRCTN,MOTION,NDIAG,LDIAG,
!                    KTER,IDLBL)
!
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           UNIT 10 - trajectory input file
!   OUTPUT FILES:          NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE TRJNDX(NTRAJ,TFMT,OLAT,OLON,OLVL,IBYR,IBMO,IBDA,IBHR,KV,TKOL,     &
                  JFYR,JFMO,JFDA,JFHR,MODEL,DIRCTN,MOTION,NDIAG,LDIAG,KTER,  &
                  IDLBL)

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER,      INTENT(INOUT) :: NTRAJ                ! numb of traj in file 
  INTEGER,      INTENT(OUT)   :: TFMT                 ! trajectory file format
  REAL,         INTENT(OUT)   :: OLAT (:)             ! starting location 
  REAL,         INTENT(OUT)   :: OLON (:)             ! starting location 
  REAL,         INTENT(OUT)   :: OLVL (:)             ! starting location 
  INTEGER,      INTENT(OUT)   :: KV                   ! numb of diff levels 
  REAL,         INTENT(OUT)   :: TKOL (:)             ! height color index  
  INTEGER,      INTENT(OUT)   :: IBYR,IBMO,IBDA,IBHR  ! trajectory start   
  INTEGER,      INTENT(OUT)   :: JFYR,JFMO,JFDA,JFHR  ! meteo start time      
  CHARACTER(8), INTENT(OUT)   :: MODEL                ! meteo data file used
  CHARACTER(8), INTENT(OUT)   :: DIRCTN               ! direction
  CHARACTER(8), INTENT(OUT)   :: MOTION               ! vertical motion method
  INTEGER,      INTENT(OUT)   :: NDIAG                ! number of diagnostics
  CHARACTER(8), INTENT(OUT)   :: LDIAG (:)            ! diagnostic label
  INTEGER,      INTENT(OUT)   :: KTER                 ! terrain height index 
  CHARACTER(8), INTENT(OUT)   :: IDLBL                ! identifier label

!-------------------------------------------------------------------------------

  INTEGER :: nfhr, ngrd, kp, kg, knum, kret

!-------------------------------------------------------------------------------
! initialization section

! number of meteorological grids used in calculation
  READ(10,'(2I6)',IOSTAT=kret)NGRD,TFMT
  IF(kret.NE.0)TFMT=0  ! missing variable indicates old format specification

! for each grid load the model and data file starting time
  DO KG=1,NGRD
     IF(KG.EQ.1)THEN
!       only save information from first grid for plotting
        READ(10,'(A8,5I6)') MODEL,JFYR,JFMO,JFDA,JFHR,NFHR
     ELSE
        READ(10,*)
     END IF
  END DO
  JFYR=MOD(JFYR,100)

! number of different trajectories in the file 
!    (IDLBL only in new-format; written by cluster program trajmean, merglist)
  IF(tfmt.EQ.0)THEN
     READ(10,'(I6,2A8)')KNUM,DIRCTN,MOTION
  ELSE
     READ(10,'(I6,3(1X,A8))',IOSTAT=kret)KNUM,DIRCTN,MOTION,IDLBL
  END IF
 !IF(kret.NE.0)... there is no IDLBL meaning it's a regular tdump file

  IF(NTRAJ.EQ.0)THEN
!    first entry returns traj number for array allocation
     NTRAJ=KNUM
!    skip remaining starting location records
     DO KP=1,NTRAJ
        READ(10,'(A)')
     END DO
!    skip the diagnostics information
     READ(10,'(I6)') NDIAG
!    return later to continue reading
     RETURN
  END IF

!-------------------------------------------------------------------------------
! load starting point for each different trajectory

  kploop : DO KP=1,KNUM  
     IF(KP.EQ.1)THEN
        IF(tfmt.EQ.0)THEN
           READ(10,'(4I6,2F8.3,F8.1)')IBYR,IBMO,IBDA,IBHR,     &
                                      OLAT(KP),OLON(KP),OLVL(KP)
        ELSE
           READ(10,'(4I6,2F9.3,F8.1)')IBYR,IBMO,IBDA,IBHR,     &
                                      OLAT(KP),OLON(KP),OLVL(KP)
        END IF
     ELSEIF(KP.GT.KNUM)THEN
!       dimension based upon first file
        CYCLE kploop
     ELSE
        IF(tfmt.EQ.0)THEN
           READ(10,'(24X,2F8.3,F8.1)')OLAT(KP),OLON(KP),OLVL(KP)
        ELSE
           READ(10,'(24X,2F9.3,F8.1)')OLAT(KP),OLON(KP),OLVL(KP)
        END IF
     END IF
  END DO kploop

! number of diagnostic meteo variables and their label
  IF(tfmt.EQ.0)THEN
     READ(10,'(I6,99A8)')NDIAG,LDIAG
  ELSE
     READ(10,'(I6,99(1X,A8))')NDIAG,LDIAG
  END IF

! save the terrain height index number if it is available
  KTER=0
  DO KP=1,NDIAG
     IF(LDIAG(KP)(1:4).EQ.'TERR')KTER=KP
  END DO

!-------------------------------------------------------------------------------
! determine height index number for trajectory color

  KV=1
  TKOL(1)=OLVL(1)
  kvloop : DO KP=2,KNUM 
     DO KG=1,KV
!       match in tlevel array means this is not a new level
        IF(TKOL(KG).EQ.OLVL(KP)) CYCLE kvloop
     END DO
!    add new level to the array
     KV=KV+1
     TKOL(KV)=OLVL(KP)
  END DO kvloop

! replace height with index
  DO KP=KNUM,1,-1
     DO KG=1,KV
        IF(TKOL(KG).EQ.OLVL(KP)) TKOL(KP)=KG
     END DO
  END DO

END SUBROUTINE trjndx
