!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSINI           EMiSsion INItialization of gridded data
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:00-09-04
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION INITIALIZATION OPENS THE GRIDDED EMISSION INVENTORY FILE.
!   DETERMINES ARRAY DIMENSIONS FOR INTERNAL EMISSION GRID BASED UPON
!   THE INDEX RECORD OF THE EMISSION FILE: emission.txt
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 04 Sep 2000 (RRD) - initial version
!                 16 Mar 2001 (RRD) - qfile initialization
!                 23 Oct 2001 (RRD) - fix format error 
!                 05 Sep 2002 (RRD) - error trap
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 02 Apr 2004 (RRD) - generic file unit definitions
!                 12 Jul 2004 (RRD) - free format file structure
!
! USAGE:  CALL EMSINI(NLOC,OLAT,OLON,QFILE,NQLON,NQLAT,NPVAL)
!
!   INPUT ARGUMENT LIST:       see below
!   OUTPUT ARGUMENT LIST:      see below
!   INPUT FILES:               unit KF31 defines emission input data file
!   OUTPUT FILES:              unit KF21 for diagnostic messages to MESSAGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE EMSINI(NLOC,OLAT,OLON,QFILE,NQLON,NQLAT,NPVAL)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)    :: nloc      ! total number of source locations
  REAL,      INTENT(IN)    :: olat(:)   ! source subgrid domain limits 
  REAL,      INTENT(IN)    :: olon(:)   ! source subgrid domain limits 
  LOGICAL,   INTENT(OUT)   :: qfile     ! indicates area source emission file
  INTEGER,   INTENT(OUT)   :: nqlon     ! number of long in subgrid
  INTEGER,   INTENT(OUT)   :: nqlat     ! number of lats in subgrid
  INTEGER,   INTENT(OUT)   :: npval     ! number of pollutants in file

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  CHARACTER(80)            :: fname
  REAL                     :: units, dlat, dlon
  INTEGER                  :: nqval, kret

!-------------------------------------------------------------------------------
! input the emissions from the inventory file
! flag set to point source emission if file is not found
!-------------------------------------------------------------------------------

  FNAME='emission.txt'
  QFILE=.FALSE.

  INQUIRE(FILE=FNAME,EXIST=QFILE)
  IF(QFILE)THEN
     OPEN(KF31,FILE=FNAME,ACTION='READ')
  ELSE
     RETURN
  END IF

!-------------------------------------------------------------------------------
! check inputs 
!-------------------------------------------------------------------------------

  IF(NLOC.NE.2)THEN
     WRITE(KF21,*)'*ERROR* emsini: 2 source locations required to define grid'
     WRITE(*,*)   '*ERROR* emsini: see message file for more information'     
     STOP 900
  END IF

!-------------------------------------------------------------------------------
! set emission file default values from index record
!-------------------------------------------------------------------------------
! NPVAL - the number of pollutants in file
! NQVAL - number of emission values per 24 hour period 
!         can be on multiple lines (limit 12 per line)
! UNITS - defines conversion from file units to internal units/hour
! DLAT  - selected resolution for input file summation, which should
! DLON    be set to some a value equal to or coarser than file
!-------------------------------------------------------------------------------

! READ(KF31,'(2I4,3F10.4)',IOSTAT=KRET) NPVAL,NQVAL,UNITS,DLAT,DLON
  READ(KF31,*,IOSTAT=KRET)              NPVAL,NQVAL,UNITS,DLAT,DLON

  IF(KRET.NE.0)THEN
     WRITE(KF21,*)'*ERROR* emsini: header information not found in emission.txt'
     WRITE(*,*)   '*ERROR* emsini: see message file for more information'     
     STOP 900
  END IF
  REWIND (KF31)

! number of grid points requested
  NQLAT=INT((OLAT(2)-OLAT(1))/DLAT)+1
  NQLON=INT((OLON(2)-OLON(1))/DLON)+1
  WRITE(KF21,*)' NOTICE emsini: internal emission grid defined - ', &
                 NPVAL,NQVAL,UNITS,DLAT,DLON,NQLAT,NQLON

END SUBROUTINE emsini
