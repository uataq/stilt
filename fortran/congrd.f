!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONGRD           SET CONCENTRATION GRID ARRAY VALUES  
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            READS THE HEADER INFORMATION IN THE CONCENTRATION FILE
!            AND SETS RELATED PROGRAM VARIABLES TO DETERMINE THE ARRAY
!            SIZES FOR DYNAMIC MEMORY ALLOCATION
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 20 Nov 2000 (RRD) - initial version
!                 22 Dec 2004 (RRD) - cpack test
!
! USAGE: CALL CONGRD(NLOC,NLAT,NLON,NTYP,NLVL,CPACK)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000    
!
!$$$

SUBROUTINE CONGRD(NLOC,NLAT,NLON,NTYP,NLVL,CPACK)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,     INTENT(OUT)   :: nloc        ! actual numb of source locations
  INTEGER,     INTENT(OUT)   :: nlat,nlon   ! number of actual pnts in grid
  INTEGER,     INTENT(OUT)   :: ntyp        ! number of pollutants in file
  INTEGER,     INTENT(OUT)   :: nlvl        ! number of concen levels in file
  INTEGER,     INTENT(OUT)   :: cpack       ! concentration data packing method

!-------------------------------------------------------------------------------
! internally defined variables
!-------------------------------------------------------------------------------

  CHARACTER(4) :: PCHAR          
  INTEGER      :: KK 
  INTEGER      :: IDATA(10)
  REAL         :: RDATA(10)

!-------------------------------------------------------------------------------

! meteo file type and initial position date
  READ(10,IOSTAT=KK)PCHAR,IDATA(1:5),NLOC,CPACK
  IF(KK.NE.0)CPACK=0  ! old style file with all grid points

  IF(CPACK.GT.1)THEN
     WRITE(*,*)'*ERROR* coninp: unsupported packing option -',CPACK
     STOP 900
  END IF

! release starting information
  DO KK=1,NLOC
     READ(10)  IDATA(1:4),RDATA(1:3)
  END DO

! number of points
  READ(10)   NLAT, NLON, RDATA(1:4)                

! vertical grid index record
  READ(10)   NLVL, (IDATA(1),KK=1,NLVL)

! pollutant identification record
  READ(10)   NTYP, (PCHAR,KK=1,NTYP)

! start over with next pass
  REWIND(10)

END SUBROUTINE congrd
