!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSDAT           EMiSsion DATa read of special file
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION INPUT - READS A LAT/LON BASED EMISSION INVENTORY FILE.
!   THE EMISSION ARRAY BASED UPON A TWO-POINT SOURCE LIMIT DEFINED IN
!   MODEL INPUT SECTION. THE ROUTINE IS CALLED ONLY ONCE.
!   THE FILE IS ONE RECORD PER EMISSION POINT WITH EACH EMISSION POINT
!   DEFINED AT A LAT/LON POINT. EMISSIONS ARE SUMMED INTO A GRID.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Feb 1998 (RRD)
!                 22 Dec 1998 (RRD) - generalized subroutine
!                 03 Sep 2000 (RRD) - fortran90 upgrade
!                 23 Oct 2001 (RRD) - format error correction
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 22 Jul 2003 (RRD) - independent emissions file
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 12 Jul 2004 (RRD) - space delimited formats / accumulation
!
! USAGE:  CALL EMSDAT(QDLON,QDLAT,NQLON,NQLAT,OLON,OLAT,QAREA,POLID)
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

SUBROUTINE EMSDAT(QDLON,QDLAT,NQLON,NQLAT,OLON,OLAT,QAREA,POLID)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,        INTENT(OUT) :: qdlon, qdlat             ! grid spacing
  INTEGER,     INTENT(IN)  :: nqlon, nqlat             ! number of points
  REAL,        INTENT(IN)  :: olon(:), olat(:)         ! sub grid corners
  REAL,        INTENT(OUT) :: qarea ( : , : , : , : )  ! houly emissions
  CHARACTER(4),INTENT(OUT) :: polid ( : )              ! pollutant ids 

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL                     :: qval(24)         ! max daily input is hourly
  REAL,   ALLOCATABLE      :: qsum(:)          ! diagnostic variable
  LOGICAL                  :: qfile,skip       ! record test
  INTEGER                  :: n,i,j,k,l,ii,jj,nn,ih
  INTEGER                  :: kret,kpts,nqval,npval,nqhrs
  REAL                     :: units,xlat,xlon
  CHARACTER(80)            :: fname

!-------------------------------------------------------------------------------

! first index record (repeated from emsini)
! READ(KF31,'(2I4,3F10.4)') NPVAL,NQVAL,UNITS,QDLAT,QDLON 
  READ(KF31,*)              NPVAL,NQVAL,UNITS,QDLAT,QDLON 

! diagnostic variables
  ALLOCATE (qsum(npval), STAT=kret)
  IF(kret.NE.0)THEN
     WRITE(*,*)'*ERROR* emsdat: memory allocation' 
     STOP 900
  END IF
  QAREA = 0.0
  QSUM  = 0.0

! hours per q value
  NQHRS=24/NQVAL
  IF(NQVAL.GT.24)THEN
     WRITE(*,*)'*ERROR* emsdat: more than 24 emissions defined - ',NQVAL
     STOP 900
  END IF

! second index record with polluant id field (should match control file)
  READ(KF31,'(255A4)')(POLID(K),K=1,NPVAL) 

! input the emissions from the inventory file
  READ(KF31,'(A)')FNAME
  CLOSE(KF31)
  WRITE(KF21,*)' NOTICE emsdat: emission file header - ',      &
                 NPVAL,NQVAL,UNITS,QDLAT,QDLON

!-------------------------------------------------------------------------------
! open the actual emissions file 
!-------------------------------------------------------------------------------

  QFILE=.FALSE.
  INQUIRE(FILE=FNAME,EXIST=QFILE)

  IF(QFILE)THEN
     OPEN(KF31,FILE=FNAME,ACTION='READ')
  ELSE
     WRITE(KF21,*)'*ERROR* emsdat: gridded emissions file not found - ', FNAME
     WRITE(*,*)   '*ERROR* emsdat: see message file for more information'      
     STOP 900
  END IF

!-------------------------------------------------------------------------------
! loop thru entire file
!-------------------------------------------------------------------------------

  inloop: DO

! READ(KF31,'(2I4,2F10.4)',END=900)II,JJ,XLON,XLAT
  READ(KF31,*,END=900)II,JJ,XLON,XLAT

! check if record is within selected limits
  SKIP=.FALSE.
  IF(XLAT.LT.OLAT(1).OR.XLAT+QDLAT.GT.OLAT(2))SKIP=.TRUE.
  IF(XLON.LT.OLON(1).OR.XLON+QDLON.GT.OLON(2))SKIP=.TRUE.

! compute index on internal emission grid
! when SW corner falls within cell then accumulate in that cell
  II=INT((XLON-OLON(1))/QDLON)+1
  JJ=INT((XLAT-OLAT(1))/QDLAT)+1

! check limits again
  IF(II.LT.1.OR.II.GT.NQLON.OR.JJ.LT.1.OR.JJ.GT.NQLAT)SKIP=.TRUE.

! Assume that emissions defined as units/interval regardless of how many
! values define the diurnal cycle. For instance, 8 input values per
! day would define the emission for a period of three hours. These
! values are then converted to an hourly rate.

  DO NN=1,NPVAL
!    READ(KF31,'(12E10.3)') QVAL(1:NQVAL)
     READ(KF31,*)           QVAL(1:NQVAL)

     IF(.NOT.SKIP)THEN
        DO IH=1,24
!          convert hour to input file index value
           K=1+(IH-1)/NQHRS
           
!          convert emissions to units/hr
           QAREA(II,JJ,NN,IH)=QAREA(II,JJ,NN,IH)+QVAL(K)*UNITS/NQHRS

!          sum average hourly rate for diagnostics
           QSUM(NN)=QSUM(NN)+QVAL(K)*UNITS/NQHRS
        END DO
     END IF
  END DO

  END DO inloop

!-------------------------------------------------------------------------------
! diagnostic dump of generated emissions
!-------------------------------------------------------------------------------

900 KPTS=0

  DO J=1,NQLAT
  DO I=1,NQLON
     SKIP=.FALSE.
     DO L=1,24
     DO K=1,NPVAL
!       set flag if any hour or pollutant is non-zero
        IF(QAREA(I,J,K,L).GT.0.0)SKIP=.TRUE.
     END DO
     END DO
!    count the number of non-zero grid points
     IF(SKIP)KPTS=KPTS+1
  END DO
  END DO

!-------------------------------------------------------------------------------
! terminate and print message according to status
!-------------------------------------------------------------------------------

  IF(KPTS.EQ.0)THEN
     WRITE(KF21,*)'*ERROR* emsdat: lat/lon emission selection'
     WRITE(KF21,*)' There are no emissions in the region'
     WRITE(*,*)   '*ERROR* emsdat: see message file for more information'
     STOP 900
  END IF

  QSUM = QSUM / NQVAL
  WRITE(KF21,*)' NOTICE emsdat: gridded emissions initialized'
  WRITE(KF21,'(10X,5A10)')(POLID(N),N=1,NPVAL)
  WRITE(KF21,'(A10,5E10.2)')' Kg/hr ',(QSUM(N),N=1,NPVAL)
  WRITE(KF21,*)'Number of emission lat/lon points: ',NQLAT,NQLON
  WRITE(KF21,*)'Total number of non-zero points  : ',KPTS
  WRITE(KF21,*)' '

  CLOSE (KF31)
  DEALLOCATE (qsum)

END SUBROUTINE emsdat
