!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONAVG           LAYER AVERAGE THE CONCENTRATION ARRAY
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!    vertically average each i,j element of array conc over element
!    k and place results in array temp.  all k elements with heights
!    between the range level1 and level2 are averaged by the layer
!    depth.  The routine is called by most plotting routines.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD) 
!                 20 Nov 2000 (RRD) - fortran90 upgrade 
!                 19 Dec 2001 (RRD) - vertical averaging options
!
! USAGE: CALL CONAVG (CONC,TEMP,HEIGHT,LEVEL1,LEVEL2,NLAT,NLON,NLVL,KMAX)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CONAVG(CONC,TEMP,HEIGHT,LEVEL1,LEVEL2,NLAT,NLON,NLVL,KMAX)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,    INTENT(IN)   :: conc (:,:,:)   ! three dimensional input array
  INTEGER, INTENT(IN)   :: height (:)     ! height of with each vertical index
  INTEGER, INTENT(IN)   :: level1         ! bottom depth for averaging
  INTEGER, INTENT(IN)   :: level2         ! top depth for averaging
  INTEGER, INTENT(IN)   :: nlat           ! number of j values to process
  INTEGER, INTENT(IN)   :: nlon           ! number of i values to process
  INTEGER, INTENT(IN)   :: nlvl           ! number of k values to process
  INTEGER, INTENT(IN)   :: kmax           ! vertical averaging options 

  REAL,    INTENT(OUT)  ::  temp (:,:)    ! two dimensional output array

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER :: kk,jj,ii,last,delz
  REAL    :: wght,vsum

!-------------------------------------------------------------------------------

  DO II=1,NLON
  DO JJ=1,NLAT

     TEMP(II,JJ)=0.0
     VSUM=0.0
     WGHT=0.0
     LAST=0

     DO KK=1,NLVL
!       only vertically average non-deposition
        IF(HEIGHT(KK).GT.0.AND.HEIGHT(KK).GE.LEVEL1.AND.                   &
           HEIGHT(KK).LE.LEVEL2)THEN

           IF(KMAX.EQ.0)THEN
!             kmax=0 standard layer weighted vertical average
              DELZ=HEIGHT(KK)-LAST
              VSUM=VSUM+CONC(II,JJ,KK)*DELZ
              WGHT=WGHT+DELZ
              LAST=HEIGHT(KK)
           ELSE
!             kmax=1 maximum value in layer
              VSUM=MAX(VSUM,CONC(II,JJ,KK))
              WGHT=1.0
           END IF

        END IF
     END DO

!    save average
     IF(WGHT.GT.0.0)TEMP(II,JJ)=VSUM/WGHT

  END DO
  END DO

END SUBROUTINE conavg
