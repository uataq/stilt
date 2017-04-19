!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONSMT           SMOOTHS THE ELEMENTS OF AN ARRAY FOR PLOTTING
!   PRGMMR:    DALE HESS        ORG: BMRC/BoM    DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE BUREAU OF METEOROLOGY, MELBOURNE
!            AND IT SMOOTHS AN X-Y ARRAY ACCORDING TO SET WEIGHTING FACTORS
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD)
!                 21 Nov 2000 (RRD) - fortran90 upgrade
!                 05 Aug 2008 (RRD) - complete rewrite
!
! USAGE:  CALL CONSMT( RVAL,SCAN )
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

SUBROUTINE CONSMT( RVAL,DIST )

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,    INTENT(INOUT) :: rval (:,:)      ! array to be smoothed
  INTEGER, INTENT(IN)    :: dist            ! scan index for smoothing

  INTEGER :: i,j,k,ii,jj
  INTEGER :: imin,imax,jmin,jmax,scan,kret 

  REAL    :: tsum, wght(0:99)               ! scaling weights by grid point 
  REAL,      ALLOCATABLE :: work (:,:)      ! working array

!-------------------------------------------------------------------------------

  IF(dist.EQ.0)RETURN

  IMIN=1
  JMIN=1
  IMAX=SIZE(rval,1)
  JMAX=SIZE(rval,2)

  ALLOCATE (work(imax,jmax), STAT=kret)
  IF(kret.NE.0)THEN
     WRITE(*,*)'ERROR consmt allocation: ',imax,jmax,kret
     STOP 900
  END IF

! size of scan box for smoothing
  scan=MIN(99,dist,imax,jmax)

! reduce scaling weight by grid point
  wght(0)=2.0
  DO k=1,scan
     wght(k)=wght(k-1)/2.0   
  END DO

!-------------------------------------------------------------------------------

  DO J=JMIN,JMAX
  DO I=IMIN,IMAX
     tsum=wght(0)
     work(i,j)=wght(0)*rval(i,j) 
     DO k=1,scan

!       left edge
        IF(j-k.GE.jmin)THEN
        DO ii=MAX(imin,i-k),MIN(imax,i+k-1)
           work(i,j)=work(i,j)+wght(k)*rval(ii,j-k)  
           tsum=tsum+wght(k)
        END DO
        END IF

!       top edge
        IF(i+k.LE.imax)THEN
        DO jj=MAX(jmin,j-k),MIN(jmax,j+k-1)
           work(i,j)=work(i,j)+wght(k)*rval(i+k,jj)  
           tsum=tsum+wght(k)
        END DO
        END IF

!       right edge
        IF(j+k.LE.jmax)THEN
        DO ii=MIN(imax,i+k),MAX(imin,i-k+1),-1
           work(i,j)=work(i,j)+wght(k)*rval(ii,j+k)  
           tsum=tsum+wght(k)
        END DO
        END IF

!       bottom edge
        IF(i-k.GE.imin)THEN
        DO jj=MIN(jmax,j+k),MAX(jmin,j-k+1),-1
           tsum=tsum+wght(k)
           work(i,j)=work(i,j)+wght(k)*rval(i-k,jj) 
        END DO
        END IF

     END DO
     work(i,j)=work(i,j)/tsum
  END DO
  END DO
 
!-------------------------------------------------------------------------------
 
  rval=work  
  DEALLOCATE (work, STAT=kret)
  IF(kret.NE.0)THEN
     WRITE(*,*)'ERROR consmt deallocation: ',kret
     STOP 900
  END IF

END SUBROUTINE consmt
