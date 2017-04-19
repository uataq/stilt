!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKOUT           PAcK OUTput converts real to character
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACK OUTPUT CONVERTS A REAL ARRAY TO CHARACTER*1 PACKED ARRAY
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 14 Feb 1997 (RRD)
!                  02 Feb 2001 (RRD) - fortran90 upgrade
!                  09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL PAKOUT(RVAR,CVAR,NX,NY,NXY,PREC,NEXP,VAR1,KSUM)
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

SUBROUTINE PAKOUT(RVAR,CVAR,NX,NY,NXY,PREC,NEXP,VAR1,KSUM)

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)  :: nx,ny,nxy   ! dimension limits
  REAL,      INTENT(IN)  :: rvar(nx,ny) ! data array to be packed  
  CHARACTER, INTENT(OUT) :: cvar(nxy)   ! packed char*1 output array
  REAL,      INTENT(OUT) :: prec        ! precision of packed data array
  INTEGER,   INTENT(OUT) :: nexp        ! packing scaling exponent
  REAL,      INTENT(OUT) :: var1        ! value of real array at position (1,1)
  INTEGER,   INTENT(OUT) :: ksum        ! rotating checksum of packed data 

!-------------------------------------------------------------------------------

  INTEGER :: icval,i,j,k 
  REAL    :: scexp,rcol,sexp,rmax,rold 

!-------------------------------------------------------------------------------

  VAR1=RVAR(1,1)

  ROLD= VAR1
  RMAX= 0.0
! find the maximum difference between adjacent elements
  DO J=1,NY
     DO I=1,NX
!       compute max difference between elements along row
        RMAX=MAX( ABS(RVAR(I,J)-ROLD), RMAX)
        ROLD=RVAR(I,J)
     END DO
!    row element 1 difference always from previous row
     ROLD=RVAR(1,J)
  END DO

  SEXP=0.0
! compute the required scaling exponent
  IF(RMAX.NE.0.0)SEXP=LOG(RMAX)/LOG(2.)
  NEXP=INT(SEXP)
! positive or whole number scaling round up for lower precision
  IF(SEXP.GE.0.0.OR.MOD(SEXP,1.0).EQ.0.0)NEXP=NEXP+1
! precision range is -127 to 127 or 254
  PREC=(2.0**NEXP)/254.0
  SCEXP=2.0**(7-NEXP)

! initialize checksum
  KSUM=0
! set column1 value
  RCOL=VAR1

  K=0
! pack the array from low to high
  DO J=1,NY
     ROLD=RCOL
     DO I=1,NX
        K=K+1

!       packed integer at element
        ICVAL=INT((RVAR(I,J)-ROLD)*SCEXP+127.5)
!       previous element as it would appear unpacked
        ROLD=FLOAT(ICVAL-127)/SCEXP+ROLD
!       save the first column element for next row
        IF(I.EQ.1)RCOL=ROLD
!       convert to character
        CVAR(K)=CHAR(ICVAL)

!       maintain rotating checksum
        KSUM=KSUM+ICVAL
!       if sum carries over the eighth bit add one
        IF(KSUM.GE.256)KSUM=KSUM-255

     END DO
  END DO

END SUBROUTINE pakout
