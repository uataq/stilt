!###############################################################################
! GEMDAY - Write concentrations to disk in binary format for initialization
!-------------------------------------------------------------------------------
! LAST REVISED: 22 May 2008 (RRD) - initial version
!               02 Jul 2008 (RRD) - eliminate normalization
!-------------------------------------------------------------------------------

SUBROUTINE gemday (IYR,IMO,IDA,IHR,IMN)

  USE gemkon  
  USE gemvar
  USE gemcfg
  USE funits

  IMPLICIT NONE

  INTEGER*4, INTENT(IN) :: iyr,imo,ida,ihr,imn

  INTEGER*4     :: nx,ny,nz,np,kgrd
  REAL*4        :: clat1,clon1,dlat,dlon

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  IF(KINIT.LT.0)RETURN

!-------------------------------------------------
! daily output of concentration field to dump file
  IF(IHR.NE.0) RETURN

  WRITE(KF28) iyr,imo,ida,ihr,imn
  WRITE(KF28) nx,ny,nz,np
  WRITE(KF28) clat1,clon1,dlat,dlon
  WRITE(KF28) sig
  WRITE(KF28) DDD
  WRITE(KF28) XXX
  REWIND(KF28)

END SUBROUTINE gemday
