!###############################################################################
! GEMSET - Allocate the meteorological data arrays
!-------------------------------------------------------------------------------
! UUU - U WIND COMPONENTS       m/s   
! VVV - V WIND COMPONENTS       m/s
! TTT - TEMPERATURE             Kelvin
! MMM - MOISTURE (RH)           percent
! HHH - HEIGHT MSL OF FIELD     meters
! WWW - DIVERGENCE & VELOCITY   1/s to m/s
! KKK - VERTICAL MIXING         m2/s
! RRR - LOCAL AIR DENSITY       Kg/m3
! XXX - MIXING RATIO            g/Kg
! CCC - TEMORARY MIXING RATIO   g/Kg
! ZZZ - TEMORARY PROFILE        g/Kg
! GSX - WE GRID DIMENSIONS      meters
! GSY - SN GRID DIMENSIONS      meters
! SFC - SURFACE PRESSURE        mb
! LUS - LANDUSE CATEGORY        index
! RGH - ROUGHNESS LENGTH        m
! FRV - FRICTION VELOCITY       m/s
! PSI - INTEGRATED STABILITY
! SWF - SHORT WAVE FLIX         w/m2
! U10 - 10 meter wind           m/s
! V10 - 10 meter wind           m/s
! T02 - 2  meter temperature    deg-K
! P6H - precipitation rate      m/min
!-------------------------------------------------------------------------------
! LAST REVISED: 30 May 2008 (RRD) - initial version
!-------------------------------------------------------------------------------
 
SUBROUTINE gemset

  USE funits
  USE gemkon
  USE gemvar 

  IMPLICIT NONE

  INTEGER*4 :: nx,ny,nz,np,kret,ktot,kgrd
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  ktot=0

  ALLOCATE (NVAR(0:nz), NGP(ny), KNDX(nz), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #01 - ',kret
  ktot=ktot+kret

  ALLOCATE (PPP(0:nz),SIG(nz),VAL(nz,7), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #02 - ',kret
  ktot=ktot+kret

  ALLOCATE (POT(nz),VB4(nz),VB8(nz), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #03 - ',kret
  ktot=ktot+kret

  ALLOCATE (UUU(nx,ny,nz),VVV(nx,ny,nz),WWW(nx,ny,nz), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #04 - ',kret
  ktot=ktot+kret

  ALLOCATE (TTT(nx,ny,nz),MMM(nx,ny,nz),HHH(nx,ny,nz), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #05 - ',kret
  ktot=ktot+kret

  ALLOCATE (KKK(nx,ny,nz),RRR(nx,ny,nz), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #06 - ',kret
  ktot=ktot+kret

  ALLOCATE (GSX(nx,ny),GSY(nx,ny),GXY(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #07 - ',kret
  ktot=ktot+kret

  ALLOCATE (SFC(nx,ny),AVG(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #08 - ',kret
  ktot=ktot+kret

  ALLOCATE (XXX(nx,ny,nz,np),CCC(nx,ny,nz),ZZZ(nz), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #09 - ',kret
  ktot=ktot+kret

  ALLOCATE (DDD(nx,ny,np),CONC(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #10 - ',kret
  ktot=ktot+kret

  ALLOCATE (LUS(nx,ny),RGH(nx,ny),FRV(nx,ny),PSI(nx,ny),SWF(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #11 - ',kret
  ktot=ktot+kret

  ALLOCATE (U10(nx,ny),V10(nx,ny),T02(nx,ny),P6H(nx,ny), STAT=kret)
  IF(kret.NE.0)WRITE(KF21,'(A,I3)')'*ERROR gemset: memory allocation #12 - ',kret
  ktot=ktot+kret

  IF(ktot.gt.0) THEN
     WRITE(KF21,*)'*ERROR gemset: unable to allocate array space'
     STOP 800
  ELSE
     WRITE(KF21,*)' NOTICE gemset: allocated array space (x,y,z) - ',nx,ny,nz
  END IF

! initialize diagnostic variable arrays (may not always be present at all levels)
  www = 0.0
  mmm = 0.0
  p6h = 0.0
  swf = 0.0
  rgh = 0.2
  lus = 2

END SUBROUTINE gemset 
