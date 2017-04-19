!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CHMPHG           Computes photolysis rate coefficient for gear method 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            uses the solar angle and surface short wave flux to compute
!            the NO2 photolysis rate coefficient.  The rate coefficient
!            is used in the integration of VOC for smog produced
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 15 Jan 1997 (RRD) - initial version
!                 24 Mar 2005 (AFS) - upgraded for cb4 chemistry
!
! USAGE:  CALL CHMPHG(DIRT,NUMTYP,EA,SWF,RK1,RKJ)
!   INPUT ARGUMENT LIST:  see below
!   OUTPUT ARGUMENT LIST: see below
!   INPUT FILES:   NONE
!   OUTPUT FILES:  NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CHMPHG(DIRT,NUMTYP,EA,SWF,RK1,RKJ)

  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
      
!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(pset),INTENT(IN)    :: dirt(:)       ! for each pollutant type 
  INTEGER,   INTENT(IN)    :: numtyp        ! number of pollutant types
  REAL,      INTENT(IN)    :: ea            ! elevation angle
  REAL,      INTENT(IN)    :: swf           ! incident short wave flux (w/m2)
  REAL,      INTENT(OUT)   :: rk1           ! NO2 photolysis rate constant
  REAL*8,    INTENT(INOUT) :: rkj(:,:)

!-------------------------------------------------------------------------------
! internal variable definitions
!-------------------------------------------------------------------------------
  real                       :: zag,zar
  integer                    :: ktime,nspecies,k,l,kz,kc

  REAL,        allocatable   :: rkjphot(:,:)
  REAL,        allocatable   :: rkjph(:)
  INTEGER,     allocatable   :: nflag(:)
  CHARACTER*4, allocatable   :: especie(:)  

  REAL,        PARAMETER     :: rpdg  = 0.017453292

  DATA KTIME/0/
  SAVE NSPECIES,KTIME,RKJPHOT,RKJPH,NFLAG,ESPECIE

  if(.not.allocated(rkjphot))then
     ALLOCATE(rkjphot(NUMTYP,10))
     ALLOCATE(rkjph(NUMTYP)) 
     ALLOCATE(nflag(NUMTYP))
     ALLOCATE(especie(NUMTYP))
  end if 
	
  IF(KTIME.EQ.0)THEN 
     OPEN(80,FILE='PHOTO.DAT',STATUS='OLD')
     READ(80,*) NSPECIES
     DO K=1,NSPECIES
        READ(80,*) ESPECIE(K),NFLAG(K),(RKJPHOT(K,KZ),KZ=1,10)
     END DO

     CLOSE(80)
     KTIME=1
  END IF

!=>zenith angle in radians from elevation angle

   ZAG=90.0-EA
   ZAR=ZAG*RPDG

!=>rate coefficient for NO2 photolysis

   IF(ZAG.GE.0.0.AND.ZAG.LT.47.0)THEN
      RK1=(4.23E-04+(1.09E-04/COS(ZAR)))*SWF    
   ELSEIF(ZAG.GE.47.0.AND.ZAG.LT.64.0)THEN
      RK1=5.82E-04*SWF    
   ELSEIF(ZAG.GE.64.0.AND.ZAG.LE.90.0)THEN
      RK1=(-0.997E-04+1.2E-03*(1.0-COS(ZAR)))*SWF    
   ELSE
      RK1=0.0
   END IF
	
   DO L=1,NSPECIES
      IF(ZAG.GE.0.0.AND.ZAG.LT.10.0)THEN
         RKJPH(L)=RKJPHOT(L,1)      
      ELSEIF(ZAG.GE.10.0.AND.ZAG.LT.20.0)THEN
         RKJPH(L)=RKJPHOT(L,2)
      ELSEIF(ZAG.GE.20.0.AND.ZAG.LT.30.0)THEN
         RKJPH(L)=RKJPHOT(L,3)
      ELSEIF(ZAG.GE.30.0.AND.ZAG.LT.40.0)THEN
         RKJPH(L)=RKJPHOT(L,4)
      ELSEIF(ZAG.GE.40.0.AND.ZAG.LT.50.0)THEN
         RKJPH(L)=RKJPHOT(L,5)
      ELSEIF(ZAG.GE.50.0.AND.ZAG.LT.60.0)THEN
         RKJPH(L)=RKJPHOT(L,6)
      ELSEIF(ZAG.GE.60.0.AND.ZAG.LT.70.0)THEN
         RKJPH(L)=RKJPHOT(L,7)
      ELSEIF(ZAG.GE.70.0.AND.ZAG.LT.76.0)THEN
         RKJPH(L)=RKJPHOT(L,8)
      ELSEIF(ZAG.GE.76.0.AND.ZAG.LT.86.0)THEN
         RKJPH(L)=RKJPHOT(L,9)
      ELSEIF(ZAG.GE.86.0.AND.ZAG.LT.90.0)THEN
         RKJPH(L)=RKJPHOT(L,10)
      END IF

      DO KC=1,NUMTYP
         IF(ESPECIE(L).EQ.DIRT(KC)%IDENT)THEN
            IF(NFLAG(L).EQ.1) RKJ(KC,1)=DBLE(RKJPH(L))
            IF(NFLAG(L).EQ.2) RKJ(KC,2)=DBLE(RKJPH(L))
         END IF  
      END DO
   END DO

   RETURN
END SUBROUTINE CHMPHG 
