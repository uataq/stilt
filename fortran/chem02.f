!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CHEM02           CHEMistry version 02 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-05-15
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   VERSION 2 OF THE CHEMISTRY TRANSFORMATION SUBROUTINE IS USED TO
!   CONVERT ONE CHEMICAL SPECIES TO ANOTHER WHEN BOTH SPECIES ARE
!   DEFINED ON THE SAME PARTICLE.  
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 15 May 2001 (RRD) - Initial version
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 18 Sep 2002 (RRD) - more generic conversion
!                 01 Apr 2004 (RRD) - generic file unit numbers
!                 17 Oct 2007 (RRD) - force positive time step
!
! USAGE:  CALL CHEM02(DIRT,DT,MASS)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             see module funits
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CHEM02(DIRT,DT,MASS)

  USE funits

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'     ! concentration and pollutant structure

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(pset), INTENT(IN)    :: dirt (:)    ! for each pollutant type 
  REAL,       INTENT(IN)    :: dt          ! time step (min)
  REAL,       INTENT(INOUT) :: mass (:)    ! particle mass array

!-------------------------------------------------------------------------------

  INTEGER,   ALLOCATABLE    :: K1(:),K2(:)  
  REAL,      ALLOCATABLE    :: RATE(:),RGMW(:)

  LOGICAL                   :: ftest
  INTEGER                   :: kt,maxdim,nrate 
  REAL                      :: total,beta

!-------------------------------------------------------------------------------

  DATA NRATE/0/
  SAVE NRATE,K1,K2,RATE,RGMW

!-------------------------------------------------------------------------------
 
! initialize the conversion routine
  IF(NRATE.EQ.0)THEN

     MAXDIM = SIZE(mass,1)  ! number of pollutants on single particle
!    routine only for single particle multiple species
     IF(MAXDIM.LE.1)THEN
        WRITE(KF21,*)'*ERROR* chem02 - multiple species definition required'
        WRITE(KF21,*)' Redefine with MAXDIM = ',MAXDIM
        WRITE(*,*)   '*ERROR* chem02 - see message file for more information'
        STOP
     END IF

     INQUIRE(FILE='CHEMRATE.TXT',EXIST=FTEST)
     IF(FTEST)THEN
        NRATE=0
        OPEN(KF33,FILE='CHEMRATE.TXT')
        loop1 : DO
           READ(KF33,*,IOSTAT=KT)
           IF(KT.NE.0) EXIT loop1
           NRATE=NRATE+1
        END DO loop1
        REWIND (KF33)
      
        ALLOCATE (K1(NRATE),K2(NRATE))        
        ALLOCATE (RATE(NRATE),RGMW(NRATE))        

        DO KT=1,NRATE
           READ(KF33,*)K1(KT),K2(KT),RATE(KT),RGMW(KT)
        END DO
        CLOSE (KF33)

     ELSE
        WRITE(KF21,*)' NOTICE chem02: CHEMRATE.TXT not found, using default'
        NRATE=1
        ALLOCATE (K1(NRATE),K2(NRATE))        
        ALLOCATE (RATE(NRATE),RGMW(NRATE))        
        K1(1)=1
        K2(1)=2
        RATE(1)=0.10
!       ratio of conversion molecular weights
!       for instance SO2/SO4 = 1.5
        RGMW(1)=1.0  
     END IF

     DO KT=1,NRATE
        WRITE(KF21,*)' NOTICE chem02: ',DIRT(K1(KT))%IDENT,' to ',  &
                                        DIRT(K2(KT))%IDENT
     END DO
  END IF

  DO KT=1,NRATE

!    change rate to per time step
     BETA=ABS(DT*RATE(KT)/60.0)

!    compute tranformation amount with the restriction that 
!    can not convert more than exists

     IF(BETA.LT.0.01)THEN
!       small values assume linear approximation
        TOTAL=MIN(MASS(K1(KT)),MASS(K1(KT))*BETA)
        MASS(K1(KT))=MASS(K1(KT))-TOTAL 

!       apply reduction to increase mass of the other species
!       according to the molecular weight ratio
        MASS(K2(KT))=MASS(K2(KT))+TOTAL*RGMW(KT) 

     ELSEIF(BETA.NE.0.0)THEN
!       apply exponential for removal > 1%
        TOTAL=MIN(MASS(K1(KT)),MASS(K1(KT))*(1.0-EXP(-BETA)))
        MASS(K1(KT))=MASS(K1(KT))-TOTAL 
        MASS(K2(KT))=MASS(K2(KT))+TOTAL*RGMW(KT) 

     ELSE 
!       no conversion
        CONTINUE    
     END IF

  END DO

END SUBROUTINE chem02
