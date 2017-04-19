!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFMIN           determines PUFf MINimum mass
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF MINIMUM DETERMINES MINIMUM PUFF/PARTICLE MASS TO ELIMINATE
!   SMALL MASS PARTICLES THROUGH BINARY SORT.  RETURNS THE PARTICLE
!   MASS   BELOW WHICH ONLY CONTAINS FRACM OF THE TOTAL MASS
!   NEEDS TO BE CALLED IMMEDIATELY AFTER PUFDEL
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL PUFMIN(KPM,MASS,NSORT,TMASS,FRAC1,FMASS1,FRAC2,FMASS2)
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

SUBROUTINE PUFMIN(KPM,MASS,NSORT,TMASS,FRAC1,FMASS1,FRAC2,FMASS2)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER,  INTENT(IN)    :: kpm         ! total number of puffs or particles
  REAL,     INTENT(IN)    :: mass (:,:)  ! particle mass array
  INTEGER,  INTENT(INOUT) :: nsort (:)   ! sortted array by position
  REAL,     INTENT(IN)    :: tmass       ! total system mass
  REAL,     INTENT(IN)    :: frac1       ! fraction of total mass tmass
  REAL,     INTENT(IN)    :: frac2       ! fraction of total mass tmass
  REAL,     INTENT(OUT)   :: fmass1      ! single particle mass at frac1*tmass
  REAL,     INTENT(OUT)   :: fmass2      ! single particle mass at frac2*tmass

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  INTEGER,  ALLOCATABLE   :: var   (:)   ! sortted array by position

  REAL      :: tmass1,tmass2,fmass,xmass
  INTEGER   :: kret,numtyp,index0,index2,index1 
  INTEGER   :: mm,i,kpd,ig1,ig2,kp,kn,len1,l2,len2,l1 

!-------------------------------------------------------------------------------

! need at least two elements to sort
  IF(KPM.LT.2)RETURN

! number of pollutants is the first dimension of mass
  NUMTYP=SIZE(mass,1)

! allocate temporary sort index array
  ALLOCATE (var(kpm),STAT=kret)
  IF(kret.NE.0)THEN
     WRITE(*,*)'*ERROR* pufmin: memory allocation' 
     STOP 900
  END IF

! length of primary sort group
  KP=1

! there are 2^n passes through n data points
  DO WHILE (KP.LT.KPM)
!    start index of sortted output
     KN=1

!    pass through the data by group
     DO WHILE (KN.LE.KPM)
!       first group starting pointer
        IG1=KN
!       second group starting pointer
        IG2=KN+KP

!       number of elements in group being tested
        LEN1=MAX(MIN(KPM-IG1+1,KP),0)
        LEN2=MAX(MIN(KPM-IG2+1,KP),0)

!       pointer to elements in each group
        L1=0
        L2=0

!       loop through the elements within a group
        DO WHILE (L1.LT.LEN1.AND.L2.LT.LEN2)
!          convert pointer to array index number
           INDEX1=NSORT(IG1+L1)
           INDEX2=NSORT(IG2+L2)

           TMASS1=MASS(1,INDEX1)
           TMASS2=MASS(1,INDEX2)
           MM=NUMTYP
           DO WHILE(MM.GT.1)
              TMASS1=TMASS1+MASS(MM,INDEX1)
              TMASS2=TMASS2+MASS(MM,INDEX2)
              MM=MM-1
           END DO

           IF(TMASS2.LT.TMASS1)THEN
              VAR(KN)=NSORT(IG2+L2)
              L2=L2+1
              KN=KN+1
           ELSE
              VAR(KN)=NSORT(IG1+L1)
              L1=L1+1
              KN=KN+1
           END IF
        END DO

!       when one group is finished ... copy the other
        DO WHILE (L1.LT.LEN1.AND.L2.EQ.LEN2)
           VAR(KN)=NSORT(IG1+L1)
           L1=L1+1
           KN=KN+1
        END DO

        DO WHILE (L2.LT.LEN2.AND.L1.EQ.LEN1)
           VAR(KN)=NSORT(IG2+L2)
           L2=L2+1
           KN=KN+1
        END DO
     END DO

!    save sortted index structure for next pass
     DO I=1,KPM
        NSORT(I)=VAR(I)
     END DO

!    double size of sort group
     KP=KP*2
  END DO

! sum mass until equals fraction-1 of total
  FMASS=0.0
  XMASS=0.0
  KPD=0
  DO WHILE (XMASS.LE.FRAC1*TMASS.AND.KPD.LT.KPM)
     FMASS1=FMASS
     KPD=KPD+1
     INDEX0=NSORT(KPD)
     FMASS=MASS(1,INDEX0)
     MM=NUMTYP
     DO WHILE(MM.GT.1)
        FMASS=FMASS+MASS(MM,INDEX0)
        MM=MM-1
     END DO
     XMASS=XMASS+FMASS
  END DO

! sum mass until equals fraction-2 of total
  FMASS=0.0
  XMASS=0.0
  KPD=0
  DO WHILE (XMASS.LE.FRAC2*TMASS.AND.KPD.LT.KPM)
     FMASS2=FMASS
     KPD=KPD+1
     INDEX0=NSORT(KPD)
     FMASS=MASS(1,INDEX0)
     MM=NUMTYP
     DO WHILE(MM.GT.1)
        FMASS=FMASS+MASS(MM,INDEX0)
        MM=MM-1
     END DO
     XMASS=XMASS+FMASS
  END DO

  DEALLOCATE (var)

END SUBROUTINE pufmin
