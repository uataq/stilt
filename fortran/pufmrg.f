!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFMRG           PUFf MeRGes puffs together at same position
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF MERGE MERGES PUFFS TOGETHER WHEN AT SAME POSITION AND
!   OF THE SAME HORIZONATL DISTRIBUTION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 (RRD)
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!                 21 Oct 2001 (RRD) - real/integer conflict
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 15 Sep 2003 (RRD) - more focused test on hdwp
!                 12 Aug 2004 (RRD) - variable name change
!                 18 Jan 2008 (RRD) - absolute value page
!                 13 Nov 2008 (RRD) - enhanced test on hdwp
!
! USAGE:  CALL PUFMRG(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,MASS,XPOS,YPOS,ZPOS,  
!                     SIGH,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,HGKM)
!
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           none
!   OUTPUT FILES:          none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PUFMRG(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,MASS,XPOS,YPOS,ZPOS,  &
                  SIGH,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,HGKM)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,    INTENT(IN)    :: frach      ! horizontal position rounding fraction
  REAL,    INTENT(IN)    :: fracv      ! vertical position rounding fraction
  REAL,    INTENT(IN)    :: fract      ! travel-time (sigma) rounding fraction
  REAL,    INTENT(IN)    :: fmass      ! particle mass sort cutoff value
  INTEGER, INTENT(IN)    :: kpm        ! number of puffs or particles
  REAL,    INTENT(IN)    :: zmdl       ! model top scaling height
  REAL,    INTENT(INOUT) :: mass (:,:) ! pollutant mass
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horiz sigma (sigma)
  REAL,    INTENT(INOUT) :: sigw (:)   ! vert sigma (sigma)
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution for pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age from release 
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteo grid index for puff
  INTEGER, INTENT(IN)    :: nsort(:)   ! sortted array by position
  REAL,    INTENT(IN)    :: hgkm (:)   ! horizontal grid distance (km)

!-------------------------------------------------------------------------------

  REAL,    ALLOCATABLE   :: MASS0(:)   ! pollutant mass array (mass species)

  LOGICAL  :: sortf, equal
  INTEGER  :: kret,numtyp,mm,kp,kk,kn,index0,index1,index2,hdwpx
  REAL     :: hgd,page0,sigw0,sigh0,xpos0,ypos0,zpos0,tmass,tmass1,tmass2

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE PUFRND( ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,                         &
             XPOS1,YPOS1,ZPOS1,SIGH1,SIGW1,HDWP1,PTYP1,MASS1,TMASS1,           &
             XPOS2,YPOS2,ZPOS2,SIGH2,SIGW2,HDWP2,PTYP2,MASS2,TMASS2,SORTF,EQUAL)
  IMPLICIT NONE
  REAL,    INTENT(IN)  :: zmdl          ! model top scaling height
  REAL,    INTENT(IN)  :: frach         ! horizontal position rounding fraction
  REAL,    INTENT(IN)  :: fracv         ! vertical position rounding fraction
  REAL,    INTENT(IN)  :: fract         ! travel-time (sigma) rounding fraction
  REAL,    INTENT(IN)  :: fmass         ! particle mass sort cutoff value
  REAL,    INTENT(IN)  :: hgd           ! horizontal grid distance (m)
  REAL,    INTENT(IN)  :: xpos1, xpos2  ! puff center positions (grid units)
  REAL,    INTENT(IN)  :: ypos1, ypos2  ! puff center positions (grid units)
  REAL,    INTENT(IN)  :: zpos1, zpos2  ! puff center positions (grid units)
  REAL,    INTENT(IN)  :: sigh1, sigh2  ! horiz sigma (sigma)
  REAL,    INTENT(IN)  :: sigw1, sigw2  ! vert sigma (sigma)
  INTEGER, INTENT(IN)  :: hdwp1, hdwp2  ! distribution
  INTEGER, INTENT(IN)  :: ptyp1, ptyp2  ! pollutant type
  REAL,    INTENT(IN)  :: mass1  (:)    ! pollutant mass
  REAL,    INTENT(IN)  :: mass2  (:)    ! pollutant mass
  REAL,    INTENT(OUT) :: tmass1,tmass2 ! mass summation over elements
  LOGICAL, INTENT(OUT) :: sortf         ! ascending sort flag for indicies 1,2
  LOGICAL, INTENT(OUT) :: equal         ! all test variables equal
  END SUBROUTINE pufrnd
  END INTERFACE
!-------------------------------------------------------------------------------

! need at least two elements to merge
  IF(KPM.LT.2)RETURN
  KP=1
  INDEX1=NSORT(KP)

! find pointer to first puff (non-particle 1,2,3,4) index
  HDWPX=HDWP(INDEX1) 
  IF(HDWPX.GE.100)HDWPX=MOD(HDWP(INDEX1)/10,10)
  DO WHILE ((HDWPX.LT.1.OR.HDWPX.GT.4).AND.KP.LE.KPM)
     KP=KP+1
     INDEX1=NSORT(KP)
  END DO

! none to merge in the array (only particles)
  IF(KP.GE.KPM)RETURN

! temporary pollutant mass array
  NUMTYP=SIZE(mass,1)
  ALLOCATE (mass0(numtyp),STAT=kret)
  IF(kret.ne.0)THEN
     WRITE(*,*)'*ERROR* pufmrg: memory allocation' 
     STOP 900
  END IF

! point to first element to test
  DO WHILE (KP.LT.KPM)
     INDEX1=NSORT(KP)
     EQUAL=.TRUE.

!    horizontal grid distance (m) for rounding
     HGD=HGKM(PGRD(INDEX1))*1000.0

     KN=KP
!    loop indicies until true (no sort required)
     DO WHILE (EQUAL.AND.KN.LT.KPM)
        KN=KN+1
        INDEX2=NSORT(KN)
        CALL PUFRND(ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,                      &
           XPOS(INDEX1),YPOS(INDEX1),ZPOS(INDEX1),SIGH(INDEX1),            &
           SIGW(INDEX1),HDWP(INDEX1),PTYP(INDEX1),                         &
           MASS(:,INDEX1),TMASS1,                                          &
           XPOS(INDEX2),YPOS(INDEX2),ZPOS(INDEX2),SIGH(INDEX2),            &
           SIGW(INDEX2),HDWP(INDEX2),PTYP(INDEX2),                         &
           MASS(:,INDEX2),TMASS2,SORTF,EQUAL)
     END DO

!    shift pointer back if last required a sort
     IF(.NOT.EQUAL)KN=KN-1

!    low mass merge option skip remaining high mass particles
     IF(TMASS1.GT.0.0.OR.TMASS2.GT.0.0)THEN
        DEALLOCATE (mass0)
        RETURN
     END IF

!    more than one observation in merge set
     IF(KN.GT.KP)THEN
!       zero out summation variables
        XPOS0=0.0
        YPOS0=0.0
        ZPOS0=0.0
        SIGH0=0.0
        SIGW0=0.0
        PAGE0=0.0

        MASS0(1)=0.0
        MM=NUMTYP
        DO WHILE(MM.GT.1)
           MASS0(MM)=0.0
           MM=MM-1
        END DO

!       sum each element weighted by total element mass
        DO KK=KP,KN
           INDEX0=NSORT(KK)
           TMASS=MASS(1,INDEX0)
           MASS0(1)=MASS0(1)+MASS(1,INDEX0)
           MM=NUMTYP
           DO WHILE(MM.GT.1)
              TMASS=TMASS+MASS(MM,INDEX0)
              MASS0(MM)=MASS0(MM)+MASS(MM,INDEX0)
              MM=MM-1
           END DO

           XPOS0=XPOS0+XPOS(INDEX0)*TMASS
           YPOS0=YPOS0+YPOS(INDEX0)*TMASS
           ZPOS0=ZPOS0+ZPOS(INDEX0)*TMASS
           PAGE0=ABS(PAGE0)+FLOAT(ABS(PAGE(INDEX0)))*TMASS
           SIGH0=SIGH0+SIGH(INDEX0)*TMASS
           SIGW0=SIGW0+SIGW(INDEX0)*TMASS
        END DO

!       compute total mass from elemental sums
!       and replace sums into first element
        INDEX0=NSORT(KP)
        TMASS=MASS0(1)
        MASS(1,INDEX0)=MASS0(1)
        MM=NUMTYP
        DO WHILE(MM.GT.1)
           MASS(MM,INDEX0)=MASS0(MM)
           TMASS=TMASS+MASS0(MM)
           MM=MM-1
        END DO

!       compute mass weighted real variables
!       integer marker variables are unchanged
        XPOS(INDEX0)=XPOS0/TMASS
        YPOS(INDEX0)=YPOS0/TMASS
        ZPOS(INDEX0)=ZPOS0/TMASS
        PAGE(INDEX0)=NINT(ABS(PAGE0)/TMASS)
        SIGH(INDEX0)=SIGH0/TMASS
        SIGW(INDEX0)=SIGW0/TMASS

!       grid identification to zero for deleted elements
        DO KK=(KP+1),KN
           INDEX0=NSORT(KK)
           PGRD(INDEX0)=0
        END DO
     END IF
!    start next search cycle
     KP=KN+1

  END DO

  DEALLOCATE (mass0)

END SUBROUTINE pufmrg
