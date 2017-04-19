!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFSRT           binary PUFf SoRTting by position
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF SORTTING IS THE BINARY SORT OF PUFFS BY ROUNDED POSITION
!   NEEDS TO BE CALLED IMMEDIATELY AFTER PUFDEL
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 (RRD)
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 12 Aug 2004 (RRD) - variable name change
!
! USAGE:  CALL PUFSRT(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,XPOS,YPOS,ZPOS,MASS,  
!                     SIGH,SIGW,HDWP,PTYP,PGRD,NSORT,HGKM)
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

SUBROUTINE PUFSRT(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,XPOS,YPOS,ZPOS,MASS,  &
                  SIGH,SIGW,HDWP,PTYP,PGRD,NSORT,HGKM)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,    INTENT(IN)    :: frach    ! horizontal position rounding fraction
  REAL,    INTENT(IN)    :: fracv    ! vertical position rounding fraction
  REAL,    INTENT(IN)    :: fract    ! travel-time rounding fraction
  REAL,    INTENT(IN)    :: fmass    ! mass sort cutoff value (0 for none)
  INTEGER, INTENT(IN)    :: kpm      ! total number of puffs or particles
  REAL,    INTENT(IN)    :: zmdl     ! vertical model scaling height
  REAL,    INTENT(IN)    :: xpos (:) ! puff center positions (grid units)
  REAL,    INTENT(IN)    :: ypos (:) ! puff center positions (grid units)
  REAL,    INTENT(IN)    :: zpos (:) ! puff center height (sigma)
  REAL,    INTENT(IN)    :: mass (:,:) ! pollutant mass
  REAL,    INTENT(IN)    :: sigh (:) ! horiz sigma (sigma)
  REAL,    INTENT(IN)    :: sigw (:) ! vert sigma (sigma)
  INTEGER, INTENT(IN)    :: hdwp (:) ! Horizontal distribution within pollutant
  INTEGER, INTENT(IN)    :: ptyp (:) ! pollutant type index number
  INTEGER, INTENT(IN)    :: pgrd (:) ! meteo grid index for puff
  INTEGER, INTENT(INOUT) :: nsort(:) ! sortted array by position
  REAL,    INTENT(IN)    :: hgkm (:) ! horizontal grid size (km)

!-------------------------------------------------------------------------------

  INTEGER, ALLOCATABLE   :: VAR(:)   ! temporary sort index
  LOGICAL                :: SORTF    ! sort flag 
  LOGICAL                :: EQUAL    ! equality flag
  REAL                   :: hgd,tmass2,tmass1 
  INTEGER                :: kret,l2,i,ig1,kn,kp,ig2,l1,len2,len1,index1,index2 

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

! need at least two elements to sort
  IF(KPM.LT.2)RETURN

! temporary sort variable
  ALLOCATE (var(kpm),STAT=kret)
  IF(kret.NE.0)THEN
     WRITE(*,*)'*ERROR* pufsrt: memory allocation' 
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

!          horizontal meteo grid distance (m) for rounding
           HGD=HGKM(PGRD(INDEX1))*1000.0

!          round position according to sigma and determine sort
           CALL PUFRND(ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,XPOS(INDEX1),         &
                YPOS(INDEX1),ZPOS(INDEX1),SIGH(INDEX1),SIGW(INDEX1),          &
                HDWP(INDEX1),PTYP(INDEX1),MASS(:,INDEX1),TMASS1,XPOS(INDEX2), &
                YPOS(INDEX2),ZPOS(INDEX2),SIGH(INDEX2),SIGW(INDEX2),          &
                HDWP(INDEX2),PTYP(INDEX2),MASS(:,INDEX2),TMASS2,SORTF,EQUAL)

           IF(SORTF)THEN
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

  DEALLOCATE (var)

END SUBROUTINE pufsrt
