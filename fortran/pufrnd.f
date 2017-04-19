!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFRND           PUFf RouNDs the puff position for sortting
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF ROUND ROUNDS THE PUFF POSITION ACCORDING TO DISPERSION SO
!   IF TWO POSITIONS ARE WITHIN THE FRACTION*SIGMA THEN
!   THEY ARE ASSUMED TO BE IN THE SAME POSITION FOR SORTING
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 12 Aug 2004 (RRD) - variable name change
!                 13 Nov 2008 (RRD) - enhanced test on hdwp
!
! USAGE:  CALL PUFRND(
!      ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,XPOS1,YPOS1,
!      ZPOS1,SIGH1,SIGW1,HDWP1,PTYP1,MASS1,TMASS1,XPOS2,YPOS2,
!      ZPOS2,SIGH2,SIGW2,HDWP2,PTYP2,MASS2,TMASS2,SORTF,EQUAL)
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

SUBROUTINE PUFRND( ZMDL,FRACH,FRACV,FRACT,FMASS,HGD,                         &
           XPOS1,YPOS1,ZPOS1,SIGH1,SIGW1,HDWPX,PTYP1,MASS1,TMASS1,           &
           XPOS2,YPOS2,ZPOS2,SIGH2,SIGW2,HDWPY,PTYP2,MASS2,TMASS2,SORTF,EQUAL)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

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
  INTEGER, INTENT(IN)  :: hdwpx, hdwpy  ! distribution
  INTEGER, INTENT(IN)  :: ptyp1, ptyp2  ! pollutant type
  REAL,    INTENT(IN)  :: mass1  (:)    ! pollutant mass
  REAL,    INTENT(IN)  :: mass2  (:)    ! pollutant mass
  REAL,    INTENT(OUT) :: tmass1,tmass2 ! mass summation over elements
  LOGICAL, INTENT(OUT) :: sortf         ! ascending sort flag for indicies 1,2
  LOGICAL, INTENT(OUT) :: equal         ! all test variables equal

!-------------------------------------------------------------------------------

  REAL     :: sprec,hprec,vprec,xmass1,xmass2,hdwp1,hdwp2
  INTEGER  :: mm 

!-------------------------------------------------------------------------------

! determine any mass sorting options
  IF(FMASS.GT.0.0)THEN
     XMASS1=MASS1(1)
     XMASS2=MASS2(1)

     MM=SIZE(MASS1,1)
     DO WHILE(MM.GT.1)
        XMASS1=XMASS1+MASS1(MM)
        XMASS2=XMASS2+MASS2(MM)
        MM=MM-1
     END DO

!    set mass values to 0 or 1 depending on fmass
     TMASS1=0.0
     TMASS2=0.0
     IF(XMASS1.GT.FMASS)TMASS1=1.0
     IF(XMASS2.GT.FMASS)TMASS2=1.0
  ELSE
     TMASS1=0.0
     TMASS2=0.0
  END IF

! horizontal sigma precision
  SPREC=FRACT*MIN(SIGH1,SIGH2)
! horizontal position precision
  HPREC=FRACH*MIN(SIGH1,SIGH2)/HGD

! horizontal distribution (simple or complex)
  HDWP1=HDWPX
  IF(HDWP1.GE.100)HDWP1=MOD(HDWPX/10,10)
  HDWP2=HDWPY
  IF(HDWP2.GE.100)HDWP2=MOD(HDWPY/10,10)

! vertical precision depends upon distribution
  VPREC=0.0
  IF(HDWP1.EQ.HDWP2)THEN
     IF(HDWP1.EQ.1.OR.HDWP1.EQ.2)THEN
!       vertical top-hat
        VPREC=FRACV*MIN(SIGW1,SIGW2)
     ELSE
!       vertical particle
        IF(SIGW1.EQ.0.0.OR.SIGW2.EQ.0.0)THEN
!          if either puff has just split then merge
!          for particles assume sigma-v to be 10% of sigma-h
           VPREC=0.10*MIN(SIGH1,SIGH2)/ZMDL
!          with a limit of 10% of the model depth (0-1)
           VPREC=FRACV*MIN(0.10,VPREC)
        ELSE
!          no merge until just after horizontal split
           VPREC=0.0
        END IF
     END IF
  END IF

! ascending order sort rules
  SORTF=.FALSE.
  EQUAL=.FALSE.

! horizontal distribution type
  IF(HDWP2.LT.HDWP1)THEN
     SORTF=.TRUE.
  ELSEIF(HDWP2.EQ.HDWP1)THEN

! pollutant type
  IF(PTYP2.LT.PTYP1)THEN
     SORTF=.TRUE.
  ELSEIF(PTYP2.EQ.PTYP1)THEN

! horizontal sigma same as age
  IF(SIGH2-SIGH1+SPREC.LT.0.0)THEN
     SORTF=.TRUE.
  ELSEIF(ABS(SIGH2-SIGH1).LT.SPREC)THEN

! particle mass
  IF(TMASS2.LT.TMASS1)THEN
     SORTF=.TRUE.
  ELSEIF(TMASS2.EQ.TMASS1)THEN

! vertical position
  IF(ZPOS2-ZPOS1+VPREC.LT.0.0)THEN
     SORTF=.TRUE.
  ELSEIF(ABS(ZPOS2-ZPOS1).LT.VPREC)THEN

! east-west position
  IF(XPOS2-XPOS1+HPREC.LT.0.0)THEN
     SORTF=.TRUE.
  ELSEIF(ABS(XPOS2-XPOS1).LT.HPREC)THEN

! north-south position
  IF(YPOS2-YPOS1+HPREC.LT.0.0)THEN
     SORTF=.TRUE.
  ELSEIF(ABS(YPOS2-YPOS1).LT.HPREC)THEN
!    total equality
     EQUAL=.TRUE.
     RETURN
  END IF

  END IF
  END IF
  END IF
  END IF
  END IF
  END IF

END SUBROUTINE pufrnd
