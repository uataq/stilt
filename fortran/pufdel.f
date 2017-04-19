!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFDEL           PUFf DELete eliminate unused puffs/particles
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF DELETE ELIMINATE UNUSED PUFFS/PARTICLES FROM ARRAY BASED
!   UPON EITHER OFF-GRID, EXCEEDING AGE, OR BELOW MINIMUM MASS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 05 Sep 2000 (RRD) - fortran90 upgrade
!                 08 May 2003 (RRD) - added alongwind sigx
!                 10 Aug 2004 (RRD) - argument list change
!                 13 Oct 2005 (RRD) - lagrangian sampling test
!                 18 Jan 2008 (RRD) - absolute value page
!
! USAGE:  CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,
!              SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,KHMAX,FMASS)
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

SUBROUTINE PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,  &
                  PGRD,NSORT,KHMAX,FMASS)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(INOUT) :: kpm        ! total number of puffs or particles
  INTEGER, INTENT(IN)    :: khmax      ! maximum hours limit
  REAL,    INTENT(IN)    :: fmass      ! lower mass limit per particle/puff
  REAL,    INTENT(INOUT) :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horizontal puff sigma 
  REAL,    INTENT(INOUT) :: sigu (:)   ! turbulence u'2    
  REAL,    INTENT(INOUT) :: sigv (:)   ! turbulence v'2
  REAL,    INTENT(INOUT) :: sigw (:)   ! turbulence w'2 or vertical puff sigma
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age since release (min)
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteorological grid of puff position
  INTEGER, INTENT(INOUT) :: nsort(:)   ! sort index array by position

!-------------------------------------------------------------------------------

  LOGICAL  :: kill  
  INTEGER  :: numtyp,khrs,mm,kp,kn
  REAL     :: tmass 

!-------------------------------------------------------------------------------

! set pollutant dimension
  NUMTYP=SIZE(mass,1)

! need at least one element
  IF(KPM.LT.1)RETURN

! pointer to existing puff index for testing
  KP=1
! pointer to new index
  KN=0

  DO WHILE (KP.LE.KPM)
!    check age against limits
     KHRS=ABS(PAGE(KP))/60

!    check mass total
     TMASS=MASS(1,KP)
     MM=NUMTYP
     DO WHILE(MM.GT.1)
        TMASS=TMASS+MASS(MM,KP)
        MM=MM-1
     END DO

!    test for puff/particle delete conditions  
     KILL=.FALSE.
     IF(PGRD(KP).EQ.0.OR.KHRS.GE.KHMAX)THEN
        KILL=.TRUE.
     ELSEIF(HDWP(KP).NE.6)THEN
        IF(TMASS.LE.FMASS) KILL=.TRUE.
     END IF

     IF(KILL)THEN
!       puff/particle to delete
        KP=KP+1

     ELSE
!       puff/particle to keep
        KN=KN+1
        IF(KP.NE.KN)THEN
!          shift position in array
           XPOS(KN)=XPOS(KP)
           YPOS(KN)=YPOS(KP)
           ZPOS(KN)=ZPOS(KP)
           PAGE(KN)=PAGE(KP)
           SIGH(KN)=SIGH(KP)
           SIGU(KN)=SIGU(KP)
           SIGV(KN)=SIGV(KP)
           SIGW(KN)=SIGW(KP)
           HDWP(KN)=HDWP(KP)
           PTYP(KN)=PTYP(KP)
           PGRD(KN)=PGRD(KP)
           MASS(1,KN)=MASS(1,KP)
           MM=NUMTYP
           DO WHILE(MM.GT.1)
              MASS(MM,KN)=MASS(MM,KP)
              MM=MM-1
           END DO
        END IF

!       save index position in dummy array
        NSORT(KN)=KN

        KP=KP+1
     END IF
  END DO

! reset counter to new total
  KPM=KN

END SUBROUTINE pufdel
