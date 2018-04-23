!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVSFC           ADVection on a SurFaCe
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:03-09-16
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION VIA IMPROVED EULER-CAUCHY COMPUTES THE 2D ADVECTION OF A
!   POINT IN SPACE AND TIME USING A TWO-STEP PROCESS.  THE ADVECTION
!   THE RESULTS OF AN AVERAGE VELOCITY AT THE INITIAL POINT IN SPACE
!   TIME AND THE VELOCITY AT THE FIRST GUESS POSITION (USING THE INITIAL
!   VELOCITY).  ALL INTERPOLATION IS LINEAR.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 16 Sep 2003 (RRD) - initial version from adviec
!                 14 Apr 2006 (RRD) - staggered velocity fields
!
! USAGE: CALL ADVSFC(U,V,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,DT,BACK,
!                    GLOBAL,NXP,NYP)
!
!   INPUT ARGUMENT LIST:       see below
!   OUTPUT ARGUMENT LIST:      see below
!   INPUT FILES:               none
!   OUTPUT FILES:              none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE ADVSFC(U,V,K1,K2,NLVL,MTIME,JET,XX,YY,DT,BACK,GLOBAL,NXP,NYP)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  REAL,    INTENT(IN)    :: u(:,:,:)        ! u wind
  REAL,    INTENT(IN)    :: v(:,:,:)        ! v wind
  INTEGER, INTENT(IN)    :: k1,k2           ! last and next time index    
  INTEGER, INTENT(IN)    :: nlvl            ! number of levels to process
  INTEGER, INTENT(IN)    :: mtime(2)        ! time of meteo observations   
  INTEGER, INTENT(IN)    :: jet             ! current elapsed time (minutes)
  REAL,    INTENT(IN)    :: dt              ! integration step (minutes)
  LOGICAL, INTENT(IN)    :: back            ! flag to indicate direction
  CHARACTER(2), INTENT(IN)    :: global          ! global cyclic boundary conditions
  INTEGER, INTENT(IN)    :: nxp,nyp         ! global boundaries           
  REAL,    INTENT(INOUT) :: xx,yy           ! old (t) and new (t+dt) position 

!-------------------------------------------------------------------------------
! internal variables

  REAL                  :: uu,uu1,uu2,vv,vv1,vv2
  REAL                  :: xt,yt,tf
  INTEGER               :: ipass
  REAL                  :: SUX,SUY,SVX,SVY

!-------------------------------------------------------------------------------
! external variables

! velocity offsets for staggered grids
  COMMON /STAGGER/ SUX,SUY,SVX,SVY

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE ADV2NT(S,XP,YP,SS,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: s(:,:)        ! field for interpolation
  REAL,      INTENT(IN)    :: xp,yp         ! position of interpolated value
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  CHARACTER(2),   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundary values
  END SUBROUTINE adv2nt
  END INTERFACE
!-------------------------------------------------------------------------------

! need to save initial position in two-pass integration
  XT=XX
  YT=YY

! determine interpolation time factor
  TF=FLOAT(ABS(JET-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))

! two pass to obtain average from velocity at old and new
  DO IPASS=1,2

!    interpolate to position at last time
     CALL ADV2NT(U(:,:,K1),XT+SUX,YT+SUY,UU1,GLOBAL,NXP,NYP)
     CALL ADV2NT(V(:,:,K1),XT+SVX,YT+SVY,VV1,GLOBAL,NXP,NYP)

!    interpolate to position at next time
     CALL ADV2NT(U(:,:,K2),XT+SUX,YT+SUY,UU2,GLOBAL,NXP,NYP)
     CALL ADV2NT(V(:,:,K2),XT+SVX,YT+SVY,VV2,GLOBAL,NXP,NYP)

!    interpolate to current time
     UU=(UU2-UU1)*TF+UU1
     VV=(VV2-VV1)*TF+VV1

     IF(IPASS.EQ.1)THEN
!       first pass position simple integration
        XT=XX+UU*DT
        YT=YY+VV*DT

!       code section added 13 Jan 2003
        TF=FLOAT(ABS(JET+NINT(DT)-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))

     ELSE
!       final pass average of first and second positions
        XX=0.5*(XT+XX+UU*DT)
        YY=0.5*(YT+YY+VV*DT)
     END IF

  END DO

END SUBROUTINE advsfc
