!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVRKM           ADVection via Runge-Kutta Method
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION VIA THE RUNGE-KUTTA METHOD (4TH ORDER ACCURACY).
!   REPLACES ADVIEC (2ND ORDER METHOD) IN SUBROUTINE ADVPNT. ARGUMENTS
!   ARE IDENTICAL FOR THE TWO ROUTINES.  AN ADDITIONAL OPTION IS TO
!   REPLACE THE LINEAR INTERPOLATION ADV3NT USED IN THIS ROUTINE WITH THE
!   POLYNOMIAL INTERPOLATION SUBROUTINE ADVRNT. ARGUMENTS ARE IDENTICAL.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 21 Jan 2003 (RRD) - initial version from ADVIEC
!                 14 Apr 2006 (RRD) - staggered velocity fields
!
! USAGE: CALL ADVRKM(U,V,W,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,DT,BACK,
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

SUBROUTINE ADVRKM(U,V,W,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,DT,BACK,  &
                  GLOBAL,NXP,NYP)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  REAL,    INTENT(IN)    :: u(:,:,:,:)      ! u wind
  REAL,    INTENT(IN)    :: v(:,:,:,:)      ! v wind
  REAL,    INTENT(IN)    :: w(:,:,:,:)      ! w wind
  INTEGER, INTENT(IN)    :: k1,k2           ! last and next time index    
  INTEGER, INTENT(IN)    :: nlvl            ! number of levels to process
  INTEGER, INTENT(IN)    :: mtime(2)        ! time of meteo observation    
  INTEGER, INTENT(IN)    :: jet             ! current elapsed time (minutes)
  REAL,    INTENT(IN)    :: zmdl            ! vertical model domain top (m)
  REAL,    INTENT(IN)    :: dt              ! integration step (minutes)
  LOGICAL, INTENT(IN)    :: back            ! flag to indicate direction
  CHARACTER(2), INTENT(IN)    :: global          ! global cyclic boundary conditions
  INTEGER, INTENT(IN)    :: nxp,nyp         ! global boundaries           

  REAL,    INTENT(INOUT) :: xx,yy,zz        ! old (t) and new (t+dt) position 
  REAL,    INTENT(OUT)   :: zx              ! last estimate of vertical index

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL      :: RTF(4) = (/0.0,0.5,0.5,1.0/) ! interpolation factors
  REAL      :: x(0:4),y(0:4),z(0:4)         ! guess position array

  INTEGER   :: ip                 
  REAL      :: uu,uu1,uu2,vv,vv1,vv2,ww,ww1,ww2
  REAL      :: xt,yt,zt,zagl,tf,dist,aa,bb,cc
  REAL      :: SUX,SUY,SVX,SVY

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------

! common block to pass parameters for vertical grid
  COMMON /ZZTOKK/ AA,BB,CC

! velocity offsets for staggered grids
  COMMON /STAGGER/ SUX,SUY,SVX,SVY

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE ADVRNT(S,XT,YT,ZX,SS,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: xt,yt         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  CHARACTER(2),   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries
  END SUBROUTINE advrnt
  END INTERFACE
!-------------------------------------------------------------------------------

  X(0)=0.0
  Y(0)=0.0
  Z(0)=0.0

! fourth order accuracy requires multiple iterations
  DO IP=1,4

!    guess position
     XT=XX+X(IP-1)*RTF(IP)
     YT=YY+Y(IP-1)*RTF(IP)
     ZT=ZZ+Z(IP-1)*RTF(IP)

!    compute vertical interpolation factor
     ZAGL=ZMDL*(1.0-MIN(1.0,ZT))
     DIST=(BB*BB-4.0*AA*(CC-ZAGL))
     IF(DIST.GE.0.0)THEN
        ZX=(-BB+SQRT(DIST))/(2.0*AA)
     ELSE
        ZX=1.0
     END IF
     ZX=MIN(MAX(1.0,ZX),FLOAT(NLVL))

!    interpolate to position
     CALL ADVRNT(U(:,:,:,K1),XT+SUX,YT+SUY,ZX,UU1,GLOBAL,NXP,NYP)
     CALL ADVRNT(V(:,:,:,K1),XT+SVX,YT+SVY,ZX,VV1,GLOBAL,NXP,NYP)
     CALL ADVRNT(W(:,:,:,K1),XT    ,YT    ,ZX,WW1,GLOBAL,NXP,NYP)
     CALL ADVRNT(U(:,:,:,K2),XT+SUX,YT+SUY,ZX,UU2,GLOBAL,NXP,NYP)
     CALL ADVRNT(V(:,:,:,K2),XT+SVX,YT+SVY,ZX,VV2,GLOBAL,NXP,NYP)
     CALL ADVRNT(W(:,:,:,K2),XT    ,YT    ,ZX,WW2,GLOBAL,NXP,NYP)

!    interpolate to integration time
     TF=FLOAT(ABS(JET+NINT(DT*RTF(IP))-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))
     IF(IP.EQ.1)THEN
        IF(BACK)THEN
           IF(TF.EQ.0.0)TF=1.0
           TF=1.0-TF
        END IF
     ELSE
        IF(BACK)THEN
           TF=1.0-TF
        ELSE
           IF(TF.EQ.0.0)TF=1.0
        END IF
     END IF
     UU=(UU2-UU1)*TF+UU1
     VV=(VV2-VV1)*TF+VV1
     WW=(WW2-WW1)*TF+WW1

!    updated guess position
     X(IP)=UU*DT
     Y(IP)=VV*DT
     Z(IP)=WW*DT

  END DO

! Final Averaged position
  XX=XX+(X(1)+2.0*X(2)+2.0*X(3)+X(4))/6.0
  YY=YY+(Y(1)+2.0*Y(2)+2.0*Y(3)+Y(4))/6.0
  ZZ=ZZ+(Z(1)+2.0*Z(2)+2.0*Z(3)+Z(4))/6.0

END SUBROUTINE advrkm
