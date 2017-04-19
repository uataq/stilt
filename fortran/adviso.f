!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVISO           ADVection on an ISO surface  
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION VIA IMPROVED EULER-CAUCHY COMPUTES THE 3D ADVECTION OF
!   POINT IN SPACE AND TIME USING A TWO-STEP PROCESS.  THE ADVECTION
!   THE RESULTS OF AN AVERAGE VELOCITY AT THE INITIAL POINT IN SPACE
!   TIME AND THE VELOCITY AT THE FIRST GUESS POSITION (USING THE INITIAL
!   VELOCITY).  ALL INTERPOLATION IS LINEAR. THIS SPECIAL VERSION FIRST
!   INTERPOLATES THE VELOCITY VECTORS TO A SELECTED ISOSURFACE PRIOR
!   TO COMPUTING THE ADVECTION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED:  12 Oct 2005 (RRD) - initial version from adviec
!                  14 Apr 2006 (RRD) - staggered velocity fields
!
! PROGRAM NOTES:
! Results will differ from standard integration due to internal precision is
! maintained while it appears in the standard approach precision is lost when
! the results are passed back to the main routine. The internal time step
! adjustment is triggered automatically when the advection distance exceeds
! 0.75 of the grid distance. The 0.75 is set in the namelist through
! the tratio variable. The automated integration feature should usually be
! invoked when the integration time step is forced to a value larger than the
! stability parameter and grid resolution permits.
!
! USAGE: CALL ADVISO(U,V,P,ZSG,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,
!                    DT,TRATIO,BACK,GLOBAL,NXP,NYP,UVEL,VVEL)
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

SUBROUTINE ADVISO(U,V,P,ZSG,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,   &
                  DT,TRATIO,BACK,GLOBAL,NXP,NYP,UVEL,VVEL,RAMSFLG)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  REAL,    INTENT(IN)    :: u(:,:,:,:)      ! u wind
  REAL,    INTENT(IN)    :: v(:,:,:,:)      ! v wind
  REAL,    INTENT(IN)    :: p(:,:,:,:)      ! pressure
  REAL,    INTENT(IN)    :: zsg(:)          ! sigma levels
  INTEGER, INTENT(IN)    :: k1,k2           ! last and next time index    
  INTEGER, INTENT(IN)    :: nlvl            ! number of levels to process
  INTEGER, INTENT(IN)    :: mtime(2)        ! time of meteo observations   
  INTEGER, INTENT(IN)    :: jet             ! current elapsed time (minutes)
  REAL,    INTENT(IN)    :: zmdl            ! vertical model domain top (m)
  REAL,    INTENT(IN)    :: dt              ! integration step (minutes)
  REAL,    INTENT(IN)    :: tratio          ! time step stability criterion
  LOGICAL, INTENT(IN)    :: back            ! flag to indicate direction
  LOGICAL, INTENT(IN)    :: global          ! global cyclic boundary conditions
  INTEGER, INTENT(IN)    :: nxp,nyp         ! global boundaries           

  REAL,    INTENT(INOUT) :: xx,yy,zz        ! old (t) and new (t+dt) position 
  REAL,    INTENT(OUT)   :: zx              ! last estimate of vertical index
  REAL,    INTENT(IN)    :: uvel,vvel       ! forced velocity vector
  LOGICAL, INTENT(IN)    :: RAMSFLG         ! RAMS flag

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  LOGICAL               :: step
  INTEGER               :: i,j,k,l,ipass,nxs,nys,delt,tsum
  REAL                  :: delx,xt,yt,zt,zagl,tf,zf,dist,aa,bb,cc
  REAL                  :: uu,uu1,uu2,vv,vv1,vv2,pp,pp1,pp2,zs,zs1,zs2
  REAL                  :: SUX,SUY,SVX,SVY

! interpolated single surface velocity vectors
  REAL, ALLOCATABLE     :: us(:,:,:),vs(:,:,:),ws(:,:,:)

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------

! common block to pass parameters for vertical grid
  COMMON /ZZTOKK/ AA,BB,CC

! velocity offsets for staggered grids
  COMMON /STAGGER/ SUX,SUY,SVX,SVY

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE ADV2NT(S,XP,YP,SS,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: s(:,:)        ! field for interpolation
  REAL,      INTENT(IN)    :: xp,yp         ! position of interpolated value
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  LOGICAL,   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundary values
  END SUBROUTINE adv2nt
!-------------------------------------------------------------------------------
  SUBROUTINE ADV3NT(S,XT,YT,ZX,SS,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: xt,yt         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  LOGICAL,   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries
  END SUBROUTINE adv3nt
  END INTERFACE
!-------------------------------------------------------------------------------

! dimension size for meteo subgrid
  NXS=SIZE(U,1)
  NYS=SIZE(U,2)

! determine initial interpolation time factor for input meteo data times
  TF=FLOAT(ABS(JET-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))

! initial time step and total integration time within this routine
  DELT=INT(DT)
  TSUM=0
  STEP=.FALSE.

! internal integration loop (usually only one interation required)
  DO WHILE (ABS(TSUM).LT.ABS(INT(DT)))

! need to save initial position in two-pass integration
  XT=XX
  YT=YY
  ZT=ZZ

! two pass to obtain average from velocity at old and new
  DO IPASS=1,2

!    compute vertical interpolation factor as fraction
!    from integer index based upon quadratic relation between
!    height/sigma (for Zsfc=0) and array index position

!    height agl when Zsfc=0, particle position in same ratio
!    to layer depth (Zmdl-Zsfc) as when Zsfc=0
!!$     ZAGL=ZMDL*(1.0-MIN(1.0,ZT))

!    relation between index and sigma when zsfc=0
!    note height not correct Zagl for particle, but it is 
!    only needed to get the correct vertical index number
!!$     DIST=(BB*BB-4.0*AA*(CC-ZAGL))
!!$     IF(DIST.GE.0.0)THEN
!!$        ZX=(-BB+SQRT(DIST))/(2.0*AA)
!!$     ELSE
!!$        ZX=1.0
!!$     END IF
     IF (.NOT.RAMSFLG) CALL ind_zsg (zmdl,zsg,nlvl,zt,zx,aa,bb,cc)
     ZX=MIN(MAX(1.0,ZX),FLOAT(NLVL))

!---------------------------------------------------------
!    compute two-dimensional velocity fields
!---------------------------------------------------------

     IF(IPASS.EQ.1)THEN
!       compute pressure at initial point for field interpolation
        CALL ADV3NT(P(:,:,:,K1),XT,YT,ZX,PP1,GLOBAL,NXP,NYP)
        CALL ADV3NT(P(:,:,:,K2),XT,YT,ZX,PP2,GLOBAL,NXP,NYP)
        PP=(PP2-PP1)*TF+PP1

!       two dimensional velocity vectors
        ALLOCATE (US(nxs,nys,2),VS(nxs,nys,2),WS(nxs,nys,2))

!       interpolate winds to PP surface 
        DO K=1,2
           DO J=1,NYS
           DO I=1,NXS

              WS(I,J,K)=ZSG(1)
              IF(UVEL.EQ.0.0.AND.VVEL.EQ.0.0)THEN
!                compute velocity field from meteorology 
                 US(I,J,K)=U(I,J,1,K)
                 VS(I,J,K)=V(I,J,1,K)
              ELSE
!                force velocity vector
                 US(I,J,K)=UVEL
                 VS(I,J,K)=VVEL
              END IF

              tloop : DO L=2,NLVL
                 IF(P(I,J,L,K).LE.PP.AND.P(I,J,L-1,K).GT.PP)THEN

!                   vertical interpolation factor
                    ZF=(PP-P(I,J,L-1,K))/(P(I,J,L,K)-P(I,J,L-1,K))
                    WS(I,J,K)=ZF*(ZSG(L)-ZSG(L-1))+ZSG(L-1)

                    IF(UVEL.EQ.0.0.AND.VVEL.EQ.0.0)THEN
                       US(I,J,K)=ZF*(U(I,J,L,K)-U(I,J,L-1,K))+U(I,J,L-1,K)
                       VS(I,J,K)=ZF*(V(I,J,L,K)-V(I,J,L-1,K))+V(I,J,L-1,K)
                    ELSE
                       US(I,J,K)=UVEL
                       VS(I,J,K)=VVEL
                    END IF
                    EXIT tloop
                 END IF

              END DO tloop

           END DO
           END DO
        END DO
     END IF

!---------------------------------------------------------

!    interpolate to position at last time
     CALL ADV2NT(US(:,:,K1),XT+SUX,YT+SUY,UU1,GLOBAL,NXP,NYP)
     CALL ADV2NT(VS(:,:,K1),XT+SVX,YT+SVY,VV1,GLOBAL,NXP,NYP)
     CALL ADV2NT(WS(:,:,K1),XT    ,YT    ,ZS1,GLOBAL,NXP,NYP)

!    interpolate to position at next time
     CALL ADV2NT(US(:,:,K2),XT+SUX,YT+SUY,UU2,GLOBAL,NXP,NYP)
     CALL ADV2NT(VS(:,:,K2),XT+SVX,YT+SVY,VV2,GLOBAL,NXP,NYP)
     CALL ADV2NT(WS(:,:,K2),XT    ,YT    ,ZS2,GLOBAL,NXP,NYP)

!    interpolate to current time
     UU=(UU2-UU1)*TF+UU1
     VV=(VV2-VV1)*TF+VV1
     ZS=(ZS2-ZS1)*TF+ZS1

     IF(IPASS.EQ.1)THEN

        IF(TSUM.EQ.0)THEN
           DELX=MAX(ABS(UU*DELT),ABS(VV*DELT))
!          test if maximum step exceeds 3/4 grid cell
           IF(DELX.GT.TRATIO)THEN
!             compute required internal time step
              DELT=MAX(1,INT(TRATIO*ABS(DELT)/DELX))
              DO WHILE (MOD(INT(ABS(DT)),DELT).NE.0.AND.DELT.GT.1)
                 DELT=DELT-1
              END DO
              DELT=SIGN(DELT,INT(DT))
              STEP=.TRUE.
           END IF
        END IF

!       first pass position simple integration
        XT=XX+UU*DELT
        YT=YY+VV*DELT
        ZT=ZS             

!       off grid test for particles that may approach the limits of the subgrid
!       most are terminated outside of this routine if within 2 grid pts of edge
        IF(STEP.AND..NOT.GLOBAL)THEN
           IF(XT.LT.1.0.OR.XT.GT.FLOAT(NXS).OR.YT.LT.1.0.OR.YT.GT.FLOAT(NYS))THEN
!             advect distance for remaining time and then terminate
              XX=XX+UU*(DT-TSUM)
              YY=YY+VV*(DT-TSUM)
              ZZ=0.5*(ZT+ZS)
              RETURN
           END IF
        END IF

!       code modified to include [delt/tsum] 13 Jan 2003
        TSUM=TSUM+DELT
        TF=FLOAT(ABS(JET+TSUM-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))

     ELSE
!       final pass average of first and second positions
        XX=0.5*(XT+XX+UU*DELT)
        YY=0.5*(YT+YY+VV*DELT)
        ZZ=0.5*(ZT+ZS)
     END IF

  END DO
  END DO

  DEALLOCATE (US,VS,WS)
END SUBROUTINE adviso
