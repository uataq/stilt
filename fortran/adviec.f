!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVIEC           ADVection via Improved Euler-Cauchy
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION VIA IMPROVED EULER-CAUCHY COMPUTES THE 3D ADVECTION OF
!   POINT IN SPACE AND TIME USING A TWO-STEP PROCESS.  THE ADVECTION
!   THE RESULTS OF AN AVERAGE VELOCITY AT THE INITIAL POINT IN SPACE
!   TIME AND THE VELOCITY AT THE FIRST GUESS POSITION (USING THE INITIAL
!   VELOCITY).  ALL INTERPOLATION IS LINEAR.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Nov 1997 (RRD)
!                 26 Aug 1999 (RRD) - common block to define vertical grid
!                 14 Aug 2000 (RRD) - check solution on quadratic equation
!                 28 Sep 2000 (RRD) - fortran90 upgrade
!                 12 Mar 2001 (RRD) - global cyclic boundary conditions
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 16 Dec 2002 (RRD) - temporal interpolation code standard
!                 13 Jan 2003 (RRD) - temporal interpolation 2nd pass
!                 10 Apr 2003 (RRD) - replaced DM with MTIME 
!                 16 Sep 2004 (RRD) - iterative time step if needed ...
!                 25 Oct 2004 (RRD) - tratio passed through argument list
!
! Results will differ from standard integration due to internal precision is
! maintained while it appears in the standard approach precision is lost when
! the results are passed back to the main routine. The internal time step
! adjustment is triggered automatically when the advection distance exceeds
! 0.75 of the grid distance. The 0.75 is set in the namelist through
! the tratio variable. The automated integration feature should usually be
! invoked when the integration time step is forced to a value larger than the
! stability parameter and grid resolution permits.
!
! USAGE: CALL ADVIEC(U,V,W,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,
!                    DT,TRATIO,BACK,GLOBAL,NXP,NYP)
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

!dwen(20090810) SUBROUTINE ADVIEC(U,V,W,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,   &
!dwen(20090810)                  DT,TRATIO,BACK,GLOBAL,NXP,NYP)
!dwen(20090810) ****************

! JCL:(8/13/01) add grd ht (ZTER) & 1st sigma level (ZSG1) as input to calculate AGL for accurate interpolation
! JCL:(3/1/01) add WWOUT used to output the Wbar [sigma/min]
! JCL: 'GD' is array of grid size (m)--needed to convert windspeed from (grid/min)=>(m/s)
! CHG(09/11/03) Pass on RAMSFLG
! CHG(09/11/03) Pass on density fields D1,D2
! JCL:(11/03/03) remove winderror arguments--do all calculations in HYMODELC
 
SUBROUTINE ADVIEC(U,V,W,d,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,   &
     DT,TRATIO,BACK,GLOBAL,NXP,NYP,                    &
     zter,zsg1,wwout,awrfflg,fluxflg,zsg,ramsflg,dead)

  IMPLICIT NONE

  !-------------------------------------------------------------------------------
  ! argument list definitions
  !-------------------------------------------------------------------------------

  REAL,    INTENT(IN)    :: u(:,:,:,:)      ! u wind
  REAL,    INTENT(IN)    :: v(:,:,:,:)      ! v wind
  REAL,    INTENT(IN)    :: w(:,:,:,:)      ! w wind
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
  !dwen(20090810) *******************
  real,    intent(in)    :: d(:,:,:,:)      !air density
  real,    intent(in)    :: zter            !terrain height
  real,    intent(in)    :: zsg1            !1st sigma level
  real,    intent(out)   :: wwout           !wbar
  real,    intent(in)    :: zsg(:)          
  logical, intent(in)    :: awrfflg,fluxflg,ramsflg
  logical, intent(out)   :: dead
  !dwen ******************************

  !-------------------------------------------------------------------------------
  ! internal variables
  !-------------------------------------------------------------------------------

  LOGICAL               :: step
  INTEGER               :: ipass,nxs,nys,delt,tsum
  REAL                  :: delx,xt,yt,zt,zagl,tf,dist,aa,bb,cc
  REAL                  :: uu,uu1,uu2,vv,vv1,vv2,ww,ww1,ww2
  REAL                  :: SUX,SUY,SVX,SVY

  !dwen(20090819)***********************
  real                  :: zagl1,xxt,yyt,zzz,  &
       dd1,dd2,zxw,xxw,yyw
  !-------------------------------------------------------------------------------
  ! external variables
  !-------------------------------------------------------------------------------

  ! common block to pass parameters for vertical grid
  COMMON /ZZTOKK/ AA,BB,CC

  ! velocity offsets for staggered grids
  COMMON /STAGGER/ SUX,SUY,SVX,SVY

  !-------------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE ADV3NT(S,XT,YT,ZX,SS,GLOBAL,NXP,NYP)
       IMPLICIT NONE
       REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
       REAL,      INTENT(IN)    :: xt,yt         ! position of interpolated value
       REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
       REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
       LOGICAL,   INTENT(IN)    :: global        ! global cyclic boundary conditions
       INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries
     END SUBROUTINE adv3nt

     SUBROUTINE ADV3NTWIND(zagl,zagl1,S,X1,Y1,ZX,SS,GLOBAL,NXP,NYP)
       IMPLICIT NONE
       REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
       REAL,      INTENT(IN)    :: x1,y1         ! position of interpolated value
       REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
       REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
       LOGICAL,   INTENT(IN)    :: global        ! cyclic boundary condition flag
       INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries
       real,      intent(in)    :: zagl,zagl1
     end subroutine adv3ntwind

  END INTERFACE
  !-------------------------------------------------------------------------------

  ! dimension size for meteo subgrid (really just the size of the subgrid array passed in -
  ! the actual subgrid could be smaller)
  ! Now only used for array bounds checking, and deadening of particles
  NXS=SIZE(U,1)
  NYS=SIZE(U,2)
  dead=.FALSE.

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
        !     ZAGL=ZMDL*(1.0-MIN(1.0,ZT))

        !    relation between index and sigma when zsfc=0
        !    note height not correct Zagl for particle, but it is 
        !    only needed to get the correct vertical index number
        !     DIST=(BB*BB-4.0*AA*(CC-ZAGL))
        !     IF(DIST.GE.0.0)THEN
        !        ZX=(-BB+SQRT(DIST))/(2.0*AA)
        !     ELSE
        !        ZX=1.0
        !     END IF
        !     ZX=MIN(MAX(1.0,ZX),FLOAT(NLVL))
        !----old STILT way

        IF (.NOT.RAMSFLG) CALL ind_zsg (zmdl,zsg,nlvl,zt,zx,aa,bb,cc)
        ! NOTE: This follows old stilt: 
        !        JCL:(8/13/01) enable vertical index to be < 1.0 for interpolation
        !        between 1st level and ground surface
        ZX = MIN(MAX(0.0,ZX),FLOAT(NLVL))

        !dwen(20090810) ******************
        ! JCL:(11/1/02) use Draxler formulation of sigma-coordinate, w/ 'terrain compression factor'
        ! JCL:(8/13/01) calculate ht of first sigma level and ht of current particle [m AGL]
        ! CHG(09/10/03) correct transformation between sigma and agl

        ZAGL1 = (1.0-ZSG1)*(ZMDL-ZTER) 
        ZAGL = (1.0-ZT)*(ZMDL-ZTER) 

        IF (RAMSFLG) THEN

           ! CHG(09/11/03) use this for RAMS
           XXT = XT-0.5
           YYT = YT-0.5
           ! check if off grid since sometimes OFFG is not set
           IF (XXT < 1.0 .OR. YYT < 1.0 .OR. INT(XT) > NXS .OR. INT(YT) > NYS) THEN
              write (*,*) 'subroutine adviec: xxt/xt and/or yyt/yt exceeds allowable limit:', &
                   & XXT,XT,YYT,YT,NXS,NYS
              dead = .TRUE.
              RETURN
           END IF

           !dwen(20090824)            XX = DNINT(XT)                            !DNINT for scalars; for x-fluxes use X1, not XX
           !dwen(20090824)            YY = DNINT(YT)                            !DNINT for scalars; for y-fluxes use YT, not YY
           XXw = aNINT(XT)                            !DNINT for scalars; for x-fluxes use X1, not XX
           YYw = aNINT(YT)                            !DNINT for scalars; for y-fluxes use YT, not YY

           ! vertical grid: 1st level = 1st flux level;
           ! so assume ZX = 0.3, means need scalars at 1, and fluxes between 0 and 1
           ! so assume ZX = 0.7, means need scalars at 1, and fluxes between 0 and 1
           !to get scalar from proper level
           !dwen(20090824)            ZZZ = DINT(ZX+1.0)
           ZZZ = aINT(ZX+1.0)
           !dwen(20090824)            ZZZ = DMIN1(DMAX1(DBLE(1.0),ZZZ),DBLE(NLVL))
           ZZZ = AMIN1(AMAX1(1.0,ZZZ),float(NLVL))
           !dwen(20090824)            IF (XXT < 1d0  .OR.  YYT < 1d0  .OR.  XX > NXS  .OR.  YY > NYS) THEN
           IF (XXT < 1.0  .OR.  YYT < 1.0  .OR.  XXw > NXS  .OR.  YYw > NYS) THEN
              write (*,*) 'subroutine adviec: xxt/xxw and/or yyt/yyw exceeds allowable limit:', &
                   & XXT,xxw,YYT,yyw,NXS,NYS
              dead = .TRUE.
              RETURN
           END IF

           !  interpolate to position at last time
           !dwen(20090810)  CALL ADVINT (U1,NXS,NYS,NZM,XXT,YY ,ZZZ,GLOBAL,NXP,NYP,UU1)
           CALL ADV3NT(U(:,:,:,K1),XXT,YYw,ZZZ,UU1,GLOBAL,NXP,NYP)
           !dwen(20090810)  CALL ADVINT (V1,NXS,NYS,NZM,XX ,YYT,ZZZ,GLOBAL,NXP,NYP,VV1)
           CALL ADV3NT(V(:,:,:,K1),XXw,YYT,ZZZ,VV1,GLOBAL,NXP,NYP)
           !  below 1st level using 'zero wind' boundary-condition at ground
           !dwen(20090810)  CALL ADVINTWIND (ZAGL,ZAGL1,W1,NXS,NYS,NZM,XX,YY,ZX,GLOBAL,NXP,NYP,WW1)
           CALL ADV3NTWIND (ZAGL,ZAGL1,W(:,:,:,K1),XXw,YYw,ZX,WW1,GLOBAL,NXP,NYP)

           !dwen(20090810)  CALL ADVINT (D1,NXS,NYS,NZM,XX,YY,ZZZ,GLOBAL,NXP,NYP,DD1)
           CALL ADV3NT(D(:,:,:,K1),XXw,YYw,ZZZ,DD1,GLOBAL,NXP,NYP)

           !  interpolate to position at next time
           !dwen(20090810)     CALL ADVINT (U2,NXS,NYS,NZM,XXT,YY ,ZZZ,GLOBAL,NXP,NYP,UU2)
           CALL ADV3NT(U(:,:,:,K2),XXT,YYw,ZZZ,UU2,GLOBAL,NXP,NYP)
           !dwen(20090810)     CALL ADVINT (V2,NXS,NYS,NZM,XX ,YYT,ZZZ,GLOBAL,NXP,NYP,VV2)
           CALL ADV3NT(V(:,:,:,K2),XXw,YYT,ZZZ,VV2,GLOBAL,NXP,NYP)
           !dwen(20090810)  CALL ADVINTWIND (ZAGL,ZAGL1,W2,NXS,NYS,NZM,XX,YY,ZX,GLOBAL,NXP,NYP,WW2)
           CALL ADV3NTWIND (ZAGL,ZAGL1,W(:,:,:,K2),XXw,YYw,ZX,WW2,GLOBAL,NXP,NYP)
           !dwen(20090810)   CALL ADVINT (D2,NXS,NYS,NZM,XX ,YY ,ZZZ,GLOBAL,NXP,NYP,DD2)
           CALL ADV3NT(D(:,:,:,K2),XXw,YYw,ZZZ,DD2,GLOBAL,NXP,NYP)

           ! CHG divide by density
           UU1 = UU1/DD1
           VV1 = VV1/DD1
           WW1 = WW1/DD1
           UU2 = UU2/DD2
           VV2 = VV2/DD2
           WW2 = WW2/DD2
           !  interpolate to current time
           IF (BACK) THEN
              UU = UU1
              VV = VV1
              WW = WW1
           ELSE
              UU = UU2
              VV = VV2
              WW = WW2
           END IF

           !WRITE(*,*)'adviec: TF',TF,JET
           !WRITE(45,*) ZX,WW

        ELSE IF (awrfflg) THEN

           !U(V) is staggered in x(y)-direction (C-grid) in WRF
           ! Note that x1 corresponds to position on mass grid, but staggered direction
           ! is off 0.5 (first staggered gridpoint is at position -0.5, so x1 on mass
           ! grid corresponds to x1+0.5 on staggered grid - this is different from RAMS)
           ! check if off grid since sometimes OFFG is not set
           IF (XT < 1.0  .OR.  YT < 1.0  .OR.  CEILING(XT+SUX) > NXS  .OR.  CEILING(YT+SVY) > NYS) THEN
              write (*,*) 'subroutine adviec: xt/xt+sux and/or yt/yt+svy exceeds allowable limit:', &
                   & XT,xt+sux,YT,YT+svy,NXS,NYS
              dead = .TRUE.
              RETURN
           END IF
           zxw = min(zx+0.5,real(nlvl))
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,U1,NXS,NYS,NZM,XXT,YT ,ZX ,GLOBAL,NXP,NYP,UU1)
           CALL ADV3NTWIND (ZAGL,ZAGL1,U(:,:,:,K1),XT+SUX,YT+SUY,ZX,UU1,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,V1,NXS,NYS,NZM,XT ,YYT,ZX ,GLOBAL,NXP,NYP,VV1)
           CALL ADV3NTWIND (ZAGL,ZAGL1,V(:,:,:,K1),XT+SVX,YT+SVY,ZX,VV1,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,W1,NXS,NYS,NZM,XT ,YT ,ZXW,GLOBAL,NXP,NYP,WW1)
           CALL ADV3NTWIND (ZAGL,ZAGL1,W(:,:,:,K1),XT,YT,ZXW,WW1,GLOBAL,NXP,NYP)

           !  interpolate to position at next time
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,U2,NXS,NYS,NZM,XXT,YT ,ZX ,GLOBAL,NXP,NYP,UU2)
           CALL ADV3NTWIND (ZAGL,ZAGL1,U(:,:,:,K2),XT+SUX,YT+SUY,ZX,UU2,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,V2,NXS,NYS,NZM,XT ,YYT,ZX ,GLOBAL,NXP,NYP,VV2)
           CALL ADV3NTWIND (ZAGL,ZAGL1,V(:,:,:,K2),XT+SVX,YT+SVY,ZX,VV2,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,W2,NXS,NYS,NZM,XT ,YT ,ZXW,GLOBAL,NXP,NYP,WW2)
           CALL ADV3NTWIND (ZAGL,ZAGL1,W(:,:,:,K2),XT,YT,ZXW,WW2,GLOBAL,NXP,NYP)

           IF (fluxflg) THEN
              !  use the time-averaged values, do not time-interpolate
              IF (BACK) THEN
                 UU = UU1
                 VV = VV1
                 WW = WW1
              ELSE
                 UU = UU2
                 VV = VV2
                 WW = WW2
              END IF
           ELSE
              !  interpolate to current time
              UU = (UU2-UU1)*TF+UU1
              VV = (VV2-VV1)*TF+VV1
              WW = (WW2-WW1)*TF+WW1
           END IF

        ELSE
           ! neither rams nor awrf:

           ! check if off grid since sometimes OFFG is not set
           IF (XT < 1.  .OR.  YT < 1.  .OR.  CEILING(XT) > NXS  .OR.  CEILING(YT) > NYS) THEN
              write (*,*) 'subroutine adviec: xt and/or yt exceeds allowable limit:', &
                   & XT,YT,NXS,NYS
              dead = .TRUE.
              RETURN
           END IF
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,U1,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,UU1)
           CALL ADV3NTWIND (ZAGL,ZAGL1,U(:,:,:,K1),XT,YT,ZX,UU1,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,V1,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,VV1)
           CALL ADV3NTWIND (ZAGL,ZAGL1,V(:,:,:,K1),XT,YT,ZX,VV1,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,W1,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,WW1)
           CALL ADV3NTWIND (ZAGL,ZAGL1,W(:,:,:,K1),XT,YT,ZX,WW1,GLOBAL,NXP,NYP)
           !dwen(20090810)
           !  interpolate to position at next time
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,U2,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,UU2)
           CALL ADV3NTWIND (ZAGL,ZAGL1,U(:,:,:,K2),XT,YT,ZX,UU2,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,V2,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,VV2)
           CALL ADV3NTWIND (ZAGL,ZAGL1,V(:,:,:,K2),XT,YT,ZX,VV2,GLOBAL,NXP,NYP)
           !dwen(20090810)            CALL ADVINTWIND (ZAGL,ZAGL1,W2,NXS,NYS,NZM,XT,YT,ZX,GLOBAL,NXP,NYP,WW2)
           CALL ADV3NTWIND (ZAGL,ZAGL1,W(:,:,:,K2),XT,YT,ZX,WW2,GLOBAL,NXP,NYP)

           !  interpolate to current time
           UU = (UU2-UU1)*TF+UU1
           VV = (VV2-VV1)*TF+VV1
           WW = (WW2-WW1)*TF+WW1

        END IF
        !dwen ******************************************

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
           ZT=ZZ+WW*DELT
           ! JCL:(3/1/01) store the Wbar from the 1st pass
           WWOUT = WW

           !       off grid test for particles that may approach the limits of the subgrid
           !       most are terminated outside of this routine if within 2 grid pts of edge
           IF(STEP.AND..NOT.GLOBAL)THEN
              IF(XT.LT.1.0.OR.XT.GT.FLOAT(NXS).OR.YT.LT.1.0.OR.YT.GT.FLOAT(NYS))THEN
                 !             advect distance for remaining time and then terminate
                 XT=XX+UU*(DT-TSUM)
                 YT=YY+VV*(DT-TSUM)
                 ZT=ZZ+WW*(DT-TSUM)
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
           ZZ=0.5*(ZT+ZZ+WW*DELT)
           ! JCL:(3/1/01) average the Wbar from the two passes
           WWOUT = 0.5*(WWOUT+WW)
        END IF

     END DO
  END DO

END SUBROUTINE adviec
