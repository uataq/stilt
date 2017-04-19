!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVMET           ADVection step returns local METeorology
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION STEP RETURNS LOCAL METEOROLOGY AT END POINT AFTER
!   INTEPOLATION OF METEOROLOGICAL INFORMATION FROM GRID TO POSITION
!   IN BOTH SPACE AND TIME.  NOT USED IN ADVECTION BUT BY OTHER
!   ROUTINES AS DISPERSION AND DEPOSITION.  THIS ROUTINE PROVIDES TH
!   ONLY INTERFACE OF METEOROLOGICAL INFORMATION TO NON-ADVECTION
!   SUBROUTINES.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 18 Aug 1998 (RRD) - added tension stress to deformation
!                 21 Dec 1998 (RRD) - generalized vertical index and
!                                     simplified non-interpolation option
!                                     added rdep to argument list
!                                     compile option for spatial interpolation
!                 04 Mar 1999 (RRD) - eliminate negative precip
!                 06 Apr 1999 (RRD) - created pressure profile, pass to METO
!                 30 Apr 1999 (RRD) - added terrain and height computation
!                 05 Jul 2000 (RRD) - dimensioned TMRK to return multiple
!                 28 Sep 2000 (RRD) - fortran90 upgrade
!                 12 Mar 2001 (RRD) - global lat lon grid option
!                 21 Jun 2001 (RRD) - ambient temperature array
!                 17 Oct 2001 (RRD) - new trajectory markers + mixing depth
!                 26 Feb 2002 (RRD) - downward shortwave flux
!                 23 Jul 2002 (RRD) - added rel humidity as marker variable
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 16 Dec 2002 (RRD) - temporal interpolation factor
!                 10 Apr 2003 (RRD) - replaced DM with MTIME
!                 15 Sep 2003 (RRD) - moved landuse and roughness length
!                                   - return friction velocity as scalar
!                 14 Oct 2003 (RRD) - local density computation
!                 05 Nov 2003 (RRD) - turbulent velocity components
!                 03 Dec 2004 (BS)  - precipitation rate (PRT6)
!                 11 May 2005 (RRD) - global cyclic boundary condition
!                 31 May 2005 (RRD) - terrain and solar flux along traj
!                 08 Mar 2006 (RRD) - static stability parameter
!                 13 Sep 2007 (RRD) - terrain interpolation
!
! USAGE:  CALL ADVMET(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,XP,YP,JET,MTIME,
!                     KCYCLE,NLVL,FHOUR,IFHR,K1,K2,GX,GY,Z0,LU,ZT,A,T,Q,P,
!                     E,X,ZI,H,U0,V0,RT,UF,VF,SF,SS,DS,DSWF,GLOBAL,NXP,NYP)
!       
!   INPUT ARGUMENT LIST:        see below
!   OUTPUT ARGUMENT LIST:       see below
!   INPUT FILES:                none
!   OUTPUT FILES:               none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

!dwen(20090811) SUBROUTINE ADVMET(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,XP,YP,JET,MTIME,      &
!dwen(20090811)                   KCYCLE,NLVL,FHOUR,IFHR,K1,K2,GX,GY,Z0,LU,ZT,A,T,Q,P,E,X, &
!dwen(20090811)                   ZI,H,U0,V0,RT,UF,VF,SF,SS,DS,DSWF,GLOBAL,NXP,NYP)

!dwen(20090811) ****************************
! JCL:add two arrays of Lagrangian timescales (TL1,TL2) and
!         two arrays of stddev of vertical velocity (SIGW1,SIGW2)
! JCL:(5/9/01)added arrays of horizontal velocity (U1,U2) and (V1,V2)
! CHG:(11/20/01) added conv. precip. flag CFLG to arguments
! CHG:(11/20/01)added conv. precip. rates (RC1,RC2) to arguments
! CHG:(11/20/01) added tot.cld and radiation flag (TCLF,RADF)
! CHG:(11/20/01)added tot. cld and radiation (TC1,TC2,SW1,SW2)
! CHG:(11/20/01) added energy fluxes (LF1,LF2,HF1,HF2) to arguments
! CHG:(11/20/01)added low.cld flag LCLF
! CHG:(11/20/01) added low.cld LC1, LC2
! JCL:(4/3/02)added arrays of mass violation (DMASS1,DMASS2)
! CHG:(22/01/03)added soil moisture flag SLMF
! CHG:(22/01/03) added soil moisture SM1 SM2
! CHG(09/16/03)changed input position XP,YP to XPI,YPI
! CHG(09/16/03) added RAMSFLG
! dwen(20090811) added vertical and horizontal mixing coefficent XM,HM
! dwen(20090819) added WMXC to replace VMIX and output vertical mixing profile,
!                      VMIX in new HYSPLIT represent v-component turbulence

 SUBROUTINE ADVMET(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,XPi,YPi,JET,MTIME,      &
                   KCYCLE,NLVL,FHOUR,IFHR,K1,K2,GX,GY,Z0,LU,ZT,A,T,Q,P,E,X, &
                   ZI,H,U0,V0,RT,UF,VF,SF,SS,DS,DSWF,GLOBAL,NXP,NYP,        &
                   cflg,tclf,lclf,radf,slmf,                                &
                   u,v,d,xm,hm,dmass,tl,sigw,rc,hf,lf,tc,lc,sm,          &
                   ramsflg,ecmflg,awrfflg,fluxflg)



  IMPLICIT NONE

  INCLUDE 'DEFMETO.INC'                ! meteo summary at last advection point

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(bset),INTENT(OUT)   :: metz (:)   ! profile advection variables
  TYPE(aset),INTENT(INOUT) :: meto       ! surface advection variables
  LOGICAL,   INTENT(IN)    :: back,vmix,cdep,rdep,traj  
!dwen(20090812)  REAL,      INTENT(IN)    :: xp,yp 
  REAL,      INTENT(IN)    :: xpi,ypi 
  INTEGER,   INTENT(IN)    :: jet 
  INTEGER,   INTENT(IN)    :: mtime(2)   ! time of meteo observations
  INTEGER,   INTENT(IN)    :: kcycle,nlvl,k1,k2,fhour(2)
  INTEGER,   INTENT(INOUT) :: ifhr
  REAL,      INTENT(IN)    :: gx(:,:), gy(:,:), z0(:,:), zt(:,:)
  INTEGER,   INTENT(IN)    :: lu(:,:)
  REAL,      INTENT(IN)    :: t(:,:,:,:), a(:,:,:,:) 
  REAL,      INTENT(IN)    :: q(:,:,:,:), p(:,:,:,:)
  REAL,      INTENT(IN)    :: e(:,:,:,:), x(:,:,:,:), h(:,:,:,:)
  REAL,      INTENT(IN)    :: u0(:,:,:),  v0(:,:,:),  ss(:,:,:)
  REAL,      INTENT(IN)    :: uf(:,:,:),  vf(:,:,:),  sf(:,:,:)
  REAL,      INTENT(IN)    :: zi(:,:,:),  rt(:,:,:),  ds(:,:,:)
  LOGICAL,   INTENT(IN)    :: dswf          ! downward shortwave flag 
  LOGICAL,   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundary values

!dwen(20090811) ******************
 logical,    intent(in)    :: cflg,tclf,lclf,radf,slmf
 logical,    intent(in)    :: ramsflg,ecmflg,awrfflg,fluxflg
 real,       intent(in)    :: u(:,:,:,:),v(:,:,:,:),d(:,:,:,:)
 real,       intent(in)    :: xm(:,:,:,:),hm(:,:,:,:)
 real,       intent(in)    :: dmass(:,:,:,:),tl(:,:,:,:)
 real,       intent(in)    :: sigw(:,:,:,:),rc(:,:,:)
 real,       intent(in)    :: hf(:,:,:),lf(:,:,:),tc(:,:,:)
 real,       intent(in)    :: lc(:,:,:),sm(:,:,:)
! ********************************

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL,          PARAMETER :: rdry  = 287.04  ! dry air (J/Kg-K)
  REAL,          PARAMETER :: p2jm  = 100.0   ! mb to j/m3

  INTEGER                  :: kl,ii,jj
  REAL                     :: xc,yc,tf,tr,ea,sea,rfhr,zx,zk,var1,var2
!dwen(20090812) ******************
  real                     :: xp,yp,xx2,yy2,xx,yy,dens1,dens2,zz,ztmp,   &
                              swfc1,swfc2,swfc
  integer                  :: jet1,jet2,jetdel,jetmax
!dwen ****************************

  ! velocity offsets for staggered grids
  REAL                  :: SUX,SUY,SVX,SVY
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
  SUBROUTINE ADV3NT(S,XP,YP,ZX,SS,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: xp,yp         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  LOGICAL,   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundary values
  END SUBROUTINE adv3nt
!-------------------------------------------------------------------------------
  SUBROUTINE SUNFLX(NLVL,SEA,QQ,SWF,TR)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: nlvl      ! number of vertical levels in profile
  REAL,     INTENT(IN)   :: sea       ! sine of the solar elevation angle 
  REAL,     INTENT(IN)   :: qq   (:)  ! RH fraction profile 
  REAL,     INTENT(OUT)  :: swf       ! incident short wave flux (w/m2)
  REAL,     INTENT(OUT)  :: tr        ! transmissivity
  END SUBROUTINE sunflx
!-------------------------------------------------------------------------------
  SUBROUTINE SUNANG(JET,OLAT,OLON,EA,SEA)
  INTEGER,  INTENT(IN)   :: jet       ! elapsed minutes since January 1st 1970
  REAL,     INTENT(IN)   :: olat      ! latitude (+ = north)
  REAL,     INTENT(IN)   :: olon      ! longitude (- = west)
  REAL,     INTENT(OUT)  :: ea        ! solar elevation angle in deg at time
  REAL,     INTENT(OUT)  :: sea       ! sine of EA
  END SUBROUTINE sunang

!dwen(20090917)********************
      SUBROUTINE SUNAVE(JET,JETDEL,OLAT,OLON,SWF)
      IMPLICIT none
      integer, intent(in) :: jet,jetdel
      real,    intent(in) :: olat,olon
      real,    intent(out):: swf
      end subroutine sunave
!**********************************
  END INTERFACE
!-------------------------------------------------------------------------------

!dwen(20090812) **************
! CHG(09/16/03) use proper grid stagger, see comments in advpnt.f
      xp=xpi
      yp=ypi
      IF (RAMSFLG) THEN
!dwen(20090824)        XP = DNINT(XPI)
!dwen(20090824)        YP = DNINT(YPI)
        XP = aNINT(XPI)
        YP = aNINT(YPI)
      end if
      if (awrfflg) then
           XX2 = Xpi+0.5
           YY2 = Ypi+0.5
      end if
!dwen *************************

! interpolation factor for current time
  TF=FLOAT(ABS(JET-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))

! compute current forecast hour
  RFHR=(FHOUR(K2)-FHOUR(K1))*TF+FHOUR(K1)
  IFHR=NINT(RFHR)

! save position indicies within meteo subgrid
  II=NINT(XP)
  JJ=NINT(YP)
  ZX=METO%ZNDX
  ZX = aMIN1(aMAX1(1.0,ZX),float(NLVL))

! cyclic boundary conditions
  IF(GLOBAL)THEN
     IF(II.GT.NXP)II=1
     IF(II.LT.1)  II=NXP
     IF(JJ.GT.NYP)JJ=NYP
     IF(JJ.LT.1)  JJ=1
  END IF

! grid distance at point
  METO%GDISX=GX(II,JJ)
  METO%GDISY=GY(II,JJ)

! terrain surface elevation
  CALL ADV2NT(ZT(:,:),XP,YP,VAR1,GLOBAL,NXP,NYP)
  METO%ZTER=VAR1

! roughness length 
  METO%AERO=Z0(II,JJ)

! land-use category 
  METO%LAND=LU(II,JJ)


!dwen(20090812) *************************
! CHG:(11/20/01) energy fluxes
!        interpolate sensible heat flux to point and time
!dwen         CALL ADVINT2D (HF1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADV2NT (HF(:,:,k1),XP,YP,var1,GLOBAL,NXP,NYP)
!dwen         CALL ADVINT2D (HF2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         CALL ADV2NT (HF(:,:,k2),XP,YP,var2,GLOBAL,NXP,NYP)
         METO%SHTF = (VAR2-VAR1)*TF+VAR1
!        interpolate latent heat flux to point and time
!dwen         CALL ADVINT2D (LF1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADV2NT (LF(:,:,k1),XP,YP,var1,GLOBAL,NXP,NYP)
!dwen         CALL ADVINT2D (LF2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         CALL ADV2NT (LF(:,:,k2),XP,YP,var2,GLOBAL,NXP,NYP)
         METO%LHTF = (VAR2-VAR1)*TF+VAR1

! CHG:(11/20/01) total cloud cover and downw. radiation
      IF (TCLF) THEN
!        interpolate total cloud cover to point and time
!dwen         CALL ADVINT2D (TC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADV2NT (TC(:,:,k1),XP,YP,var1,GLOBAL,NXP,NYP)
!dwen         CALL ADVINT2D (TC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         CALL ADV2NT (TC(:,:,k2),XP,YP,var2,GLOBAL,NXP,NYP)
         METO%TCLD = (VAR2-VAR1)*TF+VAR1
      ELSE
         METO%TCLD = -99.0
      END IF

! CHG:(22/01/03) soil moisture
      IF (SLMF) THEN
!        interpolate soil moisture to point and time
!dwen         CALL ADVINT2D (SM1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADV2NT (SM(:,:,k1),XP,YP,var1,GLOBAL,NXP,NYP)
!dwen         CALL ADVINT2D (SM2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         CALL ADV2NT (SM(:,:,k2),XP,YP,var2,GLOBAL,NXP,NYP)
         METO%SOLW = (VAR2-VAR1)*TF+VAR1
      ELSE
         METO%SOLW = -99.0
      END IF

! CHG:(12/04/01) low cloud cover
      IF (LCLF) THEN
!        interpolate low cloud cover to point and time
!dwen         CALL ADVINT2D (LC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
         CALL ADV2NT (LC(:,:,k1),XP,YP,var1,GLOBAL,NXP,NYP)
!dwen         CALL ADVINT2D (LC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
         CALL ADV2NT (LC(:,:,k2),XP,YP,var2,GLOBAL,NXP,NYP)
         METO%LCLD = (VAR2-VAR1)*TF+VAR1
      ELSE
         METO%LCLD = -99.0
      END IF
!dwen ****************************************

!-------------------------------------------------------------------------------
! trajectory option only returns simple marker variable
!-------------------------------------------------------------------------------

  IF(TRAJ)THEN
     IF(METO%FLAG%PRES.GT.0)THEN
!       all trajectories return pressure
        CALL ADV3NT(P(:,:,:,K1),XP,YP,ZX,VAR1,GLOBAL,NXP,NYP)
        CALL ADV3NT(P(:,:,:,K2),XP,YP,ZX,VAR2,GLOBAL,NXP,NYP)
        METO%TMRK(METO%FLAG%PRES)=(VAR2-VAR1)*TF+VAR1
     END IF

     IF(METO%FLAG%TPOT.GT.0)THEN
!       isentropic -> save potential temperature
        CALL ADV3NT(T(:,:,:,K1),XP,YP,ZX,VAR1,GLOBAL,NXP,NYP)
        CALL ADV3NT(T(:,:,:,K2),XP,YP,ZX,VAR2,GLOBAL,NXP,NYP)
        METO%TMRK(METO%FLAG%TPOT)=(VAR2-VAR1)*TF+VAR1
     END IF

     IF(METO%FLAG%TAMB.GT.0)THEN
!       ambient temperature
        CALL ADV3NT(A(:,:,:,K1),XP,YP,ZX,VAR1,GLOBAL,NXP,NYP)
        CALL ADV3NT(A(:,:,:,K2),XP,YP,ZX,VAR2,GLOBAL,NXP,NYP)
        METO%TMRK(METO%FLAG%TAMB)=(VAR2-VAR1)*TF+VAR1
     END IF

     IF(METO%FLAG%RELH.GT.0)THEN
!       relative humidity  
        CALL ADV3NT(Q(:,:,:,K1),XP,YP,ZX,VAR1,GLOBAL,NXP,NYP)
        CALL ADV3NT(Q(:,:,:,K2),XP,YP,ZX,VAR2,GLOBAL,NXP,NYP)
        METO%TMRK(METO%FLAG%RELH)=((VAR2-VAR1)*TF+VAR1)*100.0
     END IF

     IF(METO%FLAG%TERR.GT.0)THEN 
!       terrain height
        METO%TMRK(METO%FLAG%TERR)=METO%ZTER
     END IF
  END IF

!-------------------------------------------------------------------------------
! parameter profiles used in deposition and mixing calculations
!-------------------------------------------------------------------------------

  IF(VMIX.AND.METO%FLAG%MIXD.GT.0)THEN
!    mixing depth                   
     CALL ADV2NT(ZI(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(ZI(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
     METO%MIXD=(VAR2-VAR1)*TF+VAR1
     METO%TMRK(METO%FLAG%MIXD)=METO%MIXD 
  END IF

  IF(.NOT.TRAJ.AND.VMIX)THEN

!    to obtain the u-component turbulence (m2/s2)
     CALL ADV3NT(H(:,:,:,K1),XP,YP,ZX,VAR1,GLOBAL,NXP,NYP)
     CALL ADV3NT(H(:,:,:,K2),XP,YP,ZX,VAR2,GLOBAL,NXP,NYP)
     METO%UMIX=(VAR2-VAR1)*TF+VAR1

!    to obtain the v-component turbulence (m2/s2)
     CALL ADV3NT(E(:,:,:,K1),XP,YP,ZX,VAR1,GLOBAL,NXP,NYP)
     CALL ADV3NT(E(:,:,:,K2),XP,YP,ZX,VAR2,GLOBAL,NXP,NYP)
     METO%VMIX=(VAR2-VAR1)*TF+VAR1

!dwen(20090812) *************************
!        to obtain the horizontal diffusivity in m2/sec
! JCL:(4/26/02) since ZX could be less than 1.0, prevent ADVINT from accessing 0th element in array
         ZTMP = ZX
         IF (ZX < 1.0) ZTMP = 1.0
! CHG(09/16/03) use next higher level for RAMS vertical (see adviec comments CHG)
!dwen(20090824)         IF (RAMSFLG) ZTMP = DINT(ZX+1.0)
!dwen(20090824)         ZTMP = DMIN1(DMAX1(DBLE(1.0),ZTMP),DBLE(NLVL))
         IF (RAMSFLG) ZTMP = aINT(ZX+1.0)
         ZTMP = aMIN1(aMAX1(1.0,ZTMP),float(NLVL))
!dwen         CALL ADVINT (H1,NXS,NYS,NZM,XP,YP,ZTMP,GLOBAL,NXP,NYP,VAR1)
         CALL ADV3NT (hm(:,:,:,k1),XP,YP,ZTMP,var1,GLOBAL,NXP,NYP)
!dwen         CALL ADVINT (H2,NXS,NYS,NZM,XP,YP,ZTMP,GLOBAL,NXP,NYP,VAR2)
         CALL ADV3NT (hm(:,:,:,k2),XP,YP,ZTMP,var2,GLOBAL,NXP,NYP)
         METO%HMIX = (VAR2-VAR1)*TF+VAR1
!dwen ***********************************

     DO KL=1,NLVL
!       set vertical interpolation point to index position
        ZK=KL

!       to obtain the w-component turbulence (m2/s2)
        CALL ADV3NT(X(:,:,:,K1),XP,YP,ZK,VAR1,GLOBAL,NXP,NYP)
        CALL ADV3NT(X(:,:,:,K2),XP,YP,ZK,VAR2,GLOBAL,NXP,NYP)
        METZ(KL)%WMIX=(VAR2-VAR1)*TF+VAR1

!dwen(20090812) *******************
! JCL:(5/9/01)interpolate horizontal winds to position
! CHG(09/16/03) don't use this for RAMS
            IF (.NOT.RAMSFLG) THEN
               if (.not. awrfflg) then
!dwen                  CALL ADVINT (U1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                  CALL ADV3NT (U(:,:,:,k1),XP,YP,ZK,var1,GLOBAL,NXP,NYP)
!dwen                  CALL ADVINT (U2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
                  CALL ADV3NT (U(:,:,:,k2),XP,YP,ZK,var2,GLOBAL,NXP,NYP)
!dwen(20090826)                  METz%UUNEXT(KL) = (VAR2-VAR1)*TF+VAR1
                  METz(kl)%UUNEXT = (VAR2-VAR1)*TF+VAR1
!dwen                  CALL ADVINT (V1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                  CALL ADV3NT (V(:,:,:,k1),XP,YP,ZK,var1,GLOBAL,NXP,NYP)
!dwen                  CALL ADVINT (V2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
                  CALL ADV3NT (V(:,:,:,k2),XP,YP,ZK,var2,GLOBAL,NXP,NYP)
!dwen(20090826)                  METz%VVNEXT(KL) = (VAR2-VAR1)*TF+VAR1
                  METz(kl)%VVNEXT = (VAR2-VAR1)*TF+VAR1
               else
! for awrf:
                  if (.not. fluxflg) then
! instantantenous velocities, on staggered C-grid
!dwen                     CALL ADVINT (U1,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                     CALL ADV3NT (U(:,:,:,k1),XP+SUX,YP+SUY,ZK,var1,GLOBAL,NXP,NYP)
!dwen                     CALL ADVINT (U2,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR2)
                     CALL ADV3NT (U(:,:,:,k2),XP+SUX,YP+SUY,ZK,var2,GLOBAL,NXP,NYP)
!dwen(20090826)                     METz%UUNEXT(KL) = (VAR2-VAR1)*TF+VAR1
                     METz(kl)%UUNEXT = (VAR2-VAR1)*TF+VAR1
!dwen                     CALL ADVINT (V1,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR1)
                     CALL ADV3NT (V(:,:,:,k1),XP+SVX,YP+SVY,ZK,var1,GLOBAL,NXP,NYP)
!dwen                     CALL ADVINT (V2,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR2)
                     CALL ADV3NT (V(:,:,:,k2),XP+SVX,YP+SVY,ZK,var2,GLOBAL,NXP,NYP)
!dwen(20090826)                     METz%VVNEXT(KL) = (VAR2-VAR1)*TF+VAR1
                     METz(kl)%VVNEXT = (VAR2-VAR1)*TF+VAR1
                  else
! time-averaged coupled u,v (decoupled in prfcom): (on C-grid)
                     IF (BACK) THEN
!dwen                        CALL ADVINT (U1,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                        CALL ADV3NT (U(:,:,:,k1),XP+SUX,YP+SUY,ZK,var1,GLOBAL,NXP,NYP)
!dwen                        CALL ADVINT (V1,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR2)
                        CALL ADV3NT (V(:,:,:,k1),XP+SVX,YP+SVY,ZK,var2,GLOBAL,NXP,NYP)
                     else
!dwen                        CALL ADVINT (U2,NXS,NYS,NZM,Xx2,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                        CALL ADV3NT (U(:,:,:,k2),XP+SUX,YP+SUY,ZK,var1,GLOBAL,NXP,NYP)
!dwen                        CALL ADVINT (V2,NXS,NYS,NZM,XP,Yy2,ZK,GLOBAL,NXP,NYP,VAR2)
                        CALL ADV3NT (V(:,:,:,k2),XP+SVX,YP+SVY,ZK,var2,GLOBAL,NXP,NYP)
                     end IF !back
!dwen(20090826)                     METz%UUNEXT(KL) = VAR1
                     METz(kl)%UUNEXT = VAR1
!dwen(20090826)                     METz%VVNEXT(KL) = VAR2
                     METz(kl)%VVNEXT = VAR2
                  end if !fluxflg
               end if !awrfflg
                 !of IF(.NOT.RAMSFLG)THEN
            ELSE
! CHG(09/16/03) use this for RAMS
! use later time for fluxes
              IF (BACK) THEN
!dwen                CALL ADVINT (U1,NXS,NYS,NZM,XPI,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                CALL ADV3NT (U(:,:,:,k1),XPI,YP,ZK,var1,GLOBAL,NXP,NYP)
!dwen                CALL ADVINT (V1,NXS,NYS,NZM,XP,YPI,ZK,GLOBAL,NXP,NYP,VAR2)
                CALL ADV3NT (V(:,:,:,k1),XP,YPI,ZK,var2,GLOBAL,NXP,NYP)
              ELSE
!dwen                CALL ADVINT (U2,NXS,NYS,NZM,XPI,YP,ZK,GLOBAL,NXP,NYP,VAR1)
                CALL ADV3NT (U(:,:,:,k2),XPI,YP,ZK,var1,GLOBAL,NXP,NYP)
!dwen                CALL ADVINT (V2,NXS,NYS,NZM,XP,YPI,ZK,GLOBAL,NXP,NYP,VAR2)
                CALL ADV3NT (V(:,:,:,k2),XP,YPI,ZK,var2,GLOBAL,NXP,NYP)
              END IF
!dwen(20090826)              METz%UUNEXT(KL) = VAR1
!dwen(20090826)              METz%VVNEXT(KL) = VAR2
              METz(kl)%UUNEXT = VAR1
              METz(kl)%VVNEXT = VAR2
                   !of IF(.NOT.RAMSFLG)THEN ... ELSE ...
            END IF

!           interpolate mixing to position
!dwen            CALL ADVINT (X1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
            CALL ADV3NT (xm(:,:,:,k1),XP,YP,ZK,var1,GLOBAL,NXP,NYP)
!dwen            CALL ADVINT (X2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
            CALL ADV3NT (xm(:,:,:,k2),XP,YP,ZK,var2,GLOBAL,NXP,NYP)

!****************************************
!dwen(20090819): VMIX in the newest HYSPLIT stand for v-component turbulence(m2/s2), instead of vertical mixing porfile (m2/s). To keep consistent with the current STILT, try to add VMXC to represent vertical mixing profile(m2/s).
!            METO%VMIX(KL) = (VAR2-VAR1)*TF+VAR1
            METz(kl)%wmxc = (VAR2-VAR1)*TF+VAR1
!****************************************

! JCL:(2/23/01) need density for correction in PARDSP, so should extract
!                 density profile even when CDEP = F
!           interpolate density to position
!dwen            CALL ADVINT (D1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,DENS1)
            CALL ADV3NT (D(:,:,:,k1),XP,YP,ZK,dens1,GLOBAL,NXP,NYP)
!dwen            CALL ADVINT (D2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,DENS2)
            CALL ADV3NT (D(:,:,:,k2),XP,YP,ZK,dens2,GLOBAL,NXP,NYP)
!dwen(20090826)            METz%DENS(KL) = (DENS2-DENS1)*TF+DENS1
            METz(kl)%DENS = (DENS2-DENS1)*TF+DENS1

! JCL:(4/3/02)should NOT have spatial interpolation of mass violation b/c mass violation
!           was calculated for a gridcell, so should keep same value
! CHG(10/10/03) special for RAMS
            IF (.NOT.RAMSFLG) THEN
!dwen(20090824)              XX = DINT(XP)
!dwen(20090824)              YY = DINT(YP)
!dwen(20090824)              ZZ = DINT(ZK)
              XX = aINT(XP)
              YY = aINT(YP)
              ZZ = aINT(ZK)
                 !for RAMS
            ELSE
              XX = XP
              YY = YP
              ZZ = DBLE(ZK)
            ENDIF
!dwen            CALL ADVINT (DMASS1,NXS,NYS,NZM,XX,YY,ZZ,GLOBAL,NXP,NYP,VAR1)
            CALL ADV3NT (DMASS(:,:,:,k1),XX,YY,ZZ,var1,GLOBAL,NXP,NYP)
!dwen            CALL ADVINT (DMASS2,NXS,NYS,NZM,XX,YY,ZZ,GLOBAL,NXP,NYP,VAR2)
            CALL ADV3NT (DMASS(:,:,:,k2),XX,YY,ZZ,var2,GLOBAL,NXP,NYP)
            IF (.NOT. (RAMSFLG .OR. fluxflg) .OR. ECMFLG) THEN
!dwen(20090826)               METz%DMASSNEXT(KL) = (VAR2-VAR1)*TF+VAR1
               METz(kl)%DMASSNEXT = (VAR2-VAR1)*TF+VAR1
            ELSE IF (BACK) THEN
!dwen(20090826)                  METz%DMASSNEXT(KL) = VAR1
                  METz(kl)%DMASSNEXT = VAR1
               ELSE
!dwen(20090826)                  METz%DMASSNEXT(KL) = VAR2
                  METz(kl)%DMASSNEXT = VAR2
            END IF
! JCL:(4/3/02)the density change term is extremely small, so not worry
!          also need the density change term (kg/m3/min) as part of mass budget!
!dwen            METO%DMASSNEXT(KL) = METO%DMASSNEXT(KL)*2.0/(DENS2+DENS1)     &
!dwen                        +((DENS2-DENS1)/((DENS2+DENS1)/2.0))/DM
!dwen(20090826)           METz%DMASSNEXT(KL) = METz%DMASSNEXT(KL)*2.0/(DENS2+DENS1)     &
           METz(kl)%DMASSNEXT = METz(kl)%DMASSNEXT*2.0/(DENS2+DENS1)     &
               +((DENS2-DENS1)/((DENS2+DENS1)/2.0))/ABS(MTIME(2)-MTIME(1))
!        take into account the sign associated with direction in TIME
!dwen(20090826)            IF (BACK) METz%DMASSNEXT(KL) = -1.0*METz%DMASSNEXT(KL)
            IF (BACK) METz(kl)%DMASSNEXT = -1.0*METz(kl)%DMASSNEXT


! JCL:(2/28/01) want to output surface temperature, which uses pressure, so should
!                extract pressure profile even when CDEP = F
!              interpolate pressure to position (currently not saved)
!dwen            CALL ADVINT (P1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
            CALL ADV3NT (P(:,:,:,k1),XP,YP,ZK,var1,GLOBAL,NXP,NYP)
!dwen            CALL ADVINT (P2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
            CALL ADV3NT (P(:,:,:,k2),XP,YP,ZK,var2,GLOBAL,NXP,NYP)
!dwen(20090826)            METz%PRES(KL) = (VAR2-VAR1)*TF+VAR1
            METz(kl)%PRES = (VAR2-VAR1)*TF+VAR1

! JCL:(2/28/01) want to output surface temperature, so should extract
!                temperature profile even when CDEP = F
!           compute temperature from pressure and density
!dwen(20090826)            METz%TEMP(KL) = 100.0*METz%PRES(KL)/METz%DENS(KL)/287.04
            METz(kl)%TEMP = 100.0*METz(kl)%PRES/METz(kl)%DENS/287.04

! JCL:(2/28/01) relative humidity profile needed to calculate incident
!                solar radiation in SUNFLX, so extract even when CDEP = F
!           interpolate humidity to position
!dwen            CALL ADVINT (Q1,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR1)
            CALL ADV3NT (Q(:,:,:,k1),XP,YP,ZK,var1,GLOBAL,NXP,NYP)
!dwen            CALL ADVINT (Q2,NXS,NYS,NZM,XP,YP,ZK,GLOBAL,NXP,NYP,VAR2)
            CALL ADV3NT (Q(:,:,:,k2),XP,YP,ZK,var2,GLOBAL,NXP,NYP)
!dwen(20090826)            METz%RHFR(KL) = (VAR2-VAR1)*TF+VAR1
            METz(kl)%RHFR = (VAR2-VAR1)*TF+VAR1
!dwen **************************************

        IF(CDEP)THEN
!          interpolate pressure to position 
           CALL ADV3NT(P(:,:,:,K1),XP,YP,ZK,VAR1,GLOBAL,NXP,NYP)
           CALL ADV3NT(P(:,:,:,K2),XP,YP,ZK,VAR2,GLOBAL,NXP,NYP)
           METZ(KL)%PRES=(VAR2-VAR1)*TF+VAR1

!          interpolate humidity to position
           CALL ADV3NT(Q(:,:,:,K1),XP,YP,ZK,VAR1,GLOBAL,NXP,NYP)
           CALL ADV3NT(Q(:,:,:,K2),XP,YP,ZK,VAR2,GLOBAL,NXP,NYP)
           METZ(KL)%RHFR=(VAR2-VAR1)*TF+VAR1

!          ambient temperature
           CALL ADV3NT(A(:,:,:,K1),XP,YP,ZK,VAR1,GLOBAL,NXP,NYP)
           CALL ADV3NT(A(:,:,:,K2),XP,YP,ZK,VAR2,GLOBAL,NXP,NYP)
           METZ(KL)%TEMP=(VAR2-VAR1)*TF+VAR1

!          use pressure and temperature to find density
           METZ(KL)%DENS=(P2JM*METZ(KL)%PRES)/(METZ(KL)%TEMP*RDRY)
        END IF

     END DO
  END IF

!-------------------------------------------------------------------------------
! surface level parameters for resistance and deposition calculations
!-------------------------------------------------------------------------------

  IF(RDEP)THEN
!    interpolate friction velocity to point and time
     CALL ADV2NT(UF(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(UF(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
!    convert from gp/min to m/sec
     XC=((VAR2-VAR1)*TF+VAR1)*GX(II,JJ)/60.0

     CALL ADV2NT(VF(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(VF(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
!    convert from gp/min to m/sec
     YC=((VAR2-VAR1)*TF+VAR1)*GY(II,JJ)/60.0
!    scalar friction velocity
     METO%USTR=SQRT(XC*XC+YC*YC)

!    interpolate stability function to point and time
     CALL ADV2NT(SF(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(SF(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
     METO%PSI=(VAR2-VAR1)*TF+VAR1

!    static stability parameter 
     CALL ADV2NT(SS(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(SS(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
     METO%SSP=(VAR2-VAR1)*TF+VAR1

!    low-level scalar wind      
     CALL ADV2NT(U0(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(U0(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
     XC=((VAR2-VAR1)*TF+VAR1)

     CALL ADV2NT(V0(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(V0(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
     YC=((VAR2-VAR1)*TF+VAR1)
     METO%UBAR=SQRT(XC*XC+YC*YC)

!    mixing depth                   
     CALL ADV2NT(ZI(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
     CALL ADV2NT(ZI(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
     METO%MIXD=(VAR2-VAR1)*TF+VAR1
  END IF
 
!dwen(20090916)  IF(RDEP.OR.(TRAJ.AND.METO%FLAG%DSWF.GT.0))THEN
  IF(radf.or.RDEP.OR.(TRAJ.AND.METO%FLAG%DSWF.GT.0))THEN
!    interpolate downward shortwave flux
     IF(DSWF)THEN
        CALL ADV2NT(DS(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
        CALL ADV2NT(DS(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)
!dwen(20090915)*******************
! USE STILT CODES
! CHG:(9/24/02) better interpolate "cloudiness factor"
! i.e. DSWF/DSFW(clear)
! so first get clear sky DSWF
!dwen(20090915):replace DM with MTIME
!         IF (BACK) THEN
!! CHG:(1/28/03) fixed bug in JET1 and JET2
!            JET1 = JET+TF*DM
!            JET2 = JET-(1.0-TF)*DM
!         ELSE
!            JET1 = JET-TF*DM
!            JET2 = JET+(1.0-TF)*DM
!         END IF

         jet1 = jet-tf*float((mtime(2)-mtime(1)))
         jet2 = jet+(1.0-tf)*float((mtime(2)-mtime(1)))

!dwen(20090915): replace tlat and tlon with meto%plat, meto%plon
!         CALL SUNAVE(JETMAX,JETDEL,TLAT,TLON,SWFC1)
!         CALL SUNCLR (JET,TLAT,TLON,SWFC)

         JETDEL=IABS(JET1-JET2)
         JETMAX=MAX0(JET1,JET2)
         CALL SUNAVE(JETMAX,JETDEL,meto%pLAT,meto%pLON,SWFC1)
         CALL SUNCLR (JET,meto%PLAT,meto%PLON,SWFC)
         !for other fields than ECMWF use instantaneous clear sky radiation,
         !since radiation is not time integrated over met. output time step
         if(.not. ecmflg) SWFC1=SWFC
         
! Ratio DSWF(cld) to DSWF(clear); VAR? is now "cloudiness factor"
         IF (VAR1 < 10d0*SWFC1 .AND. SWFC1 > 0d0) THEN
            VAR1 = VAR1/SWFC1
         ELSE
            VAR1 = -1d0
         END IF
         IF (VAR2 < 10d0*SWFC1 .AND. SWFC1 > 0d0) THEN
            VAR2=VAR2/SWFC1
         ELSE
            VAR2 = -1d0
         END IF
         IF(JET1 > JET2) THEN
            IF(VAR1 < 0d0) THEN
               METO%DSWF=0d0
            ELSE
               METO%DSWF=VAR1*SWFC
            ENDIF
         ELSE
            IF(VAR2 < 0d0) THEN
              METO%DSWF=0d0
            ELSE
               METO%DSWF=VAR2*SWFC
            ENDIF
         ENDIF
!dwen****************************************
     ELSE
        CONTINUE
!       computation of solar angle only required for gaseous dry deposition
!       or certain specialized chemistry applications
        CALL SUNANG(JET,METO%PLAT,METO%PLON,EA,SEA)

!       compute solar flux from solar angle and humidity
!       required for resistance gaseous dry deposition
        CALL SUNFLX(NLVL,SEA,METZ%RHFR,meto%dswf,TR)

    END IF
     IF(METO%FLAG%DSWF.GT.0) METO%TMRK(METO%FLAG%DSWF)=meto%dswf 
!**************************************
!dwen(20090831) set default to meto%dswf
     else
     meto%dswf=-999.0
!dwen ********************************
 END IF

!dwen(20090812) **************************
! CHG:(11/20/01) get rain anyway
!   get rain anyway according to STILT codes
!dwen(20090812)  IF(CDEP.OR.(TRAJ.AND.METO%FLAG%RAIN.GT.0))THEN
!    interpolate 2D precipitation total (m) to position
     IF(KCYCLE.GT.0)THEN
        CALL ADV2NT(RT(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
        CALL ADV2NT(RT(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)

!       check if previous time requires zero accumulation
!       it is an even multiple rain accumulation initialization
        IF(MOD(FHOUR(K1)*60,KCYCLE).EQ.0)VAR1=0.0
!       convert accumulation to precipitation rate (m/min)
        METO%RAIN=MAX(0.0,(VAR2-VAR1)/ABS(MTIME(2)-MTIME(1)))

     ELSEIF(KCYCLE.LT.0)THEN
        CALL ADV2NT(RT(:,:,K1),XP,YP,VAR1,GLOBAL,NXP,NYP)
        CALL ADV2NT(RT(:,:,K2),XP,YP,VAR2,GLOBAL,NXP,NYP)

!       it is an even multiple rain accumulation initialization
        IF(MOD(FHOUR(K1)*60,ABS(KCYCLE)).EQ.0)THEN
!          use precipitation rate (m/min) during 0 to +3 (rate const over 3-h)
           METO%RAIN=MAX(0.0,VAR2)
        ELSE
!          calc rate for times +3 to +6 
!               [ (6-h rate * 6 hours) - (3-h rate * 3 hours) ] / 3 hours
           METO%RAIN=MAX(0.0,(VAR2*2.0-VAR1))
        END IF

     ELSE
!       accumlation cycle not set implies no precip field
        METO%RAIN=0.0
     END IF
!    convert from m/min to mm/hr
     IF(TRAJ)METO%TMRK(METO%FLAG%RAIN)=METO%RAIN*60000.0
!CHG:(11/20/01) get rain anyway    
!dwen(20090812)  END IF

!dwen(20090812) **************************
      IF (CFLG) THEN
         !  interpolate 2D precipitation total (m) to position
         IF (KCYCLE > 0) THEN
!dwen            CALL ADVINT2D (RC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
            CALL ADV2NT (RC(:,:,k1),XP,YP,var1,GLOBAL,NXP,NYP)
!dwen            CALL ADVINT2D (RC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
            CALL ADV2NT (RC(:,:,k2),XP,YP,var2,GLOBAL,NXP,NYP)
            ! optional statement to simplify previous code
            !VAR2 = RC2(II,JJ)
            !VAR1=RC1(II,JJ)

            ! check if previous time requires zero accumulation
            ! it is an even multiple rain accumulation initialization
            IF (MOD(FHOUR(K1)*60,KCYCLE) == 0) VAR1 = 0.0
!           convert accumulation to precipitation rate (m/min)
!dwen            METO%CRAI = MAX(DBLE(0.0),(VAR2-VAR1)/DM)
           METO%CRAI = MAX(0.0,(VAR2-VAR1)/ABS(MTIME(2)-MTIME(1)))
         ELSE IF (ECMFLG) THEN                  ! instantaneous fields (averaged of 3 hours before)
            ! interpolate convective rain (in m/s) to point and time
!dwen            CALL ADVINT2D (RC1,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR1)
            CALL ADV2NT (RC(:,:,k1),XP,YP,var1,GLOBAL,NXP,NYP)
!dwen            CALL ADVINT2D (RC2,NXS,NYS,XP,YP,GLOBAL,NXP,NYP,VAR2)
            CALL ADV2NT (RC(:,:,k2),XP,YP,var2,GLOBAL,NXP,NYP)
            METO%CRAI = ((VAR2-VAR1)*TF+VAR1)*60d0 ! output is in m/min
         ELSE
!           accumlation cycle not set implies no precip field
            METO%CRAI = -0.9
         END IF
      ELSE
! CHG:(11/20/01) if no conv. precip avail., set to -999
         METO%CRAI = -0.9
      END IF
!dwen *****************************************************

END SUBROUTINE advmet
