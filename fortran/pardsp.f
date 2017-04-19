!$$$  SUBPROGRAM DOCUMENTATION BLOCKMIN1MIN1
!
! SUBPROGRAM:  PARDSP           PARticle DISpersion turbulent component
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PARTICLE DISPERSION ADDS THE TURBULENT VELOCITY COMPONENT TO THE
!   PREVIOUSLY CALCULATED POSITION.  BASED UPON DIFFUSIVITY PROVIDING
!   THE BASIS FOR SIGMA U,V,W.  LOW TURBULENCE CORRECTION ONLY APPLIES
!   TO THE VERTICAL COMPONENT.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Nov 1997 (RRD)
!                 18 Aug 1998 (RRD) - isotropic turbulence option
!                 20 Apr 1999 (RRD) - terrain compression adjustment
!                 26 Aug 1999 (RRD) - pass vertical grid parameters in common
!                 01 Mar 2000 (RRD) - correction to vertical autocorrelation
!                 06 Jul 2000 (RRD) - convert time step to seconds
!                                     reverse sign of velocity after reflection
!                                     correction to vertical index algorithm
!                 14 Aug 2000 (RRD) - quadratic equation test
!                 20 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Mar 2001 (RRD) - global lat lon grid option
!                 23 May 2002 (RRD) - time step based upon entire profile
!                 16 Aug 2002 (AS)  - chock & winkler non-homogeneous turb
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 15 Sep 2003 (RRD) - more focused test on HDWP
!                 10 Nov 2003 (RRD) - convert from mixing to turbulence
!                 12 Oct 2005 (RRD) - lagrangian sampling option test
!                 25 May 2006 (AS)  - variable seed for initialization
!                 17 Oct 2007 (RRD) - forced time step positive
!                 08 Aug 2008 (RRD) - time scales become variables
!
! USAGE:  CALL PARDSP(UMIX,VMIX,GDISX,GDISY,DTM,ZMDL,ZSFC,NLVL,WMIX,ZSG,XPOS,
!              YPOS,ZPOS,VPRIME,WRATIO,UPRIME,HDWP,ZNDX,ISEED)
!
!   INPUT ARGUMENT LIST:   see definitions below
!   OUTPUT ARGUMENT LIST:  see definitions below
!   INPUT FILES:           none
!   OUTPUT FILES:          none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$
!************************************************ 
! JCL:(5/16/2000)Conduct well-mixed test, but implement Thomson [1997] scheme for
!   treatment of interfaces between step changes in high & low turbulence
! JCL:in this version of PARDSP, switch off vertical interpolation
! JCL:add 'SMIN'--refers to number of minutes since start of model
! JCL:add 'SEEVEG'--flag that gets set to true if particle sees vegetation
!     add 'TLFRAC'--fraction of Lagrangian timescale that determines internal timestep
!     add 'RSEED'--the random seed that will be fed to PARVAR
!     add 'SIGMAW'--OUTPUT: std dev of vertical velocity
!     add 'SIGMAU'--OUTPUT: std dev of horizontal velocity
!     add 'SAMPTT'--OUTPUT: amount of time [min.] that particle 'sees' the ground
!     add 'TLPREV'--vertical profile of Lagrangian timescales [s] @ t
!     add 'STDWPREV'--vertical profile of std dev of vert velocity @ t
!     add 'ZML1'--mixed-layer ht @ t
!     add 'DENSPREV'--air density profile [kg/m3] @ t
!     add 'TLNEXT'--vertical profile of Lagrangian timescales [s] @ t+dt
!     add 'STDWNEXT'--vertical profile of std dev of vert velocity @ t+dt
!     add 'ZML2'--mixed-layer ht @ t+dt
!     add 'DENSNEXT'--air density profile [kg/m3] @ t+dt
! CHG:Add 'SAMPTT2'--output for Tlagrange vertical
!            difference betw STDW & SIGMAW: SIGW is INPUT--determined in STBSND
! JCL:unlike 'chgpardsp.f', have CONSTANT profiles of TLPREV,TLNEXT,STDWPREV,&STDWNEXT
!            so not take these as arguments;also remove ZSG
! JCL:(5/18/00)add as 1st argument--# of min. since start
! JCL:(5/18/00)'WWPREV' is output variable that stores WPRIME/SIGMAW
! JCL:(5/18/00)Arguments when have TIME-VARYING vertical profiles of TL & SIGMAW
! JCL:(6/29/00)'VEGHT' is ht [m] below which a particle would be counted as 'seeing' grd vegetation
! JCL:(5/9/01)added horizontal position before mean advection (XPOSPREV & YPOSPREV)
!                  + vertical profiles of horizontal winds (UUPREV,UUNEXT,VVPREV,VVNEXT)
!             -these are used to implement interaction between vertical turbulence and wind shear
!             -added the flag 'BACK' to tell whether running model backwards in time or not
! JCL:(4/3/02)added vertical profile of mass violation before & after mean advection
! JCL:(4/3/02)'DMWEIGHT' is total weight of particle changed by mass violation [fraction/grid] experienced by particle in PARDSP
! CHG:(12/04/01)add 'CONVDUR' duration for conv. redistribution, used as FLAG (and DURATION in future?)
! CHG:(12/04/01)add 'ZLOC' limit of convection altitude (m)
! set to 0 when conv. redistribution was done (input & output)
! CHG(09/11/03) Pass on RAMSFLG
! CHG:(03/17/2004) pass on rel. humidity and temperature profile to get dry air density column integrated
!  also provide output for specific humidity, and output for "foot" (footprint in ppm per micro-moles/m2/s)
! JCL:(07/14/2004) split up GDIS=>GDISX&GDISY for global grids (non-conformal)
!************************************************ 

!************************************************ 
!dwen(20090318): add some arguments according to STILT
!dwen(20090816): remove SMIN,useless
!dwen(20090817): remove seeveg,isot and ecmflg from argument list, useless in PARDSP
!dwen(20090816): Starting from v4.6, HYSPLIT used WRATIO, instead of WPRIME, to stand for last time w-compoment(z)
!                  turbulent velocity. still use WPRIME here
!dwen(20090816): Starting from v4.9, HYSPLIT used ISEED,instead of RSEED, for the random seed that will be
!                  fed to PARVAR
!dwen(20090816): Starting from v4.7, HYSPLIT put WMIX(vertical turbulence profile) METZ, instead of METO
!dwen(20090817): remove KP from argument list, useless in PARDSP


!SUBROUTINE pardsp(umix,vmix,gdisx,gdisy,dtm,zmdl,zsfc,nlvl,wmix,zsg,xpos,  &
!                  ypos,zpos,vprime,wratio,uprime,hdwp,zndx,iseed)
SUBROUTINE pardsp(umix,vmix,gdisx,gdisy,dtm,zmdl,zsfc,nlvl,wmix,zsg,xpos,&
                  ypos,zpos,vprime,wprime,uprime,hdwp,zndx,iseed,        &
                  tlprev,stdwprev,zml1,densprev,tlnext,stdwnext,zml2,    &
!dwen(20090906)                  densnext,tlfrac,xposprev,yposprev,uuprev,uunext,       &
                  densnext,tlfrac,uuprev,uunext,       &
                  vvprev,vvnext,dmassprev,dmassnext,dmweight,rseed,      &
                  sigmaw,sigmau,samptt,samptt2,wwprev,veght,             &
                 convdur,zloc,ramsflg,tempnext,tempprev,rhfrnext,        &
                  rhfrprev,sphu,foot)
!************************************************ 

  IMPLICIT NONE

!dwen(20090822) ********************
   INCLUDE 'DEFARG2.INC'
!*************************************

!-------------------------------------------------------------------------------
! argument list variables 
!-------------------------------------------------------------------------------

  REAL,         INTENT(IN)     :: umix        ! u-component turbulence (m2/s2)
  REAL,         INTENT(IN)     :: vmix        ! v-component turbulence (m2/s2)
  REAL,         INTENT(IN)     :: gdisx       ! horizontal grid distance (m) 
  REAL,         INTENT(IN)     :: gdisy       ! horizontal grid distance (m) 
  REAL,         INTENT(IN)     :: dtm         ! advection time step (min)  
  REAL,         INTENT(IN)     :: zmdl        ! top of domain (m)          
  REAL,         INTENT(IN)     :: zsfc        ! terrain height (m)         
  INTEGER,      INTENT(IN)     :: nlvl        ! number of levels in subgrid
  REAL,         INTENT(IN)     :: wmix (:)    ! vertical turbulence (m2/s2)
  REAL,         INTENT(IN)     :: zsg  (:)    ! internal model sigma levels
  REAL,         INTENT(INOUT)  :: xpos        ! EW horizontal particle position
  REAL,         INTENT(INOUT)  :: ypos        ! NS horizontal particle position
  REAL,         INTENT(INOUT)  :: zpos        ! vertical particle position (sig)
  REAL,         INTENT(INOUT)  :: vprime      ! last v turbulent velocity (y)
!dwen(20090817)  REAL,         INTENT(INOUT)  :: wratio      ! last w turbulent velocity (z)
  REAL,         INTENT(INOUT)  :: wprime      ! last w turbulent velocity (z)
  REAL,         INTENT(INOUT)  :: uprime      ! last u turbulent velocity (x)
  INTEGER,      INTENT(IN)     :: hdwp        ! horizontal distribution index
  REAL,         INTENT(IN)     :: zndx        ! fractional vertical index  
  INTEGER,      INTENT(INOUT)  :: iseed       ! seed number for dispersion

!-------------------------------------------------------------------------------
! internally defined variables
!-------------------------------------------------------------------------------

  REAL                         :: aa,bb,cc          ! vertical grid polynomial
  REAL                         :: vscale            ! vert time scale (sec)
  REAL                         :: hscale            ! horz time scale (sec)
  REAL                         :: dt                ! internal time step (sec)
  REAL                         :: zx                ! vertical index
  REAL                         :: sigu,sigv,sigw    ! velocity sigma
  INTEGER                      :: kz,kb,kt          ! vertical indicies
  REAL                         :: uu,vv,ww          ! new velocity components
  REAL                         :: rauto             ! velocity autocorrelation
!dwen(20090817)  REAL                         :: wprime            ! vertical velocity turb  
  REAL                         :: delz              ! min vertical layer depth
  REAL                         :: delt              ! vertical time step 
  INTEGER                      :: knum              ! numb integration steps
  REAL                         :: fact              ! interp factor
  REAL                         :: vm                ! interp vertical mixing
  REAL                         :: sig2t,sig2b,sigdz ! velocity variances
  REAL                         :: zz,dist           ! temporary height variables

!********************************************
!dwen(20090318):adopt the following lines from pardsp of STILT version
!               change definition format
! JCL:(4/13/00)include to know how large to make arrays

! JCL:(6/13/00)logical for whether write out or not
      LOGICAL                  :: OUTYN

! JCL:(3/3/2000) needed to calculate diff in sigmaw between two internal timesteps
      REAL                     :: SIGWOLD

! JCL:vertical profile of Lagrangian timescale, std dev of W, air density, horizontal velocities
!      REAL*8 TL(NZM), STDW(NZM), DENS(NZM), UBAR(NZM), VBAR(NZM)
      real,allocatable         :: tl(:),stdw(:),dens(:),ubar(:),vbar(:)
! JCL:(11/14/02) temporary sigmal levels--needed if add a small 'sublayer' of low turbulence at level of zi
!      REAL*8 ZSGADD(NZM)
      REAL,allocatable         :: ZSGADD(:)

! JCL:(4/13/2000) vertical profile of Lagrangian timescale, std dev of W
!     before and after the mean advection timestep
      REAL, intent(inout)         :: TLPREV(:), STDWPREV(:)
      REAL, intent(inout)         :: TLNEXT(:), STDWNEXT(:)

! JCL:(2/23/2001) vertical profile of air density [kg/m3]
!     before and after the mean advection timestep
      REAL,intent(in)          :: DENSPREV(:), DENSNEXT(:)

! JCL:(5/9/01) vertical profile of horizontal wind [grid/min]
      REAL,intent(in)          :: UUPREV(:), VVPREV(:)
      REAL,intent(in)          :: UUNEXT(:), VVNEXT(:)
! JCL:(5/9/01) flag saying whether running backwards in time or not
!dwen(20090318)      LOGICAL BACK

! CHG:(12/04/01)add 'CONVDUR' set to 0 after conv. redistribution was done (input & output)
!dwen(20090824)      INTEGER,intent(inout)    :: CONVDUR
      real,   intent(inout)    :: CONVDUR

! JCL:(4/3/02)vertical profile of mass violation before & after mean advection
      REAL,intent(in)          :: DMASSPREV(:),DMASSNEXT(:)
!      REAL*8 DMASS(NZM)
      REAL,allocatable         :: DMASS(:)
! JCL:(4/3/02)weighting of particles changed by mass violation (1.0=no change)
      REAL,intent(inout)          :: DMWEIGHT

! CHG:(03/17/2004) vertical profile of virtual pot. temperature [K],
!     and relative humidity fraction
!     before and after the mean advection timestep
      REAL,intent(in)          :: TEMPPREV(:), TEMPNEXT(:)
      REAL,intent(in)          :: RHFRPREV(:), RHFRNEXT(:)
!      REAL,intent(in)          :: xposPREV
!      REAL,intent(in)          :: yposPREV

! CHG:(03/17/2004) vertical profile of dry air density [kg/m3], pressure [mbar],
!     Temperature [K], vapor pressure (mbar), sat. vap. pr. (mbar), and spec. hum. (g/g)
!      REAL*8 DENSDRY(NLVL),PRS(NLVL),TK(NLVL),E(NLVL),ES(NLVL),R(NLVL)  &
!     & ,ZADD(NLVL)
      REAL,allocatable :: DENSDRY(:),PRS(:),TK(:),E(:),ES(:),R(:),ZADD(:) 

! JCL:(4/13/00) interpolation factor in TIME
!dwen(20090318)      REAL*8 TF

! JCL:needed to test for whether particle sees vegetation
!dwen(20090318)      INTEGER SEEVEG

! JCL:fraction of Lagrangian timescale that will determine internal timestep
      REAL,intent(in)           :: TLFRAC
      REAL,intent(out)           :: sphu

! JCL:the random seed
      INTEGER,intent(in)        :: RSEED
!dwen(20090822): RAN3 is already defined in DEFARG1.INC 
!     REAL(KIND(1d0)), EXTERNAL :: RAN3


! JCL:standard deviation of vert&hor velocity that will be output to file
      REAL,intent(inout)          :: SIGMAW, SIGMAU

! JCL:(4/30/00) mixed-layer ht (interpolated in time betw ZML1&ZML2)
!         needed to implement particle reflection at top of boundary-layer
      REAL,intent(in)           :: ZML1,ZML2
!dwen(20090318):zmlprev,zmlnext,zml are useless
!      real                      :: ZMLPREV,ZMLNEXT,ZML

! JCL:(4/30/00) used to store vertical position before vertical dispersion
      REAL                      ::ZPOSOLD

! JCL:amount of time [min.] that particle 'sees' the ground
!      SAMPTT=(# particle touchdowns)*(deltat); deltat is timestep of subloop in PARDSP
      REAL,intent(out)          :: SAMPTT, FOOT
      REAL,intent(out)          :: SAMPTT2

! JCL:(5/3/00)variables to store altitudes [m] of lower & upper levels, and current altitude
!      'IFACT' is another variable for the interpolation factor
!dwen(20090318):zzb,zzt,zznow and ifact are useless
!      REAL*8 ZZB, ZZT, ZZNOW, IFACT

! JCL:(5/3/00)indices for levels (bottom,middle,top) used to determine smoothed gradient
!dwen(20090318):ks2b,ks2m and ks2t are useless
!      INTEGER KS2B,KS2M,KS2T

! JCL:vertical gradient in std dev of W--used to set max internal timestep DELT
!dwen(20090318):dsigwdz,dsigwdzn1 and dsigwdzn2 are useless
!      REAL*8 DSIGWDZ,DSIGWDZN1,DSIGWDZN2

! JCL:(5/16/00)transmission velocity thru interface and its squared value
!      REAL*8 WT,WT2
      REAL             :: WT

! JCL:(5/16/00)time [min] required for particle to reach interface
!dwen(20090318):tt is useless
!      REAL*8 TT

! JCL:(5/18/00) minutes since starting model run--required for prescribing time-varying
!        vertical profiles of SIGMAW & TL
!dwen(20090318):smin is uselss
!      INTEGER SMIN

! JCL:(5/18/00)'WWPREV' is output variable that stores WPRIME/SIGMAW
      REAL,intent(inout)   :: WWPREV

! JCL:(6/29/00)  the ht [m] below which a particle would be counted as 'seeing' grd vegetation
!         'VEGHTSIGMA' is simply conversion of VEGHT to sigma coordinate
      REAL,intent(in)   :: VEGHT
      real              :: VEGHTSIGMA

! CHG:(12/04/01)add 'ZLOC' limit of convection altitude (m)
      REAL,intent(in)   :: ZLOC

! CHG(09/16/03):add flag specifying whether data from RAMS or not
!dwen(20090318):ecmflg is useless
      LOGICAL, INTENT(IN) :: RAMSFLG
      integer             :: i,kzz,ktt,kk,kzm,nlvladd,kv,kublflag,     &
                             kupflag
      real                :: tdirsign,zmlsigma,veghtz,densint,zaddold, &
                             cumtt,prob,pp,zf,ubarint,vbarint,zzold,   &
                             zznorm,zscl
      real,    parameter  :: rdry = 287.04 ! dry air (J/Kg-k)
      real,    parameter  :: p2jm = 100  ! mb to j/m3
!      LOGICAL, INTENT(IN) :: RAMSFLG, ECMFLG
!********************************************

  COMMON /ZZTOKK/ AA,BB,CC
  COMMON /stblen/ vscale,hscale

!-------------------------------------------------------------------------------
! Configure subroutine interface argumment lists
!-------------------------------------------------------------------------------

  INTERFACE
    SUBROUTINE parvar (sigma,veloc,iseed) 
      IMPLICIT NONE
      REAL,         INTENT(IN )    :: sigma
      REAL,         INTENT(OUT)    :: veloc
      INTEGER,      INTENT(INOUT)  :: iseed
    END SUBROUTINE parvar
  END INTERFACE

!-------------------------------------------------------------------------------
      allocate(tl(nlvl+1),stdw(nlvl+1),dens(nlvl+1),ubar(nlvl+1),vbar(nlvl+1))
      allocate(ZSGADD(nlvl+1),dmass(nlvl+1))
      allocate(DENSDRY(nlvl),PRS(nlvl),TK(nlvl),E(nlvl),ES(nlvl), & 
               R(nlvl),ZADD(nlvl)) 
!******************************************
!dwen(20090318):adopted from PARDSP of STILT
! JCL:initialize VSCALE as 100 s (will be modified in PARDSP)
      VSCALE=100.0

!      xpos=xposprev
!      ypos=yposprev
  
! JCL(03/31/03):make sure that TL not = 0 when sigmaw is 0
      DO I=1,NLVL
         IF(STDWPREV(I).EQ.0.0)STDWPREV(I)=0.001
         IF(STDWNEXT(I).EQ.0.0)STDWNEXT(I)=0.001
         IF(TLPREV(I).EQ.0.0)TLPREV(I)=100.0
         IF(TLNEXT(I).EQ.0.0)TLNEXT(I)=100.0
      END DO
!*****************************************

! only compute dispersion for appropriate particle distributions (9/15/03)

  IF(HDWP.EQ.1.OR.HDWP.EQ.2.OR.HDWP.EQ.5.OR.HDWP.EQ.6) RETURN 
  
  DT=ABS(DTM*60.0)                    ! advection time step from min to sec
  ZX=ZNDX                             ! initial particle position index

!-------------------------------------------------------------------------------
! horizontal turbulence (only with 3D particle when hdwp=0)
!-------------------------------------------------------------------------------

  IF(HDWP.EQ.0)THEN 
!    horizontal velocity standard deviations
     SIGU=SQRT(UMIX)
     SIGV=SQRT(VMIX)

!*************************************
! JCL:   assign standard deviation of hor velocity to SIGMAU
         SIGMAU=SIGU

! CHG(09/16/03) normalize by density for RAMS
         IF(RAMSFLG)THEN
! get vertical index from ZNDX
!dwen(20090825)           KZZ=IDINT(ZNDX+1.0)
           KZZ=INT(ZNDX+1.0)
           SIGMAU=SIGMAU/DENSPREV(KZZ)
           SIGU=SIGU/DENSPREV(KZZ)
           SIGV=SIGV/DENSPREV(KZZ)
         END IF

! JCL:(7/1/02) no horizontal dispersion in simple 2-D (x&z) system
!        call random number generator for each component
!        where prime defines the velocity standard deviation
! JCL:(5/9/00)call function GASDEV to generate random variable instead
         UU=GASDEV(RSEED,SIGU)
         VV=GASDEV(RSEED,SIGV)

!*************************************
!    call random number generator for each component
!    where prime defines the velocity standard deviation
!dwen(20090318):use GASDEV instead 
!     CALL PARVAR(SIGU,UU,ISEED)
!     CALL PARVAR(SIGV,VV,ISEED)
!*************************************

!    autocorrelation function based upon time step and TL
     RAUTO=EXP(-DT/HSCALE)

!    compute new turbulent velocity component
     UPRIME=RAUTO*UPRIME+SQRT(1.0-RAUTO*RAUTO)*UU
     VPRIME=RAUTO*VPRIME+SQRT(1.0-RAUTO*RAUTO)*VV

!    adjust horizontal positions (m/s)(s)(grid/m)
     XPOS=XPOS+UPRIME*DT/GDISX
     YPOS=YPOS+VPRIME*DT/GDISY
  END IF

!-------------------------------------------------------------------------------
! vertical turbulence (only vertical particles hdwp = 0,3,4)
!-------------------------------------------------------------------------------
!******************************************************
!dwen(20090319):adopt the following codes according STILT
! JCL:(11/14/02)initialization of vertical index changes if vertical indices have changed; initialization done below
!     find the puff vertical index number for center position
! CHG&JCL change NINT(ZNDX) to INT() b/c if not, then
!     ZNDX may not be in between KB & KT
!     Ex: if ZNDX=4.7, then KZ=5 => KT=6 => KB=5
!         KZ=MAX0(1, MIN0(NLVL, INT(ZNDX) ))
!         KT=MIN0(KZ+1, NLVL)
!         KB=KT-1

! JCL:initialize SAMPTT to be 0
      SAMPTT=0.0
      FOOT=0.0
      SAMPTT2=0.0

! JCL:(3/3/00)initialize SIGWOLD before going into internal time loop
      SIGWOLD=-999.0

! JCL:(5/9/01)set direction of time, used in updating horizontal position below
      TDIRSIGN=1.0
!dwen(20090319):remove back and use DTM 
!      IF(BACK)TDIRSIGN=-1.0
      IF(int(DTM).LT.0)TDIRSIGN=-1.0

! JCL:(11/1/02)split internal loop into TWO parts, b/c averaging could remove the minimum in turbulence vars needed to prevent parti
!     1) using the "PREV" turbulence parameters before advection & 2) the other using the "NEXT" after advection
!************************** KTT=1 *******************************************
      DO KTT=1,2

! JCL:(11/1/02)no longer average turbulence vars
         IF(KTT.EQ.1)THEN
!           use VEGHT as fraction of mixed-layer ht if VEGHT <= 1.0; convert VEGHT from [m]=>modified sigma
            IF(VEGHT.LE.1.0)THEN
               VEGHTSIGMA=1.0-VEGHT*ZML1/(ZMDL-ZSFC)
            ELSE
               VEGHTSIGMA=1.0-VEGHT/(ZMDL-ZSFC)
            END IF

! JCL:(11/14/02)mixed-layer ht in SIGMA coordinates
            ZMLSIGMA=1.0-ZML1/(ZMDL-ZSFC)
! JCL:(11/14/02)determine vertical index of mixed-layer ht
! CHG(09/22/03) don't do this for RAMS
            IF(.NOT.RAMSFLG)THEN
            KK=1
! the following to cover the case that mixed layer height is below first model level
! this can happen even though zml>zlvls(1) in advpnt because of independent interpolations
! of zml, zsfc
            kzm=1  
!           DNINT used to prevent numerical errors
!dwen(20090825)            DO WHILE(DNINT(ZSG(KK)*1E4).GT.DNINT(ZMLSIGMA*1E4))
            DO WHILE(aNINT(ZSG(KK)*1E4).GT.aNINT(ZMLSIGMA*1E4))
               KZM=KK+1
               KK=KK+1
            END DO

!*************************************************
!dwen(20090903): avoid the dim of zsg out of bound
            kzm=min(kzm,nlvl-1)
!*************************************************

! JCL:(11/14/02)creating ADDTIONAL VERTICAL LEVEL right above mixed-layer
            NLVLADD=NLVL+1

!           turbulence profiles BEFORE advection
! JCL:(11/14/02)below level of mixed-layer ht, could just assign the previous profiles
            DO I=1,KZM
!               WRITE(45,*)I,ZSG(I),(1.0-ZSG(I))*ZMDL*ZMDL/(ZMDL-ZSFC)
               STDW(I)=STDWPREV(I)
               TL(I)=TLPREV(I)
               DENS(I)=DENSPREV(I)
               UBAR(I)=UUPREV(I)
               VBAR(I)=VVPREV(I)
               DMASS(I)=DMASSPREV(I)
               ZSGADD(I)=ZSG(I)
            END DO
! JCL:(11/14/02)create ADDITIONAL VERTICAL LEVEL that's very thin, so that particles don't get trapped in
!           low-turbulence layer above mixed-layer ht that's too thick
            ZSGADD(KZM+1)=ZSG(KZM)-0.05*(ZSG(KZM)-ZSG(KZM+1))
            STDW(KZM+1)=STDWPREV(KZM+1)
            TL(KZM+1)=TLPREV(KZM+1)
            DENS(KZM+1)=DENSPREV(KZM+1)
            UBAR(KZM+1)=UUPREV(KZM+1)
            VBAR(KZM+1)=VVPREV(KZM+1)
            DMASS(KZM+1)=DMASSPREV(KZM+1)
            DO I=KZM+2,NLVLADD
!               WRITE(45,*)I,ZSG(I),(1.0-ZSG(I))*ZMDL*ZMDL/(ZMDL-ZSFC)
               STDW(I)=STDWPREV(I-1)
               TL(I)=TLPREV(I-1)
               DENS(I)=DENSPREV(I-1)
               UBAR(I)=UUPREV(I-1)
               VBAR(I)=VVPREV(I-1)
               DMASS(I)=DMASSPREV(I-1)
               ZSGADD(I)=ZSG(I-1)
            END DO
                 !IF(.NOT.RAMSFLG)
            ELSE
            NLVLADD=NLVL
            ZSGADD(1:NLVLADD)=ZSG(1:NLVL)
            STDW(1:NLVLADD)=STDWPREV(1:NLVL)
            TL(1:NLVLADD)=TLPREV(1:NLVL)
            DENS(1:NLVLADD)=DENSPREV(1:NLVL)
            UBAR(1:NLVLADD)=UUPREV(1:NLVL)
            VBAR(1:NLVLADD)=VVPREV(1:NLVL)
            DMASS(1:NLVLADD)=DMASSPREV(1:NLVL)
            ENDIF

! CHG: (03/17/2004) get dry air density
!     to relate emissions in micro-mol/m2/s to changes in dry air mixing ratio in ppm
! HAVE: TEMP as virt. pot. temp, DENS as "wet" air density, and RHFR as relative humidity fraction
! WANT: dry air density
!     get p, T(K), e, esat, and rsat
!     get pressure (mbar), from virt. pot. temp (K) and DEN (kg/m3)
!     p=DEN*TK*RDRY/P2JM; TK=TT*(p/1000)**0.286
!     i.e. p=DEN*TT*(p/1000)**0.286*RDRY/P2JM
         PRS(1:NLVL)=(DENSPREV(1:NLVL)*TEMPPREV(1:NLVL)*(1.0/1000.0)    &
     &    **0.286 *RDRY/P2JM)**(1/(1.0-0.286))
!     get T (1:NLVL) from pot. temp. (1:NLVL)
         TK(1:NLVL)=TEMPPREV(1:NLVL)*(PRS(1:NLVL)/1000)**0.286

!     get saturation vap. pres. (mbar) from T (1:NLVL)
         WHERE(TK(1:NLVL).GE.273.15)
!     liquid
            ES(1:NLVL)=6.1078*exp(17.29694*(TK(1:NLVL)-273.16)          &
     &        /(TK(1:NLVL)-35.86))
         ELSEWHERE
!     ice
            ES(1:NLVL)=exp(23.33086-6111.72784/TK(1:NLVL)               &
     &       +0.15215*log(TK(1:NLVL)))
         END WHERE

!     from RH (rel.humidity fraction) to vapor pressure (mbar)
         E(1:NLVL)=RHFRPREV(1:NLVL)*ES(1:NLVL)

!     get spec. hum. r (g/g) from e (mbar) and p (mbar)
         R(1:NLVL)=0.622*E(1:NLVL)/(PRS(1:NLVL)-E(1:NLVL))

!     get dry air density
         DENSDRY(1:NLVL)=DENS(1:NLVL)*(1-R(1:NLVL))

! DONE getting dry air density profile

!****************  KTT is = to 2 (determine profiles for AFTER advection) ***************************
         ELSE
!           use VEGHT as fraction of mixed-layer ht if VEGHT <= 1.0; convert VEGHT from [m]=>modified sigma
            IF(VEGHT.LE.1.0)THEN
               VEGHTSIGMA=1.0-VEGHT*ZML2/(ZMDL-ZSFC)
            ELSE
               VEGHTSIGMA=1.0-VEGHT/(ZMDL-ZSFC)
            END IF

! JCL:(11/14/02)mixed-layer ht in SIGMA coordinates
            ZMLSIGMA=1.0-ZML2/(ZMDL-ZSFC)
! JCL:(11/14/02)determine vertical index of mixed-layer ht
! CHG(09/22/03) don't do this for RAMS
            IF(.NOT.RAMSFLG)THEN
            KK=1
! the following to cover the case that mixed layer height is below first model level
! this can happen even though zml>zlvls(1) in advpnt because of independent interpolations
! of zml, zsfc
            kzm=1  
!           DNINT used to prevent numerical errors
!dwen(20090825)            DO WHILE(DNINT(ZSG(KK)*1E4).GT.DNINT(ZMLSIGMA*1E4))
            DO WHILE(aNINT(ZSG(KK)*1E4).GT.aNINT(ZMLSIGMA*1E4))
               KZM=KK+1
               KK=KK+1
            END DO

!*************************************************
!dwen(20090903): avoid the dim of zsg out of bound
            kzm=min(kzm,nlvl-1)
!*************************************************

! JCL:(11/14/02)creating ADDTIONAL VERTICAL LEVEL right above mixed-layer
            NLVLADD=NLVL+1

!           turbulence profiles BEFORE advection
! JCL:(11/14/02)below level of mixed-layer ht, could just assign the previous profiles
            DO I=1,KZM
!              WRITE(45,*)I,ZSG(I),(1.0-ZSG(I))*ZMDL*ZMDL/(ZMDL-ZSFC)
               STDW(I)=STDWNEXT(I)
               TL(I)=TLNEXT(I)
               DENS(I)=DENSNEXT(I)
               UBAR(I)=UUNEXT(I)
               VBAR(I)=VVNEXT(I)
               DMASS(I)=DMASSNEXT(I)
               ZSGADD(I)=ZSG(I)
            END DO
! JCL:(11/14/02)create ADDITIONAL VERTICAL LEVEL that's very thin, so that particles don't get trapped in
!           low-turbulence layer above mixed-layer ht that's too thick
            ZSGADD(KZM+1)=ZSG(KZM)-0.05*(ZSG(KZM)-ZSG(KZM+1))
            STDW(KZM+1)=STDWNEXT(KZM+1)
            TL(KZM+1)=TLNEXT(KZM+1)
            DENS(KZM+1)=DENSNEXT(KZM+1)
            UBAR(KZM+1)=UUNEXT(KZM+1)
            VBAR(KZM+1)=VVNEXT(KZM+1)
            DMASS(KZM+1)=DMASSNEXT(KZM+1)
            DO I=KZM+2,NLVLADD
!              WRITE(45,*)I,ZSG(I),(1.0-ZSG(I))*ZMDL*ZMDL/(ZMDL-ZSFC)
               STDW(I)=STDWNEXT(I-1)
               TL(I)=TLNEXT(I-1)
               DENS(I)=DENSNEXT(I-1)
               UBAR(I)=UUNEXT(I-1)
               VBAR(I)=VVNEXT(I-1)
               DMASS(I)=DMASSNEXT(I-1)
               ZSGADD(I)=ZSG(I-1)
            END DO
                 !IF(.NOT.RAMSFLG)
            ELSE
            NLVLADD=NLVL
            ZSGADD(1:NLVLADD)=ZSG(1:NLVL)
            STDW(1:NLVLADD)=STDWNEXT(1:NLVL)
            TL(1:NLVLADD)=TLNEXT(1:NLVL)
            DENS(1:NLVLADD)=DENSNEXT(1:NLVL)
            UBAR(1:NLVLADD)=UUNEXT(1:NLVL)
            VBAR(1:NLVLADD)=VVNEXT(1:NLVL)
            DMASS(1:NLVLADD)=DMASSNEXT(1:NLVL)
            ENDIF

! CHG: (03/17/2004) get dry air density
!     to relate emissions in micro-mol/m2/s to changes in dry air mixing ratio in ppm
! HAVE: TEMP as virt. pot. temp, DENS as "wet" air density, and RHFR as relative humidity fraction
! WANT: dry air density
!     get p, T(K), e, esat, and rsat
!     get pressure (mbar), from virt. pot. temp (K) and DEN (kg/m3)
!     p=DEN*TK*RDRY/P2JM; TK=TT*(p/1000)**0.286
!     i.e. p=DEN*TT*(p/1000)**0.286*RDRY/P2JM
         PRS(1:NLVL)=(DENSNEXT(1:NLVL)*TEMPNEXT(1:NLVL)*(1.0/1000.0)    &
     &   **0.286*RDRY/P2JM)**(1/(1.0-0.286))
!     get T (1:NLVL) from pot. temp. (1:NLVL)
         TK(1:NLVL)=TEMPNEXT(1:NLVL)*(PRS(1:NLVL)/1000)**0.286

!     get saturation vap. pres. (mbar) from T (1:NLVL)
         WHERE(TK(1:NLVL).GE.273.15)
!     liquid
            ES(1:NLVL)=6.1078*exp(17.29694*(TK(1:NLVL)-273.16)          &
     &        /(TK(1:NLVL)-35.86))
         ELSEWHERE
!     ice
            ES(1:NLVL)=exp(23.33086-6111.72784/TK(1:NLVL)               &
     &       +0.15215*log(TK(1:NLVL)))
         END WHERE

!     from RH (rel.humidity fraction) to vapor pressure (mbar)
         E(1:NLVL)=RHFRNEXT(1:NLVL)*ES(1:NLVL)

!     get spec. hum. r (g/g) from e (mbar) and p (mbar)
         R(1:NLVL)=0.622*E(1:NLVL)/(PRS(1:NLVL)-E(1:NLVL))

!     get dry air density
         DENSDRY(1:NLVL)=DENS(1:NLVL)*(1-R(1:NLVL))

! DONE getting dry air density profile


!        IF(KTT.EQ.1)THEN
         END IF

! JCL:(11/14/02)initialize vertical index
! JEL:(10/01/07)restrict the search to valid range 
!     This prevents errors for particles above the top level.
!     Could be done prior to calling pardsp (using is_off_grid?)
         KZ=1
         DO WHILE(ZSGADD(KZ+1).GT.ZPOS .AND. KZ.LE.NLVL)
           !CHG i.e. if particle below HGT(2), KB=1 and KT=2; also if particle below HGT(1)
           KZ=KZ+1
         END DO
         KT=KZ+1
         IF(KZ.EQ.NLVLADD)KT=NLVLADD
         KB=KT-1
! CHG(09/16/03) for RAMS use also different vertical index to extract scalars
         IF(RAMSFLG)THEN
           KZ=1
           DO WHILE(ZSGADD(KZ).GT.ZPOS .AND. KZ.LE.NLVL)
             KZ=KZ+1
           END DO
           KT=KZ
           IF(KT.GT.NLVLADD)KT=NLVLADD
                   !CHG i.e. if particle below HGT(1), KB=0 and KT=1
           KB=KT-1
         END IF

! CHG:(03/17/2004) get dry air density between surface and veght
! integrate density up to vertical index corresponding to VEGHTSIGMA, then add remaining fraction of layer
! get Z profile and veght in meters above ground
         ZADD(1:NLVL)=(1.0-ZSGADD(1:NLVL))*(ZMDL-ZSFC)
         VEGHTZ=(1.0-VEGHTSIGMA)*(ZMDL-ZSFC)
         KV=1
         DENSINT=0.0
         ZADDOLD=0.0
         IF(RAMSFLG)THEN
           DO WHILE(ZADD(KV).LT.VEGHTZ)
             DENSINT=DENSINT+DENSDRY(KV)*(ZADD(KV)-ZADDOLD)
             ZADDOLD=ZADD(KV)
             KV=KV+1
           END DO
              ! for other than RAMS
         ELSE
           DO WHILE(ZADD(KV+1).LT.VEGHTZ)
             DENSINT=DENSINT+DENSDRY(KV)*(ZADD(KV)-ZADDOLD)
             ZADDOLD=ZADD(KV)
             KV=KV+1
           END DO
         END IF
         DENSINT=DENSINT+DENSDRY(KV)*(VEGHTZ-ZADDOLD)
! JCL:(5/18/00)initialize WPRIME with value stored for this particle from last PARDSP,
!       multiplied by new SIGMAW--required to have proper treatment of step changes in
!       turbulence statistics from one timestep to the next (Thomson et al.[1997])
! CHG(09/16/03) not for RAMS, there different vertical index required
      IF(.NOT.RAMSFLG)WPRIME=STDW(KB)*WWPREV
      IF(RAMSFLG)WPRIME=STDW(KT)*WWPREV

! JCL:(9/5/02) initialize cumulative time elapsed
      CUMTT=0.0

! JCL:(11/1/02)split internal loop into TWO parts, b/c averaging could remove the minimum in turbulence vars needed to prevent parti
! JCL:(9/5/02) new timestep implementation
!     DO WHILE(CUMTT.LT.DABS(DT))
!dwen(20090817) use DTM instead of DT
!dwen(20090825)      DO WHILE(CUMTT.LT.DABS(DTm)/2.0)

      DO WHILE(CUMTT.LT.ABS(DTm)/2.0)
! JCL:(5/16/00)no longer conduct vertical interpolation--simply take value on lower level
! CHG(09/16/03) not for RAMS, there different vertical index required
         IF(.NOT.RAMSFLG)THEN
           SIGW=STDW(KB)
           VSCALE=TL(KB)
         ELSE
           SIGW=STDW(KT)
           VSCALE=TL(KT)
         END IF

! JCL:(9/5/02)new timestep scheme
! dwen(20090817) used DTM instead of DT
!dwen(20090825)         DELT=DMIN1(DABS(DTm)/2.0-CUMTT,TLFRAC*(VSCALE/60.0))

         DELT=MIN(ABS(DTm)/2.0-CUMTT,TLFRAC*(VSCALE/60.0))

! JCL:   (5/9/00)call new function GASDEV for calculating random fluctuations
         WW=GASDEV(RSEED,SIGW)

! JCL:   store result in 'SIGMAW', which will be written to output file
         SIGMAW=SIGW

! CHG:   store result in 'SAMPTT2', which will be written to output file
         SAMPTT2=VSCALE

!        autocorrelation function based upon time step and TL
! JCL:   shouldn't RAUTO be determined by the INTERNAL timestep, not the
!        external timestep?; move RAUTO calculation into internal loop
!dwen(20090825)         RAUTO=DEXP(-60.0*DELT/VSCALE)
         RAUTO=EXP(-60.0*DELT/VSCALE)

!        new turbulent vertical velocity
!dwen(20090825)         WPRIME=RAUTO*WPRIME+DSQRT(1.0-RAUTO*RAUTO)*WW
         WPRIME=RAUTO*WPRIME+SQRT(1.0-RAUTO*RAUTO)*WW

! JCL:(4/30/00)store position before vertical dispersion
         ZPOSOLD=ZPOS

! JCL:(9/5/02)whether or not particle is in the Unresolved Boundary Layer (UBL)
         KUBLFLAG=0
! CHG(09/16/03) use different index for RAMS
         IF(KB.EQ.1.AND.(.NOT.RAMSFLG))KUBLFLAG=1
         IF(KB.EQ.0.AND.RAMSFLG)KUBLFLAG=1
! JCL:(9/5/02) time necessary to transport particle to model level interface, given the wprime
         IF(WPRIME.GT.0.0)THEN
            KUPFLAG=1
!           TAU=(ZPOS-ZSG(KT))/(60.0*(WPRIME/(ZMDL-ZSFC)))
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!dwen(20090825)            DIST=DABS((ZPOS-ZSGADD(KT))*(ZMDL-ZSFC)/60.0)
            DIST=ABS((ZPOS-ZSGADD(KT))*(ZMDL-ZSFC)/60.0)
         END IF
! CHG:(9/17/02) use LE instead of LT, to include 0.0
         IF(WPRIME.LE.0.0)THEN
            KUPFLAG=0
! JCL:      if particle is near the ground, in the Unresolved Boundary Layer (UBL)
            IF(KUBLFLAG.EQ.1)THEN
!              TAU=(ZPOS-1.0)/(60.0*(WPRIME/(ZMDL-ZSFC)))
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!dwen(20090825)               DIST=DABS((ZPOS-1.0)*(ZMDL-ZSFC)/60.0)
               DIST=ABS((ZPOS-1.0)*(ZMDL-ZSFC)/60.0)
            ELSE
!              TAU=(ZPOS-ZSG(KB))/(60.0*(WPRIME/(ZMDL-ZSFC)))
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!dwen(20090825)               DIST=DABS((ZPOS-ZSGADD(KB))*(ZMDL-ZSFC)/60.0)
               DIST=ABS((ZPOS-ZSGADD(KB))*(ZMDL-ZSFC)/60.0)
            END IF
         END IF

! JCL:(9/5/02)when particle near ground or model top, not have to implement Thomson [1997] reflection/transmission
         IF((KUBLFLAG.EQ.1.AND.WPRIME.LT.0.0).OR.                       &
     &                   (KT.EQ.NLVLADD.AND.WPRIME.GT.0.0))THEN
! JCL:(9/5/02) particle would NOT hit ground or model top, so simply advect
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!           IF(TAU.GT.DELT)THEN
!dwen(20090825)            IF(DIST.GE.DABS(DELT*WPRIME))THEN
            IF(DIST.GE.ABS(DELT*WPRIME))THEN
!              vertical displacement (m/s => sigma/min => sigma)
               ZPOS=ZPOS-(WPRIME/(ZMDL-ZSFC))*DELT*60.0

! JCL:      particle WOULD reach ground or top: implement reflection, but move particle only to ground or model top
            ELSE
! JCL:         sign of turbulent velocity reversed when reflection occurs
               WPRIME=(-1.0*WPRIME)
! JCL:         particle spent timestep TAU moving to ground or model top
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!              DELT=TAU
!dwen(20090825)               DELT=DABS(DIST/WPRIME)
               DELT=ABS(DIST/WPRIME)
! JCL:         particle position moved to either ground or model top
               IF(KUBLFLAG.EQ.1)ZPOS=1.0
               IF(KUBLFLAG.EQ.0)ZPOS=ZSGADD(NLVLADD)
            END IF

! JCL:   particle NOT approaching ground or model top
         ELSE

! JCL:(9/5/02) particle wouldn't pass through interface, so simply advect
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!           IF(TAU.GT.DELT)THEN
!dwen(20090825)            IF(DIST.GE.DABS(DELT*WPRIME))THEN
            IF(DIST.GE.ABS(DELT*WPRIME))THEN
!              vertical displacement (m/s => sigma/min => sigma)
               ZPOS=ZPOS-(WPRIME/(ZMDL-ZSFC))*DELT*60.0

!           particle WOULD reach interface: implement Thomson [1997] reflection/transmission scheme
            ELSE
!              ########################################################
!              particle crossed interface from *******BOTTOM=>UP*******
               IF(KUPFLAG.EQ.1)THEN
! JCL:(8/5/02)!!! new reflection/transmission scheme in which particles are transmitted at random
!                 probability of transmission; note that could be > 1.0 => then particle ALWAYS transmitted
! CHG(09/22/03) for RAMS get scalars from proper level
                  IF(.NOT.RAMSFLG)PROB=DENS(KT)*STDW(KT)/               &
     &                                 (DENS(KB)*STDW(KB))
                  IF(RAMSFLG)PROB=DENS(KT+1)*STDW(KT+1)/                &
     &                            (DENS(KB+1)*STDW(KB+1))
!                 uniform random deviate between 0.0 and 1.0
                  PP=RAN3(RSEED)
! JCL:(8/5/02)    particle TRANSMITTED
                  IF(PP.LE.PROB)THEN
! JCL:               move particle to location of interface
                     ZPOS=ZSGADD(KT)
! JCL:               only used up amount of time it takes for particle to get to interface, so update timestep
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!                    DELT=TAU
!dwen(20090825)                     DELT=DABS(DIST/WPRIME)
                     DELT=ABS(DIST/WPRIME)
! JCL:               rescaled velocity
! CHG(09/22/03) for RAMS get scalars from proper level
                     IF(.NOT.RAMSFLG)WT=WPRIME*STDW(KT)/STDW(KB)
                     IF(RAMSFLG)WT=WPRIME*STDW(KT+1)/STDW(KB+1)
! JCL:               update velocity
                     WPRIME=WT
! JCL:               update vertical index
                     KT=KT+1
                     KB=KB+1
! JCL:(9/5/02)    particle REFLECTED, but only move to location of interface
                  ELSE
! JCL:               move to location of interface
                     ZPOS=ZSGADD(KT)
! JCL:               particle spent timestep TAU moving to interface
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME  unit:[min m/s]  !!!!!
!                    DELT=TAU
!dwen(20090825)                     DELT=DABS(DIST/WPRIME)
                     DELT=ABS(DIST/WPRIME)
! JCL:               sign of turbulent velocity reversed when reflection occurs
                     WPRIME=(-1.0*WPRIME)
                  END IF
               END IF

!              particle crossed interface from *******TOP=>BOTTOM******
               IF(KUPFLAG.EQ.0)THEN
! JCL:(8/5/02)!!! new reflection/transmission scheme in which particles are transmitted at random
!                 probability of transmission; note that could be > 1.0 => then particle ALWAYS transmitted
! CHG(09/22/03) for RAMS get scalars from proper level
                  IF(.NOT.RAMSFLG)PROB=DENS(KB-1)*STDW(KB-1)/           &
     &                                 (DENS(KB)*STDW(KB))
                  IF(RAMSFLG)PROB=DENS(KB)*STDW(KB)/                    &
     &                                 (DENS(KT)*STDW(KT))
!                 uniform random deviate between 0.0 and 1.0
                  PP=RAN3(RSEED)
! JCL:(8/5/02)    particle TRANSMITTED
                  IF(PP.LE.PROB)THEN
! JCL:               move particle to location of interface
                     ZPOS=ZSGADD(KB)
! JCL:               only used up amount of time it takes for particle to get to interface, so update timestep
! CHG(11/04/02) avoid dividing by 0 (when WPRIME==0), use DIST =TAU*WPRIME
!                    DELT=TAU
!dwen(20090825)                     DELT=DABS(DIST/WPRIME)
                     DELT=ABS(DIST/WPRIME)
! JCL:               rescaled velocity
! CHG(09/22/03) for RAMS get scalars from proper level
                     IF(.NOT.RAMSFLG)WT=WPRIME*STDW(KB-1)/STDW(KB)
                     IF(RAMSFLG)WT=WPRIME*STDW(KB)/STDW(KT)
! JCL:               update velocity
                     WPRIME=WT
! JCL:               update vertical index
                     KT=KT-1
                     KB=KB-1

! JCL:(8/5/02)    particle REFLECTED, but only move to location of interface
                  ELSE
! JCL:               move to location of interface
                     ZPOS=ZSGADD(KB)
! JCL:               particle spent timestep TAU moving to interface
! CHG(11/04/02) avoid dividing by 0 when WPRIME==0; use DIST=TAU*WPRIME
!                    DELT=TAU
!dwen(20090825)                     DELT=DABS(DIST/WPRIME)
                     DELT=ABS(DIST/WPRIME)
! JCL:               sign of turbulent velocity reversed when reflection occurs
                     WPRIME=(-1.0*WPRIME)

                  END IF
               END IF

! JCL:      whether to implement reflection/transmission scheme
            END IF
! JCL:   whether particle near model bottom/top
         END IF

! JCL:(5/9/01) update horizontal position with horizontal wind as particle is transported by vertical turbulence
!         XPOS=XPOS+UBAR(KB)*TDIRSIGN*DELT
!         YPOS=YPOS+VBAR(KB)*TDIRSIGN*DELT
! JCL:(4/4/02) previously had DISCRETE wind profiles to advect particle, but now
!              LINEARLY INTERPOLATE vertically to get wind profile
! CHG(09/16/03) not for RAMS
         IF(.NOT.RAMSFLG)THEN
!        interpolation factor
         ZF=(ZSGADD(KB)-ZPOS)/(ZSGADD(KB)-ZSGADD(KT))
         UBARINT=UBAR(KB)+ZF*(UBAR(KT)-UBAR(KB))
         VBARINT=VBAR(KB)+ZF*(VBAR(KT)-VBAR(KB))
! JCL:(4/4/02)impose NO-SLIP boundary-condition if below 1st level
!        note that if particle is below 1st level, it is still assigned KB=1
         IF(ZPOS.GT.ZSGADD(1))THEN
            ZF=(1.0-ZPOS)/(1.0-ZSGADD(KB))
            UBARINT=0.0+ZF*(UBAR(KB)-0.0)
            VBARINT=0.0+ZF*(VBAR(KB)-0.0)
         END IF
         END IF
         IF(RAMSFLG)THEN
           !no interpolation, but convert fluxes to velocities
           UBARINT=UBAR(KT)
           VBARINT=VBAR(KT)
         END IF

! JCL:(7/2/02) previously linearly interpolated in SIGMA coordinates, but really should linearly
!        interpolate in INDEX COORDINATES (to be consistent with interpolation routine in ADVIEC)
!        ZX=(-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)

         XPOS=XPOS+UBARINT*TDIRSIGN*DELT
         YPOS=YPOS+VBARINT*TDIRSIGN*DELT

! JCL:(4/5/02)calculate the total mass violation experienced by particle (fraction of mass in gridcell)
! CHG(09/22/03) for RAMS get scalars from proper level
         IF(.NOT.RAMSFLG)THEN
         IF(ZPOS.LT.ZSGADD(1))THEN
!           if particle is above lowest model level
            DMWEIGHT=DMWEIGHT*(1.0+DMASS(KT)*DELT)
         ELSE
!           particle is near ground
            DMWEIGHT=DMWEIGHT*(1.0+DMASS(1)*DELT)
         END IF
              !for RAMS
         ELSE
           DMWEIGHT=DMWEIGHT*(1.0+DMASS(KT)*DELT)
         END IF

! JCL:(6/29/00)calculate how much time particle spent within specified ht above ground
!                so to calculate vegetation influence; if within specified ht, then increment SAMPTT
! CHG:(03/17/2004) get residence time multiplied with 1/(#moles dry air density) between surface and veght
!                to have FOOT as "sensitivity" of mixing ratio changes (in ppm) to fluxes (in micro-moles/m2/s)
!                multiplying FOOT with fluxes gives then mixing ratio increments in PPM

         IF(ZPOS.GT.VEGHTSIGMA)THEN
           SAMPTT=SAMPTT+DELT
           FOOT=FOOT+DELT*0.02884*60.0/DENSINT
         END IF
!           28.84e-3 is mass of a mol of air in kg, 60 is to get from minuts to seconds

! JCL:   (3/3/00)update the 'old' sigmaw at end of internal time loop
         SIGWOLD=SIGW

! JCL:   (9/5/02) increment accumulated time
         CUMTT=CUMTT+DELT

!     time step loop
      END DO   !end while CUMTT

! JCL:(5/18/00)store normalized turbulent velocity
! CHG(09/22/03) for RAMS get scalars from proper level
      IF(.NOT.RAMSFLG)WWPREV=WPRIME/STDW(KB)
      IF(RAMSFLG)WWPREV=WPRIME/STDW(KT)

! CHG (03/17/2004) provide SPHU as output
      IF(.NOT.RAMSFLG)SPHU=R(KB)
      IF(RAMSFLG)SPHU=R(KT)

      END DO  !end KTT

! CHG:(12/05/01) redistribute particles every analysis time
!     into well mixed profile up to ZLOC
!     do it only when ZLOC above zi, and when particle is below ZLOC
!     Get old altitude in m
      ZZOLD=(1.0-ZPOS)*(ZMDL-ZSFC)

!dwen(20090825)      IF (CONVDUR > 0 .AND. ZLOC > DMAX1(ZML1,ZML2) .AND. ZLOC > ZZOLD) THEN
      IF (CONVDUR > 0 .AND. ZLOC > MAX1(ZML1,ZML2) .AND. ZLOC > ZZOLD) THEN
!        get random number between 0,1
         ZZNORM=RAN3(RSEED)
!        transform into altitude between surface height ZSNC and ZLOC
!        integr(ro.normalized) from ZSFC to ZZ is set equal to ZZNORM
!        use scale height ZSCL of 9200 m, good to 10 km
         ZSCL=9200.0
!dwen(20090825)         ZZ=ZSFC-ZSCL*DLOG(1-ZZNORM*(1-DEXP((ZSFC-ZLOC)/ZSCL)))
         ZZ=ZSFC-ZSCL*LOG(1-ZZNORM*(1-EXP((ZSFC-ZLOC)/ZSCL)))

!        get altitude in sigma from AGL m
         ZPOSOLD=1.0-(ZZ/(ZMDL-ZSFC))
        ZPOS=1.0-(ZZ/(ZMDL-ZSFC))
      END IF

      CONVDUR=0

!******************************************************
!dwen(20090319):comment out the following lines because of using  
!                  the above lines
! determine the required vertical dispersion time step (sec)
! by computing the minimum time step for vertical stability
!  DELZ=(ZMDL-ZSFC)*(ZSG(1)-ZSG(2))
!  DELT=MIN(DT,0.125*DELZ*DELZ/MAXVAL(WMIX)/VSCALE)
!  DELT=MAX(1,NINT(DELT))
!
!! round down time step to even multiple of DT
!  DO WHILE (MOD(INT(DT),INT(DELT)).NE.0.AND.INT(DELT).GT.1)
!     DELT=DELT-1.0
!  END DO
!
!! autocorrelation function based upon time step and TL
!! delt replaced dt - RRD:1/3/00
!  RAUTO=EXP(-DELT/VSCALE)
!
!! go through iterations to match external DT
!
!tstep : DO KNUM=1,NINT(DT/DELT)
!
!!    layer indicies
!     KZ=MAX(1, MIN(NLVL, INT(ZX)))
!     KT=MIN(KZ+1, NLVL)
!     KB=KT-1
!
!!    mixing at particle position
!     FACT=MAX(0.0, MIN(1.0,(ZSG(KB)-ZPOS)/(ZSG(KB)-ZSG(KT))))
!     VM=(WMIX(KT)-WMIX(KB))*FACT+WMIX(KB)
!
!!    vertical velocity standard deviation (m/s)
!     SIGW=SQRT(VM)
!     CALL PARVAR(SIGW,WW,ISEED)
!
!!    gradient of velocity variance
!     SIG2T=SQRT(WMIX(KT))
!     SIG2B=SQRT(WMIX(KB))
!     SIGDZ=(SIG2T-SIG2B) / ((ZSG(KB)-ZSG(KT)) * (ZMDL-ZSFC))
!
!!    new turbulent vertical velocity
!!    WPRIME=RAUTO*WPRIME+SQRT(1.0-RAUTO*RAUTO)*WW
!!    low turbulence convergence correction
!!    WPRIME=WPRIME+VSCALE*(1.0-RAUTO)*SIGDZ
!
!!    correction for non-homogeneous turbulence: Chock and Winkler (1994)
!!    J. Geophys. Res., Vol 99, No D1, Pages 1033-1041
!
!     WPRIME=WRATIO*SIGW
!     WRATIO=RAUTO*WRATIO+SQRT(1.0-RAUTO*RAUTO)*WW/SIGW  &
!           +VSCALE*(1.0-RAUTO)*SIGDZ
!     WPRIME=WRATIO*(SIGW+(WPRIME*DELT*SIGDZ))
!
!!    vertical displacement (m/s)(1/m)(s) => sigma
!     ZPOS=ZPOS-(WPRIME/(ZMDL-ZSFC))*DELT
!
!!    full reflection at the bottom
!     IF(ZPOS.GT.1.0)THEN
!        ZPOS=2.0-ZPOS
!!       reverse turbulent component
!        WRATIO=-WRATIO
!     END IF
!
!!    full reflection at the top
!     IF(ZPOS.LT.ZSG(NLVL))THEN
!        ZPOS=2.0*ZSG(NLVL)-ZPOS
!        WRATIO=-WRATIO
!     END IF
!
!!    recompute the particle vertical index number
!!    from integer index based upon quadratic relation between
!!    height (or sigma) and array index position
!
!!     height agl when zsfc=0
!      ZZ=ZMDL*(1.0-MIN(1.0,ZPOS))
!
!!     convert to vertical index using equation for zsfc=0
!      DIST=(BB*BB-4.0*AA*(CC-ZZ))
!      IF(DIST.GE.0.0)THEN
!         ZX=(-BB+SQRT(DIST))/(2.0*AA)
!      ELSE
!         ZX=1.0
!      END IF
!
!  END DO tstep

      deallocate(tl,stdw,dens,ubar,vbar)
      deallocate(ZSGADD,dmass)
      deallocate(DENSDRY,PRS,TK,E,ES,R,ZADD) 

END SUBROUTINE pardsp
