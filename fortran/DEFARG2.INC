!-------------------------------------------------------------------------------
! Configure subroutine interface argument lists for subroutine ADVPNT     
! Last Revised: 07 Sep 2001 (RRD) - simultaneous multiple meteorology
!               26 Feb 2002 (RRD) - downward shortwave flux: metinp, advmet
!               21 May 2002 (RRD) - divergence options: metwnd
!               23 Jul 2002 (RRD) - terrain correction eta: metwnd
!               10 Apr 2003 (RRD) - advmet, adviec temporal interpolation
!               16 Sep 2003 (RRD) - added routine advsfc
!               10 Nov 2003 (RRD) - velocity variances through metinp
!               25 Oct 2004 (RRD) - pass through tratio
!               12 Oct 2005 (RRD) - added routine adviso
!               08 Mar 2006 (RRD) - static stability (prfcom, advmet)
!               22 May 2006 (RRD) - mixed layer depth (metinp)
!               21 Nov 2006 (RRD) - day night tke partition
!               04 Jun 2008 (RRD) - additional mixing options prfcom
!-------------------------------------------------------------------------------

INTERFACE

      SUBROUTINE ADVMETGRELL(METz,X,Y,CFLXUP1,CFLXUP2,CFLXDN1,DFLXUP1,DFLXUP2,    &
                       EFLXUP1,EFLXUP2,DFLXDN1,EFLXDN1,TKEG,NLVL,GLOBAL,NXP,NYP)

      IMPLICIT NONE 

      include 'DEFMETO.INC'
!------------------------------------------
!dwen(20090824)  TYPE(aset),INTENT(OUT)   :: meto       ! surface advection variables
  TYPE(bset),INTENT(OUT)   :: metz(:)       ! surface advection variables
  real,      intent(in)    :: x,y        
  real,      intent(in)    :: cflxup1(:,:,:),cflxup2(:,:,:)
  real,      intent(in)    :: cflxdn1(:,:,:),dflxup1(:,:,:)
  real,      intent(in)    :: dflxup2(:,:,:),dflxdn1(:,:,:)
  real,      intent(in)    :: eflxup1(:,:,:),eflxup2(:,:,:)
  real,      intent(in)    :: eflxdn1(:,:,:),tkeg(:,:,:)
  integer,   intent(in)    :: nlvl,nxp,nyp 
  logical,   intent(in)    :: global

   end  SUBROUTINE ADVMETGRELL

       SUBROUTINE CGRELL(CFXUP1,CFXUP2,CFXDN1,DFXUP1,DFXUP2,EFXUP1,     &
     &                   EFXUP2,DFXDN1,EFXDN1,RAUP1,RAUP2,RADN1,AREA,   &
     &                   AREAPRU,AREAPRD,DENS,NLVL,ZPROFM,Z1,Z2,        &
     &                   ICNDX1,ICNDX2,ZX,DT,BACK)

      implicit  NONE
     
      real,   intent(in)    :: cfxup1(:),cfxup2(:),cfxdn1(:)
      real,   intent(in)    :: dfxup1(:),dfxup2(:),dfxdn1(:)
      real,   intent(in)    :: efxup1(:),efxup2(:),efxdn1(:)
      real,   intent(in)    :: raup1(:),raup2(:),radn1(:)
      real,   intent(in)    :: area 
      real,   intent(inout) :: areapru,areaprd 
      real,   intent(in)    :: dens(:)
      integer,intent(in)    :: nlvl
      real,   intent(in)    :: zprofm(:) 
      real,   intent(in)    :: z1
      real,   intent(out)   :: z2
      integer,   intent(inout) :: icndx1,zx
      integer,   intent(out)   :: icndx2
      real,   intent(in)    :: dt
      logical, intent(in)   :: back
      
   end  SUBROUTINE CGRELL

!-------------------------------------------------------------------------------
  SUBROUTINE METPOS(BACK,XP,YP,JET,NGRD,NTIM,FTIME,MTIME,POINT,OFFG,   &
                    KGC,KGX,KT1,KT2)
  IMPLICIT NONE
  LOGICAL,   INTENT(IN)    :: back           ! defines backward integration
  REAL,      INTENT(INOUT) :: xp,yp          ! x,y grid particle position
  INTEGER,   INTENT(IN)    :: jet            ! elapsed time (minutes)
  INTEGER,   INTENT(IN)    :: ngrd           ! number of meteo grids   
  INTEGER,   INTENT(IN)    :: ntim           ! number of meteo times   
  INTEGER,   INTENT(IN)    :: ftime (:,:)    ! acc min of current array data
  INTEGER,   INTENT(OUT)   :: mtime (2)      ! acc min of requested input data
  INTEGER,   INTENT(OUT)   :: point (2)      ! array index pointer 
  LOGICAL,   INTENT(OUT)   :: offg           ! off entire grid flag
  INTEGER,   INTENT(IN)    :: kgc            ! current particle grid number
  INTEGER,   INTENT(OUT)   :: kgx            ! new particle grid number
  INTEGER,   INTENT(OUT)   :: kt1            ! grid number time last    
  INTEGER,   INTENT(OUT)   :: kt2            ! grid number time next    
  END SUBROUTINE metpos
!-------------------------------------------------------------------------------
!dwen(20090820)  SUBROUTINE METINP(BACK,KG,KT,KUNIT,KREC,LX1,LY1,NXS,NYS,NZS,MC,KEND,     &
!dwen(20090820)                    IFHR,ZT,DS,P0,T0,U0,V0,UF,VF,HF,RT,ZI,U,V,W,A,Q,P,E,H,X)

 subroutine metinp(BACK,KG,KT,KUNIT,KREC,LX1,LY1,NXS,NYS,NZS,MC,KEND,     &
                  IFHR,ZT,DS,P0,T0,U0,V0,UF,VF,HF,RT,ZI,U,V,W,A,Q,P,E,H,X,&
                  W0,fluxflg, deepflg, shallflg,              &
                  muu,muv,mu,msfu,msfv,msft,RC,LF,TC,LC,SW,SM,  &
                  TLRAMS,SIGWRAMS,cfu1,CFU2,CFD1,DFU1,       &
                  DFU2,EFU1,EFU2,DFD1,EFD1,TKE)

  IMPLICIT NONE
  LOGICAL,      INTENT(IN)    :: back      ! integration direction
  INTEGER,      INTENT(IN)    :: kg        ! number of active grid
  INTEGER,      INTENT(IN)    :: kt        ! number of active time
  INTEGER,      INTENT(IN)    :: kunit     ! input device unit number
  INTEGER,      INTENT(IN)    :: krec      ! record number of index record
  INTEGER,      INTENT(IN)    :: lx1,ly1   ! lower left  of subgrid FG unit
  INTEGER,      INTENT(IN)    :: nxs,nys   ! dimensions of sub-grid
  INTEGER,      INTENT(IN)    :: nzs       ! number of data levels to read
  INTEGER,      INTENT(INOUT) :: mc        ! accumulated minutes of data read
  INTEGER,      INTENT(IN)    :: kend      ! last valid record number of file
  INTEGER,      INTENT(OUT)   :: ifhr      ! current forecast hour
  REAL,         INTENT(OUT)   :: zt (:,:)  ! terrain height
  REAL,         INTENT(OUT)   :: ds (:,:)  ! downward shortwave flux
  REAL,         INTENT(OUT)   :: p0 (:,:)  ! surface pressure 
  REAL,         INTENT(OUT)   :: u0 (:,:)  ! low level u wind
  REAL,         INTENT(OUT)   :: v0 (:,:)  ! low level v wind
  REAL,         INTENT(OUT)   :: t0 (:,:)  ! low level temperature
  REAL,         INTENT(OUT)   :: uf (:,:)  ! u momentum flux
  REAL,         INTENT(OUT)   :: vf (:,:)  ! v momentum flux
  REAL,         INTENT(OUT)   :: hf (:,:)  ! sensible heat flux
  REAL,         INTENT(OUT)   :: rt (:,:)  ! rainfall total
  REAL,         INTENT(OUT)   :: zi (:,:)  ! mixed layer depth
  REAL,         INTENT(OUT)   :: u (:,:,:) ! upper level u wind
  REAL,         INTENT(OUT)   :: v (:,:,:) ! upper level v wind
  REAL,         INTENT(OUT)   :: w (:,:,:) ! upper level w wind
  REAL,         INTENT(OUT)   :: a (:,:,:) ! upper level temperature
  REAL,         INTENT(OUT)   :: q (:,:,:) ! upper level moisture
  REAL,         INTENT(OUT)   :: p (:,:,:) ! upper level pressure
  REAL,         INTENT(OUT)   :: e (:,:,:) ! turbulent kinetic energy or ...
                                           ! velocity variance v'2
  REAL,         INTENT(OUT)   :: h (:,:,:) ! velocity variance u'2
  REAL,         INTENT(OUT)   :: x (:,:,:) ! velocity variance w'2
!dwen(20090730)  ***************************
  real,         intent(out)   :: w0(:,:)
  real,         intent(out)   :: muu(:,:)
  real,         intent(out)   :: muv(:,:)
  real,         intent(out)   :: mu(:,:)
  Logical,      intent(inout) :: fluxflg,deepflg,shallflg
  real,         intent(out)   :: msfu(:,:)
  real,         intent(out)   :: msfv(:,:)
  real,         intent(out)   :: msft(:,:)
  real,         intent(out)   :: rc(:,:)
  real,         intent(out)   :: lf(:,:)
  real,         intent(out)   :: tc(:,:)
  real,         intent(out)   :: lc(:,:)
  real,         intent(out)   :: sw(:,:)
  real,         intent(out)   :: sm(:,:)
  real,         intent(out)   :: tlrams(:,:,:)
  real,         intent(out)   :: sigwrams(:,:,:)
  real,         intent(out)   :: cfu1(:,:,:)
  real,         intent(out)   :: cfu2(:,:,:)
  real,         intent(out)   :: cfd1(:,:,:)
  real,         intent(out)   :: dfu1(:,:,:)
  real,         intent(out)   :: dfu2(:,:,:)
  real,         intent(out)   :: dfd1(:,:,:)
  real,         intent(out)   :: efu1(:,:,:)
  real,         intent(out)   :: efu2(:,:,:)
  real,         intent(out)   :: efd1(:,:,:)
  real,         intent(out)   :: tke(:,:,:)
!dwen ***************************************
  END SUBROUTINE metinp
!-------------------------------------------------------------------------------
  SUBROUTINE ADVRKM(U,V,W,K1,K2,NLVL,DM,JET,ZMDL,XX,YY,ZZ,ZX,DT,BACK, &
                    GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,    INTENT(IN)    :: u(:,:,:,:)      ! u wind
  REAL,    INTENT(IN)    :: v(:,:,:,:)      ! v wind
  REAL,    INTENT(IN)    :: w(:,:,:,:)      ! w wind
  INTEGER, INTENT(IN)    :: k1,k2           ! last and next time index    
  INTEGER, INTENT(IN)    :: nlvl            ! number of levels to process
  REAL,    INTENT(IN)    :: dm              ! minutes between winds 1 and 2
  INTEGER, INTENT(IN)    :: jet             ! current elapsed time (minutes)
  REAL,    INTENT(IN)    :: zmdl            ! vertical model domain top (m)
  REAL,    INTENT(INOUT) :: xx,yy,zz        ! old position at time t
  REAL,    INTENT(IN)    :: dt              ! integration step (minutes)
  REAL,    INTENT(OUT)   :: zx              ! last estimate of vertical index
  LOGICAL, INTENT(IN)    :: back            ! flag to indicate direction
  LOGICAL, INTENT(IN)    :: global          ! global cyclic boundary conditions
  INTEGER, INTENT(IN)    :: nxp,nyp         ! grid dimensions for global  
  END SUBROUTINE advrkm
!-------------------------------------------------------------------------------
!dwen(2009082)  SUBROUTINE ADVIEC(U,V,W,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,DT,TRATIO, &
!dwen(2009082)                    BACK,GLOBAL,NXP,NYP)

 SUBROUTINE ADVIEC(U,V,W,d,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,   &
                  DT,TRATIO,BACK,GLOBAL,NXP,NYP,                    &
                  zter,zsg1,wwout,awrfflg,fluxflg,zsg,ramsflg,dead)

  IMPLICIT NONE
  REAL,    INTENT(IN)    :: u(:,:,:,:)      ! u wind
  REAL,    INTENT(IN)    :: v(:,:,:,:)      ! v wind
  REAL,    INTENT(IN)    :: w(:,:,:,:)      ! w wind
  INTEGER, INTENT(IN)    :: k1,k2           ! last and next time index    
  INTEGER, INTENT(IN)    :: nlvl            ! number of levels to process
  INTEGER, INTENT(IN)    :: mtime(2)        ! time of meteo observation 
  INTEGER, INTENT(IN)    :: jet             ! current elapsed time (minutes)
  REAL,    INTENT(IN)    :: zmdl            ! vertical model domain top (m)
  REAL,    INTENT(INOUT) :: xx,yy,zz        ! old position at time t
  REAL,    INTENT(IN)    :: dt              ! integration step (minutes)
  REAL,    INTENT(IN)    :: tratio          ! time step stability criterion
  REAL,    INTENT(OUT)   :: zx              ! last estimate of vertical index
  LOGICAL, INTENT(IN)    :: back            ! flag to indicate direction
  LOGICAL, INTENT(IN)    :: global          ! global cyclic boundary conditions
  INTEGER, INTENT(IN)    :: nxp,nyp         ! grid dimensions for global  
!dwen(20090810) *******************
  real,    intent(in)    :: d(:,:,:,:)      !air density
  real,    intent(in)    :: zter            !terrain height
  real,    intent(in)    :: zsg1            !1st sigma level
  real,    intent(out)   :: wwout           !wbar
  real,    intent(in)    :: zsg(:)
  logical, intent(in)    :: awrfflg,fluxflg,ramsflg
  logical, intent(out)   :: dead
!dwen ******************************
  END SUBROUTINE adviec
!-------------------------------------------------------------------------------
  SUBROUTINE ADVISO(U,V,P,ZSG,K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,ZX,DT,TRATIO, &
                    BACK,GLOBAL,NXP,NYP,UVEL,VVEL,RAMSFLG)
  IMPLICIT NONE
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
  END SUBROUTINE adviso
!-------------------------------------------------------------------------------
  SUBROUTINE ADVSFC(U,V,K1,K2,NLVL,MTIME,JET,XX,YY,DT,BACK,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,    INTENT(IN)    :: u(:,:,:)        ! u wind
  REAL,    INTENT(IN)    :: v(:,:,:)        ! v wind
  INTEGER, INTENT(IN)    :: k1,k2           ! last and next time index
  INTEGER, INTENT(IN)    :: nlvl            ! number of levels to process
  INTEGER, INTENT(IN)    :: mtime(2)        ! time of meteo observations
  INTEGER, INTENT(IN)    :: jet             ! current elapsed time (minutes)
  REAL,    INTENT(IN)    :: dt              ! integration step (minutes)
  LOGICAL, INTENT(IN)    :: back            ! flag to indicate direction
  LOGICAL, INTENT(IN)    :: global          ! global cyclic boundary conditions
  INTEGER, INTENT(IN)    :: nxp,nyp         ! global boundaries
  REAL,    INTENT(INOUT) :: xx,yy           ! old (t) and new (t+dt) position
  END SUBROUTINE advsfc
!-------------------------------------------------------------------------------
!dwen(20090820)  SUBROUTINE ADVMET(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,XP,YP,JET,MTIME,       &
!dwen(20090820)                    KCYCLE,NLVL,FHOUR,IFHR,K1,K2,GX,GY,Z0,LU,ZT,A,T,Q,P,E,X,  &
!dwen(20090820)                   ZI,H,U0,V0,RT,UF,VF,SF,SS,DS,DSWF,GLOBAL,NXP,NYP)


 SUBROUTINE ADVMET(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,XPi,YPi,JET,MTIME,      &
                   KCYCLE,NLVL,FHOUR,IFHR,K1,K2,GX,GY,Z0,LU,ZT,A,T,Q,P,E,X, &
                   ZI,H,U0,V0,RT,UF,VF,SF,SS,DS,DSWF,GLOBAL,NXP,NYP,        &
                   cflg,tclf,lclf,radf,slmf,                                &
                   u,v,d,xm,hm,dmass,tl,sigw,rc,hf,lf,tc,lc,sm,          &
                   ramsflg,ecmflg,awrfflg,fluxflg)

  IMPLICIT NONE
  INCLUDE 'DEFMETO.INC'
  TYPE(bset),INTENT(OUT)   :: metz (:)  ! profile advection variables
  TYPE(aset),INTENT(INOUT) :: meto      ! surface advection variables
  LOGICAL,   INTENT(IN)    :: back,vmix,cdep,rdep,traj  
  REAL,      INTENT(IN)    :: xpi,ypi     ! horizontal particle position
  INTEGER,   INTENT(IN)    :: jet       ! elapsed time
  INTEGER,   INTENT(IN)    :: mtime(2)  ! time of meteo observations 
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
  LOGICAL,   INTENT(IN)    :: dswf       ! downward shortwave flag 
  LOGICAL,   INTENT(IN)    :: global     ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp    ! grid dimensions for global  
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
  END SUBROUTINE advmet
!-------------------------------------------------------------------------------
  SUBROUTINE METGRD(KG,KT,LX1,LY1,NXS,NYS,GX,GY,Z0,LU,ZT)
  IMPLICIT NONE
  INTEGER,    INTENT(IN)     :: kg            ! grid selection index
  INTEGER,    INTENT(IN)     :: kt            ! time selection index
  INTEGER,    INTENT(IN)     :: lx1,ly1       ! subgrid lower left position
  INTEGER,    INTENT(IN)     :: nxs,nys       ! subgrid dimensions
  REAL,       INTENT(OUT)    :: gx (:,:)      ! grid size array (m)
  REAL,       INTENT(OUT)    :: gy (:,:)      ! grid size array (m)
  REAL,       INTENT(OUT)    :: z0 (:,:)      ! aerodynamic roughness length (m)
  INTEGER,    INTENT(OUT)    :: lu (:,:)      ! land-use category (1-11)
  REAL,       INTENT(INOUT)  :: zt (:,:)      ! terrain height (m)
  END SUBROUTINE metgrd
!-------------------------------------------------------------------------------
  SUBROUTINE METWND(K1,K2,KVEL,NXS,NYS,NLVL,DM,ZMDL,ZSG,ZT,U,V,W,P,T,A)
  IMPLICIT NONE
  INTEGER,     INTENT(IN)    :: k1,k2         ! time index pointers
  INTEGER,     INTENT(IN)    :: kvel          ! vertical motion flag
  INTEGER,     INTENT(IN)    :: nxs,nys       ! dimensions of sub-grid
  INTEGER,     INTENT(IN)    :: nlvl          ! number of data levels
  REAL,        INTENT(IN)    :: dm            ! time between obs
  REAL,        INTENT(IN)    :: zmdl          ! vertical model domain top (m)
  REAL,        INTENT(IN)    :: zsg (:)       ! sigma levels
  REAL,        INTENT(IN)    :: zt(:,:)       ! terrain height
  REAL,        INTENT(IN)    :: u (:,:,:)     ! velocity components
  REAL,        INTENT(IN)    :: v (:,:,:)     ! velocity components
  REAL,        INTENT(INOUT) :: w (:,:,:)     ! velocity components
  REAL,        INTENT(IN)    :: p (:,:,:,:)   ! pressure             
  REAL,        INTENT(IN)    :: t (:,:,:,:)   ! potential temperature          
  REAL,        INTENT(IN)    :: a (:,:,:,:)   ! ambient temperature  
  END SUBROUTINE metwnd
!-------------------------------------------------------------------------------
!dwen(20090820)  SUBROUTINE PRFCOM(TKERD,TKERN,KG,KT,KSFC,GX,GY,Z0,ZT,NXS,NYS,NZS,ZMDL,      &
!dwen(20090820)                    ZSG,NLVL,VMIX,KMIXD,KMIX0,ZI,P0,T0,U0,V0,UF,VF,HF,SF,     &
!dwen(20090820)                   SS,U,V,W,A,T,Q,P,E,H,X)

     subroutine PRFCOM(TKERD,TKERN,KG,KT,KSFC,GX,GY,Z0,              &
         ZT,NXS,NYS,NZS,ZMDL,ZSG,NLVL,VMIX,KMIXD,KMIX0,         &
         ZI,P0,T0,U0,V0,UF,VF,HF,SF,SS,U,V,W,A,T,Q,P,E,H,X,     &
         iconvect,w0,tl,d,sigw,zloc,delmass,           &
         ziscale,alt0,mu,muu,muv,msfu,msfv,msft,              &
         fluxflg, deepflg, shallflg, cfxup1,cfxup2,cfxdn1,      &
         dfxup1,dfxup2,efxup1,efxup2,dfxdn1,efxdn1,tke,xm,hm)     

  IMPLICIT NONE
  REAL,    INTENT(IN)    :: tkerd      ! day turbulent kinetic energy ratio    
  REAL,    INTENT(IN)    :: tkern      ! night turbulent kinetic energy ratio    
  INTEGER, INTENT(IN)    :: kg         ! grid indicator
  INTEGER, INTENT(IN)    :: kt         ! time indicator
  INTEGER, INTENT(IN)    :: ksfc       ! index top of the sfc layer
  REAL,    INTENT(IN)    :: gx (:,:)   ! grid size (m)
  REAL,    INTENT(IN)    :: gy (:,:)   ! grid size (m)
  REAL,    INTENT(IN)    :: z0 (:,:)   ! aerodynamic roughness length (m)
  REAL,    INTENT(INOUT) :: zt (:,:)   ! terrain height elevations (m)
  INTEGER, INTENT(IN)    :: nxs,nys    ! subgrid dimensions
  INTEGER, INTENT(IN)    :: nzs        ! number of input data levels
  REAL,    INTENT(IN)    :: zmdl       ! internal model top (meters)
  REAL,    INTENT(IN)    :: zsg (:)    ! array internal sigma levels
  INTEGER, INTENT(IN)    :: nlvl       ! number of output levels
  LOGICAL, INTENT(IN)    :: vmix       ! indicator for mixing computation
  INTEGER, INTENT(IN)    :: kmixd      ! mixed layer depth options
  INTEGER, INTENT(IN)    :: kmix0      ! minimum mixing depth
  REAL,    INTENT(INOUT) :: zi (:,:)   ! mixing depth  
  REAL,    INTENT(INOUT) :: p0 (:,:)   ! pressue surface variable 
  REAL,    INTENT(IN)    :: t0 (:,:)   ! temperature 
  REAL,    INTENT(IN)    :: u0 (:,:)   ! u velocity 
  REAL,    INTENT(IN)    :: v0 (:,:)   ! v velocity 
  REAL,    INTENT(INOUT) :: uf (:,:)   ! u momentum flux and u*
  REAL,    INTENT(INOUT) :: vf (:,:)   ! v momentum flux and v* 
  REAL,    INTENT(IN)    :: hf (:,:)   ! sensible heat flux 
  REAL,    INTENT(OUT)   :: sf (:,:)   ! stability function 
  REAL,    INTENT(OUT)   :: ss (:,:)   ! static stability   
  REAL,    INTENT(INOUT) :: u  (:,:,:) ! wind 3d variables       
  REAL,    INTENT(INOUT) :: v  (:,:,:) ! v component       
  REAL,    INTENT(INOUT) :: w  (:,:,:) ! w component        
  REAL,    INTENT(INOUT) :: a  (:,:,:) ! ambient temperature        
  REAL,    INTENT(OUT)   :: t  (:,:,:) ! potential temperature        
  REAL,    INTENT(INOUT) :: q  (:,:,:) ! moisture           
  REAL,    INTENT(INOUT) :: p  (:,:,:) ! pressure           
  REAL,    INTENT(INOUT) :: e (:,:,:)  ! turbulent kinetic energy or ...
                                       ! velocity variance v'2
  REAL,    INTENT(INOUT) :: h (:,:,:)  ! velocity variance u'2
  REAL,    INTENT(INOUT) :: x (:,:,:)  ! velocity variance w'2
!dwen(20090803) *******************
  real,      intent(in)       :: w0(:,:)
  real,      intent(inout)    :: d(:,:,:)
  real,      intent(out)      :: tl(:,:,:)
  real,      intent(out)      :: sigw(:,:,:)
  real,      intent(out)      :: hm(:,:,:) !horizontal mixing coefficient
! CHG:(12/04/2001)array to store lim. of conv. heights
  real,      intent(out)         :: zloc(:,:)
! JCL:(4/3/02)mass violation grid [fraction of gridcell/min]
  real,      intent(out)   :: delmass(:,:,:)
  real,      intent(in)         :: ziscale
  real,      intent(in)         :: alt0(:,:,:)
  real,      intent(in)         :: mu(:,:)
  real,      intent(in)         :: muu(:,:)
  real,      intent(in)         :: muv(:,:)
  real,      intent(in)         :: msfu(:,:)
  real,      intent(in)         :: msfv(:,:)
  real,      intent(in)         :: msft(:,:)
  logical,   intent(inout) :: fluxflg  ! flag for WRF momentum flux input
  logical,   intent(inout) :: deepflg, shallflg
  real,      intent(inout)         :: cfxup1(:,:,:)
  real,      intent(inout)         :: cfxup2(:,:,:)
  real,      intent(inout)         :: cfxdn1(:,:,:)
  real,      intent(inout)         :: dfxup1(:,:,:)
  real,      intent(inout)         :: dfxup2(:,:,:)
  real,      intent(inout)         :: efxup1(:,:,:)
  real,      intent(inout)         :: efxup2(:,:,:)
  real,      intent(inout)         :: dfxdn1(:,:,:)
  real,      intent(inout)         :: efxdn1(:,:,:)
  real,      intent(inout)         :: tke(:,:,:)  !turbulence kinetic energy
  real,      intent(out)         :: xm(:,:,:)
  integer,   intent(in)            :: iconvect
!************************************
  END SUBROUTINE prfcom
!-------------------------------------------------------------------------------
  SUBROUTINE METSUB(XP,YP,KG,KT,OFFG,LX1,LY1,LXR,LYR,KRET)
  use module_defgrid         ! meteorology grid and file
  IMPLICIT NONE
!  INCLUDE 'DEFGRID.INC' ! meteorology grid and file
  REAL,      INTENT(IN)    :: xp,yp          ! x,y grid particle position
  INTEGER,   INTENT(IN)    :: kg             ! grid number being loaded
  INTEGER,   INTENT(IN)    :: kt             ! grid number being loaded
  LOGICAL,   INTENT(OUT)   :: offg           ! load new subgrid flag  
  INTEGER,   INTENT(INOUT) :: lx1,ly1        ! subgrid corner point
  INTEGER,   INTENT(INOUT) :: lxr,lyr        ! subgrid range        
  INTEGER,   INTENT(OUT)   :: kret           ! particle termination
  END SUBROUTINE metsub
!-------------------------------------------------------------------------------
!dwen(20090315): add gasdev
      FUNCTION GASDEV(IDUM,SIGMA)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER          :: IDUM
      REAL             :: GASDEV
      REAL             :: SIGMA
      INTEGER          :: ISET
      REAL             :: FAC,GSET,RSQ,V1,V2
      REAL             :: RAN3
      END FUNCTION gasdev
!dwen(20090315):add ran3
      FUNCTION RAN3(IDUMI)
!
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER, INTENT(IN) :: IDUMI
      real                :: ran3
end function ran3

END INTERFACE
