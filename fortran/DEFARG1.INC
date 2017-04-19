!-------------------------------------------------------------------------------
! Configure subroutine interface argument lists for main programs         
! Last Revised: 17 Nov 2000 (RRD) - conini, condsk
!               14 May 2001 (RRD) - mpi implementation of emspnt
!               21 Jun 2001 (RRD) - temperature array added
!               29 Aug 2001 (RRD) - simultaneous multiple meteorology
!               24 Oct 2001 (RRD) - simplified backward dispersion options
!               21 Mar 2002 (RRD) - pm10 dust emission routine
!               09 Sep 2002 (RRD) - fortran coding standards
!               13 Feb 2002 (RRD) - added subroutine metlvl
!               08 May 2003 (RRD) - added parinp and parout 
!               07 Aug 2003 (RRD) - added metsum and mod to condsk
!               15 Sep 2003 (RRD) - deposition probability to depelm
!               10 Aug 2004 (RRD) - enhanced particle-puff conversions
!               13 Oct 2004 (RRD) - horizontal puff dispersion options
!               22 Dec 2004 (RRD) - non regular concentration grids
!               25 May 2005 (RRD) - added massum 
!               13 Oct 2005 (RRD) - added lagrangian sampling emslag
!               25 May 2006 (AS)  - random generator seed in pardsp
!               21 Nov 2006 (RRD) - day night tke partition
!               19 Mar 2008 (AS)  - emission duration to emrise
!               04 Jun 2008 (RRD) - mixing options to advpnt 
!               01 Jul 2008 (RRD) - changed conversion test to age
!               15 Aug 2008 (RRD) - more mixing options metini
!-------------------------------------------------------------------------------

INTERFACE

!-------------------------------------------------------------------------------
!dwen(20090823)  SUBROUTINE PARDSP(umix,vmix,gdisx,gdisy,dtm,zmdl,zsfc,nlvl,wmix,zsg,xpos,  &
!dwen(20090823)                    ypos,zpos,vprime,wprime,uprime,hdwp,zndx,iseed)

SUBROUTINE pardsp(umix,vmix,gdisx,gdisy,dtm,zmdl,zsfc,nlvl,wmix,zsg,xpos,&
                  ypos,zpos,vprime,wprime,uprime,hdwp,zndx,iseed,        &
                  tlprev,stdwprev,zml1,densprev,tlnext,stdwnext,zml2,    &
!dwen(20090906)                  densnext,tlfrac,xposprev,yposprev,uuprev,uunext,       &
                  densnext,tlfrac,uuprev,uunext,       &
                  vvprev,vvnext,dmassprev,dmassnext,dmweight,rseed,      &
                  sigmaw,sigmau,samptt,samptt2,wwprev,veght,             &
                 convdur,zloc,ramsflg,tempnext,tempprev,rhfrnext,        &
                  rhfrprev,sphu,foot)

  IMPLICIT NONE
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
  REAL,         INTENT(INOUT)  :: wprime      ! last w turbulent velocity (z)
  REAL,         INTENT(INOUT)  :: uprime      ! last u turbulent velocity (x)
  INTEGER,      INTENT(IN)     :: hdwp        ! horizontal distribution index
  REAL,         INTENT(IN)     :: zndx        ! fractional vertical index  
  INTEGER,      INTENT(INOUT)  :: iseed       ! random generator seed

!dwen(20090823)****************************************
      REAL, intent(inout)         :: TLPREV(:), STDWPREV(:)
      REAL, intent(inout)         :: TLNEXT(:), STDWNEXT(:)
      REAL,intent(in)          :: DENSPREV(:), DENSNEXT(:)
      REAL,intent(in)          :: UUPREV(:), VVPREV(:)
      REAL,intent(in)          :: UUNEXT(:), VVNEXT(:)
      real,intent(inout)    :: CONVDUR
      REAL,intent(inout)          :: DMWEIGHT
      REAL,intent(in)          :: TEMPPREV(:), TEMPNEXT(:)
      REAL,intent(in)          :: RHFRPREV(:), RHFRNEXT(:)
!      REAL,intent(in)          :: xposPREV
!      REAL,intent(in)          :: yposPREV
      real,intent(in)          :: dmassprev(:),dmassnext(:)
      REAL,intent(in)           :: TLFRAC
      REAL,intent(out)           :: sphu
      INTEGER,intent(in)        :: RSEED
      REAL,intent(inout)          :: SIGMAW, SIGMAU
      REAL,intent(in)           :: ZML1,ZML2
      REAL,intent(out)          :: SAMPTT, FOOT
      REAL,intent(out)          :: SAMPTT2
      REAL,intent(inout)   :: WWPREV
      REAL,intent(in)   :: VEGHT
      REAL,intent(in)   :: ZLOC
      LOGICAL, INTENT(IN) :: RAMSFLG
!********************************************

  END SUBROUTINE pardsp
!-------------------------------------------------------------------------------
  SUBROUTINE CONSUM(conc,numgrd,plat,plon,dt,jet,zmdl,zsfc,kmsl,     &
                  cgsize,mass,dept,zpos,sigh,sigv,hdwp,ptyp,csum)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(cset), INTENT(IN) :: conc(:)           ! for each concentration grid 
  INTEGER, INTENT(IN)    :: numgrd            ! number of concentration grids
  REAL,    INTENT(IN)    :: plat,plon         ! particle position lat & long
  REAL,    INTENT(IN)    :: dt                ! time step (min)
  INTEGER, INTENT(IN)    :: jet               ! elapsed time (min)
  REAL,    INTENT(IN)    :: zmdl              ! model domain top
  REAL,    INTENT(IN)    :: zsfc              ! height of ground surface (m)
  INTEGER, INTENT(IN)    :: kmsl              ! agl=0 or msl=1 height units
  REAL,    INTENT(INOUT) :: cgsize            ! minimum sampling grid size (km)
  REAL,    INTENT(IN)    :: mass (:)          ! mass of pollutant 
  REAL,    INTENT(IN)    :: dept (:)          ! deposition amount of mass
  REAL,    INTENT(IN)    :: zpos              ! puff center height (sigma)
  REAL,    INTENT(IN)    :: sigh,sigv         ! horiz (m) and vert sigma (sigma)
  INTEGER, INTENT(IN)    :: hdwp              ! horiz distribution in pollutant
  INTEGER, INTENT(IN)    :: ptyp              ! pollutant type index number
  REAL,    INTENT(INOUT) :: csum (:,:,:,:,:)  ! summation matrix             
  END SUBROUTINE consum
!-------------------------------------------------------------------------------
  SUBROUTINE MASSUM(conc,numgrd,plat,plon,dt,jet,zmdl,zsfc,kmsl,   &
                    cgsize,mass,zpos,ptyp,csum)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC' 
  TYPE(cset), INTENT(IN) :: conc(:)           ! for each concentration grid
  INTEGER, INTENT(IN)    :: numgrd            ! number of concentration grids
  REAL,    INTENT(IN)    :: plat,plon         ! particle position lat & long
  REAL,    INTENT(IN)    :: dt                ! time step (min)
  INTEGER, INTENT(IN)    :: jet               ! elapsed time (min)
  REAL,    INTENT(IN)    :: zmdl              ! model domain top
  REAL,    INTENT(IN)    :: zsfc              ! height of ground surface (m)
  INTEGER, INTENT(IN)    :: kmsl              ! agl=0 or msl=1 height units
  REAL,    INTENT(INOUT) :: cgsize            ! minimum sampling grid size (km)
  REAL,    INTENT(IN)    :: mass (:)          ! mass of pollutant
  REAL,    INTENT(IN)    :: zpos              ! puff center height (sigma)
  INTEGER, INTENT(IN)    :: ptyp              ! pollutant type index number
  REAL,    INTENT(INOUT) :: csum (:,:,:,:,:)  ! concentration summation matrix
  END SUBROUTINE massum
!-------------------------------------------------------------------------------
  SUBROUTINE METSUM(conc,numgrd,xp,yp,gdx,gdy,dt,jet,zmdl,zsfc,kmsl,        &
                    cgsize,mass,dept,zpos,sigh,sigv,hdwp,ptyp,csum)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC' ! pollutant and concentration grid
  TYPE(cset), INTENT(IN) :: conc(:)           ! for each concentration grid
  INTEGER, INTENT(IN)    :: numgrd            ! number of concentration grids
  REAL,    INTENT(IN)    :: xp,yp             ! particle position lat & long
  REAL,    INTENT(IN)    :: gdx,gdy           ! grid distance at part position
  REAL,    INTENT(IN)    :: dt                ! time step (min)
  INTEGER, INTENT(IN)    :: jet               ! elapsed time (min)
  REAL,    INTENT(IN)    :: zmdl              ! model domain top
  REAL,    INTENT(IN)    :: zsfc              ! height of ground surface (m)
  INTEGER, INTENT(IN)    :: kmsl              ! agl=0 or msl=1 height units
  REAL,    INTENT(INOUT) :: cgsize            ! minimum sampling grid size (km)
  REAL,    INTENT(IN)    :: mass (:)          ! mass of pollutant
  REAL,    INTENT(IN)    :: dept (:)          ! deposition amount of mass
  REAL,    INTENT(IN)    :: zpos              ! puff center height (sigma)
  REAL,    INTENT(IN)    :: sigh,sigv         ! horiz (m) and vert sigma (sigma)
  INTEGER, INTENT(IN)    :: hdwp              ! horiz distribution in pollutant
  INTEGER, INTENT(IN)    :: ptyp              ! pollutant type index number
  REAL,    INTENT(INOUT) :: csum (:,:,:,:,:)  ! concentration summation matrix
  END SUBROUTINE metsum
!-------------------------------------------------------------------------------
  SUBROUTINE CONINI(spot,conc,dirt,nloc,numgrd,numtyp,cpack)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  INCLUDE 'DEFSPOT.INC'
  TYPE(rset), INTENT(IN)    :: spot(:)   ! source location characteristics
  TYPE(cset), INTENT(INOUT) :: conc(:)   ! for each concentration grid 
  TYPE(pset), INTENT(IN)    :: dirt(:)   ! for each pollutant type 
  INTEGER,    INTENT(IN)    :: nloc      ! number of sources
  INTEGER,    INTENT(IN)    :: numgrd    ! number of concentration grids
  INTEGER,    INTENT(IN)    :: numtyp    ! number of pollutants
  INTEGER,    INTENT(IN)    :: cpack     ! concentration packing flag
  END SUBROUTINE conini
!-------------------------------------------------------------------------------
  SUBROUTINE CONDSK(conc,dirt,ichem,kgm,ktm,numgrd,numtyp,dt,jet,ifhr,     &
                    cpack,csum)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(cset),  INTENT(IN)     :: conc(:)  ! for each concentration grid 
  TYPE(pset),  INTENT(IN)     :: dirt(:)  ! for each pollutant type 
  INTEGER,     INTENT(IN)     :: ichem               ! chemistry parameter
  INTEGER,     INTENT(IN)     :: kgm                 ! current meteo grid number
  INTEGER,     INTENT(IN)     :: ktm                 ! current meteo time number
  INTEGER,     INTENT(IN)     :: numgrd              ! number of concen grids
  INTEGER,     INTENT(IN)     :: numtyp              ! number of pollutants
  REAL,        INTENT(IN)     :: dt                  ! time step
  INTEGER,     INTENT(IN)     :: jet                 ! current elapsed time
  INTEGER,     INTENT(IN)     :: ifhr                ! current forecast hour
  INTEGER,     INTENT(IN)     :: cpack               ! concen packing flag
  REAL,        INTENT(IN)     :: csum (:,:,:,:,:)    ! summation matrix
  END SUBROUTINE condsk
!-------------------------------------------------------------------------------
  SUBROUTINE CONSET(CONC,ZMDL,NUMGRD,IUNIT,CGSIZE,CPACK,OLAT,OLON,IBYR,IBMO, &
                    IBDA,IBHR,BACK)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(cset), INTENT(OUT) :: conc(:)              ! for each concentration grid 
  REAL,    INTENT(IN)     :: zmdl                 ! model domain top meters AGL
  INTEGER, INTENT(IN)     :: numgrd               ! number concentration grids
  INTEGER, INTENT(IN)     :: iunit                ! unit number for input
  REAL,    INTENT(INOUT)  :: cgsize               ! minimum grid spacing (km) 
  INTEGER, INTENT(IN)     :: cpack                ! concentration grid packing
  REAL,    INTENT(IN)     :: olat, olon           ! starting point
  INTEGER, INTENT(IN)     :: ibyr,ibmo,ibda,ibhr  ! starting time
  LOGICAL, INTENT(IN)     :: back                 ! integration direction
  END SUBROUTINE conset 
!-------------------------------------------------------------------------------
  SUBROUTINE CONZRO(conc,numgrd,dt,jet,ifhr,csum)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(cset), INTENT(INOUT)   :: conc(:)             ! each concentration grid 
  INTEGER,     INTENT(IN)     :: numgrd              ! number of concen grids
  REAL,        INTENT(IN)     :: dt                  ! time step (min)
  INTEGER,     INTENT(IN)     :: jet                 ! current elapsed time
  INTEGER,     INTENT(IN)     :: ifhr                ! current forecast hour
  REAL,        INTENT(INOUT)  :: csum (:,:,:,:,:)    ! summation matrix
  END SUBROUTINE conzro
!-------------------------------------------------------------------------------
  SUBROUTINE DEPSUS(CONC,DIRT,INITD,KGM,KTM,NUMGRD,NUMTYP,DT,ICHEM,KPM,  &
             CSUM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,     &
             PTYP,PGRD,NSORT,MAXPAR)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(cset), INTENT(IN)    :: conc(:)      ! for each concentration grid 
  TYPE(pset), INTENT(IN)    :: dirt(:)      ! for each pollutant type 
  INTEGER,    INTENT(IN)    :: initd        ! initial puff/particle distribution
  INTEGER,    INTENT(IN)    :: kgm          ! current meteorological grid number
  INTEGER,    INTENT(IN)    :: ktm          ! current meteorological time number
  INTEGER,    INTENT(IN)    :: numgrd       ! number of concentration grids
  INTEGER,    INTENT(IN)    :: numtyp       ! number of pollutants
  REAL,       INTENT(IN)    :: dt           ! time step (min)
  INTEGER,    INTENT(IN)    :: ichem        ! chemistry parameter
  INTEGER,    INTENT(INOUT) :: kpm          ! total number of puffs or particles
  REAL,       INTENT(INOUT) :: csum  (:,:,:,:,:)                
  REAL,       INTENT(INOUT) :: mass  (:,:)  ! (species, particles)
  REAL,       INTENT(INOUT) :: xpos  (:)    ! x position (numb part)
  REAL,       INTENT(INOUT) :: ypos  (:)    ! y position (numb part)
  REAL,       INTENT(INOUT) :: zpos  (:)    ! z position (numb part)
  REAL,       INTENT(INOUT) :: sigh  (:)    ! horizontal sigma    
  REAL,       INTENT(INOUT) :: sigu  (:)    ! turbulence
  REAL,       INTENT(INOUT) :: sigv  (:)    ! turbulence
  REAL,       INTENT(INOUT) :: sigw  (:)    ! vertical sigma
  INTEGER,    INTENT(INOUT) :: hdwp  (:)    ! Horizontal distribution
  INTEGER,    INTENT(INOUT) :: page  (:)    ! pollutant age (min)     
  INTEGER,    INTENT(INOUT) :: ptyp  (:)    ! pollutant type index numb 
  INTEGER,    INTENT(INOUT) :: pgrd  (:)    ! particle on this grid numb
  INTEGER,    INTENT(INOUT) :: nsort (:)    ! sorted array index values
  INTEGER,    INTENT(IN)    :: maxpar       ! maximum particle count
  END SUBROUTINE depsus
!-------------------------------------------------------------------------------
  SUBROUTINE DEPRAD(CONC,DIRT,NUMGRD,NUMTYP,DT,CSUM)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(cset),   INTENT(IN)    :: conc(:)   ! for each concentration grid 
  TYPE(pset),   INTENT(IN)    :: dirt(:)   ! for each pollutant type 
  INTEGER,      INTENT(IN)    :: numgrd    ! number of concentration grids
  INTEGER,      INTENT(IN)    :: numtyp    ! number of pollutants
  REAL,         INTENT(IN)    :: dt        ! intgration time step (min)
  REAL,         INTENT(INOUT) :: csum (:,:,:,:,:) 
  END SUBROUTINE deprad
!-------------------------------------------------------------------------------
!dwen(20090823)  SUBROUTINE ADVPNT(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,KSFC,TKERD,TKERN,HDWP,  &
!dwen(20090823)                    XP,YP,ZP,JET,DT,TRATIO,KG,KT,NGRD,NTIM,ZSG,NLVL,ZMDL,KVEL, &
!dwen(20090823)                    KMIXD,KMIX0,UBAR,IFHR,KRET)

 SUBROUTINE ADVPNT(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,KSFC,TKERD,TKERN,HDWP,   &
                   XP,YP,ZP,JET,DT,TRATIO,KG,KT,NGRD,NTIM,ZSG,NLVL,ZMDL,KVEL,  &
                   KMIXD,KMIX0,UBAR,IFHR,KRET,                                 &
                   conc,wwout,zicontroltf,nhrszi,zipresc,convdur1,iconvect,       &
                   ramsflg,ecmflg,xposprev,yposprev)

  IMPLICIT NONE
!dwen(20090823)
  include 'DEFCONC.INC'
  INCLUDE 'DEFMETO.INC'
  TYPE(bset),INTENT(OUT)   :: metz (:)     ! profile advection variables
  TYPE(aset),INTENT(INOUT) :: meto         ! surface advection variables
  LOGICAL,   INTENT(IN)    :: back         ! integration direction flag
  LOGICAL,   INTENT(IN)    :: vmix         ! return mixing profile flag
  LOGICAL,   INTENT(IN)    :: cdep         ! return deposition variable flag
  LOGICAL,   INTENT(IN)    :: rdep         ! resistance deposition flag
  LOGICAL,   INTENT(IN)    :: traj         ! return trajectory variables flag
  INTEGER,   INTENT(IN)    :: ksfc         ! top of surface layer index
  REAL,      INTENT(IN)    :: tkerd        ! day turbulent kinetic energy ratio
  REAL,      INTENT(IN)    :: tkern        ! night turbulent kinetic energy ratio
  INTEGER,   INTENT(IN)    :: hdwp         ! puff/particle distribution type 
  REAL,      INTENT(INOUT) :: xp,yp,zp     ! particle position for advection
  INTEGER,   INTENT(IN)    :: jet          ! elapsed time (min)
  REAL,      INTENT(IN)    :: dt           ! advection time step (min)
  REAL,      INTENT(IN)    :: tratio       ! time step stability criterion
  INTEGER,   INTENT(INOUT) :: kg           ! grid index for calculation
  INTEGER,   INTENT(INOUT) :: kt           ! index to indicate current time
  INTEGER,   INTENT(IN)    :: ngrd         ! number of meteo grids
  INTEGER,   INTENT(IN)    :: ntim         ! number of meteo times
  REAL,      INTENT(IN)    :: zsg (:)      ! vertical sigma levels
  INTEGER,   INTENT(IN)    :: nlvl         ! number of vertical levels
  REAL,      INTENT(IN)    :: zmdl         ! model top m agl
  INTEGER,   INTENT(IN)    :: kvel         ! vertical velocity remapping option
  INTEGER,   INTENT(IN)    :: kmixd        ! mixed layer depth options
  INTEGER,   INTENT(IN)    :: kmix0        ! minimum mixing depth
  REAL,      INTENT(OUT)   :: ubar         ! advection velocity for component
  INTEGER,   INTENT(OUT)   :: ifhr         ! current forecast hour
  INTEGER,   INTENT(OUT)   :: kret         ! return for point off grid
!dwen(20090813) ******************
  type(cset),intent(in)    :: conc(:)
  real,      intent(out)   :: wwout 
  integer,   intent(in)    :: zicontroltf
  real,      intent(in)    :: zipresc(:) 
  integer,   intent(in)    :: nhrszi
  real,      intent(out)   :: convdur1
  integer,   intent(in)    :: iconvect
  logical,   intent(in)    :: ramsflg,ecmflg
  real,      intent(out)   :: xposprev,yposprev
!dwen *****************************
  END SUBROUTINE advpnt
!-------------------------------------------------------------------------------
  SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
  INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time
  END SUBROUTINE tm2day
!-------------------------------------------------------------------------------
  SUBROUTINE DEPSET(DIRT,NUMTYP,IUNIT,CDEP,RDEP,SDEP)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(pset), INTENT(INOUT) :: dirt(:)   ! for each pollutant type 
  INTEGER,    INTENT(IN)    :: numtyp    ! number of pollutant types
  INTEGER,    INTENT(IN)    :: iunit     ! unit number for input data
  LOGICAL,    INTENT(OUT)   :: cdep      ! indicate wet or dry deposition
  LOGICAL,    INTENT(OUT)   :: rdep      ! resistance deposition
  LOGICAL,    INTENT(OUT)   :: sdep      ! resuspension option
  END SUBROUTINE depset
!-------------------------------------------------------------------------------
  SUBROUTINE DEPELM(DIRT,OLAT,IBMO,NLVL,ZX,DT,ZMDL,ZSFC,ZSG,MASS,DEPT,ZPOS, &
             SIGV,ICHEM,KTYP,LAND,ROUG,SFCL,USTR,PSI,SFLX,HDWP,RAIN,        &
             DD,TT,QQ,DRYV)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(pset), INTENT(IN)   :: dirt (:) ! for each pollutant type 
  REAL,     INTENT(IN)     :: olat     ! origin location
  INTEGER,  INTENT(IN)     :: ibmo     ! computational month
  INTEGER,  INTENT(IN)     :: nlvl     ! number of vertical levels
  REAL,     INTENT(IN)     :: zx       ! vertical array position (index)
  REAL,     INTENT(IN)     :: dt       ! time step (min)
  REAL,     INTENT(IN)     :: zmdl     ! model domain top (meters)
  REAL,     INTENT(IN)     :: zsfc     ! height of terrain surface (m)
  REAL,     INTENT(IN)     :: zsg  (:) ! model sigma levels
  REAL,     INTENT(INOUT)  :: mass (:) ! mass (arbitrary units)
  REAL,     INTENT(INOUT)  :: zpos     ! puff center height (sigma)
  REAL,     INTENT(INOUT)  :: sigv     ! vert sigma (sigma)
  INTEGER,  INTENT(IN)     :: ichem    ! special depositon options   
  INTEGER,  INTENT(IN)     :: ktyp     ! pollutant type index number
  INTEGER,  INTENT(IN)     :: land     ! land use category (1-11)
  REAL,     INTENT(IN)     :: roug     ! aerodynamic roughness length (m)
  REAL,     INTENT(IN)     :: sfcl     ! height of the surface layer (m)
  REAL,     INTENT(IN)     :: ustr     ! friction velocity (m/s)
  REAL,     INTENT(IN)     :: psi      ! integrated stability function for heat
  REAL,     INTENT(IN)     :: sflx     ! incident short wave flux (w/m2)
  INTEGER,  INTENT(INOUT)  :: hdwp     ! pollutant distribution type (index)
  REAL,     INTENT(IN)     :: rain     ! precipitation value (m/min)
  REAL,     INTENT(IN)     :: dd   (:) ! air density profile (kg/m3)
  REAL,     INTENT(IN)     :: tt   (:) ! temperature profile
  REAL,     INTENT(IN)     :: qq   (:) ! humidity profile (fraction 0-1)
  REAL,     INTENT(OUT)    :: dept (:) ! deposition total (mass units)
  REAL,     INTENT(OUT)    :: dryv (:) ! dry deposition velocity as computed
  END SUBROUTINE depelm
!-------------------------------------------------------------------------------
  SUBROUTINE DEPDRY(KPOL,LAND,ROUG,SFCL,USTR,PSI,SFLX,AIRD,TEMP,PDIA,VG,VD)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)    :: kpol       ! polluant index number
  INTEGER,  INTENT(IN)    :: land       ! land-use category
  REAL,     INTENT(IN)    :: roug       ! aerodynamic roughness length (m)
  REAL,     INTENT(IN)    :: sfcl       ! height of constant flux layer (m)
  REAL,     INTENT(IN)    :: ustr       ! friction velocity (m/s)
  REAL,     INTENT(IN)    :: psi        ! integrated stability function heat
  REAL,     INTENT(IN)    :: sflx       ! solar irradiation at sfc (watts/m2)
  REAL,     INTENT(IN)    :: aird       ! ambient air density (g/m3)
  REAL,     INTENT(IN)    :: temp       ! canopy air temperature (deg K)
  REAL,     INTENT(IN)    :: pdia       ! particle diameter (meters)
  REAL,     INTENT(IN)    :: vg         ! gravitational settling velocity (m/s)
  REAL,     INTENT(OUT)   :: vd         ! total deposition velocity (m/s)
  END SUBROUTINE depdry
!-------------------------------------------------------------------------------
  SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
  INTEGER,  INTENT(OUT)  :: macc                 ! minutes since 1 Jan 1970
  END SUBROUTINE tm2min
!-------------------------------------------------------------------------------
  SUBROUTINE METINI(KDEF,KZMIX,TVMIX,KBLS,KBLT,NGRD,NTIM,OLAT,IBYR,     &
                    IBMO,IBDA,IBHR,IBMN,BACK,KVEL,KFOR)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: kdef                 ! deformation horizontal mix
  INTEGER,   INTENT(IN)    :: kzmix                ! vertical mixing averaging  
  REAL,      INTENT(IN)    :: tvmix                ! trop vert mix scale factor
  INTEGER,   INTENT(IN)    :: kbls                 ! pbl stability derived from
  INTEGER,   INTENT(IN)    :: kblt                 ! pbl scheme for mixing
  INTEGER,   INTENT(INOUT) :: ngrd                 ! number of meteo grids 
  INTEGER,   INTENT(INOUT) :: ntim                 ! number of meteo times 
  REAL,      INTENT(IN)    :: olat                 ! source latitude
  INTEGER,   INTENT(INOUT) :: ibyr,ibmo,ibda       ! calculation starting date
  INTEGER,   INTENT(INOUT) :: ibhr,ibmn            ! calculation starting time
  LOGICAL,   INTENT(IN)    :: back                 ! integration direction 
  INTEGER,   INTENT(INOUT) :: kvel                 ! vertical motion method
  INTEGER,   INTENT(IN)    :: kfor                 ! 0=NoForecasts 1=AllData
  END SUBROUTINE metini
!-------------------------------------------------------------------------------
  SUBROUTINE METLVL(NGRD,ZDATA,NLVL,KSFC,SFCL)
  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: NGRD  ! number of meteo grids       
  REAL,    INTENT(IN)   :: ZDATA ! maximum height for input met data
  INTEGER, INTENT(OUT)  :: NLVL  ! internal grid number levels
  INTEGER, INTENT(OUT)  :: KSFC  ! index of the surface layer
  REAL,    INTENT(OUT)  :: SFCL  ! height of the surface layer
  END SUBROUTINE metlvl
!-------------------------------------------------------------------------------
  SUBROUTINE EMSSET(DIRT,NUMTYP,IUNIT,IBYR,IBMO,IBDA,IBHR,OLAT,OLON,BACK)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  TYPE(pset),   INTENT(OUT)   :: dirt(:)             ! for each pollutant type 
  INTEGER,      INTENT(IN)    :: numtyp              ! numb of pollutant types
  INTEGER,      INTENT(IN)    :: iunit               ! unit numb for input data
  INTEGER,      INTENT(IN)    :: ibyr,ibmo,ibda,ibhr ! starting time
  REAL,         INTENT(INOUT) :: olat,olon           ! emission location
  LOGICAL,      INTENT(IN)    :: back                ! integration direction
  END SUBROUTINE emsset
!-------------------------------------------------------------------------------
  SUBROUTINE EMSPNT(SPOT,DIRT,NLOC,NUMTYP,KPM,INITD,DT,JET,NSORT,MASS,        &
             XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,QCYCLE,   &
             NUMPAR,MAXPAR,job_id,num_job,ichem,kret)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  INCLUDE 'DEFSPOT.INC'
  TYPE(rset),INTENT(INOUT) :: spot(:) ! source location characteristics
  TYPE(pset),INTENT(INOUT) :: dirt(:) ! for each pollutant type 
  INTEGER,   INTENT(IN)    :: nloc    ! total number of source locations
  INTEGER,   INTENT(IN)    :: numtyp  ! number of pollutant types
  INTEGER,   INTENT(INOUT) :: kpm     ! number of puffs or particles
  INTEGER,   INTENT(IN)    :: initd   ! initial distribution type
  REAL,      INTENT(IN)    :: dt      ! time step (min)
  INTEGER,   INTENT(IN)    :: jet     ! current elapsed time (min)
  REAL,      INTENT(IN)    :: qcycle  ! optional emission cycle time in hours
  INTEGER,   INTENT(IN)    :: numpar  ! maximum number of particles permitted
  INTEGER,   INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,      INTENT(INOUT) :: mass  (:,:) ! mass of pollutant (arbitrary units)
  REAL,      INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,      INTENT(INOUT) :: sigh  (:)   ! horizontal sigma    
  REAL,      INTENT(INOUT) :: sigu  (:)   ! turbulence
  REAL,      INTENT(INOUT) :: sigv  (:)   ! turbulence
  REAL,      INTENT(INOUT) :: sigw  (:)   ! vertical sigma
  INTEGER,   INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,   INTENT(INOUT) :: page  (:)   ! pollutant age since release (min)
  INTEGER,   INTENT(INOUT) :: ptyp  (:)   ! pollutant type index number
  INTEGER,   INTENT(INOUT) :: pgrd  (:)   ! meteorological grid of puff position
  INTEGER,   INTENT(IN)    :: maxpar      ! maximum particle number
  INTEGER,   INTENT(IN)    :: job_id      ! mpi implementation 
  INTEGER,   INTENT(IN)    :: num_job     ! mpi implementation
  INTEGER,   INTENT(IN)    :: ichem       ! chemistry options index
  INTEGER,   INTENT(OUT)   :: kret        ! emit denial return code
  END SUBROUTINE emspnt
!-------------------------------------------------------------------------------
  SUBROUTINE EMSTMP(SPRT,KG,NLOC,NUMTYP,KPM,INITD,DT,JET,NUMPAR,MAXPAR,      &
                    NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,      &
                    PAGE,PTYP,PGRD,job_id,num_job,kret)
  IMPLICIT NONE
  INCLUDE 'DEFSPRT.INC' ! source emissions matrix
  TYPE(qset),INTENT(IN)    :: sprt(:,:) ! source location characteristics
  INTEGER,   INTENT(IN)    :: kg        ! active grid number
  INTEGER,   INTENT(IN)    :: nloc      ! total number of source locations
  INTEGER,   INTENT(IN)    :: numtyp    ! number of pollutant types
  INTEGER,   INTENT(INOUT) :: kpm       ! number of puffs or particles
  INTEGER,   INTENT(IN)    :: initd     ! initial distribution type
  REAL,      INTENT(IN)    :: dt        ! time step (min)
  INTEGER,   INTENT(IN)    :: jet       ! current elapsed time (min)
  INTEGER,   INTENT(IN)    :: numpar    ! maximum number of particles permitted
  INTEGER,   INTENT(IN)    :: maxpar    ! maximum particle number
  INTEGER,   INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,      INTENT(INOUT) :: mass  (:,:) ! mass of pollutant (arbitrary units)
  REAL,      INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,      INTENT(INOUT) :: sigh  (:)   ! horizontal puff sigma
  REAL,      INTENT(INOUT) :: sigu  (:)   ! turbulence u'2
  REAL,      INTENT(INOUT) :: sigv  (:)   ! turbulence v'2
  REAL,      INTENT(INOUT) :: sigw  (:)   ! turbulence w'2 or vertical puff sigma
  INTEGER,   INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,   INTENT(INOUT) :: page  (:)   ! pollutant age since release (min)
  INTEGER,   INTENT(INOUT) :: ptyp  (:)   ! pollutant type index number
  INTEGER,   INTENT(INOUT) :: pgrd  (:)   ! meteorological grid of puff position
  INTEGER,   INTENT(IN)    :: job_id      ! mpi implementation
  INTEGER,   INTENT(IN)    :: num_job     ! mpi implementation
  INTEGER,   INTENT(OUT)   :: kret        ! emit denial return code
  END SUBROUTINE emstmp
!-------------------------------------------------------------------------------
  SUBROUTINE EMSMAT(SPRT,NGRD,JET,KG,KT,NLOC,NUMTYP,EFILE,QTEMP,BACK)
  use module_defgrid         ! meteorology grid and file
  IMPLICIT NONE
!  INCLUDE 'DEFGRID.INC'         ! meteorology grid and file
  INCLUDE 'DEFSPRT.INC'         ! special source-emission matrix
  TYPE(qset),    INTENT(INOUT) :: sprt(:,:) ! for each source and pollutant
  INTEGER,       INTENT(IN)    :: ngrd      ! number of meteo grids
  INTEGER,       INTENT(IN)    :: jet       ! current time
  INTEGER,       INTENT(IN)    :: kg,kt     ! active grid number
  INTEGER,       INTENT(IN)    :: nloc      ! number of source locations
  INTEGER,       INTENT(IN)    :: numtyp    ! number of pollutant types
  CHARACTER(80), INTENT(INOUT) :: efile     ! temporal emission file name
  LOGICAL,       INTENT(INOUT) :: qtemp     ! temporal file flag
  LOGICAL,       INTENT(IN)    :: back      ! integration direction
  END SUBROUTINE emsmat
!-------------------------------------------------------------------------------
  SUBROUTINE EMRISE(KDUR,HEAT,STAB,UBAR,USTR,MIXD,RISE)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: KDUR   ! emission duration (min)
  REAL(4), INTENT(IN)  :: HEAT   ! heat release (watts)
  REAL(4), INTENT(IN)  :: STAB   ! static stability (1/s2)
  REAL(4), INTENT(IN)  :: UBAR   ! wind speed (m/s)
  REAL(4), INTENT(IN)  :: USTR   ! friction velocity (m/s)
  REAL(4), INTENT(IN)  :: MIXD   ! mixed layer depth (m)
  REAL(4), INTENT(OUT) :: RISE   ! plume rise (m)
  END SUBROUTINE emrise
!-------------------------------------------------------------------------------
  SUBROUTINE LAGEMS(NUMLAG,LAGS,MAXPAR,JET,KPM,NSORT,MASS,XPOS,YPOS,ZPOS, &
                    SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD)
  IMPLICIT NONE
  INCLUDE 'DEFLAGS.INC'                    ! meteorology grid and file
  INTEGER,    INTENT(IN)    :: numlag      ! total number of source locations
  TYPE(lset), INTENT(INOUT) :: lags(:)     ! source location characteristics
  INTEGER,    INTENT(IN)    :: maxpar      ! maximum particle number
  INTEGER,    INTENT(IN)    :: jet         ! current elapsed time (min)
  INTEGER,    INTENT(INOUT) :: kpm         ! number of puffs or particles
  INTEGER,    INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,       INTENT(INOUT) :: mass  (:,:) ! mass of pollutant (arbitrary units)
  REAL,       INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,       INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,       INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,       INTENT(INOUT) :: sigh  (:)   ! horizontal puff sigma
  REAL,       INTENT(INOUT) :: sigu  (:)   ! turbulence u'2
  REAL,       INTENT(INOUT) :: sigv  (:)   ! turbulence v'2
  REAL,       INTENT(INOUT) :: sigw  (:)   ! turbulence w'2 or vertical puff sigma
  INTEGER,    INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,    INTENT(INOUT) :: page  (:)   ! pollutant age since release (min)
  INTEGER,    INTENT(INOUT) :: ptyp  (:)   ! pollutant type index number
  INTEGER,    INTENT(INOUT) :: pgrd  (:)   ! meteorological grid of puff position
  END SUBROUTINE lagems
!-------------------------------------------------------------------------------
  SUBROUTINE EMSINI(NLOC,OLAT,OLON,QFILE,NQLON,NQLAT,NPVAL)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: nloc      ! total number of source locations
  REAL,      INTENT(IN)    :: olat(:)   ! source subgrid domain limits 
  REAL,      INTENT(IN)    :: olon(:)   ! source subgrid domain limits 
  LOGICAL,   INTENT(OUT)   :: qfile     ! indicates area source emission file
  INTEGER,   INTENT(OUT)   :: nqlon     ! number of long in subgrid
  INTEGER,   INTENT(OUT)   :: nqlat     ! number of lats in subgrid
  INTEGER,   INTENT(OUT)   :: npval     ! number of pollutants in file
  END SUBROUTINE emsini
!-------------------------------------------------------------------------------
  SUBROUTINE EMSDAT(QDLON,QDLAT,NQLON,NQLAT,OLON,OLAT,QAREA,POLID)
  IMPLICIT NONE
  REAL,        INTENT(OUT) :: qdlon, qdlat             ! grid spacing
  INTEGER,     INTENT(IN)  :: nqlon, nqlat             ! number of points
  REAL,        INTENT(IN)  :: olon(:), olat(:)         ! sub grid corners
  REAL,        INTENT(OUT) :: qarea ( : , : , : , : )  ! houly emissions
  CHARACTER(4),INTENT(OUT) :: polid ( : )              ! pollutant ids 
  END SUBROUTINE emsdat
!-------------------------------------------------------------------------------
  SUBROUTINE EMSGRD(SPOT,DIRT,NUMTYP,KPM,INITD,DT,JET,NSORT,MASS,XPOS,         &
                    YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,QCYCLE,  &
                    NUMPAR,NPVAL,NQLON,NQLAT,QDLON,QDLAT,QAREA,POLID,MAXPAR)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'
  INCLUDE 'DEFSPOT.INC'
  TYPE(rset), INTENT(IN)   :: spot(:) ! source location characteristics
  TYPE(pset),INTENT(INOUT) :: dirt(:) ! for each pollutant type 
  INTEGER,   INTENT(IN)    :: numtyp  ! number of pollutant types
  INTEGER,   INTENT(INOUT) :: kpm     ! number of puffs or particles
  INTEGER,   INTENT(IN)    :: initd   ! initial distribution type
  REAL,      INTENT(IN)    :: dt      ! time step (min)
  INTEGER,   INTENT(IN)    :: jet     ! current elapsed time (min)
  REAL,      INTENT(IN)    :: qcycle  ! optional emission cycle time in hours
  INTEGER,   INTENT(IN)    :: numpar  ! maximum number of particles permitted
  INTEGER,   INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,      INTENT(INOUT) :: mass  (:,:) ! mass of pollutant (arbitrary units)
  REAL,      INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,      INTENT(INOUT) :: sigh  (:)   ! horizontal sigma   
  REAL,      INTENT(INOUT) :: sigu  (:)   ! turbulence
  REAL,      INTENT(INOUT) :: sigv  (:)   ! turbulence    
  REAL,      INTENT(INOUT) :: sigw  (:)   ! vertical sigma
  INTEGER,   INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,   INTENT(INOUT) :: page  (:)   ! pollutant age since release (min)
  INTEGER,   INTENT(INOUT) :: ptyp  (:)   ! pollutant type index number
  INTEGER,   INTENT(INOUT) :: pgrd  (:)   ! meteorological grid of puff position
  INTEGER,      INTENT(IN) :: npval            ! numb pollutants
  INTEGER,      INTENT(IN) :: nqlon, nqlat     ! number of grid points
  REAL,         INTENT(IN) :: qdlon, qdlat     ! grid spacing
  REAL,         INTENT(IN) :: qarea (:,:,:,:)  ! gridded emissions array 
  CHARACTER(4), INTENT(IN) :: polid (:)        ! emissions pollutant id
  INTEGER,      INTENT(IN) :: maxpar           ! maximum particle number
  END SUBROUTINE emsgrd
!-------------------------------------------------------------------------------
  SUBROUTINE PUFMIN(KPM,MASS,NSORT,TMASS,FRAC1,FMASS1,FRAC2,FMASS2)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)    :: kpm         ! total number of puffs or particles
  REAL,     INTENT(IN)    :: mass (:,:)  ! particle mass array
  INTEGER,  INTENT(INOUT) :: nsort (:)   ! sortted array by position
  REAL,     INTENT(IN)    :: tmass       ! total system mass
  REAL,     INTENT(IN)    :: frac1       ! fraction of total mass tmass
  REAL,     INTENT(IN)    :: frac2       ! fraction of total mass tmass
  REAL,     INTENT(OUT)   :: fmass1      ! single particle mass at frac1*tmass
  REAL,     INTENT(OUT)   :: fmass2      ! single particle mass at frac2*tmass
  END SUBROUTINE pufmin
!-------------------------------------------------------------------------------
  SUBROUTINE TRJSET(SPOT,NGRD,NTIM,NLOC,BACK,NDIAG,KVEL,IUNIT,LABEL,   &
                    job_id,num_job)
  IMPLICIT NONE
  INCLUDE 'DEFSPOT.INC'
  TYPE(rset),   INTENT(OUT):: spot(:)      ! source location characteristics
  INTEGER,      INTENT(IN) :: ngrd         ! number of grids defined
  INTEGER,      INTENT(IN) :: ntim         ! number of times defined
  INTEGER,      INTENT(IN) :: nloc         ! number of trajectories in file
  LOGICAL,      INTENT(IN) :: back         ! direction indicator
  INTEGER,      INTENT(IN) :: ndiag        ! number of diagnostic variables
  INTEGER,      INTENT(IN) :: kvel         ! vertical velocity method
  INTEGER,      INTENT(IN) :: iunit        ! model input
  CHARACTER(8), INTENT(IN) :: label(ndiag) ! diagnostic variable names
  INTEGER,      INTENT(IN) :: job_id       ! MPI job id (0 to proc-1)
  INTEGER,      INTENT(IN) :: num_job      ! MPI number of jobs
  END SUBROUTINE trjset
!-------------------------------------------------------------------------------
  SUBROUTINE TRJDSK(KP,KG,JET,MC0,IFHR,TLAT,TLON,TZHT,NDIAG,TDIAG)
  IMPLICIT NONE
  INTEGER,  INTENT(IN) :: kp            ! trajectory id sequence number
  INTEGER,  INTENT(IN) :: kg            ! current meteo grid number
  INTEGER,  INTENT(IN) :: jet           ! elapsed minutes
  INTEGER,  INTENT(IN) :: mc0           ! initial minutes
  INTEGER,  INTENT(IN) :: ifhr          ! current forecast hour
  REAL,     INTENT(IN) :: tlat,tlon     ! trajectory position
  REAL,     INTENT(IN) :: tzht          ! trajectory height AGL
  INTEGER,  INTENT(IN) :: ndiag         ! number of diagnostic meteo variables
  REAL,     INTENT(IN) :: tdiag (ndiag) ! diagnostic meteo array
  END SUBROUTINE trjdsk
!-------------------------------------------------------------------------------
  SUBROUTINE TMPLUS(IY,IM,ID,IH,IC)
  IMPLICIT NONE
  INTEGER,  INTENT(INOUT) :: iy,im,id      ! current date
  INTEGER,  INTENT(INOUT) :: ih            ! current hour
  INTEGER,  INTENT(IN)    :: ic            ! increment in hours
  END SUBROUTINE tmplus
!-------------------------------------------------------------------------------
  SUBROUTINE RUNSET(SPOT,NLOC,NHRS,NGRD,NTIM,ZMDL,ZDATA,BACK,IUNIT,KVEL)
  IMPLICIT NONE
  INCLUDE 'DEFSPOT.INC'
  TYPE(rset), INTENT(INOUT) :: spot(:) ! source location characteristics
  INTEGER,    INTENT(IN )   :: nloc    ! number of starting locations
  INTEGER,    INTENT(OUT)   :: nhrs    ! simulation duration in hours
  INTEGER,    INTENT(OUT)   :: ngrd    ! number of data grids for this run
  INTEGER,    INTENT(OUT)   :: ntim    ! number of data times for this run
  REAL,       INTENT(OUT)   :: zmdl    ! max height for scaling coord system
  REAL,       INTENT(OUT)   :: zdata   ! max height for scaling input data   
  LOGICAL,    INTENT(OUT)   :: back    ! integration direction flag
  INTEGER,    INTENT(IN )   :: iunit   ! unit for file of input parameters
  INTEGER,    INTENT(OUT)   :: kvel    ! vertical velocity remapping
  END SUBROUTINE runset
!-------------------------------------------------------------------------------
  SUBROUTINE PARPUF(KPM,PAGE,SIGU,SIGV,SIGW,HDWP,PGRD,CONAGE)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)    :: kpm        ! total number of puffs or particles
  INTEGER,  INTENT(IN)    :: page (:)   ! particle age (minutes)        
  REAL,     INTENT(INOUT) :: sigu (:)   ! horizontal u-component turbulence
  REAL,     INTENT(INOUT) :: sigv (:)   ! horizontal v-component turbulence
  REAL,     INTENT(INOUT) :: sigw (:)   ! vertical   w-component turbulence
  INTEGER,  INTENT(INOUT) :: hdwp (:)   ! horizontal distribution of pollutant
  INTEGER,  INTENT(IN)    :: pgrd (:)   ! meteorological grid index
  INTEGER,  INTENT(IN)    :: conage     ! conversion age (minutes)        
  END SUBROUTINE parpuf
!-------------------------------------------------------------------------------
  SUBROUTINE PUFPAR(KPM,PAGE,SIGU,SIGV,SIGW,HDWP,PGRD,CONAGE)
  IMPLICIT NONE
  INTEGER,  INTENT(IN)    :: kpm        ! total number of puffs or particles
  INTEGER,  INTENT(IN)    :: page (:)   ! particle age (minutes)        
  REAL,     INTENT(INOUT) :: sigu (:)   ! horizontal u-component turbulence
  REAL,     INTENT(INOUT) :: sigv (:)   ! horizontal v-component turbulence
  REAL,     INTENT(INOUT) :: sigw (:)   ! vertical   w-component turbulence
  INTEGER,  INTENT(INOUT) :: hdwp (:)   ! horizontal distribution of pollutant
  INTEGER,  INTENT(IN)    :: pgrd (:)   ! meteorological grid index
  INTEGER,  INTENT(IN)    :: conage     ! conversion age (minutes)        
  END SUBROUTINE pufpar
!-------------------------------------------------------------------------------
  SUBROUTINE PUFDSP(KPUFF,UMIX,VMIX,DTM,ZMDL,ZSFC,NLVL,WMIX,ZSG,PAGE,ZPOS,    &
                    SIGH,SIGV,HDWP,ZNDX)
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: kpuff      ! linear (0) or sqr root (1) dispersion
  REAL,    INTENT(IN)    :: umix       ! u-component turbulence (m2/s2)
  REAL,    INTENT(IN)    :: vmix       ! v-component turbulence (m2/s2)
  REAL,    INTENT(IN)    :: dtm        ! horizontal time step from advection
  REAL,    INTENT(IN)    :: zmdl       ! top of computational domain (meters)
  REAL,    INTENT(IN)    :: zsfc       ! surface terrain elevation (m)
  INTEGER, INTENT(IN)    :: nlvl       ! number of levels in subgrid
  REAL,    INTENT(IN)    :: wmix (:)   ! w-vertical turbulence profile (m2/s2)
  REAL,    INTENT(IN)    :: zsg  (:)   ! internal model sigma levels
  INTEGER, INTENT(IN)    :: hdwp       ! Horizontal distribution pollutant
  INTEGER, INTENT(IN)    :: page       ! Age of puff in minutes
  REAL,    INTENT(IN)    :: zndx       ! fraction vertical index position
  REAL,    INTENT(INOUT) :: zpos       ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh       ! horiz (meters) sigma
  REAL,    INTENT(INOUT) :: sigv       ! vert sigma (sigma)
  END SUBROUTINE pufdsp
!-------------------------------------------------------------------------------

  SUBROUTINE PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE, &
                    PTYP,PGRD,NSORT,KHMAX,FMASS)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: kpm        ! total number of puffs or particles
  INTEGER, INTENT(IN)    :: khmax      ! maximum hours limit
  REAL,    INTENT(IN)    :: fmass      ! lower mass limit per particle/puff
  REAL,    INTENT(INOUT) :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horizontal sigma   
  REAL,    INTENT(INOUT) :: sigu (:)   ! turbulence
  REAL,    INTENT(INOUT) :: sigv (:)   ! turbulence
  REAL,    INTENT(INOUT) :: sigw (:)   ! vertical sigma   
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age since release (min)
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteorological grid of puff position
  INTEGER, INTENT(INOUT) :: nsort(:)   ! sort index array by position
  END SUBROUTINE pufdel
!-------------------------------------------------------------------------------
  SUBROUTINE PUFSRT(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,XPOS,YPOS,ZPOS,MASS,  &
                    SIGH,SIGW,HDWP,PTYP,PGRD,NSORT,HGKM)
  IMPLICIT NONE
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
  REAL,    INTENT(IN)    :: hgkm (:) ! horizontal grid in km
  END SUBROUTINE pufsrt
!-------------------------------------------------------------------------------
  SUBROUTINE PUFMRG(FRACH,FRACV,FRACT,FMASS,KPM,ZMDL,MASS,XPOS,YPOS,ZPOS,  &
                    SIGH,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,HGKM)
  IMPLICIT NONE
  REAL,    INTENT(IN)    :: frach      ! horizontal position rounding fraction
  REAL,    INTENT(IN)    :: fracv      ! vertical position rounding fraction
  REAL,    INTENT(IN)    :: fract      ! travel-time (sigma) rounding fraction
  REAL,    INTENT(IN)    :: fmass      ! particle mass sort cutoff value
  INTEGER, INTENT(IN)    :: kpm        ! number of puffs or particles
  REAL,    INTENT(IN)    :: zmdl       ! model top scaling height
  REAL,    INTENT(INOUT) :: mass (:,:) ! pollutant mass
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horiz sigma (sigma)
  REAL,    INTENT(INOUT) :: sigw (:)   ! vert sigma (sigma)
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution for pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age from release 
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteo grid index for puff
  INTEGER, INTENT(IN)    :: nsort(:)   ! sortted array by position
  REAL,    INTENT(IN)    :: hgkm (:)   ! horizontal grid in km
  END SUBROUTINE pufmrg
!-------------------------------------------------------------------------------
  SUBROUTINE PUFSPV(KPM,NLVL,ZSG,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,HDWP,PAGE,  &
                    PTYP,PGRD,NSORT,NUMPAR,KRET)
  IMPLICIT NONE
  INTEGER,  INTENT(INOUT) :: kpm          ! total number of puffs or particles
  INTEGER,  INTENT(IN)    :: nlvl         ! number of levels in subgrid
  REAL,     INTENT(IN)    :: zsg   (:)    ! internal model sigma levels
  REAL,     INTENT(INOUT) :: mass  (:,:)  ! mass of pollutant (arbitrary units)
  REAL,     INTENT(INOUT) :: xpos  (:)    ! center positions (grid units)
  REAL,     INTENT(INOUT) :: ypos  (:)    ! center positions (grid units)
  REAL,     INTENT(INOUT) :: zpos  (:)    ! puff center height (sigma)
  REAL,     INTENT(INOUT) :: sigh  (:)    ! horiz sigma (sigma)
  REAL,     INTENT(INOUT) :: sigw  (:)    ! vert sigma (sigma)
  INTEGER,  INTENT(INOUT) :: hdwp  (:)    ! Horizontal distribution pollutant
  INTEGER,  INTENT(INOUT) :: page  (:)    ! pollutant age since release (min)
  INTEGER,  INTENT(INOUT) :: ptyp  (:)    ! pollutant type index number
  INTEGER,  INTENT(INOUT) :: pgrd  (:)    ! meteorological grid puff position
  INTEGER,  INTENT(INOUT) :: nsort (:)    ! sorted array index values
  INTEGER,  INTENT(IN)    :: numpar       ! number of particles in calculation
  INTEGER,  INTENT(OUT)   :: kret         ! status of splitting 0=yes 1=no
  END SUBROUTINE pufspv
!-------------------------------------------------------------------------------
  SUBROUTINE PUFSPH(KPM,HGD,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,HDWP,PAGE,PTYP,  &
                    PGRD,NSORT,NUMPAR,CGSIZE,SPLITF,KRET)
  IMPLICIT NONE
  INTEGER,  INTENT(INOUT) :: kpm          ! total number of puffs or particles
  REAL,     INTENT(IN)    :: hgd   (:)    ! meteo grid spacing (km)
  REAL,     INTENT(INOUT) :: mass  (:,:)  ! mass of pollutant (arbitrary units)
  REAL,     INTENT(INOUT) :: xpos  (:)    ! center positions (grid units)
  REAL,     INTENT(INOUT) :: ypos  (:)    ! center positions (grid units)
  REAL,     INTENT(INOUT) :: zpos  (:)    ! puff center height (sigma)
  REAL,     INTENT(INOUT) :: sigh  (:)    ! horiz sigma (sigma)
  REAL,     INTENT(INOUT) :: sigw  (:)    ! vert sigma (sigma)
  INTEGER,  INTENT(INOUT) :: hdwp  (:)    ! Horizontal distribution pollutant
  INTEGER,  INTENT(INOUT) :: page  (:)    ! pollutant age since release (min)
  INTEGER,  INTENT(INOUT) :: ptyp  (:)    ! pollutant type index number
  INTEGER,  INTENT(INOUT) :: pgrd  (:)    ! meteorological grid puff position
  INTEGER,  INTENT(INOUT) :: nsort (:)    ! sorted array index values
  INTEGER,  INTENT(IN)    :: numpar       ! number of particles in calculation
  REAL,     INTENT(IN)    :: cgsize       ! minimum size (conc/meteo) in km
  REAL,     INTENT(IN)    :: splitf       ! horizontal splitting adjustment
  INTEGER,  INTENT(OUT)   :: kret         ! status of splitting 0=yes 1=no
  END SUBROUTINE pufsph 
!-------------------------------------------------------------------------------
  SUBROUTINE DATSET(NLOC,IBYR,IBMO,IBDA,IBHR,IBMN,IUNIT,KLEN,FNAME)
  IMPLICIT NONE
  INTEGER,      INTENT(OUT)   :: nloc   ! number of starting locations
  INTEGER,      INTENT(OUT)   :: ibyr,ibmo,ibda,ibhr,ibmn
  INTEGER,      INTENT(OUT)   :: iunit  ! unit for file of input parameters
  INTEGER,      INTENT(IN)    :: klen   ! length of unique file name string
  CHARACTER(*) ,INTENT(IN)    :: fname  ! unique file name to input data
  END SUBROUTINE datset 
!-------------------------------------------------------------------------------
  SUBROUTINE ADVRNG(NGRD,MGMIN,UMAX,KPM,XPOS,YPOS,PGRD)
  use module_defgrid         ! meteorology grid and file
  IMPLICIT NONE
!  INCLUDE 'DEFGRID.INC' ! meteorology grid and file
  INTEGER,    INTENT(IN)    :: ngrd     ! minimum subgrid size (from namelist)
  INTEGER,    INTENT(IN)    :: mgmin    ! minimum subgrid size (from namelist)
  REAL,       INTENT(IN)    :: umax     ! maximum wind speed (grid pts / min)
  INTEGER,    INTENT(IN)    :: kpm      ! number of particles
  REAL,       INTENT(IN)    :: xpos (:) ! particle center positions (grid units)
  REAL,       INTENT(IN)    :: ypos (:) ! particle center positions (grid units)
  INTEGER,    INTENT(IN)    :: pgrd (:) ! particle meteorological grid 
  END SUBROUTINE advrng
!-------------------------------------------------------------------------------
  SUBROUTINE P10ADJ(P10F,ZWND,UBAR,RAIN,AREA,MASS,PGRD)
  IMPLICIT NONE
  REAL,   INTENT(IN)    :: P10F ! threshold emission sensitivity factor
  REAL,   INTENT(IN)    :: ZWND ! release height or ref meteo height
  REAL,   INTENT(IN)    :: UBAR ! wind speed in meters / sec
  REAL,   INTENT(IN)    :: RAIN ! precip rate
  REAL,   INTENT(IN)    :: AREA ! emission area in sq meters
  REAL,   INTENT(INOUT) :: MASS ! particle/puff mass
  INTEGER,INTENT(INOUT) :: PGRD ! grid identification index
  END SUBROUTINE p10adj
!-------------------------------------------------------------------------------
  SUBROUTINE PARINP(JET,KG,KT,KPM,MASS,TLAT,TLON,XPOS,YPOS,ZPOS,SIGH,SIGU,    &
                    SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,ZMDL,MAXPAR,NINIT)
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: jet        ! current elapsed time
  INTEGER, INTENT(IN)    :: kg         ! current meteo grid in space
  INTEGER, INTENT(IN)    :: kt         ! current meteo grid in time
  INTEGER, INTENT(INOUT) :: kpm        ! total number of puffs or particles
  REAL,    INTENT(INOUT) :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,    INTENT(INOUT) :: tlat (:)   ! puff center positions (lat)
  REAL,    INTENT(INOUT) :: tlon (:)   ! puff center positions (lon)
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horizontal sigma
  REAL,    INTENT(INOUT) :: sigu (:)   ! turbulence         
  REAL,    INTENT(INOUT) :: sigv (:)   ! turbulence          
  REAL,    INTENT(INOUT) :: sigw (:)   ! vertical sigma
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age since release (min)
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteorological grid of puff position
  INTEGER, INTENT(INOUT) :: nsort(:)   ! sort index array by position
  REAL,    INTENT(IN)    :: zmdl       ! vertical index scaling height
  INTEGER, INTENT(IN)    :: maxpar     ! maximum number of particles   
  INTEGER, INTENT(IN)    :: ninit      ! particle initialization type
  END SUBROUTINE parinp
!-------------------------------------------------------------------------------
  SUBROUTINE PAROUT(JET,KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,    &
                    SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,ZMDL)
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: jet        ! current elapsed time
  INTEGER, INTENT(INOUT) :: kpm        ! total number of puffs or particles
  REAL,    INTENT(INOUT) :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horizontal sigma   
  REAL,    INTENT(INOUT) :: sigu (:)   ! turbulence
  REAL,    INTENT(INOUT) :: sigv (:)   ! turbulence
  REAL,    INTENT(INOUT) :: sigw (:)   ! vertical sigma
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age since release (min)
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteorological grid of puff position
  INTEGER, INTENT(INOUT) :: nsort(:)   ! sort index array by position
  REAL,    INTENT(IN)    :: zmdl       ! vertical index scaling height
  END SUBROUTINE parout
!-------------------------------------------------------------------------------
  SUBROUTINE LAGOUT(LAGS,JET,KPM,DELT,MASS,XPOS,YPOS,ZPOS,HDWP,PAGE,   &
                    PTYP,PGRD,ZMDL)
  IMPLICIT NONE
  INCLUDE 'DEFLAGS.INC'                   ! lagrangian sampling configuration
  TYPE(lset), INTENT(INOUT) :: lags(:)    ! lagrangian sampler confuration
  INTEGER,    INTENT(IN)    :: jet        ! current elapsed time
  INTEGER,    INTENT(IN)    :: kpm        ! particle count
  REAL,       INTENT(IN)    :: delt       ! integration time step
  REAL,       INTENT(IN)    :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,       INTENT(IN)    :: xpos (:)   ! center positions (grid units)
  REAL,       INTENT(IN)    :: ypos (:)   ! center positions (grid units)
  REAL,       INTENT(IN)    :: zpos (:)   ! center height (sigma)
  INTEGER,    INTENT(IN)    :: hdwp (:)   ! always=6 as lagrangian marker
  INTEGER,    INTENT(IN)    :: page (:)   ! averaging counter
  INTEGER,    INTENT(IN)    :: ptyp (:)   ! lagrangian sampler number
  INTEGER,    INTENT(IN)    :: pgrd (:)   ! meteorological grid of puff position
  REAL,       INTENT(IN)    :: zmdl       ! vertical index scaling height
  END SUBROUTINE lagout
!-------------------------------------------------------------------------------
  SUBROUTINE LAGSUM(LAGS,CONC,NUMGRD,PLAT,PLON,ZPOS,DELT,JET,ZMDL,ZSFC,  &
                    KMSL,MASS,PAGE,PTYP,CSUM)
  IMPLICIT NONE
  INCLUDE 'DEFLAGS.INC'                   ! lagrangian sampling configuration
  INCLUDE 'DEFCONC.INC'                   ! pollutant and concentration grid
  TYPE(lset), INTENT(INOUT) :: lags(:)    ! lagrangian sampler configuration
  TYPE(cset), INTENT(IN)    :: conc(:)    ! for each concentration grid
  INTEGER,    INTENT(IN)    :: numgrd     ! number of concentration grids
  REAL,       INTENT(IN)    :: plat,plon  ! particle position lat & long
  REAL,       INTENT(IN)    :: zpos       ! center height (sigma)
  REAL,       INTENT(IN)    :: delt       ! current integration time step
  INTEGER,    INTENT(IN)    :: jet        ! current elapsed time
  REAL,       INTENT(IN)    :: zmdl       ! model domain top
  REAL,       INTENT(IN)    :: zsfc       ! height of ground surface (m)
  INTEGER,    INTENT(IN)    :: kmsl       ! agl=0 or msl=1 height units
  REAL,       INTENT(INOUT) :: mass (:)   ! mass of pollutant (arbitrary units)
  INTEGER,    INTENT(INOUT) :: page       ! averaging counter
  INTEGER,    INTENT(IN)    :: ptyp       ! lagrangian sampler number
  REAL,       INTENT(IN)    :: csum (:,:,:,:,:)  ! concentration summation
  END SUBROUTINE lagsum
!-------------------------------------------------------------------------------
  SUBROUTINE GEMPAR(KPM,HDWP,PTYP,PAGE,XPOS,YPOS,ZPOS,PGRD,MASS,CONAGE,ZMDL)
  use module_defgrid         ! meteorology grid and file
  IMPLICIT NONE
!  INCLUDE 'DEFGRID.INC'
  INTEGER,  INTENT(IN)    :: kpm        ! total number of puffs or particles
  INTEGER,  INTENT(IN)    :: hdwp (:)   ! horizontal distribution of pollutant
  INTEGER,  INTENT(IN)    :: ptyp (:)   ! pollutant type 
  INTEGER,  INTENT(IN)    :: page (:)   ! particle age (minutes)         
  REAL,     INTENT(IN)    :: xpos (:)
  REAL,     INTENT(IN)    :: ypos (:)
  REAL,     INTENT(IN)    :: zpos (:)
  INTEGER,  INTENT(INOUT) :: pgrd (:)   ! meteorological grid index
  REAL,     INTENT(INOUT) :: mass (:,:)
  INTEGER,  INTENT(IN)    :: conage     ! conversion age (minutes)         
  REAL,     INTENT(IN)    :: zmdl       ! vertical index scaling height
  END SUBROUTINE gempar
!-------------------------------------------------------------------------------
  SUBROUTINE GEMSUM (DIRT,DT,BACK,CDEP,IMO)
  IMPLICIT NONE
  INCLUDE     'DEFCONC.INC'           ! concentration and pollutant structure
  TYPE(pset),  INTENT(IN) :: DIRT (:) ! for each pollutant type
  INTEGER*4,   INTENT(IN) :: DT       ! external time step
  LOGICAL,     INTENT(IN) :: BACK     ! integration direction
  LOGICAL,     INTENT(IN) :: CDEP     ! deposition calculation
  INTEGER*4,   INTENT(IN) :: IMO      ! current month      
  END SUBROUTINE gemsum
!-------------------------------------------------------------------------------
  SUBROUTINE GEMTST (METO) 
  IMPLICIT NONE
  INCLUDE 'DEFMETO.INC'
  TYPE(aset),INTENT(IN) :: meto  ! surface and advection variables
  END SUBROUTINE gemtst
!-------------------------------------------------------------------------------
!dwen(20090315):add ext_zsg
subroutine read_zsg(fname,iu,zmdl,zsg,nzm,nlvl,ksfc,sfcl,aa,bb,cc,ret_code)
  implicit none

! Interface:
  CHARACTER(LEN=*), INTENT(IN)    :: fname               !input filename
  INTEGER         , INTENT(IN)    :: iu, nzm             !input unit number, max array dimension for zsg
  INTEGER         , INTENT(INOUT) :: nlvl, ksfc          !number of levels, sfc layer index
  REAL            , INTENT(INOUT) :: zmdl,zsg(nzm),sfcl  !model top height, model z-sigma levels, top of sfc layer
  REAL            , INTENT(INOUT) :: aa,bb,cc            !parameters for default quadratic relationship
  INTEGER, INTENT(OUT) :: ret_code !Return code; =0: ok; <0: nlvl>nzm; >0: other errors

! Automatic arrays:
  REAL                            :: hgt(nzm)
! Other local variables
  REAL                            :: zmdl_new
  INTEGER                         :: ierr, iline, k, nlvl_new, ksfc_new
  CHARACTER(LEN=256)              :: reason, line
  LOGICAL                         :: ftest
end subroutine read_zsg
!dwen(20090315):add output
SUBROUTINE OUTPUT (IFILN,IVMAX,VARSIWANT,NTURB,IETIM,KP,ICNDX,XOUT, &
                   YOUT,ZOUT,ZSFC,SIGMAW,TEMP,SAMPTTCUM,FOOTCUM,    &
                   SWF,WWOUT,ZMLAVG,RAIN,CRAI,ZFXCUM,SHTF,WHTF,     &
                   TCLD,DMASSWT,TLGR,DENSLOCAL,RHFR,SPHU,SOLW,LCLD, &
                   ZLOCNEXT,PRESLOCAL,TEMPLOCAL)

   IMPLICIT NONE

! INTENT(IN) VARIABLES
   INTEGER        , INTENT(IN)  :: IFILN,IVMAX,NTURB,IETIM,KP,ICNDX
   REAL           , INTENT(IN)  :: XOUT,YOUT,ZOUT,ZSFC,SIGMAW,TEMP,       &
                                   SAMPTTCUM,FOOTCUM,SWF,WWOUT,ZMLAVG,    &
                                   RAIN,CRAI,ZFXCUM,SHTF,WHTF,TCLD,       &
                                   DMASSWT,TLGR,DENSLOCAL,RHFR,SPHU,SOLW, &
                                   LCLD,ZLOCNEXT,PRESLOCAL,TEMPLOCAL
   CHARACTER(LEN=*), INTENT(IN) :: VARSIWANT(IVMAX)
   REAL                         :: WRITING(IVMAX)
   CHARACTER(LEN=270)           :: FORMATO
   CHARACTER(LEN=12)            :: FORMATO2
   INTEGER                      :: IV,E2,EF
end subroutine output
!dwen(20080831):add sunflx
SUBROUTINE SUNFLX(NLVL,SEA,QQ,SWF,TR)

  IMPLICIT NONE

  INTEGER,  INTENT(IN)   :: nlvl      ! number of vertical levels in profile
  REAL,     INTENT(IN)   :: sea       ! sine of the solar elevation angle 
  REAL,     INTENT(IN)   :: qq(:)     ! RH fraction profile 
  REAL,     INTENT(OUT)  :: swf       ! incident short wave flux (w/m2)
  REAL,     INTENT(OUT)  :: tr        ! transmissivity
end subroutine sunflx

!dwen(20090831):add sunang
SUBROUTINE SUNANG(JET,OLAT,OLON,EA,SEA)

  IMPLICIT NONE

  INTEGER,  INTENT(IN)   :: jet       ! elapsed minutes 
  REAL,     INTENT(IN)   :: olat      ! latitude (+ = north)
  REAL,     INTENT(IN)   :: olon      ! longitude (- = west)
  REAL,     INTENT(OUT)  :: ea        ! solar elevation angle in deg at time
  REAL,     INTENT(OUT)  :: sea       ! sine of EA

end subroutine sunang

END INTERFACE
