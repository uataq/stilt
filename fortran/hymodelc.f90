
!
! MAIN PROGRAM: HYMODELC     MAIN HYSPLIT PROGRAM FOR AIR CONCENTRATIONS
!   PRGMMR: DRAXLER          ORG: R/ARL      DATE: 1998-08-26
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   IS THE MAIN PROGRAM FOR THE TRANSPORT AND DISPERSION MODEL HYSPLIT. GIVEN
!   A TIME, SOURCE LOCATION, AND RELEASE AMOUNT, THE MODEL COMPUTES AIR
!   CONCENTRATIONS OVER PRESPECIFIED SAMPLING PERIODS, ON PRESELECTED
!   LATITUDE-LONGITUDE GRIDS AT VARIOUS HEIGHTS ABOVE GROUND.  REQUIRED
!   METEOROLOGICAL DATA ARE PRESUMED TO HAVE ALREADY BEEN PREPARED FOR MODEL
!   INPUT. SEE ROUTINES AVN2ARL FOR A DISCUSSION OF METEO DATA PREPARATION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 May 1998 (RRD)
!                 17 Aug 1998 (RRD) - added delt & isot option to namelist
!                                   - particle dump and initialization
!                 22 Dec 1998 (RRD) - included area source emission file
!                 12 Feb 1999 (RRD) - pardump option fixed for grid change
!                 04 Mar 1999 (RRD) - subgrid optimiztion
!                 20 Apr 1999 (RRD) - terrain compression factor
!                 14 May 1999 (RRD) - pgrd=0 test prior to advpnt
!                 09 Jun 1999 (RRD) - changed structure of PARDUMP file
!                 15 Jun 1999 (RRD) - added tratio to namelist
!                 26 Aug 1999 (RRD) - pass parameters for vertical grid
!                                   - puff to particle conversion
!                 08 Nov 1999 (RRD) - correction to PARDUMP read
!                 02 Dec 1999 (RRD) - added meteo subgrid min to namelist
!                 03 Mar 2000 (RRD) - advpnt subgrid fix, pardsp delt for R
!                 05 Jul 2000 (RRD) - pardump initialization generalized
!                                     passed kmsl variable to consum
!                                     relative emission start time
!                 21 Sep 2000 (RRD) - fortran90 upgrade
!                 17 Nov 2000 (RRD) - concentration output packing
!                 12 Dec 2000 (RRD) - ensemble -DENS and multiple -DMPI 
!                 02 Feb 2001 (RRD) - expanded calendar function
!                 22 Feb 2001 (RRD) - ensemble initial pressure correction
!                 01 Mar 2001 (RRD) - ensemble grid limits test
!                 09 Mar 2001 (RRD) - global lat lon grid option
!                 24 Mar 2001 (RRD) - particle initialization correction
!                 18 Jun 2001 (RRD) - temporal emission file definintion
!                 02 Oct 2001 (RRD) - simultaneous multiple meteorology
!                 24 Oct 2001 (RRD) - simplified backward dispersion
!                                   - ensemble options to namelist
!                 15 Nov 2001 (RRD) - particle dump file name
!                 17 Jan 2002 (RRD) - standardize umax units
!                 26 Feb 2002 (RRD) - downward shortwave flux for chemistry
!                 14 Mar 2002 (RRD) - mpi/ensemble single processor 
!                 22 Mar 2002 (RRD) - incorporated pm10 dust emission module
!                 18 Apr 2002 (RRD) - ensemble & mpi file label name mods 
!                 20 Jun 2002 (RRD) - starting point can be off one grid
!                 09 Sep 2002 (RRD) - initial switch to global
!                 21 Oct 2002 (RRD) - ensemble offgrid test modified
!                 03 Dec 2002 (RRD) - height initialization fix for back
!                 06 Jan 2003 (RRD) - message file time step information
!                 13 Feb 2003 (RRD) - simplified vertical grid definition
!                 10 Mar 2003 (RRD) - condsk after pardump initialization
!                 09 May 2003 (RRD) - multiple time periods in PARDUMP file
!                 22 Jul 2003 (RRD) - improved linkages with gridded emissions
!                 20 Aug 2003 (RRD) - added IER chemistry as a feature 
!                 27 Aug 2003 (RRD) - added GRS chemistry as a feature 
!                 28 Aug 2003 (RRD) - added SO2 chemistry as a feature 
!                 15 Sep 2003 (RRD) - deposition probability flag
!                 10 Nov 2003 (RRD) - support for meteo data with TKE or Vvar
!                 02 Apr 2004 (RRD) - generic file unit numbers (use funits)
!                 02 Jul 2004 (RRD) - meteo analysis to message file
!                 10 Aug 2004 (RRD) - enhanced particle-puff options
!                 13 Oct 2004 (RRD) - linear and root puff dispersion options
!                 22 Dec 2004 (RRD) - one-point concentration grids (cpack=3)
!                 03 Mar 2005 (RRD) - global boundary test for ensemble
!                 24 Mar 2005 (AFS) - incorporated cb4 chemistry options
!                 25 May 2005 (RRD) - added mass summation option (cmass=1)
!                 15 Jun 2005 (RRD) - diagnostic message file open earlier
!                 12 Oct 2005 (RRD) - Lagrangian sampling option
!                 25 Oct 2005 (RRD) - p10 emission sensitivity factor
!                 23 Feb 2006 (RRD) - initial source test now multi-grid
!                 07 Mar 2006 (RRD) - complex point source emissions
!                 12 May 2006 (RRD) - variable deallocation at termination
!                 22 May 2006 (RRD) - isotropic turbulence runs <= 12 hrs
!                 25 May 2006 (AS)  - turbulence ensembles
!                 21 Nov 2006 (RRD) - split tker into day and night components
!                 18 Jan 2007 (RRD) - maximum concentration option
!                                   - particle/puff conversion factor
!                 02 Feb 2007 (RRD) - namelist structure update correction
!                 19 Mar 2007 (RRD) - added option to compute max average 
!                 17 Apr 2007 (RRD) - multiple concentration grids MPI
!                 21 May 2007 (RRD) - updated namelist save
!                 19 Jul 2007 (RRD) - ISOT short-range default = 2
!                 17 Sep 2007 (RRD) - revised MPI start/stop test
!                 22 Oct 2007 (RRD) - standardized time step variable
!                                   - moved merge parameters to namelist
!                 14 Nov 2007 (RRD) - minor correction to puff split limits
!                 17 Jan 2008 (RRD) - plume rise diagnostics / MPI test
!                 28 Jan 2008 (RRD) - dynamic puff split/merge parameters
!                 19 Mar 2008 (RRD) - emission duration to emrise
!                 30 May 2008 (RRD) - incorporated GEM routines
!                 04 Jun 2008 (RRD) - mixing depth computation options
!                 16 Jun 2008 (BS)  - per change in condsk, affect hymodelm
!                 27 Jun 2008 (RRD) - initialized SPRT array from CONTROL
!                 01 Jul 2008 (RRD) - changed conversion test to age
!                 15 Aug 2008 (RRD) - split horizontal and vertical mixing
!                 01 Oct 2008 (RRD) - revised default namelist values
!                 10 Dec 2008 (RRD) - registered user flag kfor
!                 14 Jan 2009 (RRD) - modification to gem output initialization
!                 09 Feb 2009 (RRD) - trial or registered version response
!
! USAGE:  HYMODELC [optional process ID]
!
!   INPUT PARAMETERS: 
!         PROMPTED ON STANDARD INPUT UNLESS FILE NAMED "CONTROL" EXISTS
!   INPUT FILES:   see module funits
!   OUTPUT FILES:  see module funits
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

PROGRAM HYMODELC

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE

  !-------------------------------------------------------------------------------
  ! The following list explains all preprocessor directives that have been
  ! assigned. Subroutines with directive statements ending with .F rather than .f
  ! and should be compiled with the following flags: -WF,-Dnull, where null
  ! defines one of the directives listed below. Null skips all directives.
  !-------------------------------------------------------------------------------
  ! Executable File Name Convention - hy{c|t}{m|s}_{xxx}
  ! {c|t} where c=concentration and t=trajectory
  ! {m|s} where m=multi-processor and s=single-processor
  ! {xxx} where xxx=compilation-variation
  ! Flag  Old_name   New_name	Description
  !-------------------------------------------------------------------------------
  ! STD - hymodelc - hycs_std 	concentration single-processor
  ! MPI - hymodelm - hycm_std	concentration  multi-processor 
  ! ENS - hymodele - hycm_ens	concentration  multi-processor ensemble  
  !     - hymodels - hycs_ens	concentration single-processor ensemble
  ! VAR - hymodelv - hycs_var	concentration single-processor variance 
  ! GEM - hymodelg - hycs_gem	concentration single-processor with GEM 
  ! IER - hysp_ier - hycs_ier	concentration single-processor IER ozone 
  ! GRS - hysp_grs - hycs_grs	concentration single-processor GRS ozone
  ! SO2 - hysp_so2 - hycs_so2	concentration single-processor SO2 to SO4 
  ! CB4 - hysp_cb4 - hycs_cb4	concentration single-processor CB4 ozone
  ! PUB -          - hycp_???	public version unregistered no forecast data
  !-------------------------------------------------------------------------------

  ! registered user flag to access forecast data

  INTEGER :: KFOR = 1          ! 0=NoForecast  1=AllData

  INTEGER :: ens_id   = 1      ! first member to process
  INTEGER :: job_id   = 0      ! root job number always zero
  INTEGER :: num_job  = 1      ! defaults to single processor

  !-------------------------------------------------------------------------------
  ! special variable definitions for conditional compilation
  !-------------------------------------------------------------------------------

  INCLUDE 'DEFARG1.INC' ! main program subroutine interface
  INCLUDE 'DEFARG2.INC' ! main program subroutine interface
  INCLUDE 'DEFARG3.INC' ! global subroutine interface
  INCLUDE 'DEFCHEM.INC' ! chemistry subprogram interface
  INCLUDE 'DEFCONC.INC' ! pollutants and concentration grid
  INCLUDE 'DEFMETO.INC' ! meteo variables returned after advection
  INCLUDE 'DEFSPOT.INC' ! multiple source structure
  INCLUDE 'DEFLAGS.INC' ! lagrangian sampling structure
  INCLUDE 'DEFSPRT.INC' ! source-pollutant definition matrix

  !-------------------------------------------------------------------------------
  ! variable allocations
  !-------------------------------------------------------------------------------
  !********************************************************
  !dwen(20090305): add following variables declaration 
  !      IMPLICIT REAL (A-H,O-Z)    
  INTEGER, PARAMETER :: dp=KIND(1d0) 
  ! JCL:(11/03/03) store profile of wind errors
  TYPE WINDPROF
     !        REAL*8 U(NZM)
     REAL          :: U(mlvl)
     !   REAL, ALLOCATABLE :: u(:)         ! model sigma levels
     !         REAL*8 V(NZM)
     REAL         :: V(mlvl)
     !   REAL, ALLOCATABLE :: v(:)         ! model sigma levels
  END TYPE WINDPROF
  TYPE (WINDPROF), ALLOCATABLE :: UVERR(:)
  !      REAL*8 UERR(NZM),VERR(NZM),UERR2(NZM),VERR2(NZM)
  !      REAL*8 UUERR_T(NZM),VVERR_T(NZM)
  !      REAL*8 SIGUVERR,TLUVERR,ZCORLEN
  !      REAL*8 SIGZIERR,TLZIERR,HORCORZILEN,RELZIERR,ZIERR(MAXPAR)
  REAL, ALLOCATABLE :: uerr(:)      ! model sigma levels
  REAL, ALLOCATABLE :: verr(:)         ! model sigma levels
  REAL, ALLOCATABLE :: uerr2(:)         ! model sigma levels
  REAL, ALLOCATABLE :: verr2(:)         ! model sigma levels
  REAL, ALLOCATABLE :: uuerr_t(:)         ! model sigma levels
  REAL, ALLOCATABLE :: vverr_t(:)         ! model sigma levels
  REAL, ALLOCATABLE :: zierr(:)         ! model sigma levels
  REAL              :: siguverr,tluverr,zcorlen,sigzierr 
  REAL              :: tlzierr,horcorzilen,relzierr      
  !*********************************************************


  ! particle and puff arrays
  !                                                       particle    puff
  REAL,    ALLOCATABLE :: tlat (:)     ! position
  REAL,    ALLOCATABLE :: tlon (:)     ! position
  REAL,    ALLOCATABLE :: xpos (:)     ! position
  REAL,    ALLOCATABLE :: ypos (:)     ! position
  REAL,    ALLOCATABLE :: zpos (:)     ! height
  REAL,    ALLOCATABLE :: sigh (:)     !                Sigma_h (m) Sigma_h (m) 
  REAL,    ALLOCATABLE :: sigu (:)     !                u'2
  REAL,    ALLOCATABLE :: sigv (:)     !                v'2
  REAL,    ALLOCATABLE :: sigw (:)     !                w'2         Sigma_z (m)
  REAL,    ALLOCATABLE :: mass (:,:)   ! source mass
  INTEGER, ALLOCATABLE :: page (:)     ! age
  INTEGER, ALLOCATABLE :: hdwp (:)     ! distribution
  INTEGER, ALLOCATABLE :: ptyp (:)     ! pollutant
  INTEGER, ALLOCATABLE :: pgrd (:)     ! meteo grid
  INTEGER, ALLOCATABLE :: nsort(:)     ! sort array

  ! gridded emissions variables

  REAL,         ALLOCATABLE :: qarea (:,:,:,:)          
  CHARACTER(4), ALLOCATABLE :: polid (:)  
  REAL                      :: qdlon    ! grid spacing
  REAL                      :: qdlat    ! grid spacing
  INTEGER                   :: nqlon    ! number of long in subgrid
  INTEGER                   :: nqlat    ! number of lats in subgrid
  INTEGER                   :: npval    ! number of pollutants in file
  INTEGER                   :: nqval    ! number of time periods per day

  CHARACTER(7)              :: eflag     ! grid projection flag
  REAL,         ALLOCATABLE :: xllat (:) ! grid latitudes          
  REAL,         ALLOCATABLE :: xllon (:) ! grid longitudes         

  ! concentration output array

  REAL, ALLOCATABLE :: csum (:,:,:,:,:) ! conc (x,y,z,species,grids)
  REAL, ALLOCATABLE :: tnum (:,:,:)     ! particle total by cell
  REAL, ALLOCATABLE :: znum (:,:,:,:)   ! particles number by species

  REAL, ALLOCATABLE :: dept (:)         ! deposition accumulation
  REAL, ALLOCATABLE :: dryd (:)         ! deposition velocity 
  REAL, ALLOCATABLE :: zmass(:)         ! vertical mass diagnostic
  INTEGER           :: maxdim           ! pollutants per particle
  TYPE(cset), ALLOCATABLE :: conc(:)    ! for each concentration grid 
  TYPE(pset), ALLOCATABLE :: dirt(:)    ! for each pollutant type 

  !*****************************************************
  !dwen(20090305): add declaration of the following variables
  ! JCL:(5/18/00)array to hold WPRIME/SIGMAW for each particle exiting PARDSP
  !dwen      REAL*8 WWPREV(MAXPAR)
  REAL, ALLOCATABLE :: wwprev(:)        
  ! JCL:(09/01/03)array to store UPRIME/UERR & VPRIME/VERR for each particle--to store wind err fluctuations
  !      REAL*8 UERRPREV(MAXPAR),VERRPREV(MAXPAR)
  ! CHG(09/24/03)array to store rel. area coverage (updraft/downdraft) at particle location
  !dwen      REAL*8 AREAPRU(MAXPAR),AREAPRD(MAXPAR)
  REAL, ALLOCATABLE :: areapru(:)        
  REAL, ALLOCATABLE :: areaprd(:)        
  ! CHG(09/24/03)array to store cloud index (1=updraft, 2=env., 3=downdraft)
  !dwen     INTEGER ICNDX(MAXPAR)
  integer, ALLOCATABLE :: icndx(:)        

  ! JCL:(6/29/00)array to hold total SAMPTT (calculated in PARDSP) for each particle before results are written out
  !dwen      REAL*8 SAMPTTCUM(MAXPAR)
  !dwen      REAL*8 FOOTCUM(MAXPAR)
  REAL, ALLOCATABLE :: sampttcum(:)        
  REAL, ALLOCATABLE :: footcum(:)        

  ! CHG(24/06/04)array to store vertical displacement due to convection
  !dwen      REAL*8 ZFXCUM(MAXPAR)
  REAL, ALLOCATABLE :: zfxcum(:)        

  ! CHG:(12/05/01)array to hold 'CONVDUR' (calc. in ADVPNT, reset to 0 after conv. redistribution was done in PARDSP)
  !dwen      INTEGER CONVDUR(MAXPAR), CONVTMP
  REAL, ALLOCATABLE :: convdur(:)        
  INTEGER           :: convtmp
  !dwen      logical is_off_grid !replaces alternate return mechanism
  LOGICAL           :: is_off_grid ! replaces alternate return mechanism

  ! JCL:(4/3/02)weighting of particles due to mass violation (calculated in ADVPNT and tallied in PARDSP) for each particle
  !dwen      REAL*8 DMASSWT(MAXPAR)
  REAL, ALLOCATABLE :: dmasswt(:)        
  ! JCL:(4/3/02)the total mass violation experienced by particle in PARDSP during each timestep
  !dwen      REAL*8 DMASSTT
  REAL, ALLOCATABLE              :: dmasstt(:)        
  !dwen      real*8 zprofm(nzm) !flux-level zsg expressed as height AGL (without terrain compression)
  REAL,  ALLOCATABLE :: zprofm(:) ! flux-level zsg expressed as height AGL (without terrain compression)

  !dwen      real*8 zsigw, zsigwdn !flux-level zsg used in zprofm computation
  REAL               :: zsigw, zsigwdn !flux-level zsg used in zprofm computation

  ! JCL:variables to denote particles' positions that would be written out
  !dwen      REAL*8 XOUT,YOUT,ZOUT
  REAL               :: xout,yout,zout

  ! JCL:standard deviation of vertical & horizontal velocity that is calculated in PARDSP
  !dwen      REAL*8 SIGMAW
  !dwen      REAL*8 SIGMAU
  REAL               :: sigmaw, sigmau

  ! JCL:amount of time [min.] that particle 'sees' the ground--calculated in PARDSP:
  !      SAMPTT=(# particle touchdowns)*(deltat); deltat is timestep of subloop in PARDSP
  !dwen      REAL*8 SAMPTT
  !dwen      REAL*8 FOOT
  REAL               :: samptt, foot

  ! CHG:used to store Lagrangian timescale--calculated in PARDSP
  !dwen      REAL*8 SAMPTT2
  REAL               :: samptt2

  ! JCL:variables used to store wind shear [(m/s)/m]
  !dwen      REAL*8 DUDZ, DVDZ
  REAL               :: dudz, dvdz

  !dwen      INTEGER SMIN
  INTEGER            :: smin
  Logical            :: preset

  ! JCL:logical flag to specify whether data should be written out to PARTICLE.DAT or not
  !dwen      LOGICAL WRITEOUT
  LOGICAL             :: writeout

  ! JCL:flag to tell whether particle has 'seen' vegetation on ground
  !dwen      INTEGER SEEVEG
  INTEGER            :: seeveg

  ! JCL:random seed for random number generator--would always be
  !     '5' if the flag 'RANDOM' is selected to be '0'; otherwise
  !     the random seed itself is 'randomized' using the current time & RANDOM_NUMBER
  INTEGER              :: RSEED, dtt(8) 
  REAL                 :: harvest

  ! JCL:the number of seconds since midnight--used to generate random seed
  !dwen      REAL*8 SECS
  REAL                :: secs

  ! JCL:declre the variables read in as NAMELIST
  !dwen      INTEGER INITD,KHMAX,NUMPAR,QCYCLE,KRND,ISOT,NDUMP,RANDOM
  INTEGER             :: random
  !dwen      REAL*8 FRMR,DELT

  ! JCL:(2/28/2004) enable user to choose DYNAMICALLY the output variables
  CHARACTER(4), DIMENSION(100) :: VARSIWANT
  INTEGER                          :: IVMAX

  ! JCL:(9/16/02)flag to say whether PBL height was prescribed or not
  !dwen      INTEGER ZICONTROLTF
  INTEGER                          :: zicontroltf
  ! JCL:(9/16/02)vector to store the scaling factors to prescribe PBL height
  !dwen(20090823)      REAL*8 ZIPRESC(150)
  REAL                             ::  ZIPRESC(150)

  ! JCL:(09/01/03)flag specifying whether to include wind errors as Markov process
  !dwen      INTEGER WINDERRTF
  INTEGER                           :: winderrtf

  ! JCL:counter to keep track of # of [min] since last time output written to PARTICLE.DAT
  !dwen      REAL*8 COUNTOUT
  REAL                               :: countout

  !     assumed maximum wind speed in km/min for calculation of subgrid size
  !dwen      REAL(dp) :: UMAXI=1.2
  REAL                               :: umaxi=1.2

  ! JCL:(2/13/2001)variable to store Wbar (vertical velocity)
  !dwen      REAL*8 WWOUT
  REAL                              :: wwout

  ! JCL:(4/27/01)counter to keep track of how many particles have left model area
  !dwen      INTEGER COUNTNPAROUT
  INTEGER                           :: countnparout

  !dwen      REAL(dp) :: TLON=-HUGE(TLON), TLAT=-HUGE(TLON)

  ! CHG(09/16/03):add flag specifying whether data from RAMS or not
  !dwen      LOGICAL RAMSFLG, awrfflg, ECMFLG
  LOGICAL                           :: ramsflg,awrfflg,ecmflg
  !dwen       INTEGER :: dummy
  INTEGER                           :: dummy,nturb,iconvect
  REAL                              :: tlfrac,outdt,veght,outfrac
  !dwen(20090315): define variables                                            
  integer                           :: iconvectrams,nhrszi,    &
       knzi,kpp,kzz,kk,kcb1,kcb2,      &
       kdb1,kct1,kct2,kdt1,icndx1,icndx2,     &
       ktt,kbb
  real                              :: vmax,ea,sea,max_hgt,horcorlen,deltaz,  &
       rauto,uu,vv,z1,raup,area,z2,   &
       convfacx,convfacy,uu2,vv2,     &
       uuerrprev,vverrprev,cfacx,     &
       cfacy,xposprev,yposprev,       &
       rautoh,uuerrnext,vverrnext,    &
       dmass,swf,tr,zf,preslocal,     &
       templocal,denslocal,rhfrlocal, &
       sphu


  !******************************************************
  ! meteorological data file and arrays

  CHARACTER(80)     :: fname, nbase     ! file information
  REAL, ALLOCATABLE :: zsg  (:)         ! model sigma levels

  TYPE(qset), ALLOCATABLE :: sprt (:,:) ! source-pollutant definition matrix
  TYPE(lset), ALLOCATABLE :: lags (:)   ! lagrangian sampling information
  TYPE(rset), ALLOCATABLE :: spot (:)   ! source location characteristics
  TYPE(bset), ALLOCATABLE :: metz (:)   ! profile advection variables
  TYPE(aset)              :: meto       ! surface advection variables

  ! simulation flags

  LOGICAL           :: back             ! integration direction
  LOGICAL           :: cdep             ! deposition process
  LOGICAL           :: rdep             ! resistance method
  LOGICAL           :: sdep             ! resuspension
  !dwen(20090315)  LOGICAL           :: ftest1,ftest2    ! generic file test
  LOGICAL           :: ftest,ftest1,ftest2    ! generic file test
  LOGICAL           :: qtemp = .false.  ! temporal emission file
  LOGICAL           :: qfile            ! gridded emission file
  LOGICAL           :: lagsam           ! lagrangian sampling file
  LOGICAL           :: spliton = .true. ! puff split status

  ! rounding and merging parameters
  REAL     :: frhs,frvs,frts,frhe,frve,frte,frhmax,splitf

  REAL     :: tkerd,tkern,ubar,tmass,dt,pmass,frac,p10f  
  REAL     :: pcnt,rmass,height,fmass,sgb,zsfc,cgsize,sfcl,zmdl,zdata 
  REAL     :: zx,zz,sgt,umax,tratio,delt,cc,aa,bb,frmr,frme,qcycle 
  REAL     :: dist, xtmp, ytmp, ztmp, dxf, dyf, dzf
  REAL     :: cntr, avgrise, avgmixd, vscale, hscale, tvmix  

  INTEGER  :: kret,kspk,kemit,kmixd,kmix0,conage,kzmix,kbls,kblt,initk
  INTEGER  :: isot,kagl,nver,hdwpx,kpuff,mhrs,numlag,kavrg,kmaxc,ksnap,kspl
  INTEGER  :: ktime,nstep,kz,ks,kt,kh,mm,maxdt,ifhr,ida,ihr,imo,numpol
  INTEGER  :: iyr,imn,kp,kpt,i,j,k,k1,k2,iunit,ibda,ibhr,ibmn,nstr,ichem,cmass 
  INTEGER  :: kdef,mgmin,nloc,ibyr,numspl,ktyp,kmsl,klen,initd,khmax,cpack 
  INTEGER  :: n,jet,kg,numtyp,nzp,kpm,nyp,kgrid,numgrd,nxp,nhrs,ngrd,ntim
  INTEGER  :: ninit,ncycl,ndump,numpar,krnd,maxpar,ksfc,ibmo,nlvl,kvel,mc0 
  INTEGER  :: TOUT,TM_PRES,TM_TPOT,TM_TAMB,TM_RAIN,TM_MIXD,TM_RELH 
  INTEGER  :: TM_TERR,TM_DSWF 
  integer :: max_nlvl

  CHARACTER(80) :: ecode ! memory allocation error message
  CHARACTER(80) :: efile ! temporal emission file name
  CHARACTER(80) :: pinpf ! particle input file
  CHARACTER(80) :: poutf ! particle output file

  INTEGER, PARAMETER :: NUMVAR=27
  INTEGER            :: VX(NUMVAR)

  INTEGER :: ret_code

  !-------------------------------------------------------------------------------
  ! external definitions
  !-------------------------------------------------------------------------------
  ! meteorolgical grid, record, and file information

  ! parameters that define vertical grid polynomial
  COMMON /ZZTOKK/ AA,BB,CC

  ! mixing length scales
  COMMON /stblen/ vscale,hscale

  ! special simulation setup parameters read from namelist
  !  NAMELIST/SETUP/INITD,KHMAX,NUMPAR,MAXPAR,MAXDIM,QCYCLE,P10F,KMIXD,KMIX0,     &
  !                 FRME,FRMR,KRND,DELT,KDEF,TKERD,TKERN,NDUMP,NCYCL,TRATIO,ISOT, &
  !                 MGMIN,KMSL,NSTR,MHRS,NVER,CPACK,CMASS,ICHEM,EFILE,VSCALE,     &
  !                 TOUT,TM_PRES,TM_TPOT,TM_TAMB,TM_RAIN,TM_MIXD,TM_TERR,HSCALE,  &
  !                 TM_DSWF,TM_RELH,DXF,DYF,DZF,PINPF,POUTF,NINIT,KAGL,KPUFF,     &
  !                 KBLT,KZMIX,TVMIX,KBLS,CONAGE,KSPL,FRHS,FRVS,FRTS,FRHMAX,SPLITF
  !
  !***************************************************************
  !dwen(20090305): modify namelist according to the requirement of STILT
  ! special simulation setup parameters read from namelist
  !  NAMELIST/SETUP/INITD,KHMAX,NUMPAR,MAXPAR,MAXDIM,QCYCLE,P10F,KMIXD,KMIX0,     &
  !                 FRME,FRMR,KRND,DELT,KDEF,TKERD,TKERN,NDUMP,NCYCL,TRATIO,ISOT, &
  !                 MGMIN,KMSL,NSTR,MHRS,NVER,CPACK,CMASS,ICHEM,EFILE,VSCALE,     &
  !                 TOUT,TM_PRES,TM_TPOT,TM_TAMB,TM_RAIN,TM_MIXD,TM_TERR,HSCALE,  &
  !                 TM_DSWF,TM_RELH,DXF,DYF,DZF,PINPF,POUTF,NINIT,KAGL,KPUFF,     &
  !                 KBLT,KZMIX,TVMIX,KBLS,CONAGE,KSPL,FRHS,FRVS,FRTS,FRHMAX,SPLITF &
  !                 random,aa,bb,cc,tlfrac,outdt,nturb,veght,outfrac,zicontroltf, &
  !                 iconvect,winderrtf,ivmax,umaxi,varsiwant

  !dwen(20090310):remove aa,bb,cc from setup namelist,because their values are assigned 
  !               in metlvl.f in hysplit4.9
  NAMELIST/SETUP/INITD,KHMAX,NUMPAR,MAXPAR,MAXDIM,QCYCLE,P10F,KMIXD,KMIX0,     &
       FRME,FRMR,KRND,DELT,KDEF,TKERD,TKERN,NDUMP,NCYCL,TRATIO,ISOT, &
       MGMIN,KMSL,NSTR,MHRS,NVER,CPACK,CMASS,ICHEM,EFILE,VSCALE,     &
       TOUT,TM_PRES,TM_TPOT,TM_TAMB,TM_RAIN,TM_MIXD,TM_TERR,HSCALE,  &
       TM_DSWF,TM_RELH,DXF,DYF,DZF,PINPF,POUTF,NINIT,KAGL,KPUFF,     &
       KBLT,KZMIX,TVMIX,KBLS,CONAGE,KSPL,FRHS,FRVS,FRTS,FRHMAX,SPLITF, &
       random,tlfrac,outdt,nturb,veght,outfrac,zicontroltf, &
       iconvect,winderrtf,ivmax,umaxi,varsiwant

  !******************************************************************


  !-------------------------------------------------------------------------------
  ! required for NCEP operational implementation
  !-------------------------------------------------------------------------------

  ! CALL W3TAGB('HYMODELC',1998,0238,0068,'R/ARL  ')

  ! unique input file information added as suffix to CONTROL.{FNAME}
  ! standard output files get same extension
  CALL GETARG(1,FNAME)

  ! identify as registered or trial version if argument = version
  IF(FNAME(1:7).EQ.'version')THEN
     IF(KFOR.EQ.0)THEN
        !       trial version
        WRITE(*,'(A)')'trial'
        STOP
     ELSE
        !       registered version
        WRITE(*,'(A)')'registered'
        STOP
     END IF
  END IF

  ! constant seed value
  VX=-1

  KLEN=INDEX(FNAME,' ')-1

  ! open diagnostic message file with unique name if required
  IF(KLEN.LE.0)THEN
     OPEN(KF21,FILE='MESSAGE')
  ELSE
     OPEN(KF21,FILE='MESSAGE.'//FNAME(1:KLEN))
  END IF

  ! title line to standard output
  WRITE(*,*)'HYSPLIT49 (Feb 2009) - Initialization ',FNAME(:KLEN)

  !***************************************************************
  !dwen(20090306):add this line
  ! JCL:initalize output flag to FALSE
  WRITEOUT=.FALSE.
  !***************************************************************

  !-------------------------------------------------------------------------------
  ! trajectory model namelist parameters (not used in this code)
  !-------------------------------------------------------------------------------

  METO%FLAG%PRES=0 
  METO%FLAG%TPOT=0
  METO%FLAG%TAMB=0
  METO%FLAG%RAIN=0
  METO%FLAG%MIXD=0
  METO%FLAG%RELH=0
  METO%FLAG%TERR=0
  METO%FLAG%DSWF=0
  TM_PRES=0
  TM_TPOT=0
  TM_TAMB=0
  TM_RAIN=0
  TM_MIXD=0
  TM_TERR=0
  TM_DSWF=0
  TM_RELH=0
  NSTR=0
  MHRS=0
  NVER=0
  TOUT=0

  !-------------------------------------------------------------------------------
  ! concentration model default namelist parameters read from file: SETUP.CFG  
  !-------------------------------------------------------------------------------
  !dwen(20090306):set default to 0
  !     isotropic turbulence option (0 - off; 1 - on)
  ISOT=-99         ! no longer used
  !
  !-------------------------------------------------------------------------------

  !dwen(20090306):set default to 4
  !   INITD=0      ! initial distribution:
  INITD=4      ! initial distribution:
  !                  0 - 3D particle model
  !                  1 - 3D puff model: horizontal Gaussian, vertical top-hat
  !                  2 - 3D puff model: horizontal top-hat,  vertical top-hat 
  !                  3 - Puff/Particle: horizontal Gaussian, vertical particle
  !                  4 - Puff/Particle: horizontal top-hat,  vertical particle
  !                  5 - deposited particles subject to 2D sfc-water transport
  !                  6 - defines particle as a Lagrangian sampler
  !                103 - converts 3D particle #0 to #3
  !                104 - converts 3D particle #0 to #4
  !                130 - converts #3 to #0 3D particle
  !                140 - converts #4 to #0 3D particle

  CONAGE=48    ! particle/puff conversions at conage (hours)

  KPUFF=0      ! linear (0) or square root (1) horizontal puff dispersion
  KHMAX=9999   ! maximum duration in hours of any particle/puff
  NUMPAR=2500  ! number of puffs in simulation or particles to release
  MAXPAR=100000! maximum number of particles carried in simulation
  MAXDIM=1     ! maximum number of pollutants to carry on one mass particle
  QCYCLE=0.0   ! optional cycling of emissions (hours) anatex=60

  FRMR=0.0     ! mass removal fraction during enhanced merging
  FRME=0.10    ! mass rounding fraction for enhanced merging

  KSPL=1       ! standard splitting interval (hours)
  FRHS=1.00    ! standard horizontal puff rounding fraction for merge
  FRVS=0.01    ! vertical puff rounding
  FRTS=0.10    ! temporal puff rounding

  KRND=6       ! enhanced merge interval (hours)
  FRHE=1.5*FRHS
  FRVE=1.5*FRVS
  FRTE=1.5*FRTS

  SPLITF=1.0   ! =1 automatic size adjustment factor for horizontal splitting
  ! <0 disable 
  ! >1 use value

  FRHMAX=3.0   ! maximum value for the horizontal rounding parameter

  TKERD=0.18   ! day (unstable) turbulent kinetic energy ratio 
  TKERN=0.18   ! night (stable) turbulent kinetic energy ratio
  ! ratio = w'2/(u'2+v'2); night>day sets 40% urban enhancement
  ! setting a non-zero value forces that ratio at all levels

  KDEF=1       ! horizontal turbulence computation; default as for ISOT=0
  ! 0 - in proportion to the vertical turbulence
  ! 1 - computed from the velocity deformation

  KZMIX=1      ! vertical mixing adjustments; default as for ISOT=0
  ! 0 - NONE vertical diffusivity in PBL varies with height  
  ! 1 - Vertical diffusivity in PBL single average value
  ! 2 - scale boundary-layer values by TVMIX
  ! 3 - scale free-troposphere values by TVMIX 

  TVMIX=1.0    ! vertical mixing scale factor

  KBLS=1       ! boundary layer stability derived from ; default as for ISOT=0
  ! 1 - heat and momentum fluxes
  ! 2 - wind and temperature profiles

  KBLT=1       ! boundary layer turbulence parameterizations; default as for ISOT=0
  ! 1 - Beljaars/Holtslag and Betchov/Yaglom
  ! 2 - Kanthar/Clayson
  ! 3 - TKE field from input meteorology data file 
  ! 4 - velocity variances from input meteorology

  VSCALE=200.0   ! vertical Lagrangian time scale (sec)
  HSCALE=10800.0 ! horizontal Lagrangian time scale (sec)

  ! Note on defaults used if needed input data are missing:
  !   if ZICONTROLTF=FALSE: will use temperature profile if data are missing for
  !                         KMIXD=0 or KMIXD=2
  !   if ZICONTROLTF=TRUE with a single negative ziscale:
  !      will use whatever is specified by KMIXD, then overwrite with model if available
  KMIXD=3      ! 0 = use input data MIXD if available, otherwise compute; default=3 as in STILT
  ! 1 = compute from temperature profile (also used as default for KMIXD=0 or 2 if data is missing)
  ! 2 = compute from TKE profile
  ! 3 = compute from bulk Ri profile - added for STILT
  ! > = 10 use this value as a constant

  KMIX0=250    ! minimum mixing depth (abs(kmix0) is used as the minimum mixing depth,
               ! negative values are used to force mixing heights coincident with model levels)

  NINIT=1      ! particle initialization (0-none; 1-once; 2-add; 3-replace)
  NDUMP=0      ! dump particles to/from file  0-none or nhrs-output interval 
  NCYCL=0      ! pardump output cycle time (added to initial output)
  TRATIO=0.75  ! advection stability ratio
  DELT=0.0     ! integration time step (0 - autoset; >0 - constant minutes)
  MGMIN=10     ! minimum meteorological subgrid size
  KMSL=0       ! starting heights default to AGL=0 or MSL=1

  CPACK=1      ! binary concentration packing 
  ! 0 = all grid points written to file
  ! 1 = only nonzero points written to file
  ! 2 = special non-regular grid 

  CMASS=0      ! compute grid concentrations (0) or grid mass (1)

  ICHEM=0      ! special chemistry or conversion modules:
  ! 1 = concentration grid treated as source-receptor matrix format
  ! 2 = convert pollutant from species #1 to species #2 
  ! 3 = enable pm10 dust storm emission module 
  ! 4 = configure concentration grid similar to meteorology grid 
  ! 5 = treat 3D particle deposition using probability function
  ! 6 = (not used)
  ! 7 = enable water surface transport of particle deposition

  P10F=1.0     ! dust threshold velocity sensitivity factor
  EFILE=' '    ! temporal emission file name for input on unit=52
  DXF=1.0      ! horizontal grid adjustment factor for ensemble]
  DYF=1.0        
  DZF=0.01     ! vertical factor for ensemble (0.01 ~ 250m)

  PINPF='PARINIT' ! particle dump input file
  POUTF='PARDUMP' ! particle dump output file

  !******************************************************************
  !dwen(20090306):add default values for the following variables
  ! JCL:default random seed is 5
  RSEED=5
  ! JCL:default of RANDOM is TRUE  (0-FALSE; 1-TRUE)
  RANDOM=1
  ! JCL:fraction of TL (Lagrangian timescale) to set as timestep in dispersion subroutine
  TLFRAC=0.1
  ! JCL:default value of OUTDT is 0.0--means output is written EVERY timestep
  OUTDT=0.0
  ! JCL:ht [fraction of PBL ht or m] below which a particle's time spent is tallied;
  !        If <=1.0, then specifies fraction of PBL ht
  VEGHT=0.5
  ! CHG:default value of NTURB is FALSE means TURBULENCE is on
  NTURB=0
  ! JCL:(4/27/01)default value of OUTFRAC--specifies fraction of particle number over which
  !        model would stop if fraction of particles that leaves model area exceeds this number
  !        e.g., if OUTFRAC=0.9, once over 90% of particles leave model area, then model stops
  OUTFRAC=0.9
  ! JCL:(9/16/02) initialize PBL height prescription flag to FALSE (0-FALSE; 1-TRUE)
  ZICONTROLTF=0
  ! CHG:(9/17/02)default value of ICONVECT is FALSE ~means CONVECTION is off
  ICONVECT=0
  ICONVECTRAMS=0
  ! JCL:(09/01/03) wind err flag to FALSE (0-FALSE; 1-TRUE)
  ! JCL:(09/01/03) wind err flag to FALSE (0-FALSE; 1-TRUE, 2: zi error, 3: zi and wind error)
  WINDERRTF=0
  ! JCL:(02/28/2004) default output variables to be written to 'PARTICLE.DAT'
  !     number of variables to output
  IVMAX=5
  !     4-lettered codes for the different variables
  VARSIWANT(:)=' '
  VARSIWANT(1:5)=(/'TIME','INDX','LONG','LATI','ZAGL'/)
  !     vertical grid polynomial defaults Z= aa*k^2 + bb*k + cc
  !     mechanism: if set in SETUP.CFG - no further change, otherwise set in sbr runset depending
  !        on meteo model ID

  !dwen(20090309):comment out AA,BB,CC because they are assigned in subroutine metlvl
  !      AA = -HUGE(AA)
  !      BB = AA
  !      CC = AA
  !
  !******************************************************************

  ! open namelist file (use special file suffix if required)

  IF(KLEN.LE.0)THEN
     NBASE='SETUP.CFG'
  ELSE
     NBASE='SETUP.'//FNAME(1:KLEN)
  END IF

  ! read namelist file if found in /working directory
  INQUIRE(FILE=NBASE,EXIST=FTEST1)
  IF(FTEST1)THEN
     OPEN(KF26,FILE=NBASE)
     READ(KF26,SETUP)
     CLOSE(KF26)
     WRITE(*,*)' NOTICE   main: using namelist file - ',NBASE
  END IF

  !*****************************************************************
  !dwen(20090306):add the following line
  ! JCL:(02/28/2004) make sure that not select too many output variables
  IF(IVMAX.GT.100) STOP 'Too many output vars (IVMAX) in SETUP.CFG!'
  !*****************************************************************

  IF(ICHEM.EQ.5.OR.ICHEM.EQ.7)THEN
     IF(INITD.NE.0)THEN
        WRITE(*,*)' NOTICE   main: Probability deposition requires INITD=0'
        WRITE(*,*)' Configuration corrected, continuing ...'
        INITD=0
     END IF

     IF(MAXDIM.NE.1)THEN
        WRITE(*,*)' *ERROR*  main: Probability deposition requires MAXDIM=1'
        WRITE(*,*)' Configuration not corrected ...'
        STOP 900   
     END IF
  END IF

  ! trajectory option (kmsl=2) not valid
  IF(KMSL.LT.0.OR.KMSL.GT.1)THEN
     WRITE(*,*)'WARNING main: BL fraction height not valid in concentration run'
     KMSL=0
  END IF

  IF(num_job.GT.1.AND.(INITD.EQ.1.OR.INITD.EQ.2))THEN
     WRITE(*,*)'WARNING main: Puff dispersion (initd=1,2) not optimized for MPI'
     WRITE(*,*)'Simulation will continue on one processor!'
  END IF

  ! horizontal split factor, if turned off, insure value set to one
  IF(SPLITF.LT.0.0) SPLITF=-1.0  

  ! mixing depth options limits
  KMIX0=MAX(-999,MIN( 999,KMIX0))
  if (KMIX0 .eq. 0) KMIX0=1
  KMIXD=MAX(0,MIN(9999,KMIXD))

  ! vertical mixing adjustments
  KZMIX=MAX(0,MIN(3,KZMIX))

  ! legacy ISOT option
  IF(ISOT.NE.-99)THEN
     IF(ISOT.EQ.0)THEN
        KBLS=1
        KDEF=1
        KZMIX=1
        KBLT=1
     ELSEIF(ISOT.EQ.1)THEN
        KBLS=1
        KDEF=0
        KZMIX=0
        KBLT=2
     ELSEIF(ISOT.EQ.2)THEN
        KBLS=2
        KDEF=0
        KZMIX=0
        KBLT=2
     ELSEIF(ISOT.EQ.3)THEN
        KBLS=1
        KDEF=1
        KZMIX=0
        KBLT=1
     ELSEIF(ISOT.EQ.4)THEN
        KBLT=3
     ELSEIF(ISOT.EQ.5)THEN
        KBLT=3
        KMIXD=2
     ELSEIF(ISOT.EQ.6)THEN
        KBLT=4 
     END IF
     ISOT=-99
  END IF

  !-------------------------------------------------------------------------------
  ! allocate particle arrays
  !-------------------------------------------------------------------------------

  ALLOCATE (xpos(maxpar),ypos(maxpar),zpos(maxpar),STAT=kret)
  IF(kret.ne.0)THEN 
     ECODE='Particle position'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  ALLOCATE (sigh(maxpar),sigu(maxpar),sigv(maxpar),sigw(maxpar),STAT=kret)    
  IF(kret.ne.0)THEN      
     ECODE='Particle sigma'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  ALLOCATE (mass(maxdim,maxpar), STAT=kret)
  IF(kret.ne.0)THEN     
     ECODE='Particle mass'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  ALLOCATE (page(maxpar),hdwp(maxpar),ptyp(maxpar),STAT=kret)
  IF(kret.ne.0)THEN      
     ECODE='Particle characteristics'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  ALLOCATE (pgrd(maxpar),nsort(maxpar),STAT=kret)
  IF(kret.ne.0)THEN      
     ECODE='Particle grid and sort'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  !dwen********************************************************
  !dwen(20090305):allocate variable
  ALLOCATE (zierr(maxpar),STAT=kret)
  IF(kret.ne.0)THEN      
     ECODE='Particle sigma err profile'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF
  !********************************************************

  !dwen********************************************************
  !dwen(20090305):allocate variables
  ALLOCATE (wwprev(maxpar),uverr(maxpar),                    &
       areapru(maxpar),areaprd(maxpar),   &
       icndx(maxpar),                     &
       sampttcum(maxpar),footcum(maxpar),zfxcum(maxpar),  &
       convdur(maxpar),dmasswt(maxpar),STAT=kret)
  IF(kret.ne.0)THEN      
     ECODE='Particle sigma err profile'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF
  !********************************************************



  !-------------------------------------------------------------------------------
  ! standard meteorological model initialization
  !-------------------------------------------------------------------------------

  ! time function initialization
  CALL TMINIT

  ! set number of locations and starting time
  CALL DATSET(NLOC,IBYR,IBMO,IBDA,IBHR,IBMN,IUNIT,KLEN,FNAME)

  ALLOCATE (spot(nloc), STAT=kret)    ! defined for each source location
  IF(kret.ne.0)THEN     
     ECODE='Source location array'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  SPOT(1)%IBYR=IBYR
  SPOT(1)%IBMO=IBMO
  SPOT(1)%IBDA=IBDA
  SPOT(1)%IBHR=IBHR
  SPOT(1)%IBMN=IBMN

  ! set the basic simulation control parameters

  !*********************************************************************
  !dwen(20090309):use updated RUNSET
  !               HYSPLIT4.5 added SPOT to argument list, 
  !                  instead of using COMMON bloc in old verison of EMSGRD
  !               HYSPLIT4.6 added NTIM (number of meteo data times defined in CONTROL)
  !                  to argument list and initialized it in RUNSET
  !               HYSPLIT4.7 moved KSFC,SFCL and NLVL to METLVL
  !               ZSG is calculated in main program after HYSPLIT4.6, not in RUNSET any more
  !               HYSPLIT4.5 removed FNAME and KLEN in RUNSET
  !                      and IUNIT is assigned in DATSET 
  !               AA,BB,CC are assigned in subroutine METLVL after HYSPLIT4.6

  CALL RUNSET(SPOT,NLOC,NHRS,NGRD,NTIM,ZMDL,ZDATA,BACK,IUNIT,KVEL)
  !               
  !      CALL RUNSET(NLOC,NHRS,NGRD,NLVL,KSFC,ZSG,ZMDL,SFCL,               &
  !                  BACK,IUNIT,KVEL,KLEN,FNAME, AA, BB, CC)
  !*********************************************************************

  ! initialize information about each meteo data grid
  !*********************************************************************
  !dwen(20090314): added CFLG,TCLF,LCLF,RADF,SLMF,
  !                     CPPT,CPPD,CPP3,CPP6,CPRC in METINIT according to STILT
  !dwen(20090317): HYSPLIT4.5 moved HEADER to subroutine METSET that is called by METINI 
  !dwen(20090317): HYSPLIT4.9 added KDEF,KZMIX,TVMIX,KBLS,KBLT to argument list, they are initialized by default or 
  !                   assigned in SETUP.CFG file, used here for consistent check of meteo data
  !                   with turbulence options in SETUP namelist
  !dwen(20090317): HYSPLIT4.6 added NTIM (number of meteo data times defined in CONTROL) to argument list
  !
  CALL METINI(KDEF,KZMIX,TVMIX,KBLS,KBLT,NGRD,NTIM,SPOT(1)%OLAT,SPOT(1)%IBYR, &
       SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR,SPOT(1)%IBMN,BACK,KVEL,KFOR)

  !       CALL METINI(HEADER,NGRD,SPOT(1)%OLAT,SPOT(1)%IBYR,                &
  !     &   SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR,BACK,KVEL)
  !*********************************************************************
  WRITE(KF21,*)' NOTICE   main: number meteo grids and times - ',NGRD,NTIM

  !-------------------------------------------------------------------------------
  ! define the resolution of the internal vertical meteo grid
  !-------------------------------------------------------------------------------

  CALL METLVL(NGRD,ZDATA,NLVL,KSFC,SFCL)

  if (nlvl .gt. mlvl) then
     ! Add this here because loop limits below based on nlvl, for arrays dimensioned by mlvl
     write (*,*) 'Diagnosed nlvl=',nlvl,' greater than array limit mlvl=',mlvl
     write (*,*) 'Rerun with smaller ZDATA in CONTROL, or recompile with larger mlvl'
     stop 900
  end if

  ALLOCATE (zsg(nlvl),STAT=kret)
  IF(kret.ne.0)THEN      
     ECODE='Sigma profile - initial'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  DO N=1,NLVL
     !    compute levels as sigma coordinate for terrain=0
     !    correction applied in each subroutine as needed
     ZSG(N)=1.0-(AA*FLOAT(N*N)+BB*FLOAT(N)+CC)/ZMDL
  END DO

  !-------------------------------------------------------------------------------
  ! configure the starting location, grid, and time
  !-------------------------------------------------------------------------------
  !************************************************************
  !dwen(20090306):add the following lines
  ! CHG(09/18/03) moved from after setting JET=MC0
  !     set default meteorological grid number for starting particles
  KG=1
  !     set current active meteorological grid ( when no particles = 0 )
  !dwen(20090824)      KGRID=0
  Kt=1

  ! CHG(09/09/03) Need to use excact flux levels from RAMS
  ! Assume no change in vertical grid between nests
  ! i.e. set NLVL (number of internal levels),
  ! ZSG (pseudo sigma internal coordinate),
  ! ZMDL (to corresponding value from RAMS at (NZRAMS-1)),
  ! KSFC (vertical index for top of sfc layer, was 2, now 1), SFCL (=HGT(KSFC)),
  IF(GRID(KG,kt)%MODEL_ID.EQ.'RAMS')THEN
     ! CHG&JCL (03/10/2004) set convection flag for RAMS
     !     so that other routines are not affected
     ICONVECTRAMS=ICONVECT
     ! CHG&JCL (03/10/2004) no excess convection for RAMS
     IF(ICONVECT.EQ.1)ICONVECT=0
     RAMSFLG = .TRUE.
     NLVL=GRID(KG,kt)%NZ-2
     deallocate(zsg)
     ALLOCATE (zsg(nlvl),STAT=kret)
     IF(kret.ne.0)THEN      
        ECODE='Sigma profile - BRAMS branch'
        WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
        STOP 900
     END IF
     ZMDL=DREC(KG,kt)%HEIGHT(NLVL+1)
     ZSG(1:NLVL)=1.0-DREC(KG,kt)%HEIGHT(2:(NLVL+1))/ZMDL
     KSFC=1
     SFCL=DREC(KG,kt)%HEIGHT(KSFC+1)
  ELSE
     RAMSFLG = .FALSE.
  ENDIF
  awrfflg = GRID(KG,kt)%MODEL_ID(2:4) .EQ. 'WRF'
  ECMFLG  = GRID(KG,kt)%MODEL_ID(1:2) == 'EC'
  if (awrfflg .OR. ECMFLG) then
     if (iconvect .eq. 1) then !turn on cgrell, turn off extreme conv. (as for RAMS)
        iconvectrams = 1
        iconvect = 0
     elseif (iconvect .eq. -1) then !turn off cgrell, turn on extreme conv.
        iconvectrams = 0
        iconvect = 1
     endif
  end if

  ! Optionally read in AGL heights of model layers from external file:
  if (.not. RAMSFLG) then
     write (*,*) 'Initial call to read_zsg to retrieve nlvl'
     max_nlvl=1
     call read_zsg ('ZSG_LEVS.IN',KFJCLMSG,zmdl,zsg,max_nlvl,nlvl,ksfc,sfcl,aa,bb,cc,ret_code)
     if (ret_code .eq. 0) then
        write (*,*) 'ERROR hymodelc: Unexpected successfull return for max_nlvl=',max_nlvl
        STOP 900
     elseif (ret_code .eq. -1) then
        ! Retrieved nlvl, reallocate array and read them in:
        deallocate(zsg)
        ALLOCATE (zsg(nlvl),STAT=kret)
        IF(kret.ne.0)THEN      
           ECODE='Sigma profile - read_zsg branch'
           WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
           STOP 900
        END IF
        write (*,*) 'Second call to read_zsg to read zsg levels'
        max_nlvl=nlvl
        call read_zsg ('ZSG_LEVS.IN',KFJCLMSG,zmdl,zsg,max_nlvl,nlvl,ksfc,sfcl,aa,bb,cc,ret_code)
     else
        ! Cannot use external zsg levels, leave things be
     end if
     !         call read_zsg ('ZSG_LEVS.IN',45,zmdl,zsg,nzm,nlvl,ksfc,sfcl,aa,bb,cc)
  endif

  ALLOCATE (metz(nlvl),STAT=kret)  
  IF(kret.ne.0)THEN      
     ECODE='Meteo return profile'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  ALLOCATE (zmass(nlvl),STAT=kret) 
  IF(kret.ne.0)THEN     
     ECODE='Level diagnostic information'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  !*******************************************************
  !dwen(20090305):allocate variables
  ALLOCATE (uerr(nlvl),verr(nlvl),uerr2(nlvl),                    &
       verr2(nlvl),uuerr_t(nlvl),vverr_t(nlvl),zprofm(nlvl), &
       STAT=kret) 
  IF(kret.ne.0)THEN     
     ECODE='wind sigma err profile'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF
  !*******************************************************


  ! CHG(09/16/03) moved from runset
  ! check starting height limit against model top, using ZSG and ZMDL
  max_hgt = zmdl * (1.-zsg(nlvl))
!!$      IF(RAMSFLG)THEN
  DO N=1,NLOC
     !         starting height limit
!!$        IF(SPOT(N)%OLVL.GT.ZMDL)THEN
     IF(SPOT(N)%OLVL.GT.max_hgt)THEN
        WRITE(*,*)'ERROR hymodelc: Start height above mdl domain'
        WRITE(*,*)'   Source - ',N,'    Height - ',SPOT(N)%OLVL
!!$          WRITE(*,*)'   Model domain - ',INT(DREC(KG)%HEIGHT(NLVL+1))
        WRITE(*,*)'   Model domain - ',max_hgt
        WRITE(*,*)'   Model top ht - ',ZMDL
        STOP
     END IF
  END DO
  !************************************************************


  ! convert starting time to accumulated minutes
  CALL TM2MIN(SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR,  &
       SPOT(1)%IBMN,MC0)

  ! elapsed time
  JET=MC0

  ! same time applies to all source locations
  IF(NLOC.GT.1)THEN
     DO N=2,NLOC
        spot(n)%ibyr=spot(1)%ibyr
        spot(n)%ibmo=spot(1)%ibmo
        spot(n)%ibda=spot(1)%ibda
        spot(n)%ibhr=spot(1)%ibhr 
        spot(n)%ibmn=spot(1)%ibmn 
     END DO
  END IF

  ! set initial meteorological grid (time & grid index)
  KT=1
  KG=1

  ! default initial meteo grid (km) for time-step calculation
  ! the first grid defined always has the finest resolution
  CGSIZE=GRID(KG,KT)%SIZE

  !***************************************************************
  !dwen(20090306):moved from STILT
  ! JCL:
  OPEN(KFJCLMSG,FILE='JCLmessage',FORM='FORMATTED')

  ! JCL:open file to write particle information
  OPEN(KFPARDAT,FILE='PARTICLE.DAT',FORM='FORMATTED')

  ! get a random number RSEED to initialize RAN3, if RANDOM=1
  ! JCL(050803): implemented new random seed initialization suggested by Stefan Koerner
  IF(RANDOM.EQ.1)THEN
     CALL RANDOM_SEED                       ! initializes the build-in random number generator
     CALL DATE_AND_TIME (values=dtt)
     DO i=1,MAX(dtt(8),1)                   ! dtt(8) contains the milliseconds of the wall-clock
        CALL RANDOM_NUMBER (harvest)          ! the build-in random number generator for uniform
     END DO                                 !    distribution between 0 and 1
     RSEED = ABS(NINT((harvest-0.5)*1e4))
  END IF

  ! JCL:
  WRITE(KFJCLMSG,*) 'RSEED: ',RSEED

  ! JCL:(3/16/01) Instead of setting initial WWPREV (WPRIME normalized by stddev) to 0.0,
  !               set it to a GAUSSIAN DISTRIBUTION
  ! JCL:(6/29/00) Initialized the cumulative SAMPTT to 0.0
  ! JCL:(4/3/02)  Initialized the weights of particles from mass violation [fraction of gridcell]
  ! CHG(09/24/03) Initialize rel. area coverage up/downdraft at particle location
  ! CHG(09/24/03) Initialize cloud index to environment
  DO I=1,MAXPAR
     !dwen(20090315):
     !         WWPREV(I)=GASDEV(RSEED,DBLE(1.0))
     WWPREV(I)=GASDEV(RSEED,1.0)
     SAMPTTCUM(I)=0.0
     FOOTCUM(I)=0.0
     DMASSWT(I)=1.0
     AREAPRU(I)=0.0
     AREAPRD(I)=0.0
     !start in environment
     ICNDX(I)=2
     ! CHG(24/06/04) add ZFXCUM as vertical displacement due to convective flux
     !  (deep or shallow, up or downdraft) along trajectory [m]
     ZFXCUM(I)=0.0
  END DO
  !***************************************************************


  ! set initial meteorological grid index for each starting location
  DO N=1,NLOC

     KG=1       
     SPOT(N)%KG=0 
     newkg : DO WHILE (KG.LE.NGRD)

        IF(GRID(KG,KT)%LATLON)THEN
           !****************************************************************
           !dwen(20090309): HYSPLIT4.6 added KT (active time number) to argument list to 
           !              considere simultaneous multiple meteo 

           CALL GBL2XY(KG,KT,SPOT(N)%OLAT,SPOT(N)%OLON,                        &
                SPOT(N)%XP,SPOT(N)%YP)
           !            CALL GBL2XY(KG,SPOT(N)%OLAT,SPOT(N)%OLON,                   &
           !     &                     SPOT(N)%XP,SPOT(N)%YP)
           !****************************************************************
        ELSE
           !****************************************************************
           !dwen(20090309): HYSPLIT4.6 added KT (active time number) to argument list to 
           !              considere simultaneous multiple meteo 

           CALL CLL2XY_wps(GRID(KG,KT)%GBASE,SPOT(N)%OLAT,SPOT(N)%OLON,            &
                SPOT(N)%XP,SPOT(N)%YP,GRID(KG,KT)%proj)
           !            CALL CLL2XY(GRID(KG)%GBASE,SPOT(N)%OLAT,SPOT(N)%OLON,       &
           !     &                  SPOT(N)%XP, SPOT(N)%YP, GRID(KG)%proj)
           !****************************************************************
        END IF

        IF(SPOT(N)%XP.LT.2.0.OR.SPOT(N)%XP.GT.FLOAT(GRID(KG,KT)%NX)-1.0.OR.    &
             SPOT(N)%YP.LT.2.0.OR.SPOT(N)%YP.GT.FLOAT(GRID(KG,KT)%NY)-1.0)THEN
           IF(GRID(KG,KT)%GBLDAT)THEN
              MGMIN=MAX(GRID(KG,KT)%NX,GRID(KG,KT)%NY)
              SPOT(N)%ZP=(ZMDL-SPOT(N)%OLVL)/ZMDL
              SPOT(N)%KG=KG
              EXIT newkg
           ELSE
              WRITE(*,*)'WARNING main: source point off grid ',KG,N
              WRITE(*,*)'Position:',SPOT(N)%OLAT,SPOT(N)%OLON
              WRITE(*,*)'Grid Loc:',SPOT(N)%XP,  SPOT(N)%YP
           END IF
           KG=KG+1     

        ELSE
           !          default without terrain adjustment
           SPOT(N)%ZP=(ZMDL-SPOT(N)%OLVL)/ZMDL
           !          meteorological grid
           SPOT(N)%KG=KG
           EXIT newkg
        END IF

     END DO newkg

  END DO

  ! eliminate off-grid starting locations
  N=1
  K=0  
  DO WHILE (N.LE.NLOC)
     IF(SPOT(N)%KG.NE.0)THEN
        K=K+1
        IF(K.LT.N) spot(k)=spot(n)
     END IF
     N=N+1
  END DO

  IF(K.EQ.0)THEN     
     WRITE(*,*)'*ERROR* main: all source points off all meteo grids!'
     STOP 900
  END IF
  IF(K.LT.NLOC)THEN 
     NLOC=K
     WRITE(*,*)'WARNING main: number of sources reduced to - ',NLOC
  END IF
  KG=1  ! reset to initial file 

  !**********************************************************************
  !dwen(20090306):copied from STILT
  ! JCL:write starting location to PARTICLE.DAT
  !     use any old variable--ZOUT in this case, to get formatting correct
  ! JCL:(11/6/02) not needed anymore
  WRITE(KFPARDAT,*)SPOT(1)%OLAT,SPOT(1)%OLON,SPOT(1)%OLVL

  ! JCL:(9/16/02) read in prescribed scaling factors for mixed-layer height
  IF(ZICONTROLTF.EQ.1)THEN
     INQUIRE(FILE='ZICONTROL',EXIST=FTEST)
     IF(FTEST)THEN
        OPEN(KF26,FILE='ZICONTROL',ACTION='READ')
        READ(KF26,*)NHRSZI
        !         WRITE(45,*)'NHRSZI=',NHRSZI
        DO KNZI=1,NHRSZI
           READ(KF26,*)ZIPRESC(KNZI)
           !            WRITE(45,*)KNZI,' ZIPRESC=',ZIPRESC(KNZI)
        END DO
        WRITE(KFJCLMSG,'(a,i4,a,g15.6)') ' ZI prescribed from ZICONTROL: nhrszi= ', &
             & nhrszi,' zipresc(1)= ',ZIPRESC(1)
        CLOSE(KF26)
     ELSE
        WRITE(*,*)'File ZICONTROL not found!'
        STOP
     END IF
  END IF
  !**********************************************************************

  !-------------------------------------------------------------------------------
  ! emission (source term) and pollutant initializaton
  !-------------------------------------------------------------------------------

  NUMTYP = 1
  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Number of different pollutants'
     WRITE(*,*)NUMTYP
  END IF
  READ(IUNIT,*)NUMTYP
  IF(IUNIT.EQ.5)WRITE(KF22,*)NUMTYP

  ! matrix modification tranforms multiple locations to multiple pollutants
  NUMPOL=NUMTYP
  IF(ICHEM.EQ.1)NUMPOL=NLOC

  ALLOCATE (dirt(numpol),STAT=kret)  
  IF(kret.ne.0)THEN      
     ECODE='Pollutant description'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  ALLOCATE (dept(numpol),dryd(numpol),STAT=kret)
  IF(kret.ne.0)THEN     
     ECODE='Deposition description'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  !***************************************************************
  !dwen(20090306):use updated subroutine EMSSET (updated after HYSPLIT4.5)
  !dwen(20090306): HYSPLIT4.5 removed DIRT from COMMON bloc and added to argument list
  !       HYSPLIT4.6 added BACK(backward option) to argument list
  !       HYSPLIT4.6 added OLAT and OLON in argument list

  CALL EMSSET(DIRT,NUMTYP,IUNIT,SPOT(1)%IBYR,SPOT(1)%IBMO, SPOT(1)%IBDA,  &
       SPOT(1)%IBHR,SPOT(1)%OLAT,SPOT(1)%OLON,BACK)
  !     emission (source term) initializaton
  !      CALL EMSSET(NUMTYP,IUNIT,                                         &
  !     &   SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR)
  !***************************************************************

  ! check dimensions for special applications
  IF(MAXDIM.GT.NUMTYP) WRITE(*,*)    &
       'WARNING: namelist file (maxdim) > control file (numtyp)'

  IF(ICHEM.GT.0)THEN
     IF(ICHEM.EQ.1) &
          WRITE(*,*)' NOTICE   main: concentration file treated as matrix'
     IF(ICHEM.EQ.2) &
          WRITE(*,*)' NOTICE   main: simple particle species conversion enabled'
     IF(ICHEM.EQ.3) &
          WRITE(*,*)' NOTICE   main: PM10 dust emission algorithm enabled'
     IF(ICHEM.EQ.4) &
          WRITE(*,*)' NOTICE   main: concentration grid set to meteorology grid'
     IF(ICHEM.EQ.5) &
          WRITE(*,*)' NOTICE   main: particle dry deposition using probability' 
     IF(ICHEM.EQ.6) &
          WRITE(*,*)' NOTICE   main: invalid ichem=6 option ignored' 
     IF(ICHEM.EQ.7) THEN
        WRITE(*,*)' NOTICE   main: particle transport on water surfaces' 
        IF(.NOT.DREC(KG,KT)%UFLX)THEN
           WRITE(*,*)' *ERROR*  main: input data with momentum fluxes required'
           STOP 900
        END IF
     END IF
  END IF

  ! input file for point source emissions
  IF(EFILE.NE.'')THEN
     ALLOCATE (sprt(nloc,numtyp),STAT=kret)  
     IF(kret.ne.0)THEN      
        ECODE='Source-Pollutant Matrix (from EFILE)'
        WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
        STOP 900
     END IF

     !    first guess default fill in from CONTROL file information
     DO N=1,NLOC
        SPRT(N,:)%START=DIRT(:)%START%MACC
        SPRT(N,:)%STOP =DIRT(:)%START%MACC+DIRT(:)%QHRS*60
        SPRT(N,:)%KG   =SPOT(N)%KG
        SPRT(N,:)%XP   =SPOT(N)%XP
        SPRT(N,:)%YP   =SPOT(N)%YP
        SPRT(N,:)%QLVL =SPOT(N)%OLVL
        IF(SPOT(N)%QTRM.GT.0.0)THEN
           SPRT(N,:)%RATE=SPOT(N)%QTRM   ! from source record  
        ELSE
           SPRT(N,:)%RATE=DIRT(:)%QRATE  ! from pollutant record 
        END IF
        SPRT(N,:)%AREA =SPOT(N)%AREA
        SPRT(N,:)%HEAT =0.0           
     END DO
  END IF

  !-------------------------------------------------------------------------------
  ! sampling (concentration) grid initialization (update CGSIZE if required)
  !-------------------------------------------------------------------------------

  NUMGRD = 1
  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Number of simultaneous concentration grids'
     WRITE(*,*)NUMGRD
  END IF
  READ(IUNIT,*)NUMGRD
  IF(IUNIT.EQ.5)WRITE(KF22,*)NUMGRD

  ALLOCATE (conc(numgrd),STAT=kret)  
  IF(kret.ne.0)THEN     
     ECODE='Concentration grid structure'
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF

  !****************************************************************
  ! JCL: added BACK as argument to be able to assign concentration
  !      sampling starting & stopping times accurately
  !     sampling grid initialization (update CGSIZE if required)
  !dwen(20090309): HYSPLIT4.5 added CONC to argument list, removed from COMMON bloc
  !dwen(20090814): HYSPLIT4.9 added CPACk for concentration packing flag
  !               

  CALL CONSET(CONC,ZMDL,NUMGRD,IUNIT,CGSIZE,CPACK,SPOT(1)%OLAT,SPOT(1)%OLON,   &
       SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR,BACK)

  !      CALL CONSET(ZMDL,NUMGRD,IUNIT,CGSIZE,                             &
  !     &   SPOT(1)%OLAT,SPOT(1)%OLON,                                     &
  !     &   SPOT(1)%IBYR,SPOT(1)%IBMO,SPOT(1)%IBDA,SPOT(1)%IBHR,BACK)
  !*****************************************************

  ! set the maximum and snapshot grid numbers
  ! definition of a maximum concentration grid requires a snapshot grid defined
  KMAXC=0
  KSNAP=0
  KAVRG=0
  DO K=1,NUMGRD
     IF(CONC(K)%SNAP.EQ.0) KAVRG=K
     IF(CONC(K)%SNAP.EQ.1) KSNAP=K
     IF(CONC(K)%SNAP.EQ.2) KMAXC=K
  END DO

  ! set loop indicies and allocate array space
  NXP=MAXVAL(CONC(:)%NUMB_LON)
  NYP=MAXVAL(CONC(:)%NUMB_LAT)
  NZP=MAXVAL(CONC(:)%LEVELS)

  ! chem=4 sets conc array to the same size as the meteo data array
  IF(ICHEM.EQ.4)THEN
     NXP=MAX(NXP,GRID(KG,KT)%NX)
     NYP=MAX(NYP,GRID(KG,KT)%NY)
     IF(NGRD.GT.1)THEN
        WRITE(*,*)' *ERROR* main: ichem=4 option only supports one meteo grid'
        STOP 900
     END IF
  END IF

  ECODE='STD Concentration array'
  ALLOCATE (csum (nxp,nyp,nzp,numpol,numgrd),STAT=kret)

  IF(kret.ne.0)THEN       
     WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
     STOP 900
  END IF
  csum = 0.0

  ! set the deposition parameters for each pollutant type

  !*************************************************
  !dwen(20090306):  HYSPLIT4.5 removed DIRT from COMMON bloc and added to argument list

  CALL DEPSET(DIRT,NUMTYP,IUNIT,CDEP,RDEP,SDEP)

  !      CALL DEPSET(NUMTYP,IUNIT,CDEP,RDEP,SDEP)
  !*************************************************

  ! matrix modification tranforms multiple locations to multiple pollutants
  IF(ICHEM.EQ.1)THEN
     !     reset pollutant type number as one for each source location 
     NUMTYP=NLOC
     KTYP=1
     WRITE(DIRT(1)%IDENT,'(I4.4)') KTYP

     !     define polluant parameter array for remaining source locations
     DO KTYP=2,NUMTYP
        DIRT(KTYP)=DIRT(1)
        WRITE(DIRT(KTYP)%IDENT,'(I4.4)') KTYP
     END DO
  END IF

  ! dust storm module requires deposition to be turned on for precip field
  IF(ICHEM.EQ.3) CDEP=.TRUE.

  ! open concentration output files (units 21,...)

  !***********************************************
  !dwen(20090306):  HYSPLIT4.5 removed SPOT,DIRT and CONC from COMMON bloc and added to argument list
  !dwen(20090814):  HYSPLIT4.5 removed CSUM 
  !dwen(20090814):  HYSPLIT4.9 added CPACk for concentration packing flag
  !dwen(20090306):  open conc output file as FORMATTED according to STILT code, instead of UNFROMATTED

  CALL CONINI(SPOT,CONC,DIRT,NLOC,NUMGRD,NUMTYP,CPACK)

  !     open concentration output files (units 21,...), zero variables
  !      CALL CONINI(NLOC,NUMGRD,NUMTYP,CSUM)
  !***********************************************

  !-------------------------------------------------------------------------------
  ! concentration grid sensitive switches 
  !-------------------------------------------------------------------------------

  ! special limit for splitting puffs to leave room in the array for subsequent
  ! emissions and further splitting

  INITK=INITD
  IF(INITK.GE.100) INITK=MOD(INITD/10,10)

  IF(INITK.EQ.1.OR.INITK.EQ.2)THEN
     !    leave room in the array for one emission per minute 
     NUMSPL=MAX(1,MAXPAR-INT(60.0*NLOC*DIRT(1)%QHRS))
  ELSEIF(INITK.EQ.3.OR.INITK.EQ.4)THEN
     NUMSPL=MAXPAR/2
  ELSE
     NUMSPL=MAXPAR
  END IF

  !-------------------------------------------------------------------------------
  ! diagnostic file initialization
  !-------------------------------------------------------------------------------

  ! create backup namelist for reference
  OPEN(KF26,FILE='CONC.CFG')
  WRITE(KF26,'(A)')' &SETUP'
  WRITE(KF26,'(A,F4.2,A)')' tratio = ',tratio,','

  IF(delt.LT.10.0)THEN
     WRITE(KF26,'(A,F3.1,A)')' delt = ',delt,','
  ELSE
     WRITE(KF26,'(A,F4.1,A)')' delt = ',delt,','
  END IF

  IF(initd.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' initd = ',initd,','
  ELSE
     WRITE(KF26,'(A,I3,A)')' initd = ',initd,','
  END IF

  WRITE(KF26,'(A,I1,A)')' kpuff = ',kpuff,','

  IF(khmax.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' khmax = ',khmax,','
  ELSEIF(khmax.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' khmax = ',khmax,','
  ELSEIF(khmax.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' khmax = ',khmax,','
  ELSE
     WRITE(KF26,'(A,I4,A)')' khmax = ',khmax,','
  END IF

  IF(numpar.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' numpar = ',numpar,','
  ELSEIF(numpar.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' numpar = ',numpar,','
  ELSEIF(numpar.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' numpar = ',numpar,','
  ELSEIF(numpar.LT.10000)THEN
     WRITE(KF26,'(A,I4,A)')' numpar = ',numpar,','
  ELSEIF(numpar.LT.100000)THEN
     WRITE(KF26,'(A,I5,A)')' numpar = ',numpar,','
  ELSEIF(numpar.LT.1000000)THEN
     WRITE(KF26,'(A,I6,A)')' numpar = ',numpar,','
  ELSE
     WRITE(KF26,'(A,I7,A)')' numpar = ',numpar,','
  END IF

  IF(maxpar.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' maxpar = ',maxpar,','
  ELSEIF(maxpar.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' maxpar = ',maxpar,','
  ELSEIF(maxpar.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' maxpar = ',maxpar,','
  ELSEIF(maxpar.LT.10000)THEN
     WRITE(KF26,'(A,I4,A)')' maxpar = ',maxpar,','
  ELSEIF(maxpar.LT.100000)THEN
     WRITE(KF26,'(A,I5,A)')' maxpar = ',maxpar,','
  ELSEIF(maxpar.LT.1000000)THEN
     WRITE(KF26,'(A,I6,A)')' maxpar = ',maxpar,','
  ELSE
     WRITE(KF26,'(A,I7,A)')' maxpar = ',maxpar,','
  END IF

  IF(qcycle.LT.10.0)THEN
     WRITE(KF26,'(A,F3.1,A)')' qcycle = ',qcycle,','
  ELSEIF(qcycle.LT.100.0)THEN
     WRITE(KF26,'(A,F4.1,A)')' qcycle = ',qcycle,','
  ELSE
     WRITE(KF26,'(A,F5.1,A)')' qcycle = ',qcycle,','
  END IF

  WRITE(KF26,'(3A)')' efile = ',''''//trim(efile),''','

  WRITE(KF26,'(A,I1,A)')' kdef = ',kdef,','
  WRITE(KF26,'(A,I1,A)')' kzmix = ',kzmix,','
  WRITE(KF26,'(A,I1,A)')' kbls = ',kbls,','
  WRITE(KF26,'(A,I1,A)')' kblt = ',kblt,','
  WRITE(KF26,'(A,I3,A)')' isot = ',isot,','

  WRITE(KF26,'(A,F5.1,A)')' vscale = ',vscale,','
  WRITE(KF26,'(A,F7.1,A)')' hscale = ',hscale,','

  WRITE(KF26,'(A,F4.2,A)')' tvmix = ',tvmix,','
  WRITE(KF26,'(A,F4.2,A)')' tkerd = ',tkerd,','
  WRITE(KF26,'(A,F4.2,A)')' tkern = ',tkern,','

  IF(kmix0 .ge. 0 .and. kmix0 .LT. 10)THEN
     WRITE(KF26,'(A,I1,A)')' kmix0 = ',kmix0,','
  ELSEIF(kmix0 .gt. -10 .and. kmix0.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' kmix0 = ',kmix0,','
  ELSEIF(kmix0 .gt. -100 .and. kmix0.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' kmix0 = ',kmix0,','
  ELSE
     WRITE(KF26,'(A,I4,A)')' kmix0 = ',kmix0,','
  END IF

  IF(kmixd.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' kmixd = ',kmixd,','
  ELSEIF(kmixd.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' kmixd = ',kmixd,','
  ELSEIF(kmixd.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' kmixd = ',kmixd,','
  ELSE
     WRITE(KF26,'(A,I4,A)')' kmixd = ',kmixd,','
  END IF

  WRITE(KF26,'(A,I1,A)')' ninit = ',ninit,','
  IF(ndump.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' ndump = ',ndump,','
  ELSEIF(ndump.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' ndump = ',ndump,','
  ELSEIF(ndump.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' ndump = ',ndump,','
  ELSE
     WRITE(KF26,'(A,I4,A)')' ndump = ',ndump,','
  END IF

  IF(ncycl.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' ncycl = ',ncycl,','
  ELSEIF(ncycl.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' ncycl = ',ncycl,','
  ELSEIF(ncycl.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' ncycl = ',ncycl,','
  ELSE
     WRITE(KF26,'(A,I4,A)')' ncycl = ',ncycl,','
  END IF

  WRITE(KF26,'(3A      )')' pinpf = ',''''//trim(pinpf),''','
  WRITE(KF26,'(3A)')' poutf = ',''''//trim(poutf),''','

  IF(mgmin.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' mgmin = ',mgmin,','
  ELSEIF(mgmin.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' mgmin = ',mgmin,','
  ELSEIF(mgmin.LT.1000)THEN
     WRITE(KF26,'(A,I3,A)')' mgmin = ',mgmin,','
  ELSE
     WRITE(KF26,'(A,I4,A)')' mgmin = ',mgmin,','
  END IF

  IF(conage.LT.10)THEN
     WRITE(KF26,'(A,I1,A)')' conage = ',conage,','
  ELSEIF(conage.LT.100)THEN
     WRITE(KF26,'(A,I2,A)')' conage = ',conage,','
  ELSE
     WRITE(KF26,'(A,I3,A)')' conage = ',conage,','
  END IF

  WRITE(KF26,'(A,I1,A)')' kmsl = ',kmsl,','
  WRITE(KF26,'(A,I1,A)')' ichem = ',ichem,','
  WRITE(KF26,'(A,I1,A)')' cpack = ',cpack,','
  WRITE(KF26,'(A,I1,A)')' cmass = ',cmass,','
  WRITE(KF26,'(A,I1,A)')' kspl = ',kspl,','    
  WRITE(KF26,'(A,I1,A)')' krnd = ',krnd,',' 
  WRITE(KF26,'(A,F4.2,A)')' frhmax = ',frhmax,',' 
  WRITE(KF26,'(A,F4.2,A)')' splitf = ',splitf,',' 
  WRITE(KF26,'(A,F4.2,A)')' frhs = ',frhs,',' 
  WRITE(KF26,'(A,F4.2,A)')' frvs = ',frvs,','
  WRITE(KF26,'(A,F4.2,A)')' frts = ',frts,','
  WRITE(KF26,'(A,F4.2,A)')' dxf = ',dxf,','
  WRITE(KF26,'(A,F4.2,A)')' dyf = ',dyf,','
  WRITE(KF26,'(A,F4.2,A)')' dzf = ',dzf,','
  WRITE(KF26,'(A)')' /'
  CLOSE(KF26)

  ! close all initialization files
  IF(IUNIT.EQ.5)THEN
     CLOSE(KF22)
  ELSE
     CLOSE(KF25)
  END IF

  WRITE(KF21,*)'HYSPLIT49 (Feb 2009)'
  WRITE(KF21,*)' '
  WRITE(KF21,*)'------------- Start Namelist configuration -------------------'
  WRITE(KF21,*)'Internal grid parameters (nlvl,aa,bb,cc): ',NLVL,AA,BB,CC
  WRITE(KF21,SETUP)
  WRITE(KF21,*)'------------- End Namelist configuration -------------------'
  WRITE(KF21,*)' '

  WRITE(*,*)'Calculation Started ... please be patient'

  ! diagnostic message on status of deposition flags

  WRITE(KF21,*)' NOTICE   main: pollutant initialization flags'
  IF(DIRT(1)%DOGAS)WRITE(KF21,*)'      Gas pollutant - ',DIRT(1)%DOGAS
  IF(DIRT(1)%DODRY)WRITE(KF21,*)'     Dry deposition - ',DIRT(1)%DODRY
  IF(DIRT(1)%DOWET)WRITE(KF21,*)'     Wet deposition - ',DIRT(1)%DOWET
  IF(DIRT(1)%DORES)WRITE(KF21,*)'        Dynamic dry - ',DIRT(1)%DORES
  IF(DIRT(1)%DOGRV)WRITE(KF21,*)'      Grav settling - ',DIRT(1)%DOGRV
  IF(DIRT(1)%DOSUS)WRITE(KF21,*)'       Resuspension - ',DIRT(1)%DOSUS
  IF(DIRT(1)%DORAD)WRITE(KF21,*)'        Radioactive - ',DIRT(1)%DORAD
  IF(DIRT(1)%DOVOL)WRITE(KF21,*)'Volume Unit Convert - ',DIRT(1)%DOVOL

  WRITE(KF21,*)' NOTICE   main: meteorological data flags'
  WRITE(KF21,*)'    Ten meter wind - ',DREC(KG,KT)%UFLG
  WRITE(KF21,*)'    Two meter temp - ',DREC(KG,KT)%TFLG
  WRITE(KF21,*)'    Exchange coeff - ',DREC(KG,KT)%EFLX
  WRITE(KF21,*)'         Heat flux - ',DREC(KG,KT)%HFLX
  WRITE(KF21,*)'     Momentum flux - ',DREC(KG,KT)%UFLX
  WRITE(KF21,*)' Velocity variance - ',DREC(KG,KT)%VELV
  WRITE(KF21,*)'      3D TKE field - ',DREC(KG,KT)%TKEN
  WRITE(KF21,*)'    Down shortwave - ',DREC(KG,KT)%DSWF
  WRITE(KF21,*)' Friction velocity - ',DREC(KG,KT)%USTR
  WRITE(KF21,*)'     Friction temp - ',DREC(KG,KT)%TSTR
  WRITE(KF21,*)'    Terrain height - ',DREC(KG,KT)%SHGT
  WRITE(KF21,*)'  Surface pressure - ',DREC(KG,KT)%PRSS
  WRITE(KF21,*)' Mixed layer depth - ',DREC(KG,KT)%MIXD 
  WRITE(KF21,*)' Profile averaging - ',DREC(KG,KT)%KZMIX
  WRITE(KF21,*)'  Stability method - ',DREC(KG,KT)%KBLS 
  WRITE(KF21,*)' Horizontal mixing - ',DREC(KG,KT)%KDEF 
  WRITE(KF21,*)' PBL mixing scheme - ',DREC(KG,KT)%KBLT 
  WRITE(KF21,*)'  Free Trop Mixing - ',DREC(KG,KT)%TVMIX  
  WRITE(KF21,*)' Precip accumulate - ',DREC(KG,KT)%ACYCLE
  WRITE(KF21,*)'------------- Start computation messages -------------------'

  !-------------------------------------------------------------------------------
  ! Initialize gridded emissions array, which enables gridded emissions from 
  ! an ascii file (emission.txt) used in conjunction with call emsgrd.
  !-------------------------------------------------------------------------------


  !***********************
  !dwen(20090814) HSYPLIT4.5 splitted old EMSDAT into EMSINI and EMSDAT
  CALL EMSINI(NLOC,SPOT(:)%OLAT,SPOT(:)%OLON,QFILE,NQLON,NQLAT,NPVAL)

  IF(QFILE)THEN

     ALLOCATE (polid (MAX(NPVAL,NUMTYP)),STAT=kret)
     IF(kret.ne.0)THEN       
        ECODE='Gridded emission pollutant'
        WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
        STOP 900
     END IF

     ALLOCATE (qarea (nqlon,nqlat,MAX(NPVAL,NUMTYP),24),STAT=kret)
     IF(kret.ne.0)THEN
        ECODE='Gridded emission array'
        WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
        STOP 900
     END IF

     !***************************************************************
     !dwen(20090814) HSYPLIT4.5 splitted old EMSDAT into EMSINI and EMSDAT
     !               HYSPLIT4.5 removed QDLON,QDLAT,NQLON,NQLAT,QAREA,POLID from 
     !                          COMMON bloc and added to argument list
     !               HYSPLIT4.5 moved NLOC to MESINI, and assigned from CONTROL
     !               HYSPLIT4.5 output QFILE flag in EMSINI instead of EMSDAT
     !               HYSPLIT4.5 removed SPOT from COMMON bloc and added SPOT to argument list
     !               HYSPLIT4.5 removed NUMTYP for input dimension consistency check

     CALL EMSDAT(QDLON,QDLAT,NQLON,NQLAT,SPOT(:)%OLON,SPOT(:)%OLAT,QAREA,POLID)

     !     enables gridded emissions from an ascii file (emission.txt)
     !     used in conjunction with call emsgrd
     !      CALL EMSDAT(NLOC,NUMTYP,QFILE)
     !***************************************************************

     !    the emission array must have an element for each pollutant type defined
     !    in the control file even if it is not defined in the emissions file
     NPVAL=MAX(NPVAL,NUMTYP)

  END IF

  !-------------------------------------------------------------------------------
  ! optional initialization from previous simulation dump file
  !-------------------------------------------------------------------------------

  ! initial value of particle/puff counter
  KPM=0

  ! if file exists then read all points
  !**********************************************************
  !dwen(20090306):copied from STILT
  IF(NDUMP.GT.0)THEN
     !**********************************************************
     INQUIRE(FILE=PINPF,EXIST=FTEST1)
     IF(NINIT.EQ.0)FTEST1=.FALSE.
     IF(PINPF(1:1).EQ.' ')FTEST1=.FALSE.

     IF(FTEST1)THEN
        ALLOCATE (tlat(maxpar),tlon(maxpar),STAT=kret)
        IF(kret.ne.0)THEN      
           ECODE='Particle lat/lon'
           WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
           STOP 900
        END IF

        METO%GDISX=CGSIZE*1000.0
        METO%GDISY=CGSIZE*1000.0
        DT=CONC(1)%DELTA%MACC 

        !    for LINUX add: CONVERT='BIG_ENDIAN')
        OPEN(KF23,FILE=PINPF,FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

        !    load particles if time matches
        !****************************
        !dwen(20090814): HYSPLIT4.7 supplemented PARINP
        CALL PARINP(JET,KG,KT,KPM,MASS,TLAT,TLON,XPOS,YPOS,ZPOS,SIGH,SIGU,   &
             SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,ZMDL,MAXPAR,NINIT)


        !    initialize output as snapshot concentration (3/10/2003 - RRD)  
        CALL CONZRO(CONC,NUMGRD,DT,JET,IFHR,CSUM)

        DO J=1,KPM
           HDWPX=HDWP(J)                            ! complex mode
           IF(HDWPX.GE.100)HDWPX=MOD(HDWP(J)/10,10) ! simple mode

           IF(ICHEM.EQ.4)THEN
              !          sum to meteo projection concentration grid
              !dwen(20090814): HYSPLIT4.7 supplemented METSUM
              CALL METSUM(CONC,NUMGRD,XPOS(J),YPOS(J),METO%GDISX,METO%GDISX,      &
                   DT,JET,ZMDL,ZSFC,KMSL,CGSIZE,MASS(:,J),DEPT,ZPOS(J),    &
                   SIGH(J),SIGW(J),HDWPX,PTYP(J),CSUM)
           ELSE
              IF(CMASS.EQ.0)THEN
                 !          sum to concentration array (3/10/2003 - RRD)  
                 !dwen(20090817): no CONSUM call around here in STILT code
                 CALL CONSUM(CONC,NUMGRD,TLAT(J),TLON(J),DT,JET,ZMDL,ZSFC,KMSL,      &
                      CGSIZE,MASS(:,J),DEPT,ZPOS(J),SIGH(J),SIGW(J),          &
                      HDWPX,PTYP(J),CSUM)
              ELSE
                 !          sum mass to grid as mass on standard lat/lon projection
                 !dwen(20090814): HYSPLIT4.9 supplemented MASSUM
                 CALL MASSUM(CONC,NUMGRD,TLAT(J),TLON(J),DT,JET,ZMDL,ZSFC,KMSL,      &
                      CGSIZE,MASS(:,J),ZPOS(J),PTYP(J),CSUM)
              END IF
           END IF
        END DO

        !    output snapshot concentration (3/10/2003 - RRD)  
        !    requires output start set prior to model run start time
        CALL CONDSK(CONC,DIRT,ICHEM,KG,KT,NUMGRD,NUMTYP,DT,JET,IFHR,CPACK,CSUM)
        CALL CONZRO(CONC,NUMGRD,DT,JET,IFHR,CSUM)

        !    delete any particles that may not be on new grid
        !******************************************************************
        !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
        !                           SIGU and SIGV for horizontal velocity sigma 
        !                SIGH:horizontal puff sigma (new verson), horiz sigma (sigma)(old version) 
        !                SIGV:turbulence v'2 (new version), vert sigma (sigma)(old version)

        CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,             &
             PAGE,PTYP,PGRD,NSORT,KHMAX,0.0)

        !            XXX=0.0
        !            CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,HDWP,         &
        !     &           PAGE,PTYP,PGRD,NSORT,KHMAX,XXX)
        !******************************************************************

        !    file still open
        INQUIRE(FILE=PINPF,OPENED=FTEST1)
        IF(NINIT.EQ.1)FTEST1=.FALSE.

        ! end pardump exist test
     END IF
     !**********************************************
     !dwen(20090306):
  end if   ! end ndump>0 
  !**********************************************

  !-------------------------------------------------------------------------------
  ! starting sigma height initialization
  !-------------------------------------------------------------------------------

  ! initial wind speed max (1.2 km/min = 20 m/s) for dt
  UMAX=1.2

  !**********************************************************
  !dwen(20090306): used in STILT
  !     convert to gp/min
  !     ??????vmax
  !dwen(20090824)      VMAX=UMAX/GRID(MAX0(KG,KGRID))%SIZE
  VMAX=UMAX/GRID(KG,kt)%SIZE
  !**********************************************************

  ! redefine subgrid parameters
  IF(KPM.LE.0)THEN
     !    no particles or first time initialize with source
     !**********************************************************
     !dwen(20090815): check the unit of UMAX in STILT and hysplit
     !dwen(20090815): HYSPLIT4.4 added MGMIN to argument list
     !dwen(20090815): HYSPLIT4.6 added NGRD and PGRD to argument list for multiple grids
     CALL ADVRNG(NGRD,MGMIN,UMAX,NLOC,SPOT%XP,SPOT%YP,SPOT%KG)
     !      CALL ADVRNG(MAX0(KG,KGRID),VMAX,NLOC,XARG,YARG)
     !**********************************************************
  ELSE
     !    include source for area determination
     DO N=1,NLOC
        XPOS(KPM+N)=SPOT(N)%XP
        YPOS(KPM+N)=SPOT(N)%YP
        PGRD(KPM+N)=SPOT(N)%KG
     END DO
     !    subsequent time use mean particle position to define grid
     !**********************************************************
     !dwen(20090815): check the unit of UMAX in STILT and hysplit
     !dwen(20090815): HYSPLIT4.4 added MGMIN to argument list
     !dwen(20090815): HYSPLIT4.6 added NGRD and PGRD to argument list for multiple grids
     CALL ADVRNG(NGRD,MGMIN,UMAX,(KPM+NLOC),XPOS,YPOS,PGRD)
     !      CALL ADVRNG(MAX0(KG,KGRID),VMAX,NLOC,XARG,YARG)
     !**********************************************************
  END IF


  DT=1.0

  !******************************************
  !dwen(20090306): dt = 10 used in STILT
  !      DT=10.0
  !******************************************

  IF(BACK)DT=-DT
  nloop : DO N=1,NLOC

     XTMP =SPOT(N)%XP
     YTMP =SPOT(N)%YP
     ZTMP =SPOT(N)%ZP
     KGRID=SPOT(N)%KG
     KTIME=KT

     !dummy call to advpnt to return terrain height

     !dwen(20090813)     CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,     &
     !dwen(20090813)                 INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,     &
     !dwen(20090813)                 ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET)

     !******************************************
     ! JCL:   add DUDZ&DVDZ&WWOUT as output arguments in dummy call as well
     ! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
     ! CHG:(12/05/01) add 'CONVDUR' (here only dummy, use first particles CONVDUR)
     ! CHG:(9/17/02) add 'ICONVECT' as convection flag
     ! CHG(09/18/03) pass on RAMSFLG
     !        dummy call to advpnt to return terrain height
     !dwen(20090310): ADVPNT needs more check 
     !dwen(20090319): HYSPLIT4.5 splitted METO into METZ and METO 
     !                METZ defines profile advection variables, and METO defines 
     !                advection surface variables.
     !dwen(20090319): HYSPLIT4.9 added TRATIO (advection stability ratio) in argument list  
     !dwen(20090319): HYSPLIT4.6 added KTIME and made 
     !                 simultaneous multiple meteorology available  
     !dwen(20090319): NTIM (number of meteo data times defined in CONTROL)
     !                  has been added and initialized in new version
     !dwen(20090319): turbelent kinetic energy (TKERD and TKERN ) has been added in new version
     !dwen(20090813): remove DUDZ and DVDZ. no use in ADVPNT 
     !dwen(20090813): remove ISOT. HYSPLIT removed it in version 4.9
     !dwen(20090813): remove is_off_grid from argument list
     !dwen(20090822): add CONC in argument list
     !                 
     CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,     &
          INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,     &
          ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET,                &
          conc,wwout,zicontroltf,nhrszi,zipresc,convdur(1),iconvect,    &
          ramsflg,ecmflg,xposprev,yposprev)

     !******************************************

     IF(KRET.NE.0) CYCLE nloop

     !    initial grid may not be valid, hence particles are remapped
     IF(KGRID.NE.SPOT(N)%KG)THEN  
        SPOT(N)%KG=KGRID
        IF(GRID(KGRID,KTIME)%LATLON)THEN
           CALL GBL2XY(KGRID,KTIME,SPOT(N)%OLAT,SPOT(N)%OLON,                  &
                SPOT(N)%XP,SPOT(N)%YP)
        ELSE
           CALL CLL2XY_wps(GRID(KGRID,KTIME)%GBASE,SPOT(N)%OLAT,SPOT(N)%OLON,      &
                SPOT(N)%XP,SPOT(N)%YP,GRID(KGRID,KTIME)%proj)
        END IF
     END IF

     !    if input specified MSL rather than AGL, then convert to AGL
     IF(KMSL.EQ.1)SPOT(N)%OLVL=MAX(0.0,SPOT(N)%OLVL-METO%ZTER)
     !    terrain adjusted sigma
     SPOT(N)%ZP=1.0-SPOT(N)%OLVL/(ZMDL-METO%ZTER)

     IF(KMSL.EQ.1)THEN
        !       another dummy call to get correct initial pressure
        XTMP=SPOT(N)%XP
        YTMP=SPOT(N)%YP
        ZTMP=SPOT(N)%ZP


        !dwen(20090814)        CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,  &
        !dwen(20090814)                    INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,  &
        !dwen(20090814)                    ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET)

        !******************************************
        ! JCL:   add DUDZ&DVDZ&WWOUT as output arguments in dummy call as well
        ! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
        ! CHG:(12/05/01) add 'CONVDUR' (here only dummy, use first particles CONVDUR)
        ! CHG:(9/17/02) add 'ICONVECT' as convection flag
        ! CHG(09/18/03) pass on RAMSFLG
        !        dummy call to advpnt to return terrain height
        !dwen(20090310): ADVPNT needs more check 
        !dwen(20090319): HYSPLIT4.5 splitted METO into METZ and METO 
        !                METZ defines profile advection variables, and METO defines 
        !                advection surface variables.
        !dwen(20090319): HYSPLIT4.9 added TRATIO (advection stability ratio) in argument list  
        !dwen(20090319): HYSPLIT4.6 added KTIME and made 
        !                 simultaneous multiple meteorology available  
        !dwen(20090319): NTIM (number of meteo data times defined in CONTROL)
        !                  has been added and initialized in new version
        !dwen(20090319): turbelent kinetic energy (TKERD and TKERN ) has been added in new version
        !dwen(20090813): remove DUDZ and DVDZ. no use in ADVPNT 
        !dwen(20090813): remove ISOT. HYSPLIT removed it in version 4.9
        !dwen(20090813): remove is_off_grid from argument list
        !                 
        CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,     &
             INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,     &
             ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET,                &
             conc,wwout,zicontroltf,nhrszi,zipresc,convdur(1),iconvect,    &
             ramsflg,ecmflg,xposprev,yposprev)

        !******************************************

        IF(KRET.NE.0) CYCLE nloop
     END IF
  END DO nloop


  !-------------------------------------------------------------------------------
  ! lagrangian sampling configuration
  !-------------------------------------------------------------------------------

  NUMLAG=0
  INQUIRE(FILE='LAGSET.CFG',EXIST=LAGSAM)
  IF(LAGSAM)THEN
     OPEN(KF50,FILE='LAGSET.CFG')
     READ(KF50,*)NUMLAG

     ALLOCATE (lags(numlag), STAT=kret) 
     IF(kret.ne.0)THEN     
        ECODE='Lagrangian sampling initialization'
        WRITE(*,*)'ERROR hymodelc: memory allocation - ',ECODE,KRET
        STOP 900
     END IF

     DO N=1,NUMLAG
        READ(KF50,*)LAGS(N)%RLAT,LAGS(N)%RLON,LAGS(N)%RLVL
        READ(KF50,*)LAGS(N)%WDIR,LAGS(N)%WSPD  
        READ(KF50,*)LAGS(N)%RYR, LAGS(N)%RMO, LAGS(N)%RDA,LAGS(N)%RHR,LAGS(N)%RMN
        READ(KF50,*)LAGS(N)%SYR, LAGS(N)%SMO, LAGS(N)%SDA,LAGS(N)%SHR,LAGS(N)%SMN
        READ(KF50,*)LAGS(N)%AVRG
        READ(KF50,*)LAGS(N)%DISK
        READ(KF50,*)LAGS(N)%FILE

        !    convert date/time to accumulated minutes
        CALL TM2MIN(LAGS(N)%RYR,LAGS(N)%RMO,LAGS(N)%RDA,                          &
             LAGS(N)%RHR,LAGS(N)%RMN,LAGS(N)%RACM)
        CALL TM2MIN(LAGS(N)%SYR,LAGS(N)%SMO,LAGS(N)%SDA,                          &
             LAGS(N)%SHR,LAGS(N)%SMN,LAGS(N)%SACM)

        !    open the output files
        IF(KF51+N-1.GT.KF52)THEN
           WRITE(KF21,*)'WARNING main: exceed max files for lagrangian samplers'  
           WRITE(KF21,*)'NO DISK OUTPUT for sampler numbers ',N,' and greater'   
           LAGS(N)%DISK=999999
           LAGS(N)%SACM=0 
        ELSE
           OPEN(KF51+N-1,FILE=LAGS(N)%FILE) 
           !       first output time always one interval greater
           LAGS(N)%SACM=LAGS(N)%SACM+LAGS(N)%DISK
        END IF

        !    convert starting position to meteorology grid units
        IF(GRID(KG,KT)%LATLON)THEN
           CALL GBL2XY(KG,KT,LAGS(N)%RLAT,LAGS(N)%RLON,                     &
                LAGS(N)%XP,LAGS(N)%YP)
        ELSE
           CALL CLL2XY_wps(GRID(KG,KT)%GBASE,LAGS(N)%RLAT,LAGS(N)%RLON,         &
                LAGS(N)%XP,LAGS(N)%YP,GRID(KG,KT)%proj)
        END IF

        !    test for position on grid
        LAGS(N)%KG=KG  ! initial meteorological grid
        IF(LAGS(N)%XP.LT.1.0.OR.LAGS(N)%XP.GT.FLOAT(GRID(KG,KT)%NX).OR.     &
             LAGS(N)%YP.LT.1.0.OR.LAGS(N)%YP.GT.FLOAT(GRID(KG,KT)%NY))THEN
           IF(.NOT.GRID(KG,KT)%GBLDAT)THEN
              WRITE(KF21,*)'WARNING main: lagrangian release point off grid'
              WRITE(KF21,*)'Position:',LAGS(N)%RLAT,LAGS(N)%RLON
              WRITE(KF21,*)'Grid Loc:',LAGS(N)%XP,  LAGS(N)%YP
              LAGS(N)%RACM=0  ! turn off release for this location
           END IF
        END IF

        XTMP=LAGS(N)%XP
        YTMP=LAGS(N)%YP
        ZTMP=1.0 
        KGRID=KG
        KTIME=KT

        !    dummy call to advpnt to return terrain height

        !dwen(20090814)     CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,     &
        !dwen(20090814)                 INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,     &
        !dwen(20090814)                 ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET)

        !******************************************
        ! JCL:   add DUDZ&DVDZ&WWOUT as output arguments in dummy call as well
        ! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
        ! CHG:(12/05/01) add 'CONVDUR' (here only dummy, use first particles CONVDUR)
        ! CHG:(9/17/02) add 'ICONVECT' as convection flag
        ! CHG(09/18/03) pass on RAMSFLG
        !        dummy call to advpnt to return terrain height
        !dwen(20090310): ADVPNT needs more check 
        !dwen(20090319): HYSPLIT4.5 splitted METO into METZ and METO 
        !                METZ defines profile advection variables, and METO defines 
        !                advection surface variables.
        !dwen(20090319): HYSPLIT4.9 added TRATIO (advection stability ratio) in argument list  
        !dwen(20090319): HYSPLIT4.6 added KTIME and made 
        !                 simultaneous multiple meteorology available  
        !dwen(20090319): NTIM (number of meteo data times defined in CONTROL)
        !                  has been added and initialized in new version
        !dwen(20090319): turbelent kinetic energy (TKERD and TKERN ) has been added in new version
        !dwen(20090813): remove DUDZ and DVDZ. no use in ADVPNT 
        !dwen(20090813): remove ISOT. HYSPLIT removed it in version 4.9
        !dwen(20090813): remove is_off_grid from argument list
        !                 
        CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,     &
             INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,     &
             ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET,                &
             conc,wwout,zicontroltf,nhrszi,zipresc,convdur(1),iconvect,    &
             ramsflg,ecmflg,xposprev,yposprev)

        !*************************************************

        IF(KRET.NE.0) METO%ZTER=0.0 
        LAGS(N)%KG=KGRID

        !    if input specified MSL rather than AGL, then convert to AGL
        IF(KMSL.EQ.1)LAGS(N)%RLVL=MAX(0.0,LAGS(N)%RLVL-METO%ZTER)

        !    terrain adjusted sigma
        LAGS(N)%ZP=1.0-LAGS(N)%RLVL/(ZMDL-METO%ZTER)

        !    test for optional vector to force the trajectory 
        IF(LAGS(N)%WDIR.NE.0.0.OR.LAGS(N)%WSPD.NE.0.0)THEN
           !        convert to vector m/s to grid-pts/min
           XTMP=60.0*LAGS(N)%WSPD*SIN(LAGS(N)%WDIR/57.29578)/METO%GDISX
           YTMP=60.0*LAGS(N)%WSPD*COS(LAGS(N)%WDIR/57.29578)/METO%GDISY
           !        rotate from true to grid coordinates
           IF(GRID(KGRID,KTIME)%LATLON)THEN
              METO%UVEL=XTMP
              METO%VVEL=YTMP
           ELSE
              !dwen(20090815):HYSPLIT4.7 supplemented CC2GXY
              !            CALL CC2GXY(GRID(KGRID,KTIME)%GBASE,LAGS(N)%XP,LAGS(N)%YP,         &
              !                        XTMP,YTMP,METO%UVEL,METO%VVEL,GRID(KGRID,KTIME)%proj)
              CALL CC2GXY(GRID(KGRID,KTIME)%GBASE,LAGS(N)%XP,LAGS(N)%YP,         &
                   XTMP,YTMP,METO%UVEL,METO%VVEL)
           END IF
        END IF
     END DO

     WRITE(KF21,*)' NOTICE   main: configured for lagrangian sampling - ',NUMLAG
  END IF
  CLOSE (KF50)

  !-------------------------------------------------------------------------------
  ! final check on consistency of certain options
  !-------------------------------------------------------------------------------

  ! mass removal fraction set to zero for particle simulations
  IF(INITD.EQ.0)FRMR=0.0

  IF(ICHEM.EQ.4)THEN
     IF(GRID(KG,KT)%LXR.NE.GRID(KG,KT)%NX.OR.                                  &
          GRID(KG,KT)%LYR.NE.GRID(KG,KT)%NY) THEN
        WRITE(*,*)'Configuration error, see MESSAGE file for more information'
        WRITE(KF21,*)'*ERROR* main: ICHEM=4 requires full grid meteo data'
        WRITE(KF21,*)' Set MGMIN in SETUP.CFG namelist to a value > than: ',   &
             MAX(GRID(KG,KT)%NX,GRID(KG,KT)%NY)
        STOP 900
     END IF
  END IF

  !-------------------------------------------------------------------------------
  ! optional code section to setup the Global Eulerian Model (GEM)
  !-------------------------------------------------------------------------------

  !****************************************************
  !dwen(20090306):moved from STILT
  ! JCL:initialize counter: # of [min] since last time output written to PARTICLE.DAT
  COUNTOUT=0.0
  !****************************************************

  !-------------------------------------------------------------------------------
  ! main loop over number of simulation hours
  !-------------------------------------------------------------------------------


  KGRID=KG
  KTIME=KT

  DO KH=1,NHRS

     !    time step minutes function of gridsize and max wind speed
     !    such that (uvel dt <= TRATIO * delx) and dt(max)=60, where uvel
     !    is defined as grid points per minute and delx is the grid
     !    spacing in km

     IF(UMAX.LE.0.0)UMAX=1.2
     MAXDT=MIN(60,MINVAL(CONC%DELTA%MACC),MINVAL(DREC(1:KGRID,KTIME)%DELTA))
     DT=MAX(1, MIN(MAXDT, NINT(TRATIO*CGSIZE/UMAX)))
     DO WHILE (MOD(MAXDT,INT(DT)).NE.0.AND.INT(DT).GT.1)
        DT=DT-1.0
     END DO

     !    namelist over-ride for fixed time step when delt<>0
     IF(DELT.NE.0.0)DT=ABS(DELT)

     !    when nothing to advect set to maximum value
     !### IF(KH.GT.1.AND.KPM.EQ.0)DT=60.0

     !*****************************************************
     !dwen(20090306):from STILT
     ! CHG(09/23/03) don't use fixed time step for RAMS
     IF(RAMSFLG.AND.DELT.GT.0.0)THEN
        WRITE(*,*)"ERROR: don't use fixed time step for RAMS"
        STOP
     END IF
     !*****************************************************

     !    set integration direction dependent time step 
     IF(BACK)DT=-DT

     !    redefine subgrid parameters
     IF(KPM.LE.0)THEN
        !       no particels or first time initialize with source

        !**********************************************************
        !dwen(20090815): check the unit of UMAX in STILT and hysplit
        !dwen(20090815): HYSPLIT4.4 added MGMIN to argument list
        !dwen(20090815): HYSPLIT4.6 added NGRD and PGRD to argument list for multiple grids
        CALL ADVRNG(NGRD,MGMIN,UMAX,NLOC,SPOT%XP,SPOT%YP,SPOT%KG)
        !          CALL ADVRNG(MAX0(KG,KGRID),UMAX,NLOC,XARG,YARG)
        !*******************************************************
     ELSE
        !       subsequent time use mean particle position to define grid
        !**********************************************************
        !dwen(20090815): check the unit of UMAX in STILT and hysplit
        !dwen(20090815): HYSPLIT4.4 added MGMIN to argument list
        !dwen(20090815): HYSPLIT4.6 added NGRD and PGRD to argument list for multiple grids
        CALL ADVRNG(NGRD,MGMIN,UMAX,KPM,XPOS,YPOS,PGRD)
        !         CALL ADVRNG(MAX0(KG,KGRID),UMAX,KPM,XPOS,YPOS)
        !*******************************************************
     END IF

     !    diagnostic model output
     IF(KH.EQ.1)WRITE(KF21,*)' NOTICE   main: Initial time step (min) ',INT(DT)

     !    reset time step variables
     UMAX=0.0
     CGSIZE=GRID(KGRID,KTIME)%SIZE

     !    temporal emissions variations defined by an input file    
     !dwen(20090815): HYSPLIT4.9 supplemented EMSMAT
     CALL EMSMAT(SPRT,NGRD,JET,KGRID,KTIME,NLOC,NUMTYP,EFILE,QTEMP,BACK)

     !-------------------------------------------------------------------------------
     !    sub-loop for number of time steps per hour
     !-------------------------------------------------------------------------------



     NSTEP=ABS(60/INT(DT))
     DO KS=1,NSTEP

        !*************************************************
        !dwen(20090306):the following lines used in STILT
        RAMSFLG=GRID(kgrid,ktime)%MODEL_ID.EQ.'RAMS'
        ECMFLG =GRID(kgrid,ktime)%MODEL_ID(1:2) == 'EC'
        AWRFFLG=GRID(kgrid,ktime)%MODEL_ID(2:4).EQ.'WRF'

        !!        update grid on which new particles start if required
        !         IF(KGRID.GT.KG)THEN
        !            KG=KGRID
        !! CHG(09/22/03) update RAMSFLG
        !            RAMSFLG=GRID(KG)%MODEL_ID.EQ.'RAMS'
        !            ECMFLG =GRID(KG)%MODEL_ID(1:2) == 'EC'
        !            AWRFFLG=GRID(KG)%MODEL_ID(2:4).EQ.'WRF'
        !
        !!           remap starting positions to new grid
        !            DO N=1,NLOC
        !
        !! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
        !              IF(GRID(KG)%LATLON)THEN
        !!dwen(20090315): add KT
        !!                 CALL GBL2XY(KG,SPOT(N)%OLAT,SPOT(N)%OLON,              &
        !!     &                          SPOT(N)%XP, SPOT(N)%YP)
        !                 CALL GBL2XY(KG,kt,SPOT(N)%OLAT,SPOT(N)%OLON,              &
        !     &                          SPOT(N)%XP, SPOT(N)%YP)
        !              ELSE
        !!dwen(20090315): remove GRID(KG)%proj,add KT
        !!                 CALL CLL2XY(GRID(KG)%GBASE,SPOT(N)%OLAT,SPOT(N)%OLON,  &
        !!     &                    SPOT(N)%XP,SPOT(N)%YP, GRID(KG)%proj)
        !                 CALL CLL2XY(GRID(KG,kt)%GBASE,SPOT(N)%OLAT,SPOT(N)%OLON,  &
        !     &                    SPOT(N)%XP,SPOT(N)%YP)
        !                     !IF(GRID(KG)%LATLON)THEN
        !              END IF
        !
        !                     !DO N=1,NLOC
        !            END DO
        !         END IF
        !*************************************************

        !       emissions each time step are from all grid cells or from a point
        IF(QFILE)THEN
           !          start any new particles from grid cells (from emsdat)

           !**************************************************
           ! JCL:  add BACK as argument
           !dwen(20090306):use updated EMSGRD (updated from HYSPLIT)
           !dwen(20090815): HYSPLIT4.5 removed TBAR. It is no use in EMSGRD
           !dwen(20090815): HYSPLIT4.5 removed SPOT and DIRT from COMMON bloc and added to argument list
           !dwen(20090815): HYSPLIT4.6 removed KG from argument list, assigned it from SPOT
           !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
           !                           SIGU and SIGV for horizontal velocity sigma 
           !dwen(20090815): HYSPLIT4.5 added NPVAL,NQLON,NQLAT,QDLON,QDLAT,QAREA,POLID,MAXPAR to argument list
           !dwen(20090815): HYSPLIT uses DT to determine BACK flag 

           CALL EMSGRD(SPOT,DIRT,NUMTYP,KPM,INITD,DT,JET,NSORT,MASS,           &
                XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD, &
                QCYCLE,NUMPAR,NPVAL,NQLON,NQLAT,QDLON,QDLAT,            &
                QAREA,POLID,MAXPAR)

           !            CALL EMSGRD(NUMTYP,KPM,INITD,DT,JET,KG,                     &
           !     &         NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,      &
           !     &         PTYP,PGRD,REAL(QCYCLE,dp),NUMPAR,TBAR,BACK)
           !**************************************************
        ELSEIF(QTEMP)THEN
           !          temporal emission variations from input file (emsmat)
           !          file emission option not valid with ensemble or backward mode
           KPT=KPM
           !dwen(20090815): HYSPLIT4.9 supplemented EMSTMP
           CALL EMSTMP(SPRT,KGRID,NLOC,NUMTYP,KPM,INITD,DT,JET,NUMPAR,MAXPAR,  &
                NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,     &
                PAGE,PTYP,PGRD,job_id,num_job,KEMIT)

           IF(KPM.GT.KPT)THEN
              CNTR=0.0     
              AVGMIXD=0.0                 
              AVGRISE=0.0

              !             new points have been added when KPM has increased
              DO KP=KPT+1,KPM
                 METO%ZTER=0.0
                 XTMP=XPOS(KP)
                 YTMP=YPOS(KP)  
                 ZTMP=1.0-ZPOS(KP)/(ZMDL-METO%ZTER)

                 !                dummy call to get terrain and other meteo variables
                 !                rdep is set to true to return flux variables

                 !dwen(20090814)                 CALL ADVPNT(                                                  &
                 !dwen(20090814)                 METZ,METO,BACK,.TRUE.,CDEP,.TRUE.,.FALSE.,KSFC,TKERD,TKERN,   &
                 !dwen(20090814)                 INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,     &
                 !dwen(20090814)                 ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET)

                 !******************************************
                 ! JCL:   add DUDZ&DVDZ&WWOUT as output arguments in dummy call as well
                 ! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
                 ! CHG:(12/05/01) add 'CONVDUR' (here only dummy, use first particles CONVDUR)
                 ! CHG:(9/17/02) add 'ICONVECT' as convection flag
                 ! CHG(09/18/03) pass on RAMSFLG
                 !        dummy call to advpnt to return terrain height
                 !dwen(20090310): ADVPNT needs more check 
                 !dwen(20090319): HYSPLIT4.5 splitted METO into METZ and METO 
                 !                METZ defines profile advection variables, and METO defines 
                 !                advection surface variables.
                 !dwen(20090319): HYSPLIT4.9 added TRATIO (advection stability ratio) in argument list  
                 !dwen(20090319): HYSPLIT4.6 added KTIME and made 
                 !                 simultaneous multiple meteorology available  
                 !dwen(20090319): NTIM (number of meteo data times defined in CONTROL)
                 !                  has been added and initialized in new version
                 !dwen(20090319): turbelent kinetic energy (TKERD and TKERN ) has been added in new version
                 !dwen(20090813): remove DUDZ and DVDZ. no use in ADVPNT 
                 !dwen(20090813): remove ISOT. HYSPLIT removed it in version 4.9
                 !dwen(20090813): remove is_off_grid from argument list
                 !                 
                 CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,     &
                      INITD,XTMP,YTMP,ZTMP,JET,DT,TRATIO,KGRID,KTIME,NGRD,NTIM,     &
                      ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET,                &
                      conc,wwout,zicontroltf,nhrszi,zipresc,convdur(1),iconvect,    &
                      ramsflg,ecmflg,xposprev,yposprev)

                 !**********************************************

                 IF(KRET.EQ.0)THEN 
                    !                   meteorological variables are available
                    IF(SIGW(KP).GT.0.0)THEN
                       !                      non-zero value is btu/hr heat release for plume rise
                       !dwen(20090815): HYSPLIT4.9 supplemented EMRISE
                       CALL EMRISE(PAGE(KP),SIGW(KP),METO%SSP,METO%UBAR,       &
                            METO%USTR,METO%MIXD,HEIGHT)
                       ZPOS(KP)=MAX(0.0,1.0-HEIGHT/(ZMDL-METO%ZTER))
                       SIGW(KP)=0.0
                       PAGE(KP)=0

                       CNTR=CNTR+1.0
                       AVGMIXD=AVGMIXD+METO%MIXD
                       AVGRISE=AVGRISE+HEIGHT

                    ELSE
                       !                      normally all turbulence values set to zero at start
                       ZPOS(KP)=1.0-ZPOS(KP)/(ZMDL-METO%ZTER)
                    END IF
                 ELSE
                    ZPOS(KP)=ZTMP 
                 END IF
              END DO
              IF(CNTR.GT.0.0.AND.KS.EQ.1)                                      &
                   WRITE(KF21,*)' NOTICE emrise: (mixd,rise) - ',                   &
                   NINT(AVGMIXD/CNTR),NINT(AVGRISE/CNTR)
           END IF

        ELSE
           !          start any new particles from a point (from emsset)
           !**************************************************
           ! JCL:  add BACK as argument
           !dwen(20090306): HYSPLIT4.5 removed SPOT,DIRT and CONC from COMMON bloc and added to argument list
           !dwen(20090815): HYSPLIT4.6 removed KG from argument list, assigned it from SPOT
           !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
           !                           SIGU and SIGV for horizontal velocity sigma 
           !dwen(20090816): HYSPLIT4.6 added job_id and num_job for mpi implementation
           !dwen(20090815): HYSPLIT uses DT to determine BACK flag 

           CALL EMSPNT(SPOT,DIRT,NLOC,NUMTYP,KPM,INITD,DT,JET,NSORT,           &
                MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP, &
                PGRD,QCYCLE,NUMPAR,MAXPAR,job_id,num_job,ichem,KEMIT)

           !            CALL EMSPNT(NLOC,NUMTYP,KPM,INITD,DT,JET,KG,                &
           !     &         NSORT,MASS,XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,      &
           !     &         PTYP,PGRD,REAL(QCYCLE,dp),NUMPAR,BACK)
           !*******************************************************
        END IF

        !       diagnostic variables initialization
        TMASS=0.0
        DO KZ=1,NLVL   
           ZMASS(KZ)=0.0
        END DO

        !       emission of lagrangian samplers
        IF(LAGSAM)                                                            &
             !dwen(200900816): HYSPLIT4.9 supplemented LAGEMS
             CALL LAGEMS(NUMLAG,LAGS,MAXPAR,JET,KPM,NSORT,MASS,XPOS,YPOS,ZPOS,  &
             SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD)

        !*********************************************************************
        !dwen(20090309): add the following call according to STILT
        !dwen(20090816): HYSPLIT4.6 removed IBYR from argument list and determined it with tm2day

        IF(RDEP)                                                       &
             CALL SUNANG(JET,METO%PLAT,METO%PLON,EA,SEA)      ! SUNANG was called this way in HYSPLIT4.9
        !     &      CALL SUNANG                                                 &
        !     &      (SPOT(1)%IBYR,JET,SPOT(1)%OLAT,SPOT(1)%OLON,EA,SEA)  ! called like this in STILT
        !*********************************************************************

        !       particle resuspension 

        !***************************************************************
        !dwen(20090306): HYSPLIT4.5 removed DIRT and CONC from COMMON bloc and added to argument list
        !dwen(20090815): HYSPLIT4.6 removed KG from argument list, assigned it from SPOT
        !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
        !                           SIGU and SIGV for horizontal velocity sigma 
        !dwen(20090319): HYSPLIT4.6 added KTIME and made 
        !                 simultaneous multiple meteorology available  

        IF(SDEP)                                                               &
             CALL DEPSUS(CONC,DIRT,INITD,KGRID,KTIME,NUMGRD,NUMTYP,DT,ICHEM,     &
             KPM,CSUM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,  &
             PAGE,PTYP,PGRD,NSORT,MAXPAR)

        !     &      CALL DEPSUS(INITD,KG,NUMGRD,NUMTYP,DT,KPM,CSUM,MASS,        &
        !     &      XPOS,YPOS,ZPOS,SIGH,SIGV,SIGX,HDWP,PAGE,PTYP,PGRD,NSORT)
        !***************************************************************


        !***************************************************************
        !dwen(20090306):following lines from STILT
        ! JCL:   increment # of [min] since last time output written to PARTICLE.DAT
        COUNTOUT=COUNTOUT+ABS(INT(DT))

        ! JCL:   write out results to PARTICLE.DAT only in specified intervals
        !             and OUTDT=0 means that results would be output EACH timestep

        WRITEOUT=(OUTDT.EQ.0.0).OR.(COUNTOUT.GE.OUTDT)
        ! JCL:   ensures results always written for LAST timestep
        WRITEOUT=WRITEOUT.OR.(KH.EQ.NHRS.AND.KS.EQ.NSTEP)

        ! JCL:(4/27/01)initialize counter for number of particles that left grid area
        COUNTNPAROUT=0

        ! CHG:(9/17/02) initialize temporary convtmp variable to 0
        CONVTMP=0
        !-----------------------------------------------------------
        ! JCL:(11/03/03) store profile of wind errors
        IF(KH.EQ.1.AND.KS.EQ.1.AND.(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3))THEN
           WRITE(*,*)'ENTERED INITIALIZATION!!!!'

           !           read in statistics of wind errors
           INQUIRE(FILE='WINDERR',EXIST=FTEST)
           IF(FTEST)THEN
              OPEN(KF26,FILE='WINDERR')
              !stdev of U&V errs [m/s]
              READ(KF26,*)SIGUVERR
              !correlation Lagr timescale of U&V errs [min]
              READ(KF26,*)TLUVERR
              !correlation lengthscale of U&V errs in vertical dir [m]
              READ(KF26,*)ZCORLEN
              !correlation lengthscale of U&V errs in hor. dir [km]
              READ(KF26,*)HORCORLEN
              CLOSE(KF26)
           ELSE
              WRITE(*,*)'============File "WINDERR" not found=========='
              STOP
              !IF(FTEST)THEN
           END IF

           !           initialize error profiles with appropriate vertical correlation
           DO KPP=1,KPM
              !dwen(20090315):
              UVERR(KPP)%U(1)=GASDEV(RSEED,1.0)
              UVERR(KPP)%V(1)=GASDEV(RSEED,1.0)
              !               UVERR(KPP)%U(1)=GASDEV(RSEED,DBLE(1.0))
              !               UVERR(KPP)%V(1)=GASDEV(RSEED,DBLE(1.0))
              !dwen(20090315):substitute NZM with NLVL
              !            DO KZZ=2,NZM
              DO KZZ=2,nlvl
                 DELTAZ=(ZMDL-METO%ZTER)*(ZSG(KZZ-1)-ZSG(KZZ))
                 RAUTO=EXP(-1.0*DELTAZ/ZCORLEN)
                 !dwen(20090315):
                 UU=GASDEV(RSEED,1.0)
                 !               UU=GASDEV(RSEED,DBLE(1.0))
                 UVERR(KPP)%U(KZZ)=RAUTO*UVERR(KPP)%U(KZZ-1)+             &
                      &                           SQRT(1.0-RAUTO*RAUTO)*UU
                 !dwen(20090315):
                 VV=GASDEV(RSEED,1.0)
                 !               VV=GASDEV(RSEED,DBLE(1.0))
                 UVERR(KPP)%V(KZZ)=RAUTO*UVERR(KPP)%V(KZZ-1)+             &
                      &                           SQRT(1.0-RAUTO*RAUTO)*VV
                 !               WRITE(*,*)"initial:",KZZ,UVERR(KPP)%U(KZZ)
                 !DO KZZ=2,NZM
              END DO
              !DO KPP=1,KPM
           END DO

        END IF
        !-----------------------------------------------------------
        ! CHG:(27/04/06) initialize zi errors
        IF(KH.EQ.1.AND.KS.EQ.1.AND.(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3))THEN
           WRITE(*,*)'ENTERED ZI ERROR INITIALIZATION!!!!'

           !           read in statistics of ZI errors
           INQUIRE(FILE='ZIERR',EXIST=FTEST)
           IF(FTEST)THEN
              OPEN(KF26,FILE='ZIERR')
              !stdev of ZI errs [%]
              READ(KF26,*)SIGZIERR
              !correlation Lagr timescale of U&V errs [min]
              READ(KF26,*)TLZIERR
              !correlation lengthscale of U&V errs in hor. dir [km]
              READ(KF26,*)HORCORZILEN
              CLOSE(KF26)
              WRITE(*,*)'SIGZIERR:   ',SIGZIERR
              WRITE(*,*)'TLZIERR:    ',TLZIERR
              WRITE(*,*)'HORCORZILEN:',HORCORZILEN
           ELSE
              WRITE(*,*)'============File "ZIERR" not found=========='
              STOP
              !IF(FTEST)THEN
           END IF
           DO KPP=1,KPM
              ZIERR(KPP)=GASDEV(RSEED,SIGZIERR)
           END DO
           !IF(KH.EQ.1.AND.KS.EQ.1.AND.(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3))THEN
        END IF
        !
        !****************************************************************

        !-------------------------------------------------------------------------------
        !       loop through all particles/puffs
        !-------------------------------------------------------------------------------

        ploop : DO KPT=1,KPM

           !       current computation grid index number
           KP=NSORT(KPT)

           !       skip terminated particles
           IF(PGRD(KP).EQ.0) CYCLE ploop

           !dwen(20090827)*******************
           !               from original STILT codes
           XPOSPREV=XPOS(KP)
           YPOSPREV=YPOS(KP)
           !*********************************

           !       set current horizontal distribtion mode
           HDWPX=HDWP(KP)                            ! array in complex mode
           IF(HDWPX.GE.100)HDWPX=MOD(HDWP(KP)/10,10) ! final value in simple mode

           !       advects a single point for one time step

           !dwen(20090814)        CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,  &
           !dwen(20090814)             HDWPX,XPOS(KP),YPOS(KP),ZPOS(KP),JET,DT,TRATIO,PGRD(KP),          &
           !dwen(20090814)             KTIME,NGRD,NTIM,ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET)

           !******************************************
           ! JCL:   add DUDZ&DVDZ&WWOUT as output arguments in dummy call as well
           ! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
           ! CHG:(12/05/01) add 'CONVDUR' (here only dummy, use first particles CONVDUR)
           ! CHG:(9/17/02) add 'ICONVECT' as convection flag
           ! CHG(09/18/03) pass on RAMSFLG
           !        dummy call to advpnt to return terrain height
           !dwen(20090310): ADVPNT needs more check 
           !dwen(20090319): HYSPLIT4.5 splitted METO into METZ and METO 
           !                METZ defines profile advection variables, and METO defines 
           !                advection surface variables.
           !dwen(20090319): HYSPLIT4.9 added TRATIO (advection stability ratio) in argument list  
           !dwen(20090319): HYSPLIT4.6 added KTIME and made 
           !                 simultaneous multiple meteorology available  
           !dwen(20090319): NTIM (number of meteo data times defined in CONTROL)
           !                  has been added and initialized in new version
           !dwen(20090319): turbelent kinetic energy (TKERD and TKERN ) has been added in new version
           !dwen(20090813): remove DUDZ and DVDZ. no use in ADVPNT 
           !dwen(20090813): remove ISOT. HYSPLIT removed it in version 4.9
           !dwen(20090813): remove is_off_grid from argument list
           !                 

           CALL ADVPNT(METZ,METO,BACK,.TRUE.,CDEP,RDEP,.FALSE.,KSFC,TKERD,TKERN,  &
                HDWPX,XPOS(KP),YPOS(KP),ZPOS(KP),JET,DT,TRATIO,PGRD(KP),          &
                KTIME,NGRD,NTIM,ZSG,NLVL,ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET,    &
                conc,wwout,zicontroltf,nhrszi,zipresc,convdur(kp),iconvect,       &
                ramsflg,ecmflg,xposprev,yposprev)

           !****************************************************************
           IF(KRET.NE.0) CYCLE ploop

           !****************************************************************
           !dwen(20090306):following lines used in STILT
           !CCCCCCCCCCCCCCCCCCCC  CONVECTION USING RAMS FLUXES CCCCCCCCCCCCCCCCCCCC
           ! CHG(09/23/03) call CGRELL for convective redistribution w/ RAMS fluxes
           IF ((RAMSFLG .OR. ECMFLG .or. awrfflg) .AND. NTURB .EQ. 0 .AND. ICONVECTRAMS .EQ. 1)THEN
              !convert to normalized height
              Z1=(1.0-ZPOS(KP))*ZMDL
              if (awrfflg) then
                 ! this uses the same logic as the dzsig computation in prfcom:
                 zsigwdn=1.
                 zsigw=2.*zsg(1)-zsigwdn
                 zprofm(1) = (1.0-zsigw)*zmdl
                 zsigwdn=zsigw
                 do kk = 2,nlvl
                    zsigw=2.*zsg(kk)-zsigwdn
                    zprofm(kk) = (1.-zsigw)*zmdl
                    zsigwdn=zsigw
                 enddo
              elseif (ramsflg) then
                 !dwen(20090824)                 zprofm(1:nlvl) = DREC(KG)%HEIGHT(2:(NLVL+1))
                 zprofm(1:nlvl) = DREC(pgrd(kp),ktime)%HEIGHT(2:(NLVL+1))
              ELSE IF (ECMFLG) THEN                         ! arithm. mean of mid-level z values
                 FORALL (kk=1:nlvl-1) zprofm(kk) = 0.5d0*zmdl*(2d0-zsg(kk)-zsg(kk+1))
                 zprofm(nlvl) = zmdl
              else
                 stop 'hymodelc cgrell if-struct: internal logic error'
              endif
              ! CHG(09/24/03) relative area coverage up/downdraft not yet available from RAMS
              ! area should be near zero if no flux
              ! CHG(09/25/03) try to use RAUP1 ~ CFU1(cloud base) / (DENS(CB)*(1+sqrt(TKE*2) ))
              !             Get cloud base: first level w/ non-zero CFXUP
              KCB1=1
              KCB2=1
              KDB1=1
              !cloud base deep
              !dwen(20090315): add cfxup1 to meto type
              !dwen(20090826)              DO WHILE(METz%CFXUP1(KCB1).EQ.0.0.AND.KCB1.LT.NLVL)
              DO WHILE(METz(kcb1)%CFXUP1.EQ.0.0.AND.KCB1.LT.NLVL)
                 KCB1=KCB1+1
              END DO
              !cloud base shallow
              !dwen(20090315): add cfxup2 to meto type
              !dwen(20090826)              DO WHILE(METz%CFXUP2(KCB2).EQ.0.0.AND.KCB2.LT.NLVL)
              DO WHILE(METz(kcb2)%CFXUP2.EQ.0.0.AND.KCB2.LT.NLVL)
                 KCB2=KCB2+1
              END DO
              !downdraft (should allways be to ground)
              !dwen(20090315): add cfxdn1 to meto type
              !dwen(20090826)              DO WHILE(METz%CFXDN1(KDB1).EQ.0.0.AND.KDB1.LT.NLVL)
              DO WHILE(METz(kdb1)%CFXDN1.EQ.0.0.AND.KDB1.LT.NLVL)
                 KDB1=KDB1+1
              END DO
              !             get top of convection
              KCT1=NLVL
              KCT2=NLVL
              KDT1=NLVL
              !cloud top deep
              !dwen(20090826)              DO WHILE(METz%CFXUP1(KCT1).EQ.0.0.AND.KCT1.GT.1)
              DO WHILE(METz(kct1)%CFXUP1.EQ.0.0.AND.KCT1.GT.1)
                 KCT1=KCT1-1
              END DO
              !cloud top shallow
              !dwen(20090826)              DO WHILE(METz%CFXUP2(KCT2).EQ.0.0.AND.KCT2.GT.1)
              DO WHILE(METz(kct2)%CFXUP2.EQ.0.0.AND.KCT2.GT.1)
                 KCT2=KCT2-1
              END DO
              !downdraft top
              !dwen(20090826)              DO WHILE(METz%CFXDN1(KDT1).EQ.0.0.AND.KDT1.GT.1)
              DO WHILE(METz(kdt1)%CFXDN1.EQ.0.0.AND.KDT1.GT.1)
                 KDT1=KDT1-1
              END DO
              !             Extract TKEN and relative area coverage
              RAUP=0.0
              !dwen(20090315): add TKEF to DREC 
              !dwen(20090824)              IF (DREC(KG)%TKEF) THEN
              IF (DREC(pgrd(kp),ktime)%TKEF) THEN
                 !dwen(20090315): add DENS and TKEN to METO
                 !dwen(20090826)                IF (KCB2.LT.NLVL)RAUP=METz%CFXUP2(KCB2)/                   &
                 !dwen(20090826)                   (METz%DENS(KCB2)*(1+SQRT(METz%TKEN(KCB2)*2)))
                 IF (KCB2.LT.NLVL)RAUP=METz(kcb2)%CFXUP2/                   &
                      (METz(kcb2)%DENS*(1+SQRT(METz(kcb2)%TKEN*2)))
                 !dwen(20090826)                IF (KCB1.LT.NLVL)RAUP=METz%CFXUP1(KCB1)/                   &
                 !dwen(20090826)                   (METz%DENS(KCB1)*(1+SQRT(METz%TKEN(KCB1)*2)))
                 IF (KCB1.LT.NLVL)RAUP=METz(kcb1)%CFXUP1/                   &
                      (METz(kcb1)%DENS*(1+SQRT(METz(kcb1)%TKEN*2)))
              ELSE
                 !dwen(20090315): add SIGW to METO
                 !dwen(20090826)                IF (KCB2.LT.NLVL)RAUP=METz%CFXUP2(KCB2)/                   &
                 !dwen(20090826)                   (METz%DENS(KCB2)*(1+METz%SIGW(KCB2)*2))
                 IF (KCB2.LT.NLVL)RAUP=METz(kcb2)%CFXUP2/                   &
                      (METz(kcb2)%DENS*(1+METz(kcb2)%SIGW*2))
                 !dwen(20090826)                IF (KCB1.LT.NLVL)RAUP=METz%CFXUP1(KCB1)/                   &
                 !dwen(20090826)                   (METz%DENS(KCB1)*(1+METz%SIGW(KCB1)*2))
                 IF (KCB1.LT.NLVL)RAUP=METz(kcb1)%CFXUP1/                   &
                      (METz(kcb1)%DENS*(1+METz(kcb1)%SIGW*2))
              END IF
              !reset to environment
              IF(RAUP.EQ.0.0)ICNDX2=2
              !reset area coverage here since not entering CGRELL
              IF(KCT1.LT.KCB1)AREAPRU(KP)=0.0
              !reset area coverage here since not entering CGRELL
              IF(KCT2.LT.KCB2)AREAPRU(KP)=0.0
              !reset area coverage here since not entering CGRELL
              IF(KDT1.LT.KDB1)AREAPRD(KP)=0.0

              ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
              !             AREA=METO%GDIS*METO%GDIS
              AREA=METO%GDISX*METO%GDISY
              !        IF(KP.EQ.1)WRITE(*,*)'Z before',Z1,ICNDX(KP)
              ! Only call CGRELL if convection
              IF(RAUP.GT.0.0)THEN
                 ! assign relative area coverage at proper levels
                 !initialize with small number
                 !dwen(20090315): add RAUP1, RAUP2 and RADN1 to METZ
                 !dwen(20090315): replace NZM with NLVL
                 !                METO%RAUP1(1:NZM)=1.0E-15
                 !                METO%RAUP2(1:NZM)=1.0E-15
                 !                METO%RADN1(1:NZM)=1.0E-15
                 METz%RAUP1=1.0E-15
                 METz%RAUP2=1.0E-15
                 METz%RADN1=1.0E-15

                 !dwen(20090824) reformat the following lines due to compilation problem
                 !                IF(KCT1.GE.KCB1)METz%RAUP1(KCB1:KCT1)=RAUP
                 !                IF(KCT2.GE.KCB2)METz%RAUP2(KCB2:KCT2)=RAUP
                 !                IF(KDT1.GE.KDB1)METz%RADN1(KDB1:KDT1)=RAUP
                 IF(KCT1.GE.KCB1)METz(kcb1:kct1)%RAUP1=RAUP
                 IF(KCT2.GE.KCB2)METz(kcb2:kct2)%RAUP2=RAUP
                 IF(KDT1.GE.KDB1)METz(kdb1:kdt1)%RADN1=RAUP

                 !                WRITE(*,*)'hymodelc, before CGRELL',Z1,ICNDX(KP)
                 !dwen(20090315): add DFXUP1,DFXUP2,EFXUP1,EFXUP2,DFXDN1,EFXDN1,
                 !                    to METO
                 CALL CGRELL(METz%CFXUP1,METz%CFXUP2,METz%CFXDN1,        &
                      METz%DFXUP1,METz%DFXUP2,METz%EFXUP1,METz%EFXUP2,       &
                      METz%DFXDN1,METz%EFXDN1,                               &
                      METz%RAUP1,METz%RAUP2,METz%RADN1,                      &
                      AREA,                                                  &
                      AREAPRU(KP),AREAPRD(KP),METz%DENS,NLVL,                &
                      zprofm(1:nlvl),Z1,Z2,                     &
                      ICNDX(KP),ICNDX2,dummy,DT,BACK)
                 !convert back to sigma
                 ZPOS(KP)=1.0-Z2/ZMDL
                 !                WRITE(*,*)'Z after',Z2,ICNDX2
                 ! CHG(24/06/04) add ZFXCUM as vertical displacement due to convective flux
                 !  (deep or shallow, up or downdraft) along trajectory [m]
                 ZFXCUM(KP)=ZFXCUM(KP)+ABS(Z2-Z1)

              END IF !of IF(RAUP.GT.0.0)

              !assign changed cloud index (moved by CHG(3/4/2004))
              ICNDX(KP)=ICNDX2

           END IF !of IF(RAMSFLG.AND.NTURB.EQ.0)THEN
           !CCCCCCCCCCCCCCCCC END CONVECTION USING RAMS FLUXES CCCCCCCCCCCCCCCCCCCC

           ! CHG:need to set CONVDUR only for first call of ADVPNT during current particle loop
           IF(CONVDUR(KP).LT.CONVTMP)CONVDUR(KP)=CONVTMP
           CONVTMP=CONVDUR(KP)

           !****************************************************************

           !       convert to (km/min) = (gp/min) (km/gp)
           UBAR=UBAR*GRID(KGRID,KTIME)%SIZE

           !       pm10 dust emission algorithm
           IF(ICHEM.EQ.3)THEN
              !          negative mass is flag for recent particle emission
              !dwen(20090816): HYSPLIT4.6 supplemented P10ADJ
              IF(MASS(1,KP).LT.0.0) CALL P10ADJ                                   &
                   (P10F,SPOT(1)%OLVL,UBAR,METO%RAIN,SPOT(1)%AREA,MASS(1,KP),PGRD(KP))
           END IF

           !       skip terminated particles
           IF(PGRD(KP).EQ.0) CYCLE ploop
           KGRID=PGRD(KP)

           !       surface terrain height for this position
           ZSFC=METO%ZTER

           !       save maximim advection speed (0.06 km/min = 1 m/s)
           UMAX=MAX(0.060, UMAX, UBAR)

           !       increment particle age after advection (negative for no dispersion)
           PAGE(KP)=SIGN( (ABS(PAGE(KP))+INT(DT)), PAGE(KP) )


           !*****************************************************
           !dwen(20090306):following lines used in STILT
           ! JCL:      default is .FALSE. (=0)
           SEEVEG=0

           !-----------------------------------------------------------------------
           ! JCL:(11/03/03) store profile of wind errors
           IF(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3)THEN
              ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
              !        CONVFAC=METO%GDIS/60.0   !conversion factor for [grid/min]=>[m/s]
              CONVFACX=METO%GDISX/60.0
              CONVFACY=METO%GDISY/60.0

              !        impose VERTICAL correlation
              UERR(1)=GASDEV(RSEED,SIGUVERR)
              VERR(1)=GASDEV(RSEED,SIGUVERR)
              UERR2(1)=GASDEV(RSEED,SIGUVERR)
              VERR2(1)=GASDEV(RSEED,SIGUVERR)
              DO KZZ=2,NLVL
                 DELTAZ=(ZMDL-ZSFC)*(ZSG(KZZ-1)-ZSG(KZZ))
                 RAUTO=EXP(-1.0*DELTAZ/ZCORLEN)
                 UU=GASDEV(RSEED,SIGUVERR)
                 UERR(KZZ)=RAUTO*UERR(KZZ-1)+SQRT(1.0-RAUTO*RAUTO)*UU
                 UU2=GASDEV(RSEED,SIGUVERR)
                 UERR2(KZZ)=RAUTO*UERR2(KZZ-1)+SQRT(1.0-RAUTO*RAUTO)*UU2
                 VV=GASDEV(RSEED,SIGUVERR)
                 VERR(KZZ)=RAUTO*VERR(KZZ-1)+SQRT(1.0-RAUTO*RAUTO)*VV
                 VV2=GASDEV(RSEED,SIGUVERR)
                 VERR2(KZZ)=RAUTO*VERR2(KZZ-1)+SQRT(1.0-RAUTO*RAUTO)*VV2
              END DO

              !        impose TEMPORAL correlation--assume decorrelation process INDEPENDENT to that of VERTICAL process
              RAUTO=EXP(-1.0*ABS(DT)/TLUVERR)
              DO KZZ=1,NLVL
                 UUERRPREV=UVERR(KP)%U(KZZ)*SIGUVERR
                 VVERRPREV=UVERR(KP)%V(KZZ)*SIGUVERR
                 !in [grid/min]
                 !dwen(20090315):add UUPREV to METO
                 !dwen(20090826)           METz%UUPREV(KZZ)=UUERRPREV/CONVFACX+METz%UUPREV(KZZ)
                 METz(kzz)%UUPREV=UUERRPREV/CONVFACX+METz(kzz)%UUPREV
                 UUERR_T(KZZ)=RAUTO*UUERRPREV+SQRT(1.0-RAUTO*RAUTO)*UERR(KZZ)
                 !in [grid/min]
                 !dwen(20090315):add VVPREV to METO
                 !dwen(20090826)           METz%VVPREV(KZZ)=VVERRPREV/CONVFACY+METz%VVPREV(KZZ)
                 METz(kzz)%VVPREV=VVERRPREV/CONVFACY+METz(kzz)%VVPREV
                 VVERR_T(KZZ)=RAUTO*VVERRPREV+SQRT(1.0-RAUTO*RAUTO)*VERR(KZZ)
                 !DO KZZ=1,NLVL
              END DO

              !        impose SPATIAL correlation--assume decorrelation process INDEPENDENT to that of VERTICAL process
              ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
              !        CFAC=METO%GDIS/1000.0   ![m/grid]=>[km/grid]
              ![m/grid]=>[km/grid]
              CFACX=METO%GDISX/1000.0
              ![m/grid]=>[km/grid]
              CFACY=METO%GDISY/1000.0
              !distance travelled [km]
              DIST=SQRT((CFACX*(XPOS(KP)-XPOSPREV))**2+                     &
                   &              (CFACY*(YPOS(KP)-YPOSPREV))**2)
              RAUTOH=EXP(-1.0*DIST/HORCORLEN)
              DO KZZ=1,NLVL
                 UUERRNEXT=RAUTOH*UUERR_T(KZZ)+                               &
                      &                SQRT(1.0-RAUTOH*RAUTOH)*UERR2(KZZ)
                 !in [grid/min]
                 !dwen(20090315):add UUNEXT and VVNEXT to METO
                 !dwen(20090826)           METz%UUNEXT(KZZ)=UUERRNEXT/CONVFACX+METz%UUNEXT(KZZ)
                 METz(kzz)%UUNEXT=UUERRNEXT/CONVFACX+METz(kzz)%UUNEXT

                 VVERRNEXT=RAUTOH*VVERR_T(KZZ)+                               &
                      &                SQRT(1.0-RAUTOH*RAUTOH)*VERR2(KZZ)
                 !in [grid/min]
                 !dwen(20090826)           METz%VVNEXT(KZZ)=VVERRNEXT/CONVFACY+METz%VVNEXT(KZZ)
                 METz(kzz)%VVNEXT=VVERRNEXT/CONVFACY+METz(kzz)%VVNEXT
                 !          store for next timestep
                 UVERR(KP)%U(KZZ)=UUERRNEXT/SIGUVERR
                 UVERR(KP)%V(KZZ)=VVERRNEXT/SIGUVERR
                 !           WRITE(*,*)"SPATIAL :",KZZ,UUERR_T(KZZ),
                 !     &           UUERRNEXT,METO%UUNEXT(KZZ)*CONVFAC,UERR2(KZZ),RAUTOH
                 !DO KZZ=1,NLVL
              END DO

              !IF(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3)THEN
           END IF
           !-----------------------------------------------------------------------
           ! CHG:(28/04/06) get zi errors: here relative error in percent
           IF(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3)THEN

              !        impose TEMPORAL correlation
              RAUTO=EXP(-1.0*ABS(DT)/TLZIERR)
              ZIERR(KP)=RAUTO*ZIERR(KP)+SQRT(1.0-RAUTO*RAUTO)*GASDEV(RSEED,SIGZIERR)

              !        impose SPATIAL correlation--assume decorrelation process INDEPENDENT to that of VERTICAL process
              ![m/grid]=>[km/grid]
              CFACX=METO%GDISX/1000.0
              ![m/grid]=>[km/grid]
              CFACY=METO%GDISY/1000.0
              !distance travelled [km]
              DIST=SQRT((CFACX*(XPOS(KP)-XPOSPREV))**2+                     &
                   &              (CFACY*(YPOS(KP)-YPOSPREV))**2)
              RAUTOH=EXP(-1.0*DIST/HORCORZILEN)
              ZIERR(KP)=RAUTOH*ZIERR(KP)+SQRT(1.0-RAUTOH*RAUTOH)*GASDEV(RSEED,SIGZIERR)
              !IF(WINDERRTF.EQ.1.OR.WINDERRTF.EQ.3)THEN
           END IF
           !-----------------------------------------------------------------------

           !            WRITE(*,*)'before PARDSP, XPOS, YPOS:',XPOSPREV,YPOSPREV
           !CCCCCCCCCCCFollowing comments are for PARDSP
           ! JCL:      add SEEVEG as argument--gets set to 1 if particle sees veg
           !           also add TLFRAC--fraction of TL that fixes internal timestep
           !           also add RSEED (random seed) as argument
           !           also add SIGMAW,SIGMAU (std dev of vert&hor velocity) as argument
           !           SAMPTT is amt of time [min.] that particle 'sees' the ground
           !           METO%TL is vertical profile of Lagrangian timescales [s]
           !           METO%SIGW is vertical profile of std dev of vert velocity
           !           METO%ZMLNEXT is mixed-layer ht [m]
           !           METO%DENSNEXT is vertical profile air density [kg/m3]
           ! JCL:(4/13/2000)
           !           METO%TLPREV,METO%SIGWPREV,METO%ZMLPREV,METO%DENSPREV are values from
           !                 BEFORE mean advection timestep
           ! CHG:      SAMPTT2 is variable for output from pardsp, now Tlagrange vert.
           !               (different from 'SIGMAW'--the OUTPUT from PARDSP)
           ! JCL:(4/27/00)Since want vertical profiles to be CONSTANT, not have them as argument
           ! JCL:(5/2/00) also remove ZSG as argument
           ! JCL:(5/18/00)add as 1st argument--# of min. since start
           ! JCL:(5/18/00)WWPREV is array storing NORMALIZED turbulent velocity for each particle
           !              WWPREV gets updated in PARDSP
           ! JCL:(5/18/00)Arguments when have TIME-VARYING vertical profiles of TL & SIGMAW
           ! CHG:      NTURB true switched off turbulence
           ! JCL:(6/29/00)'VEGHT' is ht [m] below which a particle would be counted as 'seeing' grd vegetation
           !                   VEGHT determines whether SAMPTT is > 0 or not
           ! JCL:(5/9/01)added horizontal position before mean advection (XPOSPREV & YPOSPREV)
           !                  + vertical profiles of horizontal winds (UUPREV,UUNEXT,VVPREV,VVNEXT)
           !             -these are used to implement interaction between vertical turbulence and wind shear
           ! JCL:(4/3/02)added vertical profile of mass violation before & after mean advection
           ! JCL:(4/3/02)added weighting of particle due to mass violation experienced by particle in PARDSP (DMASSWT)
           ! CHG:(12/05/01) add 'CONVDUR' (calc. in ADVPNT, reset to 0 after conv. redistribution was done in PARDSP)
           ! CHG:(12/04/01)add 'ZLOC' limit of convection altitude (m)
           ! CHG(09/11/03) Pass on RAMSFLG
           ! CHG:(03/17/2004) pass on rel. humidity and temperature profile to get dry air density column integrated
           ! CHG:(03/17/2004) also add output variable for specific humidity SPHU and FOOT
           ! JCL:(07/14/2004) split up GDIS=>GDISX&GDISY for global grids (non-conformal)
           SMIN=JET+INT(DT)-(CONC(1)%START%MACC)

           IF (NTURB == 0) THEN

              !       adjusts advection to simulate particle dispersion
              !******************************************************
              !dwen(20090816): HYSLPIT4.7 splitted HMIX into UMIX and VMIX 
              !dwen(20090816): HYSPLIT4.7 used WMIX,instead of VMIX, to stand for vertical turbulence
              !dwen(20090816): HYSPLIT started from v4.9 to use ISEED,instead of RSEED, for the random seed that will be
              !                  fed to PARVAR
              !dwen(20090816): remove ISOT,back,ecmflg and SEEVEG because they are useless any more in PARDSP  
              !dwen(20090319): Starting from v4.5, HYSPLIT splitted METO into METZ and METO 
              !                METZ defines profile advection variables, and METO defines 
              !                advection surface variables.
              !dwen(20090816): Starting from v4.7, HYSPLIT put WMIX(vertical turbulence profile) to METZ, instead of METO
              !                 more check is needed.
              !        CALL PARDSP(METO%UMIX,METO%VMIX,METO%GDISX,METO%GDISY,DT,ZMDL,ZSFC,    &
              !                    NLVL,METZ%WMIX,ZSG,XPOS(KP),YPOS(KP),ZPOS(KP),             &
              !                    SIGV(KP),SIGW(KP),SIGU(KP),HDWPX,METO%ZNDX,                &
              !                    VX(job_id+ens_id) )

              xpos(kp)=xposprev    
              ypos(kp)=yposprev
              CALL PARDSP(METO%UMIX,METO%VMIX,METO%GDISX,METO%GDISY,DT,ZMDL,ZSFC,    &
                   !dwen(20090906)                    NLVL,METZ%WMIX,ZSG,XPOS(KP),YPOS(KP),ZPOS(KP),             &
                   NLVL,METZ%WMIX,ZSG,xpos(kP),ypos(kp),ZPOS(KP),             &
                   SIGV(KP),SIGW(KP),SIGU(KP),HDWPX,METO%ZNDX,                &
                   VX(job_id+ens_id),                                         &
                   metz%tlprev, metz%sigwprev, meto%zmlprev, metz%densprev,       &
                   metz%tl, metz%sigw, meto%zmlnext, metz%dens, tlfrac,           &
                   !dwen(20090906)                xposprev, yposprev,                                            &
                   metz%uuprev, metz%uunext, metz%vvprev, metz%vvnext,            &
                   metz%dmassprev, metz%dmassnext, dmasswt(kp),                   &
                   rseed, sigmaw, sigmau, samptt, samptt2, wwprev(kp), veght,     &
                   convdur(kp), meto%zlocnext, ramsflg,                           &
                   metz%temp, metz%tempprev, metz%rhfr, metz%rhfrprev,            &
                   sphu, foot)


              !dwen(20090816)             CALL PARDSP(                                                      &
              !dwen(20090816)                METO%HMIX, METO%GDISX, METO%GDISY, DT, ZMDL, ZSFC, NLVL,       &
              !dwen(20090816)                METz%wMIX, ZSG, XPOS(KP), YPOS(KP), ZPOS(KP), SIGv(KP),        &
              !dwen(20090816)                SIGw(KP), SIGu(KP), HDWP(KP), METO%ZNDX,                       &
              !dwen(20090816)                METO%TLPREV, METO%SIGWPREV, METO%ZMLPREV, METO%DENSPREV,       &
              !dwen(20090816)                METO%TL, METO%SIGW, METO%ZMLNEXT, METO%DENS, TLFRAC,           &
              !dwen(20090816)                XPOSPREV, YPOSPREV,                                            &
              !dwen(20090816)                METO%UUPREV, METO%UUNEXT, METO%VVPREV, METO%VVNEXT,            &
              !dwen(20090816)                METO%DMASSPREV, METO%DMASSNEXT, DMASSWT(KP),                   &
              !dwen(20090816)                RSEED, SIGMAW, SIGMAU, SAMPTT, SAMPTT2, WWPREV(KP), KP, VEGHT, &
              !dwen(20090816)                CONVDUR(KP), METO%ZLOCNEXT, RAMSFLG,                           &
              !dwen(20090816)                METO%TEMP, METO%TEMPPREV, METO%RHFR, METO%RHFRPREV,            &
              !dwen(20090816)                SPHU, FOOT)
           ELSE
              SIGMAW  = 0d0                                ! define output values
              SAMPTT  = 0d0
              SAMPTT2 = 0d0
              FOOT    = 0d0
           END IF
           !
           !******************************************************


           !******************************************************
           !dwen(20090306):following lines used in STILT
           ! JCL:(11/03/03) apply transport error by directly changing wind profile
           ! JCL:(8/28/03)calculate total mass violation experienced by particle (fraction of mass in gridcell) for MEAN trajs
           ! CHG(10/10/03) different for RAMS
           IF (.NOT. RAMSFLG .AND.  NTURB == 1) THEN
              KZ=MAX0(1, MIN0(NLVL,INT(METO%ZNDX)))
              KT=MIN0(KZ+1,NLVL)
              IF(ZPOS(KP).LT.ZSG(1))THEN
                 !dwen(20090315):add DMASSPREV and DMASSNEXT to METO
                 !dwen(20090826)                  DMASS=(METz%DMASSPREV(KT)+METz%DMASSNEXT(KT))/2.0
                 DMASS=(METz(kt)%DMASSPREV+METz(kt)%DMASSNEXT)/2.0
              ELSE
                 !dwen(20090826)                  DMASS=(METz%DMASSPREV(1)+METz%DMASSNEXT(1))/2.0
                 DMASS=(METz(1)%DMASSPREV+METz(1)%DMASSNEXT)/2.0
              END IF
              DMASSWT(KP)=DMASSWT(KP)*(1.0+DMASS*ABS(DT))
           END IF
           IF(RAMSFLG.AND.NTURB.EQ.1)THEN
              KT=MAX0(1, MIN0(NLVL,INT(METO%ZNDX)+1))
              !dwen(20090826)              DMASS=(METz%DMASSPREV(KT)+METz%DMASSNEXT(KT))/2.0
              DMASS=(METz(kt)%DMASSPREV+METz(kt)%DMASSNEXT)/2.0
              DMASSWT(KP)=DMASSWT(KP)*(1.0+DMASS*ABS(DT))
           END IF

           ! JCL:(6/29/00)increment cumulative SAMPTT
           SAMPTTCUM(KP)=SAMPTTCUM(KP)+SAMPTT
           ! CHG:(28/04/2006) apply ZI error to FOOT
           IF(WINDERRTF.EQ.2.OR.WINDERRTF.EQ.3)THEN
              RELZIERR=(100+ZIERR(KP))/100
              RELZIERR=MAX(RELZIERR,0.0)
              FOOT=FOOT*RELZIERR
              !         increment cumulative FOOT with ZI error
              FOOTCUM(KP)=FOOTCUM(KP)+FOOT
              !              WRITE(45,*)KP,JET,RELZIERR,FOOTCUM(KP),FOOT
           ELSE
              ! CHG:(03/18/2004)increment cumulative FOOT
              FOOTCUM(KP)=FOOTCUM(KP)+FOOT
           END IF
           !******************************************************

           !       computes vertical and horizontal puff mixing
           !       requires complex mode rather than simple


           !******************************************************
           !dwen(20090817): HYSPLIR4.5 removed MASS, useless in PUFDSP
           !dwen(20090816): HYSLPIT4.7 splitted HMIX into UMIX and VMIX 
           !dwen(20090817): HYSPLIT from v4.7 uses WMIX to represent vertical turbulence(m2/s2)
           !                        UMIX u-component turbulence(m2/s2) and VMIX v-component turbulence(m2/s2)
           !dwen(20090817): VMIX in HYSPLIT earlier than v4.7 represented vertical mixing(m2/s) 
           !                        and HMIX represented horizontal diffusivity(m2/s)
           !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
           !                           SIGU and SIGV for horizontal velocity sigma 
           !dwen(20090319): HYSPLIT4.5 splitted METO into METZ and METO 
           !                METZ defines profile advection variables, and METO defines 
           !                advection surface variables.
           !dwen(20090816): HYSPLIT started from v4.7 to put WMIX(vertical turbulence profile) METZ, instead of METO
           !dwen(20090816): HYSPLIT from v4.7 removed ISOT

           CALL PUFDSP(KPUFF,METO%UMIX,METO%VMIX,DT,ZMDL,ZSFC,NLVL,METZ%WMIX,ZSG, &
                PAGE(KP),ZPOS(KP),SIGH(KP),SIGW(KP),HDWP(KP),METO%ZNDX)

           !            CALL PUFDSP(METO%HMIX,DT,ZMDL,ZSFC,NLVL,METO%VMIX,ZSG,      &
           !     &         MASS(:,KP),ZPOS(KP),SIGH(KP),SIGV(KP),HDWP(KP),          &
           !     &         METO%ZNDX,ISOT)
           !******************************************************

           !       simple chemistry conversion option
           !dwen(20090817) HYSPLIT from v4.7 supplements CHEMO2
           IF(ICHEM.EQ.2) CALL CHEM02(DIRT,DT,MASS(:,KP))

           !       optional dry, wet, and decay routines

           !***************************************************
           !dwen(20090306): HYSPLIT4.5 removed DIRT and CONC from COMMON bloc and added to argument list
           !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
           !                           SIGU and SIGV for horizontal velocity sigma 
           !               SWF, substituted by meto%dswf, is no longer used since HYSPLIT4.6 
           !dwen(20090311): HYSPLIT from v4.7 uses HDWPX instead of HDWP for complex and simple modes
           !               HDWPX is used for accounting for complex and simple modes
           !               hdwpx=kdwp(kp) for complex
           !               IF(HDWPX.GE.100)HDWPX=MOD(HDWP(KP)/10,10)  for simple mode
           !dwen(2009081&): HYSPLIT from v4.7 added DRYD to argument list for dry deposition velocity
           !dwen(20090319): HYSPLIT4.5 splitted METO into METZ and METO 
           !                METZ defines profile advection variables, and METO defines 
           !                advection surface variables.
           !dwen(20090817): From v4.5, HYSPLIT puts DENS, TEMP AND RHFR into METZ instead of METO
           !dwen(20090817): From v4.7, HYSPLIT adds ICHEM for special deposition options

           IF(CDEP) CALL DEPELM(                                                  &
                DIRT,SPOT(1)%OLAT,SPOT(1)%IBMO,NLVL,METO%ZNDX,DT,ZMDL,ZSFC,   &
                ZSG,MASS(:,KP),DEPT,ZPOS(KP),SIGW(KP),ICHEM,PTYP(KP),         &
                METO%LAND,METO%AERO,SFCL,METO%USTR,METO%PSI,meto%dswf,        &
                HDWPX,METO%RAIN,METZ%DENS,METZ%TEMP,METZ%RHFR,DRYD)

           !     &         CALL DEPELM(NLVL,METO%ZNDX,DT,ZMDL,ZSFC,ZSG,             &
           !     &         MASS(:,KP),DEPT,ZPOS(KP),SIGV(KP),PTYP(KP),              &
           !     &         METO%LAND,METO%AERO,SFCL,METO%USTR,METO%PSI,             &
           !     &         SWF,HDWP(KP),METO%RAIN,METO%DENS,METO%TEMP,METO%RHFR)
           !***************************************************

           !       surface water transport option
           IF(ICHEM.EQ.7)THEN
              IF(HDWPX.EQ.5)THEN 
                 !             hdwp=5 only if the particle has deposited over a water surface
                 !             to show movement particle redeposits all mass each time step 
                 DEPT=MASS(:,KP)

                 IF(METO%LAND.NE.7)THEN
                    !                the overwater particle cannot move over land 
                    PGRD(KP)=0
                    MASS(:,KP)=0.0
                    DEPT=0.0
                 END IF

                 !             reset in case particle has deposited in the DEPELM step
                 HDWP(KP)=HDWPX

              ELSE
                 !             show no deposition for undeposited particles 
                 DEPT=0.0
              END IF
           END IF



           IF(ICHEM.EQ.4)THEN
              !          sum mass to concentration grid on meteo data projection
              CALL METSUM(CONC,NUMGRD,XPOS(KP),YPOS(KP),METO%GDISX,METO%GDISY,    &
                   DT,JET,ZMDL,ZSFC,KMSL,CGSIZE,MASS(:,KP),DEPT,ZPOS(KP),  &
                   SIGH(KP),SIGW(KP),HDWPX,PTYP(KP),CSUM)
           ELSE
              IF(CMASS.EQ.0)THEN
                 !          sum mass to grid as concentration grid on standard lat/lon projection

                 !***************************************************
                 ! JCL:      add BACK as argument to have more complicated conditions
                 !           that would still alter conc array when JET is decreasing
                 !           sum values to array
                 !dwen(20090306): HYSPLIT4.5 removed CONC from COMMON bloc and added to argument list
                 !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
                 !                           SIGU and SIGV for horizontal velocity sigma 
                 !dwen(20090816): Starting from v4.9, HYSPLIT uses HDWPX for complex and simple modes
                 !               HDWPX is used for accounting for complex and simple modes
                 !               hdwpx=kdwp(kp) for complex
                 !               IF(HDWPX.GE.100)HDWPX=MOD(HDWP(KP)/10,10)  for simple mode
                 !dwen(20090817): From v4.4, HYSPLIT adds KMSL for msl or agl option.
                 !                KMSL is starting height default:  0 for AGL=0 or 1 for MSL(relative to mean sea level)
                 !                 starting height can be changed to meters MSL or meters AGL through KMSL 
                 !dwen(20090815): HYSPLIT uses DT to determine BACK flag 


                 CALL CONSUM(CONC,NUMGRD,METO%PLAT,METO%PLON,DT,JET,ZMDL,ZSFC,KMSL,  &
                      CGSIZE,MASS(:,KP),DEPT,ZPOS(KP),SIGH(KP),SIGW(KP),      &
                      HDWPX,PTYP(KP),CSUM)
                 !            CALL CONSUM(NUMGRD,METO%PLAT,METO%PLON,DT,JET,ZMDL,ZSFC,    &
                 !     &         CGSIZE,MASS(:,KP),DEPT,ZPOS(KP),SIGH(KP),SIGV(KP),       &
                 !     &         HDWP(KP),PTYP(KP),CSUM,BACK)
                 !*****************************************************
              ELSE
                 !          sum mass to grid as mass on standard lat/lon projection
                 CALL MASSUM(CONC,NUMGRD,METO%PLAT,METO%PLON,DT,JET,ZMDL,ZSFC,KMSL,  &
                      CGSIZE,MASS(:,KP),ZPOS(KP),PTYP(KP),CSUM)
              END IF
           END IF

           IF(LAGSAM.AND.HDWPX.EQ.6)                                              &
                !dwen(20090817): HYSPLIT v4.9 supplemented LAGSUM
                CALL LAGSUM(LAGS,CONC,NUMGRD,METO%PLAT,METO%PLON,ZPOS(KP),DT,JET,   &
                ZMDL,ZSFC,KMSL,MASS(:,KP),PAGE(KP),PTYP(KP),CSUM)

           !       sum mass for diagnostic analysis
           IF(INITK.EQ.1.OR.INITK.EQ.2)THEN
              !          vertical puffs summed into all levels within puff
              SGT=MAX(ZPOS(KP)-1.54*SIGW(KP), 0.0)
!!$              ZZ=ZMDL*(1.0-MIN(1.0,SGT))
!!$
!!$              !          convert to vertical index using equation for zsfc=0
!!$              DIST=(BB*BB-4.0*AA*(CC-ZZ))
!!$              IF(DIST.GE.0.0)THEN
!!$                 ZX=(-BB+SQRT(DIST))/(2.0*AA)
!!$              ELSE
!!$                 ZX=1.0
!!$              END IF
              IF (.NOT.RAMSFLG) CALL ind_zsg (zmdl,zsg,nlvl,sgt,zx,aa,bb,cc)
              K2=MIN(MAX(1,NINT(ZX)),NLVL)

              SGB=MIN(ZPOS(KP)+1.54*SIGW(KP), 1.0)
!!$              ZZ=ZMDL*(1.0-MIN(1.0,SGB))
!!$
!!$              !          convert to vertical index using equation for zsfc=0
!!$              DIST=(BB*BB-4.0*AA*(CC-ZZ))
!!$              IF(DIST.GE.0.0)THEN
!!$                 ZX=(-BB+SQRT(DIST))/(2.0*AA)
!!$              ELSE
!!$                 ZX=1.0
!!$              END IF
              IF (.NOT.RAMSFLG) CALL ind_zsg (zmdl,zsg,nlvl,sgb,zx,aa,bb,cc)
              K1=MIN(MAX(1,INT(ZX)),NLVL)
           ELSE
              !          particles summed into nearest meteo index
              K1=MAX(1, MIN(NLVL, NINT(METO%ZNDX) ))
              K2=K1
           END IF

           FRAC=FLOAT(K2-K1+1)
           DO KZ=K1,K2
              PMASS=MASS(1,KP)
              ZMASS(KZ)=ZMASS(KZ)+MASS(1,KP)/FRAC
              MM=MAXDIM
              DO WHILE(MM.GT.1)
                 PMASS=PMASS+MASS(MM,KP)
                 ZMASS(KZ)=ZMASS(KZ)+MASS(MM,KP)/FRAC
                 MM=MM-1
              END DO
              TMASS=TMASS+PMASS/FRAC
           END DO
           !       mark zero mass for removal (except lagrangian samplers) 
           IF(PMASS.EQ.0.0.AND.HDWP(KP).NE.6)PGRD(KP)=0

           !********************************************************************
           !dwen(20090309):following lines from STILT
           ! JCL:      skip the steps reserved for particles which have moved off grid
           GO TO 400

           ! JCL:      Following are steps reserved for particles that have moved off grid
           !CCCCCCCCCCCCCCCCCCCCCCC
           !           terminated particles skip to here
200        CONTINUE

           ! JCL:      also increment particle age even if it has moved off grid
           !           increment particle age after advection
           PAGE(KP)=PAGE(KP)+DT

           ! JCL:(4/27/01)increment counter of number of particles that have moved off grid
           COUNTNPAROUT=COUNTNPAROUT+1

           ! JCL:      not write out particle results when particle has moved off grid
           !           so branch off to end of particle loop
           GO TO 500

           !CCCCCCCCCCCCCCCCCCCCCCC
           ! JCL:
400        CONTINUE

           ! JCL:(11/1/02)use Draxler formulation, with 'terrain compression factor'
           ! JCL:      follwing lines are to write out each particle's position to file
           !           convert sigma to AGL
           ! CHG(09/10/03) correct transformation between sigma and agl
           !           ZOUT=(1.0-ZPOS(KP))*(ZMDL-ZSFC)
           !           ZOUT=(1.0-ZPOS(KP))*ZMDL*(ZMDL/(ZMDL-ZSFC))
           ZOUT=(1.0-ZPOS(KP))*(ZMDL-ZSFC)

           ! JCL:(3/1/01)convert 'WWOUT' from [sigma/min]=>[m/s]
           !           multiply by (-1.0*DT/ABS(DT)) b/c want to keep sign of direction
           WWOUT=(-1.0*DT/ABS(DT))*WWOUT*(ZMDL-ZSFC)/60.0

           ! JCL:(07/12/2004) implement global lat/lon code from HYSPLIT Ver 45
           !dwen(20090824)            IF(GRID(PGRD(KP))%LATLON)THEN
           IF(GRID(PGRD(KP),ktime)%LATLON)THEN

              !dwen(20090315): add KT,ktime?
              !              CALL GBL2LL(PGRD(KP),XPOS(KP),YPOS(KP),YOUT,XOUT)
              CALL GBL2LL(PGRD(KP),ktime,XPOS(KP),YPOS(KP),yout,xout)
           ELSE
              ! JCL:        call CXY2LL to convert positions to lat/lon before dump
              !dwen(20090315): remove GRID(KG)%proj,add KT
              !              CALL CXY2LL(GRID(PGRD(KP))%GBASE,XPOS(KP),YPOS(KP),       &
              !     &                     YOUT,XOUT, GRID(PGRD(KP))%proj)
              CALL CXY2LL_wps(GRID(PGRD(KP),ktime)%GBASE,XPOS(KP),YPOS(KP),       &
                   YOUT,XOUT,GRID(PGRD(KP),ktime)%proj)
           END IF

           ! JCL:(2/28/01)add call to SUNANG b/c need solar angle to output incident solar radiation
           !               feed current position instead of starting position, as original call does
           ! CHG:(11/20/01) use downward short.wave radiation from met data,
           !           if available
           ! CHG:(9/24/02) this is flux interpolated from nearest 2 analysis times
           ! NEW: interpolation of cloud effect rather then interpolation of radiation
           IF(meto%dswf.GE.0)THEN
              SWF=meto%dswf
           ELSE
              !dwen(20090816): HYSPLIT4.6 removed IBYR from argument list and determined it with tm2day
              !            CALL SUNANG                                                 &
              !     &         (SPOT(1)%IBYR,JET,YOUT,XOUT,EA,SEA)
              CALL SUNANG(JET,yout,xout,EA,SEA)      ! SUNANG was called this way in HYSPLIT4.9

              ! JCL:(2/28/01)add call to SUNFLX to derive solar flux
              !           compute solar flux from solar angle and humidity
              !           required for resistance gaseous dry deposition

              !dwen(20090315): add RHFR to METO
              CALL SUNFLX(NLVL,SEA,METz%RHFR,SWF,TR)
           END IF

           ! JCL:      write out results to PARTICLE.DAT only in specified intervals
           !             and OUTDT=0 means that results would be output EACH timestep
           !           'WRITEOUT' is a logical flag that gets set before particle loop

           IF(WRITEOUT)THEN

              ! CHG: (20/11/01) write also rain
              ! CHG: (08/11/02) CHANGED FOR H2O budget study
              ! JCL:      write to 'particle.dat'  1)mins since start
              !           2)particle index  3)particle LAT  4)particle LON
              !           5)particle altitude  6) terrain height
              !           7)standard deviation of vertical velocity [m/s]
              !           8)air temperature at lowest model layer [K]
              !           9)amount of time that particle 'sees' the ground [min]
              !           10)downward solar short wave radiation [W/m2]
              !           11)vertical mean wind [m/s] -->CHANGED to temp [k]
              !           12)ht of mixed-layer [m]
              !           13)total rain fall rate [m/min]
              !           14)convective rain fall rate [m/min]
              !           15)sensible heat flux [w/m2] --> CHANGED to loc. dens
              !           16)latent heat flux [w/m2] (NGM: kg/m2/s)
              !           17)total cloud cover (%)
              !           18)particle weighting due to mass violation
              !           19)Lagrangian timescale [s] -->CHANGED to RH fraction
              ! CHG:         Only write all the stuff for particle trajectories
              ! CHG: (20/11/01) write also rain

              !-------------------------------------------------------------------------
              ! JCL(060929): determine vertical indices to properly extract the values of variables at proper height
              KZZ=1
              DO WHILE(ZSG(KZZ+1).GT.ZPOS(KP))
                 !CHG i.e. if particle below HGT(2), KB=1 and KT=2; also if particle below HGT(1)
                 KZZ=KZZ+1
              END DO
              KTT=KZZ+1
              KBB=KTT-1
              !            for RAMS use also different vertical index to extract scalars
              IF(RAMSFLG)THEN
                 KZZ=1
                 DO WHILE(ZSG(KZZ).GT.ZPOS(KP))
                    KZZ=KZZ+1
                 END DO
                 KTT=KZZ
                 KBB=KTT-1
              END IF

              ! JCL(060929): vertical interpolation
              IF(KBB.EQ.0.AND.RAMSFLG)KBB=1
              ZF=(ZSG(KBB)-ZPOS(KP))/(ZSG(KBB)-ZSG(KTT))
              !            if below 1st level, near the ground
              IF(ZPOS(KP)>ZSG(1))ZF=0.0
              !            since don't have pressure profile @ t (only have at t+dt), just use the values extracted after mean advection
              ! dwen(20090315): add PRES and TEMP to METO
              !dwen(20090826)             PRESLOCAL=METz%PRES(KBB)+ZF*(METz%PRES(KTT)-METz%PRES(KBB))
              !dwen(20090826)             TEMPLOCAL=METz%TEMP(KBB)+ZF*(METz%TEMP(KTT)-METz%TEMP(KBB))
              PRESLOCAL=METz(kbb)%PRES+ZF*(METz(ktt)%PRES-METz(kbb)%PRES)
              TEMPLOCAL=METz(kbb)%TEMP+ZF*(METz(ktt)%TEMP-METz(kbb)%TEMP)
              !            NOTE: this is more appropriate for density at final transported particle location than
              !            METO%DENSLOCAL (which uses vertical position before PARDSP)
              !dwen(20090826)             DENSLOCAL=METz%DENS(KBB)+ZF*(METz%DENS(KTT)-METz%DENS(KBB))
              !dwen(20090826)             RHFRLOCAL=METz%RHFR(KBB)+ZF*(METz%RHFR(KTT)-METz%RHFR(KBB))
              DENSLOCAL=METz(kbb)%DENS+ZF*(METz(ktt)%DENS-METz(kbb)%DENS)
              RHFRLOCAL=METz(kbb)%RHFR+ZF*(METz(ktt)%RHFR-METz(kbb)%RHFR)

              !-------------------------------------------------------------------------

             !-----------------------------------------------------
              ! tk (05/17/2011)
              ! some problems exists for output format, if dmas value is gt 100000
               if (DMASSWT(kp) .gt. 100000) then 
                   DMASSWT(kp)=100000
                  endif

              !-------------------------------------------------------------------------
              ! JCL(02/28/2004):  Call output routine written by Marcos
              ! CHG:(03/17/2004) also add output variable for specific humidity SPHU and FOOT
              !dwen(20090315):add ZMLNEXT,ZMLPREV,CRAI,SHTF,LTHF,TCLD,LCLD,SOLW,ZLOCPREV AND ZLOCNEXT to METO
              CALL OUTPUT(KFPARDAT,IVMAX,VARSIWANT(1:IVMAX),NTURB,             &
                   JET+INT(DT)-(CONC(1)%START%MACC),KP,ICNDX(KP), &
                   !dwen20090826)                        XOUT,YOUT,ZOUT,ZSFC,SIGMAW,METz%TEMP(1),       &
                   XOUT,YOUT,ZOUT,ZSFC,SIGMAW,METz(1)%TEMP,       &
                   SAMPTTCUM(KP),FOOTCUM(KP),SWF,WWOUT,           &
                   (METO%ZMLNEXT+METO%ZMLPREV)/2.0,METO%RAIN,     &
                   METO%CRAI,ZFXCUM(KP),METO%SHTF,METO%LHTF,METO%TCLD, &
                   DMASSWT(KP),SAMPTT2,DENSLOCAL,                 &
                   RHFRLOCAL,SPHU,METO%SOLW,METO%LCLD,            &
                   METO%ZLOCNEXT,PRESLOCAL,TEMPLOCAL)
              ! JCL:              reinialize the cumulative SAMPTT
              SAMPTTCUM(KP)=0.0
              FOOTCUM(KP)=0.0
              ZFXCUM(KP)=0.0
              !-------------------------------------------------------------------------

              ! JCL:      end of IF(WRITEOUT) statement
           END IF
           ! JCL:      branch from particles which have moved off grid
500        CONTINUE
           !
           !********************************************************************
           !    particle/puff loop
        END DO ploop


        !********************************************************************
        !dwen(20090309):following lines from STILT
        ! JCL:   reset flag and output counter
        IF(WRITEOUT)THEN
           WRITEOUT=.FALSE.
           COUNTOUT=0.0
        END IF

        ! CHG:(12/05/01) reset temporary convtmp variable
        CONVTMP=0
        !********************************************************************

        !-------------------------------------------------------------------------------
        !    end of particle loop ... now call routines requiring global data
        !-------------------------------------------------------------------------------

        !    update maximum concentration array if requested
        IF(KMAXC.NE.0.AND.KSNAP.NE.0)  &  
             CSUM(:,:,:,:,KMAXC)=MAX(CSUM(:,:,:,:,KMAXC),CSUM(:,:,:,:,KSNAP))
        IF(KMAXC.NE.0.AND.KAVRG.NE.0)  &  
             CSUM(:,:,:,:,KMAXC)=MAX(CSUM(:,:,:,:,KMAXC),CSUM(:,:,:,:,KAVRG))

        !    decay of ground-level deposition amounts
        !*************************************************************
        !dwen(20090311):used updated DEPRAD
        !dwen(20090306): HYSPLIT4.5 removed DIRT and CONC from COMMON bloc and added to argument list
        IF(CDEP) CALL DEPRAD(CONC,DIRT,NUMGRD,NUMTYP,DT,CSUM)
        !         IF(CDEP) CALL DEPRAD(NUMGRD,NUMTYP,DT,CSUM)
        !**************************************************************

        !    increment elapsed time
        JET=JET+INT(DT)

        !---------------------------------------
        !    STANDARD SINGLE PROCESSOR SECTION

        !    check each time step for output to disk and zero-out array

        !*************************************************************
        ! JCL: add BACK as argument to subroutine to avoid some tests that
        !       would cause program to terminate when program is run backwards
        !dwen(20090306): HYSPLIT4.5 removed DIRT and CONC from COMMON bloc and added to argument list
        !dwen(20090319): HYSPLIT4.7 added KTIME and KGRID to argument lis to make  
        !                 simultaneous multiple meteorology and grids available  
        !dwen(20090817): HYSPLIT from v4.7 added ICHEM for ichem=4 option for meteo grid

        CALL CONDSK(CONC,DIRT,ICHEM,KGRID,KTIME,NUMGRD,NUMTYP,DT,JET,IFHR,    &
             CPACK,CSUM)

        ! JCL: add BACK as argument to subroutine to avoid some tests that
        !       would cause program to terminate when program is run backwards
        !         CALL CONDSK(NUMGRD,NUMTYP,JET,IFHR,CSUM,BACK)
        !***********************************************************************
        !
        IF(LAGSAM) CALL LAGOUT                                                &
             !dwen(20090817): HYSPLIT v4.9 supplemented LAGOUT
             (LAGS,JET,KPM,DT,MASS,XPOS,YPOS,ZPOS,HDWP,PAGE,PTYP,PGRD,ZMDL)

        !    zero-out array if required

        !***********************************************************************
        !dwen(20090306): HYSPLIT from v4.5 removed CONC from COMMON bloc and added to argument list
        !dwen(20090817): HYSPLIT from v4.9 used DT to determine backward or forward 
        !dwen(20090817): HYSPLIT from v4.5 removed NUMTYP, no need for pollutants loop anymore

        CALL CONZRO(CONC,NUMGRD,DT,JET,IFHR,CSUM)

        !         CALL CONZRO(NUMGRD,NUMTYP,DT,JET,IFHR,CSUM,BACK)
        !***********************************************************************


        !    diagnostics to message file
        WRITE(KF21,*)' NOTICE   main: ',KH,JET,KPM,TMASS
        !****************************************************************
        !dwen(20090309):from STILT
        ! JCL:(4/27/01) branch out of timestep loop and stop program if fraction of particles
        !           leaving model area exceeds specified fraction
        !         WRITE(45,*)KS,COUNTNPAROUT,IDINT(OUTFRAC*KPM),
        !     &                    MAX(IDINT(OUTFRAC*KPM),1)
        IF(COUNTNPAROUT > MAX(INT(OUTFRAC*KPM),1))GO TO 300
        !*************************************************************************
        !
        !-------------------------------------------------------------------------------
        ! optional code section for GEM integration 
        !-------------------------------------------------------------------------------

        ! terminate number of time steps per hour loop
     END DO

     !-------------------------------------------------------------------------------
     ! after end of time steps per hour loop call more infrequent routines 
     !-------------------------------------------------------------------------------

     ! special diagnostics dump
     IF(MOD(KH, 6).EQ.0.AND.TMASS.GT.0.0)THEN
        WRITE(KF21,'(55X,A6,A8,A6)')'Index','Height','%Mass'
        DO KZ=NLVL,1,-1
           PCNT=100.0*ZMASS(KZ)/TMASS
           IF(PCNT.GT.0.0)THEN
              HEIGHT=(1.0-ZSG(KZ))*ZMDL
              WRITE(KF21,'(55X,I6,F8.1,F6.2)')KZ,HEIGHT,PCNT
           END IF
        END DO
     END IF

     ! puff/particle conversion routines (call before merge)
     !-------------------------------------------------------------------------------

     ! optional conversion of puffs to particles
     IF(INITD.EQ.103.OR.INITD.EQ.104)                                             &
          !dwen(20090817): HYSPLIT v4.9 supplemented LAGOUT
          CALL PARPUF(KPM,PAGE,SIGU,SIGV,SIGW,HDWP,PGRD,ABS(CONAGE*60))
     IF(INITD.EQ.130.OR.INITD.EQ.140)                                             &
          !******************************************************************
          !dwen(20090309):use updated PUFPAR 
          !               when PAGE > CONAGE, puff is converted into particles in new version,
          !               whereas SIGH > PDIST(25000m) in old version  
          !dwen(20090815): HYSPLIT4.9 used SIGW, instead of SIGV, to stand for vertical velocity sigma, and 
          !                           SIGU and SIGV for horizontal velocity sigma 
          !               NUMPAR is an argument, but never used in old PUFPAR
          
          CALL PUFPAR(KPM,PAGE,SIGU,SIGV,SIGW,HDWP,PGRD,ABS(CONAGE*60))
     !      IF(PRESET)CALL PUFPAR(KPM,SIGH,SIGV,SIGX,HDWP,PGRD,NUMPAR)
     !******************************************************************

     !JCL: comment out the following PUFF splitting & management statements, since
     !     working mainly with PARTICLES
     !-------------------------------------------------------------------------------
     ! puff splitting routines called each hour (default KSPL=1)
     !-------------------------------------------------------------------------------

     !! skip routines for 3D particle simulations
     !  IF(INITD.GT.0.AND.MOD(KH,KSPL).EQ.0)THEN
     !!    vertical puff splitting
     !     CALL PUFSPV(KPM,NLVL,ZSG,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,                   &
     !                 HDWP,PAGE,PTYP,PGRD,NSORT,NUMSPL,KSPK)
     !
     !!    horizontal puff splitting
     !     CALL PUFSPH(KPM,GRID(1:,KTIME)%SIZE,MASS,XPOS,YPOS,ZPOS,SIGH,SIGW,        &
     !                 HDWP,PAGE,PTYP,PGRD,NSORT,NUMSPL,CGSIZE,SPLITF,KSPK)
     !
     !!    modify horizonal merge parameters when array limit is reached 
     !     IF(SPLITON)THEN
     !        IF(KSPK.GT.0)THEN
     !!          splitting is on and the array limit is reached, then 
     !           SPLITON=.FALSE.
     !        WRITE(KF21,*)' NOTICE   main: split stopped array limit reached'
     ! 
     !           IF(FRHS.LT.FRHMAX)THEN
     !!             increase horizontal merging first few times
     !              FRHS=MIN(FRHMAX,FRHS+0.5)
     !              FRHE=MIN(FRHMAX,FRHE+0.5)
     !              WRITE(KF21,*)' NOTICE   main: merge increased FRHS=',FRHS
     !
     !!             Horizontal splitting distance adjustment based upon how many
     !!             particles are required to cover the concentration grid (assume 
     !!             20% in PBL). If there are insufficient number of particles, 
     !!             then the puffs are allowed to grow larger before they split, 
     !!             providing for smoother patterns. SQRT converts area to distance. 
     !              IF(SPLITF.EQ.1.0) THEN
     !                 SPLITF=MAX(1.1,SQRT(5.0*NXP*NYP/MAXPAR))
     !                 WRITE(KF21,*)' NOTICE   main: split changed SPLITF=',SPLITF
     !              END IF
     !           END IF
     !        END IF
     !
     !     ELSE
     !        IF(KSPK.EQ.0)THEN
     !!          previously splitting was off but last call to split routine has
     !!          a normal exit which means all puffs meeting criteria split
     !           SPLITON=.TRUE.
     !           WRITE(KF21,*)' NOTICE   main: split turned on; FRHS=',FRHS
     !
     !        ELSEIF(KSPK.EQ.2)THEN
     !!          previously splitting was off but was briefly turned on and off
     !!          inside the last call to the splitting routines
     !           IF(FRHS.LT.FRHMAX)THEN
     !              WRITE(KF21,*)' NOTICE   main: split turned on/off; FRHS=',FRHS
     !!             increase horizontal merging first few times
     !              FRHS=MIN(FRHMAX,FRHS+0.5)
     !              FRHE=MIN(FRHMAX,FRHE+0.5)
     !           END IF
     !        END IF
     !     END IF
     !  END IF

     !-------------------------------------------------------------------------------
     ! puff management routines called each hour
     !-------------------------------------------------------------------------------

     !! skip routines for 3D particle simulations
     !  IF(INITD.GT.0.AND.KPM.GT.MAXPAR/10.AND.INITD.NE.109)THEN
     !!    eliminate unused puffs (pgrd=0)
     !     CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,   &
     !                 PGRD,NSORT,KHMAX,0.0)
     !
     !!    sort puffs by position before merge (fractions: h,v,t)
     !     CALL PUFSRT(FRHS,FRVS,FRTS,0.0,KPM,ZMDL,XPOS,YPOS,ZPOS,MASS,SIGH,         &
     !                 SIGW,HDWP,PTYP,PGRD,NSORT,GRID(1:,KTIME)%SIZE)
     !
     !!    merge puffs together if they are in the same position
     !     CALL PUFMRG(FRHS,FRVS,FRTS,0.0,KPM,ZMDL,MASS,XPOS,YPOS,ZPOS,              &
     !                 SIGH,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,GRID(1:,KTIME)%SIZE)
     !  END IF
     !
     !! eliminate unused puffs (pgrd=0)
     !  CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,      &
     !              PGRD,NSORT,KHMAX,0.0)
     !
     !!-------------------------------------------------------------------------------
     !! less frequent but less restictive merge using enhanced merge parameters
     !! and the option to revove low mass particles if mass revoval turned on  
     !!-------------------------------------------------------------------------------
     !
     !  IF(INITD.GT.0.AND.KPM.GT.MAXPAR/4.AND.MOD(KH,KRND).EQ.0.AND.INITD.NE.109)THEN
     !!    find the frme mass level for puff merging  - fmass
     !!    find the frmr mass level for puff deletion - rmass
     !     CALL PUFMIN(KPM,MASS,NSORT,TMASS,FRME,FMASS,FRMR,RMASS)
     !
     !!    sort puffs by position before merge (fractions: h,v,t)
     !     CALL PUFSRT(FRHE,FRVE,FRTE,FMASS,KPM,ZMDL,XPOS,YPOS,ZPOS,MASS,SIGH,       &
     !                 SIGW,HDWP,PTYP,PGRD,NSORT,GRID(1:,KTIME)%SIZE)
     !
     !!    merge puffs (mass<fmass) if they are in the same position
     !     CALL PUFMRG(FRHE,FRVE,FRTE,FMASS,KPM,ZMDL,MASS,XPOS,YPOS,ZPOS,            &
     !                 SIGH,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,GRID(1:,KTIME)%SIZE)
     !
     !!    eliminate unused puffs (pgrd=0) and those with mass <rmass
     !     CALL PUFDEL(KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,   &
     !                 PGRD,NSORT,KHMAX,RMASS)
     !  END IF

     ! percent completion
     PCNT=FLOAT(KH)*100.0/FLOAT(NHRS)

     WRITE(*,'(1X,A,F5.1)')'Percent complete: ',PCNT

     !-------------------------------------------------------------------------------
     ! particle input and output for model reinitialization
     !-------------------------------------------------------------------------------

     ! dump particle positions
     IF(NDUMP.GT.0.AND.KH.EQ.NDUMP)THEN
        !    open file once (5/7/2003-multiple times) 
        INQUIRE(FILE=POUTF,OPENED=FTEST2)
        !    for LINUX: CONVERT='BIG_ENDIAN'
        IF(.NOT.FTEST2)OPEN(KF24,FILE=POUTF,FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

        CALL PAROUT(JET,KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,    &
             SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,ZMDL)

        !    check for pardump file cycling
        IF(NCYCL.EQ.0)THEN
           !       no cycling then shutdown after first time
           NDUMP=0
        ELSE
           !       set next dump time according to cycle time
           NDUMP=NDUMP+NCYCL
        END IF
     END IF

     ! load more particles if time matches
     IF(FTEST1)THEN
        CALL PARINP(JET,KG,KT,KPM,MASS,TLAT,TLON,XPOS,YPOS,ZPOS,SIGH,SIGU,  &
             SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,ZMDL,MAXPAR,NINIT)
        INQUIRE(FILE=PINPF,OPENED=FTEST1)
     END IF

     !-------------------------------------------------------------------------------
     ! optional code section to update GEM meteorology and output         
     !-------------------------------------------------------------------------------

     !-------------------------------------------------------------------------------
     ! termination
     !-------------------------------------------------------------------------------

     ! terminate hours to run loop
  END DO

  !***********************************************************
  !dwen(20090309):following lines from STILT
  ! JCL:(4/27/01) branch here if fraction of particles which have moved off grid exceeds OUTFRAC
300 CONTINUE

  ! JCL:close 'JCLMESSAGE'
  CLOSE(KFJCLMSG)

  ! JCL:close 'PARTICLE.DAT'
  CLOSE(KFPARDAT)

  ! JCL:(11/03/03)
  IF(WINDERRTF.EQ.1)DEALLOCATE(UVERR)

  !***********************************************************

  ! special termination messages
  IF(.NOT.SPLITON)THEN
     WRITE(*,*)   'WARNING main: puff splitting turned off ... check MESSAGE'
     WRITE(KF21,*)'WARNING main: puff splitting turned off ... expand MAXPAR'
     WRITE(KF21,*)'        or increase merge parameters in CFG namelist file'
  END IF
  IF(KEMIT.EQ.1)THEN
     WRITE(*,*)   'WARNING particle array limit reached; emissions turned off'
     WRITE(KF21,*)'WARNING particle array limit reached; expand MAXPAR'
  END IF
  WRITE(*,*)'Complete Hysplit'

  DEALLOCATE (xpos,ypos,zpos,sigh,sigu,sigv,sigw,mass)
  DEALLOCATE (page,hdwp,ptyp,pgrd,nsort,spot,zsg,metz)  
  DEALLOCATE (zmass,dirt,dept,dryd,conc,csum)
  !****************************************************
  !dwen(20090309):
  deALLOCATE (zierr,wwprev,uverr,areapru,areaprd,icndx,      &
       sampttcum,footcum,zfxcum,convdur,dmasswt,      &
       uerr,verr,uerr2,verr2,uuerr_t,             &
       vverr_t,zprofm)
  !****************************************************

  IF(ALLOCATED(sprt))     DEALLOCATE (sprt)  
  IF(QFILE)               DEALLOCATE (polid,qarea)
  IF(ALLOCATED(tlat))     DEALLOCATE (tlat,tlon)
  IF(LAGSAM)              DEALLOCATE (lags) 

  CALL DALLOC

  ! required for NCEP operational implementation
  ! CALL W3TAGE('HYMODELC')
  ! STOP

END PROGRAM hymodelc
