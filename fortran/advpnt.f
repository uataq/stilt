!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ADVPNT           ADVection of one PoiNT in space
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ADVECTION OF ONE POINT IN SPACE IS THE PRIMARY ROUTINE THAT IS USED BY
!   THE TRAJECTORY AND DISPERSION SIMULATIONS.  IT IS CALLED EACH TIME STEP.
!   THE ROUTINE CHECKS METEO DATA IN ARRAY, IF POINT FITS WITHIN THE TIME AND
!   SPACE LIMITS CALCULATION CONTINUES, OTHERWISE NEW DATA ARE INPUT.
!   IN ADDITION,  METEO VARIABLES AT THE END OF THE STEP ARE PLACED IN THE
!   METO STRUCTURE.  THOSE VALUES ARE USED IN SUBSEQUENT DISPERSION AND
!   TRAJECTORY OUTPUT ROUTINES.  NOTE THAT ARRAY DIMENSIONS ARE AT THE
!   COMPILED MAXIMUM IN THIS UPPER LEVEL ROUTINE.  ALL ROUTINES CALL
!   HERE REQUIRE SUB-GRID DIMENSIONS, EXCEPT IN THE VERTICAL DIMENSION.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 10 Apr 1998 (RRD)
!                 18 Aug 1998 (RRD) - added isotropic turbulence option
!                 27 Oct 1998 (RRD) - Clean return (KG=0) on no met data
!                 21 Dec 1998 (RRD) - added endrec to metinp call
!                                     new subroutine to compute horiz mixing
!                                     added rdep to argument list
!                 25 Feb 1999 (RRD) - corrected ftime initialization
!                 04 Mar 1999 (RRD) - changed argument list to metpos
!                 19 Apr 1999 (RRD) - added terrain height array
!                 01 Mar 2000 (RRD) - fixed subgrid problem between grids
!                 22 Sep 2000 (RRD) - fortran90 upgrade
!                 16 Mar 2001 (RRD) - optional global lat lon grids
!                 21 Jun 2001 (RRD) - added temperature array
!                 04 Oct 2001 (RRD) - simultaneous multiple meteorology
!                 17 Oct 2001 (RRD) - added mixed layer depth as variable
!                 05 Dec 2001 (RRD) - option to read terrain file
!                 18 Dec 2001 (RRD) - correction to global cyclic BC
!                 26 Feb 2002 (RRD) - downward shortwave flux
!                 21 May 2002 (RRD) - velocity divergence option
!                 23 Jul 2002 (RRD) - eta terrain correction
!                 13 Aug 2002 (RRD) - limit subgrid to data grid size
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 13 Feb 2003 (RRD) - corrected problem with subgrid size
!                 26 Mar 2003 (RRD) - subgrid size requires grid # dimension
!                 10 Apr 2003 (RRD) - replaced DM with MTIME in related subs
!                 24 Jul 2003 (RRD) - moved meteo data array to module metval.f
!                 16 Sep 2003 (RRD) - surface advection option (chng arg list)
!                 14 Oct 2003 (RRD) - removed density as a metwnd argument
!                 17 Oct 2003 (RRD) - added turbulent kinetic energy 
!                 05 Nov 2003 (RRD) - reconfigured for turbulence 
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 01 Oct 2004 (RRD) - subgrid limits test after advection
!                 25 Oct 2004 (RRD) - pass through tratio argument
!                 08 Mar 2006 (RRD) - static stability parameter
!                 18 May 2006 (RRD) - grid dimension to forecast hour
!                 22 May 2006 (RRD) - patch test for complex HDWP/INITD 
!                                   - mixed layer depth to metinp
!                 21 Nov 2006 (RRD) - day night tke partition
!                 21 May 2007 (RRD) - preserve sign of DM to METWND
!                 04 Jun 2008 (RRD) - additional mixing depth options
!
! USAGE:  CALL ADVPNT(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,KSFC,TKERD,TKERN,
!                     HDWP,XP,YP,ZP,JET,DT,TRATIO,KG,KT,NGRD,NTIM,ZSG,NLVL,
!                     ZMDL,KVEL,KMIXD,KMIX0,UBAR,IFHR,KRET)
!     
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none 
!   OUTPUT FILES:           unit KF21 for diagnostic messages
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

!dwen(20090813) SUBROUTINE ADVPNT(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,KSFC,TKERD,TKERN,HDWP,   &
!dwen(20090813)                   XP,YP,ZP,JET,DT,TRATIO,KG,KT,NGRD,NTIM,ZSG,NLVL,ZMDL,KVEL,  &
!dwen(20090813)                   KMIXD,KMIX0,UBAR,IFHR,KRET)

!dwen(20090813) ****************
! JCL:   add DUDZ&DVDZ&WWOUT as output arguments in dummy call as well
! JCL:(9/16/02) add ZICONTROLTF, NHRSZI, ZIPRESC to prescribe mixed-layer height
! CHG:(12/05/01) add 'convdur' (here only dummy, use first particles convdur)
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! CHG(09/18/03) pass on RAMSFLG
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
!dwen(20090822): add conc in argument list.

 SUBROUTINE ADVPNT(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,KSFC,TKERD,TKERN,HDWP,   &
                   XP,YP,ZP,JET,DT,TRATIO,KG,KT,NGRD,NTIM,ZSG,NLVL,ZMDL,KVEL,  &
                   KMIXD,KMIX0,UBAR,IFHR,KRET,                                 &
                   conc,wwout,zicontroltf,nhrszi,zipresc,convdur1,iconvect,       &
                   ramsflg,ecmflg,xposprev,yposprev)

!dwen ******************************************
  USE funits
  USE metval
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE

  INCLUDE 'DEFARG2.INC' ! subroutine interfaces
  INCLUDE 'DEFARG3.INC' ! global routines subroutine interfaces
  INCLUDE 'DEFMETO.INC' ! meteo variables returned at advection point
!dwen(20090819) ***************
  INCLUDE 'DEFCONC.INC' ! meteo variables returned at advection point
!******************************

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(bset),INTENT(OUT)   :: metz (:)     ! profile advection variables
  TYPE(aset),INTENT(OUT)   :: meto         ! surface advection variables
  LOGICAL,   INTENT(IN)    :: back         ! integration direction flag
  LOGICAL,   INTENT(IN)    :: vmix         ! return mixing profile flag
  LOGICAL,   INTENT(IN)    :: cdep         ! return deposition variable flag
  LOGICAL,   INTENT(IN)    :: rdep         ! resistance deposition flag
  LOGICAL,   INTENT(IN)    :: traj         ! return trajectory variables flag
  INTEGER,   INTENT(IN)    :: ksfc         ! top of surface layer index
  REAL,      INTENT(IN)    :: tkerd        ! day turbulent kinetic eneregy ratio  
  REAL,      INTENT(IN)    :: tkern        ! night turbulent kinetic eneregy ratio  
  INTEGER,   INTENT(IN)    :: hdwp         ! puff/particle distribution type  
  REAL,      INTENT(INOUT) :: xp,yp,zp     ! particle position for advection
  INTEGER,   INTENT(IN)    :: jet          ! elapsed time (min)
  REAL,      INTENT(IN)    :: dt           ! advection time step (min)
  REAL,      INTENT(IN)    :: tratio       ! time step stability criterion
  INTEGER,   INTENT(INOUT) :: kg           ! grid index for calculation
  INTEGER,   INTENT(INOUT) :: kt           ! current grid time index      
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
 
  LOGICAL                  :: new1,new2    ! time read flags    
  LOGICAL                  :: offg         ! off-grid, sub-grid
  INTEGER                  :: mtime(2)     ! input data times requested 
  INTEGER,     ALLOCATABLE :: fhour(:,:)   ! forecast hour associated with data
  INTEGER,     ALLOCATABLE :: ftime(:,:)   ! time of meteo data in array    

  INTEGER,     ALLOCATABLE :: lx1(:)       ! subgrid corner point
  INTEGER,     ALLOCATABLE :: ly1(:)
  INTEGER,     ALLOCATABLE :: nxs(:)       ! subgrid size 
  INTEGER,     ALLOCATABLE :: nys(:)


  INTEGER                  :: nxt,nyt,nzs
  INTEGER                  :: kend1,kend2,kk,k1,k2,kgx,kt1,kt2
  INTEGER                  :: jtime,krec1,krec2,kunit1,kunit2,hdwpx
  REAL                     :: xx,yy,zz,dm
  CHARACTER(80)            :: ecode

!dwen(20090813) *********************
  real                     :: dens1,dens2
  logical                  :: awrfflg,fluxflg,deepflg,shallflg
  real,        allocatable :: tlk1(:),tlk2(:),sigwk1(:),sigwk2(:),zlvls(:)
  real                     :: ziscale,gridhr,var1,var2,temp1,temp2,rhfr1, &
                              tf,zk,rhfr2,zzz,zmlk1,zmlk2,zgrd,zmin,      &
                              zmixnew,zagl,zmix,znowsc,zza,               &
                              z1z,z1,frac,ztmp,ztmpr
  integer                  :: kgnow,kgold,kk1,kk2,kl,kz
  integer                  :: NZMAX ! (greater of NZS and NLVL, for allocation purposes)
!dwen ******************************
  integer :: dummy_kret
  real :: xx_min, yy_min, xx_max, yy_max
  LOGICAL                  :: dead    ! time read flags    

!-------------------------------------------------------------------------------
! external variable definitons
!-------------------------------------------------------------------------------


  SAVE LX1,LY1,NXS,NYS,MTIME,FTIME,FHOUR

  ! velocity offsets for staggered grids
  REAL       :: SUX,SUY,SVX,SVY
  COMMON /STAGGER/ SUX,SUY,SVX,SVY

!dwen(20090826)***********************

 interface

  SUBROUTINE ADV3NT(S,XP,YP,ZX,SS,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: xp,yp         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  LOGICAL,   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundary values
  END SUBROUTINE adv3nt

  SUBROUTINE ADV3NTZML(S,X1,Y1,ZX,SS,GLOBAL,NXP,NYP)

  IMPLICIT NONE

  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: x1,y1         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  LOGICAL,   INTENT(IN)    :: global        ! cyclic boundary condition flag
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries

  end subroutine adv3ntzml

  SUBROUTINE ADV3NTZLOC(S,X1,Y1,ZX,SS,GLOBAL,NXP,NYP)

  IMPLICIT NONE

  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: x1,y1         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  LOGICAL,   INTENT(IN)    :: global        ! cyclic boundary condition flag
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries

 end subroutine adv3ntzloc
!

 SUBROUTINE ADV2NT(S,X1,Y1,SS,GLOBAL,NXP,NYP)

  IMPLICIT NONE

  REAL,      INTENT(IN)    :: s(:,:)        ! field for interpolation
  REAL,      INTENT(IN)    :: x1,y1         ! position of interpolated value
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  LOGICAL,   INTENT(IN)    :: global        ! cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundaries  

 end subroutine adv2nt
!

 end interface
!**************************************


!-------------------------------------------------------------------------------
! Initialize array data time marker variable to indicate the status of the
! meteorological data for each grid at the pointer variable time.  
! =0 - indicates data have not yet been loaded into the array
! >0 - the internal time of the meteo data in array at the pointer location
!-------------------------------------------------------------------------------


  IF(.NOT.ALLOCATED(ftime))THEN

     ECODE='Meteo array time pointer'
     ALLOCATE(ftime(ngrd,2), STAT=kret)
     IF(kret.NE.0)GOTO 9000
     FTIME=0

     ECODE='Meteo forecast time pointer'
     ALLOCATE(fhour(2,ngrd), STAT=kret)
     IF(kret.NE.0)GOTO 9000
     FHOUR=0

     ECODE='Subgrid pointers'
     ALLOCATE(lx1(ngrd),ly1(ngrd),nxs(ngrd),nys(ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

!    The current subgrid size is saved in this routine with initial values 
!    that were set in advrng. Subsquent changes to the subgrid structure
!    variables represent potential new subgrid definitions if required.
!    This is determined in the call to metsub, when offg is true.

     DO KK=1,NGRD
        LX1(KK)=GRID(KK,1)%LX1
        LY1(KK)=GRID(KK,1)%LY1
        NXS(KK)=GRID(KK,1)%LXR
        NYS(KK)=GRID(KK,1)%LYR
     END DO
  END IF

!dwen(20090820) ********************
  IF(.NOT.ALLOCATED(tlk1))THEN

     ECODE='temporary variables'
     ALLOCATE(tlk1(nlvl),tlk2(nlvl),sigwk1(nlvl),sigwk2(nlvl),zlvls(nlvl),STAT=kret)
     IF(kret.NE.0)GOTO 9000
  end if
! **********************************

!
!-------------------------------------------------------------------------------
! For the input position and current time, determine the grid from which data
! should be loaded, the appropriate record number, and subgrid requirements.
! Calculation may terminate if point OFF all Grids and OFF in Time.
!-------------------------------------------------------------------------------

  CALL METPOS(BACK,XP,YP,JET,NGRD,NTIM,FTIME,MTIME,POINT,OFFG,KG,KGX,KT1,KT2)

!
!      CALL METPOS(BACK,XP,YP,JET,KG,NGRD,NXS,NYS,                       &
!     &   METO%LX1,METO%LY1,METO%LXC,METO%LYC,METO%LXR,METO%LYR,         &
!     &   KREC1,KREC2,OFFG,NEWX,MTIME,FTIME,KGRID)
!

   kgold = 0

! location not on any meteorological grid
  IF(OFFG)THEN
     kgold=kg
     KG=0       ! zero grid id flags particles for elimination
     KRET=1
     RETURN
  END IF

!****************************
!dwen(20090906)
!add the following four lines to keep values of xp and yp when grid switch occurs
if(kg.ne.kgx) then
  xposprev=xp
  yposprev=yp
end if
!************************************

! switched grids data have been remapped, mark particle on new grid 
  IF(KG.NE.KGX) KG=KGX
  IF(KT.NE.KT1) KT=KT1

! set meteo array index for last (1) and next (2) time
  K1=POINT(1)
  K2=POINT(2)

! unit numbers for data
  KUNIT1=FILE(KG,KT1)%KUNIT
  KUNIT2=FILE(KG,KT2)%KUNIT

! last valid record in each file
  KEND1=FILE(KG,KT1)%ENDREC
  KEND2=FILE(KG,KT2)%ENDREC

! each entry the default is that no data load required
  NEW1=.FALSE.
  NEW2=.FALSE.

! time mismatch between request and array load new data
  IF(MTIME(1).NE.FTIME(KG,K1)) NEW1=.TRUE.
  IF(MTIME(2).NE.FTIME(KG,K2)) NEW2=.TRUE.

! compute record numbers to index record for indicated time 
  KREC1=DREC(KG,KT1)%REC_PER*(MTIME(1)-FILE(KG,KT1)%FIRST%MACC)  &
       /DREC(KG,KT1)%DELTA+1
  KREC2=DREC(KG,KT2)%REC_PER*(MTIME(2)-FILE(KG,KT2)%FIRST%MACC)  &
       /DREC(KG,KT2)%DELTA+1

!-------------------------------------------------------------------------------
! Check particle position on the subgrid.  When it is OFFGrid then redefine
! the subgrid domain and set the flags to load new meteorological data fields.
!-------------------------------------------------------------------------------

! modified argument (02/07/2003 & 03/26/2003) 
  CALL METSUB(XP,YP,KG,KT1,OFFG,LX1(KG),LY1(KG),NXS(KG),NYS(KG),dummy_KRET)
!
!  no subroutine METSUB in old stilt
!

  IF(OFFG)THEN
!    any change in subgrid requires data reload at both times
     NEW1=.TRUE.
     NEW2=.TRUE.
  END IF

!-------------------------------------------------------------------------------
! Determine the number of input data levels NZS that will be processed
! according to the model-top limit set in the control file.  Input data
! are read into the same array that is used for the internal vertical   
! coordinate system with NLVL grid points.
!-------------------------------------------------------------------------------

  IF(NEW1.OR.NEW2)THEN
!    nzs: input data   nlvl: STILT internal array   grid: data read
     NZS=GRID(KG,KT)%NZ-1
     NZMAX=MAX(NZS,NLVL)
  END IF

!dwen(20090813) *******************
awrfflg = GRID(kg,kt)%model_id(2:4) .eq. 'WRF'
!dwen *****************************


!-------------------------------------------------------------------------------
! determine if new subgrid allocation is required (use 'a' for testing)
!-------------------------------------------------------------------------------

! check if array size has changed
  IF(.NOT.ALLOCATED(a).OR.(SIZE(a,1).LT.NXS(KG).OR.SIZE(a,2).LT.NYS(KG)))THEN

!    subgrid size maintained as maximum of all sub-grids
     NXT=MAXVAL(NXS)
     NYT=MAXVAL(NYS)

     IF(ALLOCATED(a))THEN
!       if already allocated deallocate before redefining
        ECODE='3d variable deallocation'
!dwen(20090813)        DEALLOCATE( u,v,w,a,t,q,e,p,x,h,STAT=kret)
        DEALLOCATE( u,v,w,a,t,q,e,p,x,h,                &
             d,tlrams,sigwrams,cfxup1,cfxup2,cfxdn1,    &
             dfxup1,dfxup2,dfxdn1,efxup1,efxup2,efxdn1, &
             tke,tl,sigw,dmass,xm,hm, STAT=kret)
        IF(kret.NE.0)GOTO 9000

        ECODE='2d variable deallocation'
!dwen(20090813)        DEALLOCATE( p0,rt,u0,v0,t0,zi,STAT=kret)
        DEALLOCATE( p0,rt,u0,v0,t0,zi,                  &
             rc,tc,sw,lc,sm,w0,muu,muv,mu,msfu,    &
             msfv,msft,zloc, STAT=kret)
        IF(kret.NE.0)GOTO 9000

        ECODE='Flux variable deallocation'
        DEALLOCATE( ds,ss,uf,vf,hf,lf,sf,STAT=kret)
        IF(kret.NE.0)GOTO 9000

        ECODE='Fixed variable deallocation'
        DEALLOCATE( gx,gy,z0,zt,lu,STAT=kret)
        IF(kret.NE.0)GOTO 9000
     END IF

!    all meteo variables defined by (x,y,z,t), where t varies from
!    1 (last observation time) to 2 (next observation time)

     ECODE='Velocity variables'
     ALLOCATE(u(nxt,nyt,NZMAX,2,ngrd),v(nxt,nyt,NZMAX,2,ngrd), &
              w(nxt,nyt,NZMAX,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

!dwen(20090813) *****************
     ECODE='stilt 3d variables'
     ALLOCATE(d(nxt,nyt,NZMAX,2,ngrd),tlrams(nxt,nyt,NZMAX,2,ngrd), &
              sigwrams(nxt,nyt,NZMAX,2,ngrd),cfxup1(nxt,nyt,NZMAX,2,ngrd),&
              cfxup2(nxt,nyt,NZMAX,2,ngrd),cfxdn1(nxt,nyt,NZMAX,2,ngrd),&
              dfxup1(nxt,nyt,NZMAX,2,ngrd),dfxup2(nxt,nyt,NZMAX,2,ngrd),&
              dfxdn1(nxt,nyt,NZMAX,2,ngrd),efxup1(nxt,nyt,NZMAX,2,ngrd),&
              efxup2(nxt,nyt,NZMAX,2,ngrd),efxdn1(nxt,nyt,NZMAX,2,ngrd),&
              tke(nxt,nyt,NZMAX,2,ngrd),tl(nxt,nyt,NZMAX,2,ngrd),&
              sigw(nxt,nyt,NZMAX,2,ngrd),dmass(nxt,nyt,NZMAX,2,ngrd),&
              xm(nxt,nyt,NZMAX,2,ngrd),hm(nxt,nyt,NZMAX,2,ngrd),&
              STAT=kret)
     IF(kret.NE.0)GOTO 9000
!dwen **************************


     ECODE='Temperature variables'
     ALLOCATE(t(nxt,nyt,NZMAX,2,ngrd),a(nxt,nyt,NZMAX,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     ECODE='Moisture & TKE variables'
     ALLOCATE(q(nxt,nyt,NZMAX,2,ngrd),e(nxt,nyt,NZMAX,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     ECODE='Diagnostic variables'
     ALLOCATE(p(nxt,nyt,NZMAX,2,ngrd),x(nxt,nyt,NZMAX,2,ngrd), &
              h(nxt,nyt,NZMAX,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     ECODE='Surface variables'
     ALLOCATE(p0(nxt,nyt,2,ngrd),rt(nxt,nyt,2,ngrd),u0(nxt,nyt,2,ngrd), &
              zi(nxt,nyt,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

!(20090813) *********************
     ECODE='stilt 2d variables'
     ALLOCATE(rc(nxt,nyt,2,ngrd),tc(nxt,nyt,2,ngrd),sw(nxt,nyt,2,ngrd), &
              lc(nxt,nyt,2,ngrd),sm(nxt,nyt,2,ngrd),w0(nxt,nyt,2,ngrd), &
              muu(nxt,nyt,2,ngrd),muv(nxt,nyt,2,ngrd), &
              mu(nxt,nyt,2,ngrd),msfu(nxt,nyt,2,ngrd),msfv(nxt,nyt,2,ngrd), &
              msft(nxt,nyt,2,ngrd),lf(nxt,nyt,2,ngrd), &
              zloc(nxt,nyt,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000
!dwen ***************************

     ECODE='Flux variables'
     ALLOCATE(uf(nxt,nyt,2,ngrd),vf(nxt,nyt,2,ngrd),hf(nxt,nyt,2,ngrd), &
              sf(nxt,nyt,2,ngrd),ss(nxt,nyt,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     ECODE='Chemistry variables'
     ALLOCATE(ds(nxt,nyt,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     ECODE='Shelter variables'
     ALLOCATE(v0(nxt,nyt,2,ngrd),t0(nxt,nyt,2,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     ECODE='Grid size variables'
     ALLOCATE(gx(nxt,nyt,ngrd),gy(nxt,nyt,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     ECODE='Fixed variables'
     ALLOCATE(z0(nxt,nyt,ngrd),zt(nxt,nyt,ngrd),lu(nxt,nyt,ngrd),STAT=kret)
     IF(kret.NE.0)GOTO 9000

     WRITE(KF21,*)' NOTICE advpnt: (kg,nx,ny,nz) - ',kg,nxt,nyt,NZMAX

!    Reset time of data in array so that subgrid data from all meteo
!    grids are reload.  Reallocation destroys existing contents.
     FTIME=0
  END IF

!-------------------------------------------------------------------------------
! Change in subgrid requires reload of surface boundary data 
!-------------------------------------------------------------------------------


  IF(NEW1.AND.NEW2)  &
!dwen(2009/07/30) updated the call of METGRD according to HYSPLIT v4.9
!dwen             KT1 and ZT were added in new version of METGRID

     CALL METGRD(KG,KT1,LX1(KG),LY1(KG),NXS(KG),NYS(KG),                       &
                 GX(:,:,KG),GY(:,:,KG),Z0(:,:,KG),LU(:,:,KG),ZT(:,:,KG))
!
!  CALL METGRD(KGRID,METO%LX1,METO%LY1,NXS,NYS,GX,GY,Z0,LU)
!

!-------------------------------------------------------------------------------
! test if new data required at the last time (k1)
!-------------------------------------------------------------------------------

  IF(NEW1)THEN

!    load meteorological data according to positioning
!    reads data for NZS levels
!************************************************
!dwen     CALL METINP(BACK,KG,KT1,KUNIT1,KREC1,LX1(KG),LY1(KG),NXS(KG),NYS(KG),     &
!dwen        NZS,FTIME(KG,K1),KEND1,FHOUR(K1,KG),ZT(:,:,KG),DS(:,:,K1,KG),          &
!dwen        P0(:,:,K1,KG),T0(:,:,K1,KG),U0(:,:,K1,KG),V0(:,:,K1,KG),               &
!dwen        UF(:,:,K1,KG),VF(:,:,K1,KG),HF(:,:,K1,KG),RT(:,:,K1,KG),ZI(:,:,K1,KG), &
!dwen        U(:,:,:,K1,KG),V(:,:,:,K1,KG),W(:,:,:,K1,KG),A(:,:,:,K1,KG),           &
!dwen        Q(:,:,:,K1,KG),P(:,:,:,K1,KG),E(:,:,:,K1,KG),H(:,:,:,K1,KG),           &
!dwen        X(:,:,:,K1,KG))
!dwen(20090323):try to use updated verison of METINP according to HYSPLIT4.9
!dwen(20090731):   note that SF (sensible heat flux) was changed into HF after HYSPLIT v4.6,
!dwen(20090731):      and SF represents stability function after HYSPLIT v4.6
!dwen(20090730):   added following arguments according STILT codes
!dwen(20090730):   hpbl
!dwen(20090730)    w0                             surface level w wind
!dwen(20090730)    muu,muv                        momentum fluxes of u,v
!dwen(20090730)    mu
!dwen(20090730)    msfu,mufv,msft
!dwen(20090730)    fluxflg,deepflg,shallflg      flags
!dwen(20090730)    A                             upper level air temperature 
!dwen(20090730)    T                             upper level air temperature
!dwen(20090730)    E                             turbulent kinetic energy
!dwen(20090730)    H                             velocity variance U'2 and V'2
!dwen(20090730)    X                             velocity variance W'2
!dwen(20090730)    rc                            convective precipitation
!dwen(20090730)    lf                            latent heat flux
!dwen(20090730)    tc                            total cloud cover
!dwen(20090730)    lc                            low cloud cover
!dwen(20090730)    sw                            shortwave radiative flux
!dwen(20090730)    sm                            soil moisture
!dwen(20090730)    tlrams                        turbulence variables directly from RAMS
!dwen(20090730)    sigwrams                      sigw directly from RAMS
!dwen(20090730)    cfxup1,cfxdn1                 convective flux directly from rams, emc, wrf
!dwen(20090730)    cfxup2,
!dwen(20090730)    dfxup1,dfxdn1
!dwen(20090730)    dfxup2
!dwen(20090730)    efxup1,efxdn1
!dwen(20090730)    efxup2
!dwen(20090730)    tke                          turbulent kinetic energy
                   
                    
!dwen     CALL METINP(BACK,KG,KT1,KUNIT1,KREC1,LX1(KG),LY1(KG),NXS(KG),NYS(KG),     &
!dwen        NZS,FTIME(KG,K1),KEND1,FHOUR(K1,KG),ZT(:,:,KG),DS(:,:,K1,KG),          &
!dwen        P0(:,:,K1,KG),T0(:,:,K1,KG),U0(:,:,K1,KG),V0(:,:,K1,KG),               &
!dwen        UF(:,:,K1,KG),VF(:,:,K1,KG),HF(:,:,K1,KG),RT(:,:,K1,KG),ZI(:,:,K1,KG), &
!dwen        U(:,:,:,K1,KG),V(:,:,:,K1,KG),W(:,:,:,K1,KG),A(:,:,:,K1,KG),           &
!dwen        Q(:,:,:,K1,KG),P(:,:,:,K1,KG),E(:,:,:,K1,KG),H(:,:,:,K1,KG),           &
!dwen        X(:,:,:,K1,KG))

      CALL METINP(BACK,KG,KT1,KUNIT1,KREC1,LX1(KG),LY1(KG),NXS(KG),NYS(KG),       &
           NZS,FTIME(KG,K1),KEND1,FHOUR(K1,KG),ZT(:,:,KG),DS(:,:,K1,KG),          &
           P0(:,:,K1,KG),T0(:,:,K1,KG),U0(:,:,K1,KG),V0(:,:,K1,KG),               &
           UF(:,:,K1,KG),VF(:,:,K1,KG),HF(:,:,K1,KG),RT(:,:,K1,KG),ZI(:,:,K1,KG), &
           U(:,:,:,K1,KG),V(:,:,:,K1,KG),W(:,:,:,K1,KG),A(:,:,:,K1,KG),           &
           Q(:,:,:,K1,KG),P(:,:,:,K1,KG),E(:,:,:,K1,KG),H(:,:,:,K1,KG),           &
           X(:,:,:,K1,KG),                                                        &
           W0(:,:,K1,KG),fluxflg, deepflg, shallflg,              &
           muu(:,:,K1,KG),muv(:,:,K1,KG),mu(:,:,K1,KG),                         &
           msfu(:,:,K1,KG),msfv(:,:,K1,KG),msft(:,:,K1,KG),                     &
           RC(:,:,K1,KG),LF(:,:,K1,KG),            &
           TC(:,:,K1,KG),LC(:,:,K1,KG),SW(:,:,K1,KG),SM(:,:,K1,KG),             &
           TLRAMS(:,:,:,K1,KG),SIGWRAMS(:,:,:,K1,KG),CFXUP1(:,:,:,K1,KG),       &
           CFXUP2(:,:,:,K1,KG),CFXDN1(:,:,:,K1,KG),DFXUP1(:,:,:,K1,KG),         &
           DFXUP2(:,:,:,K1,KG),EFXUP1(:,:,:,K1,KG),EFXUP2(:,:,:,K1,KG),         &
           DFXDN1(:,:,:,K1,KG),EFXDN1(:,:,:,K1,KG),TKE(:,:,:,K1,KG))
!
! CHG:(11/20/01)added conv. precip. rates (RC) to arguments
! CHG:(11/20/01)add tot. cloud cover and shortw. flux
! CHG:(11/20/01) add sfc flux for water (w/m2)
! CHG:(12/04/01) add low cloud cover
! CHG:(22/01/03) add soil moisture
! JCL(03/27/03): add arrays to store grids of TL & SIGW from RAMS
! CHG(09/23/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
!         CALL METINP(BACK,KGRID,KUNIT1,KREC1,NXY,METO%LX1,METO%LY1,     &
!     &      NXS,NYS,NZS,FTIME(1),KEND1,FHOUR(K1),ZT,                    &
!     &      P0(:,:,K1),T0(:,:,K1),U0(:,:,K1),V0(:,:,K1),W0(:,:,K1),hpbl(:,:,K1),     &
!     &      muu(:,:,K1),muv(:,:,K1),mu(:,:,K1),                         &
!     &      msfu(:,:,K1),msfv(:,:,K1),msft(:,:,K1),fluxflg, deepflg, shallflg,             &
!     &      UF(:,:,K1),VF(:,:,K1),SF(:,:,K1),RT(:,:,K1),                &
!     &      U(:,:,:,K1),V(:,:,:,K1),W(:,:,:,K1),T(:,:,:,K1),            &
!     &      Q(:,:,:,K1),P(:,:,:,K1),RC(:,:,K1),LF(:,:,K1),              &
!     &      TC(:,:,K1),LC(:,:,K1),SW(:,:,K1),SM(:,:,K1),                &
!     &      TLRAMS(:,:,:,K1),SIGWRAMS(:,:,:,K1),CFXUP1(:,:,:,K1),       &
!     &      CFXUP2(:,:,:,K1),CFXDN1(:,:,:,K1),DFXUP1(:,:,:,K1),         &
!     &      DFXUP2(:,:,:,K1),EFXUP1(:,:,:,K1),EFXUP2(:,:,:,K1),         &
!     &      DFXDN1(:,:,:,K1),EFXDN1(:,:,:,K1),TKEN(:,:,:,K1))
!!        WRITE(45,*)'advpnt, K1,K2:',K1,K2
         ZISCALE=1.0
! JCL:(9/16/02) if PBL height is prescribed
         IF(ZICONTROLTF.EQ.1)THEN
            if (nhrszi .eq. 1 .and. zipresc(nhrszi) .lt. 0.0) then
! Special case to denote using hpbl from met file
               ziscale = -1.
            else
! JCL:      time of loaded met grid (# of hrs elapsed since model starting time)
               GRIDHR=ABS((MTIME(1)-(CONC(1)%START%MACC))/60)
!dwen(20090825)               IF((GRIDHR+1).LE.NHRSZI)ZISCALE=ZIPRESC(GRIDHR+1)
               IF((GRIDHR+1).LE.NHRSZI)ZISCALE=ZIPRESC(int(GRIDHR+1))
            endif
         END IF
!************************************************
!    vertical interpolation to terrain following coordinates
!    data from NZS levels processed to NLVL of internal grid
!***********************************************************
!dwen(20090323):try to use updated PRFCOM
!               the call of stbhor has been moved in PRFCOM since HYSPLIT4.7
!               HF: sensible heat flux
!               SF(stability function): added after HYSPLIT v4.6
!               since HYSPLIT4 version of PRFCOM has an argument called X (
!               vertical velocity variance ), The argument X (vertical mixing 
!               coefficient ) in STILT was renamed to XM
!dwen(20090806)   add argument HM to output horizontal mixing coefficient
!dwen     CALL PRFCOM(TKERD,TKERN,KG,KT1,KSFC,GX(:,:,KG),GY(:,:,KG),Z0(:,:,KG),     &
!dwen        ZT(:,:,KG),NXS(KG),NYS(KG),NZS,ZMDL,ZSG,NLVL,VMIX,KMIXD,KMIX0,         &
!dwen        ZI(:,:,K1,KG),                                                         &
!dwen        P0(:,:,K1,KG),T0(:,:,K1,KG),U0(:,:,K1,KG),V0(:,:,K1,KG),UF(:,:,K1,KG), &
!dwen        VF(:,:,K1,KG),HF(:,:,K1,KG),SF(:,:,K1,KG),SS(:,:,K1,KG),U(:,:,:,K1,KG),&
!dwen        V(:,:,:,K1,KG),W(:,:,:,K1,KG),A(:,:,:,K1,KG),T(:,:,:,K1,KG),           &
!dwen        Q(:,:,:,K1,KG),P(:,:,:,K1,KG),E(:,:,:,K1,KG),H(:,:,:,K1,KG),           &
!dwen        X(:,:,:,K1,KG))


     CALL PRFCOM(TKERD,TKERN,KG,KT1,KSFC,GX(:,:,KG),GY(:,:,KG),Z0(:,:,KG),     &
        ZT(:,:,KG),NXS(KG),NYS(KG),NZS,ZMDL,ZSG,NLVL,VMIX,KMIXD,KMIX0,         &
        ZI(:,:,K1,KG),                                                         &
        P0(:,:,K1,KG),T0(:,:,K1,KG),U0(:,:,K1,KG),V0(:,:,K1,KG),UF(:,:,K1,KG), &
        VF(:,:,K1,KG),HF(:,:,K1,KG),SF(:,:,K1,KG),SS(:,:,K1,KG),U(:,:,:,K1,KG),&
        V(:,:,:,K1,KG),W(:,:,:,K1,KG),A(:,:,:,K1,KG),T(:,:,:,K1,KG),           &
        Q(:,:,:,K1,KG),P(:,:,:,K1,KG),E(:,:,:,K1,KG),H(:,:,:,K1,KG),           &
        X(:,:,:,K1,KG),                                                        &
        iconvect,w0(:,:,K1,KG),tl(:,:,:,K1,KG),d(:,:,:,k1,kg), &
        sigw(:,:,:,K1,KG),zloc(:,:,K1,KG),dmass(:,:,:,K1,KG),   &
        ziscale,tlrams(:,:,:,K1,KG),mu(:,:,K1,KG),muu(:,:,K1,KG),              &
        muv(:,:,K1,KG),msfu(:,:,K1,KG),msfv(:,:,k1,kg),msft(:,:,k1,kg),        &
        fluxflg, deepflg, shallflg, cfxup1(:,:,:,k1,kg),cfxup2(:,:,:,k1,kg),   &
        cfxdn1(:,:,:,k1,kg),dfxup1(:,:,:,k1,kg),dfxup2(:,:,:,k1,kg),           &
        efxup1(:,:,:,k1,kg),efxup2(:,:,:,k1,kg),dfxdn1(:,:,:,k1,kg),           &
        efxdn1(:,:,:,k1,kg),tke(:,:,:,k1,kg),xm(:,:,:,k1,kg),hm(:,:,:,k1,kg))


! JCL:   add 'TL' & 'SIGW' as arguments to PRFCOM
! JCL:(4/28/00)add ZML (mixed-layer ht) as output from PRFCOM
! CHG:(12/04/04)add ZLOC (lim. of convection ht) as output from PRFCOM
! JCL:(4/3/02)add mass violation grid as output from PRFCOM
! JCL:(9/16/02) ZISCALE is scaling factor for mixed-layer height
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
!        CALL PRFCOM(ISOT,KGRID,KSFC,GD,Z0,ZT,
!         CALL PRFCOM(ISOT,KGRID,KSFC,GX,GY,Z0,ZT,                       &
!     &      NXS,NYS,NZS,ZMDL,ZSG,NLVL,VMIX,ICONVECT,                    &
!     &      P0(:,:,K1),T0(:,:,K1),U0(:,:,K1),V0(:,:,K1),W0(:,:,K1),hpbl(:,:,K1),UF(:,:,K1),     &
!     &      VF(:,:,K1),SF(:,:,K1),U(:,:,:,K1),V(:,:,:,K1),              &
!     &      W(:,:,:,K1),T(:,:,:,K1),Q(:,:,:,K1),                        &
!     &      P(:,:,:,K1),D(:,:,:,K1),X(:,:,:,K1),TL(:,:,:,K1),           &
!     &      SIGW(:,:,:,K1),ZML(:,:,K1),ZLOC(:,:,K1),DMASS(:,:,:,K1),    &
!     &      ZISCALE,TLRAMS(:,:,:,K1),mu(:,:,K1),                        &
!     &      muu(:,:,K1),muv(:,:,K1),msfu(:,:,K1),msfv(:,:,K1),msft(:,:,K1), &
!     &      fluxflg, deepflg, shallflg, &
!     &      CFXUP1(:,:,:,K1),       &
!     &      CFXUP2(:,:,:,K1),CFXDN1(:,:,:,K1),DFXUP1(:,:,:,K1),         &
!     &      DFXUP2(:,:,:,K1),EFXUP1(:,:,:,K1),EFXUP2(:,:,:,K1),         &
!     &      DFXDN1(:,:,:,K1),EFXDN1(:,:,:,K1),TKEN(:,:,:,K1))
!
!************************************************************
!    reset meteo time position if data are missing
     IF(FTIME(KG,K1).NE.MTIME(1))THEN
     WRITE(KF21,*)'WARNING advpnt: Time 1 input does not match request!'
        WRITE(KF21,*)' Internal time: ',JET
        WRITE(KF21,*)'  Time of data: ',FTIME(KG,K1)
        WRITE(KF21,*)'  Request time: ',MTIME(1)
        MTIME(1)=FTIME(KG,K1)
     END IF
  END IF

!---------------------------------------------------------------------------
! test if new data required at the next time (k2)
!--------------------------------------------------------------------------

  IF(NEW2)THEN

!    load data for subgrid
!*******************************************************
!dwen(20090323):try to use updated METINP
!dwen(20090806)     CALL METINP(BACK,KG,KT2,KUNIT2,KREC2,LX1(KG),LY1(KG),NXS(KG),NYS(KG),     &
!dwen(20090806)        NZS,FTIME(KG,K2),KEND2,FHOUR(K2,KG),ZT(:,:,KG),DS(:,:,K2,KG),          &
!dwen(20090806)        P0(:,:,K2,KG),T0(:,:,K2,KG),U0(:,:,K2,KG),V0(:,:,K2,KG),               &
!dwen(20090806)        UF(:,:,K2,KG),VF(:,:,K2,KG),HF(:,:,K2,KG),RT(:,:,K2,KG),ZI(:,:,K2,KG), &
!dwen(20090806)        U(:,:,:,K2,KG),V(:,:,:,K2,KG),W(:,:,:,K2,KG),A(:,:,:,K2,KG),           &
!dwen(20090806)        Q(:,:,:,K2,KG),P(:,:,:,K2,KG),E(:,:,:,K2,KG),H(:,:,:,K2,KG),           &
!dwen(20090806)        X(:,:,:,K2,KG))

! CHG:(11/20/01)added conv. precip. rates (RC) to arguments
! CHG:(11/20/01)add tot. cloud cover and shortw. flux
! CHG:(11/20/01) add sfc flux for water (w/m2)
! CHG:(12/04/01) add low cloud cover
! CHG:(22/01/03) add soil moisture
! JCL(03/27/03): add arrays to store grids of TL & SIGW from RAMS
! CHG(09/23/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
      
      CALL METINP(BACK,KG,KT2,KUNIT2,KREC2,LX1(KG),LY1(KG),NXS(KG),NYS(KG),       &
           NZS,FTIME(KG,K2),KEND2,FHOUR(K2,KG),ZT(:,:,KG),DS(:,:,K2,KG),          &
           P0(:,:,K2,KG),T0(:,:,K2,KG),U0(:,:,K2,KG),V0(:,:,K2,KG),               &
           UF(:,:,K2,KG),VF(:,:,K2,KG),HF(:,:,K2,KG),RT(:,:,K2,KG),ZI(:,:,K2,KG), &
           U(:,:,:,K2,KG),V(:,:,:,K2,KG),W(:,:,:,K2,KG),A(:,:,:,K2,KG),           &
           Q(:,:,:,K2,KG),P(:,:,:,K2,KG),E(:,:,:,K2,KG),H(:,:,:,K2,KG),           &
           X(:,:,:,K2,KG),                                                        &
           W0(:,:,K2,KG),fluxflg, deepflg, shallflg,              &
           muu(:,:,K2,KG),muv(:,:,K2,KG),mu(:,:,K2,KG),                         &
           msfu(:,:,K2,KG),msfv(:,:,K2,KG),msft(:,:,K2,KG),                     &
           RC(:,:,K2,KG),LF(:,:,K2,KG),            &
           TC(:,:,K2,KG),LC(:,:,K2,KG),SW(:,:,K2,KG),SM(:,:,K2,KG),             &
           TLRAMS(:,:,:,K2,KG),SIGWRAMS(:,:,:,K2,KG),CFXUP1(:,:,:,K2,KG),       &
           CFXUP2(:,:,:,K2,KG),CFXDN1(:,:,:,K2,KG),DFXUP1(:,:,:,K2,KG),         &
           DFXUP2(:,:,:,K2,KG),EFXUP1(:,:,:,K2,KG),EFXUP2(:,:,:,K2,KG),         &
           DFXDN1(:,:,:,K2,KG),EFXDN1(:,:,:,K2,KG),TKE(:,:,:,K2,KG))
!

!dwen(20090806)         CALL METINP(BACK,KGRID,KUNIT2,KREC2,NXY,METO%LX1,METO%LY1,     &
!dwen(20090806)     &      NXS,NYS,NZS,FTIME(2),KEND2,FHOUR(K2),ZT,                    &
!dwen(20090806)     &      P0(:,:,K2),T0(:,:,K2),U0(:,:,K2),V0(:,:,K2),W0(:,:,K2),hpbl(:,:,K2),                &
!dwen(20090806)     &      muu(:,:,K2),muv(:,:,K2),mu(:,:,K2),                         &
!dwen(20090806)     &      msfu(:,:,K2),msfv(:,:,K2),msft(:,:,K2),fluxflg, deepflg, shallflg,             &
!dwen(20090806)     &      UF(:,:,K2),VF(:,:,K2),SF(:,:,K2),RT(:,:,K2),                &
!dwen(20090806)     &      U(:,:,:,K2),V(:,:,:,K2),W(:,:,:,K2),T(:,:,:,K2),            &
!dwen(20090806)     &      Q(:,:,:,K2),P(:,:,:,K2),RC(:,:,K2),LF(:,:,K2),              &
!dwen(20090806)     &      TC(:,:,K2),LC(:,:,K2),SW(:,:,K2),SM(:,:,K2),                &
!dwen(20090806)     &      TLRAMS(:,:,:,K2),SIGWRAMS(:,:,:,K2),CFXUP1(:,:,:,K2),       &
!dwen(20090806)     &      CFXUP2(:,:,:,K2),CFXDN1(:,:,:,K2),DFXUP1(:,:,:,K2),         &
!dwen(20090806)     &      DFXUP2(:,:,:,K2),EFXUP1(:,:,:,K2),EFXUP2(:,:,:,K2),         &
!dwen(20090806)     &      DFXDN1(:,:,:,K2),EFXDN1(:,:,:,K2),TKEN(:,:,:,K2))

         ZISCALE=1.0
! JCL:(9/16/02) if PBL height is prescribed
         IF(ZICONTROLTF.EQ.1)THEN
            if (nhrszi .eq. 1 .and. zipresc(nhrszi) .lt. 0.0) then
! Special case to denote using hpbl from met file
               ziscale = -1.
            else
! JCL:      time of loaded met grid (# of hrs elapsed since model starting time)
               GRIDHR=ABS((MTIME(2)-(CONC(1)%START%MACC))/60)
!dwen(20090825)               IF((GRIDHR+1).LE.NHRSZI)ZISCALE=ZIPRESC(GRIDHR+1)
               IF((GRIDHR+1).LE.NHRSZI)ZISCALE=ZIPRESC(int(GRIDHR+1))
            endif
         END IF
!**********************************************************

!    vertical interpolation of profile
!*********************************************************
!dwen(20090323):try to use updated PRFCOM
!               the call of stbhor has been moved in PRFCOM since HYSPLIT4.7
!dwen(20090806)     CALL PRFCOM(TKERD,TKERN,KG,KT2,KSFC,GX(:,:,KG),GY(:,:,KG),Z0(:,:,KG),     &
!dwen(20090806)        ZT(:,:,KG),NXS(KG),NYS(KG),NZS,ZMDL,ZSG,NLVL,VMIX,KMIXD,KMIX0,         &
!dwen(20090806)        ZI(:,:,K2,KG),       &
!dwen(20090806)        P0(:,:,K2,KG),T0(:,:,K2,KG),U0(:,:,K2,KG),V0(:,:,K2,KG),UF(:,:,K2,KG), &
!dwen(20090806)        VF(:,:,K2,KG),HF(:,:,K2,KG),SF(:,:,K2,KG),SS(:,:,K2,KG),U(:,:,:,K2,KG),&
!dwen(20090806)        V(:,:,:,K2,KG),W(:,:,:,K2,KG),A(:,:,:,K2,KG),T(:,:,:,K2,KG),           &
!dwen(20090806)        Q(:,:,:,K2,KG),P(:,:,:,K2,KG),E(:,:,:,K2,KG),H(:,:,:,K2,KG),           &
!dwen(20090806)        X(:,:,:,K2,KG))

! JCL:   add 'TL' & 'SIGW' as arguments to PRFCOM
! JCL:(4/28/00)add ZML (mixed-layer ht) as output from PRFCOM
!        vertical interpolation of profile
! CHG:(12/04/02)add ZLOC (lim. of convection ht) as output from PRFCOM
! JCL:(4/3/02)add mass violation grid as output from PRFCOM
! JCL:(9/16/02) ZISCALE is scaling factor for mixed-layer height
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)

     CALL PRFCOM(TKERD,TKERN,KG,KT2,KSFC,GX(:,:,KG),GY(:,:,KG),Z0(:,:,KG),     &
        ZT(:,:,KG),NXS(KG),NYS(KG),NZS,ZMDL,ZSG,NLVL,VMIX,KMIXD,KMIX0,         &
        ZI(:,:,K2,KG),                                                         &
        P0(:,:,K2,KG),T0(:,:,K2,KG),U0(:,:,K2,KG),V0(:,:,K2,KG),UF(:,:,K2,KG), &
        VF(:,:,K2,KG),HF(:,:,K2,KG),SF(:,:,K2,KG),SS(:,:,K2,KG),U(:,:,:,K2,KG),&
        V(:,:,:,K2,KG),W(:,:,:,K2,KG),A(:,:,:,K2,KG),T(:,:,:,K2,KG),           &
        Q(:,:,:,K2,KG),P(:,:,:,K2,KG),E(:,:,:,K2,KG),H(:,:,:,K2,KG),           &
        X(:,:,:,K2,KG),                                                        &
        iconvect,w0(:,:,K2,KG),tl(:,:,:,K2,KG),d(:,:,:,k2,kg), &
        sigw(:,:,:,K2,KG),zloc(:,:,K2,KG),dmass(:,:,:,K2,KG),   &
        ziscale,tlrams(:,:,:,K2,KG),mu(:,:,K2,KG),muu(:,:,K2,KG),              &
        muv(:,:,K2,KG),msfu(:,:,K2,KG),msfv(:,:,k2,kg),msft(:,:,k2,kg),        &
        fluxflg, deepflg, shallflg, cfxup1(:,:,:,k2,kg),cfxup2(:,:,:,k2,kg),   &
        cfxdn1(:,:,:,k2,kg),dfxup1(:,:,:,k2,kg),dfxup2(:,:,:,k2,kg),           &
        efxup1(:,:,:,k2,kg),efxup2(:,:,:,k2,kg),dfxdn1(:,:,:,k2,kg),           &
        efxdn1(:,:,:,k2,kg),tke(:,:,:,k2,kg),xm(:,:,:,k2,kg),hm(:,:,:,k2,kg))


!dwen(20090806)         CALL PRFCOM(ISOT,KGRID,KSFC,GX,GY,Z0,ZT,                       &
!dwen(20090806)     &      NXS,NYS,NZS,ZMDL,ZSG,NLVL,VMIX,ICONVECT,                    &
!dwen(20090806)     &      P0(:,:,K2),T0(:,:,K2),U0(:,:,K2),V0(:,:,K2),W0(:,:,K2),hpbl(:,:,K2),UF(:,:,K2),     &
!dwen(20090806)     &      VF(:,:,K2),SF(:,:,K2),U(:,:,:,K2),                          &
!dwen(20090806)     &      V(:,:,:,K2),                                                &
!dwen(20090806)     &      W(:,:,:,K2),T(:,:,:,K2),Q(:,:,:,K2),                        &
!dwen(20090806)     &      P(:,:,:,K2),D(:,:,:,K2),X(:,:,:,K2),TL(:,:,:,K2),           &
!dwen(20090806)     &      SIGW(:,:,:,K2),ZML(:,:,K2),ZLOC(:,:,K2),DMASS(:,:,:,K2),    &
!dwen(20090806)     &      ZISCALE,TLRAMS(:,:,:,K2),mu(:,:,K2),                        &
!dwen(20090806)     &      muu(:,:,K2),muv(:,:,K2),msfu(:,:,K2),msfv(:,:,K2),msft(:,:,K1), &
!dwen(20090806)     &      fluxflg, deepflg, shallflg, &
!dwen(20090806)     &      CFXUP1(:,:,:,K2),       &
!dwen(20090806)     &      CFXUP2(:,:,:,K2),CFXDN1(:,:,:,K2),DFXUP1(:,:,:,K2),         &
!dwen(20090806)     &      DFXUP2(:,:,:,K2),EFXUP1(:,:,:,K2),EFXUP2(:,:,:,K2),         &
!dwen(20090806)     &      DFXDN1(:,:,:,K2),EFXDN1(:,:,:,K2),TKEN(:,:,:,K2))

!*********************************************************
!    missing data metpos will wait longer before reading
     IF(FTIME(KG,K2).NE.MTIME(2))THEN
     WRITE(KF21,*)'WARNING advpnt: Time 2 input does not match request!'
        WRITE(KF21,*)' Internal time: ',JET
        WRITE(KF21,*)'  Time of data: ',FTIME(KG,K2)
        WRITE(KF21,*)'  Request time: ',MTIME(2)
        MTIME(2)=FTIME(KG,K2)
     END IF
  END IF

!-------------------------------------------------------------------------------
! determine minutes between data time periods for interpolation
! 21 May 2007 - preserve sign of DM to METWND subroutine
!-------------------------------------------------------------------------------

  DM=FLOAT(MTIME(2)-MTIME(1))
  IF(ABS(DM).GT.720.0)THEN
     WRITE(KF21,*)'*ERROR* advpnt: excessive interpolation'
     WRITE(KF21,*)' Greater than 720 min: ',ABS(DM)
     STOP 900
  END IF

!-----------------------------------------------------------------------
! optional vertical velocity remapping 
!----------------------------------------------------------------------

  IF(NEW1.AND.KVEL.GT.0)                                                       &
     CALL METWND(K1,K2,KVEL,NXS(KG),NYS(KG),NLVL,                              &
     DM,ZMDL,ZSG,ZT(:,:,KG),                                                   &
     U(:,:,:,K1,KG),V(:,:,:,K1,KG),W(:,:,:,K1,KG),                             &
     P(:,:,:,:,KG),T(:,:,:,:,KG),A(:,:,:,:,KG))

!************************************************
!dwen(20090323):try to use updated METWND
!     &   CALL METWND(KVEL,NXS,NYS,NZM,NLVL,DM,ZSG,                      &
!     &   U(:,:,:,K1),V(:,:,:,K1),W(:,:,:,K1),P(:,:,:,K1),P(:,:,:,K2),   &
!     &   T(:,:,:,K1),T(:,:,:,K2),D(:,:,:,K1),D(:,:,:,K2))
!************************************************

  IF(NEW2.AND.KVEL.GT.0)                                                       &
     CALL METWND(K1,K2,KVEL,NXS(KG),NYS(KG),NLVL,                              &
     DM,ZMDL,ZSG,ZT(:,:,KG),                                                   &
     U(:,:,:,K2,KG),V(:,:,:,K2,KG),W(:,:,:,K2,KG),                             &
     P(:,:,:,:,KG),T(:,:,:,:,KG),A(:,:,:,:,KG))

!************************************************
!dwen(20090323):try to use updated METWND
!     &   CALL METWND(KVEL,NXS,NYS,NZM,NLVL,DM,ZSG,                      &
!     &   U(:,:,:,K2),V(:,:,:,K2),W(:,:,:,K2),P(:,:,:,K1),P(:,:,:,K2),   &
!     &   T(:,:,:,K1),T(:,:,:,K2),D(:,:,:,K1),D(:,:,:,K2))
!
!************************************************
!-------------------------------------------------------------------------------
! compute new position
!-------------------------------------------------------------------------------

! map position from meteo grid to sub-grid
  XX=XP-LX1(KG)+1
  YY=YP-LY1(KG)+1
  ZZ=ZP

!*****************************************************
!dwen(20090323):adopted from STILT
! JCL4/13/00: following lines calculate the vertical profiles of TL & SIGMAW
!     at location BEFORE mean advection--will be used in PARDSP, along with
!     vertical profiles AFTER mean advection, for vertical turbulent transport
! JCL: a lot of the code below are similar to lines in ADVMET

! JCL:(4/25/00)have lines diff. from those in ADVMET b/c interpolation
!        factor for before advection should have value of 0 when subgrid is reloaded,
!        instead of value of 1 for the case after advection
!     interpolation factor for current time
      IF (ABS(DM) <= TINY(DM)) STOP 'ADVPNT: DM == 0. STOP.'
!dwen(20090825)      TF = MOD(DBLE(JET),DM)/DM
!dwen(20090825)     TF = MOD(float(JET),DM)/DM
!dwen(20090902)      IF (BACK .AND. (ABS(TF) >= TINY(TF))) TF=1d0-TF
!dwen(20090825)     IF (BACK .AND. (ABS(TF) >= TINY(TF))) TF=1.0-abs(TF)
      tf=FLOAT(ABS(JET-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))

! CHG:(12/04/01) set duration for conv. redistribution,
!     to be done once for every analysis time
         IF (ABS(TF) <= TINY(TF)) convdur1=DM
         IF(JET.EQ.CONC(1)%START%MACC)convdur1=IDNINT(DM*(1d0-TF))

! CHG(09/11/03) don't use this for RAMS
      IF(.NOT.RAMSFLG)THEN
      DO KL=1,NLVL
!        set vertical interpolation point to index position
!dwen(20090826)         ZK=DBLE(KL)
         ZK=float(KL)

! JCL:   interpolate Lagrangian timescale to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
!dwen(20090820) X1,Y1 are used in STILT version, not used in HYSPLIT version. X1 in STILT is 
!               equivalent to XX in HYSPLIT version, so replace X1 with XX, and Y1 with YY.
!dwen(20090820)         XX=DNINT(X1)
!dwen(20090820)         YY=DNINT(Y1)
!******************************
!dwen(2009086)   try to use adv3nt instead of advint
!******************************
!dwen         CALL ADVINT(TL(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADV3NT(TL(:,:,:,K1,kg),anint(XX),anint(YY),ZK,                 &
             var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)         CALL ADVINT(TL(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
         CALL ADV3NT(TL(:,:,:,K2,kg),anint(XX),anint(YY),ZK,                 &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
! JCL:(10/31/02) more sophisticated temporal interpolation that involves scaling by local zi
!         METO%TLPREV(KL)=(VAR2-VAR1)*TF+VAR1
         TLK1(KL)=VAR1
         TLK2(KL)=VAR2

! JCL:   interpolate std dev of vertical velocity to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
!dwen(20090820) X1,Y1 are used in STILT version, not used in HYSPLIT version. X1 in STILT is 
!               equivalent to XX in HYSPLIT version, so replace X1 with XX, and Y1 with YY.
!dwen(20090820)         XX=DNINT(X1)
!dwen(20090820)         YY=DNINT(Y1)
!dwen(20090806)         CALL ADVINT(SIGW(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,               &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADV3NT(SIGW(:,:,:,K1,kg),anint(XX),anint(YY),ZK,               &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)
!dwen(20090806)         CALL ADVINT(SIGW(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,               &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
         CALL ADV3NT(SIGW(:,:,:,K2,kg),anint(XX),anint(YY),ZK,               &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
! JCL:(10/31/02) more sophisticated temporal interpolation that involves scaling by local zi
!         METO%SIGWPREV(KL)=(VAR2-VAR1)*TF+VAR1
         SIGWK1(KL)=VAR1
         SIGWK2(KL)=VAR2
!         WRITE(45,*)"TLSIGW",KL,TLK1(KL),TLK2(KL),SIGWK1(KL),SIGWK2(KL)

! JCL:(5/9/01)interpolate horizontal density profiles to position
         if (awrfflg) then
                   !U(V) is staggered in x(y)-direction (C-grid) in WRF
! Note that x1 corresponds to position on mass grid, but staggered direction
! is off 0.5 (first staggered gridpoint is at position -0.5, so x1 on mass
! grid corresponds to x1+0.5 on staggered grid - this is different from RAMS)
!dwen(20090820) X1,Y1 are used in STILT version, not used in HYSPLIT version. X1 in STILT is 
!               equivalent to XX in HYSPLIT version, so replace X1 with XX, and Y1 with YY.
!dwen(20090820)            XX1=X1+0.5
!dwen(20090820)            YY1=Y1+0.5
            if (.not. fluxflg) then
! instantantenous velocities
!dwen(20090806)               CALL ADVINT(U(:,:,:,K1),NXS,NYS,NZM,XX1,Y1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
               CALL ADV3NT(U(:,:,:,K1,kg),XX+SUX,YY+SUY,ZK,                  &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)               CALL ADVINT(U(:,:,:,K2),NXS,NYS,NZM,XX1,Y1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               CALL ADV3NT(U(:,:,:,K2,kg),XX+SUX,YY+SUY,ZK,                  &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)

!dwen(20090826)               METz%UUPREV(KL)=(VAR2-VAR1)*TF+VAR1
               METz(kl)%UUPREV=(VAR2-VAR1)*TF+VAR1
!dwen(20090806)               CALL ADVINT(V(:,:,:,K1),NXS,NYS,NZM,X1,YY1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
               CALL ADV3NT(V(:,:,:,K1,kg),XX+SVX,YY+SVY,ZK,                  &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)               CALL ADVINT(V(:,:,:,K2),NXS,NYS,NZM,X1,YY1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               CALL ADV3NT(V(:,:,:,K2,kg),XX+SVX,YY+SVY,ZK,                  &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)

!dwen(20090826)               METz%VVPREV(KL)=(VAR2-VAR1)*TF+VAR1
               METz(kl)%VVPREV=(VAR2-VAR1)*TF+VAR1
            else
! time-averaged coupled u,v (decoupled in prfcom):
!dwen(20090813) *********************
!      k1 is always time last and k2 time next whether back or not, refer to metpos.f
!               IF(BACK)KLATE=K1
!              IF(.NOT.BACK)KLATE=K2
!dwen ********************************
!dwen(20090806)               CALL ADVINT(U(:,:,:,KLATE),NXS,NYS,NZM,XX1,Y1,ZK,             &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               CALL ADV3NT(U(:,:,:,K2,kg),XX+SUX,YY+SUY,ZK,             &
                var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)               METz%UUPREV(KL)=VAR2
               METz(kl)%UUPREV=VAR2
!dwen(20090806)               CALL ADVINT(V(:,:,:,KLATE),NXS,NYS,NZM,X1,YY1,ZK,             &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
               CALL ADV3NT(V(:,:,:,K2,kg),XX+SVX,YY+SVY,ZK,             &
                var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)               METz%VVPREV(KL)=VAR2
               METz(kl)%VVPREV=VAR2
            endif !fluxflg
         else !not awrfflg
!dwen(20090806)            CALL ADVINT(U(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
            CALL ADV3NT(U(:,:,:,K1,kg),XX,YY,ZK,                  &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)            CALL ADVINT(U(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
            CALL ADV3NT(U(:,:,:,K2,kg),XX,YY,ZK,                  &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
         
!dwen(20090826)            METz%UUPREV(KL)=(VAR2-VAR1)*TF+VAR1
            METz(kl)%UUPREV=(VAR2-VAR1)*TF+VAR1

!dwen(20090806)            CALL ADVINT(V(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
            CALL ADV3NT(V(:,:,:,K1,kg),XX,YY,ZK,                  &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)            CALL ADVINT(V(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
            CALL ADV3NT(V(:,:,:,K2,kg),XX,YY,ZK,                  &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)

!dwen(20090826)            METz%VVPREV(KL)=(VAR2-VAR1)*TF+VAR1
            METz(kl)%VVPREV=(VAR2-VAR1)*TF+VAR1
         endif

! JCL:   interpolate air density to position
!dwen(20090806)         CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
         CALL ADV3NT(D(:,:,:,K1,kg),XX,YY,ZK,                  &
                dens1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)         CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
         CALL ADV3NT(D(:,:,:,K2,kg),XX,YY,ZK,                  &
                dens2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)         METz%DENSPREV(KL)=(DENS2-DENS1)*TF+DENS1
         METz(kl)%DENSPREV=(DENS2-DENS1)*TF+DENS1

! CHG:(03/17/2004) interpolate temperature to position (need for dry air density)
!dwen(20090806)         CALL ADVINT(T(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP1)
         CALL ADV3NT(T(:,:,:,K1,kg),XX,YY,ZK,                  &
                temp1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)         CALL ADVINT(T(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP2)
         CALL ADV3NT(T(:,:,:,K2,kg),XX,YY,ZK,                  &
                temp2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)         METz%TEMPPREV(KL)=(TEMP2-TEMP1)*TF+TEMP1
         METz(kl)%TEMPPREV=(TEMP2-TEMP1)*TF+TEMP1

! CHG:(03/17/2004) interpolate rel. humidity to position (need for dry air density)
!dwen(20090806)         CALL ADVINT(Q(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR1)
         CALL ADV3NT(Q(:,:,:,K1,kg),XX,YY,ZK,                  &
                rhfr1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)         CALL ADVINT(Q(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,                  &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR2)
         CALL ADV3NT(Q(:,:,:,K2,kg),XX,YY,ZK,                  &
                rhfr2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)         METz%RHFRPREV(KL)=(RHFR2-RHFR1)*TF+RHFR1
         METz(kl)%RHFRPREV=(RHFR2-RHFR1)*TF+RHFR1

! JCL:(4/7/02)should NOT have spatial interpolation of mass violation b/c mass violation
!         was calculated for a gridcell, so should keep same value
! JCL:(4/3/02)interpolate mass violation (fraction of mass/min) profile to position
!         CALL ADVINT(DMASS(:,:,:,K1),NXS,NYS,NZM,X1,Y1,ZK,VAR1)
!         CALL ADVINT(DMASS(:,:,:,K2),NXS,NYS,NZM,X1,Y1,ZK,VAR2)
!         So decided to use ADVINT instead and feed it REAL values that are integers
!dwen(20090820) X1,Y1 are used in STILT version, not used in HYSPLIT version. X1 in STILT is 
!               equivalent to XX in HYSPLIT version, so replace X1 with XX, and Y1 with YY.
!dwen(20090820)         XX=DINT(X1)
!dwen(20090820)         YY=DINT(Y1)
         ZZa=aINT(ZK)
!dwen(20090806)         CALL ADVINT(DMASS(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZZ,              &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADV3NT(DMASS(:,:,:,K1,kg),aint(XX),aint(YY),ZZa,              &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)         CALL ADVINT(DMASS(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZZ,              &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
         CALL ADV3NT(DMASS(:,:,:,K2,kg),aint(XX),aint(YY),ZZa,              &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)

! JCL:   temporal interpolation
!dwen(20090826)         METz%DMASSPREV(KL)=(VAR2-VAR1)*TF+VAR1
         METz(kl)%DMASSPREV=(VAR2-VAR1)*TF+VAR1
! JCL:(4/3/02)the density change term is extremely small, so not worry
!        also need the density change term (kg/m3/min) as part of mass budget!
!dwen(20090826)         METz%DMASSPREV(KL)=METz%DMASSPREV(KL)+                         &
         METz(kl)%DMASSPREV=METz(kl)%DMASSPREV+                         &
     &                     ((DENS2-DENS1)/((DENS2+DENS1)/2.0))/DM
!        take into account the sign associated with direction in TIME
!dwen(20090826)         IF(BACK)METz%DMASSPREV(KL)=-1.0*METz%DMASSPREV(KL)
         IF(BACK)METz(kl)%DMASSPREV=-1.0*METz(kl)%DMASSPREV

      END DO

!CCCCCCCCCCCCCCCRAMS specialCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CHG(09/11/03) Do following steps for RAMS to "interpolate"
! Allways use K2 for fluxes (later time, they contain averages for period K1-K2
! Use nearest time for scalars
      ELSE IF(RAMSFLG)THEN
! CHG In Rams met files fluxes are staggered 1/2 Dx aggainst scalars
! location of (i,j)=(1,1) is scalar location (center of grid cell)
! all dummies should have been skipped in the rams2arl job file
! so assume i=2.3, means need scalars at 2, and fluxes at 1.8
! so assume i=2.7, means need scalars at 3, and fluxes at 2.2
                   !staggered
!dwen(20090820) X1,Y1 are used in STILT version, not used in HYSPLIT version. X1 in STILT is 
!               equivalent to XX in HYSPLIT version, so replace X1 with XX, and Y1 with YY.
!dwen(2009082)        XX1=X1-0.5
                   !staggered
!dwen(2009082)        YY1=Y1-0.5
                     !DNINT for scalars; for x-fluxes use XX1, not XX
!dwen(2009082)        XX=DNINT(X1)
                     !DNINT for scalars; for y-fluxes use YX1, not YY
!dwen(2009082)        YY=DNINT(Y1)
!        WRITE(*,*)'advpnt: XX1,YY1',XX1,YY1
!dwen(20090822)        IF(X1.GE.DBLE(NXS))WRITE(*,*)'advpnt: X1 off grid, >NXS',XX1
!dwen(20090822)        IF(Y1.GE.DBLE(NYS))WRITE(*,*)'advpnt: Y1 off grid, >NYS',YY1
!dwen(20090822)        IF(XX1.LE.1.0)WRITE(*,*)'advpnt: XX1 off grid, <1.0',XX1
!dwen(20090822)        IF(YY1.LE.1.0)WRITE(*,*)'advpnt: YY1 off grid, <1.0',YY1
        IF(XX.GE.DBLE(NXS(kg)))WRITE(*,*)'advpnt: X1 off grid, >NXS',XX-0.5
        IF(YY.GE.DBLE(NYS(kg)))WRITE(*,*)'advpnt: Y1 off grid, >NYS',YY-0.5
        IF(XX-0.5.LE.1.0)WRITE(*,*)'advpnt: XX1 off grid, <1.0',XX-0.5
        IF(YY-0.5.LE.1.0)WRITE(*,*)'advpnt: YY1 off grid, <1.0',YY-0.5

!dwen(20090813)        IF(TF.GE.0.5)THEN
!dwen(20090813)         KNEAR=K2
!dwen(20090813)          KFAR=K1
!dwen(20090813)       ELSE
!dwen(20090813)         KNEAR=K1
!dwen(20090813)          KFAR=K2
!dwen(20090813)       END IF
!dwen(20090813) *********************
!      k1 is always time last and k2 time next whether back or not, refer to metpos.f
!        IF(BACK)KLATE=K1
!        IF(.NOT.BACK)KLATE=K2
!dwen *******************************
        DO KL=1,NLVL
!        set vertical interpolation point to index position
!dwen(20090828)          ZK=DBLE(KL)
          ZK=float(KL)

!dwen(20090806)          CALL ADVINT(TLRAMS(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,         &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
        IF(TF.GE.0.5)THEN
          CALL ADV3NT(TLRAMS(:,:,:,K2,kg),anint(XX),anint(YY),ZK,         &
                 var1,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)

!dwen(20090806)          CALL ADVINT(SIGWRAMS(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,       &
!dwen(20090806)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
          CALL ADV3NT(SIGWRAMS(:,:,:,K2,kg),anint(XX),anint(YY),ZK,       &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
        else
           CALL ADV3NT(TLRAMS(:,:,:,K1,kg),anint(XX),anint(YY),ZK,         &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

          CALL ADV3NT(SIGWRAMS(:,:,:,K1,kg),anint(XX),anint(YY),ZK,       &
                 var2,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)
        end if
 

! DMM: 11/9/2004 Negative values should not be present in TLPREV; it causes pardsp to hang
        IF(VAR1.LT.0.0) THEN
!dwen(20090826)            METz%TLPREV(KL)=0.0
            METz(kl)%TLPREV=0.0
        ELSE
!dwen(20090826)            METz%TLPREV(KL)=VAR1
            METz(kl)%TLPREV=VAR1
        ENDIF
!dwen(20090826)          METz%SIGWPREV(KL)=VAR2
          METz(kl)%SIGWPREV=VAR2

!dwen(20090806)          CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
          CALL ADV3NT(D(:,:,:,K1,kg),anint(XX),anint(YY),ZK,                 &
                dens1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)          CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
          CALL ADV3NT(D(:,:,:,K2,kg),anint(XX),anint(YY),ZK,                 &
                dens2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)          METz%DENSPREV(KL)=(DENS2-DENS1)*TF+DENS1
          METz(kl)%DENSPREV=(DENS2-DENS1)*TF+DENS1

! CHG:(03/17/2004) interpolate temperature to position (need for dry air density)
!dwen(20090806)          CALL ADVINT(T(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP1)
          CALL ADV3NT(T(:,:,:,K1,kg),anint(XX),anint(YY),ZK,                 &
                temp1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)          CALL ADVINT(T(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,TEMP2)
          CALL ADV3NT(T(:,:,:,K2,kg),anint(XX),anint(YY),ZK,                 &
                temp2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)          METz%TEMPPREV(KL)=(TEMP2-TEMP1)*TF+TEMP1
          METz(kl)%TEMPPREV=(TEMP2-TEMP1)*TF+TEMP1

! CHG:(03/17/2004) interpolate rel. humidity to position (need for dry air density)
!dwen(20090806)          CALL ADVINT(Q(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR1)
          CALL ADV3NT(Q(:,:,:,K1,kg),anint(XX),anint(YY),ZK,                 &
                rhfr1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090806)          CALL ADVINT(Q(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,RHFR2)
          CALL ADV3NT(Q(:,:,:,K2,kg),anint(XX),anint(YY),ZK,                 &
                rhfr2,GRID(kt,kt2)%GLOBAL,GRID(kt,kt2)%NX,GRID(kt,kt2)%NY)
!dwen(20090826)          METz%RHFRPREV(KL)=(RHFR2-RHFR1)*TF+RHFR1
          METz(kl)%RHFRPREV=(RHFR2-RHFR1)*TF+RHFR1

!dwen(20090806)          CALL ADVINT(U(:,:,:,KLATE),NXS,NYS,NZM,XX1,YY,ZK,             &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
          CALL ADV3NT(U(:,:,:,K2,kg),XX-0.5,anint(YY),ZK,             &
                var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)          METz%UUPREV(KL)=VAR2/DENS2
          METz(kl)%UUPREV=VAR2/DENS2
!dwen(20090806)          CALL ADVINT(V(:,:,:,KLATE),NXS,NYS,NZM,XX,YY1,ZK,             &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
          CALL ADV3NT(V(:,:,:,K2,kg),anint(XX),YY-0.5,ZK,             &
                var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)          METz%VVPREV(KL)=VAR2/DENS2
          METz(kl)%VVPREV=VAR2/DENS2

!dwen(20090806)          CALL ADVINT(DMASS(:,:,:,KLATE),NXS,NYS,NZM,XX,YY,ZK,          &
!dwen(20090806)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADV3NT(DMASS(:,:,:,K2,kg),anint(XX),anint(YY),ZK,          &
               var1,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)          METz%DMASSPREV(KL)=VAR1
          METz(kl)%DMASSPREV=VAR1
! JCL:(4/3/02)the density change term is extremely small, so not worry
!        also need the density change term (kg/m3/min) as part of mass budget!
!        WRITE(45,*)'advpnt,divu:',METO%DMASSPREV(KL)/(DENS2+DENS1),KL
!dwen(20090826)          METz%DMASSPREV(KL)=METz%DMASSPREV(KL)*2.0/(DENS2+DENS1)       &
          METz(kl)%DMASSPREV=METz(kl)%DMASSPREV*2.0/(DENS2+DENS1)       &
                          +((DENS2-DENS1)/((DENS2+DENS1)/2.0))/DM
!dwen(20090826)          IF(BACK)METz%DMASSPREV(KL)=-1.0*METz%DMASSPREV(KL)
          IF(BACK)METz(kl)%DMASSPREV=-1.0*METz(kl)%DMASSPREV
!        WRITE(45,*)'K1, K2:',K1,K2
        END DO
      END IF

! CHG(09/11/03) END Rams special
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! JCL:(4/28/00)interpolate horizontally to get mixed-layer height
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
! JCL:(10/30/02)no longer conduct horizontal spatial interpolation, since want ZML to be consistent w/ turbulence profiles
!dwen(20090820) X1,Y1 are used in STILT version, not used in HYSPLIT version. X1 in STILT is 
!               equivalent to XX in HYSPLIT version, so replace X1 with XX, and Y1 with YY.
!dwen(20090820)      XX=DNINT(X1)
!dwen(20090820)      YY=DNINT(Y1)
      ZZZ=1.0
!dwen(20090910)      CALL ADVINTZML(ZML(:,:,K1:K1),NXS,NYS,XX,YY,ZZZ,                   &
!dwen(20090910)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)

      CALL ADV3NTZML(ZI(:,:,K1:K1,kg),anint(XX),anint(YY),ZZZ,                   &
              var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090810)      CALL ADVINTZML(ZML(:,:,K2:K2),NXS,NYS,XX,YY,ZZZ,                   &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      CALL ADV3NTZML(ZI(:,:,K2:K2,kg),anint(XX),anint(YY),ZZZ,                   &
              var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
      ZMLK1=VAR1
      ZMLK2=VAR2
      METO%ZMLPREV=(VAR2-VAR1)*TF+VAR1

! CHG:(12/04/01)get NEAREST gridpoints value for ZLOC (limit of convection
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
      ZZZ=1.0
!dwen(20090810)      CALL ADVINTZLOC(ZLOC(:,:,K1:K1),NXS,NYS,X1,Y1,ZZZ,                 &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADV3NTZLOC(ZLOC(:,:,K1:K1,kg),XX,YY,ZZZ,                 &
              var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090810)      CALL ADVINTZLOC(ZLOC(:,:,K2:K2),NXS,NYS,X1,Y1,ZZZ,                 &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      CALL ADV3NTZLOC(ZLOC(:,:,K2:K2,kg),XX,YY,ZZZ,                 &
              var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
      IF(TF.LE.0.5)METO%ZLOCPREV=VAR1
      IF(TF.GT.0.5)METO%ZLOCPREV=VAR2


!CCCCCCCCCCC EXTRACT TEMPORALLY-INTERPOLATED SIGMAW & TL ***BEFORE*** ADVECTION
! JCL:(10/31/02) New way to interpolate vertical profiles of sigmaw & TL between analysis times by scaling
!         the two analysis profiles with the temporally interpolated mixed-layer height
! 1) Find nearest model level to move Zi to
!dwen(20090820) X1,Y1 are used in STILT version, not used in HYSPLIT version. X1 in STILT is 
!               equivalent to XX in HYSPLIT version, so replace X1 with XX, and Y1 with YY.
!dwen(20090820)      XX=DNINT(X1)
!dwen(20090820)      YY=DNINT(Y1)
! JCL:(10/31/02) extract ground height with new interpolation routine that deals with 2D arrays
!dwen(20090810)      CALL ADVINT2D(ZT,NXS,NYS,XX,YY,                                   &
!dwen(20090810)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADV2NT(ZT(:,:,kg),anint(XX),anint(YY),                                   &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)
      ZGRD=VAR1
! calculate true agl heights [m] of model levels:
      do kk=1,nlvl
         ZLVLS(kk)=(1.0-ZSG(kk))*(ZMDL-ZGRD)
      end do
! Force zmix to be at a model level if KMIX0 < 0 (as in stbanl,
! redone here because of horizontal interpolation of zi, zgrd)
      zmix=METO%ZMLPREV
      if (KMIX0 .lt. 0) then
         ! Find closest model level above k=1
         zmixnew=zlvls(2)
         do kk=3,nlvl
            IF(ABS(ZLVLS(kK)-ZMIX) .LT. ABS(ZMIXNEW-ZMIX))THEN
               ZMIXNEW=ZLVLS(kK)
            END IF
         enddo
         METO%ZMLPREV=zmixnew
      endif
      

! 2)  Find CLOSEST vertical level based on coordinate that's SCALED BY zi
! CHG(09/11/03) don't do this for RAMS
      IF(.NOT.RAMSFLG)THEN
         ZMIX=METO%ZMLPREV
         KK1=1
         KK2=1
         DO KK=1,NLVL
           ! altitude of current timestep (scaled by local zi)
           ZNOWSC=ZLVLS(KK)/ZMIX
           IF (KK1 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK1)/ZMLK1-ZNOWSC).GT.                      &
                            ABS(ZLVLS(KK1+1)/ZMLK1-ZNOWSC))
                 KK1=KK1+1
                 IF (KK1 == NLVL) EXIT
              END DO
           END IF
           IF (KK2 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK2)/ZMLK2-ZNOWSC).GT.                      &
                            ABS(ZLVLS(KK2+1)/ZMLK2-ZNOWSC))
                 KK2=KK2+1
                 IF (KK2 == NLVL) EXIT
              END DO
           END IF
!dwen(20090826)           METz%SIGWPREV(KK)=(SIGWK2(KK2)-SIGWK1(KK1))*TF+SIGWK1(KK1)
!dwen(20090826)           METz%TLPREV(KK)=(TLK2(KK2)-TLK1(KK1))*TF+TLK1(KK1)
           METz(kk)%SIGWPREV=(SIGWK2(KK2)-SIGWK1(KK1))*TF+SIGWK1(KK1)
           METz(kk)%TLPREV=(TLK2(KK2)-TLK1(KK1))*TF+TLK1(KK1)
         END DO
      ELSE
         ! CHG(09/11/03) assign vertical index here for RAMS, not in ADVIEC
         ! Do this consistent with other runs: ZNDX can be < 1 for below 1st level
         ! For RAMS the 1st level is a flux level...
         ! Here use DREC()%HEIGHT rather than internal heights ZSG(), since need to compare
         ! heights in magl; note there is a shift in indices: HEIGHT(2)~ZSG(1); see hymodelc
         KGNOW = MAX(KG,KGOLD)
         KZ=1
!dwen(20090825)         Z1Z=ZMDL*(1.0-DMIN1(DBLE(1.0),Z1))
         Z1Z=ZMDL*(1.0-aMIN1(1.0,ZZ))
!dwen(20090825)         DO WHILE(DREC(KGNOW)%HEIGHT(KZ).LT.Z1Z)
         DO WHILE(DREC(KGNOW,kt1)%HEIGHT(KZ).LT.Z1Z)
           KZ=KZ+1
         END DO
         KZ=KZ-1
!dwen(20090825)         IF (KZ.GT.0) FRAC=(Z1Z-DREC(KGNOW)%HEIGHT(KZ))/                   &
!dwen(20090825)                              (DREC(KGNOW)%HEIGHT(KZ+1)-DREC(KGNOW)%HEIGHT(KZ))
!dwen(20090825)         IF(KZ.EQ.0)FRAC=Z1Z/DREC(KGNOW)%HEIGHT(KZ+1)
         IF (KZ.GT.0) FRAC=(Z1Z-DREC(KGNOW,kt1)%HEIGHT(KZ))/                   &
                              (DREC(KGNOW,kt1)%HEIGHT(KZ+1)-DREC(KGNOW,kt1)%HEIGHT(KZ))
         IF(KZ.EQ.0)FRAC=Z1Z/DREC(KGNOW,kt1)%HEIGHT(KZ+1)
         METO%ZNDX=DBLE(KZ)+FRAC-1.0
      END IF
!*****************************************************

! advection for one time step
  HDWPX=HDWP
  IF(HDWPX.GE.100) HDWPX=MOD(HDWP/10,10) ! complex mode

  IF(HDWPX.LE.4)THEN
!    standard 3D atmpospheric advection
!dwen(20090810)     CALL ADVIEC(U(:,:,:,:,KG),V(:,:,:,:,KG),W(:,:,:,:,KG),                    &
!dwen(20090810)          K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,METO%ZNDX,DT,TRATIO,BACK,         &
!dwen(20090810)          GRID(KG,KT1)%GLOBAL,GRID(KG,KT1)%NX,GRID(KG,KT1)%NY)
!********************************************************
!dwen(20090323):try to use updated ADVIEC
! JCL:(3/1/01)add WWOUT to output mean vertical velocity [sigma/min]
! JCL:add GD, DUDZ, DVDZ to argument list--used by CHG for box size
!     compute new position
! JCL:(8/13/01)add grdht(METO%ZTER) & 1st sigma level to calculate AGL accurately--to enable accurate interpolation
! JCL:(09/01/03)pass on wind error flag, UPRIME/UERR and VPRIME/VERR from previous timestep
! JCL:(09/01/03) pass on random seed 'RSEED' for modeling transport error as stochastic process
! JCL:(09/01/03) DXERR and DYERR are the horizontal displacements resulting from transport error
! CHG(09/11/03) Pass on RAMSFLG
! CHG(09/11/03) Pass on Density fields for RAMS
! JCL:(11/03/03) remove winderror arguments--do all calculations in HYMODELC
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
! dwen(20090810)   DVDZ and DUDZ are useless in ADVIEC, removed them
     CALL ADVIEC(U(:,:,:,:,KG),V(:,:,:,:,KG),W(:,:,:,:,KG),d(:,:,:,:,kg),      &
          K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,METO%ZNDX,DT,TRATIO,BACK,         &
          GRID(KG,KT1)%GLOBAL,GRID(KG,KT1)%NX,GRID(KG,KT1)%NY,                 &
          meto%zter,zsg(1),wwout,awrfflg,fluxflg,zsg,ramsflg,dead)
     kgold = 0
     if (dead) then
        WRITE(KF21,*)' NOTICE advpnt: dead particle with xx,yy,zz,kg=',xx,yy,zz,kg
        kgold=kg
        kg=0
        kret=1
        RETURN
     end if

!dwen(20090810)      CALL ADVIEC                                                       &
!dwen(20090810)     &   (U(:,:,:,K1), V(:,:,:,K1), W(:,:,:,K1),                        &
!dwen(20090810)     &    U(:,:,:,K2), V(:,:,:,K2), W(:,:,:,K2),                        &
!dwen(20090810)     &    D(:,:,:,K1), D(:,:,:,K2),                                     &
!dwen(20090810)     &    NXS,NYS,NZM,NLVL,DM,JET,ZMDL,METO%ZTER,ZSG(1),                &
!dwen(20090810)     &    X1,Y1,Z1,X2,Y2,Z2,METO%ZNDX,DT,BACK,GX,GY,DUDZ,DVDZ,WWOUT,    &
!dwen(20090810)     &    awrfflg, fluxflg, zsg, &
!dwen(20090810)     &    RAMSFLG,ECMFLG,GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DEAD)
!
!********************************************************
  ELSEIF(HDWPX.EQ.5)THEN
!    test for special surface advecting particles
     CALL ADVSFC(UF(:,:,:,KG),VF(:,:,:,KG),                                    &
          K1,K2,NLVL,MTIME,JET,XX,YY,DT,BACK,                                  &
          GRID(KG,KT1)%GLOBAL,GRID(KG,KT1)%NX,GRID(KG,KT1)%NY)

  ELSEIF(HDWPX.EQ.6)THEN
!    lagrangian isobaric sampling
     CALL ADVISO(U(:,:,:,:,KG),V(:,:,:,:,KG),P(:,:,:,:,KG),ZSG,                &
          K1,K2,NLVL,MTIME,JET,ZMDL,XX,YY,ZZ,METO%ZNDX,DT,TRATIO,BACK,         &
          GRID(KG,KT1)%GLOBAL,GRID(KG,KT1)%NX,GRID(KG,KT1)%NY,                 &
          METO%UVEL,METO%VVEL,RAMSFLG)

  ELSE
     KG=0
     KRET=1
     WRITE(KF21,*)' NOTICE advpnt: invalid HDWP/INITD option - ',HDWP
     RETURN
  END IF

! save advection distance as a wind speed (grid pts / min)
  UBAR = MAX(ABS(XX-(XP-LX1(KG)+1)),ABS(YY-(YP-LY1(KG)+1)))/ABS(DT)

! check for cyclic boundary conditions
  IF(GRID(KG,KT1)%GLOBAL)THEN
     IF(XX.GE.FLOAT(GRID(KG,KT1)%NX+1)) XX=XX-FLOAT(GRID(KG,KT1)%NX)
     IF(XX.LT.1.0)                      XX=GRID(KG,KT1)%NX+XX
     IF(YY.GT.FLOAT(GRID(KG,KT1)%NY))   YY=2.0*GRID(KG,KT1)%NY-YY
     IF(YY.LT.1.0)                      YY=2.0-YY
  ELSE
!    This condition should never occur after the advection step because
!    limits have been set to switch particles when too close to edge of
!    subgrid prior to calling the advection step. However, some multi-
!    meteo grid runs may have occasional problems with particle switch.
!    In this situation the particle is just terminated.
     xx_min = min(min(xx,xx+sux),xx+svx)
     xx_max = max(max(xx,xx+sux),xx+svx)
     yy_min = min(min(yy,yy+suy),yy+svy)
     yy_max = max(max(yy,yy+suy),yy+svy)
     IF(xx_min .LT. 1.0 .OR. yy_min .LT. 1.0 .OR. xx_max .GT. FLOAT(NXS(KG)) .OR.                     &
        yy_max .GT. FLOAT(NYS(KG)) )THEN
        WRITE(KF21,*)' NOTICE advpnt: particle with xx_min, xx_max, yy_min, yy_max ', &
             & xx_min, xx_max, yy_min, yy_max,' off subgrid for kg,nxs,nys ', &
             & kg,nxs(kg),nys(kg)
        KG=0
        KRET=1
        RETURN
     END IF
  END IF

! map position back to meteo grid
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2.
  XP=XX+LX1(KG)-1
  YP=YY+LY1(KG)-1
  ZP=MIN(1.0,ZZ)

! convert end-point grid position to true coordinates
  IF(GRID(KG,KT1)%LATLON)THEN
     CALL GBL2LL(KG,KT1,XP,YP,METO%PLAT,METO%PLON)
  ELSE
     CALL CXY2LL_wps(GRID(KG,KT1)%GBASE,XP,YP,METO%PLAT,METO%PLON,GRID(KG,KT1)%proj)
  END IF

! meteo variables interpolated to last advection point
  JTIME=JET+NINT(DT)

!dwen(20090811)  CALL ADVMET(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,XX,YY,JTIME,MTIME,            &
!dwen(20090811)              DREC(KG,KT1)%ACYCLE,NLVL,FHOUR(:,KG),IFHR,K1,K2,GX(:,:,KG),      &
!dwen(20090811)              GY(:,:,KG),Z0(:,:,KG),LU(:,:,KG),ZT(:,:,KG),                     &
!dwen(20090811)              A(:,:,:,:,KG),T(:,:,:,:,KG),Q(:,:,:,:,KG),P(:,:,:,:,KG),         &
!dwen(20090811)              E(:,:,:,:,KG),X(:,:,:,:,KG),ZI(:,:,:,KG), H(:,:,:,:,KG),         &
!dwen(20090811)              U0(:,:,:,KG),V0(:,:,:,KG),RT(:,:,:,KG),UF(:,:,:,KG),VF(:,:,:,KG),&
!dwen(20090811)              SF(:,:,:,KG), SS(:,:,:,KG),DS(:,:,:,KG),DREC(KG,KT1)%DSWF,       &
!dwen(20090811)              GRID(KG,KT1)%GLOBAL,GRID(KG,KT1)%NX,GRID(KG,KT1)%NY)

!************************************************
!dwen(20090323):adopted from STILT


! CHG(09/16/03)pass on RAMSFLG
!dwen(20090731)  note that SF is sensible heat flux, not stability function in ADVMET in STILT
!               x: w-component turbulence
!dwen(20090811) xm: vertical mixing coefficient, and hm: horizontal mixing coefficient
!dwen(20090812) remove radf, use rdep instead
!dwen(20090812) remove sw, use ds instead
!               hf: sensible heat flux ;   sf: stability function

  CALL ADVMET(METZ,METO,BACK,VMIX,CDEP,RDEP,TRAJ,XX,YY,JTIME,MTIME,            &
              DREC(KG,KT1)%ACYCLE,NLVL,FHOUR(:,KG),IFHR,K1,K2,GX(:,:,KG),      &
              GY(:,:,KG),Z0(:,:,KG),LU(:,:,KG),ZT(:,:,KG),                     &
              A(:,:,:,:,KG),T(:,:,:,:,KG),Q(:,:,:,:,KG),P(:,:,:,:,KG),         &
              E(:,:,:,:,KG),X(:,:,:,:,KG),ZI(:,:,:,KG), H(:,:,:,:,KG),         &
              U0(:,:,:,KG),V0(:,:,:,KG),RT(:,:,:,KG),UF(:,:,:,KG),VF(:,:,:,KG),&
              SF(:,:,:,KG), SS(:,:,:,KG),DS(:,:,:,KG),DREC(KG,KT1)%DSWF,       &
              GRID(KG,KT1)%GLOBAL,GRID(KG,KT1)%NX,GRID(KG,KT1)%NY,             &
              drec(kg,kt1)%cflg,drec(kg,kt1)%tclf,drec(kg,kt1)%lclf,           &
              drec(kg,kt1)%radf,drec(kg,kt1)%slmf,                             &
              u(:,:,:,:,kg),v(:,:,:,:,kg),d(:,:,:,:,kg),                       &
              xm(:,:,:,:,kg),hm(:,:,:,:,kg),                                   &
              dmass(:,:,:,:,kg),tl(:,:,:,:,kg),sigw(:,:,:,:,kg),rc(:,:,:,kg),  &
              hf(:,:,:,kg),lf(:,:,:,kg),tc(:,:,:,kg),lc(:,:,:,kg),             &
              sm(:,:,:,kg),ramsflg,ecmflg,awrfflg,fluxflg)
!
!dwen(20090811)     CALL ADVMET(BACK,VMIX,CDEP,RDEP,TRAJ,X2,Y2,TLON,TLAT,JTIME,DM,      &
!dwen(20090811)     &   KVEL,DREC(KGRID)%ACYCLE,NXS,NYS,NLVL,FHOUR,IFHR,K1,K2,           &
!dwen(20090811)     &   GX,GY,Z0,LU,ZT,DREC(KGRID)%CFLG,DREC(KGRID)%TCLF,                &
!dwen(20090811)     &   DREC(KGRID)%LCLF,DREC(KGRID)%RADF,DREC(KGRID)%SLMF,              &
!dwen(20090811)     &   T(:,:,:,K1),Q(:,:,:,K1),P(:,:,:,K1),D(:,:,:,K1),                 &
!dwen(20090811)     &   X(:,:,:,K1),H(:,:,:,K1),RT(:,:,K1),UF(:,:,K1),VF(:,:,K1),        &
!dwen(20090811)     &   U(:,:,:,K1),V(:,:,:,K1),DMASS(:,:,:,K1),                         &
!dwen(20090811)     &   T(:,:,:,K2),Q(:,:,:,K2),P(:,:,:,K2),D(:,:,:,K2),                 &
!dwen(20090811)     &   X(:,:,:,K2),H(:,:,:,K2),RT(:,:,K2),UF(:,:,K2),VF(:,:,K2),        &
!dwen(20090811)     &   U(:,:,:,K2),V(:,:,:,K2),DMASS(:,:,:,K2),                         &
!dwen(20090811)     &   TL(:,:,:,K1),TL(:,:,:,K2),SIGW(:,:,:,K1),SIGW(:,:,:,K2),         &
!dwen(20090811)     &   RC(:,:,K1),RC(:,:,K2),                                           &
!dwen(20090811)     &   SF(:,:,K1),SF(:,:,K2),LF(:,:,K1),LF(:,:,K2),                     &
!dwen(20090811)     &   TC(:,:,K1),TC(:,:,K2),LC(:,:,K1),LC(:,:,K2),                     &
!dwen(20090811)     &   SW(:,:,K1),SW(:,:,K2),SM(:,:,K1),SM(:,:,K2),                     &
!dwen(20090811)     &   RAMSFLG,ECMFLG,GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY, &
!dwen(20090811)     &   awrfflg, fluxflg)

! CHG(09/23/03) call ADVMETGRELL to prepare profiles of convective fluxes
      IF (RAMSFLG .OR. GRID(kg,kt1)%model_id(1:3) == 'ECX' .or. (awrfflg .and. (deepflg .or. shallflg))) then
         IF (BACK) THEN
! CHG(09/23/03) for back, use K1 (absolute later time)
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
!dwen(20090812) *******************
!dwen(20090812)            CALL ADVMETGRELL (X2,Y2,CFXUP1(:,:,:,K1),CFXUP2(:,:,:,K1),            &
!dwen(20090812)                              CFXDN1(:,:,:,K1),DFXUP1(:,:,:,K1),DFXUP2(:,:,:,K1), &
!dwen(20090812)                              EFXUP1(:,:,:,K1),EFXUP2(:,:,:,K1),DFXDN1(:,:,:,K1), &
!dwen(20090812)                              EFXDN1(:,:,:,K1),TKEN(:,:,:,K1),NXS,NYS,NLVL,       &
!dwen(20090812)                              GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY)
!dwen(20090812)           add  METO
!dwen(20090812)           remove NXS,NYS          
!dwen(20090812)           use TKE instead of TKEN to represent turbulent kinetic energy  
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2. so replace X2 with 
!               XX, Y2 with YY and Z2 with ZZ
            CALL ADVMETGRELL (METz,XX,YY,CFXUP1(:,:,:,K1,kg),CFXUP2(:,:,:,K1,kg),          &
                              CFXDN1(:,:,:,K1,kg),DFXUP1(:,:,:,K1,kg),DFXUP2(:,:,:,K1,kg), &
                              EFXUP1(:,:,:,K1,kg),EFXUP2(:,:,:,K1,kg),DFXDN1(:,:,:,K1,kg), &
                              EFXDN1(:,:,:,K1,kg),TKE(:,:,:,K1,kg),NLVL,                  &
                              GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)
 
         ELSE
! CHG(09/23/03) for forward, use K2 (absolute later time)
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
!dwen(20090812)            CALL ADVMETGRELL (X2,Y2,CFXUP1(:,:,:,K2),CFXUP2(:,:,:,K2),            &
!dwen(20090812)                              CFXDN1(:,:,:,K2),DFXUP1(:,:,:,K2),DFXUP2(:,:,:,K2), &
!dwen(20090812)                              EFXUP1(:,:,:,K2),EFXUP2(:,:,:,K2),DFXDN1(:,:,:,K2), &
!dwen(20090812)                              EFXDN1(:,:,:,K2),TKEN(:,:,:,K2),NXS,NYS,NLVL,       &
!dwen(20090812)                              GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY)
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2. so replace X2 with 
!               XX, Y2 with YY and Z2 with ZZ
            CALL ADVMETGRELL (METz,XX,YY,CFXUP1(:,:,:,K2,kg),CFXUP2(:,:,:,K2,kg),          &
                              CFXDN1(:,:,:,K2,kg),DFXUP1(:,:,:,K2,kg),DFXUP2(:,:,:,K2,kg), &
                              EFXUP1(:,:,:,K2,kg),EFXUP2(:,:,:,K2,kg),DFXDN1(:,:,:,K2,kg), &
                              EFXDN1(:,:,:,K2,kg),TKE(:,:,:,K2,kg),NLVL,                  &
                              GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
 
         END IF
      END IF

! JCL:(4/28/00)interpolate horizontally to get mixed-layer height
!     interpolation factor for time after advection
!dwen(20090825)      TF=DMOD(DBLE(JTIME),DM)/DM
!dwen(20090902)      TF=aMOD(float(JTIME),DM)/DM
!dwen(20090902)      IF(BACK)THEN
!dwen(20090902)         TF=1.0-TF
!dwen(20090902)           !make sure after advection use K2 when TF=0.0 (JTIME updated)
!dwen(20090902)      ELSE
!dwen(20090902)         IF(TF.EQ.0.0)TF=1.0
!dwen(20090902)      END IF

      tf=FLOAT(ABS(Jtime-MTIME(1)))/FLOAT(ABS(MTIME(2)-MTIME(1)))

! JCL:(10/31/02)extract profiles of TL & SIGMAW at the advected particle location for the two analysis times
! CHG(09/11/03) don't use this for RAMS
      IF(.NOT.RAMSFLG)THEN
      DO KL=1,NLVL
!        set vertical interpolation point to index position
         ZK=float(KL)

! JCL:   interpolate Lagrangian timescale to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2. so replace X2 with 
!               XX, Y2 with YY and Z2 with ZZ
!dwen(20090820)         XX=DNINT(X2)
!dwen(20090820)         YY=DNINT(Y2)
!dwen(20090807)         CALL ADVINT(TL(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090807)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADV3NT(TL(:,:,:,K1,kg),anint(XX),anint(YY),ZK,                 &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090807)         CALL ADVINT(TL(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,                 &
!dwen(20090807)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
         CALL ADV3NT(TL(:,:,:,K2,kg),anint(XX),anint(YY),ZK,                 &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
         TLK1(KL)=VAR1
         TLK2(KL)=VAR2

! JCL:   interpolate std dev of vertical velocity to position
! JCL:(10/30/02)no longer conduct linear interpolation, since screws up turbulence profiles
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2. so replace X2 with 
!               XX, Y2 with YY and Z2 with ZZ
!dwen(20090820)         XX=DNINT(X2)
!dwen(20090820)         YY=DNINT(Y2)
!dwen(20090807)         CALL ADVINT(SIGW(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZK,               &
!dwen(20090807)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
         CALL ADV3NT(SIGW(:,:,:,K1,kg),anint(XX),anint(YY),ZK,               &
                 var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090807)         CALL ADVINT(SIGW(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZK,               &
!dwen(20090807)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
         CALL ADV3NT(SIGW(:,:,:,K2,kg),anint(XX),anint(YY),ZK,               &
                 var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
         SIGWK1(KL)=VAR1
         SIGWK2(KL)=VAR2

      END DO
             !of IF(.NOT.RAMSFLG)THEN
      END IF

! JCL:spatial interpolation to get local ML ht with 'ADVINTZML'
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
! JCL:(10/30/02)no longer conduct horizontal spatial interpolation, since want ZML to be consistent w/ turbulence profiles
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2. so replace X2 with 
!               XX, Y2 with YY and Z2 with ZZ
!dwen(20090820)      XX=DNINT(X2)
!dwen(20090820)      YY=DNINT(Y2)
      ZZZ=1.0
!dwen(20090810)      CALL ADVINTZML(ZML(:,:,K1:K1),NXS,NYS,XX,YY,ZZZ,                   &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADV3NTZML(ZI(:,:,K1:K1,kg),anint(XX),anint(YY),ZZZ,                   &
              var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090810)      CALL ADVINTZML(ZML(:,:,K2:K2),NXS,NYS,XX,YY,ZZZ,                   &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      CALL ADV3NTZML(ZI(:,:,K2:K2,kg),anint(XX),anint(YY),ZZZ,                   &
              var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
      ZMLK1=VAR1
      ZMLK2=VAR2
      METO%ZMLNEXT=(VAR2-VAR1)*TF+VAR1

! CHG:(12/04/01)get nearest gridpoints value for ZLOC (limit of convection
!     note that put in values of '1' for Z b/c only doing 2-D interpolation
      ZZZ=1.0
!dwen(20090810)      CALL ADVINTZLOC(ZLOC(:,:,K1:K1),NXS,NYS,X2,Y2,ZZZ,                 &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADV3NTZLOC(ZLOC(:,:,K1:K1,kg),XX,YY,ZZZ,                 &
              var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090810)      CALL ADVINTZLOC(ZLOC(:,:,K2:K2),NXS,NYS,X2,Y2,ZZZ,                 &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
      CALL ADV3NTZLOC(ZLOC(:,:,K2:K2,kg),XX,YY,ZZZ,                 &
              var2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
      IF(TF.LE.0.5)METO%ZLOCNEXT=VAR1
      IF(TF.GT.0.5)METO%ZLOCNEXT=VAR2

! JCL:(6/28/02)interpolate air density to particle position
      ZTMP=METO%ZNDX
      IF(ZTMP.LT.1.0)ZTMP=1.0
                           !see comments in advmet CHG(09/16/03)
      ZTMPR=aINT(ZTMP+1.0)
!dwen(20090825)      ZTMPR=DMIN1(DMAX1(DBLE(1.0),ZTMPR),DBLE(NLVL))
      ZTMPR=aMIN1(aMAX1(1.0,ZTMPR),float(NLVL))
! CHG(09/16/03) for rams no spatial interpolation in density
      IF (.NOT. RAMSFLG) THEN
!dwen(20090807)        CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,X2,Y2,ZTMP,                 &
!dwen(20090807)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
        CALL ADV3NT(D(:,:,:,K1,kg),XX,YY,ZTMP,                 &
                dens1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090807)        CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,X2,Y2,ZTMP,                 &
!dwen(20090807)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
        CALL ADV3NT(D(:,:,:,K2,kg),XX,YY,ZTMP,                 &
                dens2,GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
      ELSE
!dwen(20090807)        CALL ADVINT(D(:,:,:,K1),NXS,NYS,NZM,XX,YY,ZTMPR,                &
!dwen(20090807)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS1)
        CALL ADV3NT(D(:,:,:,K1,kg),anint(XX),anint(YY),ZTMPR,                &
                dens1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)

!dwen(20090807)        CALL ADVINT(D(:,:,:,K2),NXS,NYS,NZM,XX,YY,ZTMPR,                &
!dwen(20090807)     &           GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,DENS2)
        CALL ADV3NT(D(:,:,:,K2,kg),anint(XX),anint(YY),ZTMPR,                &
               dens2, GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
      END IF
      METO%DENSLOCAL=(DENS2-DENS1)*TF+DENS1

!CCCCCCCCCCC EXTRACT TEMPORALLY-INTERPOLATED SIGMAW & TL ***AFTER*** ADVECTION
! JCL:(10/31/02) New way to interpolate vertical profiles of sigmaw & TL between analysis times by scaling
!         the two analysis profiles with the temporally interpolated mixed-layer height
! 1) Find nearest model level to move Zi to
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2. so replace X2 with 
!               XX, Y2 with YY and Z2 with ZZ
!dwen(20090820)      XX=DNINT(X2)
!dwen(20090820)      YY=DNINT(Y2)
! JCL:(10/31/02) extract ground height with new interpolation routine that deals with 2D arrays
!dwen(20090810)      CALL ADVINT2D(ZT,NXS,NYS,XX,YY,                                   &
!dwen(20090810)     &         GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
      CALL ADV2NT(ZT(:,:,kg),anint(XX),anint(YY),                                   &
              var1,GRID(kg,kt1)%GLOBAL,GRID(kg,kt1)%NX,GRID(kg,kt1)%NY)
      ZGRD=VAR1
! calculate true agl heights [m] of model levels:
      do kk=1,nlvl
         ZLVLS(kk)=(1.0-ZSG(kk))*(ZMDL-ZGRD)
      end do
! Force zmix to be at a model level if KMIX0 < 0 (as in stbanl,
! redone here because of horizontal interpolation of zi, zgrd)
      zmix=METO%ZMLNEXT
      if (KMIX0 .lt. 0) then
         ! Find closest model level above k=1
         zmixnew=zlvls(2)
         do kk=3,nlvl
            IF(ABS(ZLVLS(kK)-ZMIX) .LT. ABS(ZMIXNEW-ZMIX))THEN
               ZMIXNEW=ZLVLS(kK)
            END IF
         enddo
         METO%ZMLNEXT=zmixnew
      endif

! 2)  Find CLOSEST vertical level based on coordinate that's SCALED BY zi
! CHG(09/11/03) don't use this for RAMS
      IF(.NOT.RAMSFLG)THEN
         ZMIX=METO%ZMLNEXT
         KK1=1
         KK2=1
         DO KK=1,NLVL
           ! altitude of current timestep (scaled by local zi)
           ZNOWSC=ZLVLS(KK)/ZMIX
           IF (KK1 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK1)/ZMLK1-ZNOWSC).GT.                      &
                            ABS(ZLVLS(KK1+1)/ZMLK1-ZNOWSC))
                 KK1=KK1+1
                 IF (KK1 == NLVL) EXIT
              END DO
           END IF
           IF (KK2 < NLVL) THEN
              DO WHILE(ABS(ZLVLS(KK2)/ZMLK2-ZNOWSC).GT.                      &
                           ABS(ZLVLS(KK2+1)/ZMLK2-ZNOWSC))
                 KK2=KK2+1
                 IF (KK2 == NLVL) EXIT
              END DO
           END IF
!dwen(20090826)           METz%SIGW(KK)=(SIGWK2(KK2)-SIGWK1(KK1))*TF+SIGWK1(KK1)
!dwen(20090826)           METz%TL(KK)=(TLK2(KK2)-TLK1(KK1))*TF+TLK1(KK1)
           METz(kk)%SIGW=(SIGWK2(KK2)-SIGWK1(KK1))*TF+SIGWK1(KK1)
           METz(kk)%TL=(TLK2(KK2)-TLK1(KK1))*TF+TLK1(KK1)
         END DO
           !do this for RAMS
      ELSE
                     !DNINT for scalars
!dwen(20090820) X2,Y2 are used in STILT version, not used in HYSPLIT version. XX here is 
!               equivalent to X2 in STILT version, YY to Y2 and ZZ to Z2. so replace X2 with 
!               XX, Y2 with YY and Z2 with ZZ
!dwen(20090820)        XX=DNINT(X2)
                     !DNINT for scalars
!dwen(20090820)        YY=DNINT(Y2)
! CHG(3/3/2004) Now get time averaged sigw and tl from RAMS, stored at later time
!dwen(20090813)        KNEAR=K2
        DO KK=1,NLVL
          ZK=DBLE(KK)
!dwen(20090807)          CALL ADVINT(TL(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,             &
!dwen(20090807)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR1)
          CALL ADV3NT(TL(:,:,:,k2,kg),anint(XX),anint(YY),ZK,             &
                var1, GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)

!dwen(20090807)          CALL ADVINT(SIGWRAMS(:,:,:,KNEAR),NXS,NYS,NZM,XX,YY,ZK,       &
!dwen(20090807)     &            GRID(KGRID)%GLOBAL,GRID(KGRID)%NX,GRID(KGRID)%NY,VAR2)
          CALL ADV3NT(SIGWRAMS(:,:,:,K2,kg),anint(XX),anint(YY),ZK,       &
                var2, GRID(kg,kt2)%GLOBAL,GRID(kg,kt2)%NX,GRID(kg,kt2)%NY)
!dwen(20090826)          METz%TL(KK)=VAR1
          METz(kk)%TL=VAR1
!dwen(20090826)          METz%SIGW(KK)=VAR2
          METz(kk)%SIGW=VAR2
        END DO
            !of IF(.NOT.RAMSFLG)THEN ... ELSE
      ENDIF
!************************************************

! test if particle above the model top                  
  IF(TRAJ.AND.ZP.LT.ZSG(NLVL))THEN
!    terminate trajectories
     KG=0
     KRET=1
     WRITE(KF21,*)' NOTICE advpnt: trajectory above data domain'
     WRITE(KF21,*)' Trajectory sigma pt: ',ZP 
     WRITE(KF21,*)' Top of model domain: ',ZSG(NLVL) 
     RETURN
  ELSE
!    maintain particles (full reflection turbulence assumed)
     ZP=MAX(ZSG(NLVL),ZP)
  END IF

  KRET=0
  RETURN

!-------------------------------------------------------------------------------
! memory allocation errors
!-------------------------------------------------------------------------------

9000 WRITE(*,*)'*ERROR* advpnt: memory allocation - ',KRET,ECODE 
     STOP 900

END SUBROUTINE advpnt
