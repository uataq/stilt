!-------------------------------------------------------------------------------
! Configure subroutine interface argument lists stb routines 
! Last Revised: 02 Dec 2003 (RRD) - initial version extracted from sub prfcom
!               02 Jul 2004 (RRD) - stbanl argument list for hflx
!               11 May 2005 (RRD) - revised argument list stbsnd
!               08 Mar 2006 (RRD) - static stability parameter (stbanl)
!               14 Apr 2006 (RRD) - wrf interface issues 
!               21 Nov 2006 (RRD) - tke day night partition
!               04 Jun 2008 (RRD) - minimum mixing depth into stbanl
!               15 Aug 2008 (RRD) - more mixing options through stb routines
!-------------------------------------------------------------------------------

  INTERFACE
!-------------------------------------------------------------------------------

  SUBROUTINE PRFSIG(VMIX,TKEN,VELV,ZFLG,QFLG,UFLG,TFLG,PFLG,DZDT,OFFSET,     &
             ZSFC,P0,U0,V0,T0,Z0,NZ,PSG,P,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,       &
             PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX)
  IMPLICIT NONE
  LOGICAL,    INTENT(IN)    :: vmix     ! vertical mixing flag
  LOGICAL,    INTENT(IN)    :: tken     ! turbulent kinetic energy
  LOGICAL,    INTENT(IN)    :: velv     ! velocity variance
  LOGICAL,    INTENT(IN)    :: zflg     ! pressure available (nonhydrostatic)
  LOGICAL,    INTENT(IN)    :: qflg     ! specific humidity indicator
  LOGICAL,    INTENT(IN)    :: uflg     ! low level wind
  LOGICAL,    INTENT(IN)    :: tflg     ! low level temp 
  LOGICAL,    INTENT(IN)    :: pflg     ! surface pressure
  LOGICAL,    INTENT(IN)    :: dzdt     ! vertical velocity in height units
  REAL,       INTENT(IN)    :: offset   ! pressure offset
  REAL,       INTENT(INOUT) :: zsfc     ! terrain height (m)
  REAL,       INTENT(INOUT) :: p0       ! surface pressure at data terrain (mb)
  REAL,       INTENT(IN)    :: u0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: v0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: t0       ! low level temperaure (deg K)
  REAL,       INTENT(IN)    :: z0       ! roughness length (m)
  INTEGER,    INTENT(IN)    :: nz       ! number of input levels
  REAL,       INTENT(IN)    :: psg(:)   ! data sigma-p profile
  REAL,       INTENT(IN)    :: p  (:)   ! pressure data (non-hydrostatic)
  REAL,       INTENT(IN)    :: u  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: v  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: w  (:)   ! vertical motion component (dp/dt)
  REAL,       INTENT(IN)    :: t  (:)   ! temperature profile (deg K)
  REAL,       INTENT(IN)    :: q  (:)   ! specific humidity (kg/kg)
  REAL,       INTENT(IN)    :: e  (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(IN)    :: h  (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(IN)    :: x  (:)   ! w-component velocity var (m2/s2)
  INTEGER,    INTENT(IN)    :: nl       ! number of output sigma levels
  REAL,       INTENT(IN)    :: zmdl     ! internal model top (meters)
  REAL,       INTENT(IN)    :: zsg(:)   ! internal model output sigma levels
  REAL,       INTENT(OUT)   :: pp (:)   ! pressure at sigma level (mb)
  REAL,       INTENT(OUT)   :: uu (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: vv (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: ww (:)   ! vertical motion term (sigma/time)
  REAL,       INTENT(OUT)   :: tt (:)   ! virtual potential temperature (pot K)
  REAL,       INTENT(OUT)   :: zz (:)   ! internal model sigma height (meters)
  REAL,       INTENT(OUT)   :: rh (:)   ! relative humidity fraction (0-1)
  REAL,       INTENT(OUT)   :: den(:)   ! air density (kg/m3)
  REAL,       INTENT(OUT)   :: aa (:)   ! ambient temperature (deg K)
  REAL,       INTENT(OUT)   :: ee (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(OUT)   :: hh (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(OUT)   :: xx (:)   ! w-component velocity var (m2/s2)
  END SUBROUTINE prfsig
!-------------------------------------------------------------------------------
!dwen(20090822)  SUBROUTINE PRFTER(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,ZMDLT,    &
!dwen(20090822)                    P0,U0,V0,T0,Z0,NZ,PSG,P,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,   &
!dwen(20090822)                    PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX)

 SUBROUTINE PRFTER(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,ZMDLT,    &
                   P0,U0,V0,T0,Z0,NZ,PSG,P,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,   &
                   PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX,ramsflg)



  IMPLICIT NONE
  LOGICAL,    INTENT(IN)    :: vmix     ! vertical mixing flag
  LOGICAL,    INTENT(IN)    :: tken     ! turbulent kinetic energy
  LOGICAL,    INTENT(IN)    :: velv     ! velocity variance
  LOGICAL,    INTENT(IN)    :: qflg     ! specific humidity indicator
  LOGICAL,    INTENT(IN)    :: uflg     ! low level wind
  LOGICAL,    INTENT(IN)    :: tflg     ! low level temp 
  LOGICAL,    INTENT(IN)    :: pflg     ! surface pressure
  LOGICAL,    INTENT(IN)    :: sflg     ! surface terrain 
  REAL,       INTENT(INOUT) :: zsfc     ! terrain height (m)
  REAL,       INTENT(IN)    :: zmdlt    ! scaling height for Z* coordinates
  REAL,       INTENT(INOUT) :: p0       ! surface pressure at data terrain (mb)
  REAL,       INTENT(IN)    :: u0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: v0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: t0       ! low level temperaure (deg K)
  REAL,       INTENT(IN)    :: z0       ! roughness length (m)
  INTEGER,    INTENT(IN)    :: nz       ! number of input levels
  REAL,       INTENT(IN)    :: psg(:)   ! data sigma-p profile
  REAL,       INTENT(IN)    :: p  (:)   ! pressure data (non-hydrostatic)
  REAL,       INTENT(IN)    :: u  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: v  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: w  (:)   ! vertical motion component (dp/dt)
  REAL,       INTENT(IN)    :: t  (:)   ! temperature profile (deg K)
  REAL,       INTENT(IN)    :: q  (:)   ! specific humidity (kg/kg)
  REAL,       INTENT(IN)    :: e  (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(IN)    :: h  (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(IN)    :: x  (:)   ! w-component velocity var (m2/s2)
  INTEGER,    INTENT(IN)    :: nl       ! number of output sigma levels
  REAL,       INTENT(IN)    :: zmdl     ! internal model top (meters)
  REAL,       INTENT(IN)    :: zsg(:)   ! internal model output sigma levels
  REAL,       INTENT(OUT)   :: pp (:)   ! pressure at sigma level (mb)
  REAL,       INTENT(OUT)   :: uu (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: vv (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: ww (:)   ! vertical motion term (sigma/time)
  REAL,       INTENT(OUT)   :: tt (:)   ! virtual potential temperature (pot K)
  REAL,       INTENT(OUT)   :: zz (:)   ! internal model sigma height (meters)
  REAL,       INTENT(OUT)   :: rh (:)   ! relative humidity fraction (0-1)
  REAL,       INTENT(OUT)   :: den(:)   ! air density (kg/m3)
  REAL,       INTENT(OUT)   :: aa (:)   ! ambient temperature (deg K)
  REAL,       INTENT(OUT)   :: ee (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(OUT)   :: hh (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(OUT)   :: xx (:)   ! w-component velocity var (m2/s2)
!dwen(20090805) *********************
! JCL(03/27/03): set flag: whether data from RAMS or not 
 logical,     intent(in)    :: ramsflg
!dwen *********************************
  END SUBROUTINE prfter
!-------------------------------------------------------------------------------
  SUBROUTINE PRFPRS(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,          &
                    P0,U0,V0,T0,Z0,NZ,PSG,Z,U,V,W,T,R,E,H,X,NL,ZMDL,ZSG,   &
                    PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX)
  IMPLICIT NONE
  LOGICAL,    INTENT(IN)    :: vmix     ! vertical mixing flag
  LOGICAL,    INTENT(IN)    :: tken     ! turbulent kinetic energy
  LOGICAL,    INTENT(IN)    :: velv     ! velocity variance
  LOGICAL,    INTENT(IN)    :: qflg     ! specific humidity indicator
  LOGICAL,    INTENT(IN)    :: uflg     ! low level wind
  LOGICAL,    INTENT(IN)    :: tflg     ! low level temp 
  LOGICAL,    INTENT(IN)    :: pflg     ! surface pressure
  LOGICAL,    INTENT(IN)    :: sflg     ! surface terrain 
  REAL,       INTENT(INOUT) :: zsfc     ! terrain height (m)
  REAL,       INTENT(INOUT) :: p0       ! surface pressure at data terrain (mb)
  REAL,       INTENT(IN)    :: u0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: v0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: t0       ! low level temperaure (deg K)
  REAL,       INTENT(IN)    :: z0       ! roughness length (m)
  INTEGER,    INTENT(IN)    :: nz       ! number of input levels
  REAL,       INTENT(IN)    :: psg(:)   ! data sigma-p profile
  REAL,       INTENT(IN)    :: z  (:)   ! pressure data (non-hydrostatic)
  REAL,       INTENT(IN)    :: u  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: v  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: w  (:)   ! vertical motion component (dp/dt)
  REAL,       INTENT(IN)    :: t  (:)   ! temperature profile (deg K)
  REAL,       INTENT(IN)    :: r  (:)   ! specific humidity (kg/kg)
  REAL,       INTENT(IN)    :: e  (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(IN)    :: h  (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(IN)    :: x  (:)   ! w-component velocity var (m2/s2)
  INTEGER,    INTENT(IN)    :: nl       ! number of output sigma levels
  REAL,       INTENT(IN)    :: zmdl     ! internal model top (meters)
  REAL,       INTENT(IN)    :: zsg(:)   ! internal model output sigma levels
  REAL,       INTENT(OUT)   :: pp (:)   ! pressure at sigma level (mb)
  REAL,       INTENT(OUT)   :: uu (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: vv (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: ww (:)   ! vertical motion term (sigma/time)
  REAL,       INTENT(OUT)   :: tt (:)   ! virtual potential temperature (pot K)
  REAL,       INTENT(OUT)   :: zz (:)   ! internal model sigma height (meters)
  REAL,       INTENT(OUT)   :: rh (:)   ! relative humidity fraction (0-1)
  REAL,       INTENT(OUT)   :: den(:)   ! air density (kg/m3)
  REAL,       INTENT(OUT)   :: aa (:)   ! ambient temperature (deg K)
  REAL,       INTENT(OUT)   :: ee (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(OUT)   :: hh (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(OUT)   :: xx (:)   ! w-component velocity var (m2/s2)
  END SUBROUTINE prfprs
!-------------------------------------------------------------------------------
  SUBROUTINE PRFECM(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PFLG,SFLG,ZSFC,           &
                    P0,U0,V0,T0,Z0,NZ,PSG,ESG,U,V,W,T,Q,E,H,X,NL,ZMDL,ZSG,  &
                    PP,UU,VV,WW,TT,ZZ,RH,DEN,AA,EE,HH,XX)
  IMPLICIT NONE
  LOGICAL,    INTENT(IN)    :: vmix     ! vertical mixing flag
  LOGICAL,    INTENT(IN)    :: tken     ! turbulent kinetic energy
  LOGICAL,    INTENT(IN)    :: velv     ! velocity variance
  LOGICAL,    INTENT(IN)    :: qflg     ! specific humidity indicator
  LOGICAL,    INTENT(IN)    :: uflg     ! low level wind
  LOGICAL,    INTENT(IN)    :: tflg     ! low level temp 
  LOGICAL,    INTENT(IN)    :: pflg     ! surface pressure
  LOGICAL,    INTENT(IN)    :: sflg     ! surface terrain 
  REAL,       INTENT(INOUT) :: zsfc     ! terrain height (m)
  REAL,       INTENT(INOUT) :: p0       ! surface pressure at data terrain (mb)
  REAL,       INTENT(IN)    :: u0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: v0       ! low level horizontal wind component
  REAL,       INTENT(IN)    :: t0       ! low level temperaure (deg K)
  REAL,       INTENT(IN)    :: z0       ! roughness length (m)
  INTEGER,    INTENT(IN)    :: nz       ! number of input levels
  REAL,       INTENT(IN)    :: psg(:)   ! data sigma-p profile
  REAL,       INTENT(OUT)   :: esg(:)   ! data sigma-p profile (sfc pres corr)
  REAL,       INTENT(IN)    :: u  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: v  (:)   ! horizontal wind component
  REAL,       INTENT(IN)    :: w  (:)   ! vertical motion component (dp/dt)
  REAL,       INTENT(IN)    :: t  (:)   ! temperature profile (deg K)
  REAL,       INTENT(IN)    :: q  (:)   ! specific humidity (kg/kg)
  REAL,       INTENT(IN)    :: e  (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(IN)    :: h  (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(IN)    :: x  (:)   ! w-component velocity var (m2/s2)
  INTEGER,    INTENT(IN)    :: nl       ! number of output sigma levels
  REAL,       INTENT(IN)    :: zmdl     ! internal model top (meters)
  REAL,       INTENT(IN)    :: zsg(:)   ! internal model output sigma levels
  REAL,       INTENT(OUT)   :: pp (:)   ! pressure at sigma level (mb)
  REAL,       INTENT(OUT)   :: uu (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: vv (:)   ! horizontal wind component
  REAL,       INTENT(OUT)   :: ww (:)   ! vertical motion term (sigma/time)
  REAL,       INTENT(OUT)   :: tt (:)   ! virtual potential temperature (pot K)
  REAL,       INTENT(OUT)   :: zz (:)   ! internal model sigma height (meters)
  REAL,       INTENT(OUT)   :: rh (:)   ! relative humidity fraction (0-1)
  REAL,       INTENT(OUT)   :: den(:)   ! air density (kg/m3)
  REAL,       INTENT(OUT)   :: aa (:)   ! ambient temperature (deg K)
  REAL,       INTENT(OUT)   :: ee (:)   ! turb kin energy or w-var (m2/s2)
  REAL,       INTENT(OUT)   :: hh (:)   ! u-component velocity var (m2/s2)
  REAL,       INTENT(OUT)   :: xx (:)   ! w-component velocity var (m2/s2)
  END SUBROUTINE prfecm
!-------------------------------------------------------------------------------
!dwen(20090822)  SUBROUTINE STBANL(KBLS,I,J,KS,KMIX0,KMIXD,MIXD,TKEN,HFLX,EFLX,UFLX,USTAR,  &
!dwen(20090822)                    TSTAR,Z0,FMU,FMV,FHS,NL,UU,VV,TT,ZZ,EE,DEN,PBLH,PSIR,SSP)

subroutine STBANL(KBLS,I,J,KS,KMIX0,KMIXD,MIXD,TKEN,HFLX,EFLX,UFLX,USTAR,   &
                 TSTAR,Z0,FMU,FMV,FHS,NL,UU,VV,TT,ZZ,EE,DEN,PBLH,PSIR,SSP,  &
                 rh,zloc,iconvect,ziscale)



  IMPLICIT NONE
  INTEGER,   INTENT(IN)    :: kbls         ! stability compuation method
  INTEGER,   INTENT(IN)    :: i,j          ! horizontal grid point proccesed
  INTEGER,   INTENT(IN)    :: ks           ! index top of surface layer
  INTEGER,   INTENT(IN)    :: kmix0        ! minimum mixing depth
  INTEGER,   INTENT(IN)    :: kmixd        ! mixed layer depth options
  LOGICAL,   INTENT(IN)    :: mixd         ! mixed layer depth in data file
  LOGICAL,   INTENT(IN)    :: tken         ! tke field available         
  LOGICAL,   INTENT(IN)    :: hflx         ! fluxes available (heat)     
  LOGICAL,   INTENT(IN)    :: eflx         ! momentum flux as scalar exchange
  LOGICAL,   INTENT(IN)    :: uflx         ! fluxes available (momentum)
  LOGICAL,   INTENT(IN)    :: ustar,tstar  ! friction veloc and temp available
  REAL,      INTENT(IN)    :: z0           ! roughness length (m)
  REAL,      INTENT(INOUT) :: fmu          ! u-component momentum flux (N/m2) 
  REAL,      INTENT(INOUT) :: fmv          ! v-component momentum flux (N/m2) 
  REAL,      INTENT(IN)    :: fhs          ! sensible heat flux (W/m2)
  INTEGER,   INTENT(IN)    :: nl           ! number of output sigma levels
  REAL,      INTENT(IN)    :: uu  (:)      ! horizontal wind component (m/s)
  REAL,      INTENT(IN)    :: vv  (:)      ! horizontal wind component (m/s)
  REAL,      INTENT(IN)    :: tt  (:)      ! virtual potential temp (deg K)
  REAL,      INTENT(INOUT) :: zz  (:)      ! height at levels (m)
  REAL,      INTENT(IN)    :: ee  (:)      ! turbulent kinetic energy (J/kg)
  REAL,      INTENT(IN)    :: den (:)      ! air density (kg/m3)
  REAL,      INTENT(INOUT) :: pblh         ! mixed layer depth (m)
  REAL,      INTENT(OUT)   :: psir         ! integrated stability function heat
  REAL,      INTENT(OUT)   :: ssp          ! static stability parameter
!dwen(20090805) ****************
! CHG:(11/19/01) add humidity as input argument, to allow calculation of LCL, LFC and LOC
  real,      intent(in)    :: rh(:)
! CHG:(12/04/01) add ZLOC (limit of convection) as output argument
  real,      intent(out)   :: zloc
! CHG:(9/17/02) add 'ICONVECT' as convection flag 
  integer,   intent(in)    :: iconvect
  real,   intent(in)    :: ziscale
!dwen ******************************
  END SUBROUTINE stbanl
!-------------------------------------------------------------------------------
!dwen(20090822)  SUBROUTINE STBSND(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,EE,HH,XX)
SUBROUTINE STBSND(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,EE,HH,XX,  &
                  xv,ttl,sigmaw,z0)

  IMPLICIT NONE
  REAL,      INTENT(IN)    :: tkerd     ! day ratio of vert to horiz turb
  REAL,      INTENT(IN)    :: tkern     ! noght ratio of vert to horiz turb
  INTEGER,   INTENT(IN)    :: kzmix     ! height averaged mixing   
  REAL,      INTENT(IN)    :: tvmix     ! tropospheric mixing scale factor
  INTEGER,   INTENT(IN)    :: ksfc      ! index of top of surface layer
  INTEGER,   INTENT(IN)    :: nl        ! number of sigma levels
  REAL,      INTENT(IN)    :: uu (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: vv (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: tt (:)    ! virtual potential temperature (pot K)
  REAL,      INTENT(IN)    :: zz (:)    ! height at levels (m)
  REAL,      INTENT(INOUT) :: ee (:)    ! in TKE (J/kg), output  - V'2 (m2/s2)
  REAL,      INTENT(INOUT) :: hh (:)    ! U component turbulence - U'2 (m2/s2)
  REAL,      INTENT(INOUT) :: xx (:)    ! vertical turbulence    - W'2 (m2/s2)
!dwen(20090806) ***********************
!        add xv as vertical mixing coefficient
 real,      intent(out)    :: xv(:)

 real,      intent(out)    :: ttl(:)
 real,      intent(out)    :: sigmaw(:)
 real,      intent(in)     :: z0
!dwen **********************************
  END SUBROUTINE stbsnd
!-------------------------------------------------------------------------------
  SUBROUTINE STBVAR(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,TT,ZZ,EE,HH,XX)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: tkerd     ! day ratio of vert to horiz turb
  REAL,      INTENT(IN)    :: tkern     ! noght ratio of vert to horiz turb
  INTEGER,   INTENT(IN)    :: kzmix     ! height averaged mixing   
  REAL,      INTENT(IN)    :: tvmix     ! tropospheric mixing scale factor
  INTEGER,   INTENT(IN)    :: ksfc      ! index of top of surface layer
  INTEGER,   INTENT(IN)    :: nl        ! number of sigma levels
  REAL,      INTENT(IN)    :: uu (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: vv (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: tt (:)    ! virtual potential temperature (pot K)
  REAL,      INTENT(IN)    :: zz (:)    ! height at levels (m)
  REAL,      INTENT(INOUT) :: ee (:)    ! in TKE (J/kg), output  - V'2 (m2/s2)
  REAL,      INTENT(INOUT) :: hh (:)    ! U component turbulence - U'2 (m2/s2)
  REAL,      INTENT(INOUT) :: xx (:)    ! vertical turbulence    - W'2 (m2/s2)
  END SUBROUTINE stbvar
!-------------------------------------------------------------------------------
  SUBROUTINE STBTKE(TKERD,TKERN,KZMIX,TVMIX,KSFC,NL,UU,VV,ZZ,EE,HH,XX)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: tkerd     ! day ratio of vert to horiz turb
  REAL,      INTENT(IN)    :: tkern     ! night ratio of vert to horiz turb
  INTEGER,   INTENT(IN)    :: kzmix     ! height averaged mixing   
  REAL,      INTENT(IN)    :: tvmix     ! tropospheric mixing scale factor
  INTEGER,   INTENT(IN)    :: ksfc      ! index of top of surface layer
  INTEGER,   INTENT(IN)    :: nl        ! number of sigma levels
  REAL,      INTENT(IN)    :: uu (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: vv (:)    ! horizontal wind component
  REAL,      INTENT(IN)    :: zz (:)    ! height at levels (m)
  REAL,      INTENT(INOUT) :: ee (:)    ! in TKE (J/kg), output  - V'2 (m2/s2)
  REAL,      INTENT(INOUT) :: hh (:)    ! U component turbulence - U'2 (m2/s2)
  REAL,      INTENT(INOUT) :: xx (:)    ! vertical turbulence    - W'2 (m2/s2)
  END SUBROUTINE stbtke
!-------------------------------------------------------------------------------
!dwen(20090822)  SUBROUTINE STBHOR(NXS,NYS,NLVL,GX,GY,U,V,H,E)
SUBROUTINE STBHOR(NXS,NYS,NLVL,GX,GY,U,V,H,hm,E)
  IMPLICIT NONE
  INTEGER,     INTENT(IN)    :: nxs,nys         ! horizontal subgrid dimensions
  INTEGER,     INTENT(IN)    :: nlvl            ! number of output levels
  REAL,        INTENT(IN)    :: gx(:,:)         ! horizontal grid spacing
  REAL,        INTENT(IN)    :: gy(:,:)         ! horizontal grid spacing
  REAL,        INTENT(IN)    :: u (:,:,:)       ! horizontal wind component
  REAL,        INTENT(IN)    :: v (:,:,:)       ! horizontal wind component
  REAL,        INTENT(OUT)   :: h (:,:,:)       ! u component turbulence
  REAL,        INTENT(OUT)   :: e (:,:,:)       ! v component turbulence
!dwen(20090806) ***********************
!              add hm as horizontal mixing coefficient
  real,        intent(out)   :: hm(:,:,:)
!dwen           ***********************
  END SUBROUTINE stbhor

!-------------------------------------------------------------------------------
  END INTERFACE 
