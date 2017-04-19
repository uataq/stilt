!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METINI           METeorological INItialization
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL INITIALIZATION OPENS THE METEO FILES THE FIRST TIME
!   ON EACH DEFINED GRID SYSTEM; DEFINES DEFAULT GRID STRUCTURE WHICH
!   CANNOT CHANGE WITH TIME FOR A DEFINED GRID NUMBER. MULTIPLE FILE
!   TIME ARE FOR THE SAME GRID NUMBER AND REQUIRE THE SAME STRUCTURE
!   DEFINITIONS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 31 Mar 1998 (RRD)
!                  15 Apr 1999 (RRD) - added surface height definition
!                  28 Jul 2000 (RRD) - flux test correction
!                  14 Aug 2000 (RRD) - defined 1 hour precip / zflg
!                  04 Sep 2000 (RRD) - fortran90 upgrade
!                  09 Mar 2001 (RRD) - option to handle global grids
!                  01 Oct 2001 (RRD) - simultaneous multiple meteorology
!                  05 Dec 2001 (RRD) - additional surface terrain name
!                  26 Feb 2002 (RRD) - downward shortwave radiation
!                  23 Jul 2002 (RRD) - test for eta correction option
!                  03 Mar 2003 (RRD) - ichar test
!                  07 Aug 2003 (RRD) - sync test for lat/lon grids
!                  16 Sep 2003 (RRD) - do not set uflx with ngm data
!                  17 Oct 2003 (RRD) - added turbulent kinetic energy
!                  10 Nov 2003 (RRD) - check for velocity variances
!                  02 Apr 2004 (RRD) - generic file unit numbers
!                  03 Dec 2004 (BS)  - precipitation rate (PRT6)
!                  20 Dec 2004 (RRD) - disable checksum test for 2m/10m  
!                  14 Apr 2006 (RRD) - test for wrf versions
!                  22 May 2006 (RRD) - added test for mixed layer depth
!                  22 Mar 2007 (RRD) - sort same grid by forecast time
!                  05 Mar 2008 (RRD) - test level dependent variables TKE
!                  10 Jun 2008 (RRD) - starting minutes and isot option
!                  15 Aug 2008 (RRD) - split horizontal and vertical mixing
!                  14 Oct 2008 (RRD) - some warnings to message file
!                  10 Dec 2008 (RRD) - common code for unregistered version
!
! USAGE:  CALL METINI(KDEF,KZMIX,TVMIX,KBLS,KBLT,NGRD,NTIM,OLAT,IBYR,
!                     IBMO,IBDA,IBHR,IBMN,BACK,KVEL,KFOR)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE METINI(KDEF,KZMIX,TVMIX,KBLS,KBLT,NGRD,NTIM,OLAT,IBYR,    &
     IBMO,IBDA,IBHR,IBMN,BACK,KVEL,KFOR)

  USE funits
  use module_defgrid ! meteorology grid and file
  use map_utils
  IMPLICIT NONE

  !-------------------------------------------------------------------------------
  ! argument list variables
  !-------------------------------------------------------------------------------

  INTEGER,   INTENT(IN)    :: kdef                 ! deformation horizontal mix
  INTEGER,   INTENT(IN)    :: kzmix                ! vertical mixing averaging
  REAL,      INTENT(IN)    :: tvmix                ! tropospheric mixing scale
  INTEGER,   INTENT(IN)    :: kbls                 ! pbl stability derived from
  INTEGER,   INTENT(IN)    :: kblt                 ! pbl scheme for mixing
  INTEGER,   INTENT(INOUT) :: ngrd                 ! number of meteo grids 
  INTEGER,   INTENT(INOUT) :: ntim                 ! number of meteo times 
  REAL,      INTENT(IN)    :: olat                 ! source latitude
  INTEGER,   INTENT(INOUT) :: ibyr,ibmo,ibda       ! calculation starting date
  INTEGER,   INTENT(INOUT) :: ibhr,ibmn            ! calculation starting time
  LOGICAL,   INTENT(IN)    :: back                 ! integration direction 
  INTEGER,   INTENT(INOUT) :: kvel                 ! vertical motion method
  INTEGER,   INTENT(IN)    :: kfor                 ! 0=NoForecast 1=AllData 

  !-------------------------------------------------------------------------------
  ! internal variables
  !-------------------------------------------------------------------------------

  LOGICAL    :: fields
  INTEGER    :: i,j,nvar,kg,kt,k2,nlvl,macc,mjet,itda,ithr,itmn
  REAL       :: SUX,SUY,SVX,SVY

  !dwen(20090818) ******************************
  !  Need local real*4 variables for interface with WPS modules:
  REAL              :: dxm, true1, true2, sync_lat, synclon180, &
       reflon180, sync_xp, sync_yp
  REAL, ALLOCATABLE :: xlat_array(:,:)
  INTEGER           :: iproj, i_xlat, j_xlat
  REAL              :: xi, xj, xlon, xlat
  !*******************************************

  !-------------------------------------------------------------------------------
  ! external variables
  !-------------------------------------------------------------------------------

  COMMON /STAGGER/ SUX,SUY,SVX,SVY

  !-------------------------------------------------------------------------------
  ! Some Fortran compilers the ICHAR function is undefined for values greater
  ! Generic Fortran and C code functions for conversions can be used (ncep mova2i)
  ! In Fortran:  JCHAR=ISHFT(TRANSFER(MYCHR,JCHAR),8-BIT_SIZE(JCHAR))
  ! In C:        INT MOVA2I(UNSIGNED CHAR *A) { RETURN (INT) (*A); }
  ! If ICHAR always returns a signed integer, then ICHAR can be replaced with:
  !-------------------------------------------------------------------------------
  CHARACTER(1)  :: mychr
  INTEGER       :: jchar
  JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
  !-------------------------------------------------------------------------------

  INTERFACE
     !-------------------------------------------------------------------------------
     SUBROUTINE TM2DAY(MACM,IY,IM,ID,IH,MN)
       IMPLICIT NONE
       INTEGER,   INTENT(IN)    :: macm            ! accumulate minutes
       INTEGER,   INTENT(OUT)   :: iy,im,id,ih,mn  ! current date/time
     END SUBROUTINE tm2day
     !-------------------------------------------------------------------------------
     SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
       IMPLICIT NONE
       INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
       INTEGER,  INTENT(OUT)  :: macc                 ! minutes since 1 Jan 1970
     END SUBROUTINE tm2min
     !-------------------------------------------------------------------------------
     SUBROUTINE METSET(KG,KT,KFOR)
       IMPLICIT NONE
       INTEGER,       INTENT(IN)    :: kg           ! grid identfication number
       INTEGER,       INTENT(IN)    :: kt           ! data time period (1 or 2)
       INTEGER,       INTENT(IN)    :: kfor         ! 0=NoForecasts 1=AllData    
     END SUBROUTINE metset
     !-------------------------------------------------------------------------------
     SUBROUTINE GBLSET(KG,KT)
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: kg
       INTEGER,  INTENT(IN) :: kt
     END SUBROUTINE gblset
     !-------------------------------------------------------------------------------
  END INTERFACE

  ! set all grid numbers to missing
  DO KG=1,MGRD
     DO KT=1,MTIM
        GRID(KG,KT)%NUMBER=-99
     END DO
  END DO

  DO KG=1,NGRD
     DO KT=1,NTIM

        !-------------------------------------------------------------------------------

        MYCHR=CHAR(128)
        IF(JCHAR(MYCHR).NE.128)THEN
           WRITE(*,*)'*ERROR*: metini - ASCII character conversion not working'
           WRITE(*,*)'Contact support - modifications required to subroutine PAKINP'
           STOP 900
        END IF

        !-------------------------------------------------------------------------------
        ! file parameters
        !-------------------------------------------------------------------------------

        ! record counter offset usually set to zero (see metold)
        ! however old style NH/SH data sets require the offset to
        ! position to the SH data if the starting lat indicates SH

        DREC(KG,KT)%OFFSET=0
        IF(OLAT.LT.0.0)DREC(KG,KT)%OFFSET=INT(OLAT)

        ! define unit numbers for each meteorology grid     
        FILE(KG,KT)%KUNIT=KF01+(KG-1)*NTIM+(KT-1)
        IF(FILE(KG,KT)%KUNIT.GT.KF02)THEN 
           WRITE(*,*)'*ERROR* metini: exceeding maximum number of meteorology files'
           WRITE(*,*)'Number permitted: ',(KF02-KF01+1)
           STOP 900
        END IF
        ! initializes time period 
        !dwen(20090817): HYSPLIT4.5 removed HEADER from argument list
        !dwen(20090817): HYSPLIT4.9 added KFOR for forecast data option
        CALL METSET(KG,KT,KFOR)

        !-------------------------------------------------------------------------------
        ! surface parameters
        ! set flags for the type of surface data available
        !-------------------------------------------------------------------------------

        ! ten meter wind and two meter temp
        DREC(KG,KT)%UFLG=.FALSE.
        DREC(KG,KT)%TFLG=.FALSE.

        ! surface exchange coefficient, heat, momentum, shortwave flux
        DREC(KG,KT)%EFLX=.FALSE.
        DREC(KG,KT)%HFLX=.FALSE.
        DREC(KG,KT)%UFLX=.FALSE.
        DREC(KG,KT)%DSWF=.FALSE.

        ! friction velocity and temperature
        DREC(KG,KT)%USTR=.FALSE.
        DREC(KG,KT)%TSTR=.FALSE.

        ! surface terrain height and pressure
        DREC(KG,KT)%SHGT=.FALSE.
        DREC(KG,KT)%PRSS=.FALSE.

        ! set default precipitation accumulation time (default = none)
        DREC(KG,KT)%ACYCLE=0

        !*************************************************
        !dwen(20090316):adopted from metini of STILT, add kt
        ! CHG:(20/11/01) set default conv. prec. flag
        DREC(KG,kt)%CFLG = .FALSE.
        ! CHG:(20/11/01) set default total cloud flag
        DREC(KG,kt)%TCLF = .FALSE.
        ! CHG:(12/04/01) set default low cloud flag
        DREC(KG,kt)%LCLF = .FALSE.
        ! CHG:(20/11/01) set default downw. rad. flag
        DREC(KG,kt)%RADF = .FALSE.
        ! CHG:(22/01/03) set default soil moisture flag
        DREC(KG,kt)%SLMF = .FALSE.
        !*************************************************

        ! input vertical velocity units m/s rather than hPa/s
        DREC(KG,KT)%DZDT=.FALSE.

        ! mixed layer depth available 
        DREC(KG,KT)%MIXD=.FALSE.

        ! all WRF versions have u,v grids staggered with respect to mass points
        ! WRF=1 - all versions assume instantaneous velocities (from arw2arl)
        ! WRF=2 - not used; nearest neighbor interpolation (NMM-WRF from wrf2arl) 
        ! WRF=3 - associated with AER-WRF for STILT model (averaged fields)

        IF(GRID(KG,KT)%MODEL_ID.EQ.'AWRF') THEN
           DREC(KG,KT)%WRF=1  
           !    defines u,v points relative to mass points
           SUX=0.5
           SUY=0.0
           SVX=0.0
           SVY=0.5
        ELSE
           DREC(KG,KT)%WRF=0    
           SUX=0.0
           SUY=0.0
           SVX=0.0
           SVY=0.0
        END IF

        ! scan index record settings for certain variable ids
        ! if one component there assume the other is as well

        NVAR=DREC(KG,KT)%NUM_VARB(1)
        DO I=1,NVAR

           !    mixed layer depth available
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'MXHT'.OR.     & 
                DREC(KG,KT)%VARB_ID(I,1).EQ.'HPBL'.OR.     &
                DREC(KG,KT)%VARB_ID(I,1).EQ.'PBLH')        &
                DREC(KG,KT)%MIXD=.TRUE.

           !    when wrf stagger pressure fields are available then WRF=3
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'MUU0'.OR.     & 
                DREC(KG,KT)%VARB_ID(I,1).EQ.'MUV0') DREC(KG,KT)%WRF=3
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'DZDT') DREC(KG,KT)%DZDT=.TRUE.

           !    surface height (shgt,hgts) or pressure (prss) is available
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'SHGT'.AND.    &
                DREC(KG,KT)%CHK_SUM(I,1).GE.0)             &
                DREC(KG,KT)%SHGT=.TRUE.
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'HGTS'.AND.    &
                DREC(KG,KT)%CHK_SUM(I,1).GE.0)             &
                DREC(KG,KT)%SHGT=.TRUE.
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'PRSS'.AND.    &
                DREC(KG,KT)%CHK_SUM(I,1).GE.0)             &
                DREC(KG,KT)%PRSS=.TRUE.

           !    two meter temperatures 
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TMPS'.OR.     & 
                DREC(KG,KT)%VARB_ID(I,1).EQ.'T02M')THEN
              DREC(KG,KT)%TFLG=.TRUE.
              !       IF(DREC(KG,KT)%CHK_SUM(I,1).EQ.0)          &
              WRITE(KF21,*)'WARNING: metini - initial 2 m TEMP missing' 
           END IF

           !    ten meter winds 
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'U10M')THEN
              DREC(KG,KT)%UFLG=.TRUE.
              !       IF(DREC(KG,KT)%CHK_SUM(I,1).EQ.0)          &
              WRITE(KF21,*)'WARNING: metini - initial 10 m WIND missing' 
           END IF

           !    shortwave flux 
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'DSWF')DREC(KG,KT)%DSWF=.TRUE.

           !    momentum flux (checksum positive for valid)
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'UMOF'.AND.    &
                DREC(KG,KT)%CHK_SUM(I,1).GE.0)             &
                DREC(KG,KT)%UFLX=.TRUE.

           !    heat flux (checksum positive for valid)
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'SHTF'.AND.    &
                DREC(KG,KT)%CHK_SUM(I,1).GE.0)             &
                DREC(KG,KT)%HFLX=.TRUE.
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'HFLX'.AND.    &
                DREC(KG,KT)%CHK_SUM(I,1).GE.0)             &
                DREC(KG,KT)%HFLX=.TRUE.

           !    special scalar momentum (exchange coefficient) from NGM
           !dwen(20090316):stilt verison of metini includes DREC(KG)%UFLX = .TRUE.
           !               uflx is no longer set for NGM data
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'EXCO')        &
                DREC(KG,KT)%EFLX=.TRUE.

           !    check for normalized momentum and temperature profiles
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'USTR')DREC(KG,KT)%USTR=.TRUE.
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TSTR')DREC(KG,KT)%TSTR=.TRUE.

           !    cycle time that the precip bucket is emptied (min)
           IF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TPPA'.OR.   &
                DREC(KG,KT)%VARB_ID(I,1).EQ.'TPPS')THEN
              !       precip accumulation or summed over entire data set
              DREC(KG,KT)%ACYCLE=99999
           ELSEIF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TPPD')THEN
              !       accumulation daily over 24 hours
              DREC(KG,KT)%ACYCLE=1440
           ELSEIF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TPPT')THEN
              !       accumulation over 12 hours
              DREC(KG,KT)%ACYCLE=720
           ELSEIF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TPP6')THEN
              !       accumulation over 6 hours
              DREC(KG,KT)%ACYCLE=360
           ELSEIF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TPP3')THEN
              !       accumlation over 3 hours
              DREC(KG,KT)%ACYCLE=180
           ELSEIF(DREC(KG,KT)%VARB_ID(I,1).EQ.'TPP1')THEN
              !       accumlation over 3 hours
              DREC(KG,KT)%ACYCLE=60
           ELSEIF(DREC(KG,KT)%VARB_ID(I,1).EQ.'PRT6')THEN
              !       precipitation rate over 6 hours
              DREC(KG,KT)%ACYCLE=-360			! negative flags rate
           END IF
           !****************************************************
           !dwen(20090315):copy the following lines from STILT, add KT
           ! convective precip available?
           IF (ANY(DREC(KG,kt)%VARB_ID(I,1) == (/'CPPT','CPPD','CPP3','CPP6','CPRC'/))) &
                DREC(KG,kt)%CFLG = .TRUE.

           ! CHG:(11/20/01) total cloud cover available
           IF (DREC(KG,kt)%VARB_ID(I,1) == 'TCLD') DREC(KG,kt)%TCLF = .TRUE.

           ! CHG:(12/04/01) low cloud cover available
           IF (DREC(KG,kt)%VARB_ID(I,1) == 'LCLD') DREC(KG,kt)%LCLF = .TRUE.

           ! CHG:(11/20/01) shortwave radiative flux available
           IF (DREC(KG,kt)%VARB_ID(I,1) == 'DSWF') DREC(KG,kt)%RADF = .TRUE.

           ! CHG:(11/20/01) soil moisture available
           IF (DREC(KG,kt)%VARB_ID(I,1) == 'SOLW') DREC(KG,kt)%SLMF = .TRUE.
           !*****************************************************

        END DO

        !-------------------------------------------------------------------------------
        ! upper level parameters (start at level 2)
        !-------------------------------------------------------------------------------

        ! some 3D variables may not be available at all levels
        DREC(KG,KT)%TKEN=.FALSE.  
        DREC(KG,KT)%QFLG=.FALSE.
        DREC(KG,KT)%VELV=.FALSE.
        DREC(KG,KT)%ZFLG=.FALSE.

        NLVL=GRID(KG,KT)%NZ
        DO J=2,NLVL

           !    default flags for other meteorological data by level
           DREC(KG,KT)%RFLG(J-1)=.FALSE.
           DREC(KG,KT)%WFLG(J-1)=.FALSE.

           !    shift index down to corresspond with upper level variables
           NVAR=DREC(KG,KT)%NUM_VARB(J)
           DO I=1,NVAR
              !       upper level pressure data indicates non-hydrostatic
              IF(DREC(KG,KT)%VARB_ID(I,J).EQ.'PRES')DREC(KG,KT)%ZFLG=.TRUE.

              !       set moisture as specific humidity ... false means relative humidity
              IF(DREC(KG,KT)%VARB_ID(I,J).EQ.'SPHU')DREC(KG,KT)%QFLG=.TRUE.

              !       set moisture availability by level
              IF(DREC(KG,KT)%VARB_ID(I,J).EQ.'SPHU'.OR.    &
                   DREC(KG,KT)%VARB_ID(I,J).EQ.'RELH')       &
                   DREC(KG,KT)%RFLG(J-1)=.TRUE.

              !       set vertical motion flag if field found
              IF(DREC(KG,KT)%VARB_ID(I,J).EQ.'WWND'.OR.    &
                   DREC(KG,KT)%VARB_ID(I,J).EQ.'DZDT') DREC(KG,KT)%WFLG(J-1)=.TRUE.

              !       set turbulent kinetic energy if field found
              IF(DREC(KG,KT)%VARB_ID(I,J).EQ.'TKEN')DREC(KG,KT)%TKEN=.TRUE.

              !       set velocity variance flag [ note: TKE = 0.5 (UVAR+VVAR+WVAR) ]
              !       assume if the u-component is present, then all three available
              IF(DREC(KG,KT)%VARB_ID(I,J).EQ.'UVAR')DREC(KG,KT)%VELV=.TRUE.
           END DO
        END DO

        ! check selected vertical motion method versus available field
        ! no selection and no data default to isobaric   (kvel=1)
        ! no selection and no data default to divergence (kvel=5)
        !dwen(20090316):starting from HYSPLIT4.9, KVEL=5 instead of 1
        IF((.NOT.DREC(KG,KT)%WFLG(2)).AND.KVEL.EQ.0)KVEL=5

        ! only apply eta correction to ETA coordinate input data
        IF(KVEL.EQ.6)THEN 
           IF(.NOT.(GRID(KG,KT)%MODEL_ID.EQ.' ETA'.OR.   &
                GRID(KG,KT)%MODEL_ID.EQ.'ETA '.OR.   & 
                GRID(KG,KT)%MODEL_ID.EQ.'EDAS'))THEN 
              WRITE(*,*)'*ERROR*: metini - WVEL option 6 only valid with ETA data'
              WRITE(*,*)'Input data file: ',GRID(KG,KT)%MODEL_ID 
              STOP 900
           END IF
        END IF

        ! surface pressure or terrain required for calculations
        IF(.NOT.(DREC(KG,KT)%SHGT.OR.DREC(KG,KT)%PRSS))THEN
           WRITE(*,*)'*ERROR*: metini - input meteorology requires either'
           WRITE(*,*)'   Surface pressure: ',DREC(KG,KT)%PRSS
           WRITE(*,*)'   Terrain height  : ',DREC(KG,KT)%SHGT
           STOP 900
        END IF

        !-------------------------------------------------------------------------------
        ! initialize grid conversion subroutines
        !-------------------------------------------------------------------------------

        !!***************************************************
        !!dwen(20090316):copy the following codes from metini of STILT,add KT
        WRITE (*,*) "MODEL_ID: ",GRID(KG,kt)%MODEL_ID
        WRITE (*,*) "TANG_LAT: ",GRID(KG,kt)%TANG_LAT
        WRITE (*,*) "REF_LAT:  ",GRID(KG,kt)%REF_LAT
        WRITE (*,*) "REF_LON:  ",GRID(KG,kt)%REF_LON

        if (GRID(KG,kt)%pole_lat .le. vmissle .and. GRID(KG,kt)%pole_lon .le. vmissle &
             .and. GRID(KG,kt)%ref_lat .le. vmissle) then
           GRID(KG,kt)%LATLON = .FALSE.
           GRID(KG,kt)%GLOBAL = .FALSE.
           ! using WPS geolocation routines
           dxm = GRID(KG,kt)%size*1000.
           true1 = GRID(KG,kt)%tang_lat
           true2 = GRID(KG,kt)%dummy
           iproj = nint(GRID(KG,kt)%orient)
           synclon180 = GRID(KG,kt)%sync_lon
           if (synclon180 > 180) synclon180 = synclon180-360
           reflon180 = GRID(KG,kt)%ref_lon
           if (reflon180 > 180) reflon180 = reflon180-360
           sync_lat = GRID(KG,kt)%sync_lat
           sync_xp = GRID(KG,kt)%sync_xp
           sync_yp = GRID(KG,kt)%sync_yp
           call map_set(proj_code = iproj,proj = GRID(KG,kt)%proj,&
                lat1 = sync_lat,lon1 = synclon180,knowni=sync_xp,knownj=sync_yp, &
                dx = dxm, stdlon = reflon180, truelat1=true1, truelat2=true2)
           GRID(KG,kt)%gbase(:) = vmiss
           write(*,*) "metini: Using WPS geolocation routines with proj_code: ",iproj
           ! Initialize GRID(KG,kt)%mapfactor
           allocate(xlat_array(grid(kg,kt)%nx,grid(kg,kt)%ny))
           allocate(grid(kg,kt)%mapfactor(grid(kg,kt)%nx,grid(kg,kt)%ny))
           do j_xlat=1,grid(kg,kt)%ny
              xj = j_xlat
              do i_xlat=1,grid(kg,kt)%nx
                 xi = i_xlat
                 call ij_to_latlon(grid(kg,kt)%proj, xi, xj, xlat, xlon)
                 xlat_array(i_xlat,j_xlat) = xlat
              end do
           end do
           call get_map_factor(grid(kg,kt)%proj,xlat_array,grid(kg,kt)%mapfactor, &
                1, 1, grid(kg,kt)%nx, grid(kg,kt)%ny)
           deallocate(xlat_array)

        ELSE
           ! test for global latlon grid system
           CALL GBLSET(KG,KT)

           IF(GRID(KG,KT)%MODEL_ID.EQ.'RAMS') THEN
              WRITE (*,*) "yes RAMS"
!             initialize grid conversion variable (into gbase)
              CALL SOBSTR (GRID(KG,KT)%GBASE,dble(GRID(KG,KT)%REF_LAT),dble(GRID(KG,KT)%REF_LON))
           ELSE
              WRITE (*,*) "not RAMS"
! JCL:(07/06/2004) not on global latlon grid system (flag set in GBLSET)
              IF(.NOT.GRID(KG,KT)%LATLON)THEN
!             initialize grid conversion variable (into gbase)
                 CALL STLMBR_STILT(GRID(KG,KT)%GBASE,dble(GRID(KG,KT)%TANG_LAT),dble(GRID(KG,KT)%REF_LON))
              ENDIF
           ENDIF
           IF (.NOT.GRID(KG,kt)%LATLON) THEN

              !    use single point grid definition
              CALL STCM1P_STILT(GRID(KG,KT))
           END IF
        ENDIF
        WRITE (*,*) "SYNC_XP,SYNC_YP:",GRID(KG,KT)%SYNC_XP,GRID(KG,KT)%SYNC_YP
        WRITE (*,*) "SYNC_LON,SYNC_LAT:",GRID(KG,KT)%SYNC_LON,GRID(KG,KT)%SYNC_LAT
        WRITE (*,*) "SIZE,NX,NY:",GRID(KG,KT)%SIZE,GRID(KG,KT)%NX,GRID(KG,KT)%NY
        !*********************************************
        !
        !-------------------------------------------------------------------------------
        ! check for consistency of meteorological data with turbulence namelist options
        ! then set internal values for each input meteorological data set
        !-------------------------------------------------------------------------------

        DREC(KG,KT)%KDEF=MAX(0,MIN(1,KDEF))
        DREC(KG,KT)%KBLS=MAX(1,MIN(2,KBLS))
        DREC(KG,KT)%KBLT=MAX(1,MIN(4,KBLT))
        DREC(KG,KT)%KZMIX=MAX(0,MIN(1,KZMIX))
        DREC(KG,KT)%TVMIX=MAX(0.01,MIN(100.0,TVMIX))

        IF(KBLS.EQ.1)THEN
           !    needs flux field
           FIELDS=(DREC(KG,KT)%EFLX.OR.DREC(KG,KT)%USTR.OR.DREC(KG,KT)%UFLX).AND.  &
                (DREC(KG,KT)%TSTR.OR.DREC(KG,KT)%HFLX)
           IF(.NOT.FIELDS)THEN
              WRITE(KF21,*)'WARNING: metini - FLUXES not found in data (kg,kt)',KG,KT
              WRITE(KF21,*)'Setting KBLS=2 to use the wind/temp profiles!'
              DREC(KG,KT)%KBLS=2
           END IF
        ELSEIF(KBLS.EQ.2)THEN
           !    wind and temperature profile
           FIELDS=(DREC(KG,KT)%UFLG.AND.DREC(KG,KT)%TFLG)
           IF(.NOT.FIELDS)THEN
              WRITE(KF21,*)'WARNING: metini - SFC wind/temp not found in data (kg,kt)',KG,KT
              WRITE(KF21,*)'KBLS setting unchanged but compuatation may not be as accurate!'
           END IF
        ELSE
           CONTINUE
        END IF

        IF(KBLT.EQ.3)THEN
           !    tke input field
           IF(.NOT.DREC(KG,KT)%TKEN)THEN
              WRITE(KF21,*)'WARNING: metini - TKE not found in data (kg,kt)',KG,KT
              WRITE(KF21,*)'Setting KBLT=1 to use Beljaaars-Betchov parameterization!'
              DREC(KG,KT)%KBLT=1
           END IF
        ELSEIF(KBLT.EQ.4)THEN
           !    velocity variance input field
           IF(.NOT.DREC(KG,KT)%VELV)THEN
              WRITE(KF21,*)'WARNING: metini - VARIANCE not found in data (kg,kt)',KG,KT
              WRITE(KF21,*)'Setting KBLT=1 to use Beljaaars-Betchov parameterization!'
              DREC(KG,KT)%KBLT=1
           END IF
        ELSE
           CONTINUE
        END IF

     END DO
  END DO

  !-------------------------------------------------------------------------------
  ! After all meteorology files have been loaded into multiple grids for time
  ! period one, check for duplicate grids, which, by definition, will be assigned
  ! to different time periods. 
  !-------------------------------------------------------------------------------

  IF(NGRD.GT.1.AND.NTIM.EQ.1)THEN

     !    sort file information by grid size

     DO KG=1,NGRD
        !       initialize dummy array for sorting
        GRID(KG,1)%DUMMY=GRID(KG,1)%SIZE
     END DO

     DO KG=1,NGRD
        DO K2=(KG+1),NGRD
           IF(GRID(K2,1)%DUMMY.LT.GRID(KG,1)%DUMMY)THEN
              !             temporary holding array
              GRID(0,1)=GRID(K2,1)
              DREC(0,1)=DREC(K2,1)
              FILE(0,1)=FILE(K2,1)

              !             exchange with larger grid
              GRID(K2,1)=GRID(KG,1)
              DREC(K2,1)=DREC(KG,1)
              FILE(K2,1)=FILE(KG,1)
              GRID(KG,1)=GRID(0,1)
              DREC(KG,1)=DREC(0,1)
              FILE(KG,1)=FILE(0,1)
           END IF
        END DO
        GRID(KG,1)%DUMMY=0.0
     END DO

     !    combine multiple time periods for the same grid

     DO KG=1,NGRD
        DO K2=(KG+1),NGRD

           !       valid data file defined at both locations
           IF(GRID(KG,1)%NUMBER.GE.0.AND.GRID(K2,1)%NUMBER.GE.0)THEN

              !          grid definition and dimensions the same (sync test 8/7/03)
              IF(ALL(GRID(K2,1)%GBASE.EQ.GRID(KG,1)%GBASE).AND.  &
                   GRID(K2,1)%NX.EQ.GRID(KG,1)%NX.AND.             &
                   GRID(K2,1)%NY.EQ.GRID(KG,1)%NY.AND.             &
                   GRID(K2,1)%SYNC_LAT.EQ.GRID(KG,1)%SYNC_LAT.AND. &
                   GRID(K2,1)%SYNC_LON.EQ.GRID(KG,1)%SYNC_LON)THEN

                 KT=1
                 !             find the next available time index
                 DO WHILE (KT.LT.MTIM.AND.GRID(KG,KT)%NUMBER.GE.0)
                    KT=KT+1
                    NTIM=MAX(NTIM,KT)
                 END DO

                 !             transfer information to new time location
                 GRID(KG,KT)=GRID(K2,1)
                 DREC(KG,KT)=DREC(K2,1)
                 FILE(KG,KT)=FILE(K2,1)

                 !             disable the old location
                 GRID(K2,1)%NUMBER=-99
              END IF

           END IF

        END DO
     END DO

     !    remove deleted grid definitions

     KG=1
     K2=0   
     DO WHILE (KG.LE.NGRD)

        IF(GRID(KG,1)%NUMBER.LT.0)THEN
           !          delete grid entry
           KG=KG+1
        ELSE
           !          valid entry
           K2=K2+1

           IF(KG.NE.K2)THEN
              !             transfer information to new grid location
              GRID(K2,:)=GRID(KG,:)
              DREC(K2,:)=DREC(KG,:)
              FILE(K2,:)=FILE(KG,:)

              !             disable old location
              GRID(KG,:)%NUMBER=-99
           END IF
           KG=KG+1
        END IF
     END DO
     NGRD=K2

  END IF

  !-------------------------------------------------------------------------------
  ! If the same meteorological data grid is defined for multiple time periods,   
  ! perhaps because each represents the next forcast cycle, then sort the files
  ! from the oldest to the most recent forecast - added 22 March 2007 (RRD)
  !-------------------------------------------------------------------------------

  IF(NTIM.GT.1)THEN

     DO KG=1,NGRD
        DO KT=1,NTIM

           IF(GRID(KG,KT)%NUMBER.GE.0.AND.GRID(K2,KT)%NUMBER.GE.0)THEN

              DO K2=(KT+1),NTIM
                 IF(FILE(KG,K2)%FIRST%MACC.LT.FILE(KG,KT)%FIRST%MACC)THEN
                    !                temporary holding array
                    GRID(KG,0)=GRID(KG,K2)
                    DREC(KG,0)=DREC(KG,K2)
                    FILE(KG,0)=FILE(KG,K2)

                    !                exchange with more recent data
                    GRID(KG,K2)=GRID(KG,KT)
                    DREC(KG,K2)=DREC(KG,KT)
                    FILE(KG,K2)=FILE(KG,KT)
                    GRID(KG,KT)=GRID(KG,0)
                    DREC(KG,KT)=DREC(KG,0)
                    FILE(KG,KT)=FILE(KG,0)
                 END IF
              END DO

           END IF

        END DO
     END DO

  END IF

  !-------------------------------------------------------------------------------
  ! Set proper starting time if mo=0. The initial model start time is determined
  ! from the initial or last meteorological file time.  In this situation,
  ! day and hour are assumed to be relative to the file time.
  !-------------------------------------------------------------------------------

  IF(IBMO.EQ.0)THEN

     CALL TM2MIN(39,12,31,23,0,MACC) ! max date is 2039
     IF(BACK) MACC=0 

     DO KG=1,NGRD
        DO KT=1,NTIM

           IF(GRID(KG,KT)%NUMBER.GE.0)THEN
              IF(BACK)THEN
                 IBYR=FILE(KG,KT)%LAST%YR
                 IBMO=FILE(KG,KT)%LAST%MO
                 ITDA=FILE(KG,KT)%LAST%DA-IBDA
                 ITHR=FILE(KG,KT)%LAST%HR-IBHR
                 ITMN=FILE(KG,KT)%LAST%MN-IBMN
                 CALL TM2MIN(IBYR,IBMO,ITDA,ITHR,ITMN,MJET)
                 MACC=MAX(MACC,MJET)
              ELSE
                 IBYR=FILE(KG,KT)%FIRST%YR
                 IBMO=FILE(KG,KT)%FIRST%MO
                 ITDA=FILE(KG,KT)%FIRST%DA+IBDA
                 ITHR=FILE(KG,KT)%FIRST%HR+IBHR
                 ITMN=FILE(KG,KT)%FIRST%MN+IBMN
                 CALL TM2MIN(IBYR,IBMO,ITDA,ITHR,ITMN,MJET)
                 MACC=MIN(MACC,MJET)
              END IF
           END IF

        END DO
     END DO

     !    adjust relative date for potential month crossing error
     CALL TM2DAY(MACC,IBYR,IBMO,IBDA,IBHR,IBMN)
  END IF

END SUBROUTINE metini
