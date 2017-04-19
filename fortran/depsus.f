!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPSUS           DEPosition reSUSpension of a pollutant
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION RESUSPENSION OF A POLLUTANT - ASSUME THAT RATIO (K)
!   OF POLLUTANT IN AIR (C=MASS/M3) TO SURFACE VALUE (S=MASS/M2) =
!   10^-6 DIVIDED BY DURATION OF DEPOSITION IN DAYS.  FOR SIMPLICITY
!   WE ASSUME DAYS ALWAYS = 1.  K ALSO DEFINED AS R/S DS/DT, WHERE R
!   IS RELATED TO THE ATMOSPHERIC RESISTENCE = 1/KU*   THEREFORE
!   THE RESUSPENSION FLUX = S K / R.  NOTE THAT MULTIPLE SPECIES
!   CONCENTRATION FILES WILL RESULTS IN INDEPENDENT PARTICLES FOR EACH
!   DEFINED POLLUTANT FROM ONLY THE FIRST CONCENTRATION GRID.  MULTI
!   GRID DEFINITIONS ARE NOT SUPPORTED IN THIS APPLICATION. NOTE THAT
!   THE METEO VARIABLE FRICTION VELOCITY IS PASSED THROUGH AS THE
!   VALUE COMPUTED IN ADVPNT, THEREFORE ONLY VALID IN A SMALL DOMAIN
!   STANDARD MODEL CONFIGURATION DOES NOT SUPPORT FULL-GRID METEO.
! Revised (05 Aug 2003) such that subgrid meteorology is not supported when
! running this subroutine. The meteo subgrid must be set to full grid by
! setting the MGMIN parameter in the namelist to a large value. Required
! to determine U* at each emission grid point.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 11 Jun 1997 (RRD)  
!                 21 Sep 2000 (RRD) - fortran90 upgrade
!                 16 Mar 2001 (RRD) - global lat long grid option
!                 29 Aug 2001 (RRD) - simultaneous multiple meteorology
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 05 Aug 2003 (RRD) - now requires full meteo grid
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 10 Aug 2004 (RRD) - enhanced particle-puff conversion
!                 17 Oct 2007 (RRD) - absolute value on time step
!
! USAGE:  CALL DEPSUS(CONC,DIRT,INITD,KGM,KTM,NUMGRD,NUMTYP,DT,ICHEM,KPM,
!                     CSUM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,
!                     PAGE,PTYP,PGRD,NSORT,MAXPAR)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE DEPSUS(CONC,DIRT,INITD,KGM,KTM,NUMGRD,NUMTYP,DT,ICHEM,KPM,  &
                  CSUM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGWW,HDWP,   &
                  PAGE,PTYP,PGRD,NSORT,MAXPAR)

  USE funits
  USE metval
  use module_defgrid ! meteorological grid and file definitions

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC' ! concentration grid and pollutant definitions

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(cset), INTENT(IN)    :: conc(:)      ! for each concentration grid 
  TYPE(pset), INTENT(IN)    :: dirt(:)      ! for each pollutant type 
  INTEGER,    INTENT(IN)    :: initd        ! initial puff/particle distribution
  INTEGER,    INTENT(IN)    :: kgm          ! current meteorological grid number
  INTEGER,    INTENT(IN)    :: ktm          ! current meteorological time number
  INTEGER,    INTENT(IN)    :: numgrd       ! number of concentration grids
  INTEGER,    INTENT(IN)    :: numtyp       ! number of pollutants
  REAL,       INTENT(IN)    :: dt           ! time step (min)
  INTEGER,    INTENT(IN)    :: ichem        ! chemistry option parameter
  INTEGER,    INTENT(INOUT) :: kpm          ! total number of puffs or particles
  REAL,       INTENT(INOUT) :: csum  (:,:,:,:,:)                
  REAL,       INTENT(INOUT) :: mass  (:,:)  ! (species, particles)
  REAL,       INTENT(INOUT) :: xpos  (:)    ! x position (numb part)
  REAL,       INTENT(INOUT) :: ypos  (:)    ! y position (numb part)
  REAL,       INTENT(INOUT) :: zpos  (:)    ! z position (numb part)
  REAL,       INTENT(INOUT) :: sigh  (:)    ! horiz sigma (meters)
  REAL,       INTENT(INOUT) :: sigu  (:)    ! u'2 turbulence  
  REAL,       INTENT(INOUT) :: sigv  (:)    ! v'2 turbulence  
  REAL,       INTENT(INOUT) :: sigww (:)    ! vertical sigma
  INTEGER,    INTENT(INOUT) :: hdwp  (:)    ! Horizontal distribution
  INTEGER,    INTENT(INOUT) :: page  (:)    ! pollutant age (min)     
  INTEGER,    INTENT(INOUT) :: ptyp  (:)    ! pollutant type index numb 
  INTEGER,    INTENT(INOUT) :: pgrd  (:)    ! particle on this grid numb
  INTEGER,    INTENT(INOUT) :: nsort (:)    ! sorted array index values
  INTEGER,    INTENT(IN)    :: maxpar       ! maximum number of particles

!-------------------------------------------------------------------------------
! local variables definitions
!-------------------------------------------------------------------------------

  REAL,    PARAMETER   :: vonk    = 0.40      ! Von Karman's
  REAL,    PARAMETER   :: sigr    = 1.54      ! Radius definition
  REAL,    PARAMETER   :: PI      = 3.14159265358979
  REAL,    PARAMETER   :: DEGPRD  = 180.0/PI  ! deg per radian
 
  INTEGER              :: ii,jj,kt,kl,kg,kgc,kgl,nxp,nyp
  REAL                 :: xp,yp,depamt,qtot,qrate,area,plat,plon,ustr

!-------------------------------------------------------------------------------
! external variables definitions
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE GBL2XY(KG,KT,CLAT,CLON,X,Y)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: kg          ! active grid number    
  INTEGER, INTENT(IN)  :: kt          ! active time number    
  REAL,    INTENT(IN)  :: clat,clon   ! latlon location       
  REAL,    INTENT(OUT) :: x,y         ! grid position         
  END SUBROUTINE GBL2XY
  END INTERFACE
!-------------------------------------------------------------------------------

  KGC=0
! find the first grid and level which has deposition defined
  DO KG=1,NUMGRD
     IF(KGC.EQ.0)THEN
        DO KL=1,CONC(KG)%LEVELS
           IF(CONC(KG)%HEIGHT(KL).EQ.0)THEN
              KGC=KG
              KGL=KL
           END IF
        END DO
     END IF
  END DO
  IF(KGC.EQ.0)THEN
     WRITE(KF21,*)'*ERROR* depsus: requires a deposition grid'
     WRITE(*,*)   '*ERROR* depsus: see message file for more information'
     STOP 900
  END IF

! determine the concentration grid size
  IF(ICHEM.EQ.4)THEN
     NXP=GRID(KGM,KTM)%NX
     NYP=GRID(KGM,KTM)%NY
  ELSE
     NXP=CONC(KG)%NUMB_LON
     NYP=CONC(KG)%NUMB_LAT
  END IF

  DO KT=1,NUMTYP
!    resuspension must be defined for this pollutant
     IF(DIRT(KT)%DOSUS)THEN

     DO JJ=1,NYP
     DO II=1,NXP

        DEPAMT=CSUM(II,JJ,KGL,KT,KGC)
        IF(DEPAMT.GT.0.0)THEN

!          increment particle/puff counter
           KPM=KPM+1
           IF(KPM.GT.MAXPAR)THEN
              KPM=MAXPAR
              WRITE(KF21,*)'WARNING depsus: exceeding puff limit'
              RETURN
           END IF
           NSORT(KPM)=KPM

           IF(ICHEM.EQ.4)THEN
!             internal concentration grid matches meteorology grid
              XP=FLOAT(II)
              YP=FLOAT(JJ)
           ELSE
!             convert position to meteorological grid units
              PLON=FLOAT(II-1)*CONC(KGC)%DELT_LON+CONC(KGC)%X1Y1_LON
              PLAT=FLOAT(JJ-1)*CONC(KGC)%DELT_LAT+CONC(KGC)%X1Y1_LAT
              IF(GRID(KGM,KTM)%LATLON)THEN
                 CALL GBL2XY(KGM,KTM,PLAT,PLON,XP,YP)
              ELSE
                 CALL CLL2XY_wps(GRID(KGM,KTM)%GBASE,PLAT,PLON, XP, YP,GRID(KGM,KTM)%proj)
              END IF
           END IF

!          initial position always at ground
           XPOS(KPM)=XP
           YPOS(KPM)=YP
           ZPOS(KPM)=1.0

!          alongwind and vertical variances start at zero
           SIGU(KPM)=0.0
           SIGV(KPM)=0.0
           SIGWW(KPM)=0.0

!          initial distribution (see main for definitions)
           HDWP(KPM)=INITD
!          initial age at zero
           PAGE(KPM)=0
!          pollutant type definition
           PTYP(KPM)=KT
!          initial grid is the default startup grid from main
           PGRD(KPM)=KGM

!          assume grid-cell area source (111000 m / deg - squared)
           AREA=1.2E+10*CONC(KGC)%DELT_LAT                              &
                       *CONC(KGC)%DELT_LON*COS(PLAT/DEGPRD)
!          compute sigma for uniform radius
           SIGH(KPM)=SQRT(AREA/PI)/SIGR

!          pollutant flux - mass/m2-s
           USTR=UF(NINT(XP),NINT(YP),POINT(2),KGM)
           QRATE=DEPAMT*VONK*USTR*DIRT(KT)%SRATE

!          determine mass lost from surface
           QTOT=MIN(DEPAMT, 60.0*ABS(DT)*QRATE)
           CSUM(II,JJ,KGL,KT,KGC)=DEPAMT-QTOT
           MASS(KT,KPM)=QTOT*AREA

        END IF

!    horizontal grid loop
     END DO
     END DO

! pollutant types
  END IF
  END DO

END SUBROUTINE depsus
