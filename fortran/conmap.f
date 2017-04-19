!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONMAP           CONCENTRATION MAPPING ROUTINE FOR PLOTS
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            IS THE MAIN SUBROUTINE THAT CONTROLS THE PLOTTING OF AIR
!            CONCENTRATIONS.  IT DRAWS THE MAP BACKGROUND,
!            CALLS ROUTINES GRID THE DATA, DRAW CONTOURS,
!            PROVIDE FOR COLOR FILL, AND LABEL THE COMPONENTS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 06 Feb 1998 (RRD)
!                 18 Aug 1998 (RRD) - maphi option to contouring
!                 09 Dec 1998 (RRD) - convert NCAR grapics to psplot
!                 28 May 1999 (RRD) - argument list change mapbdy
!                 28 Mar 2000 (RRD) - additional temp variable
!                 20 Nov 2000 (RRD) - fortran90 upgrade
!                 12 Dec 2000 (RRD) - optional smoothing
!                 20 Apr 2001 (RRD) - mapbdy argument list
!                 02 Nov 2001 (RRD) - additional label information
!                 16 Nov 2001 (RRD) - revised fixed contour method
!                 19 Dec 2001 (RRD) - source plotting symbol variable
!                 11 Jun 2002 (RRD) - station file plotting option (-q)
!                 21 Jun 2002 (RRD) - concentration minimum for plot
!                 09 Sep 2002 (RRD) - test circle for label within plot
!                 07 Jan 2003 (RRD) - multiple source labels
!                 01 Jul 2003 (RRD) - no station file plot values <0
!                 20 Nov 2003 (RRD) - added minutes to emission time
!                 04 Dec 2003 (RRD) - converted llint to aintr
!                 26 Oct 2004 (RJP) - added dbf capability
!                 20 Dec 2004 (RRD) - g95 compatibility
!                 07 Mar 2006 (GDR) - added source lat/lon to contur
!                 16 Mar 2006 (GDR) - added current map number to contur
!                 10 Jan 2007 (RRD) - more lat-lon labelling options
!                 12 Jan 2007 (RRD) - added cylindrical equidistant
!                 26 Jan 2007 (RRD) - minutes field pass to contxt & contur
!                 09 Oct 2007 (GDR) - passed IYR through to contur
!                 22 Oct 2007 (BS)  - -k3 option
!                 30 Nov 2007 (GDR) - passed IYR through to contxt
!                 12 Aug 2008 (RRD) - shapefile mapping option
!
! USAGE: CALL CONMAP(SCAN,KONTUR,KMAP,KAVG,CONC,LEVEL1,LEVEL2,NVAL,TFACT,MODEL,
!        PTYPE,KOLOR,PMODE,QUNIT,NLOC,IBYR,IBMO,IBDA,IBHR,IBMN,OLAT,OLON,
!        OLVL1,OLVL2,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH,NLAT,NLON,DLAT,DLON,
!        CLAT,CLON,GFILE,PARMAP,XCON,NXP,NYP,ALONL,ALONR,ALATT,ALATB,CMAX,
!        CMIN,CVAL,NRAD,DIST,SYM,MAPN,QCHAR,QFILE,PROCESS,GOPEN,MAPDEL,KPROJ,CLBL)
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

SUBROUTINE CONMAP(SCAN,KONTUR,KMAP,KAVG,CONC,LEVEL1,LEVEL2,NVAL,TFACT,MODEL,   &
           PTYPE,KOLOR,PMODE,QUNIT,NLOC,IBYR,IBMO,IBDA,IBHR,IBMN,OLAT,OLON,    &
           OLVL1,OLVL2,IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH,NLAT,NLON,  &
           DLAT,DLON,CLAT,CLON,GFILE,PARMAP,XCON,NXP,NYP,ALONL,ALONR,ALATT,    &
           ALATB,CMAX,CMIN,CVAL,NRAD,DIST,SYM,MAPN,QCHAR,QFILE,PROCESS,GOPEN,  &
           MAPDEL,KPROJ,CLBL)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)    :: scan          ! scan distance for smoothing    
  INTEGER, INTENT(IN)    :: kontur        ! contour flag for auto or fixed
  INTEGER, INTENT(IN)    :: kmap,kavg     ! flag to indicate concen, exposure, deposit or AEGLs
  REAL,    INTENT(IN)    :: conc (:,:)    ! input concentration grid 
  INTEGER, INTENT(IN)    :: level1,level2 ! range of display levels or layer
  INTEGER, INTENT(IN)    :: nval          ! number of contour values
  REAL,    INTENT(IN)    :: tfact         ! units multiplication factor
  INTEGER, INTENT(IN)    :: kolor         ! B&W=0 or color=1
  INTEGER, INTENT(IN)    :: nloc          ! actual number of source locations
  INTEGER, INTENT(IN)    :: ibyr (:)      ! emission start year 
  INTEGER, INTENT(IN)    :: ibmo (:)      ! emission start month 
  INTEGER, INTENT(IN)    :: ibda (:)      ! emission start day
  INTEGER, INTENT(IN)    :: ibhr (:)      ! emission start hour
  INTEGER, INTENT(IN)    :: ibmn (:)      ! emission start minute 
  REAL,    INTENT(IN)    :: olat (:)      ! starting location 
  REAL,    INTENT(IN)    :: olon (:)      ! starting location 
  REAL,    INTENT(IN)    :: olvl1,olvl2   ! release height range
  INTEGER, INTENT(IN)    :: imo,ida,ihr   ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: imn,iyr       ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: jmo,jda,jhr   ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jmn,jyr       ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jfh           ! stop time of the sample
  INTEGER, INTENT(IN)    :: nlat,nlon     ! number of points in concen grid
  REAL,    INTENT(IN)    :: dlat,dlon     ! resolution of concen grid
  REAL,    INTENT(IN)    :: clat,clon     ! center of the concentration grid
  REAL,    INTENT(IN)    :: parmap (9)    ! describes the conformal map
  REAL,    INTENT(OUT)   :: xcon (:,:)    ! remapped concentration array
  INTEGER, INTENT(IN)    :: nxp,nyp       ! remapped array dim limits
  REAL,    INTENT(IN)    :: alonl,alatb   ! limits of the display map
  REAL,    INTENT(IN)    :: alonr,alatt   ! limits of the display map
  REAL,    INTENT(INOUT) :: cmax,cmin     ! maximum and minimum concen values
  REAL,    INTENT(INOUT) :: cval(:)       ! concentration contour values

  CHARACTER(4),  INTENT(IN) :: model      ! ident of meteorological data
  CHARACTER(4),  INTENT(IN) :: ptype      ! ident of pollutant displayed
  CHARACTER(80), INTENT(IN) :: qunit      ! pollutant mass units
  CHARACTER(80), INTENT(IN) :: gfile      ! graphics map background file
  LOGICAL,       INTENT(IN) :: pmode      ! true indicates portrait display
  INTEGER,       INTENT(IN) :: nrad       ! number of concentric circles
  REAL,          INTENT(IN) :: dist       ! distance between circles (km)
  REAL,          INTENT(IN) :: sym        ! plotting symbol threshold    
  INTEGER,       INTENT(IN) :: mapn       ! current plot map number        
  CHARACTER(1),  INTENT(IN) :: qchar      ! source plot symbol
  CHARACTER(80), INTENT(IN) :: qfile      ! data plot file name          
  CHARACTER(8),  INTENT(IN) :: process    ! process ID suffix
  CHARACTER(8),  INTENT(IN) :: clbl(:)    ! concentration contour labels

  LOGICAL,       INTENT(INOUT) :: gopen   ! Google Earth output file open
  INTEGER,       INTENT(IN)    :: mapdel  ! advanced labelling options
  INTEGER,       INTENT(IN)    :: kproj   ! projection flag           

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  LOGICAL      :: ftest, ltest, ptest1, ptest2
  REAL         :: XU1,XU2,YU1,YU2
  REAL         :: qid,qlat,qlon,qvalue,delt
  REAL         :: tmax,tmin,alon,xmet,ymet,xdis,ydis,xp,yp,cgszll,cylzll
  INTEGER      :: n,kl,ngpt,kerr
  CHARACTER(8) :: label

! potential lat/lon label intervals
  INTEGER, PARAMETER :: NLAB = 10
  REAL               :: AINTR(NLAB)
  DATA AINTR/45.0,30.0,20.0,15.0,10.0,5.0,2.0,1.0,0.5,0.2/

! map boundary testing
  COMMON /MAPBND/ XU1,XU2,YU1,YU2
  PTEST1(XP,YP)=(XP.GE.XU1.AND.XP.LE.XU2.AND.YP.GE.YU1.AND.YP.LE.YU2)
  PTEST2(XP,YP)=(XP.GE.XU2.AND.XP.LE.XU1.AND.YP.GE.YU2.AND.YP.LE.YU1)

!-------------------------------------------------------------------------------
  INTERFACE
!-------------------------------------------------------------------------------
  SUBROUTINE CON2XY(KPROJ,NXP,NYP,XCON,TFACT,SCAN,CONC,PARMAP,NLAT,NLON,       &
                    DLAT,DLON,CLAT,CLON)
  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: kproj       ! projection
  INTEGER, INTENT(IN)   :: nxp,nyp     ! output grid dimensions
  REAL,    INTENT(IN)   :: tfact       ! concentration units multiplier
  INTEGER, INTENT(IN)   :: scan        ! map smoothing distance scan        
  REAL,    INTENT(IN)   :: conc (:,:)  ! concentration grid on lat/lon system
  REAL,    INTENT(IN)   :: parmap (9)  ! concentration grid on lat/lon system
  INTEGER, INTENT(IN)   :: nlat,nlon   ! size of input concentration grid
  REAL,    INTENT(IN)   :: dlat,dlon   ! grid resolution in lat/lon 
  REAL,    INTENT(IN)   :: clat,clon   ! center of input concentration grid
  REAL,    INTENT(OUT)  :: xcon (:,:)  ! concentration array on conformal grid 
  END SUBROUTINE con2xy
!-------------------------------------------------------------------------------
  SUBROUTINE CONTUR(KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,KONTUR,NXP,NYP,XCON,CMIN,   &
                    CMAX,TMAX,TMIN,CVAL,NVAL,KOLOR,NGPT,SYM,IYR,IMO,IDA,IHR,IMN, &
                    JYR,JMO,JDA,JHR,JMN,JFH,OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,    &
                    KAVG,MAPN,PROCESS)
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: kproj         ! projection                   
  INTEGER, INTENT(IN)    :: level1,level2 ! range of display levels
  CHARACTER(4),INTENT(IN):: ptype         ! ident of pollutant displayed
  REAL,    INTENT(IN)    :: parmap (9)    ! conformal map projection parameters
  INTEGER, INTENT(IN)    :: kontur        ! flag auto or fixed contours
  INTEGER, INTENT(IN)    :: nxp,nyp       ! size of the concentration array
  REAL,    INTENT(IN)    :: xcon (:,:)    ! concentration array
  INTEGER, INTENT(IN)    :: nval          ! number of contours
  INTEGER, INTENT(IN)    :: kolor         ! flag for B&W=0 or color=1
  REAL,    INTENT(IN)    :: cmin,cmax     ! minimum and maximum concentrations
  REAL,    INTENT(OUT)   :: tmax,tmin     ! local time max for display 
  REAL,    INTENT(INOUT) :: cval(:)       ! contour array
  INTEGER, INTENT(IN)    :: ngpt          ! number of grid points per input
  REAL,    INTENT(IN)    :: sym           ! plotting symbol threshold (>=0)
  INTEGER, INTENT(IN)    :: imo,ida,ihr   ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: imn,iyr       ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: jmo,jda,jhr   ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jmn,jyr       ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jfh           ! stop time of the sample
  INTEGER, INTENT(IN)    :: mapn          ! current plotting map number
  REAL,    INTENT(IN)    :: olat (:)      ! starting location 
  REAL,    INTENT(IN)    :: olon (:)      ! starting location
  REAL,    INTENT(IN)    :: olvl1,olvl2   ! release height range
  INTEGER, INTENT(IN)    :: kmap,kavg     ! flag to indicate concen, exposure, deposit or AEGLs
  LOGICAL, INTENT(INOUT) :: gopen         ! Google Earth output file open
  CHARACTER(8),INTENT(IN):: process       ! process ID suffix
  END SUBROUTINE contur
!-------------------------------------------------------------------------------
  SUBROUTINE CONTXT(LEVEL1,LEVEL2,KMAP,KAVG,MODEL,PTYPE,NVAL,CMAX,CMIN,CVAL,   &
                    KOLOR,PMODE,QUNIT,IBYR,IBMO,IBDA,IBHR,IBMN,OLAT,OLON,OLVL1,&
                    OLVL2,IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH,MAPN,QCHAR,  &
                    PROCESS,CLBL)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: level1,level2     ! range of display levels 
  INTEGER, INTENT(IN)  :: kmap,kavg         ! map type (concen/exposure/deposit/AEGLs)
  INTEGER, INTENT(IN)  :: nval              ! number of contour values
  REAL,    INTENT(IN)  :: cmax,cmin         ! max and min concentrations  
  REAL,    INTENT(IN)  :: cval(:)           ! concentration contour values
  INTEGER, INTENT(IN)  :: kolor             ! flag to set B&W=0 or color=1
  LOGICAL, INTENT(IN)  :: pmode             ! true indicates portrait display
  INTEGER, INTENT(IN)  :: ibyr (:)          ! emission start year 
  INTEGER, INTENT(IN)  :: ibmo (:)          ! emission start month
  INTEGER, INTENT(IN)  :: ibda (:)          ! emission start day
  INTEGER, INTENT(IN)  :: ibhr (:)          ! emission start hour
  INTEGER, INTENT(IN)  :: ibmn (:)          ! emission start minute
  REAL,    INTENT(IN)  :: olat (:)          ! starting location 
  REAL,    INTENT(IN)  :: olon (:)          ! starting location
  REAL,    INTENT(IN)  :: olvl1,olvl2       ! release height range
  INTEGER, INTENT(IN)  :: imn,jmn           ! minutes field    
  INTEGER, INTENT(IN)  :: imo,ida,ihr       ! starting date/time of the sample
  INTEGER, INTENT(IN)  :: iyr               ! starting date/time of the sample
  INTEGER, INTENT(IN)  :: jmo,jda,jhr       ! stop date/time of the sample
  INTEGER, INTENT(IN)  :: jyr,jfh           ! stop date/time of the sample
  INTEGER, INTENT(IN)  :: mapn              ! current plotting map number
  CHARACTER(4),  INTENT(IN) :: model        ! ident of meteorological data
  CHARACTER(4),  INTENT(IN) :: ptype        ! ident of pollutant displayed
  CHARACTER(80), INTENT(IN) :: qunit        ! pollutant mass units
  CHARACTER(1),  INTENT(IN) :: qchar        ! source plotting symbol
  CHARACTER(8),  INTENT(IN) :: process      ! process ID suffix
  CHARACTER(8),  INTENT(IN) :: clbl(:)      ! concentration contour labels
  END SUBROUTINE contxt  
!-------------------------------------------------------------------------------
  END INTERFACE

!-------------------------------------------------------------------------------
! contour data array
!-------------------------------------------------------------------------------

! determine the number of interpolated grid points per input grid point
  IF(KPROJ.EQ.4)THEN
     NGPT=(DLAT+DLON)*55.5/CYLZLL(CLAT+DLAT*NLAT*0.5,CLON+DLON*NLON*0.5)
  ELSE
     NGPT=(DLAT+DLON)*55.5/CGSZLL(PARMAP,CLAT+DLAT*NLAT*0.5,CLON+DLON*NLON*0.5)
  END IF

! map lat/lon concentration array to conformal grid
  CALL CON2XY(KPROJ,NXP,NYP,XCON,TFACT,SCAN,CONC,PARMAP,NLAT,NLON,             &
              DLAT,DLON,CLAT,CLON)

! put symbols on map corresponding to concentrations
  CALL CONTUR(KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,KONTUR,NXP,NYP,XCON,CMIN,CMAX,  &
           TMAX,TMIN,CVAL,NVAL,KOLOR,NGPT,SYM,IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA, &
           JHR,JMN,JFH,OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,KAVG,MAPN,PROCESS)

!-------------------------------------------------------------------------------
! draw map background
!-------------------------------------------------------------------------------

  ALON=ALONR
! check for dateline/prime meridian switch
  IF(SIGN(1.0,ALONL).NE.SIGN(1.0,ALONR).AND.ALONR.LT.0.0) ALON=360.0+ALONR

! set interval to have at least 4 lat/lon lines on a map

  DELT=MAX(ABS(ALATT-ALATB)/4.0,ABS(ALONL-ALON)/4.0,AINTR(NLAB))
  KL=1
  DO WHILE (KL.LT.NLAB.AND.AINTR(KL).GT.DELT)
     KL=KL+1
  END DO

! enhanced labelling options (10 Jan 2007)
  IF(MAPDEL.EQ.0)THEN
!    labels turned off
     DELT=0.0
  ELSEIF(MAPDEL.LT.0)THEN
!    interval set on command line integer tenths
     DELT=ABS(FLOAT(MAPDEL)/10.0)
  ELSE
!    default auto label option mapdel=1
     DELT=AINTR(KL)
  END IF

  IF(KOLOR.EQ.1.OR.KOLOR.EQ.2)CALL SETCOLR(0.2,0.6,0.8)
! map background
  KL=LEN_TRIM(GFILE)
  IF(KL.GT.13)THEN
     IF(GFILE(kl-13:kl).EQ.'shapefiles.txt')THEN
        CALL SHPBDY(KPROJ,.FALSE.,PARMAP,DELT,1.0,GFILE)
     ELSE
        CALL MAPBDY(KPROJ,.FALSE.,PARMAP,DELT,1.0,GFILE)
     END IF
  ELSE
     CALL MAPBDY(KPROJ,.FALSE.,PARMAP,DELT,1.0,GFILE)
  END IF

!-------------------------------------------------------------------------------
! label plot
!-------------------------------------------------------------------------------

! plot source location position
  DO N=1,NLOC
     IF(KPROJ.EQ.4)THEN
        CALL CYL2XY(OLAT(N),OLON(N),XMET,YMET)
     ELSE
        CALL CLL2XY(PARMAP,OLAT(N),OLON(N),XMET,YMET)
     END IF
     IF(PTEST1(XMET,YMET).OR.PTEST2(XMET,YMET))THEN
        CALL MAP2XY(XMET,YMET,XP,YP)
        CALL SETFNT(35)
        CALL KEKSYMC(XP,(YP-0.06),0.12,QCHAR,0.0,1,1)
        CALL SETFNT(20)
     END IF
  END DO

! optional concentric circles
  IF(NRAD.GT.0)THEN
     IF(KPROJ.EQ.4)THEN
        CALL CYL2XY(OLAT(1),OLON(1),XMET,YMET)
     ELSE
        CALL CLL2XY(PARMAP,OLAT(1),OLON(1),XMET,YMET)
     END IF
     CALL MAP2XY(XMET,YMET,XP,YP)

     IF(KPROJ.EQ.4)THEN
        CALL CYL2XY(OLAT(1)+DIST/111.0,OLON(1),XMET,YMET)
     ELSE
        CALL CLL2XY(PARMAP,OLAT(1)+DIST/111.0,OLON(1),XMET,YMET)
     END IF
     CALL MAP2XY(XMET,YMET,XDIS,YDIS)
     DO N=1,NRAD
!       CALL CIRCLE(XP,YP,(ABS(YDIS-YP)*N),.FALSE.)
        CALL XCIRCLE(XP,YP,(ABS(YDIS-YP)*N),LTEST)
        IF(LTEST)THEN
           WRITE(LABEL,'(I4,A)')NINT(DIST*N),' km '
           CALL KEKSYMC(XP,(YP-ABS(YDIS-YP)*N),0.12,LABEL,0.0,8,1)
        END IF
     END DO
  END IF

! station file data plotting
  INQUIRE(FILE=QFILE,EXIST=FTEST)
  IF(FTEST)THEN 
     OPEN(20,FILE=QFILE,STATUS='OLD')
     KERR=0
     DO WHILE (KERR.EQ.0)
        READ(20,*,IOSTAT=KERR)QID,QLAT,QLON,QVALUE
        IF(KPROJ.EQ.4)THEN
           CALL CYL2XY(QLAT,QLON,XMET,YMET)
        ELSE
           CALL CLL2XY(PARMAP,QLAT,QLON,XMET,YMET)
        END IF
        CALL MAP2XY(XMET,YMET,XP,YP)
        LABEL='        ' 
!       x,y,height,symb,ang,nchar,{0=ll,1=cntr,2=lr)
        CALL KEKSYMC(XP,YP,0.08,'+',0.0,1,1)
        IF(QVALUE.GT.0.0)THEN
           WRITE(LABEL(1:7),'(I7)')NINT(QVALUE) 
           CALL KEKSYMC(XP,YP,0.08,LABEL,0.0,8,2)
        END IF
     END DO
     CLOSE (20)
  END IF

! place identification labels
  CALL CONTXT(LEVEL1,LEVEL2,KMAP,KAVG,MODEL,PTYPE,NVAL,TMAX,TMIN,CVAL,   &
              KOLOR,PMODE,QUNIT,IBYR,IBMO,IBDA,IBHR,IBMN,OLAT,OLON,OLVL1,&
              OLVL2,IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,JFH,MAPN,QCHAR,  &
              PROCESS,CLBL)

END SUBROUTINE conmap

SUBROUTINE XCIRCLE(XP,YP,DIST,LTEST)

   LOGICAL PTEST, LTEST
   REAL XARR(500), YARR(500)
   COMMON /MAPCLP/ XC1,XC2,YC1,YC2

   PTEST(XP,YP)=(XP.GE.XC1.AND.XP.LE.XC2.AND.YP.GE.YC1.AND.YP.LE.YC2)    &
            .OR.(XP.GE.XC2.AND.XP.LE.XC1.AND.YP.GE.YC2.AND.YP.LE.YC1)

   LTEST=.FALSE.  ! flag to write distance label

   NPTS=0
   DO IANG=0,360,2
      ANGINR=2.0*FLOAT(IANG)*8.726646E-03
      XA=XP+DIST*SIN(ANGINR)
      YA=YP+DIST*COS(ANGINR)

      IF(PTEST(XA,YA))THEN
         NPTS=NPTS+1
         XARR(NPTS)=XA
         YARR(NPTS)=YA
      ELSE
         IF(NPTS.GT.1)CALL SLDCRV(XARR,YARR,NPTS,0.01)
         NPTS=0
      END IF

!     if bottom of circle within plot then write label
      IF(IANG.EQ.180.AND.PTEST(XA,YA)) LTEST=.TRUE.

   END DO
   IF(NPTS.GT.1)CALL SLDCRV(XARR,YARR,NPTS,0.01)

END SUBROUTINE xcircle
