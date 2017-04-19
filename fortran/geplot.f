!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GEPLOT           MAP CONTOUR DRAW
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-12-11
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            CREATES THE TEXT NEEDED TO CREATE GOOGLE EATH kml OUTPUT
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 23 Feb 2007 (GDR)
!                 14 Sep 2007 (GDR) - agl/msl issues
!                 09 Oct 2007 (GDR) - pass in IYR from contur
!                 09 Oct 2007 (GDR) - added TimeSpan to Google Earth
!                 09 Oct 2007 (GDR) - made more transparent contours
!                 30 Nov 2007 (GDR) - Correction for 4 digit year
!                 03 Jan 2008 (GDR) - Modification to position information 
!                                     for ESRI Explorer
!                 28 Jan 2009 (GDR) - Improved display by draping plume over terrain
!                                   - Added NOAA Logo and improved text on graphic
!
! USAGE: CALL GEPLOT( see below )
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            GISTMP, GISOUT, GISATT
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!-------------------------------------------------------------------------------
! KMAP defines type of graphic to display
!    KMAP=1 is concentration
!    KMAP=2 is exposure
!    KMAP=3 is deposition
!    KMAP=4 is chemical threshold
!    KMAP=5 is volcanic ash (not used here)
!-------------------------------------------------------------------------------
! KGE defines type of kml to display
!    KGE=1 : write kml header information (first time only)
!    KGE=2 : zero concentrations still needs a frame
!    KGE=3 : normal contour output of data
!-------------------------------------------------------------------------------
!$$$

SUBROUTINE GEPLOT(KGE,KPROJ,LEVEL1,LEVEL2,PTYPE,PARMAP,NXP,NYP,   &
                  TMAX,TMIN,XMAX,YMAX,CVAL,NVAL,NGPT,             &
                  IYR,IMO,IDA,IHR,IMN,JYR,JMO,JDA,JHR,JMN,        &
                  OLAT,OLON,OLVL1,OLVL2,GOPEN,KMAP,               &
                  MAPN,PROCESS,XLEN,YLEN)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)     :: kge           ! 0=normal, 1=initialize, 2=zero concentrations
  INTEGER, INTENT(IN)     :: kproj         ! projection                    
  INTEGER, INTENT(IN)     :: level1,level2 ! range of display levels
  CHARACTER(4),INTENT(IN) :: ptype         ! ident of pollutant displayed

  REAL,    INTENT(IN)    :: parmap (9)    ! conformal map projection parameters
  INTEGER, INTENT(IN)    :: nxp,nyp       ! size of the concentration array
  INTEGER, INTENT(IN)    :: nval          ! number of contours
  INTEGER, INTENT(IN)    :: ngpt          ! number of grid points per input 

  REAL,    INTENT(IN)    :: xmax,ymax     ! global min max for contours
  REAL,    INTENT(IN)    :: tmax,tmin     ! local max min for display
  REAL,    INTENT(INOUT) :: cval(:)       ! contour array

  CHARACTER(8),INTENT(IN):: process       ! process ID suffix
  INTEGER, INTENT(IN)    :: mapn          ! current plot map number        
  INTEGER, INTENT(IN)    :: imo,ida,ihr   ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: imn,iyr       ! starting date/time of the sample
  INTEGER, INTENT(IN)    :: jmo,jda,jhr   ! stop date/time of the sample
  INTEGER, INTENT(IN)    :: jmn,jyr       ! stop date/time of the sample
  REAL,    INTENT(IN)    :: olat(:)       ! source latitude
  REAL,    INTENT(IN)    :: olon(:)       ! source longitude
  REAL,    INTENT(IN)    :: olvl1,olvl2   ! release height range
  INTEGER, INTENT(IN)    :: kmap          ! flag to indicate concen or deposit
  LOGICAL, INTENT(INOUT) :: gopen         ! Google Earth output file open
  REAL,    INTENT(IN)    :: xlen,ylen     ! size of plot in paper units

!-------------------------------------------------------------------------------

  REAL         :: XARR(5), YARR(5)
  INTEGER      :: n,ii,jj,nexp,kret,k 
  REAL         :: ctmax,ctmin 
  REAL         :: x1,y1,x2,y2 

! Rolph added for GIS
  REAL          :: xp,yp,xr,yr,clat,clon,cv,xra,yra,cvl,clvl(5)
  REAL          :: cvi(400),xpt1,ypt1
  CHARACTER(80) :: gisout,gisatt,gistmp,gislbl
  CHARACTER(8)  :: COLOR(5), COLOR1(5), COLOR2(5)
  INTEGER       :: igis,icnt,nn,nc,np,maxct,maxpt,iclvl,iprv,jlen
  INTEGER       :: npi(400),ilps,iclcnt(99)
  INTEGER       :: iyear,jyear,lenlvl,lvldat

! number of contour points
  REAL     :: XPT(9999,400), YPT(9999,400)
  REAL     :: XLAT(9999), XLON(9999),AREA
  LOGICAL  :: CLSPLY

!-------------------------------------------------------------------------------
! colors used for KMAP<4
! cyan, green, dark blue, yellow, red 
  DATA COLOR1/'c8ffff03','c800ff00','c8ff0000','c800ffff','c80000ff'/

! colors used for KMAP=4
! black, gray, yellow, orange, red
  DATA COLOR2/'c8000000','c8CCCCCC','c800ffff','c80080ff','c80000ff'/

!-------------------------------------------------------------------------------
! GIS option common block
  COMMON /GIS/GISOUT,GISATT,GISTMP,IGIS

  IF(IGIS.NE.3)RETURN

! set the contour colors based on type of output
  IF(KMAP.EQ.4)THEN
     COLOR=COLOR2
!    check to see if some contour levels are missing
     DO N=1,3
       IF(CVAL(N).EQ.0.0)THEN
         COLOR(N)='c8CCCCCC'
       END IF
     END DO
  ELSE
     COLOR=COLOR1
  END IF

! Open Google Earth format XML file (.kml) and write defaults at initial time
  IF(KGE.EQ.1)THEN
     OPEN(17,FILE=GISOUT)
     OPEN(18,FILE=GISATT)
     WRITE(17,'(A)')'<?xml version="1.0" encoding="UTF-8"?>'
     WRITE(17,'(A)')'<kml xmlns="http://opengis.net/kml/2.2">'
     WRITE(17,'(A)')'<Document>'
     WRITE(17,'(A)')'<name>NOAA HYSPLIT RESULTS</name>'
     WRITE(17,'(A)')'      <open>1</open>'
     WRITE(17,'(A)')'  <Style id="conc1">'
     WRITE(17,'(A)')'    <PolyStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(1),'</color>'
     WRITE(17,'(A,A8,A)')'      <outline>0</outline>'
     WRITE(17,'(A,A8,A)')'      <fill>1</fill>'
     WRITE(17,'(A)')'    </PolyStyle>'
     WRITE(17,'(A)')'    <LineStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(1),'</color>'
     WRITE(17,'(A)')'    </LineStyle>'
     WRITE(17,'(A)')'  </Style>'
     WRITE(17,'(A)')'  <Style id="conc2">'
     WRITE(17,'(A)')'    <PolyStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(2),'</color>'
     WRITE(17,'(A,A8,A)')'      <outline>0</outline>'
     WRITE(17,'(A,A8,A)')'      <fill>1</fill>'
     WRITE(17,'(A)')'    </PolyStyle>'
     WRITE(17,'(A)')'    <LineStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(2),'</color>'
     WRITE(17,'(A)')'    </LineStyle>'
     WRITE(17,'(A)')'  </Style>'
     WRITE(17,'(A)')'  <Style id="conc3">'
     WRITE(17,'(A)')'    <PolyStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(3),'</color>'
     WRITE(17,'(A,A8,A)')'      <outline>0</outline>'
     WRITE(17,'(A,A8,A)')'      <fill>1</fill>'
     WRITE(17,'(A)')'    </PolyStyle>'
     WRITE(17,'(A)')'    <LineStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(3),'</color>'
     WRITE(17,'(A)')'    </LineStyle>'
     WRITE(17,'(A)')'  </Style>'
     WRITE(17,'(A)')'  <Style id="conc4">'
     WRITE(17,'(A)')'    <PolyStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(4),'</color>'
     WRITE(17,'(A,A8,A)')'      <outline>0</outline>'
     WRITE(17,'(A,A8,A)')'      <fill>1</fill>'
     WRITE(17,'(A)')'    </PolyStyle>'
     WRITE(17,'(A)')'    <LineStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(4),'</color>'
     WRITE(17,'(A)')'    </LineStyle>'
     WRITE(17,'(A)')'  </Style>'
     WRITE(17,'(A)')'  <Style id="conc5">'
     WRITE(17,'(A)')'    <PolyStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(5),'</color>'
     WRITE(17,'(A,A8,A)')'      <outline>0</outline>'
     WRITE(17,'(A,A8,A)')'      <fill>1</fill>'
     WRITE(17,'(A)')'    </PolyStyle>'
     WRITE(17,'(A)')'    <LineStyle>'
     WRITE(17,'(A,A8,A)')'      <color>',COLOR(5),'</color>'
     WRITE(17,'(A)')'    </LineStyle>'
     WRITE(17,'(A)')'  </Style>'
! max square
     WRITE(17,'(A)')'  <Style id="maxv">'
     WRITE(17,'(A)')'    <LineStyle>'
     WRITE(17,'(A)')'      <color>FFFFFFFF</color>'
     WRITE(17,'(A)')'      <width>3</width>'
     WRITE(17,'(A)')'    </LineStyle>'
     WRITE(17,'(A)')'    <PolyStyle>'
     WRITE(17,'(A,A8,A)')'      <fill>0</fill>'
     WRITE(17,'(A)')'    </PolyStyle>'
     WRITE(17,'(A)')'  </Style>'

     WRITE(17,'(A)')'    <Placemark>'
     WRITE(17,'(A)')'      <name>Source Location</name>'
     WRITE(17,'(A)')'      <description><![CDATA[<pre>'
     WRITE(17,'(2(A,F9.4))')'LAT: ',OLAT(1),' LON: ',OLON(1)
     WRITE(17,'(2(A,F5.0),A)')'Released between ',OLVL1,' and ',OLVL2,'m AGL'
     WRITE(17,'(A)')'      </pre>]]>'
     WRITE(17,'(A)')'      </description>'
!    white line to source height
     WRITE(17,'(A)')'     <Style id="sorc">'
     WRITE(17,'(A)')'      <IconStyle>'
     WRITE(17,'(A)')'      <scale>0.8</scale>'
     WRITE(17,'(A)')'        <color>ff0000ff</color>'
     WRITE(17,'(A)')'       <Icon>'
!    WRITE(17,'(A)')'         <href>icon63.png</href>'
     WRITE(17,'(A)')'         <href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>'
     WRITE(17,'(A)')'       </Icon>'
     WRITE(17,'(A)')'      <hotSpot x="0.5" y="0.5" xunits="fraction" yunits="fraction"></hotSpot>'
     WRITE(17,'(A)')'      </IconStyle>'
     WRITE(17,'(A)')'      <LabelStyle>'
     WRITE(17,'(A)')'        <color>ff0000ff</color>'
     WRITE(17,'(A)')'      </LabelStyle>'
     WRITE(17,'(A)')'      <LineStyle>'
     WRITE(17,'(A)')'        <color>c8ffffff</color>'
     WRITE(17,'(A)')'        <width>2</width>'
     WRITE(17,'(A)')'      </LineStyle>'
     WRITE(17,'(A)')'     </Style>'
     WRITE(17,'(A)')'     <Point>'
     WRITE(17,'(A)')'         <extrude>1</extrude>'
     WRITE(17,'(A)')'         <altitudeMode>relativeToGround</altitudeMode>'

     IF(OLVL2.GE.10000.0)THEN
       LENLVL=5
     ELSE IF(OLVL2.GE.1000.0)THEN
       LENLVL=4
     ELSE IF(OLVL2.GE.100.0)THEN
       LENLVL=3
     ELSE IF(OLVL2.GE.10.0)THEN
       LENLVL=2
     ELSE
       LENLVL=1
     END IF

     IF(LENLVL.EQ.5)THEN
      IF(OLAT(1).LT.0.0)THEN
        WRITE(17,'(A,F9.4,A1,F8.4,A1,F7.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      ELSE
        WRITE(17,'(A,F9.4,A1,F7.4,A1,F7.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      END IF
     ELSE IF(LENLVL.EQ.4)THEN
      IF(OLAT(1).LT.0.0)THEN
        WRITE(17,'(A,F9.4,A1,F8.4,A1,F6.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      ELSE
        WRITE(17,'(A,F9.4,A1,F7.4,A1,F6.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      END IF
     ELSE IF(LENLVL.EQ.3)THEN
      IF(OLAT(1).LT.0.0)THEN
        WRITE(17,'(A,F9.4,A1,F8.4,A1,F5.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      ELSE
        WRITE(17,'(A,F9.4,A1,F7.4,A1,F5.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      END IF
     ELSE IF(LENLVL.EQ.2)THEN
      IF(OLAT(1).LT.0.0)THEN
        WRITE(17,'(A,F9.4,A1,F8.4,A1,F4.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      ELSE
        WRITE(17,'(A,F9.4,A1,F7.4,A1,F4.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      END IF
     ELSE
      IF(OLAT(1).LT.0.0)THEN
        WRITE(17,'(A,F9.4,A1,F8.4,A1,F3.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      ELSE
        WRITE(17,'(A,F9.4,A1,F7.4,A1,F3.1,A)')'         <coordinates>',OLON(1),',',  &
                   OLAT(1),',',OLVL2,'</coordinates>'
      END IF
     END IF
     WRITE(17,'(A)')'      </Point>'
     WRITE(17,'(A)')'    </Placemark>'
     WRITE(17,'(A)')'    <ScreenOverlay>'
     WRITE(17,'(A)')'      <name>HYSPLIT Information</name>'
     WRITE(17,'(A)')'      <Snippet maxLines="0"></Snippet>'
     WRITE(17,'(A)')'      <description>NOAA ARL HYSPLIT Model  http://www.arl.noaa.gov/HYSPLIT_info.php</description>'
     WRITE(17,'(A)')'      <Icon>'
     WRITE(17,'(A)')'        <href>logocon.gif</href>'
     WRITE(17,'(A)')'      </Icon>'
     WRITE(17,'(A)')'      <overlayXY x="1" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'      <screenXY x="1" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'      <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'      <size x="0" y="0" xunits="pixels" yunits="pixels"/>'
     WRITE(17,'(A)')'    </ScreenOverlay>'
     WRITE(17,'(A)')'    <ScreenOverlay>'
     WRITE(17,'(A)')'      <name>NOAA</name>'
     WRITE(17,'(A)')'      <Snippet maxLines="0"></Snippet>'
     WRITE(17,'(A)')'      <description>National Oceanic and Atomospheric Administration  http://www.noaa.gov</description>'
     WRITE(17,'(A)')'      <Icon>'
     WRITE(17,'(A)')'        <href>http://www.arl.noaa.gov/images/noaa_google.gif</href>'
     WRITE(17,'(A)')'      </Icon>'
     WRITE(17,'(A)')'      <overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'      <screenXY x="0.2" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'      <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'      <size x="0" y="0" xunits="pixels" yunits="pixels"/>'
     WRITE(17,'(A)')'    </ScreenOverlay>'
!    add a link to NOAA NWS kml weather data overlays
     WRITE(17,'(A)')'<Folder>'
     WRITE(17,'(A)')'<name>NOAA NWS kml Weather Data</name>'
     WRITE(17,'(A)')'<visibility>0</visibility>'
     WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
     WRITE(17,'(A)')'<description>'
     WRITE(17,'(A)')'http://weather.gov/gis/<br/>Access weather related Google Earth overlays from the National Weather Service.'
     WRITE(17,'(A)')'</description>'
     WRITE(17,'(A)')'</Folder>'
!    add a link to NOAA NESDIS kml smoke/fire data overlays
     WRITE(17,'(A)')'<Folder>'
     WRITE(17,'(A)')'<name>NOAA NESDIS kml Smoke/Fire Data</name>'
     WRITE(17,'(A)')'<visibility>0</visibility>'
     WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
     WRITE(17,'(A)')'<description>'
     WRITE(17,'(A)')'http://www.ssd.noaa.gov/PS/FIRE/hms.html/<br/>Access wildfire smoke Google Earth overlays from NOAA NESDIS.'
     WRITE(17,'(A)')'</description>'
     WRITE(17,'(A)')'</Folder>'
!    add a network link to EPA AIRnow kml Air Quality Index (AQI)
     WRITE(17,'(A)')'<Folder>'
     WRITE(17,'(A)')'<name>EPA AIRNow Air Quality Index (AQI)</name>'
     WRITE(17,'(A)')'<visibility>0</visibility>'
     WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
     WRITE(17,'(A)')'<description>'
     WRITE(17,'(A)')'http://www.epa.gov/airnow/today/airnow.kml<br/>Access Air Quality Index (AQI) data from EPAs AIRNow.gov.'
     WRITE(17,'(A)')'</description>'
     WRITE(17,'(A)')'</Folder>'
     GOPEN=.TRUE.
     RETURN
  END IF

! set the 4 digit year
  IYEAR=IYR+2000
  JYEAR=JYR+2000
  IF(IYR.GT.40)IYEAR=IYR+1900
  IF(JYR.GT.40)JYEAR=JYR+1900

! avoid contouring nothing, however create Google Earth label before exiting
  IF(KGE.EQ.2)THEN
     WRITE(17,'(A)')'<Folder>'
     IF(KMAP.LT.3.OR.KMAP.EQ.4)THEN
       WRITE(17,'(A)')'<name><![CDATA[<pre>Concentration '
       WRITE(17,'(A,I4,2I2.2,A,2I2.2,A)')'(Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,'UTC)</pre>]]>'
       WRITE(17,'(A)')'</name>'
       WRITE(17,'(A)')'      <description><![CDATA[<pre>'
       WRITE(17,'(2(A,I5.5),A)')'Averaged from ',LEVEL1,' to ',LEVEL2,'m'
       WRITE(17,'(A,I4,2I2.2,A1,2I2.2,A)')'Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,' UTC'
       WRITE(17,'(A)')'      </pre>]]>'
       WRITE(17,'(A)')'    </description>'
       WRITE(17,'(A)')'<ScreenOverlay>'
       WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
       WRITE(17,'(A)')'<TimeSpan>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
       WRITE(17,'(A)')'</TimeSpan>'
     ELSE
       WRITE(17,'(A)')'<name><![CDATA[<pre>Deposition '
       WRITE(17,'(A,I4,2I2.2,A,2I2.2,A)')'(Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,'UTC)</pre>]]>' 
       WRITE(17,'(A)')'</name>'
       WRITE(17,'(A)')'      <description><![CDATA[<pre>'
       WRITE(17,'(A,I4,2I2.2,A1,2I2.2,A)')'Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,' UTC'
       WRITE(17,'(A)')'      </pre>]]>'
       WRITE(17,'(A)')'    </description>'
       WRITE(17,'(A)')'<ScreenOverlay>'
       IF(GOPEN.AND.MAPN.GT.1)WRITE(17,'(A)')'    <visibility>0</visibility>'
       WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
       WRITE(17,'(A)')'<TimeSpan>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
       WRITE(17,'(A)')'</TimeSpan>'
     END IF
     WRITE(17,'(A)')'<name>Legend</name>'
     WRITE(17,'(A)')'<Icon>'
     JLEN=INDEX(PROCESS,' ')-1
     WRITE(17,'(A,I2.2,A1,A,A)')'<href>GELABEL_',MAPN,'_',PROCESS(:JLEN),'.gif</href>'
     WRITE(17,'(A)')'</Icon>'
     WRITE(17,'(A)')'<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'<rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'<size x="0" y="0" xunits="pixels" yunits="pixels"/>'
     WRITE(17,'(A)')'</ScreenOverlay>'
!    plot grid size box at maximum value location
     IF(KPROJ.EQ.4)THEN
       CALL CYL2LL(XMAX,YMAX,CLAT,CLON)
     ELSE
       CALL CXY2LL(PARMAP,XMAX,YMAX,CLAT,CLON)
     END IF
     WRITE(17,'(A)')'    <Placemark>'
     IF(KMAP.LT.3.OR.KMAP.EQ.4)THEN
        WRITE(17,'(A)')'      <name>Maximum Value Grid Cell</name>'
        WRITE(17,'(A)')'      <description><![CDATA[<pre>'
        WRITE(17,'(2(A,F9.4))')'LAT: ',CLAT,' LON: ',CLON
        WRITE(17,'(A,(1P,E7.1))')'Value: ',TMAX
        WRITE(17,'(A)')'      The square represents the location'
        WRITE(17,'(A)')'      of maximum concentration and the'
        WRITE(17,'(A)')'      size of the square represents the'
        WRITE(17,'(A)')'      concentration grid cell size.'
     ELSE
        IF(GOPEN.AND.MAPN.GT.1)WRITE(17,'(A)')'    <visibility>0</visibility>'
        WRITE(17,'(A)')'      <name>Maximum Value Grid Cell</name>'
        WRITE(17,'(A)')'      <description><![CDATA[<pre>'
        WRITE(17,'(2(A,F9.4))')'LAT: ',CLAT,' LON: ',CLON
        WRITE(17,'(A,(1P,E7.1))')'Value: ',TMAX
        WRITE(17,'(A)')'      The square represents the location'
        WRITE(17,'(A)')'      of maximum deposition and the'
        WRITE(17,'(A)')'      size of the square represents the'
        WRITE(17,'(A)')'      deposition grid cell size.'
     END IF
     WRITE(17,'(A)')'      </pre>]]>'
     WRITE(17,'(A)')'      </description>'
     WRITE(17,'(A)')'      <styleUrl>#maxv</styleUrl>'
     WRITE(17,'(A)')'      <Polygon>'
     WRITE(17,'(A)')'       <extrude>1</extrude>'
     WRITE(17,'(A)')'       <altitudeMode>relativeToGround</altitudeMode>'
     WRITE(17,'(A)')'       <outerBoundaryIs>'
     WRITE(17,'(A)')'       <LinearRing>'
     WRITE(17,'(A)')'       <coordinates>'
     IF(KPROJ.EQ.4)THEN
        CALL CYL2LL((XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(1),YARR(1))
        CALL CYL2LL((XMAX+0.5*NGPT),(YMAX-0.5*NGPT),XARR(2),YARR(2))
        CALL CYL2LL((XMAX+0.5*NGPT),(YMAX+0.5*NGPT),XARR(3),YARR(3))
        CALL CYL2LL((XMAX-0.5*NGPT),(YMAX+0.5*NGPT),XARR(4),YARR(4))
        CALL CYL2LL((XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(5),YARR(5))
     ELSE
        CALL CXY2LL(PARMAP,(XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(1),YARR(1))
        CALL CXY2LL(PARMAP,(XMAX+0.5*NGPT),(YMAX-0.5*NGPT),XARR(2),YARR(2))
        CALL CXY2LL(PARMAP,(XMAX+0.5*NGPT),(YMAX+0.5*NGPT),XARR(3),YARR(3))
        CALL CXY2LL(PARMAP,(XMAX-0.5*NGPT),(YMAX+0.5*NGPT),XARR(4),YARR(4))
        CALL CXY2LL(PARMAP,(XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(5),YARR(5))
     END IF
     DO II=1,5
        IF(XARR(II).LT.0.0)THEN
          WRITE(17,'(A,F10.5,A1,F9.5,A1,I5.5)')'        ',YARR(II),',',XARR(II),',',LEVEL2
        ELSE
          WRITE(17,'(A,F10.5,A1,F8.5,A1,I5.5)')'        ',YARR(II),',',XARR(II),',',LEVEL2
        END IF
     END DO
     WRITE(17,'(A)')'       </coordinates>'
     WRITE(17,'(A)')'       </LinearRing>'
     WRITE(17,'(A)')'       </outerBoundaryIs>'
     WRITE(17,'(A)')'      </Polygon>'
     WRITE(17,'(A)')'    </Placemark>'
     WRITE(17,'(A)')'</Folder>'
   RETURN
  END IF

  IF(KGE.EQ.0)THEN
     REWIND(16)
     WRITE(17,'(A)')'<Folder>'
     IF(GOPEN.AND.MAPN.GT.1)THEN
         WRITE(17,'(A)')'    <visibility>0</visibility>'
     ELSE
         WRITE(17,'(A)')'    <visibility>1</visibility>'
         WRITE(17,'(A)')'      <open>1</open>'
     END IF
     IF(KMAP.LT.3.OR.KMAP.EQ.4)THEN
       WRITE(17,'(A)')'<name><![CDATA[<pre>Concentration '
       WRITE(17,'(A,I4,2I2.2,A,2I2.2,A)')'(Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,'UTC)</pre>]]>'
       WRITE(17,'(A)')'</name>'
       WRITE(17,'(A)')'      <description><![CDATA[<pre>'
       WRITE(17,'(2(A,I5.5),A)')'Averaged from ',LEVEL1,' to ',LEVEL2,'m'
       WRITE(17,'(A,I4,2I2.2,A1,2I2.2,A)')'Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,' UTC'
       WRITE(17,'(A)')'      </pre>]]>'
       WRITE(17,'(A)')'    </description>'
       WRITE(17,'(A)')'<ScreenOverlay>'
       WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
       WRITE(17,'(A)')'<TimeSpan>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
       WRITE(17,'(A)')'</TimeSpan>'
     ELSE
       WRITE(17,'(A)')'<name><![CDATA[<pre>Deposition '
       WRITE(17,'(A,I4,2I2.2,A,2I2.2,A)')'(Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,'UTC)</pre>]]>'
       WRITE(17,'(A)')'</name>'
       WRITE(17,'(A)')'      <description><![CDATA[<pre>'
       WRITE(17,'(A,I4,2I2.2,A1,2I2.2,A)')'Valid:',JYEAR,JMO,JDA,' ',JHR,JMN,' UTC'
       WRITE(17,'(A)')'      </pre>]]>'
       WRITE(17,'(A)')'    </description>'
       WRITE(17,'(A)')'<ScreenOverlay>'
       IF(GOPEN.AND.MAPN.GT.1)WRITE(17,'(A)')'    <visibility>0</visibility>'
       WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
       WRITE(17,'(A)')'<TimeSpan>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
       WRITE(17,'(A)')'</TimeSpan>'
     END IF

     WRITE(17,'(A)')'<name>Legend</name>'
     WRITE(17,'(A)')'<Icon>'
     JLEN=INDEX(PROCESS,' ')-1
     WRITE(17,'(A,I2.2,A1,A,A)')'<href>GELABEL_',MAPN,'_',PROCESS(:JLEN),'.gif</href>'
     WRITE(17,'(A)')'</Icon>'
     WRITE(17,'(A)')'<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'<rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>'
     WRITE(17,'(A)')'<size x="0" y="0" xunits="pixels" yunits="pixels"/>'
     WRITE(17,'(A)')'</ScreenOverlay>'

!    read in contour levels, reorder them highest conc to lowest conc, and 
!    write them to Google Earth file
     ILPS=0
     KRET=0
     ICNT=0
     ICLCNT=0
     xloop : DO WHILE (kret.EQ.0)
        ILPS=ILPS+1
        IF(ILPS.EQ.400)THEN
          WRITE(*,*)'CONCPLOT ERROR: More than 400 closed contours.  Try increasing the number of particles.'
          STOP
        END IF
        READ(16,'(E14.7,I6)',IOSTAT=kret) CVI(ILPS),NPI(ILPS)
        IF(kret.NE.0) EXIT xloop
        IF(NPI(ILPS).GT.9999)THEN
          WRITE(*,*)'CONCPLOT ERROR: Increase 1st dimension of XPT,YPT in geplot.f - NPII(ILPS)=',NPI(ILPS)
          STOP
        END IF

        ICLVL=0
        IF(ICNT.GT.0)THEN
          DO NN=1,ICNT
            IF(CLVL(NN).EQ.CVI(ILPS))THEN
              ICLVL=NN
              ICLCNT(NN)=ICLCNT(NN)+1
            END IF
          END DO
        END IF
        IF(ICLVL.EQ.0)THEN
           ICNT=ICNT+1
           CLVL(ICNT)=CVI(ILPS)
           ICLVL=ICNT
           ICLCNT(ICLVL)=ICLCNT(ICLVL)+1
        END IF

        DO NN=1,NPI(ILPS)
           IF(ICLCNT(ICLVL).LT.99.AND.NN.LT.9999)THEN
               READ(16,'(2F10.7)') XPT(NN,ILPS),YPT(NN,ILPS)
           ELSE
               READ(16,'(2F10.7)') XPT1,YPT1
           END IF
        END DO
        IF(ICLCNT(ICLVL).GE.99.OR.NPI(ILPS).GE.9999)ILPS=ILPS-1
     END DO xloop
     ILPS=ILPS-1

     ICNT=0
     IPRV=0
     CLVL=0
     CLSPLY=.TRUE.
     DO II=1,ILPS 
        CV=CVI(II)
        NP=NPI(II)
        ICLVL=0
        IF(ICNT.GT.0)THEN
          DO NN=1,ICNT
            IF(CLVL(NN).EQ.CV)ICLVL=NN
          END DO
        END IF
        IF(ICLVL.EQ.0)THEN
           ICNT=ICNT+1
           CLVL(ICNT)=CV
           ICLVL=ICNT
        END IF
!       check for multiple polygons with same contour level
        IF(IPRV.NE.ICLVL)THEN
           IF(IPRV.GT.0)THEN
             IF(.NOT.CLSPLY)THEN
                WRITE(17,'(A)')'      </Polygon>'
               CLSPLY=.TRUE.
             END IF
             WRITE(17,'(A,I1,A)')'    </MultiGeometry>'
             WRITE(17,'(A)')'    </Placemark>'
           END IF
           IPRV=ICLVL
           WRITE(17,'(A)')'    <Placemark>'
           IF(GOPEN.AND.KMAP.EQ.3.AND.MAPN.GT.1)WRITE(17,'(A)')'    <visibility>0</visibility>'
           WRITE(17,'(A,(1P,E7.1),A)')'      <name>Contour Level: ',CV,'</name>'
           IF(KMAP.LT.3.OR.KMAP.EQ.4)THEN
             WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
             WRITE(17,'(A)')'<TimeSpan>'
             WRITE(17,'(A,I4,4(A,I2.2),A)')'<begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
             WRITE(17,'(A,I4,4(A,I2.2),A)')'<end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
             WRITE(17,'(A)')'</TimeSpan>'
           ELSE
             IF(GOPEN.AND.MAPN.GT.1)WRITE(17,'(A)')'    <visibility>0</visibility>'
             WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
             WRITE(17,'(A)')'<TimeSpan>'
             WRITE(17,'(A,I4,4(A,I2.2),A)')'<begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
             WRITE(17,'(A,I4,4(A,I2.2),A)')'<end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
             WRITE(17,'(A)')'</TimeSpan>'
           END IF
           WRITE(17,'(A)')'   '
           IF(KMAP.EQ.4)THEN
             ICLVL=ICLVL+1
           END IF
           WRITE(17,'(A,I1,A)')'       <styleUrl>#conc',ICLVL,'</styleUrl>'
           WRITE(17,'(A,I1,A)')'       <MultiGeometry>'
!          set the output level to an arbitrary height above ground in order of increasing concentration
           LVLDAT=5000*ICLVL
        END IF

        XRA=0
        YRA=0
        DO NN=1,NP
           XR=XPT(NN,II)/XLEN
           YR=YPT(NN,II)/YLEN
           XP=XR*(NXP-1.0)+1.0
           YP=YR*(NYP-1.0)+1.0
           IF(KPROJ.EQ.4)THEN
              CALL CYL2LL(XP,YP,XLAT(NN),XLON(NN))
           ELSE
              CALL CXY2LL(PARMAP,XP,YP,XLAT(NN),XLON(NN))
           END IF
        END DO
!       determine the area of polygon to see if it is inner(area>0,clockwise) or outer (area<0,counterclockwise)
        DO NN=1,NP
           AREA = .5 * (XLON(1) + XLON(NN) ) * (XLAT(1) - XLAT(NN))
           DO K = 2,NP
              AREA = AREA + .5 * (XLON(K) + XLON(K-1)) * (XLAT(K) - XLAT(K-1))
           END DO
        END DO
        IF(AREA.LE.0.0)THEN
          IF(.NOT.CLSPLY)THEN
             WRITE(17,'(A)')'      </Polygon>'
             CLSPLY=.TRUE. 
          END IF
          WRITE(17,'(A)')'      <Polygon>'
          CLSPLY=.FALSE. 
          WRITE(17,'(A)')'       <extrude>1</extrude>'
          WRITE(17,'(A)')'       <tessellate>1</tessellate>'
        END IF
        IF(AREA.LE.0.0)THEN
           WRITE(17,'(A)')'       <outerBoundaryIs>'
        ELSE
           WRITE(17,'(A)')'       <innerBoundaryIs>'
        END IF
        WRITE(17,'(A)')'       <LinearRing>'
        WRITE(17,'(A)')'       <coordinates>'
        DO NN=NP,1,-1
          IF(XLAT(NN).LT.0.0)THEN
            WRITE(17,'(F10.5,A1,F9.5,A1,I5.5)')XLON(NN),',',XLAT(NN),',',LVLDAT
          ELSE
            WRITE(17,'(F10.5,A1,F8.5,A1,I5.5)')XLON(NN),',',XLAT(NN),',',LVLDAT
          END IF
        END DO
        WRITE(17,'(A)')'       </coordinates>'
        WRITE(17,'(A)')'       </LinearRing>'
        IF(AREA.LE.0.0)THEN
           WRITE(17,'(A)')'       </outerBoundaryIs>'
        ELSE
           WRITE(17,'(A)')'       </innerBoundaryIs>'
        END IF
     END DO

     IF(.NOT.CLSPLY)THEN
         WRITE(17,'(A)')'      </Polygon>'
         CLSPLY=.TRUE. 
     END IF
     IF(ILPS.GT.0)THEN
       WRITE(17,'(A,I1,A)')'    </MultiGeometry>'
       WRITE(17,'(A)')'    </Placemark>'
     END IF
!    plot grid size box at maximum value location
     IF(KPROJ.EQ.4)THEN
       CALL CYL2LL(XMAX,YMAX,CLAT,CLON)
     ELSE
       CALL CXY2LL(PARMAP,XMAX,YMAX,CLAT,CLON)
     END IF
     WRITE(17,'(A)')'    <Placemark>'
     WRITE(17,'(A)')'      <name>Maximum Value Grid Cell</name>'
     WRITE(17,'(A)')'      <description><![CDATA[<pre>'
     WRITE(17,'(2(A,F9.4))')'LAT: ',CLAT,' LON: ',CLON
     WRITE(17,'(A,(1P,E7.1))')'Value: ',TMAX
     IF(KMAP.LT.3.OR.KMAP.EQ.4)THEN
        WRITE(17,'(A)')'      The square represents the location'
        WRITE(17,'(A)')'      of maximum concentration and the'
        WRITE(17,'(A)')'      size of the square represents the'
        WRITE(17,'(A)')'      concentration grid cell size.'
     ELSE
        WRITE(17,'(A)')'      The square represents the location'
        WRITE(17,'(A)')'      of maximum deposition and the'
        WRITE(17,'(A)')'      size of the square represents the'
        WRITE(17,'(A)')'      deposition grid cell size.'
     END IF
     WRITE(17,'(A)')'      </pre>]]>'
     WRITE(17,'(A)')'      </description>'
     IF(KMAP.LT.3.OR.KMAP.EQ.4)THEN
       WRITE(17,'(A)')'      <Snippet maxLines="0"></Snippet>'
       WRITE(17,'(A)')'      <TimeSpan>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'           <begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'           <end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
       WRITE(17,'(A)')'      </TimeSpan>'
     ELSE
       IF(GOPEN.AND.MAPN.GT.1)WRITE(17,'(A)')'    <visibility>0</visibility>'
       WRITE(17,'(A)')'<Snippet maxLines="0"></Snippet>'
       WRITE(17,'(A)')'<TimeSpan>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<begin>',IYEAR,'-',IMO,'-',IDA,'T',IHR,':',IMN,':00Z</begin>'
       WRITE(17,'(A,I4,4(A,I2.2),A)')'<end>',JYEAR,'-',JMO,'-',JDA,'T',JHR,':',JMN,':00Z</end>'
       WRITE(17,'(A)')'</TimeSpan>'
     END IF
     WRITE(17,'(A)')'      <styleUrl>#maxv</styleUrl>'
     WRITE(17,'(A)')'      <Polygon>'
     WRITE(17,'(A)')'       <extrude>1</extrude>'
     WRITE(17,'(A)')'       <tessellate>1</tessellate>'
     WRITE(17,'(A)')'       <altitudeMode>relativeToGround</altitudeMode>'
     WRITE(17,'(A)')'       <outerBoundaryIs>'
     WRITE(17,'(A)')'       <LinearRing>'
     WRITE(17,'(A)')'       <coordinates>'
     IF(KPROJ.EQ.4)THEN
        CALL CYL2LL((XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(1),YARR(1))
        CALL CYL2LL((XMAX+0.5*NGPT),(YMAX-0.5*NGPT),XARR(2),YARR(2))
        CALL CYL2LL((XMAX+0.5*NGPT),(YMAX+0.5*NGPT),XARR(3),YARR(3))
        CALL CYL2LL((XMAX-0.5*NGPT),(YMAX+0.5*NGPT),XARR(4),YARR(4))
        CALL CYL2LL((XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(5),YARR(5))
     ELSE
        CALL CXY2LL(PARMAP,(XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(1),YARR(1))
        CALL CXY2LL(PARMAP,(XMAX+0.5*NGPT),(YMAX-0.5*NGPT),XARR(2),YARR(2))
        CALL CXY2LL(PARMAP,(XMAX+0.5*NGPT),(YMAX+0.5*NGPT),XARR(3),YARR(3))
        CALL CXY2LL(PARMAP,(XMAX-0.5*NGPT),(YMAX+0.5*NGPT),XARR(4),YARR(4))
        CALL CXY2LL(PARMAP,(XMAX-0.5*NGPT),(YMAX-0.5*NGPT),XARR(5),YARR(5))
     END IF
     DO II=1,5
        IF(XARR(II).LT.0.0)THEN
          WRITE(17,'(A,F10.5,A1,F9.5,A1,I5.5)')'        ',YARR(II),',',XARR(II),',',LEVEL2
        ELSE
          WRITE(17,'(A,F10.5,A1,F8.5,A1,I5.5)')'        ',YARR(II),',',XARR(II),',',LEVEL2
        END IF
     END DO
     WRITE(17,'(A)')'       </coordinates>'
     WRITE(17,'(A)')'       </LinearRing>'
     WRITE(17,'(A)')'       </outerBoundaryIs>'
     WRITE(17,'(A)')'      </Polygon>'
     WRITE(17,'(A)')'    </Placemark>'

     WRITE(17,'(A)')'</Folder>'
     CLOSE (16)
     IF(.NOT.GOPEN)THEN
        GOPEN=.TRUE.
     END IF
  END IF
  RETURN

END SUBROUTINE geplot
