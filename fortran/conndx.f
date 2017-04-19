!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONNDX           SET CONCENTRATION VALUES FOR PLOTTING
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            READS THE HEADER INFORMATION IN THE CONCENTRATION FILE
!            AND SETS RELATED PROGRAM VARIABLES.  THE ROUTINE THEN
!            READS ALL THE CONCENTRATION RECORDS TO DETERMINE THE MAXIMUM
!            SIZE OF THE CLOUD FOR THE DURATION OF THE SIMULATION
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 09 Dec 1998 (RRD)
!                 28 Mar 2000 (RRD) - added map centroid computation
!                 03 Nov 2000 (RRD) - grid limits tests and dateline tests
!                 20 Nov 2000 (RRD) - fortran90 upgrade
!                 01 Dec 2000 (RRD) - concentration packing
!                 21 Feb 2001 (RRD) - improved longitude span check
!                 29 Oct 2001 (RRD) - save release height ranges
!                 02 Nov 2001 (RRD) - standardization
!                 16 Nov 2001 (RRD) - max min contour by level
!                 14 Jan 2002 (RRD) - revised map limits test
!                 28 Feb 2002 (RRD) - map center based upon mean position
!                 17 Jul 2002 (RRD) - major revison
!                 20 Nov 2003 (RRD) - added minutes field to header
!                 04 Dec 2003 (RRD) - include to module
!                 22 Dec 2004 (RRD) - cpack=2 test 
!                 02 Jun 2008 (BS)  - make map even if all concentrations ZERO
!
! USAGE: CALL CONNDX(CONC,TEMP,HEIGHT,MODEL,IBYR,IBMO,IBDA,IBHR,IBMN,
!             OLAT,OLON,OLVL1,OLVL2,DLAT,DLON,CLAT,CLON,IDENT,NMAPS,
!             CMAX,CMIN,CPACK,NMAP1,NMAP2,IZRO)
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

SUBROUTINE CONNDX(CONC,TEMP,HEIGHT,MODEL,IBYR,IBMO,IBDA,IBHR,IBMN,         &
                  OLAT,OLON,OLVL1,OLVL2,DLAT,DLON,CLAT,CLON,IDENT,NMAPS,   &
                  CMAX,CMIN,CPACK,NMAP1,NMAP2,IZRO)

  USE mapbox

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,        INTENT(OUT)   :: conc(:,:,:) ! concen grid on lat/lons
  REAL,        INTENT(OUT)   :: temp(:,:)   ! working array space
  INTEGER,     INTENT(OUT)   :: height(:)   ! levels of concentration data
  CHARACTER(4),INTENT(OUT)   :: model       ! identification of meteo data
  INTEGER,     INTENT(OUT)   :: ibyr(:)     ! starting year of emissions
  INTEGER,     INTENT(OUT)   :: ibmo(:)     ! starting month 
  INTEGER,     INTENT(OUT)   :: ibda(:)     ! starting day 
  INTEGER,     INTENT(OUT)   :: ibhr(:)     ! starting hour 
  INTEGER,     INTENT(OUT)   :: ibmn(:)     ! starting minutes   
  REAL,        INTENT(OUT)   :: olat(:)     ! starting location 
  REAL,        INTENT(OUT)   :: olon(:)     ! starting location
  REAL,        INTENT(OUT)   :: olvl1,olvl2 ! pollutant release height range
  REAL,        INTENT(OUT)   :: dlat,dlon   ! resolution of grid
  REAL,        INTENT(OUT)   :: clat,clon   ! center of the concen grid
  CHARACTER(4),INTENT(OUT)   :: ident(:)    ! identification of each pollutant
  INTEGER,     INTENT(OUT)   :: nmaps       ! number of time periods in file
  REAL,        INTENT(OUT)   :: cmin(:)     ! max and min concen values
  REAL,        INTENT(OUT)   :: cmax(:)     ! max and min concen values
  INTEGER,     INTENT(OUT)   :: cpack       ! concentration data packing method
  INTEGER,     INTENT(INOUT) :: nmap1,nmap2 ! range of time periods to plot  
  INTEGER,     INTENT(IN)    :: izro        ! izro=1 make map if all conc ZERO
!-------------------------------------------------------------------------------
! internally defined variables
!-------------------------------------------------------------------------------

  INTEGER      :: nlat,nlon         ! number of actual pnts in grid
  INTEGER      :: nloc              ! actual numb of source locations
  INTEGER      :: ntyp              ! number of pollutants in file
  INTEGER      :: nlvl              ! number of concen levels in file

  CHARACTER(4) :: PTYPE*4
  INTEGER(2)   :: IP,JP
  INTEGER      :: KRET

  REAL         :: rcon,olvl,plon,plat
  INTEGER      :: n,np,ii,jj,kp,kl,ig,kr,jg,kk 
  INTEGER      :: level,nxyp,numb,nhit
  INTEGER      :: myr,mmo,mhr,mda,mfh 
  INTEGER      :: jmo,jyr,jda,jmn,jhr,jfh 
  INTEGER      :: iyr,ida,imo,ihr,ifh,imn 

!-------------------------------------------------------------------------------
! determine input concentration grid/data characteristics
!-------------------------------------------------------------------------------

  100 CONTINUE ! special restart if intial time limits incorrect

! meteo file type and initial position date
  READ(10,IOSTAT=KRET)MODEL,MYR,MMO,MDA,MHR,MFH,NLOC,CPACK
  IF(KRET.NE.0)CPACK=0  ! old style file with all grid points

  IF(CPACK.GT.1)THEN
     WRITE(*,*)'*ERROR* conndx: unsupported packing option -',CPACK
     STOP 900
  END IF

! release starting information (added minutes 20 Nov 2003)
  OLVL1=99999.0
  OLVL2=0.0
  DO N=1,NLOC
     READ(10,IOSTAT=KRET) IBYR(N),IBMO(N),IBDA(N),IBHR(N),    &
                          OLAT(N),OLON(N),OLVL,IBMN(N)
     IF(KRET.NE.0)IBMN(N)=0  ! old style default to zero minutes
     OLVL1=MIN(OLVL1,OLVL)
     OLVL2=MAX(OLVL2,OLVL)
  END DO

! number of points, delta increment, lower left corner
  READ(10)   NLAT, NLON, DLAT, DLON, CLAT, CLON

! vertical grid index record, number of levels and heights
  READ(10) NLVL, (HEIGHT(KK),KK=1,NLVL)

! pollutant identification record
  READ(10) NTYP, (IDENT(KK),KK=1,NTYP)

!-------------------------------------------------------------------------------
! read through concentration gridded data records to determine limits
!-------------------------------------------------------------------------------

! save max min for all time periods on file
  CMAX=0.0
  CMIN=1.0E+25

  NMAPS=0
  DO WHILE (NMAPS.LT.NMAP2)

!    sample start and stop times
     READ(10,END=500,ERR=500)IYR,IMO,IDA,IHR,IMN,IFH
     READ(10,END=500,ERR=500)JYR,JMO,JDA,JHR,JMN,JFH
     NMAPS=NMAPS+1

     DO KP=1,NTYP
     DO KL=1,NLVL

        IF(CPACK.EQ.1)THEN
!          condensed style with only nonzero grid points written to file 
           TEMP = 0.0
           READ(10,IOSTAT=kret)PTYPE,LEVEL,NXYP,(IP,JP,TEMP(IP,JP),NP=1,NXYP)

        ELSE
!          old style record with all grid points written to file
           READ(10,END=500,ERR=500)PTYPE,LEVEL,     &
               ((TEMP(II,JJ),II=1,NLON),JJ=1,NLAT)
        END IF

!       fill array to determine mapping limits 
        IF(NMAPS.GE.NMAP1)THEN
           DO II=1,NLON
           DO JJ=1,NLAT
              RCON=TEMP(II,JJ)
              IF(RCON.GT.0.0)CMIN(KL)=AMIN1(RCON,CMIN(KL))
              CMAX(KL)=AMAX1(RCON,CMAX(KL))
              CONC(II,JJ,1)=CONC(II,JJ,1)+RCON
           END DO
           END DO
        END IF

     END DO
     END DO

  END DO

  500 REWIND(10)

! check time limits for mapping
  NMAP2=NMAPS
  IF(NMAP2.GE.1.AND.NMAP1.GT.NMAP2)THEN
!    special case where fewer maps in file than selected
     NMAP1=1
     GOTO 100
  END IF

! reposition file to start of sequential gridded output past index records
  DO KR=1,(NLOC+4)
     READ(10,END=600,ERR=600)
  END DO
  600 CONTINUE

!-------------------------------------------------------------------------------
! define internal plume position counter grid for map optimization
!-------------------------------------------------------------------------------

! for small maps set up a finer grid, otherwise use module defaults
  IF((NLAT-1)*DLAT.LE.2.0.AND.(NLON-1)*DLON.LE.2.0)THEN
     GINC=0.10
     GCLAT=CLAT
     GCLON=CLON
     KLAT=NLAT*DLAT/GINC
     KLON=NLON*DLON/GINC
  ELSEIF((NLAT-1)*DLAT.LE.5.0.AND.(NLON-1)*DLON.LE.5.0)THEN
     GINC=0.20
     GCLAT=CLAT
     GCLON=CLON
     KLAT=NLAT*DLAT/GINC
     KLON=NLON*DLON/GINC
  END IF

  ALLOCATE(LATLON(KLAT,KLON),STAT=KRET)
  IF(KRET.NE.0)THEN
     WRITE(*,*)'Memory allocation error for internal grid:',KLAT,KLON
     WRITE(*,*)'Rerun model with smaller concentration grid'  
     STOP 900
  END IF

! zero global map plume position counter
  LATLON = 0

!-------------------------------------------------------------------------------
! sum non-zero plume hits to cells
!-------------------------------------------------------------------------------

  NHIT=0
  DO II=1,NLON
     PLON=(II-1)*DLON+CLON
     IF(PLON.LT.0.0.AND.GCLON.GE.0.0)PLON=PLON+360.0

     DO JJ=1,NLAT
        PLAT=(JJ-1)*DLAT+CLAT

        IF(CONC(II,JJ,1).GT.0.0)THEN
!          convert to global coordinates
           IG=INT((PLAT-GCLAT)/GINC)+1
           JG=INT((PLON-GCLON)/GINC)+1
!          limits check (rrd 2/21/01)
           IG=MAX(1,MIN(IG,KLAT))
           JG=MAX(1,MIN(JG,KLON))
!          count hits in cells
           LATLON(IG,JG)=LATLON(IG,JG)+1
           NHIT=NHIT+1
        END IF
     END DO
  END DO

  IF(NHIT.EQ.0)THEN
     IF(IZRO.EQ.1)THEN
        WRITE(*,*)'All concentrations ZERO'
     ELSE 
        WRITE(*,*)'All concentrations ZERO - no maps'
        STOP 900
     END IF
  END IF

END SUBROUTINE conndx
