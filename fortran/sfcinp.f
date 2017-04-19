!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SFCINP           SurFaCe data INPut reads lat/lon files
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SURFACE DATA INPUT READS A LAT/LON BASED SFC CHARACTERISTICS FILE
!   CURRENT COMPILED VERSION ASSUMES A FILE WITH 1-DEG RESOLUTION WITH
!   THE FIRST RECORD STARTING AT THE NORTHWEST CORNER (CNTR 179.5W,89.5N)
!   RETURNS THE ROUGHNESS LENGTH AND LAND-USE AT SPECIFIED POINT.
!
!   Land-Use Values are defined as follows:
!      1-urban          2-agricultural      3-dry range      4-deciduous
!      5-coniferous     6-mixed forest      7-water          8-desert
!      9-wetlands      10-mixed range      11-rocky
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 Apr 1997 (RRD)
!                 14 Jun 1999 (RRD) - generalized data file input
!                 21 Oct 1999 (RRD) - LREC added to save statement
!                 26 Jul 2000 (RRD) - open statements as read only
!                 29 Sep 2000 (RRD) - fortran90 upgrade
!                 05 Dec 2001 (RRD) - added terrain height file read
!                 24 Jun 2002 (RRD) - initial default dir is local
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 03 Apr 2008 (RRD) - dynamic array allocation
!                 03 Jun 2008 (RRD) - embedded blanks dir/file
!                 07 Aug 2008 (RRD) - terminate if no ASCDATA.CFG
!
! USAGE:  CALL SFCINP(CLAT,CLON,ZNOT,LUSE,HSET,HGTS)
!
!   INPUT ARGUMENT LIST:       see below
!   OUTPUT ARGUMENT LIST:      see below
!   INPUT FILES:               unit 60 landuse categories (LANDUSE.ASC)
!                              unit 62 roughness length (ROUGLEN.ASC)
!                              unit 64 terrain height (TERRAIN.ASC)
!   OUTPUT FILES:              none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE SFCINP(CLAT,CLON,ZNOT,LUSE,HSET,HGTS)

  USE funits

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,     INTENT(IN)    :: clat,clon        ! Lat/Lon of required point
  REAL,     INTENT(OUT)   :: znot             ! Aerodynamic rougness length (m)
  INTEGER,  INTENT(OUT)   :: luse             ! Land-use categories (1-11)
  LOGICAL,  INTENT(IN)    :: hset             ! Read terrain file flag
  REAL,     INTENT(OUT)   :: hgts             ! terrain height value (m)

!-------------------------------------------------------------------------------
! internal variables     
!-------------------------------------------------------------------------------

  LOGICAL                 :: tfile = .FALSE.  
  LOGICAL                 :: rfile = .FALSE.  
  LOGICAL                 :: qfile = .FALSE.  
  LOGICAL                 :: hfile = .FALSE.
  INTEGER                 :: krec  = 0
  REAL                    :: alatb, alonl, dlat, dlon
  INTEGER                 :: nlat, nlon, lrec, luse0
  INTEGER                 :: jrec, klen, k, kk 
  REAL                    :: znot0, hgts0
  CHARACTER(4)            :: dval
  CHARACTER(80)           :: fdir

  CHARACTER(1),  ALLOCATABLE :: tdat(:),rdat(:),hdat(:)

!-------------------------------------------------------------------------------

  SAVE QFILE,TFILE,RFILE,HFILE,KREC,ALATB,ALONL,DLAT,DLON, &
       NLAT,NLON,LUSE0,ZNOT0,HGTS0,LREC,TDAT,RDAT,HDAT

!-------------------------------------------------------------------------------
! open the files only once
!-------------------------------------------------------------------------------

  IF(.NOT.QFILE)THEN

!    the configuration file

     INQUIRE(FILE='ASCDATA.CFG',EXIST=QFILE)
     IF(QFILE)THEN
        OPEN(KF41,FILE='ASCDATA.CFG',ACTION='READ',recl=1)
     ELSE
        INQUIRE(FILE='../bdyfiles/ASCDATA.CFG',EXIST=QFILE)
        IF(QFILE)OPEN(KF41,FILE='../bdyfiles/ASCDATA.CFG',ACTION='READ')
     END IF

     IF(QFILE)THEN
        WRITE(KF21,*)' NOTICE sfcinp: reading ASCDATA.CFG'
        READ(KF41,*)ALATB,ALONL
        READ(KF41,*)DLAT,DLON
        READ(KF41,*)NLAT,NLON
        READ(KF41,*)LUSE0
        READ(KF41,*)ZNOT0
        READ(KF41,*)FDIR
        CLOSE(KF41)
        FDIR=ADJUSTL(FDIR)
        KLEN=LEN_TRIM(FDIR)

!       allocate array space based upon input data resolution
!       ascii data record length = number points x 4 bytes plus one
        LREC = 1+nlon*4  
        ALLOCATE (TDAT(LREC),RDAT(LREC),HDAT(LREC))

     ELSE
        WRITE(KF21,*)' NOTICE sfcinp: no ASCDATA.CFG, using default'
!       lower left corner of the file
        ALATB=-90.0
        ALONL=-180.0
!       incremental spacing
        DLAT=1.0
        DLON=1.0
!       number of points in each direction
        NLAT=180
        NLON=360
!       content default values
        LUSE0=2
        ZNOT0=0.2
        HGTS0=0.0
!       data file directory
        FDIR='../bdyfiles/'
        WRITE(KF21,*)'*ERROR* sfcinp: ASCDATA.CFG file not found!'
        WRITE(KF21,*)'Required in local directory or ../bdyfiles/'
        WRITE(*,*)'*ERROR* sfcinp: ASCDATA.CFG file not found!'
        WRITE(*,*)'  See MESSAGE file for more information  '
        STOP 
     END IF


!    land use file

     INQUIRE(FILE='LANDUSE.ASC',EXIST=QFILE)
     IF(QFILE)THEN
        OPEN(KF41,FILE='LANDUSE.ASC',ACCESS='DIRECT',                       &
             RECL=LREC,FORM='UNFORMATTED',ACTION='READ')
        TFILE=.TRUE.
     ELSE
        INQUIRE(FILE=FDIR(1:KLEN)//'LANDUSE.ASC',EXIST=QFILE)
        IF(QFILE)THEN
           OPEN(KF41,FILE=FDIR(1:KLEN)//'LANDUSE.ASC',ACCESS='DIRECT',       &
                RECL=LREC,FORM='UNFORMATTED',ACTION='READ')
           TFILE=.TRUE.
        ELSE
           WRITE(KF21,*)' NOTICE sfcinp: LANDUSE.ASC file not found'
           WRITE(KF21,*)' On default directory:',FDIR(1:KLEN)
           WRITE(KF21,*)' Using default value category 2 agricultural'
        END IF
     END IF

!    roughness length file

     INQUIRE(FILE='ROUGLEN.ASC',EXIST=QFILE)
     IF(QFILE)THEN
        OPEN(KF42,FILE='ROUGLEN.ASC',ACCESS='DIRECT',                        &
             RECL=LREC,FORM='UNFORMATTED',ACTION='READ')
        RFILE=.TRUE.
     ELSE
        INQUIRE(FILE=FDIR(1:KLEN)//'ROUGLEN.ASC',EXIST=QFILE)
        IF(QFILE)THEN
           OPEN(KF42,FILE=FDIR(1:KLEN)//'ROUGLEN.ASC',ACCESS='DIRECT',       &
                RECL=LREC,FORM='UNFORMATTED',ACTION='READ')
           RFILE=.TRUE.
        ELSE
           WRITE(KF21,*)' NOTICE sfcinp: ROUGLEN.ASC file not found'
           WRITE(KF21,*)' On default directory:',FDIR(1:KLEN)
           WRITE(KF21,*)' Using default value 0.20 m'
        END IF
     END IF

!    terrain height file   

     IF(HSET)THEN
!       only required for pressure sigma with no terrain
        INQUIRE(FILE='TERRAIN.ASC',EXIST=QFILE)
        IF(QFILE)THEN
           OPEN(KF43,FILE='TERRAIN.ASC',ACCESS='DIRECT',                     &
                RECL=LREC,FORM='UNFORMATTED',ACTION='READ')
           HFILE=.TRUE.
        ELSE
           INQUIRE(FILE=FDIR(1:KLEN)//'TERRAIN.ASC',EXIST=QFILE)
           IF(QFILE)THEN
              OPEN(KF43,FILE=FDIR(1:KLEN)//'TERRAIN.ASC',ACCESS='DIRECT',    &
                   RECL=LREC,FORM='UNFORMATTED',ACTION='READ')
              HFILE=.TRUE.
           ELSE
              WRITE(KF21,*)' NOTICE sfcinp: TERRAIN.ASC file not found'
              WRITE(KF21,*)' On default directory:',FDIR(1:KLEN)
              WRITE(KF21,*)' Using default value 0 m'
           END IF
        END IF
     END IF

!    set logical file test so that files are not opened again
     QFILE=.TRUE.

  END IF

!-------------------------------------------------------------------------------
! load data from file as required
!-------------------------------------------------------------------------------

! default values
  LUSE=LUSE0
  ZNOT=ZNOT0
  HGTS=HGTS0

  IF(RFILE.OR.TFILE.OR.HFILE)THEN
!    determine required record number based upon latitude
     JREC=NLAT-INT((CLAT-ALATB)/DLAT)
     JREC=MIN(MAX(1,JREC),NLAT)

!    compute the byte number of first digit from longitude
     KK=INT((CLON-ALONL)/DLON)+1
     KK=1+(MIN(MAX(1,KK),NLON)-1)*4

     IF(JREC.NE.KREC)THEN
!       load new data into buffers if record different from previous
        KREC=JREC
        IF(TFILE) READ(KF41,REC=JREC) TDAT   
        IF(RFILE) READ(KF42,REC=JREC) RDAT   
        IF(HFILE) READ(KF43,REC=JREC) HDAT   
     END IF

     IF(TFILE)THEN                  
        DVAL=TDAT(KK)//TDAT(KK+1)//TDAT(KK+2)//TDAT(KK+3)
        READ(DVAL,'(I4)') LUSE 
     END IF

     IF(RFILE)THEN                  
        DVAL=RDAT(KK)//RDAT(KK+1)//RDAT(KK+2)//RDAT(KK+3)
        READ(DVAL,'(F4.0)') ZNOT  
        ZNOT=ZNOT*0.01             ! cm to meters
     END IF

     IF(HFILE)THEN                  
        DVAL=HDAT(KK)//HDAT(KK+1)//HDAT(KK+2)//HDAT(KK+3)
        READ(DVAL,'(F4.0)') HGTS 
     END IF
  END IF

END SUBROUTINE sfcinp
