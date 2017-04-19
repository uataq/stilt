!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  IERVOC           IER VOC emissions adjustment
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   NATURAL EMISSIONS OF ISOPRENE ARE COMPUTED FOR ALL FOREST LAND-USE 
!   CATEGORIES. THIS PROGRAM IS CALLED HOURLY BEFORE EMSGRD BECAUSE IT 
!   ADJUSTS THE VOC EMISSION ARRAY FOR SOLAR ANGLE AND TEMPERATURE
!   CORRECTED ISOPRENE. NO INITIAL VALUES ARE ALSO COMPUTED FROM NOX
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Feb 1998 (RRD)
!                 22 Jul 1999 (RRD) - subroutine name change
!                 24 Jul 2003 (RRD) - new version created from ierems
!                 02 Apr 2004 (RRD) - generic file unit numbers
!
! USAGE:  CALL IERVOC(SPOT,JET,NPVAL,NQLON,NQLAT,QDLON,QDLAT,QAREA,POLID)
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

SUBROUTINE IERVOC(SPOT,JET,NPVAL,NQLON,NQLAT,QDLON,QDLAT,QAREA,POLID)

  USE funits
  USE metval
  use module_defgrid ! meteorology file

  IMPLICIT NONE

  INCLUDE 'DEFSPOT.INC' ! multiple source information

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(rset),   INTENT(IN)    :: spot(:)          ! source characteristics
  INTEGER,      INTENT(IN)    :: jet              ! current minutes
  INTEGER,      INTENT(IN)    :: npval            ! number gridded pollutants
  INTEGER,      INTENT(IN)    :: nqlon, nqlat     ! number of grid points
  REAL,         INTENT(IN)    :: qdlon, qdlat     ! grid spacing
  REAL,         INTENT(INOUT) :: qarea (:,:,:,:)  ! gridded emissions array 
  CHARACTER(4), INTENT(IN)    :: polid (:)        ! emissions pollutant id

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL, PARAMETER :: PI     = 3.14159   ! deg per radian
  REAL, PARAMETER :: DEGPRD = 180.0/PI  ! deg per radian

  INTEGER         :: KISO=0 
  INTEGER         :: KNOX=0 
  INTEGER         :: KNOB=0 
  INTEGER         :: KNO2=0 

  INTEGER         :: iyr,imo,ida,ihr,imn,luse
  INTEGER         :: kp,kd,kg,kt,ii,jj,ll
  REAL            :: xp,yp,ea,sea,clat,clon,area,temp 


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
! Find the index number of required pollutants. Those that have not been defined
! in the emission.txt file must be defined in the main program prior to calling
! this routine.
!-------------------------------------------------------------------------------

  DO KD=1,NPVAL
     IF(POLID(KD).EQ.'ISOP') KISO=KD
     IF(POLID(KD).EQ.'NOx ') KNOX=KD
     IF(POLID(KD).EQ.'NO# ') KNOB=KD
     IF(POLID(KD).EQ.'NO2 ') KNO2=KD
  END DO
  IF(KISO.EQ.0.OR.KNOX.EQ.0.OR.KNOB.EQ.0)THEN
     WRITE(KF21,*)'*ERROR* iervoc: pollutants undefined'
     STOP
  END IF

! compute current hour and convert to array index (depends on #q per day)
  CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)
  LL=IHR+1               ! index 1 = average 0000 to 0100 UTC
  KP=POINT(2)            ! set temporal pointer to the next time period
  KG=1                   ! initial grid index set in main program
  KT=1                   ! file time period index

!-------------------------------------------------------------------------------
! check each grid point
!-------------------------------------------------------------------------------

  DO JJ=1,NQLAT

!    compute position from lower left corner point (1,1 point)
!    assume emissions at center of grid cell (from SW corner)
     CLAT=FLOAT(JJ-1)*QDLAT+SPOT(1)%OLAT+QDLAT/2.0

!    area source (111 km / deg - squared)
     AREA=12321.0*QDLAT*QDLON*COS(CLAT/DEGPRD)

  DO II=1,NQLON

     CLON=FLOAT(II-1)*QDLON+SPOT(1)%OLON+QDLON/2.0

!    convert to grid units
     IF(GRID(KG,KT)%LATLON)THEN
        CALL GBL2XY(KG,KT,CLAT,CLON,XP,YP)
     ELSE
        CALL CLL2XY(GRID(KG,KT)%GBASE,CLAT,CLON,XP,YP)
     END IF

     LUSE=LU(NINT(XP),NINT(YP),KG)
     TEMP=A(NINT(XP),NINT(YP),1,KP,KG)

!    select grid points with forest variation
     IF((LUSE.GE.4.AND.LUSE.LE.6).OR.LUSE.EQ.10)THEN
!       emissions adjusted by local sun angle
        CALL SUNANG(JET,CLAT,CLON,EA,SEA)

!       assume isoprene (index=4) maximum equals 2 kg/(km^2 hr)
        QAREA(II,JJ,KISO,LL)=MAX(0.0, 2.0*AREA*SEA)
     END IF

!    temperature adjustment for isoprene
!    from jacob et al (JGR,1993,14797-14813)
     QAREA(II,JJ,KISO,LL)=QAREA(II,JJ,KISO,LL)*EXP(0.096*(TEMP-298.0))

!    NOX is given in NO2 (mw=46) molar units
!    NO (30) is assumed to be 90% of NOx after correction for MW (30/46)
     QAREA(II,JJ,KNOB,LL)=0.58696*QAREA(II,JJ,KNOX,LL)

!    add NO2 emissions if required (GRS module)
     IF(KNO2.NE.0) QAREA(II,JJ,KNO2,LL)=0.10*QAREA(II,JJ,KNOX,LL)

! number of sources loop
  END DO
  END DO

END SUBROUTINE iervoc
