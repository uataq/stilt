!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METLVL           METeorological LeVeL analysis 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:03-02-05
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   ANALYZES THE METEOROLOGICAL DATA AND DETERMINES THE APPROPRIATE  
!   INTERNAL GRID TO USE SO THAT THERE ARE SUFFICIENT LEVELS TO 
!   INTERPOLATE ALL THE METEOROLOGICAL INPUT DATA LEVELS TO AN INTERNAL
!   LEVEL WITHOUT SKIPPING INPUT DATA DUE TO INSUFFICIENT VERTICAL
!   RESOLUTION OR INSUFFICIENT NUMBER OF DATA LEVELS.  THE NUMBER OF
!   LEVELS IS DETERMINED FROM THE TOP OF THE MODEL DOMAIN DEFINED IN THE
!   INPUT CONTROL FILE IN SUBROUTINE RUNSET. 
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 04 Feb 2003 (RRD) - initial version
!                 13 Feb 2003 (RRD) - added max meteo data height
!                 09 Oct 2007 (RRD) - nintg loop termination
!
! USAGE:  CALL METLVL(NGRD,ZDATA,NLVL,KSFC,SFCL)
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

SUBROUTINE METLVL(NGRD,ZDATA,NLVL,KSFC,SFCL)

  use module_defgrid

  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: NGRD  ! number of meteo grids      
  REAL,    INTENT(IN)   :: ZDATA ! maximum height of input data
  INTEGER, INTENT(OUT)  :: NLVL  ! internal grid number levels
  INTEGER, INTENT(OUT)  :: KSFC  ! index of the surface layer
  REAL,    INTENT(OUT)  :: SFCL  ! height of the surface layer

  REAL     :: SFCP = 1013.0      ! assumed surface pressure

!  LOGICAL  :: DIAG = .false.     ! diagnostic testing 
 LOGICAL  :: DIAG = .true.      ! diagnostic testing 
  INTEGER  :: I,KK,KG            ! indicies
  INTEGER  :: LL                 ! input data level index
  REAL     :: DELZ               ! height increment per hPA
  REAL     :: AA,BB,CC,DIST      ! polynomial grid coefficients
  REAL     :: PLEVEL,ZLEVEL      ! pressure and height values
  REAL     :: OFFSET,PSIGMA      ! hybrid coordinate factors


  INTEGER,PARAMETER   :: NINTG = 5
  REAL :: AZ(NINTG) = (/ 30., 20., 15., 10., 5. /)
  REAL :: BZ(NINTG) = (/-25.,-15., -5.,-10., 5. /)
  REAL :: CZ(NINTG) = (/  5.,  5.,  0., 10., 0. /)

! meters per hPa by 100 hPa intervals (100 to 900)
  REAL :: ZPH1(9) = (/ 17.98,14.73,13.09,11.98,11.15,10.52,10.04,9.75,9.88 /) 

! meters per hPa by 10 hPa intervals (10 to 90)
  REAL :: ZPH2(9) = (/ 31.37,27.02,24.59,22.92,21.65,20.66,19.83,19.13,18.51 /)

! see defgrid.inc for common variables description

! vertical coordinate system coefficients
  COMMON /ZZTOKK/ AA,BB,CC

!-------------------------------------------------------------------
! determine which grid has the most levels

  NLVL=0
  DO KK=1,NGRD
     IF(GRID(KK,1)%NZ.GT.NLVL)THEN
        NLVL=GRID(KK,1)%NZ
        KG=KK
     END IF
  END DO
  IF(DIAG) WRITE(*,*)'Max number of input lvls: ',NLVL
  IF(DIAG) WRITE(*,*)'Vertical coordinate flag: ',DREC(KG,1)%Z_FLAG

!-------------------------------------------------------------------
! Find the input data level index above the selected domain top
! The internal grid should have at least this number of levels

  LL=0
  DELZ=10.0
  ZLEVEL=0.0
  DO WHILE (ZLEVEL.LE.ZDATA.AND.LL.LT.NLVL)
     LL=LL+1

     IF(DREC(KG,1)%Z_FLAG.EQ.1)THEN
!       pressure sigma levels
        OFFSET=GRID(KG,1)%DUMMY
        PLEVEL=OFFSET+(SFCP-OFFSET)*DREC(KG,1)%HEIGHT(LL)

     ELSEIF(DREC(KG,1)%Z_FLAG.EQ.2)THEN
!       absolute pressure units
        PLEVEL=DREC(KG,1)%HEIGHT(LL)
        IF(LL.EQ.1)PLEVEL=SFCP

     ELSEIF(DREC(KG,1)%Z_FLAG.EQ.3)THEN
!       terrain following height units
        ZLEVEL=DREC(KG,1)%HEIGHT(LL)

     ELSEIF(DREC(KG,1)%Z_FLAG.EQ.4)THEN
!       ecmwf hybrid coordinate system
        OFFSET=INT(DREC(KG,1)%HEIGHT(LL))
        PSIGMA=DREC(KG,1)%HEIGHT(LL)-OFFSET
        PLEVEL=SFCP*PSIGMA+OFFSET
        IF(LL.EQ.1)PLEVEL=SFCP
     END IF

!    convert pressure to approximate height
     IF(DREC(KG,1)%Z_FLAG.NE.3)THEN
        IF(PLEVEL.GE.100.0)THEN
           DELZ=ZPH1(MIN(9,INT(PLEVEL/100.0)))
        ELSE
           DELZ=ZPH2(MAX(1,MIN(9,INT(PLEVEL/10.0))))
        END IF
        ZLEVEL=(SFCP-PLEVEL)*DELZ
     END IF
     IF(DIAG) WRITE(*,*)'Input data level (k,ppp,zzz,dz): ', &
              LL,PLEVEL,ZLEVEL,DELZ
  END DO
  LL=MIN(LL+1,NLVL)

  IF(DIAG) WRITE(*,*)'Input data required numb levels: ',LL,ZLEVEL
  IF(DIAG) WRITE(*,*)'To read all data to height: ',ZDATA

!-------------------------------------------------------------------
! Find which grid definition gives at least LL levels at ZLEVEL
! Z = AA*K^2 + BB*K + CC

  I=0
  KK=0
! revised test 10/9/2007 (rrd) 
  DO WHILE (KK.LE.LL.AND.I.LT.NINTG)
     I=I+1
     DIST=(BZ(I)*BZ(I)-4.0*AZ(I)*(CZ(I)-ZDATA))
     IF(DIST.GE.0.0)THEN
        KK=INT((-BZ(I)+SQRT(DIST))/(2.0*AZ(I)))+1
     ELSE
        KK=1    
     END IF
     IF(DIAG) WRITE(*,*)'Grid loop (kk,ll,i): ',KK,LL,I
  END DO

  AA=AZ(I)
  BB=BZ(I)
  CC=CZ(I)
  ZLEVEL=AA*KK*KK+BB*KK+CC
  NLVL=KK

  IF(DIAG) WRITE(*,*)'Coefficients (aa,bb,cc): ',AA,BB,CC
  IF(DIAG) WRITE(*,*)'Internal grid number of levels: ',KK
  IF(DIAG) WRITE(*,*)'Valid to interpolation height : ',ZLEVEL

!-------------------------------------------------------------------
! Define the index that represents the top of the sfc layer
! always defined to be at 75 m

  DIST=(BB*BB-4.0*AA*(CC-75.0))
  IF(DIST.GE.0.0)THEN
     KSFC=INT((-BB+SQRT(DIST))/(2.0*AA))
  ELSE
     KSFC=1    
  END IF
  SFCL=AA*KSFC*KSFC+BB*KSFC+CC
  IF(DIAG) WRITE(*,*)'Index and top height of sfc layer: ',KSFC, SFCL

END SUBROUTINE metlvl
