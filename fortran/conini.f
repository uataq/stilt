!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONINI           CONcentration array INItialization routine
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CONCENTRATION ARRAY INITIALIZATION ROUTINE IS A SINGLE CALL
!   TO OPEN CONCENTRATION OUTPUT FILES AND WRITE THE INITIAL
!   INDEX RECORD DATA TO THE CONCENTRATION OUTPUT FILES.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Jun 1997 (RRD) 
!                 21 Sep 2000 (RRD) - fortran90 upgrade 
!                 17 Nov 2000 (RRD) - concentration packing
!                 29 Aug 2001 (RRD) - simultaneous multiple meteo files
!                 20 Nov 2003 (RRD) - add minutes field to header
!                 01 Apr 2004 (RRD) - generic file unit definitions
!                 22 Dec 2004 (RRD) - non regular concentration grid
!                 03 Jun 2008 (RRD) - embedded blanks in dir/file
!
! USAGE:  CALL CONINI(SPOT,CONC,DIRT,NLOC,NUMGRD,NUMTYP,CPACK)
!
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           none
!   OUTPUT FILES:          see module funits  
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE CONINI(spot,conc,dirt,nloc,numgrd,numtyp,cpack)

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC' ! pollutant and concentration grid
  INCLUDE 'DEFSPOT.INC' ! multiple source structure

!-------------------------------------------------------------------------------
! argument list variable definitions
!-------------------------------------------------------------------------------

  TYPE(rset), INTENT(IN)    :: spot (:)  ! source location characteristics
  TYPE(cset), INTENT(INOUT) :: conc (:)  ! for each concentration grid 
  TYPE(pset), INTENT(IN)    :: dirt (:)  ! for each pollutant type 
  INTEGER,    INTENT(IN)    :: nloc      ! number of sources
  INTEGER,    INTENT(IN)    :: numgrd    ! number of concentration grids
  INTEGER,    INTENT(IN)    :: numtyp    ! number of pollutants
  INTEGER,    INTENT(IN)    :: cpack     ! concentration packing flag

!-------------------------------------------------------------------------------
! internal variable definitions
!-------------------------------------------------------------------------------

  CHARACTER(80)            :: label      ! dummy character string
  INTEGER                  :: n,ng,kk,kg ! loop indicies
  INTEGER                  :: nzp        ! number of levels
  INTEGER                  :: klen       ! string length  
  INTEGER                  :: kunit      ! output file unit number

!-------------------------------------------------------------------------------
! external variables or procedures  
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! special cpack=2 one grid station format
!-------------------------------------------------------------------------------

  NG=NUMGRD
  IF(CPACK.EQ.2)NG=1

!-------------------------------------------------------------------------------
! go through each grid
!-------------------------------------------------------------------------------

  gloop : DO KG=1,NG

     NZP=CONC(KG)%LEVELS               ! set loop indicies

!    open the files for output
     KUNIT=KF11+KG-1
     IF(KUNIT.GT.KF12)THEN
        WRITE(KF21,*)'*ERROR* conini: exceeding max numb of conc grids'
        WRITE(KF21,*)'Compiled maximum: ',(KF12-KF11+1)
        WRITE(*,*)   '*ERROR* conini: see message file for more information'
        STOP 900
     END IF

     CONC(KG)%UNIT=KUNIT
     LABEL=ADJUSTL(CONC(KG)%DIR)
     KLEN=LEN_TRIM(LABEL)

! for compilers that don't support the bigendian flag (Absoft)
!    OPEN(KUNIT,FILE=LABEL(1:KLEN)//CONC(KG)%FILE,FORM='UNFORMATTED',    &
!         ACCESS='SEQUENTIAL',CONVERT='BIG_ENDIAN')
!***************************************************
!dwen(20090316):adopted from conini of STILT
! JCL:  open conc output file as FORMATTED to generate text file
     OPEN(KUNIT,FILE=LABEL(1:KLEN)//CONC(KG)%FILE,FORM='UNFORMATTED',    &
          ACCESS='SEQUENTIAL')
!     OPEN(KUNIT,FILE=LABEL(1:KLEN)//CONC(KG)%FILE,FORM='FORMATTED',    &
!          ACCESS='SEQUENTIAL')
!***************************************************

!    rec#1  meteo file information for only the starting meteo file 
!dwen(20090826)     WRITE(KUNIT,*) GRID(1,1)%MODEL_ID,FILE(1,1)%FIRST%YR,                 &
     WRITE(KUNIT) GRID(1,1)%MODEL_ID,FILE(1,1)%FIRST%YR,                 &
                  FILE(1,1)%FIRST%MO,FILE(1,1)%FIRST%DA,                 &
                  FILE(1,1)%FIRST%HR,FILE(1,1)%FIRST%IC,                 &
                  NLOC, CPACK

!    rec#2->nloc  source information
!    emission start currently same for all sources
!    added minutes field at the end of record (20 Nov 2003)
     DO N=1,NLOC
!dwen(20090826)        WRITE(KUNIT,*)                                                      &
        WRITE(KUNIT)                                                      &
           DIRT(1)%START%YR,  DIRT(1)%START%MO,                           &
           DIRT(1)%START%DA,  DIRT(1)%START%HR,                           &
           SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL,                        &
           DIRT(1)%START%MN
     END DO

     IF(CPACK.NE.2)THEN
!       horizontal grid index information
!dwen(20090826)        WRITE(KUNIT,*)                                                       &
        WRITE(KUNIT)                                                       &
        CONC(KG)%NUMB_LAT,CONC(KG)%NUMB_LON,CONC(KG)%DELT_LAT,             &
        CONC(KG)%DELT_LON,CONC(KG)%X1Y1_LAT,CONC(KG)%X1Y1_LON
     ELSE
!       each grid has one point, special output to one file
!dwen(20090826)        WRITE(KUNIT,*)                                                       &
        WRITE(KUNIT)                                                       &
        CONC(KG)%NUMB_LAT,CONC(KG)%NUMB_LON,CONC(KG)%DELT_LAT,             &
        CONC(KG)%DELT_LON,  0.0,0.0    
     END IF

!    vertical grid information
!dwen(20090826)     WRITE(KUNIT,*) NZP,(CONC(KG)%HEIGHT(KK),KK=1,NZP)
     WRITE(KUNIT) NZP,(CONC(KG)%HEIGHT(KK),KK=1,NZP)

!    pollutant identification
!dwen(20090826)     WRITE(KUNIT,*) NUMTYP, (DIRT(KK)%IDENT,KK=1,NUMTYP)
     WRITE(KUNIT) NUMTYP, (DIRT(KK)%IDENT,KK=1,NUMTYP)

  END DO gloop

END SUBROUTINE conini
