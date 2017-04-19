!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  RUNSET           RUN SETup for basic model parameters
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   RUN SETUP SETS THE MOST BASIC MODEL SIMULATION PARAMETERS, BEFORE
!   ANY OTHER INFORMATION IS KNOWN. IF THERE EXISTS A FILE NAMED CONTROL
!   THEN ALL INPUT PARAMETERS WILL BE READ FROM THAT FILE.  OTHERWISE
!   INPUT IS EXPECTED ON STANDARD INPUT (UNIT 5).  IN THAT CASE AN OUTPUT
!   FILE CALLED STARTUP WILL BE CREATED THAT WILL CONTAIN ALL DATA ENTRIES
!   IT MAY BE USED IN SUBSEQUENT SIMULATIONS AS THE CONTROL FILE.
!   ENTRIES FROM STANDARD INPUT HAVE DEFAULT VALUES WHICH MAY BE SELECTED
!   BY JUST ENTERING "/".  THIS IS NOT VALID FOR CHARACTER STRINGS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 21 Jul 1998 (RRD)
!                 09 Mar 1999 (RRD) - additonal subgrid parameters added
!                 08 Apr 1999 (RRD) - location specific emissions
!                 22 Apr 1999 (RRD) - renamed surface layer depth variable
!                                   - source/area at release location
!                 26 Aug 1999 (RRD) - common block to define vertical grid
!                 09 Dec 1999 (RRD) - reconfigured source definition
!                 01 May 2000 (RRD) - zdata initial value instead of zmdl
!                 17 Nov 2000 (RRD) - fortran90 upgrade
!                 16 Mar 2001 (RRD) - nlvl initialization
!                 26 Sep 2001 (RRD) - ntim meteo input file order defined
!                 05 Oct 2001 (RRD) - search for comment to exclude
!                 11 Feb 2002 (RRD) - formatted STARTUP file
!                 21 May 2002 (RRD) - divergence option
!                 12 Jul 2002 (RRD) - input data decoder
!                 23 Jul 2002 (RRD) - eta terrain correction option
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 24 Sep 2002 (RRD) - test for number of meteo inputs
!                 13 Feb 2003 (RRD) - moved vertical definitions to metlvl 
!                 02 Apr 2004 (RRD) - generic file number definitions
!                 03 Jun 2008 (RRD) - embedded blanks dir/file
!                 09 Jun 2008 (RRD) - minutes to start time definition
!
! USAGE:  CALL RUNSET(SPOT,NLOC,NHRS,NGRD,NTIM,ZMDL,ZDATA,BACK,IUNIT,KVEL)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            units 5 or KF21 depending if CONTROL file is found
!   OUTPUT FILES:           unit KF22 when input defined on unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

SUBROUTINE RUNSET(SPOT,NLOC,NHRS,NGRD,NTIM,ZMDL,ZDATA,BACK,IUNIT,KVEL)

  USE funits
  use module_defgrid ! meteorology file and grid

  IMPLICIT NONE

  INCLUDE 'DEFSPOT.INC' ! source information

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(rset),INTENT(INOUT) :: spot (:)  ! source location characteristics
  INTEGER,   INTENT(IN )   :: nloc      ! number of starting locations
  INTEGER,   INTENT(OUT)   :: nhrs      ! simulation duration in hours
  INTEGER,   INTENT(OUT)   :: ngrd      ! number of data grids for this run
  INTEGER,   INTENT(OUT)   :: ntim      ! number of data time periods this run
  REAL,      INTENT(OUT)   :: zmdl      ! maximum height for scaling coord 
  REAL,      INTENT(OUT)   :: zdata     ! maximum height for input data    
  LOGICAL,   INTENT(OUT)   :: back      ! integration direction flag
  INTEGER,   INTENT(IN )   :: iunit     ! unit for file of input parameters
  INTEGER,   INTENT(OUT)   :: kvel      ! vertical velocity remapping

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  CHARACTER(80)     :: label 
  INTEGER           :: n,k,kv,kk,kg,kt,kl,kret

!-------------------------------------------------------------------------------
! external variables
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE DECODR(IUNIT,VAR1,VAR2,VAR3,VAR4,VAR5)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: IUNIT   ! unit number
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR1
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR2
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR3
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR4
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR5
  END SUBROUTINE decodr
  END INTERFACE

!-------------------------------------------------------------------------------
! generic model defaults

  NHRS=48
  KVEL=0
  ZDATA=10000.0
  NGRD=1
  NTIM=1
  FILE(1,1)%DIR='|main|sub|data|'
  FILE(1,1)%METEO='file_name'

!-------------------------------------------------------------------------------
! starting location information
!-------------------------------------------------------------------------------

   DO N=1,NLOC
      SPOT(N)%AREA=0.0
      SPOT(N)%QTRM=0.0
!     default starting locations
      IF(N.EQ.1)THEN
         SPOT(N)%OLAT=40.0
         SPOT(N)%OLON=-90.0
         SPOT(N)%OLVL=50.0
      ELSE
         SPOT(N)%OLAT=SPOT(N-1)%OLAT
         SPOT(N)%OLON=SPOT(N-1)%OLON
         SPOT(N)%OLVL=SPOT(N-1)%OLVL
      END IF

      IF(IUNIT.EQ.5)THEN
         WRITE(*,*)'Enter starting location (lat, lon, m-agl)'
         WRITE(*,*)SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL
      END IF

!-------------------------------------------------------------------------------
! read standard 3 parameters or variable parameters
!-------------------------------------------------------------------------------

!     standard input configuration
!     READ(IUNIT,*)SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL
      CALL DECODR(IUNIT,SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL,     &
                        SPOT(N)%QTRM,SPOT(N)%AREA)

      IF(IUNIT.EQ.5)WRITE(KF22,'(3F10.3)')SPOT(N)%OLAT,SPOT(N)%OLON,SPOT(N)%OLVL

!     for now only permit one starting time for all locations
      SPOT(N)%IBYR=SPOT(1)%IBYR
      SPOT(N)%IBMO=SPOT(1)%IBMO
      SPOT(N)%IBDA=SPOT(1)%IBDA
      SPOT(N)%IBHR=SPOT(1)%IBHR
      SPOT(N)%IBMN=SPOT(1)%IBMN
   END DO

!-------------------------------------------------------------------------------
! simulation run time, vertical coordinate
!-------------------------------------------------------------------------------

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Enter total run time (hours)'
     WRITE(*,*)NHRS
  END IF
  READ(IUNIT,*)NHRS
  IF(IUNIT.EQ.5)WRITE(KF22,'(I4.4)')NHRS

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Vertical (0:data 1:isob 2:isen 3:dens 4:sigma 5:diverg 6:eta)'
     WRITE(*,*)KVEL
  END IF
  READ(IUNIT,*)KVEL
  IF(IUNIT.EQ.5)WRITE(KF22,'(I1)')KVEL

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Top of model domain (internal coordinates m-msl)'
     WRITE(*,*)ZDATA
  END IF
  READ(IUNIT,*)ZDATA
  IF(IUNIT.EQ.5)WRITE(KF22,'(F10.1)')ZDATA

! ZMDL also is the terrain scaling parameter, the height at which the
! sigma surfaces go flat relative to terrain.

  ZMDL=MAX(ABS(ZDATA),25000.0)

! when zdata<0 then use abs(value) for scaling input data on terrain 
! following coordinate system such as rams (20000m) or coamps (34800m)
! zmdlt variable is only used for terrain coordinate systems

  IF(ZDATA.LT.0.0)THEN
     DREC(1,1)%ZMDLT=ABS(ZDATA)
  ELSE
!    assume default rams value
     DREC(1,1)%ZMDLT=20000.0
!    assume default coamps value
!    DREC(1,1)%ZMDLT=34800.0
  END IF

  DO N=1,NLOC
!    starting height limit
     IF(SPOT(N)%OLVL.GT.ABS(ZDATA))THEN
        WRITE(*,*)'*ERROR* runset: Start height above mdl domain'
        WRITE(*,*)'   Source - ',N,'    Height - ',SPOT(N)%OLVL
        WRITE(*,*)'   Model top ht - ',ABS(ZDATA)
        STOP 900
     END IF
  END DO

!-------------------------------------------------------------------------------
! meteorological grid information - number of grids and time periods
!-------------------------------------------------------------------------------

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Number of input grids'
     WRITE(*,*)NGRD
  END IF

! old style forced read required both variables
! READ(IUNIT,*)NGRD

! variable input includes optional number of time periods
  READ(IUNIT,'(A)')LABEL

! check for #label information
  K=INDEX(LABEL,'#')-1
  IF(K.LE.0)K=80

! extract each field from label
  KK=1
  KV=0
  DO WHILE (KK.LT.K.AND.KV.LT.5)
!    find next non-blank character
     DO WHILE (KK.LT.K.AND.LABEL(KK:KK).EQ.' ')
        KK=KK+1
     END DO

!    find relative end of character field
     KL=INDEX(LABEL(KK:),' ')
     IF(KL.GT.0.AND.KK.LT.K)THEN
        KV=KV+1
        IF(KV.EQ.1)READ(LABEL(KK:),*)NGRD 
        IF(KV.EQ.2)READ(LABEL(KK:),*)NTIM
!       convert end back to absolute units
        KK=KK+KL-1
     END IF
  END DO
  IF(IUNIT.EQ.5)WRITE(KF22,'(I1)')NGRD

! test number of grid limits
  IF(NGRD.GT.MGRD)THEN
     WRITE(*,*)'*ERROR* runset: Numb meteo grids exceed DEFGRID limit - ',MGRD
     STOP 900
  ELSEIF(NGRD.EQ.0)THEN
     WRITE(*,*)'*ERROR* runset: Number of input meteorology grids = zero'
     STOP 900
  ELSE
     CONTINUE
  END IF

! test number of times limits
  IF(NTIM.GT.MTIM)THEN
     WRITE(*,*)'*ERROR* runset: Numb meteo times exceed DEFGRID limit - ',MTIM
     STOP 900
  ELSEIF(NTIM.EQ.0)THEN
     WRITE(*,*)'*ERROR* runset: Number of meteorology time periods = zero'
     STOP 900
  ELSE
     CONTINUE
  END IF

! input file order should be defined as all files in time sequence for grid #1
! followed by all files in time sequence for grid #2, etc.  This structure maintains
! the maximum compatibility with older style control files.

  DO KG=1,NGRD

     IF(KG.GT.1)THEN
        FILE(KG,1)%DIR=FILE(KG-1,1)%DIR
        FILE(KG,1)%METEO=FILE(KG-1,1)%METEO
!       save the scaling height for each grid
        DREC(KG,1)%ZMDLT=DREC(KG-1,1)%ZMDLT
     END IF

  DO KT=1,NTIM

     IF(KT.GT.1)THEN
        FILE(KG,KT)%DIR=FILE(KG,KT-1)%DIR
        FILE(KG,KT)%METEO=FILE(KG,KT-1)%METEO
!       save the scaling height for each grid
        DREC(KG,KT)%ZMDLT=DREC(KG,KT-1)%ZMDLT
     END IF

     IF(IUNIT.EQ.5)THEN
        WRITE(*,'(A,I2,A)')' Enter grid #',KG,' directory (|...|)'
        K=LEN_TRIM(ADJUSTL(FILE(KG,KT)%DIR))
        WRITE(*,'(1X,A)')FILE(KG,KT)%DIR(1:K)
     END IF
     READ(IUNIT,'(A)')FILE(KG,KT)%DIR
     FILE(KG,KT)%DIR=ADJUSTL(FILE(KG,KT)%DIR)
     K=LEN_TRIM(FILE(KG,KT)%DIR)
     IF(IUNIT.EQ.5)WRITE(KF22,'(A)')FILE(KG,KT)%DIR(1:K)

     IF(IUNIT.EQ.5)THEN
        WRITE(*,'(A,I2,A)')' Enter grid #',KG,' file name (?????)'
        K=LEN_TRIM(ADJUSTL(FILE(KG,KT)%METEO))
        WRITE(*,'(1X,A)')FILE(KG,KT)%METEO(1:K)
     END IF
     READ(IUNIT,'(A)')FILE(KG,KT)%METEO
     FILE(KG,KT)%METEO=ADJUSTL(FILE(KG,KT)%METEO)
     K=LEN_TRIM(FILE(KG,KT)%METEO)
     IF(IUNIT.EQ.5)WRITE(KF22,'(A)')FILE(KG,KT)%METEO(1:K)

  END DO
  END DO

! set integration direction
  BACK=.FALSE.
  IF(NHRS.LT.0)THEN
     BACK=.TRUE.
     NHRS=-NHRS
  END IF

END SUBROUTINE runset
