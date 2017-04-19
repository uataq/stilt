!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAROUT           PARticle OUTput to file
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:05-08-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   WRITES ALL THE POLLUTANT PARTICLE INFORMATION TO A FILE WHICH
!   CAN BE USED TO INITIALIZE THE MODEL OR POSTPROCESSING GRPAHICS.
!   THE FILE MAY CONTAIN ONE OR MORE TIME PERIODS. 
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 08 May 2003 (RRD) - Initial version
!                 02 Apr 2004 (RRD) - Generic file unit numbers
!                 10 Aug 2004 (RRD) - Argument list change
!
! USAGE:  CALL PAROUT(JET,KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,
!              SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,ZMDL)
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

SUBROUTINE PAROUT(JET,KPM,MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,    &
                  SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,NSORT,ZMDL)

  USE funits
  use module_defgrid

  IMPLICIT NONE


!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)    :: jet        ! current elapsed time
  INTEGER, INTENT(INOUT) :: kpm        ! total number of puffs or particles
  REAL,    INTENT(INOUT) :: mass (:,:) ! mass of pollutant (arbitrary units)
  REAL,    INTENT(INOUT) :: xpos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: ypos (:)   ! puff center positions (grid units)
  REAL,    INTENT(INOUT) :: zpos (:)   ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh (:)   ! horizontal puff sigma 
  REAL,    INTENT(INOUT) :: sigu (:)   ! turbulence u'2    
  REAL,    INTENT(INOUT) :: sigv (:)   ! turbulence v'2
  REAL,    INTENT(INOUT) :: sigw (:)   ! turbulence w'2 or vertical puff sigma
  INTEGER, INTENT(INOUT) :: hdwp (:)   ! Horizontal distribution pollutant
  INTEGER, INTENT(INOUT) :: page (:)   ! pollutant age since release (min)
  INTEGER, INTENT(INOUT) :: ptyp (:)   ! pollutant type index number
  INTEGER, INTENT(INOUT) :: pgrd (:)   ! meteorological grid of puff position
  INTEGER, INTENT(INOUT) :: nsort(:)   ! sort index array by position
  REAL,    INTENT(IN)    :: zmdl       ! vertical index scaling height

!-------------------------------------------------------------------------------

  REAL     :: TLAT,TLON  ! particle position
  REAL     :: ZTMP       ! particle height
  INTEGER  :: NUMPOL     ! number of different pollutants
  INTEGER  :: IYR        ! file year
  INTEGER  :: IMO        ! file month
  INTEGER  :: IDA        ! file day 
  INTEGER  :: IHR        ! file hour 
  INTEGER  :: IMN        ! file minutes
  INTEGER  :: I,J        ! local indicies


!-------------------------------------------------------------------------------

! set pollutant dimension
  NUMPOL=SIZE(mass,1)

  CALL TM2DAY(JET,IYR,IMO,IDA,IHR,IMN)
  WRITE(KF24)KPM,NUMPOL,IYR,IMO,IDA,IHR

  DO J=1,KPM

!    convert positions to lat/lon before dump
     IF(GRID(PGRD(J),1)%LATLON)THEN
        CALL GBL2LL(PGRD(J),1,XPOS(J),YPOS(J),TLAT,TLON)
     ELSE
        CALL CXY2LL_wps(GRID(PGRD(J),1)%GBASE,XPOS(J),YPOS(J),TLAT,TLON,GRID(PGRD(J),1)%proj)
     END IF

!    convert sigma to agl assuming terrain height = 0
     ZTMP=(1.0-ZPOS(J))*ZMDL

     WRITE(KF24)(MASS(I,J),I=1,NUMPOL)
!    additional variables not yet supported in supplemental programs
!    WRITE(KF24)TLAT,TLON,ZTMP,SIGH(J),SIGW(J),SIGV(J),SIGU(J)
     WRITE(KF24)TLAT,TLON,ZTMP,SIGH(J),SIGW(J),SIGV(J)
     WRITE(KF24)PAGE(J),HDWP(J),PTYP(J),PGRD(J),NSORT(J)

  END DO

  WRITE(KF21,*)' NOTICE parout: created PARDUMP file at - ',IYR,IMO,IDA,IHR
  WRITE(KF21,*)' NOTICE parout: number particles output - ',KPM

END SUBROUTINE parout
