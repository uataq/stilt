!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPSET           DEPosition parameters SET from input data
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION PARAMETERS SET IS THE DATA ENTRY FOR POLLUTANT SPECIES
!   REQUIRED FOR GRAVITATIONAL SETTLING, DRY DEPOSITION, WET REMOVAL
!   AND RADIOACTIVE DECAY COMPUTATIONS. ONE SET OF ENTIRES REQUIRED
!   EACH DEFINED POLLUTANT TYPE.
! 
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 03 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 27 Aug 2003 (RRD) - volume units conversion 
!                 02 Apr 2004 (RRD) - generic file unit numbers
!
! USAGE:  CALL DEPSET(DIRT,NUMTYP,IUNIT,CDEP,RDEP,SDEP)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             unit 5 or unit KF21 if CONTROL file found
!   OUTPUT FILES:            unit KF22 to file STARTUP if input from unit 5
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE DEPSET(DIRT,NUMTYP,IUNIT,CDEP,RDEP,SDEP)

  USE funits

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'        ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(pset), INTENT(INOUT) :: dirt(:)   ! for each pollutant type 
  INTEGER,    INTENT(IN)    :: numtyp    ! number of pollutant types
  INTEGER,    INTENT(IN)    :: iunit     ! unit number for input data
  LOGICAL,    INTENT(OUT)   :: cdep      ! indicate wet or dry deposition
  LOGICAL,    INTENT(OUT)   :: rdep      ! resistance deposition
  LOGICAL,    INTENT(OUT)   :: sdep      ! resuspension option

!-------------------------------------------------------------------------------
! internally defined variables
!-------------------------------------------------------------------------------

  INTEGER    :: kk, numpol

!-------------------------------------------------------------------------------

  CDEP=.FALSE.
  RDEP=.FALSE.
  SDEP=.FALSE.

!-------------------------------------------------------------------------------
! number of pollutants already defined
!-------------------------------------------------------------------------------

  IF(IUNIT.EQ.5)WRITE(*,*)'Enter number of pollutants:', NUMTYP
  READ(IUNIT,*)NUMPOL
  IF(NUMPOL.NE.NUMTYP)THEN
     WRITE(*,*)'WARNING depset: # pollutants defined -',NUMTYP
     NUMPOL=NUMTYP
  END IF
  IF(IUNIT.EQ.5)WRITE(KF22,'(I2.2)')NUMPOL
 
  DO KK=1,NUMTYP

  IF(IUNIT.EQ.5) WRITE(*,*)'Enter data for Pollutant:',DIRT(KK)%IDENT

!-------------------------------------------------------------------------------
! define as gas or particle
!-------------------------------------------------------------------------------

  DIRT(KK)%DOGAS=.FALSE.
  DIRT(KK)%PDIAM=0.0
  DIRT(KK)%PDENS=0.0
  DIRT(KK)%SHAPE=0.0

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Define as gas (all zero) or particle'
     WRITE(*,*)'Diameter (um), density (g/cc), shape (1-2)'
     WRITE(*,*)DIRT(KK)%PDIAM, DIRT(KK)%PDENS, DIRT(KK)%SHAPE
  END IF

  READ(IUNIT,*)DIRT(KK)%PDIAM, DIRT(KK)%PDENS, DIRT(KK)%SHAPE
  IF(IUNIT.EQ.5) WRITE(KF22,'(3E10.3)')                                   &
     DIRT(KK)%PDIAM, DIRT(KK)%PDENS, DIRT(KK)%SHAPE

! all three options required to define a particle
  IF(DIRT(KK)%PDIAM.EQ.0.0 .OR. DIRT(KK)%PDENS.EQ.0.0 .OR.              &
     DIRT(KK)%SHAPE.EQ.0.0)     DIRT(KK)%DOGAS=.TRUE.

!-------------------------------------------------------------------------------
! define dry deposition velocity (over-rides grav settling)
!-------------------------------------------------------------------------------

  DIRT(KK)%DORES=.FALSE.
  DIRT(KK)%DODRY=.FALSE.
  DIRT(KK)%DOVOL=.FALSE.
  DIRT(KK)%DRYVL=0.0
  DIRT(KK)%GPMOL=0.0
  DIRT(KK)%ACVTY=0.0
  DIRT(KK)%DIFTY=0.0
  DIRT(KK)%HENRY=0.0

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Dry deposition - '
     WRITE(*,*)'---- explicit ----  --- resistence method --->'
     WRITE(*,*)'Set velocity (m/s), Gram Molecular Wt',                 &
               'Activity Ratio, Diffusivity Ratio, Henrys Constant'
     WRITE(*,*)DIRT(KK)%DRYVL, DIRT(KK)%GPMOL, DIRT(KK)%ACVTY,          &
               DIRT(KK)%DIFTY, DIRT(KK)%HENRY
  END IF

  READ(IUNIT,*)DIRT(KK)%DRYVL, DIRT(KK)%GPMOL, DIRT(KK)%ACVTY,          &
               DIRT(KK)%DIFTY, DIRT(KK)%HENRY
  IF(IUNIT.EQ.5) WRITE(KF22,'(5E10.3)')                                 &
               DIRT(KK)%DRYVL, DIRT(KK)%GPMOL, DIRT(KK)%ACVTY,          &
               DIRT(KK)%DIFTY, DIRT(KK)%HENRY
  IF(DIRT(KK)%DRYVL.GT.0.0)DIRT(KK)%DODRY=.TRUE.

! negative molecular weight indicates that emitted mass (kg) will be
! converted to volume units (ppm) only in subroutine emsgrd
  IF(DIRT(KK)%GPMOL.LT.0.0)THEN
     DIRT(KK)%DOVOL=.TRUE.
     DIRT(KK)%GPMOL=ABS(DIRT(KK)%GPMOL)
  END IF

! resistence method requires molecular weight for gases
! and only diameter for particles, however will use molecular
! weight as the turn-on resistence flag for both particles and gases
  IF(DIRT(KK)%GPMOL.GT.0.0)THEN
     DIRT(KK)%DORES=.TRUE.
     IF(DIRT(KK)%DOGAS)THEN   
        IF((DIRT(KK)%DRYVL+DIRT(KK)%ACVTY+DIRT(KK)%DIFTY+               &
            DIRT(KK)%HENRY).EQ.0.0) DIRT(KK)%DORES=.FALSE.
     END IF
  END IF

! check for consistent inputs for gravitational settling
  DIRT(KK)%DOGRV=.FALSE.
  IF(DIRT(KK)%PDIAM.GT.0.0.AND. DIRT(KK)%PDENS.GT.0.0.AND.              &
     DIRT(KK)%SHAPE.GT.0.0)     DIRT(KK)%DOGRV=.TRUE.

! specified dry deposition over-rides gravitational settling
! and resistance method specification
  IF(DIRT(KK)%DRYVL.GT.0.0)THEN
     DIRT(KK)%DOGRV=.FALSE.
     DIRT(KK)%DORES=.FALSE.
  END IF

!-------------------------------------------------------------------------------
! wet removal constants
!-------------------------------------------------------------------------------

  DIRT(KK)%DOWET=.FALSE.
  DIRT(KK)%WETGAS=0.0
  DIRT(KK)%WETIN=0.0
  DIRT(KK)%WETLO=0.0

  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Wet removal constants - '
     WRITE(*,*)'-- gasses ---   ----------- particles -----------'
     WRITE(*,*)'Henrys (M/atm), In-cloud (L/L), Below-cloud (1/s)'
     WRITE(*,*)DIRT(KK)%WETGAS, DIRT(KK)%WETIN, DIRT(KK)%WETLO
  END IF

  READ(IUNIT,*)DIRT(KK)%WETGAS, DIRT(KK)%WETIN, DIRT(KK)%WETLO
  IF(IUNIT.EQ.5) WRITE(KF22,'(3E10.3)')                                  &
     DIRT(KK)%WETGAS, DIRT(KK)%WETIN, DIRT(KK)%WETLO

! check for consistency of wet removal definitions
  IF(DIRT(KK)%DOGAS)THEN
!    wet removal of gasses only if they are soluable
     DIRT(KK)%WETIN=0.0
     DIRT(KK)%WETLO=0.0
  ELSE
     DIRT(KK)%WETGAS=0.0
  END IF

  IF(DIRT(KK)%WETGAS.GT.0.0 .OR. DIRT(KK)%WETIN.GT.0.0 .OR.             &
     DIRT(KK)%WETLO .GT.0.0)     DIRT(KK)%DOWET=.TRUE.

!-------------------------------------------------------------------------------
! radioactive decay options (one year default)
!-------------------------------------------------------------------------------

  DIRT(KK)%DORAD=.FALSE.
  DIRT(KK)%RHALF=0.0

  IF(IUNIT.EQ.5)THEN
      WRITE(*,*)'Radioactive decay half-life (days / 0=none)'
      WRITE(*,*)DIRT(KK)%RHALF
  END IF

  READ(IUNIT,*)DIRT(KK)%RHALF
  IF(IUNIT.EQ.5)WRITE(KF22,'(E10.3)')DIRT(KK)%RHALF

  IF(DIRT(KK)%RHALF.GT.0.0)DIRT(KK)%DORAD=.TRUE.

!-------------------------------------------------------------------------------
! deposition resuspension option
!-------------------------------------------------------------------------------

  DIRT(KK)%DOSUS=.FALSE.
  DIRT(KK)%SRATE=0.0
  IF(IUNIT.EQ.5)THEN
     WRITE(*,*)'Deposition resuspension constant (1/m)'
     WRITE(*,*)DIRT(KK)%SRATE
  END IF

  READ(IUNIT,*)DIRT(KK)%SRATE
  IF(IUNIT.EQ.5)WRITE(KF22,'(E10.3)')DIRT(KK)%SRATE
  IF(DIRT(KK)%SRATE.GT.0.0)DIRT(KK)%DOSUS=.TRUE.

!-------------------------------------------------------------------------------
! set species independent flags to see if any deposition options enabled
!-------------------------------------------------------------------------------

  IF(DIRT(KK)%DODRY.OR.DIRT(KK)%DORES.OR.                               &
     DIRT(KK)%DOGRV.OR.DIRT(KK)%DOWET.OR.DIRT(KK)%DORAD)                &
     CDEP=.TRUE.

  IF(DIRT(KK)%DORES)RDEP=.TRUE.
  IF(DIRT(KK)%DOSUS)SDEP=.TRUE.

  END DO 

END SUBROUTINE depset
