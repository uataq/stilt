!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DEPRAD           DEPosition RADioactive decay
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEPOSITION RADIOACTIVE COMPUTES THE RADIOACTIVE DECAY OF
!   POLLUTANTS THAT HAVE ALREADY BEEN DEPOSITED.  THE DECAY IS
!   APPLIED TO THE VALUES THAT ARE BEING SUMMED IN THE DEPOSITION
!   ARRAY EACH TIME STEP.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 21 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 05 Aug 2003 (RRD) - converted to implicit loop
!                 17 Oct 2007 (RRD) - forced positive time step
!
! USAGE:  CALL DEPRAD(CONC,DIRT,NUMGRD,NUMTYP,DT,CSUM)
! 
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE DEPRAD(CONC,DIRT,NUMGRD,NUMTYP,DT,CSUM)

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC'         ! pollutant and concentration grid

!-------------------------------------------------------------------------------
! argument list definitions
!-------------------------------------------------------------------------------

  TYPE(cset),   INTENT(IN)    :: conc(:)   ! for each concentration grid 
  TYPE(pset),   INTENT(IN)    :: dirt(:)   ! for each pollutant type 
  INTEGER,      INTENT(IN)    :: numgrd    ! number of concentration grids
  INTEGER,      INTENT(IN)    :: numtyp    ! number of pollutants
  REAL,         INTENT(IN)    :: dt        ! intgration time step (min)
                                           ! conc array (x,y,z,grids,species)
  REAL,         INTENT(INOUT) :: csum (:,:,:,:,:) 

!-------------------------------------------------------------------------------
! internal variable definitinons
!-------------------------------------------------------------------------------

  INTEGER         :: ii,jj,kl,kt,kg
  REAL            :: rfr,rtc

!-------------------------------------------------------------------------------

! go through each grid
  DO KG=1,NUMGRD

     DO KT=1,NUMTYP
        IF(DIRT(KT)%DORAD)THEN

!          convert half-life to time constant
           RTC=-LOG(0.5)/(DIRT(KT)%RHALF*86400.0)
!          convert time constant to removal fraction
           RFR=EXP(-60.0*ABS(DT)*RTC)

           DO KL=1,CONC(KG)%LEVELS
!             zero height is for deposition
              IF(CONC(KG)%HEIGHT(KL).EQ.0.0)             &
                 CSUM(:,:,KL,KT,KG)=CSUM(:,:,KL,KT,KG)*RFR
           END DO

        END IF
     END DO

  END DO

END SUBROUTINE deprad
