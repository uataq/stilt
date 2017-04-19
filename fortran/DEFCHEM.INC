!-------------------------------------------------------------------------------
! Configure subroutine interface argument lists for main program HYMODELC
! chemistry applications 
! Last Revised: 15 May 2001 (RRD) - initial version
!-------------------------------------------------------------------------------

INTERFACE

!-------------------------------------------------------------------------------
  SUBROUTINE CHEM02(DIRT,DT,MASS)
  IMPLICIT NONE
  INCLUDE 'DEFCONC.INC'     ! concentration and pollutant structure
  TYPE(pset), INTENT(IN)    :: dirt (:)    ! for each pollutant type 
  REAL,       INTENT(IN)    :: dt          ! time step (min)
  REAL,       INTENT(INOUT) :: mass (:)    ! particle mass array
  END SUBROUTINE chem02
!-------------------------------------------------------------------------------

END INTERFACE
