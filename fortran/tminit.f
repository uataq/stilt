!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  TMINIT           TiMeINITialize   sets up table for tm2 routines
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-02-02
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SETS UP THE ACCUMULATED DAYS TABLE (1900-2100) THAT IS REQUIRED FOR ALL
!   THE OTHER TM2--- CONVERSION PROGRAMS.                   
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 02 Feb 2001 (RRD) - initial version 
!
! USAGE:  CALL TMINIT
!
!   INPUT ARGUMENT LIST:    none
!   OUTPUT ARGUMENT LIST:   none  
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE TMINIT

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER  :: K, IYEAR, NDPY, NDTOT
  INTEGER  :: NADPY (200)     ! accumulated days table from Jan 1, 1900        
                              ! where the first value is for Jan 1, 1901
  COMMON /TMVALS/ NADPY

!-------------------------------------------------------------------------------

  NDTOT=0
  DO IYEAR=1900,2099

     K=IYEAR-1900+1
     NDPY=365

     IF(MOD(IYEAR,100).NE.0)THEN
!       regular year
        IF(MOD(IYEAR,4).EQ.0)   NDPY=366 ! leap year
     ELSE
!       century
        IF(MOD(IYEAR,400).EQ.0) NDPY=366 ! leap century
     END IF

     NDTOT=NDTOT+NDPY
     NADPY(K)=NDTOT 

  END DO

END SUBROUTINE tminit
