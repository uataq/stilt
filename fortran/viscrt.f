!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:    VISCRT
!   PRGMMR: STUNDER          ORG: R/ARL      DATE: 99-11-15
!
! ABSTRACT: This program written at the Air Resources Laboratory ...
!   is from subroutine c2vis in the VAFTAD vfcst program.  This program
!   defines the ash concentration at the visual ash boundary, which
!   is dependent on the height of the eruption column above the
!   summit and the ash reduction level.  This is used by volcplot and
!   visthrsh.
!
! PROGRAM HISTORY LOG:
!   97-??-??  HEFFTER
!   99-01-07  HEFFTER - added the Cloud Excessive Areal Coverage (CEAC)
!                         adjustment (ash reduction)
!   00-09-25  STUNDER - Output ASCII file to be used with GRIB output file
!                         for runs using the polar stereographic AVN grids.
!   00-11-08  STUNDER - Move ASCII file (VISUAL_ASH_TABLE) creation to
!                         new subroutine vistab.f
!   01-11-30  STUNDER - Return SYM instead of MVIS for hysplit.
!
! USAGE:    CALL VISCRT(OLVL1,OLVL2,IREDUC,MVIS)
!   INPUT ARGUMENT LIST:  see below
!   OUTPUT ARGUMENT LIST: see below
!
! REMARKS: 
!     ITAB     -  ash concentration cutoff for visual ash cloud as
!                   a function of summit and ash column heights for
!                   the default run (no ash reduction)
!     MVIS     - ash concentration cutoff for visual ash (exponent)
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM SP
!
!$$$


  SUBROUTINE VISCRT(OLVL1,OLVL2,IREDUC,SYM)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,    INTENT(IN)  :: olvl1,olvl2 ! release height range
  INTEGER, INTENT(IN)  :: ireduc      ! reduction factor (0=none,1,2,3=most)
  REAL,    INTENT(OUT) :: sym         ! plotting symbol threshold value

!-------------------------------------------------------------------------------

  INTEGER :: ITAB(10,20)		! ivsh,iact
  INTEGER :: i,j,mvis
  INTEGER :: ivsh			! volcano summit height (feet)
  INTEGER :: iact			! ash column top (feet)

!-------------------------------------------------------------------------------

!**critical concentrations

!  This table has summit height going across from 00 to 18 kft (kilo-feet)
!  and ash column height in kft from 00 (top) to 38 (bottom), both
!  in increments of 2 kft.  For example, for a large eruption with ash 
!  column top 34,000 ft or more, ash is 'visual' if the concentration is
!  greater than or equal to 1.0E-18.  For an eruption with an ash column
!  top 20,000 ft (11th row), if the summit is 16,000 ft, the threshold
!  concentration is 1.0E-15, but if the summit is 2000 ft, the threshold
!  is 1.0E-16.   

!                     summit height in kft
!                 00 02 04 06 08 10 12 14 16 18

      data itab/  15,00,00,00,00,00,00,00,00,00,                               &
                  15,15,00,00,00,00,00,00,00,00,                               &
                  15,15,15,00,00,00,00,00,00,00,                               &
                  15,15,15,15,00,00,00,00,00,00,                               &
                  15,15,15,15,15,00,00,00,00,00,                               &
                  16,15,15,15,15,15,00,00,00,00,                               &
                  16,16,15,15,15,15,15,00,00,00,                               &
                  16,16,16,16,15,15,15,15,00,00,                               &
                  16,16,16,16,16,15,15,15,15,00,                               &
                  16,16,16,16,16,16,15,15,15,15,                               &
                  17,16,16,16,16,16,16,15,15,15,                               &
                  17,17,17,16,16,16,16,16,16,15,                               &
                  17,17,17,17,17,16,16,16,16,16,                               &
                  17,17,17,17,17,17,17,16,16,16,                               &
                  17,17,17,17,17,17,17,17,17,16,                               &
                  18,18,17,17,17,17,17,17,17,17,                               &
                  18,18,18,18,18,18,17,17,17,17,                               &
                  18,18,18,18,18,18,18,18,18,18,                               &
                  18,18,18,18,18,18,18,18,18,18,                               &
                  18,18,18,18,18,18,18,18,18,18/         

!-------------------------------------------------------------------------------

      IVSH=OLVL1*3.2808+0.5    	  	  ! convert meters to feet
      IACT=OLVL2*3.2808+0.5 	          ! convert meters to feet

      if(iact.ge.50000) then
        mvis=-19+ireduc
      elseif(iact.ge.40000.and.iact.lt.50000) then
        mvis=-18+ireduc
      else
        i=min(ivsh,20000)/2000+1   ! summit height
        j=iact/2000+1              ! column height
        mvis=-itab(i,j)+ireduc
      end if

      sym=10.0**MVIS			! 1.0E(MVIS)

      return
      end
