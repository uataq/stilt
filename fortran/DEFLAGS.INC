!-------------------------------------------------------------------------------
! DEFLAGS.INC 
! DEFines structure of the LAGrangian Sampling option
! Last Revised: 12 Oct 2005 (RRD) - initial version
!-------------------------------------------------------------------------------

  TYPE LSET
     SEQUENCE

     REAL          :: rlat,rlon               ! release location
     REAL          :: rlvl                    ! release height (agl)
     REAL          :: wdir,wspd               ! force trajectory vector
     REAL          :: xp,yp,zp                ! grid position
     INTEGER       :: ryr,rmo,rda,rhr,rmn     ! release time
     INTEGER       :: syr,smo,sda,shr,smn     ! sample start time
     INTEGER       :: racm,sacm               ! accumulated minutes
     INTEGER       :: kg                      ! meteorological grid
     INTEGER       :: avrg                    ! sample averaging time (min)
     INTEGER       :: disk                    ! sample disk output (min)
     CHARACTER(80) :: file                    ! output file name

  END TYPE
