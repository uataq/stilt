!-------------------------------------------------------------------------------
! DEFSPOT.INC 
! Defines structure of starting location and time
! Last Revised: 24 Jun 1996
!               22 Apr 1999 (RRD) - added source term information
!               20 Sep 2000 (RRD) - fortran90 upgrade
!               04 Oct 2005 (RRD) - boundary layer fractional height
!               09 Jun 2008 (RRD) - added minutes
!-------------------------------------------------------------------------------

  TYPE RSET
     SEQUENCE

     REAL        :: olat,olon            ! starting location
     REAL        :: olvl                 ! starting height (agl)
     REAL        :: zlvl                 ! starting height (Zi fraction)
     REAL        :: qtrm                 ! emission rate (units/hr)
     REAL        :: area                 ! emission area (m^2)
     INTEGER     :: ibyr,ibmo,ibda       ! starting date 
     INTEGER     :: ibhr,ibmn            ! starting time
     REAL        :: xp,yp,zp             ! current x,y,x position
     REAL        :: zv                   ! previous height for line sources
     INTEGER     :: kg                   ! calculation meteo grid number

  END TYPE
