!-------------------------------------------------------------------------------
! DEFSPRT.INC 
! Defines structure of complex point source emissions input from a data file
! Last Revised: 13 Apr 2006 (RRD) - initial version
!-------------------------------------------------------------------------------

  TYPE QSET
     SEQUENCE

     INTEGER     :: start                ! starting time (acc minutes)
     INTEGER     :: stop                 ! stopping time (acc minutes)
     INTEGER     :: kg                   ! initial meteorology grid
     REAL        :: xp,yp                ! (x,y) source position
     REAL        :: qlvl                 ! release height (m agl)  
     REAL        :: rate                 ! emission rate (units/hr)
     REAL        :: area                 ! emission area (m^2)
     REAL        :: heat                 ! heat release rate (watts)

  END TYPE
