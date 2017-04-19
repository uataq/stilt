!-------------------------------------------------------------------------------
! Module for stability analysis constants   
! Primarily referenced in stb??? subroutines 
! Last Revised: 02 Dec 2003 (RRD) - initial version
!               08 Aug 2008 (RRD) - moved variables to common block
!-------------------------------------------------------------------------------

MODULE stbcon

!                             universal constants
  REAL,      PARAMETER     :: grav   = 9.80616      ! gravity (m/s2)
  REAL,      PARAMETER     :: vonk   = 0.4          ! Karman's
  REAL,      PARAMETER     :: prn    = 0.923        ! neutral Prandtl number
  REAL,      PARAMETER     :: cp     = 1005.0       ! heat (J/Kg-K)

!                             arbitrary diffusivity limits
  REAL,      PARAMETER     :: vmin   = 0.01         ! Mixing minimum
  REAL,      PARAMETER     :: vmax   = 200.0        ! Mixing maximum

!                             coefficients for Beljaars-Holtslag data (1991)
  REAL,      PARAMETER     :: a      =  1.0
  REAL,      PARAMETER     :: b      =  0.66667
  REAL,      PARAMETER     :: c      =  5.0
  REAL,      PARAMETER     :: d      =  0.35

!                             coefficients for mixing length scales
  REAL,      PARAMETER     :: a1     =  0.2828E-03
  REAL,      PARAMETER     :: a2     =  0.8049
  REAL,      PARAMETER     :: a3     =  1.6583
  REAL,      PARAMETER     :: a4     =  0.5090E-02
  REAL,      PARAMETER     :: a5     = -1.0063E-03

!                             coefficients for phi Ri => z/L at 75m over land
  REAL,      PARAMETER     :: b1     =  0.005
  REAL,      PARAMETER     :: b2     =  41.2
  REAL,      PARAMETER     :: b3     =  1.185
  REAL,      PARAMETER     :: b4     =  1.5
  REAL,      PARAMETER     :: b5     =  1.37
  REAL,      PARAMETER     :: b6     =  0.50

!                             coefficients for psi integral
  REAL,      PARAMETER     :: p1     =  0.1164E-04
  REAL,      PARAMETER     :: p2     = -2.7188
  REAL,      PARAMETER     :: p3     = -2.1551
  REAL,      PARAMETER     :: p4     = -0.9859
  REAL,      PARAMETER     :: p5     = -0.1765

END MODULE stbcon
