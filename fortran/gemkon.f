!#########################################################################
! GEMKON - Module with all the global model constant values 
!-------------------------------------------------------------------------
! Last revised: 21 May 2008 (RRD) - initial version
!-------------------------------------------------------------------------

MODULE gemkon

  INTEGER*4         :: KUMET = 10            ! default meteorology unit

                                             ! PHYSICAL CONSTANTS
  REAL*4, PARAMETER :: ROW     = 1.17        ! STP air density (g/kg)
  REAL*4, PARAMETER :: REARTH  = 6378140.0   ! radius of the earth (m)
  REAL*4, PARAMETER :: PI      = 3.141589    ! circumference over diameter
  REAL*4, PARAMETER :: RADPDEG = PI/180.0    ! radians per degree
  REAL*4, PARAMETER :: GRAV    = 9.80616     ! acceleration of gravity (m/s2)
  REAL*4, PARAMETER :: VONK    = 0.4         ! von Karman's constant
  REAL*4, PARAMETER :: P2JM    = 100.0       ! convert mb to J/m3 
  REAL*4, PARAMETER :: RDRY    = 287.04      ! dry air constant (J/Kg-K)
  REAL*4, PARAMETER :: RGAS    = 0.082057    ! gas constant (atm-liter / deg-mole)
  REAL*4, PARAMETER :: DMVC    = 1.789E-02   ! dynamic viscosity (g m-1 s-1)
  REAL*4, PARAMETER :: FREP    = 6.53E-08    ! mean free path (m at stp)
  REAL*4, PARAMETER :: prn     = 0.923       ! neutral Prandtl number

  REAL*4, PARAMETER :: PSFC    = 1013.0      ! default surface pressure
  REAL*4, PARAMETER :: PTOP    =   10.0      ! top of the model pressure
  REAL*4, PARAMETER :: VWLIM   =    0.02     ! vertical velocity limits (m/s)
  REAL*4, PARAMETER :: TWGHT   =    0.40     ! central point terrain smoothing
  REAL*8, PARAMETER :: DZMIN   =  100.00     ! minimum height difference

!                              coefficients for Beljaars-Holtslag data (1991)
  REAL*4, PARAMETER :: a       =  1.0
  REAL*4, PARAMETER :: b       =  0.66667
  REAL*4, PARAMETER :: c       =  5.0
  REAL*4, PARAMETER :: d       =  0.35

!                              coefficients for phi Ri => z/L at 75m over land
  REAL*4, PARAMETER :: b1      =  0.005
  REAL*4, PARAMETER :: b2      =  41.2
  REAL*4, PARAMETER :: b3      =  1.185
  REAL*4, PARAMETER :: b4      =  1.5
  REAL*4, PARAMETER :: b5      =  1.37
  REAL*4, PARAMETER :: b6      =  0.50

!                              coefficients for psi integral
  REAL*4, PARAMETER :: p1      =  0.1164E-04
  REAL*4, PARAMETER :: p2      = -2.7188
  REAL*4, PARAMETER :: p3      = -2.1551
  REAL*4, PARAMETER :: p4      = -0.9859
  REAL*4, PARAMETER :: p5      = -0.1765

! standard atmosphere heights at 25 hPa intervals from 1000 hpa to 500 hPa
!                             at 50 hPa intervals from  500 hPa to 100 hPa
!                             at 10 hPa intervals from  100 hPa to  10 hPa

  REAL*4, PARAMETER :: STDATM(38)  = (/                                  & 
                       111.,323.,540.,762.,988.,1220.,1457.,1700.,       &
                       1949.,2204.,2466.,2735.,3012.,3297.,3591.,3894.,  &
                       4206.,4530.,4865.,5213.,5574.,6344.,7185.,8117.,  &
                       9164.,10363.,11784.,13608.,16180.,16848.,17595.,  &
                       18442.,19419.,20576.,22000.,23849.,26481.,31055. /)

  REAL*4, PARAMETER :: STDPPP(38)  = (/                                  & 
                       1000.,975.,950.,925.,900.,875.,850.,825.,800.,    &
                       775.,750.,725.,700.,675.,650.,625.,600.,575.,     &
                       550.,525.,500.,450.,400.,350.,300.,250.,200.,     &
                       150.,100.,90.,80.,70.,60.,50.,40.,30.,20.,10.    /)

END MODULE gemkon
