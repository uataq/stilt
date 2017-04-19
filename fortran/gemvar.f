!#########################################################################
! GEMVAR - Module that defines the meteorological data arrays
! that are dependent upon the dimensions of the input data
!-------------------------------------------------------------------------
! Last revised: 30 May 2008 (RRD) - initial version
!-------------------------------------------------------------------------

MODULE gemvar

! NGP - HORIZONTAL POINTS       number
! PPP - PRESSURE                mb
! UUU - U WIND COMPONENTS       m/s
! VVV - V WIND COMPONENTS       m/s
! TTT - TEMPERATURE             Kelvin
! MMM - MOISTURE (RH)           percent
! HHH - HEIGHT MSL OF FIELD     meters
! WWW - VERTICAL VELOCITY       mb/s to m/s
! KKK - VERTICAL MIXING         m2/s
! RRR - LOCAL AIR DENSITY       Kg/m3
! QQQ - MIXING RATIO            g/Kg
! DDD - DEPOSITION TOTALS       g/m2
! ZZZ - MIXING RATIO PROFILE    g/Kg
! GSX - WE GRID DIMENSIONS      meters
! GSY - SN GRID DIMENSIONS      meters
! GXY - AREA OF GRID CELL       m2
! SFC - SURFACE PRESSURE        mb
! LUS - LANDUSE CATEGORY        index
! RGH - ROUGHNESS LENGTH        m  
! FRV - FRICTION VELOCITY       m/s 
! PSI - INTEGRATED STABILITY   
! SWF - SHORT WAVE FLIX         w/m2 
! U10 - 10 meter wind           m/s
! V10 - 10 meter wind           m/s
! T02 - 2  meter temperature    deg-K   
! P6H - precipitation rate      m/min

! meteorological file variables
  INTEGER*4, ALLOCATABLE :: NVAR(:)            ! number variables per level
  INTEGER*4, ALLOCATABLE :: KNDX(:)            ! standard atm index for input

! one-dimensional vertical grid variables
  REAL*4,    ALLOCATABLE :: SIG(:)             ! pressure-sigma coordinate
  REAL*4,    ALLOCATABLE :: PPP(:),DPP(:)      ! pressure and delta-pressure
  REAL*4,    ALLOCATABLE :: TCF(:),POT(:)      ! potential temp conversion and temp

! two-dimensional horizontal grid variables
  REAL*8,    ALLOCATABLE :: SFC(:,:), AVG(:,:) ! surface pressure, average pressure
  REAL*8,    ALLOCATABLE :: GSX(:,:), GSY(:,:) ! grid spacing
  REAL*8,    ALLOCATABLE :: GXY(:,:)           ! grid cell area

! other variables
  INTEGER*4, ALLOCATABLE :: NGP(:)             ! number of X grid points by latitude
  REAL*4,    ALLOCATABLE :: VB4(:)             ! temporary profile holding variable
  REAL*8,    ALLOCATABLE :: VB8(:)             ! temporary profile holding variable

! two-dimensional surface and diagnostic   
  INTEGER*4, ALLOCATABLE :: LUS(:,:)           ! landuse category  
  REAL*4,    ALLOCATABLE :: RGH(:,:)           ! aerodynamic roughness length 
  REAL*4,    ALLOCATABLE :: FRV(:,:)           ! friction velocity               
  REAL*4,    ALLOCATABLE :: PSI(:,:)           ! integrated stability function     
  REAL*4,    ALLOCATABLE :: SWF(:,:)           ! incident short-wave flux          

! two-dimensional meteorological variables
  REAL*4,    ALLOCATABLE :: U10(:,:), V10(:,:) ! 10 meter wind components
  REAL*4,    ALLOCATABLE :: T02(:,:)           ! 2 m temperature
  REAL*4,    ALLOCATABLE :: P6H(:,:)           ! 6 hour precipitation 

! three-dimensional meteorological variables
  REAL*8,    ALLOCATABLE :: VAL(:,:)           ! temporary variables
  REAL*8,    ALLOCATABLE :: UUU(:,:,:), VVV(:,:,:), TTT(:,:,:)
  REAL*8,    ALLOCATABLE :: MMM(:,:,:), HHH(:,:,:), WWW(:,:,:)
  REAL*8,    ALLOCATABLE :: KKK(:,:,:), RRR(:,:,:)

! three-dimensional pollutant mass   (x,y,z,p)  
  REAL*8, ALLOCATABLE, TARGET  :: XXX(:,:,:,:) ! master array with multiple species
  REAL*8,              POINTER :: QQQ(:,:,:)   ! computational mixing ratio array
  REAL*8, ALLOCATABLE          :: CCC(:,:,:)   ! temporary mixing ratio array
  REAL*8, ALLOCATABLE          :: ZZZ(    :)   ! temporary mixing ratio profile
  REAL*8, ALLOCATABLE          :: DDD(:,:,  :) ! deposition totals 

! two-dimensional output concentration
  REAL*4, ALLOCATABLE          :: CONC(:,:)    ! mixing ratio converted to concentration

  SAVE

END MODULE gemvar 
