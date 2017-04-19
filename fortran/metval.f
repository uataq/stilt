!-------------------------------------------------------------------------------
! Module for meteorological data arrays     
! Primarily referenced in subroutine advpnt and conditional compilation
! chemistry routines.
! Last Revised: 16 Sep 2003 (RRD) - added integrated stability as variable
!               17 Oct 2003 (RRD) - converted density to tubulent kinetic energy
!               05 Nov 2003 (RRD) - mixing variables converted to turbulent vel
!               08 Mar 2006 (RRD) - static stability and velocity scalar
!-------------------------------------------------------------------------------

MODULE metval

  REAL,    ALLOCATABLE :: u  (:,:,:,:,:) ! x wind (inp m/s; out grid/min)
  REAL,    ALLOCATABLE :: v  (:,:,:,:,:) ! v wind (inp m/s; out grid/min)
  REAL,    ALLOCATABLE :: w  (:,:,:,:,:) ! z wind (inp dp/dt; out sigma/min)
  REAL,    ALLOCATABLE :: a  (:,:,:,:,:) ! ambient temp (deg-K)
  REAL,    ALLOCATABLE :: t  (:,:,:,:,:) ! potential temp (deg-K)
  REAL,    ALLOCATABLE :: q  (:,:,:,:,:) ! humid (inp rh or kg/kg; out rh/100)
  REAL,    ALLOCATABLE :: p  (:,:,:,:,:) ! local air pressure (mb)
  REAL,    ALLOCATABLE :: e  (:,:,:,:,:) ! input TKE (J/kg), output V'2 (m2/s2)
  REAL,    ALLOCATABLE :: x  (:,:,:,:,:) ! vertical turbulence      W'2 (m2/s2)
  REAL,    ALLOCATABLE :: h  (:,:,:,:,:) ! horizontal turbulence    U'2 (m2/s2)

!dwen(20090819)***************************************
  real,    allocatable :: d  (:,:,:,:,:)
  real,    allocatable :: tlrams  (:,:,:,:,:)
  real,    allocatable :: sigwrams  (:,:,:,:,:)
  real,    allocatable :: cfxup1  (:,:,:,:,:)
  real,    allocatable :: cfxup2  (:,:,:,:,:)
  real,    allocatable :: cfxdn1  (:,:,:,:,:)
  real,    allocatable :: dfxup1  (:,:,:,:,:)
  real,    allocatable :: dfxup2  (:,:,:,:,:)
  real,    allocatable :: dfxdn1  (:,:,:,:,:)
  real,    allocatable :: efxup1  (:,:,:,:,:)
  real,    allocatable :: efxup2  (:,:,:,:,:)
  real,    allocatable :: efxdn1  (:,:,:,:,:)
  real,    allocatable :: tke  (:,:,:,:,:)
  real,    allocatable :: tl  (:,:,:,:,:)
  real,    allocatable :: sigw  (:,:,:,:,:)
  real,    allocatable :: dmass  (:,:,:,:,:)
  real,    allocatable :: xm  (:,:,:,:,:)
  real,    allocatable :: hm  (:,:,:,:,:)
!*******************************************

  REAL,    ALLOCATABLE :: p0 (:,:  ,:,:) ! model sfc press (mb)
  REAL,    ALLOCATABLE :: rt (:,:  ,:,:) ! rainfall total (m)
  REAL,    ALLOCATABLE :: u0 (:,:  ,:,:) ! low-level u horizontal wind (m/s)
  REAL,    ALLOCATABLE :: v0 (:,:  ,:,:) ! low-level v horizontal wind (m/s)
  REAL,    ALLOCATABLE :: t0 (:,:  ,:,:) ! low-level ambient temp (deg-K)

                                         ! u,v flux replaced by U* in prfcom
  REAL,    ALLOCATABLE :: uf (:,:  ,:,:) ! u momentum flux (n/m2) or (kg/m2-s)
  REAL,    ALLOCATABLE :: vf (:,:  ,:,:) ! v momentum flux (n/m2) 
  REAL,    ALLOCATABLE :: hf (:,:  ,:,:) ! sensible heat flux (w/m2)
  REAL,    ALLOCATABLE :: sf (:,:  ,:,:) ! normalized stability function 
  REAL,    ALLOCATABLE :: zi (:,:  ,:,:) ! mixed layer depth  (m)     
  REAL,    ALLOCATABLE :: ds (:,:  ,:,:) ! downward shortwave flux (w/m2)
  REAL,    ALLOCATABLE :: ss (:,:  ,:,:) ! static stability (1/s2)       

!dwen(20090819) *******************************
  real,    allocatable :: rc (:,:  ,:,:)
  real,    allocatable :: tc (:,:  ,:,:)
  real,    allocatable :: sw (:,:  ,:,:)
  real,    allocatable :: lc (:,:  ,:,:)
  real,    allocatable :: sm (:,:  ,:,:)
  real,    allocatable :: w0 (:,:  ,:,:)
  real,    allocatable :: muu (:,:  ,:,:)
  real,    allocatable :: muv (:,:  ,:,:)
  real,    allocatable :: mu (:,:  ,:,:)
  real,    allocatable :: msfu (:,:  ,:,:)
  real,    allocatable :: msfv (:,:  ,:,:)
  real,    allocatable :: msft (:,:  ,:,:)
  real,    allocatable :: zloc (:,:  ,:,:)
  real,    allocatable :: lf (:,:  ,:,:)
!******************************************

  REAL,    ALLOCATABLE :: gx (:,:    ,:) ! grid dist (m)
  REAL,    ALLOCATABLE :: gy (:,:    ,:) ! grid dist (m)
  REAL,    ALLOCATABLE :: z0 (:,:    ,:) ! roughness length (m)
  REAL,    ALLOCATABLE :: zt (:,:    ,:) ! terrain elevation (m)
  INTEGER, ALLOCATABLE :: lu (:,:    ,:) ! land use category

  INTEGER              :: point(2)       ! index pointer for each grid

!dwen(20090824) PGF: An identifier appearing in a SAVE statement must be a local variable or array.
!dwen(20090824)  SAVE point,u,v,w,a,t,q,e,p,x,h,  p0,rt,u0,v0,t0,  &
!dwen(20090824)       uf,vf,hf,sf,zi,ds,ss,       gx,gy,z0,zt,lu

END MODULE metval
