!-------------------------------------------------------------------------------
! DEFMETO.INC 
!-------------------------------------------------------------------------------
! Defines meteorological varaibles that are returned from the advection step at
! the last particle or puff position. 
!-------------------------------------------------------------------------------
! LAST REVISED: 13 Jan 1997 (RRD)
!               04 Mar 1999 (RRD) - included more subgrid information
!               06 Apr 1999 (RRD) - added pressure profile
!               19 Apr 1999 (RRD) - ground surface
!               05 Jul 2000 (RRD) - dimensioned TMRK variable
!               20 Sep 2000 (RRD) - fortran90 upgrade
!               12 Mar 2001 (RRD) - global lat lon grid option
!               16 Oct 2001 (RRD) - added marker variable flags
!               26 Feb 2002 (RRD) - downward short-wave flux for chemistry
!               23 Jul 2002 (RRD) - additional trajectory marker variables
!               31 May 2005 (RRD) - more trajectory marker variables (& Oct 05)
!               14 Oct 2005 (RRD) - aircraft sampling velocity vectors
!               08 Mar 2006 (RRD) - static stability parameter
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! define trajectory marker variables
!-------------------------------------------------------------------------------

  TYPE MARK
     SEQUENCE

     INTEGER    :: pres            ! pressure 
     INTEGER    :: tamb            ! ambient temperature
     INTEGER    :: tpot            ! potential temperature
     INTEGER    :: rain            ! precipitation
     INTEGER    :: mixd            ! mixed layer depth
     INTEGER    :: relh            ! relative humidity 
     INTEGER    :: terr            ! terrain height 
     INTEGER    :: dswf            ! downward shortwave flux

  END TYPE

!-------------------------------------------------------------------------------
! define surface advection variables for a particle
!-------------------------------------------------------------------------------

  TYPE ASET
     SEQUENCE

     REAL       :: zndx            ! vertical fractional index position
     REAL       :: zter            ! height of ground surface
     REAL       :: mixd            ! mixed layer depth
     REAL       :: gdisx           ! grid spacing (m)
     REAL       :: gdisy           ! grid spacing (m)
     REAL       :: rain            ! precipitation rate (meters/min)
     REAL       :: umix            ! u-component turbulence (m2/s2)
     REAL       :: vmix            ! v-component turbulence (m2/s2)
     REAL       :: aero            ! aerodynamic roughness length (m)
     REAL       :: ustr            ! friction velocity (m/s)
     REAL       :: psi             ! integrated stability profile
     REAL       :: ssp             ! static stability parameter (1/s2)
     REAL       :: dswf            ! downward short-wave flux (w/m2)
     REAL       :: plat,plon       ! true position of endpoint
     REAL       :: uvel,vvel       ! constant velocity vector for sampling
     REAL       :: ubar            ! low-level scalar velocity 
     REAL       :: tmrk (8)        ! trajectory endpoint meteo marker variable
     INTEGER    :: land            ! land-use category
     TYPE(mark) :: flag            ! flag to indicate marker on (1) or off (0)
     real       :: crai
     real       :: shtf
     real       :: lhtf
     real       :: tcld
     real       :: lcld
     real       :: solw
     real       :: zlocprev
     real       :: zlocnext
     real       :: denslocal
     real       :: zmlprev
     real       :: zmlnext
     real       :: hmix            ! horizontal diffusivity (m2/s)
 
  END TYPE

!-------------------------------------------------------------------------------
! define advection profile variables
!-------------------------------------------------------------------------------

  TYPE BSET
     SEQUENCE

     REAL       :: pres            ! pressure profile (mb)
     REAL       :: temp            ! temperature profile (deg K)
     REAL       :: rhfr            ! relative humidity fraction (0-1)
     REAL       :: dens            ! density profile (kg/m3)
     REAL       :: wmix            ! w-vertical turbulence profile (m2/s2)
     real       :: wmxc            ! vertical mixing profile (m2/s)
     real       :: cfxup1
     real       :: cfxup2
     real       :: cfxdn1
     real       :: tken
     real       :: sigw
     real       :: sigwprev
     real       :: raup1
     real       :: raup2
     real       :: radn1
     real       :: dfxup1
     real       :: dfxup2
     real       :: efxup1
     real       :: efxup2
     real       :: dfxdn1
     real       :: efxdn1
     real       :: uuprev
     real       :: uunext
     real       :: vvprev
     real       :: vvnext
     real       :: dmassprev
     real       :: dmassnext
     real       :: densprev
     real       :: densnext
     real       :: tempprev
     real       :: tempnext
     real       :: rhfrprev
     real       :: rhfrnext
     real       :: tl
     real       :: tlprev
  END TYPE
