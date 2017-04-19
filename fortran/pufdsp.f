!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PUFDSP           PUFf DiSPersion subgrid puff dispersion
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PUFF DISPERSION ASSUMES LINEAR GROWTH OF HORIZONTAL DISTRIBUTION
!   SQUARE ROOT VERTICAL GROWTH. GROWTH RATES ARE COMPUTED FROM THE
!   STANDARD DEVIATION OF THE TURBULENT VELOCITY.  TURBULENT VELOCITY
!   VARIANCES MAY BE OBTAINED FROM THE METEOROLOGICAL MODEL DIRECTLY
!   PARAMETERIZED FROM THE HORIZONTAL AND VERTICAL DIFFUSIVITY,
!   SIGMA^2 = H / TL.  VERTICAL GROWTH MAY REQUIRE FINER TIME STEPS
!   MAINTAIN COMPUTATIONAL STABILITY.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 17 Nov 1997 (RRD)
!                 18 Aug 1998 (RRD) - isotroptic turbulence option
!                 20 Apr 1999 (RRD) - added terrain compression
!                 26 Aug 1999 (RRD) - pass vertical grid in common block
!                 06 Jul 2000 (RRD) - time step to seconds
!                                     correction to vertical index
!                 14 Aug 2000 (RRD) - quadratic equation test
!                 05 Sep 2000 (RRD) - fortan90 upgrade
!                 16 Mar 2001 (RRD) - argument list change
!                 23 May 2002 (RRD) - time step based upon entire profile
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 15 Sep 2003 (RRD) - more focused test on hdwp
!                 05 Nov 2003 (RRD) - convert from mixing to turbulence
!                 10 Aug 2004 (RRD) - square root horizontal dispersion
!                 13 Oct 2004 (RRD) - horizontal puff dispersion options
!                 12 Oct 2005 (RRD) - langrangian sampling test hdwp=6
!                 30 Aug 2007 (RRD) - fixed gaussian kpuff=1 
!                 17 Oct 2007 (RRD) - force positive time step
!                 18 Jan 2008 (RRD) - test for splitting turned off
!                 25 Jan 2008 (RRD) - lower limit on sqrt growth rate
!                 03 Mar 2008 (RRD) - correction for sqrt age=0 test
!                 02 Jul 2008 (RRD) - refined puff/particle distribution
!                 08 Aug 2008 (RRD) - time scales become variables
!
! USAGE:  CALL PUFDSP(KPUFF,UMIX,VMIX,DT,ZMDL,ZSFC,NLVL,WMIX,ZSG,PAGE,ZPOS,
!                     SIGH,SIGW,HDWP,ZNDX)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

!dwen(20090817) used ind_zsg to calculate ZX
SUBROUTINE PUFDSP(KPUFF,UMIX,VMIX,DTM,ZMDL,ZSFC,NLVL,WMIX,ZSG,PAGE,ZPOS, &
                  SIGH,SIGW,HDWP,ZNDX)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)    :: kpuff      ! linear (0) or sqr root (1) dispersion
  REAL,    INTENT(IN)    :: umix       ! u-component turbulence (m2/s2)
  REAL,    INTENT(IN)    :: vmix       ! v-component turbulence (m2/s2)
  REAL,    INTENT(IN)    :: dtm        ! horizontal time step from advection
  REAL,    INTENT(IN)    :: zmdl       ! top of computational domain (meters)
  REAL,    INTENT(IN)    :: zsfc       ! surface terrain elevation (m)
  INTEGER, INTENT(IN)    :: nlvl       ! number of levels in subgrid
  REAL,    INTENT(IN)    :: wmix (:)   ! vertical turbulence profile (m2/s2)
  REAL,    INTENT(IN)    :: zsg  (:)   ! internal model sigma levels
  INTEGER, INTENT(IN)    :: hdwp       ! Horizontal distribution pollutant
  INTEGER, INTENT(IN)    :: page       ! Puff age in minutes
  REAL,    INTENT(IN)    :: zndx       ! fraction vertical index position
  REAL,    INTENT(INOUT) :: zpos       ! puff center height (sigma)
  REAL,    INTENT(INOUT) :: sigh       ! horiz (meters) sigma
  REAL,    INTENT(INOUT) :: sigw       ! vert sigma (sigma)

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL             :: VSCALE            ! Vertical Lagrangian time scales 
  REAL             :: HSCALE            ! Horizontal Lagrangian time scales
  REAL, PARAMETER  :: SIGR   = 1.54     ! vertical extent
  REAL, PARAMETER  :: root2  = 1.414    ! sqrt(2)

  REAL    :: dsh,dt,zx,sgb,zz,cc,aa,bb,delt,sgt,amix,delz
  REAL    :: wvvt,wvvb,vsigt,vsigb,ds2t,ds2b,dist,age,rootx
  INTEGER :: knum,kt,kz,kb,hdwpx 

!-------------------------------------------------------------------------------

! common block to pass parameters for vertical grid
  COMMON /ZZTOKK/ AA,BB,CC
  COMMON /stblen/ vscale,hscale

!-------------------------------------------------------------------------------

! particles with puff distributions calculated in this routine
  HDWPX=HDWP                            ! simple mode  
  IF(HDWPX.GE.100)HDWPX=MOD(HDWP/10,10) ! complex mode

! exit routine for pure particles
  IF(HDWPX.EQ.0.OR.HDWPX.EQ.5.OR.HDWPX.EQ.6)RETURN

! negative age set for big puffs when splitting limit reached
! reset to positive after splitting in pufsph or pufspv
  IF(PAGE.LT.0)RETURN

! convert time step to sec
  DT=ABS(DTM*60.0)
! puff age from minutes to seconds
  AGE=PAGE*60.0
! save initial vertical index position
  ZX=ZNDX

!-------------------------------------------------------------------------------
! Horizontal puff diffusion (hdwp = 1,2,3,4) defaults (kpuff=0) to a linear
! growth rate which applies prior to splitting to represent subgrid turbulence.
! Puff splitting simulates sqrt random walk approach. In the 3D puff mode,
! splitting occurs frequently (age reset to 0) so that the kpuff=0 or =1
! options give very similar results. However in the hybrid mode only horizontal
! splitting applies and hence it becomes sensitive to the size of the        
! meteorological data grid.
!-------------------------------------------------------------------------------

! horizontal turbulent velocity sigma (m/s)
  DSH=SQRT(UMIX+VMIX)

! Dispersion equation from Hysplit3 
! linear horizontal diffusion summation (m/sec => meters)
! D(Sy) = sqrt(2) * Sv * Dt
  IF(KPUFF.EQ.0)THEN 
     SIGH=SIGH+ROOT2*DSH*DT

! Revised equation to match particle dispersion
! square root horizontal diffusion summation (m/sec => meters)
! D(Sy) = sqrt( Tl / 2t ) * Sv * Dt
  ELSE
!    linear rate (kpuff=0) occurs at about one hour, 
!    maintain minimum growth rate corresponding to age of about 6 h   
     IF(AGE.GT.0.0) THEN
        ROOTX=MAX(0.25, 0.5*HSCALE/AGE)
     ELSE
        ROOTX=ROOT2*ROOT2
     END IF
     SIGH=SIGH+SQRT(ROOTX)*DSH*DT
  END IF

! particle in vertical return now
  IF(HDWPX.EQ.3.OR.HDWPX.EQ.4)RETURN

!-------------------------------------------------------------------------------
! vertical puff diffusion (hdwp = 1,2)
!-------------------------------------------------------------------------------

! determine the required vertical dispersion time step (sec)
! by computing the minimum time step for vertical stability
  DELZ=(ZMDL-ZSFC)*(ZSG(1)-ZSG(2))
  DELT=MIN(DT,0.125*DELZ*DELZ/MAXVAL(WMIX)/VSCALE)
  DELT=MAX(1,NINT(DELT))

! round down time step to even multiple of DT
  DO WHILE (MOD(INT(DT),INT(DELT)).NE.0.AND.INT(DELT).GT.1)
     DELT=DELT-1.0
  END DO

! go through iterations to match external DT
  DO KNUM=1,NINT(DT/DELT)

!    index at bottom and top of puff
     SGT=MAX(ZPOS-SIGR*SIGW, 0.0)
     SGB=MIN(ZPOS+SIGR*SIGW, 1.0)

!    compute vertical index (at point just below and just above)
!    from integer index based upon quadratic relation between
!    height (or sigma) and array index position

!***************************************
!dwen(20090817): copy the following lines from STILT
! use the ind_zsg routine instead to allow for externally specified zsg levels:
!!$         ZZ=ZMDL*(1.0-DMIN1(DBLE(1.),SGB))
!!$         ZX=(-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)
         call ind_zsg(zmdl,zsg,nlvl,sgb,zx,aa,bb,cc)
         KB=MIN(MAX(1,INT(ZX)),NLVL)

!        height agl based on zsfc=0
!        index equation based on zsfc=0
! use the ind_zsg routine instead to allow for externally specified zsg levels:
!!$         ZZ=ZMDL*(1.0-DMIN1(DBLE(1.0),SGT))
!!$         ZX=(-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)
         call ind_zsg(zmdl,zsg,nlvl,sgt,zx,aa,bb,cc)
         KT=MIN(MAX(1,NINT(ZX)),NLVL)

!*****************************************

!dwen(20090817)!    height agl based on zsfc=0
!dwen(20090817)     ZZ=ZMDL*(1.0-MIN(1.0,SGB))
!dwen(20090817)
!dwen(20090817)!    index equation based on zsfc=0
!dwen(20090817)     DIST=(BB*BB-4.0*AA*(CC-ZZ))
!dwen(20090817)     IF(DIST.GE.0.0)THEN
!dwen(20090817)        ZX=(-BB+SQRT(DIST))/(2.0*AA)
!dwen(20090817)     ELSE
!dwen(20090817)        ZX=1.0
!dwen(20090817)     END IF
!dwen(20090817)     KB=MIN(MAX(1,INT(ZX)),NLVL)
!dwen(20090817)
!dwen(20090817)!    height agl based on zsfc=0
!dwen(20090817)     ZZ=ZMDL*(1.0-MIN(1.0,SGT))
!dwen(20090817)
!dwen(20090817)!    index equation based on zsfc=0
!dwen(20090817)     DIST=(BB*BB-4.0*AA*(CC-ZZ))
!dwen(20090817)     IF(DIST.GE.0.0)THEN
!dwen(20090817)        ZX=(-BB+SQRT(DIST))/(2.0*AA)
!dwen(20090817)     ELSE
!dwen(20090817)        ZX=1.0
!dwen(20090817)     END IF
!dwen(20090817)     KT=MIN(MAX(1,NINT(ZX)),NLVL)

!    vertical turbulent velocity variance (m2/s2)
     WVVT=WMIX(KT)*WMIX(KT)
     WVVB=WMIX(KB)*WMIX(KB)

!    square of vertical diffusion rate of at bottom and top
     DS2T=2.0*WVVT*VSCALE
     DS2B=2.0*WVVB*VSCALE

!    compute new puff sigma's normalized vertical coordinate
     VSIGT=SQRT(SIGW*SIGW+DS2T*DELT/(ZMDL-ZSFC)/(ZMDL-ZSFC))
     VSIGB=SQRT(SIGW*SIGW+DS2B*DELT/(ZMDL-ZSFC)/(ZMDL-ZSFC))

!    updated bottom and top of puff
     SGB=MIN(ZPOS+SIGR*VSIGB, 1.0)
     SGT=MAX(ZPOS-SIGR*VSIGT, ZSG(NLVL))

!    update sigma and vertical position
     SIGW=(SGB-SGT)/SIGR/2.0
     ZPOS=0.5*(SGB+SGT)

! time step loop
  END DO

END SUBROUTINE pufdsp
