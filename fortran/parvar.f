!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PARVAR           PARticle VARiance returns random component
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PARTICLE VARIANCE RETURNS RANDOM COMPONENT GIVEN VELOCITY SIGMA.
!   RESULTS WITH A GAUSSIAN DISTRIBUTION AND A STANDARD
!   DEVIATION OF THE INPUT SIGMA AND A MEAN OF ZERO.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 19 Dec 1998 (RRD) - replace iterative solution with table
!                 20 Jan 1999 (MDC) - fixed array bound test failure
!                 06 Jul 2000 (RRD) - using ran1 and gasdev functions from
!                                     Numerical Recipes 2.0 - thanks to suggesti
!                                     from Gerbig & Lin, Dept Earth & Planetary
!                                     Harvard University
!                 20 Sep 2000 (RRD) - fortran90 upgrade
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 25 May 2006 (AS)  - option for variable seed
!
! USAGE:  CALL PARVAR(SIGMA,VELOC,ISEED)  
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN-90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE parvar (sigma,veloc,iseed) 

  IMPLICIT NONE
  REAL,         INTENT(IN )    :: sigma 
  REAL,         INTENT(OUT)    :: veloc
  INTEGER,      INTENT(INOUT)  :: iseed    
  REAL                         :: gasdev

  veloc = GASDEV(iseed)*sigma

END SUBROUTINE parvar

!-------------------------------------------------------------------------------
! returns a uniform random deviate between 0.0 and 1.0. Set IDUM to any
! negative value to initialize or reinitialize the sequence
! see Numerical Recipes, The Art of Scientific Computing, Press et al.,
! Cambridge University Press
!-------------------------------------------------------------------------------

REAL FUNCTION ran1(idum)

  IMPLICIT NONE

  INTEGER             :: idum 
  INTEGER, PARAMETER  :: ia       = 16807
  INTEGER, PARAMETER  :: im       = 2147483647
  INTEGER, PARAMETER  :: iq       = 127773
  INTEGER, PARAMETER  :: ir       = 2836
  INTEGER, PARAMETER  :: ntab     = 32
  INTEGER, PARAMETER  :: ndiv     = 1+(im-1)/ntab
  REAL,    PARAMETER  :: am       =  1.0/im
  REAL,    PARAMETER  :: eps      =  1.2e-7
  REAL,    PARAMETER  :: rnmx     =  1.0-eps 
  INTEGER             :: j,k
  INTEGER             :: iv(ntab) = 0
  INTEGER             :: iy       = 0

  SAVE iv,iy

  IF (idum.LE.0.OR.iy.EQ.0) THEN
      idum=max(-idum,1)
      DO j=ntab+8,1,-1
         k=idum/iq
         idum=IA*(idum-k*iq)-ir*k
         IF (idum.lt.0) idum=idum+IM
         IF (j.le.NTAB) iv(j)=idum
      END DO      
      iy=iv(1)
  END IF     

  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  IF (idum.lt.0) idum=idum+IM
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  ran1=min(am*float(iy),rnmx)

END FUNCTION

!-------------------------------------------------------------------------------
! returns a normally distributed deviate with zero mean and unit
! variance, using RAN1(IDUM) as the source of the uniform deviates
! see Numerical Recipes, The Art of Scientific Computing, Press et al.,
! Cambridge University Press
!-------------------------------------------------------------------------------

REAL FUNCTION gasdev(idum)

  IMPLICIT NONE

  INTEGER        :: idum
  INTEGER        :: iset          = 0
  REAL           :: fac,gset,rsq
  REAL           :: v1,v2,ran1

  SAVE iset,gset

  IF (iset.eq.0) THEN 
     rsq=0.0
     DO WHILE (rsq.ge.1..or.rsq.eq.0.)
        v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
     END DO
     fac=sqrt(-2.*log(rsq)/rsq)
     gset=v1*fac
     gasdev=v2*fac
     iset=1
  ELSE    
     gasdev=gset
     iset=0
  END IF   

END FUNCTION 
