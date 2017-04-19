!###############################################################################
! GEMDIV - Convert vertical velocity from mb/s to m/s
!-------------------------------------------------------------------------------
! LAST REVISED: 16 May 2007 (RRD) - initial version
!-------------------------------------------------------------------------------

SUBROUTINE gemdiv 

  USE gemkon   
  USE gemvar  
  USE gemcfg

  IMPLICIT NONE

  REAL*8    :: dzi,delu,delv,divg,wvel
  INTEGER*4 :: i,j,k,nx,ny,nz,np,im1,ip1,kgrd

  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  DO I=1,nx

     IM1=I-1                        ! index minus one
     IP1=I+1                        ! index plus one
     IF(I.EQ.1)IM1=nx               ! cyclic boundary adjustments
     IF(I.EQ.nx)IP1=1

     DO J=2,ny-1

        DIVG=0.0
        DO K=1,nz

!          set vertical box sizes between interfaces (DZI)
           IF(K.EQ.1) THEN
              DZI=HHH(I,J,K+1)-HHH(I,J,K)
           ELSEIF (K.EQ.NZ) THEN
              DZI=HHH(I,J,K)-HHH(I,J,K-1)
           ELSE
              DZI=0.5*(HHH(I,J,K+1)-HHH(I,J,K-1))
           END IF

           IF(kdivg.EQ.1)THEN
!             divergence in sec-1
              DELU=0.5*(UUU(IP1,J,K)-UUU(IM1,J,K))/GSX(I,J)
              DELV=0.5*(VVV(I,J+1,K)-VVV(I,J-1,K))/GSY(I,J)
              DIVG=DIVG+(DELU+DELV)*DZI

!             positive divergence implies downward motion or a negative vertical velocity
              WVEL=-DIVG*WFACT

!             vertical velocity limits
              IF(wvel.LT.-vwlim) wvel=-vwlim
              IF(wvel.GT.+vwlim) wvel=+vwlim

!             only apply divergence to those levels that originally had a defined velocity
              IF(www(i,j,k).NE.0.0) www(i,j,k)=wvel

           ELSE
!             convert to m/s decreasing pressure with height is a positive vertical velocity
              WVEL=-WFACT*P2JM*WWW(I,J,K)/RRR(I,J,K)/GRAV
              IF(wvel.LT.-vwlim) wvel=-vwlim
              IF(wvel.GT.+vwlim) wvel=+vwlim
              WWW(I,J,K)=wvel
           END IF

        END DO
     END DO

!    pole values the same as adjacent "j" grid point
     WWW(I,1,:)=WWW(I,2,:)
     WWW(I,NY,:)=WWW(I,NY-1,:)

  END DO

END SUBROUTINE gemdiv 
