!###############################################################################
! GEMDAT - METeorological DATa processing
!-------------------------------------------------------------------------------
! LAST REVISED: 16 May 2007 (RRD) - initial version
!               14 Jan 2009 (RRD) - moved terrain correction after k-loop
!-------------------------------------------------------------------------------

SUBROUTINE gemdat

  USE funits
  USE gemkon   
  USE gemvar   

  IMPLICIT NONE

  REAL*8    :: dzi,dsp
  REAL*4    :: pres,fact,zsfc
  INTEGER*4 :: i,j,k,kb,kt,nx,ny,nz,np,kgrd

  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

!--------------------------------------------------------------
! average terrain by averaging surface pressure field
!--------------------------------------------------------------

! main grid 
  DO I=2,(nx-1)
  DO J=2,(ny-1)
     avg(i,j)=0.25*(1.0-twght)*(sfc(i-1,j)+sfc(i+1,j)+sfc(i,j-1)+sfc(i,j+1))+twght*sfc(i,j)
  END DO
  END DO

! east side of edge
  I=1
  DO J=2,(ny-1)
     avg(i,j)=0.25*(1.0-twght)*(sfc(nx,j)+sfc(i+1,j)+sfc(i,j-1)+sfc(i,j+1))+twght*sfc(i,j)
  END DO

! west side of edge
  I=nx
  DO J=2,(ny-1)
     avg(i,j)=0.25*(1.0-twght)*(sfc(i-1,j)+sfc(1,j)+sfc(i,j-1)+sfc(i,j+1))+twght*sfc(i,j)
  END DO

! poles
  avg(:,1) =SUM(sfc(:,1))/nx
  avg(:,ny)=SUM(sfc(:,ny))/nx

! replace surface pressure with averaged field
  sfc=avg

!--------------------------------------------------------------
! interpolation to sigma surfaces
!--------------------------------------------------------------

  DO J=1,ny
  DO I=1,nx

!---------------------------------
!    compute the terrain height

     IF(sfc(i,j).LT.ppp(1))THEN
!       interpolate height from existing points
        zloop : DO k=2,nz
           IF(sfc(i,j).GE.ppp(k))THEN
              IF(sfc(i,j).EQ.ppp(k))THEN
                 zsfc=hhh(i,j,k)
              ELSE
                 fact=(ppp(k-1)-sfc(i,j))/(ppp(k-1)-ppp(k))
                 zsfc=fact*(hhh(i,j,k)-hhh(i,j,k-1))+hhh(i,j,k-1)
              END IF
              EXIT zloop
           END IF
        END DO zloop

     ELSEIF(sfc(i,j).GT.ppp(1))THEN
!       compute height to the ground
        zsfc=HHH(I,J,1)-0.5*(HHH(I,J,2)-HHH(I,J,1))
     ELSE
!       coincidently perfect match
        zsfc=hhh(i,j,1)
     END IF

!----------------------------------
!    compute the density profile

!##  DZI=0.0        
!##  DSP=1013.0

     DO K=1,nz

!##     thermodynamic method to recompute height fields
!##     RRR(I,J,K)=P2JM*PPP(K)/(TTT(I,J,K)*RDRY)
!##     DZI=DZI+P2JM*(DSP-PPP(K))/(RRR(I,J,K)*GRAV)
!##     HHH(I,J,K)=DZI
!##     DSP=PPP(K)

!       check height field for consistency (needed for error check on external data)
        IF((K.GT.1).AND.(HHH(I,J,K).LT.HHH(I,J,K-1)))THEN
           WRITE(*,*)'Height inversion found at (i,j)',I,J
           WRITE(*,*)' for lev: ',K,  ppp(K),  HHH(I,J,K)
           WRITE(*,*)' for lev: ',K-1,ppp(K-1),HHH(I,J,K-1)
           WRITE(*,*)' std atm: ',STDATM(KNDX(K))
           WRITE(*,*)' terrain: ',ZSFC 
           STOP 800
!##        if not stop then replace with std atm, send message to KF21 
!##        HHH(I,J,K)=STDATM(KNDX(K))
        END IF

!       use vertical box sizes between interfaces (DZI) to compute density
        IF(K.EQ.1) THEN
           DZI=HHH(I,J,K+1)-HHH(I,J,K)
           DSP=PPP(K)-PPP(K+1)
        ELSEIF (K.EQ.NZ) THEN
           DZI=HHH(I,J,K)-HHH(I,J,K-1)
           DSP=PPP(K-1)-PPP(K)
        ELSE
           DZI=0.5*(HHH(I,J,K+1)-HHH(I,J,K-1))
           DSP=0.5*(PPP(K-1)-PPP(K+1))
        END IF

!       average density (g/kg) within each vertical box
        RRR(I,J,K)=P2JM*DSP/DZI/GRAV

     END DO

!    optional adjustment for  heights relative to terrain
!    does not affect subsequent computations based upon height differences
     HHH(i,j,:)=HHH(i,j,:)-zsfc

!-----------------------------------
!    interpolate to sigma surfaces

     DO K=1,nz
!       pressure at the interpolation sigma level
        PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP        

        kt=1
        DO WHILE (kt.LT.nz.AND.PPP(kt).GT.PRES)
!          find first index above sigma level
           kt=kt+1
        END DO
        kb=MAX(1,kt-1)

!       vertical interpolation factor
        IF(kt.GT.kb)THEN
           FACT=(PPP(kb)-PRES)/(PPP(kb)-PPP(kt))
        ELSE
           FACT=0.0
        END IF
  
        VAL(k,1)=HHH(i,j,kb)-FACT*(HHH(i,j,kb)-HHH(i,j,kt))
        VAL(k,2)=TTT(i,j,kb)-FACT*(TTT(i,j,kb)-TTT(i,j,kt)) 
        VAL(k,3)=UUU(i,j,kb)-FACT*(UUU(i,j,kb)-UUU(i,j,kt)) 
        VAL(k,4)=VVV(i,j,kb)-FACT*(VVV(i,j,kb)-VVV(i,j,kt)) 
        VAL(k,5)=WWW(i,j,kb)-FACT*(WWW(i,j,kb)-WWW(i,j,kt)) 
        VAL(k,6)=MMM(i,j,kb)-FACT*(MMM(i,j,kb)-MMM(i,j,kt)) 
        VAL(k,7)=RRR(i,j,kb)-FACT*(RRR(i,j,kb)-RRR(i,j,kt)) 
     END DO

!    put data back into 3D array
     HHH(i,j,:)=VAL(:,1)
     TTT(i,j,:)=VAL(:,2)
     UUU(i,j,:)=VAL(:,3)
     VVV(i,j,:)=VAL(:,4)
     WWW(i,j,:)=VAL(:,5)
     MMM(i,j,:)=VAL(:,6) 
     RRR(i,j,:)=VAL(:,7) 

  END DO
  END DO

END SUBROUTINE gemdat
