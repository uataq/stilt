!###############################################################################
! GEMEQN - Main grid finite difference equations gives the solution for all
! grid points in x except the first (south pole) and last (north pole) y points.
! The south and north pole regions are solved in subroutine bounds.
!-------------------------------------------------------------------------------
! PRIMARY METEOROLOGICAL VARIABLES
! UUU - U WIND COMPONENTS    M/S
! VVV - V WIND COMPONENTS    M/S
! TTT - TEMPERATURE          KELVIN
! MMM - MOISTURE (RH)        PERCENT
! HHH - HEIGHT MSL OF FIELD  METERS
! WWW - INPUT AS DIVERGENCE  1/SEC
! KKK - VERTICAL MIXING      M2/SEC
! RRR - LOCAL AIR DENSITY    KG/M3
! QQQ - MIXING RATIO         g/KG
!-------------------------------------------------------------------------------
! LAST REVISED: 21 May 2008 (RRD) - initial version
!               02 Feb 2009 (RRD) - restructure loop outside of routine
!-------------------------------------------------------------------------------

SUBROUTINE gemeqn 

  USE gemkon   
  USE gemvar
  USE gemcfg

  IMPLICIT NONE

  REAL*8    :: HM                     ! adjusted horizontal mixing
  REAL*8    :: XDYD                   ! ratio of (E-W)/(N-S) grid distances  
  REAL*8    :: QSUM                   ! mass summation
  REAL*8    :: DQDT                   ! rate of mass change
  REAL*8    :: DIST,DZI,DZT,DZB       ! various grid distances
  REAL*8    :: TRX,TRY,TRZ            ! transport (advection) fluxes
  REAL*8    :: DFX,DFY,DFT,DFB        ! diffusion fluxes

  INTEGER*4 :: n                      ! integration loop counter 
  INTEGER*4 :: numb                   ! flux averaging grid point counter
  INTEGER*4 :: i,j,k,m                ! basic loop indicies for 3D grid
  INTEGER*4 :: ip,ii                  ! loop indices for sub-grid summation
  INTEGER*4 :: im1,ip1                ! plus and minus one indicies
  INTEGER*4 :: nx,ny,nz,np,kgrd       ! 3D grid dimensions

  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  IF(KINIT.LT.0)RETURN

! initial system mass prior to computation
  qtot = 0.0
  qsum = 0.0              

  DO K=1,nz
  DO J=1,ny

! skip to the first horizontal grid point of each multi-meteo concentration cell
! process the first and last point for east-west fluxes

  DO I=1,nx,ngp(j)

!    ratio is just COS(latitude) - horizontal mixing reduces with latitude
     XDYD=GSX(i,j)/GSY(i,j)
     HM=HMIX*XDYD

!    tropical convection enhancement (0.87 = 30 deg lat)
     if(xdyd.ge.0.87)hm=hm*2.0

!    Except near the poles, the number of concentration grid points (ngp) per
!    meteorological grid point is one.

     TRX=0.0
     DFX=0.0

     IF(j.GT.1.AND.j.LT.ny)THEN
!       Horizontal advection use upwind gradient and for multi-meteorology 
!       concentration grid cells, advection is computed from the W-E ends.
!       The polar cell is excluded for W-E advection.

        DIST=GSX(I,J)*DBLE(NGP(J))     ! adjusted grid distance

        IM1=I-1                        ! index minus one
        IP1=I+1                        ! index plus one
        IF(I.EQ.1)IM1=nx               ! cyclic boundary adjustments
        IF(I.EQ.nx)IP1=1

!       west to east advection
        IF(UUU(I,J,K)*DELTA.GT.0.0)THEN
           TRX=UUU(I,J,K)*(QQQ(I,J,K)-QQQ(IM1,J,K))/DIST
        ELSEIF(UUU(I,J,K)*DELTA.LT.0.0)THEN
           TRX=UUU(I,J,K)*(QQQ(IP1,J,K)-QQQ(I,J,K))/DIST
        END IF

!       west to east horizontal diffusion
        DFX=HM*(QQQ(IP1,J,K)-2.0*QQQ(I,J,K)+QQQ(IM1,J,K))/DIST/DIST

!       set the number of meteo cells to average per conc cell
        numb=ngp(j)
     ELSE
!       at the poles process only one cell
        numb=1
     END IF

     TRY=0.0
     DFY=0.0
     TRZ=0.0
     DFT=0.0
     DFB=0.0

     ii=0
     DO WHILE(ii.LT.numb)

!       Near the poles where ngp>1 the fluxes are averaged for each
!       concentration grid cell.

        ii=ii+1
        ip=ii+i-1

!       south to north advection and diffusion

        IF(j.GT.1.AND.j.LT.ny)THEN
!          interior of the grid       
           IF(VVV(IP,J,K)*DELTA.GT.0.0)THEN
              TRY=TRY+VVV(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP,J-1,K))/GSY(IP,J)
           ELSEIF(VVV(I,J,K)*DELTA.LT.0.0)THEN
              TRY=TRY+VVV(IP,J,K)*(QQQ(IP,J+1,K)-QQQ(IP,J,K))/GSY(IP,J)
           END IF
           DFY=DFY+HM*(QQQ(IP,J+1,K)-2.0*QQQ(IP,J,K)+QQQ(IP,J-1,K))/GSY(IP,J)/GSY(IP,J)

!       At the poles advection and diffusion is computed only for one grid cell. At other
!       grid cells rate terms are computed as the average of all meteo cells that are in
!       that concentration grid cell. The gradients are computed across the pole such that
!       a point at (I,J+1) or (I,J-1) is the same as a point at (I+NX/2,J) 

        ELSEIF(j.EQ.1.AND.IP.EQ.1)THEN
!          south pole
           IF(VVV(IP,J,K)*DELTA.GT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP+NX/2,J+1,K))/GSY(IP,J)
           ELSEIF(VVV(I,J,K)*DELTA.LT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP,J+1,K)-QQQ(IP,J,K))/GSY(IP,J)
           END IF
           DFY=HM*(QQQ(IP,J+1,K)-2.0*QQQ(IP,J,K)+QQQ(IP+NX/2,J+1,K))/GSY(IP,J)/GSY(IP,J)

        ELSEIF(J.EQ.NY.AND.IP.EQ.1)THEN
!         north pole
           IF(VVV(IP,J,K)*DELTA.GT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP,J-1,K))/GSY(IP,J)
           ELSEIF(VVV(I,J,K)*DELTA.LT.0.0)THEN
              TRY=VVV(IP,J,K)*(QQQ(IP+NX/2,J-1,K)-QQQ(IP,J,K))/GSY(IP,J)
           END IF
           DFY=HM*(QQQ(IP,J-1,K)-2.0*QQQ(IP,J,K)+QQQ(IP+NX/2,J-1,K))/GSY(IP,J)/GSY(IP,J)
        END IF
 
!       vertical diffusion and advection

        IF (k.GT.1.AND.k.LT.nz) THEN
!          interior of the grid
           DZT=MAX(DZMIN,HHH(IP,J,K+1)-HHH(IP,J,K))
           DZI=MAX(DZMIN,0.5*(HHH(IP,J,K+1)-HHH(IP,J,K-1)))
           DZB=MAX(DZMIN,HHH(IP,J,K)-HHH(IP,J,K-1))

           IF(WWW(IP,J,K)*DELTA.GT.0.0)THEN
              TRZ=TRZ+WWW(IP,J,K)*(QQQ(IP,J,K)-QQQ(IP,J,K-1))/DZB
           ELSE
              TRZ=TRZ+WWW(IP,J,K)*(QQQ(IP,J,K+1)-QQQ(IP,J,K))/DZT
           END IF
           DFT=DFT+KKK(IP,J,K)*(QQQ(IP,J,K+1)-QQQ(IP,J,K))/DZT/DZI
           DFB=DFB+KKK(IP,J,K-1)*(QQQ(IP,J,K)-QQQ(IP,J,K-1))/DZB/DZI

        ELSEIF (k.EQ.1) THEN
!          bottom boundary
           DZI=MAX(DZMIN,HHH(IP,J,2)-HHH(IP,J,1))
           TRZ=TRZ+WWW(IP,J,2)*(QQQ(IP,J,2)-QQQ(IP,J,1))/DZI
           DFT=DFT+KKK(IP,J,2)*(QQQ(IP,J,2)-QQQ(IP,J,1))/DZI/DZI

        ELSE
!          top boundary
           DZI=MAX(DZMIN,HHH(IP,J,nz)-HHH(IP,J,nz-1))
           TRZ=TRZ+WWW(IP,J,nz)*(QQQ(IP,J,nz)-QQQ(IP,J,nz-1))/DZI
           DFB=DFB+KKK(IP,J,nz-1)*(QQQ(IP,J,nz)-QQQ(IP,J,nz-1))/DZI/DZI
        END IF

     END DO

!    add all final derivative terms
     IF(DELTA.GT.0.0)THEN
        DQDT=((DFX-TRX)+(DFY-TRY-TRZ+DFT-DFB)/DBLE(numb))*DBLE(DELTA*60.0)
     ELSE
        DQDT=((DFX+TRX)+(DFY+TRY+TRZ+DFT-DFB)/DBLE(numb))*DBLE(abs(DELTA)*60.0)
     END IF

!    and update the concentration working array
     CCC(I,J,K)=QQQ(I,J,K)+DQDT

!    any instability go back to the orignal value
     IF(CCC(I,J,K).LT.0.0) CCC(I,J,K)=QQQ(I,J,K)

!    for computational purposes multi-meteo-cell concentration grid cells 
!    near the poles all have the same concentration values
     IF(ngp(j).GT.1) CCC(I+1:I+NGP(J)-1,J,K) = CCC(i,j,k)

!    diagnostic mass summations
     IF(kmass.GT.0)THEN
        DO IP=I,I+NGP(J)-1
           IF(k.LT.nz)THEN
              DZI=HHH(ip,j,k+1)-HHH(ip,j,k)
           ELSE
              DZI=HHH(ip,j,k)-HHH(ip,j,k-1)
           END IF
           qtot=qtot+qqq(ip,j,k)*RRR(ip,j,k)*GXY(ip,j)*DZI   ! before differencing     
           qsum=qsum+ccc(ip,j,k)*RRR(ip,j,k)*GXY(ip,j)*DZI   ! after  differencing 

!          error check (effectively disabled due to previous instability test)
           IF(ccc(ip,j,k).LT.0.0)THEN
              WRITE(*,*)'At grid location: ',ip,j,k
              WRITE(*,*)'Concentration <0: ',ccc(ip,j,k),qqq(ip,j,k)
              WRITE(*,*)'Surface pressure: ',sfc(ip,j)
              WRITE(*,'(7A10)')'sig','ppp','hhh','ttt','www','mmm','rrr'
              DO m=nz,1,-1
                 WRITE(*,'(7E10.3)') sig(m),ppp(m),hhh(ip,j,m),ttt(ip,j,m),www(ip,j,m), &
                                                   mmm(ip,j,m),rrr(ip,j,m)
              END DO 
              READ(*,*)
           END IF
        END DO
     END IF

  END DO
  END DO
  END DO

  IF(kmass.EQ.2.AND.qsum.NE.0.0) THEN
!    computation mass adjustment factor to maintain system mass
     qqq = qtot*ccc/qsum
  ELSE
     qqq = ccc
  END IF

END SUBROUTINE gemeqn
