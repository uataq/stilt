!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CGRELL           Convective redistribution using GRELL scheme
!   PRGMMR:    CHRISTOPH GERBIG   DATE:03-13-03
!
! ABSTRACT:
!   Convective redistribution of particles using updraft & downdraft mass
!   fluxes, entrainment and detrainment fluxes. The vertical advection
!   is done on a fixed vertical grid, with dynamic internal timestep.
!   First, particles are assigned to updraft, downdraft or environment,
!   then they are advected to the next grid location. Then they are advected
!   using the fixed grid, until time is used up. Finally they are advected
!   from there position before time was up to the exact (off-grid) location
!   corresponding to the exact time at the end of the convective
!   redistribution.
!   Initially particles are assigned a cloud index (cloud updraft, downdraft,
!   or environment. When the area coverage changes, there is a probability
!   that they change there cloud index.
!
!
! PROGRAM HISTORY LOG:
!   LAST REVISED:
!
! USAGE:  CALL CGRELL(CFXUP1,CFXUP2,CFXDN1,DFXUP1,DFXUP2,EFXUP1,EFXUP2,
!              DFXDN1,EFXDN1,RAUP1,RAUP2,RADN1,AREA,AREAPRU,AREAPRD,DENS,NLVL,
!              ZPROFM,Z1,Z2,ICNDX1,ICNDX2,ZX,DT,BACK)
!   INPUT ARGUMENT LIST:
!     ?FX*  - real      Cloud/Entrainment/Detrainment Mass Fluxes
!                 (Updraft, downdraft; 1:deep; 2: shallow)
!                 UNITS: kg/m^2s; convention: upward=positive
!                       TO GET MASSFLUX (kg/s):
!                       MASSFLUX=?FX* * DX*DY (up/dndraft & hor. fluxes)
!                       numbers represent ro*w*a w/ a=fractional coverage
!                       to get mass flux[kg/s]: ?FX* * area (gridcell)
!     RA*   - real      updraft area coverage (fract.),
!                 (Updraft, downdraft; 1:deep; 2: shallow)
!     AREA  - real      gridcell area (m^2)
!     AREAPRU     - real      area coverage updr. (fract.) of particle during last call
!     AREAPRD     - real      area coverage downdr. (fract.) of particle during last call
!     DENS  - real      Density profile [kg/m^3]
!     NLVL  - int vertical grid dimension
!     ZPROFM      - real      profile of level altitudes above ground, FLX
!                 START: First level above ground (surface not incl.)
!     Z1    - real      old position at time t
!     ICNDX1      - int cloud index (whether in up/dndrft/environment)
!                 at time t
!     DT    - real      integration step (minutes)
!     BACK      - log   logical flag to indicate direction
!   OUTPUT ARGUMENT LIST:
!     Z2    - real      new position at time t+dt
!     ICNDX2      - real      cloud index (whether in up/dndrft/environment)
!                 at time t+dt
!     AREAPRU     - real      area coverage updr. (fract.) of particle after convection
!     AREAPRD     - real      area coverage downdr. (fract.) of particle after convection
!     ZX    - int last estimate of vertical index
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  LINUX
!
! $Id: cgrell.f,v 1.4 2009/12/11 00:03:24 jel Exp $
!
!$$$

       SUBROUTINE CGRELL(CFXUP1,CFXUP2,CFXDN1,DFXUP1,DFXUP2,EFXUP1,     &
     &                   EFXUP2,DFXDN1,EFXDN1,RAUP1,RAUP2,RADN1,AREA,   &
     &                   AREAPRU,AREAPRD,DENS,NLVL,ZPROFM,Z1,Z2,        &
     &                   ICNDX1,ICNDX2,ZX,DT,BACK)

!dwen(20090825) **************************
!      IMPLICIT REAL*8 (A-H,O-Z)
!      LOGICAL BACK
!      REAL*8   CFXUP1(NLVL),CFXUP2(NLVL),CFXDN1(NLVL)
!      REAL*8   DFXUP1(NLVL),DFXUP2(NLVL),DFXDN1(NLVL)
!      REAL*8   EFXUP1(NLVL),EFXUP2(NLVL),EFXDN1(NLVL)
!      REAL*8   RAUP1(NLVL),RAUP2(NLVL),RADN1(NLVL)
!      REAL*8   DENS(NLVL),ZPROFM(NLVL),ZPROFT(NLVL)
!      REAL*8   DZT(NLVL)
!
!      INTEGER ICNDX1,ICNDX2,ZX
!
!!    grid size arrays for Fluxes
!      REAL*8 FUP(3,NLVL),FLE(3,NLVL),FRI(3,NLVL),FDN(3,NLVL),FT(3,NLVL)
!      REAL*8 FDN1(3,NLVL),FI(3,NLVL),FBUD(3,NLVL),FUPS(3,NLVL)
!
!!    grid size arrays for probabilities
!      REAL*8 PUP(3,NLVL),PLE(3,NLVL),PRI(3,NLVL),PDN(3,NLVL),PTO(3,NLVL)
!
!!    grid size arrays for w's, density and areas
!      REAL*8 DENSM(3,NLVL),AREAA(3,NLVL)
!      REAL*8 WUP(3,NLVL),WDN(3,NLVL),WTT(3,NLVL)
!
!!    grid size arrays for MASS and DZ's
!      REAL*8 DZTA(3,NLVL)
!      DATA PI/3.1415926/
!      DATA ONE/1.0/
!********************************************

      USE funits

      implicit  NONE
     
      real,   intent(in)    :: cfxup1(:),cfxup2(:),cfxdn1(:)
      real,   intent(in)    :: dfxup1(:),dfxup2(:),dfxdn1(:)
      real,   intent(in)    :: efxup1(:),efxup2(:),efxdn1(:)
      real,   intent(in)    :: raup1(:),raup2(:),radn1(:)
      real,   intent(in)    :: area 
      real,   intent(inout) :: areapru,areaprd 
      real,   intent(in)    :: dens(:)
      integer,intent(in)    :: nlvl
      real,   intent(in)    :: zprofm(:) 
      real,   intent(in)    :: z1
      real,   intent(out)   :: z2
      integer,   intent(inout) :: icndx1,zx
      integer,   intent(out)   :: icndx2
      real,   intent(in)    :: dt
      logical, intent(in)   :: back
      

      real     :: zproft(nlvl),dzt(nlvl)
      real     :: fup(3,nlvl),fle(3,nlvl),fri(3,nlvl), &
                  fdn(3,nlvl),ft(3,nlvl),fdn1(3,nlvl), &
                  fi(3,nlvl),fbud(3,nlvl),fups(3,nlvl)

      real     :: pup(3,nlvl),ple(3,nlvl),pri(3,nlvl), &
                  pdn(3,nlvl),pto(3,nlvl)

      real     :: densm(3,nlvl),areaa(3,nlvl),         &
                  wup(3,nlvl),wdn(3,nlvl),wtt(3,nlvl)

      real     :: dzta(3,nlvl)

      real     :: etime,dtn,areanow,areapr,rand,ran3, &
                  randn,pdtrup,etime2
      integer  :: izndxt,k,izndxt2


      real, parameter :: pi = 3.1415926
      real, parameter :: one = 1.0 

!     get DZT as grid spacing
      ZPROFT(2:NLVL)=0.5*(ZPROFM(1:(NLVL-1))+ZPROFM(2:NLVL))
      ZPROFT(1)=0.5*ZPROFM(1)
      DZT(1:NLVL-1)=ZPROFT(2:NLVL)-ZPROFT(1:(NLVL-1))
      DZT(NLVL)=ZPROFT(NLVL)-ZPROFT(NLVL-1)
      DZTA(1,1:NLVL)=DZT
      DZTA(2,1:NLVL)=DZT
      DZTA(3,1:NLVL)=DZT
!     get areas in array
!     relative area coverage is defined at flux levels
      AREAA(1,1:NLVL)=RAUP1(1:NLVL)+RAUP2(1:NLVL)
      AREAA(3,1:NLVL)=RADN1(1:NLVL)
      AREAA(2,1:NLVL)=1-AREAA(1,1:NLVL)-AREAA(3,1:NLVL)
!     from relative area to m^2
      AREAA(1:3,1:NLVL)=AREAA(1:3,1:NLVL)*AREA

!     put upward fluxes in 2D array
!     take into account time direction
!     ALL FLUXES are in MASSFLUX/GRIDAREA (=dxdy)
!     column 1 is updraft, 2 is environment, 3 is downdraft
      IF(.NOT.BACK)THEN
!dwen(20090825)        FUP(1,1:NLVL)=DABS(CFXUP1)+DABS(CFXUP2)
        FUP(1,1:NLVL)=ABS(CFXUP1(1:NLVL))+ABS(CFXUP2(1:NLVL))
        FUP(3,1:NLVL)=0.0
!     put downward fluxes in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
        FDN(1,1:NLVL)=0.0
!dwen(20090825)        FDN(3,1:NLVL)=-DABS(CFXDN1)
        FDN(3,1:NLVL)=-ABS(CFXDN1(1:NLVL))
!     put fluxes to "right" in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
!dwen(20090825)        FRI(1,1:NLVL)=DABS(DFXUP1(1:NLVL))+DABS(DFXUP2(1:NLVL))
!dwen(20090825)        FRI(2,1:NLVL)=DABS(EFXDN1(1:NLVL))
        FRI(1,1:NLVL)=ABS(DFXUP1(1:NLVL))+ABS(DFXUP2(1:NLVL))
        FRI(2,1:NLVL)=ABS(EFXDN1(1:NLVL))
        FRI(3,1:NLVL)=0.0
!     put fluxes to "left" in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
!dwen(20090825)        FLE(2,1:NLVL)=DABS(EFXUP1(1:NLVL))+DABS(EFXUP2(1:NLVL))
!dwen(20090825)        FLE(3,1:NLVL)=DABS(DFXDN1(1:NLVL))
        FLE(2,1:NLVL)=ABS(EFXUP1(1:NLVL))+ABS(EFXUP2(1:NLVL))
        FLE(3,1:NLVL)=ABS(DFXDN1(1:NLVL))
        FLE(1,1:NLVL)=0.0
!     reverse flows for backward
      ELSEIF(BACK)THEN
!dwen(20090825)        FDN(1,1:NLVL)=-DABS(CFXUP1)-DABS(CFXUP2)
        FDN(1,1:NLVL)=-ABS(CFXUP1(1:NLVL))-ABS(CFXUP2(1:NLVL))
        FDN(3,1:NLVL)=0.0
!     put downward fluxes in 2D array
!     column 1 is updraft, 2 is environment, 3 is downdraft
        FUP(1,1:NLVL)=0.0
!dwen(20090825)        FUP(3,1:NLVL)=DABS(CFXDN1)
!dwen(20090825)        FLE(2,1:NLVL)=DABS(DFXUP1(1:NLVL))+DABS(DFXUP2(1:NLVL))
!dwen(20090825)        FLE(3,1:NLVL)=DABS(EFXDN1(1:NLVL))
        FUP(3,1:NLVL)=ABS(CFXDN1(1:NLVL))
        FLE(2,1:NLVL)=ABS(DFXUP1(1:NLVL))+ABS(DFXUP2(1:NLVL))
        FLE(3,1:NLVL)=ABS(EFXDN1(1:NLVL))
        FLE(1,1:NLVL)=0.0
!dwen(20090825)        FRI(1,1:NLVL)=DABS(EFXUP1(1:NLVL))+DABS(EFXUP2(1:NLVL))
!dwen(20090825)        FRI(2,1:NLVL)=DABS(DFXDN1(1:NLVL))
        FRI(1,1:NLVL)=ABS(EFXUP1(1:NLVL))+ABS(EFXUP2(1:NLVL))
        FRI(2,1:NLVL)=ABS(DFXDN1(1:NLVL))
        FRI(3,1:NLVL)=0.0
      END IF
!     get environment fluxes from
!     mass balance; what goes up must come down
      FUP(2,1:NLVL)=-FUP(1,1:NLVL)-FUP(3,1:NLVL)                        &
     &             -FDN(1,1:NLVL)-FDN(3,1:NLVL)
!     from net flux, seperate up and downward fluxes
!     no heavyside function, so use (sign(x)+1)/2
      WHERE(FUP(2,1:NLVL).LT.0.0)
        FDN(2,1:NLVL)=FUP(2,1:NLVL)
        FUP(2,1:NLVL)=0.0
      ELSEWHERE
        FDN(2,1:NLVL)=0.0
      END WHERE

!     put density in array:
!     DENSM at flux levels (mean of T neighbours)
      DENSM(1,1:NLVL-1)=0.5*(DENS(1:NLVL-1)+DENS(2:NLVL))
      DENSM(1,NLVL)=DENS(NLVL)
      DENSM(2,1:NLVL)=DENSM(1,1:NLVL)
      DENSM(3,1:NLVL)=DENSM(1,1:NLVL)

!     get w (w=F/ro) (WTT from T- to T-level)
!     set small values to 0.0
      WHERE(FUP.GE.1.0E-8)
        WUP=FUP*AREA/AREAA/DENSM
      ELSEWHERE
        WUP=0.0
      END WHERE
      WHERE(FDN.LE.-1.0E-8)
        WDN=FDN*AREA/AREAA/DENSM
      ELSEWHERE
        WDN=0.0
      END WHERE
      WTT=WUP+WDN

!     get probabilities to move to above, below, right and left
!     shift flux down so that FDN1 is flux leaving current level
!      FDN1=FDN
      FDN1(1:3,2:NLVL)=FDN(1:3,1:NLVL-1)
      FDN1(1:3,1)=0.0
!     set small values to 0.0
!     NEED TO USE CORRECT PROBABILITIES WHEN FLUXES ARE ZERO
!dwen(20090825)      FT=FUP+DABS(FDN1)+FRI+FLE
      FT=FUP+ABS(FDN1)+FRI+FLE
      WHERE(FT.NE.0.0)
        PUP=FUP/FT
!dwen(20090825)        PDN=DABS(FDN1)/FT
        PDN=ABS(FDN1)/FT
        PRI=FRI/FT
        PLE=FLE/FT
      ELSEWHERE
        PUP=0.0
        PDN=0.0
        PRI=0.0
        PLE=0.0
      END WHERE
!     TOTAL PROBABILITY LEAVING
      PTO=PUP+PDN+PRI+PLE

!     CHECK MASS BUDGET
!     FLUX IN:
!     shift FUP so that use flux entering from below
!      FI(1,2:NLVL)=FUP(1,1:NLVL-1)-FDN(1,2:NLVL)+FLE(2,2:NLVL)
!      FI(1,1)=-FDN(1,1)+FLE(2,1)
!      FI(2,2:NLVL)=FUP(2,1:NLVL-1)-FDN(2,2:NLVL)+FLE(3,2:NLVL)
!     :           +FRI(1,2:NLVL)
!      FI(2,1)=-FDN(2,1)+FLE(3,1)+FRI(1,1)
!      FI(3,2:NLVL)=FUP(3,1:NLVL-1)-FDN(3,2:NLVL)+FRI(2,2:NLVL)
!      FI(3,1)=-FDN(3,1)+FRI(2,1)
!     MASS BUDGET
!      FBUD=FT-FI

!      OPEN(35,FILE='flux.out')
!      do k=1,NLVL
!      WRITE(35,*)FT(1:3,k),FUP(1:3,k),FDN(1:3,k),FRI(1:3,k),FLE(1:3,k),
!     :  WTT(1:3,k),FI(1:3,k),FBUD(1:3,k),PUP(1:3,k),PDN(1:3,k),
!     :  PRI(1:3,k),PLE(1:3,k)
!      enddo
!      close(35)


      ETIME=0.0
!dwen(20090825)      DTN=DABS(DT)
!dwen(20090825)      IF(BACK)DTN=-DABS(DT)
      DTN=ABS(DT)
      IF(BACK)DTN=-ABS(DT)

!     position before conv.: Z1
      IZNDXT=1
      DO K=1,NLVL
       IF(ZPROFM(K).LT.Z1)IZNDXT=K+1
      ENDDO

      IF(IZNDXT.EQ.1.AND.Z1.LE.ZPROFT(1))THEN
!     BELOW FIRST LEVEL
!     NOT MOVING, LEAVE
        Z2=Z1
!        ICNDX2=ICNDX1
                                          !if zero prob. to leave cell,
        ICNDX2=2
                                          !    can't have up/downdraft
        ZX=IZNDXT
                                          !    at this level
        RETURN
      ENDIF

!     Need to allow for change in ICNDX if area
!     changed (previous area cov.: AREAPRU, AREAPRD)
!     Compare previous with current area
!     When area now is smaller, need to detrain
!     Need random variable for this
    
      AREANOW=AREAA(ICNDX1,IZNDXT)/AREA
      IF(ICNDX1.EQ.1)THEN
        AREAPR=AREAPRU
      ELSEIF(ICNDX1.EQ.3)THEN
        AREAPR=AREAPRD
          !IF(ICNDX1.EQ.2) THEN
      ELSE
        AREAPR=1.0-AREAPRD-AREAPRU
      END IF
!      WRITE(*,*)"before",ICNDX1,AREANOW,AREAPR,RSEED,IZNDXT
      IF(AREANOW.GT.0.0001.AND.AREANOW.LT.AREAPR)THEN
                     !RSEED get's reset to one anyway after 1st call
        RAND=RAN3(1)
        IF(RAND.GT.AREANOW/AREAPR)THEN
!     detrain from current ICNDX1 with probability
!     p(detr)=(areapr-areanow)/areapr=1-areanow/areapr
!           from up- and downdraft to environment
          IF(ICNDX1.EQ.1.OR.ICNDX1.EQ.3)THEN
            ICNDX1=2
          ELSEIF(ICNDX1.EQ.2)THEN
!         from environment to updraft or downdraft
!         depending on random number and prev. updr. and dndr. coverage
!         transform remaining RAND to random values between 0,1
            RANDN=(RAND-AREANOW/AREAPR)/((AREAPR-AREANOW)/AREAPR)
!         probability to go to updraft: fraction of "new" updraft
!         from total "new" area in up+downdraft
            PDTRUP=(AREAA(1,IZNDXT)/AREA-AREAPRU)/                      &
     &              ((1-AREANOW)-(1-AREAPR))
            IF(AREAA(1,IZNDXT)/AREA.LT.0.0001)PDTRUP=0.0
            IF(RANDN.LT.PDTRUP) THEN
              ICNDX1=1
                !IF(RANDN.GE.PDTRUP)
            ELSE
              ICNDX1=3
            END IF
          END IF
        END IF
      END IF
!      WRITE(*,*)"after",ICNDX1,RAND,RANDN,RSEED


!     Advect particle to next grid position (T-grid)
!     need to figure out wether moving up or down

!     if zero flux leaving
!     current level:
!     leave routine, since nothing will happen
      IF(PTO(ICNDX1,IZNDXT).EQ.0.0)THEN
      Z2=Z1
!      ICNDX2=ICNDX1
                                          !if zero prob. to leave cell,
      ICNDX2=2
                                          !    can't have up/downdraft
      ZX=IZNDXT
                                          !    at this level
      AREAPRU=AREAA(1,IZNDXT)/AREA
      AREAPRD=AREAA(3,IZNDXT)/AREA
      RETURN
      ENDIF
!      WRITE(*,*)"2nd: ",z1,IZNDXT,ICNDX1

                                   !BELOW CURRENT LEVEL
      IF(Z1.LE.ZPROFT(IZNDXT))THEN
        IF(WTT(ICNDX1,IZNDXT-1).GT.0.0)THEN
!     MOVING UP
!dwen(20090825)          ETIME=ETIME+(ZPROFT(IZNDXT)-Z1)/DABS(WTT(ICNDX1,IZNDXT-1))
!dwen(20090825)          IF(ETIME.GE.DABS(DT)*60.0)THEN
!dwen(20090825)            Z2=Z1+WTT(ICNDX1,IZNDXT-1)*DABS(DT)*60.0
          ETIME=ETIME+(ZPROFT(IZNDXT)-Z1)/ABS(WTT(ICNDX1,IZNDXT-1))
          IF(ETIME.GE.ABS(DT)*60.0)THEN
            Z2=Z1+WTT(ICNDX1,IZNDXT-1)*ABS(DT)*60.0
            ICNDX2=ICNDX1
            ZX=IZNDXT
            AREAPRU=AREAA(1,ZX)/AREA
            AREAPRD=AREAA(3,ZX)/AREA
            RETURN
          ENDIF
        ELSEIF(WTT(ICNDX1,IZNDXT-1).LT.0.0)THEN
!     MOVING DOWN
          IZNDXT=IZNDXT-1
!dwen(20090825)          ETIME=ETIME+(Z1-ZPROFT(IZNDXT))/DABS(WTT(ICNDX1,IZNDXT))
!dwen(20090825)          IF(ETIME.GE.DABS(DT)*60.0)THEN
!dwen(20090825)            Z2=Z1-DABS(WTT(ICNDX1,IZNDXT))*DABS(DT)*60.0
          ETIME=ETIME+(Z1-ZPROFT(IZNDXT))/ABS(WTT(ICNDX1,IZNDXT))
          IF(ETIME.GE.ABS(DT)*60.0)THEN
            Z2=Z1-ABS(WTT(ICNDX1,IZNDXT))*ABS(DT)*60.0
            ICNDX2=ICNDX1
            IF(Z2.LE.ZPROFM(IZNDXT)) THEN
              ZX=IZNDXT
            ELSE
              ZX=IZNDXT+1
            END IF
            AREAPRU=AREAA(1,ZX)/AREA
            AREAPRD=AREAA(3,ZX)/AREA
            RETURN
          ENDIF
        ELSE
!     NOT MOVING, LEAVE
!          WRITE(*,*)'THIS IS WEIRD...(1)',ICNDX1,IZNDXT
          Z2=Z1
!        ICNDX2=ICNDX1
                                          !if zero prob. to leave cell,
          ICNDX2=2
                                          !   can't have up/downdraft
          ZX=IZNDXT
                                          !   at this level
          AREAPRU=AREAA(1,IZNDXT)/AREA
          AREAPRD=AREAA(3,IZNDXT)/AREA
          RETURN
        ENDIF
           !ABOVE CURRENT LEVEL
      ELSE
                                          !MOVING UP
        IF(WTT(ICNDX1,IZNDXT).GT.0.0)THEN
!dwen(20090825)          ETIME=ETIME+(ZPROFT(IZNDXT+1)-Z1)/DABS(WTT(ICNDX1,IZNDXT))
          ETIME=ETIME+(ZPROFT(IZNDXT+1)-Z1)/ABS(WTT(ICNDX1,IZNDXT))
          IZNDXT=IZNDXT+1
!dwen(20090825)          IF(ETIME.GE.DABS(DT)*60.0)THEN
!dwen(20090825)            Z2=Z1+WTT(ICNDX1,IZNDXT-1)*DABS(DT)*60.0
          IF(ETIME.GE.ABS(DT)*60.0)THEN
            Z2=Z1+WTT(ICNDX1,IZNDXT-1)*ABS(DT)*60.0
            ICNDX2=ICNDX1
            IF(Z2.LE.ZPROFM(IZNDXT-1)) THEN
              ZX=IZNDXT-1
            ELSE
              ZX=IZNDXT
            END IF
            AREAPRU=AREAA(1,ZX)/AREA
            AREAPRD=AREAA(3,ZX)/AREA
            RETURN
          ENDIF
                                              !MOVING DOWN
        ELSEIF(WTT(ICNDX1,IZNDXT).LT.0.0)THEN
!dwen(20090825)          ETIME=ETIME+(Z1-ZPROFT(IZNDXT))/DABS(WTT(ICNDX1,IZNDXT))
!dwen(20090825)          IF(ETIME.GE.DABS(DT)*60.0)THEN
!dwen(20090825)            Z2=Z1-DABS(WTT(ICNDX1,IZNDXT))*DABS(DT)*60.0
          ETIME=ETIME+(Z1-ZPROFT(IZNDXT))/ABS(WTT(ICNDX1,IZNDXT))
          IF(ETIME.GE.ABS(DT)*60.0)THEN
            Z2=Z1-ABS(WTT(ICNDX1,IZNDXT))*ABS(DT)*60.0
            ICNDX2=ICNDX1
            ZX=IZNDXT
            AREAPRU=AREAA(1,ZX)/AREA
            AREAPRD=AREAA(3,ZX)/AREA
            RETURN
          ENDIF
        ELSE
!     NOT MOVING, LEAVE
!          WRITE(*,*)'THIS IS WEIRD...(2)',ICNDX1,IZNDXT
          Z2=Z1
!        ICNDX2=ICNDX1
                   !if zero prob. to leave cell, can't have up/downdraft
          ICNDX2=2
          ZX=IZNDXT
          AREAPRU=AREAA(1,IZNDXT)/AREA
          AREAPRD=AREAA(3,IZNDXT)/AREA
          RETURN
        ENDIF
      ENDIF

!      OPEN(34,FILE='mass.out')

!dwen(20090825)        DO WHILE(DABS(DT)*60.0.GE.ETIME)
    !   DO WHILE(ABS(DT)*60.0.GE.ETIME) ! old until 20100809)
       DO WHILE(ABS(DT)*60.0.GE.ETIME  .AND.  IZNDXT .LE. NLVL) ! new (tk 20100810)
          ICNDX2=ICNDX1
          IZNDXT2=IZNDXT
          ETIME2=ETIME
          Z2=ZPROFT(IZNDXT)
!      WRITE(*,*)"WTT2: ",DABS(WTT(ICNDX1,IZNDXT)),ICNDX1,IZNDXT,Z2,ETIME
!     Reassign Cloud index and vert. index
!     get random number
!          RSEED=30985
          RAND=RAN3(1)

!     check if moves to left or right, do immediately, no change in ETIME
          IF(RAND.LT.PLE(ICNDX1,IZNDXT))THEN
            ICNDX1=ICNDX1-1
          ELSEIF(RAND.LT.PLE(ICNDX1,IZNDXT)+PRI(ICNDX1,IZNDXT))THEN
            ICNDX1=ICNDX1+1
!     check if moves up, do it, change ETIME
          ELSEIF(RAND.LT.PLE(ICNDX1,IZNDXT)+PRI(ICNDX1,IZNDXT)          &
     &      +PUP(ICNDX1,IZNDXT)) THEN
!dwen(20090825)              ETIME=ETIME+DZTA(ICNDX1,IZNDXT)/DABS(WTT(ICNDX1,IZNDXT))
              ETIME=ETIME+DZTA(ICNDX1,IZNDXT)/ABS(WTT(ICNDX1,IZNDXT))
              IZNDXT=IZNDXT+1
!       WRITE(*,*)"MOVED UP"
!dwen(20090825)            IF(ETIME.GT.DABS(DT)*60.0)THEN
!dwen(20090825)             Z2=ZPROFT(IZNDXT2)+DABS(WTT(ICNDX1,IZNDXT2))*              &
!dwen(20090825)     &          (DABS(DT)*60.0-ETIME2)
            IF(ETIME.GT.ABS(DT)*60.0)THEN
             Z2=ZPROFT(IZNDXT2)+ABS(WTT(ICNDX1,IZNDXT2))*              &
     &          (ABS(DT)*60.0-ETIME2)
             IF(Z2.LT.ZPROFM(IZNDXT2)) THEN
               ZX=IZNDXT2
             ELSE
               ZX=IZNDXT
             END IF
            ENDIF

               !moves down, do it, change ETIME
          ELSE
            IZNDXT=IZNDXT-1
!            if(IZNDXT.EQ.0.AND.CFXUP2(1).GT.0.0.and.EFXUP2(1).EQ.0.0
!     :         .and.EFXUP2(2).eq.0.0)then
!DMM-eliminated reference to level 2 in the quick fix, as only level one is applicable
!cg(20110304)            IF(IZNDXT.EQ.0.AND.CFXUP2(1).GT.0.0.and.EFXUP2(1).EQ.0.0)   &
            IF(IZNDXT.EQ.0)   &
    &         THEN

              WRITE(KFJCLMSG,*)"CGRELL: quickfix shallow convection"
                       !if zero prob. to leave cell, can't have updraft
              ICNDX2=2
              ZX=IZNDXT+1
              AREAPRU=AREAA(1,ZX)/AREA
              AREAPRD=AREAA(3,ZX)/AREA
              RETURN
            ELSE
!dwen(20090825)              ETIME=ETIME+DZTA(ICNDX1,IZNDXT)/DABS(WTT(ICNDX1,IZNDXT))
              ETIME=ETIME+DZTA(ICNDX1,IZNDXT)/ABS(WTT(ICNDX1,IZNDXT))
            ENDIF

!       WRITE(*,*)"MOVED DN"
!dwen(20090825)            IF(ETIME.GT.DABS(DT)*60.0)THEN
!dwen(20090825)             Z2=ZPROFT(IZNDXT2)-DABS(WTT(ICNDX1,IZNDXT))*               &
!dwen(20090825)     &          (DABS(DT)*60.0-ETIME2)
            IF(ETIME.GT.ABS(DT)*60.0)THEN
             Z2=ZPROFT(IZNDXT2)-ABS(WTT(ICNDX1,IZNDXT))*               &
     &          (ABS(DT)*60.0-ETIME2)
             IF(Z2.LT.ZPROFM(IZNDXT)) THEN
               ZX=IZNDXT
             ELSE
               ZX=IZNDXT2
             END IF
            ENDIF
          ENDIF

!          IF(ETIME.GT.1000000.0)WRITE(*,*)'cgrell too long',ETIME

!     END OF WHILE LOOP
        END DO
!     ASSIGN VALUES AFTER CONVECTION
      ICNDX2=ICNDX1
      AREAPRU=AREAA(1,ZX)/AREA
      AREAPRD=AREAA(3,ZX)/AREA

!      close(34)

      RETURN
      END
