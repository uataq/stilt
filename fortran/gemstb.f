!###############################################################################
! STBLTY - Compute the maximum permitted time step to maintain computational
! stability according to the maximum wind speed. The routine also computes the 
! vertical mixing coefficient array. Note that the vertical coordinate system is
! defined such that the meteorological data represent values at the box centers.
! The heights of the bottom and top cells are the same as their adjacent cells.
!-------------------------------------------------------------------------------
! LAST REVISED: 28 May 2008 (RRD) - initial version
!-------------------------------------------------------------------------------

SUBROUTINE gemstb 

  USE gemcfg
  USE gemkon    
  USE gemvar   
  USE funits

  IMPLICIT NONE

  REAL*4  :: HSTEP, VSTEP   ! temporary internal time step
  REAL*4  :: FRACZ = 0.50   ! stability limit for vertical mixing
  REAL*4  :: FRACH = 0.75   ! stability limit for horizontal advection  
  REAL*4  :: UMAX  = 150.0  ! maximum permitted horizontal wind speed  

  INTEGER*4 :: i,j,k,nx,ny,nz,np,kgrd
  REAL*4    :: ri,zl,tbar,tsfc,delu,delt,pres
  REAL*4    :: s,t,v,dzc,ustr,tstr,etrm,phim,phih,zmix 

  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  HSTEP = 3600.0  !  default initial time step (in sec)
  VSTEP = 3600.0  

  DO J=1,ny
  DO I=1,nx
  DO K=1,nz

!    delta box sizes between centers
     IF(K.EQ.1)THEN
        DZC=HHH(I,J,K)
     ELSE
        DZC=MAX(DZMIN,HHH(I,J,K)-HHH(I,J,K-1))
     END IF

!    potential temperature and scalar wind speed profile
     PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP        
     POT(K)=TTT(I,J,K)*(1000.0/PRES)**0.286
     VB4(K)=SQRT(UUU(I,J,K)*UUU(I,J,K)+VVV(I,J,K)*VVV(I,J,K))

!    limit maximum scalar wind speed but keep vector orientation
     IF(VB4(K).GT.UMAX)THEN
        UUU(I,J,K)=UUU(I,J,K)*UMAX/VB4(K)
        VVV(I,J,K)=VVV(I,J,K)*UMAX/VB4(K)
     END IF

!    time step due to advection
     IF(j.GT.1.AND.j.LT.ny) &
        HSTEP=MIN(HSTEP,FRACH*REAL(GSX(I,J))*NGP(J)/MIN(MAX(VB4(K),1.0),UMAX))

!    vertical mixing
     IF(K.GT.1)THEN

!       mixing coefficient (ri method) at top of each box
        TBAR=2.0/(POT(K)+POT(K-1))
        DELT=(POT(K)-POT(K-1))/DZC
        DELU=(VB4(K)-VB4(K-1))/DZC
!       RI=SIGN(998.0,DELT)
        IF (DELU.NE.0.0) RI=GRAV*DELT*TBAR/(DELU*DELU)
        RI=MAX(-1.0,MIN(RI,998.0))

!       mixing via nmc-ngm boundary layer method
        ZMIX=VMIX/(2.0+RI)

!       time step due to vertical mixing
        VSTEP=MIN(VSTEP,FRACZ*DZC*DZC/ZMIX)

!       set mixing limits
        IF(PRES.GE.250.0)THEN
           KKK(i,j,k-1)=MAX(0.03, MIN(100.0, ZMIX))
        ELSE
           KKK(i,j,k-1)=MAX(0.03, MIN( 10.0, ZMIX))
        END IF

!    surface stability functions 
     ELSE            
!       recompute richardson number 
        TSFC=T02(I,J)*(1000.0/SFC(i,j))**0.286
        TBAR=(POT(1)+TSFC)/2.0
        DELT=POT(1)-TSFC 
        DELU=(UUU(i,j,1)-U10(i,j))**2+(VVV(i,j,1)-V10(i,j))**2
        DELU=MAX(0.0001,DELU)
        RI=GRAV*DELT*DZC/TBAR/DELU
        RI=MAX(-1.0, MIN(RI,1.0))

        S=LOG(DZC/RGH(i,j)+1.0)
        T=LOG(DZC/(0.1*RGH(i,j)) +1.0)
        V=LOG(RGH(i,j)/(0.1*RGH(i,j)))

        IF(RI.GT.0.0.AND.RI.LT.0.08)THEN
           ZL=(-T+10.0*S*RI+SQRT(T*T-20.0*S*T*RI+20.0*S*S*RI))/(10.0*(1.0-5.0*RI))
        ELSEIF(RI.GE.0.08)THEN
           ZL=(B1*S+B2)*RI*RI+(B3*S-B4*V-B5)*RI
        ELSE
           ZL=RI*(S*S/T-B6)
        END IF
        ZL=SIGN(ZL,RI)

        IF(ZL.GE.0.0)THEN
!          stable surface layer Beljaars-Holtslag
           ETRM=B*EXP(-D*ZL)*(1.0-D*ZL+C)
           PHIM=1.0+(A+ETRM)*ZL
           PHIH=PRN*(1.0+(A*SQRT(1.0+A*B*ZL)+ETRM)*ZL)
        ELSE
!          unstable surface layer Betchov-Yaglom / Kadar-Perepelkin
           PHIM=((1.0+0.625*ZL*ZL)/(1.0-7.5*ZL))**0.3333
           PHIH=0.64*((3.0-2.5*ZL)/(1.0-10.0*ZL+50.0*ZL*ZL))**0.3333
        END IF

!       compute friction terms
        USTR=VONK*DZC*SQRT(DELU)/PHIM/DZC
        TSTR=VONK*DZC*DELT/PHIH/DZC

!       recompute Z/L from friction terms to be consistent
        ZL=DZC*VONK*GRAV*TSTR/(USTR*USTR*TTT(i,j,1))

!       check limits on Z/L (-2 <= z/L <= 15 )
        ZL=MAX(-2.0, MIN(15.0, ZL))
        IF(ZL.EQ.0.0)ZL=0.001

!       compute integral (psi) required for deposition calculations
        IF(ZL.LT.-0.001)THEN
!          unstable
           PSI(i,j)=P1+ZL*(P2+ZL*(P3+ZL*(P4+P5*ZL)))
        ELSEIF(ZL.GE.-0.001.AND.ZL.LT.0.0)THEN
!          neutral
           PSI(i,j)=-2.7283*ZL
        ELSE
!          stable
           PSI(i,j)=-(1.0+A*B*ZL)**1.5-B*(ZL-C/D)*EXP(-D*ZL)-B*C/D+1.0
        END IF

!       other variables to pass through
        FRV(i,j)=USTR
     END IF

  END DO
  KKK(i,j,nz)=0.0

  END DO
  END DO

! time step in even minutes evenly divisible into 60
  DELTA=MAX(1.0,MIN(HSTEP,VSTEP)/60.0)
  DO WHILE (MOD(60,DELTA).NE.0)
     DELTA=DELTA-1
  END DO
  DELTA=MAX(1,DELTA)    

  WRITE(KF21,*)' NOTICE gemstb: minimum H,V time steps (s) - ', &
                 NINT(hstep),NINT(vstep),'  Final (min):',DELTA

END SUBROUTINE gemstb
