!-----------------------------------------------------------------------------
! This file contains subroutines to compute ozone formation using the GRS
! chemistry model.  The routines were developed by Martin Cope, Victoria EPA,
! and later modified by G.D. Hess, BMRC, BoM, Melbourne, Australia.
!-----------------------------------------------------------------------------
 
SUBROUTINE GCHEMS(DIFTM,CON,CON1,K,active1,active2,cisop)

  USE funits

! number of species and reactions
  PARAMETER (NOSPEC=6,NOREAC=7)

! rrd- added isoprene (cisop) and its activity (active2)
! SUBROUTINE INTEGRATES THE CBM-XR CHEMISTRY WITH A SULPHUR
! OXIDATION EXTENSION FOR ONE CELL FOR ONE TIME STEP
 
! REACTION RATES AND TOTAL FORMATION/LOSS TERMS
 
  REAL K(NOREAC)
  REAL A(NOREAC),AP(NOREAC),B(NOREAC),AOB(NOREAC),R(NOREAC)
  REAL ZERO(NOREAC)
 
! CONCENTRATION ARRAYS. CON= INPUT ARRAY
!                       CON1=PREDICTED CONCENTRATION
 
  DIMENSION CON(NOREAC),CON1(NOREAC)
  REAL LMBDA

  DATA ZERO /7*1.0E-12/

  SAVE AOB,A,AP,B,R,ZERO

  elpstm=0.

! ENTRY FOR INTEGRATION STAGE. (PECE MODE). CONTINUE UNTIL
! FULL TIME STEP COMPLETED
! ALSO ENTRY POINT FOR NEW PREDICTION
 
  DO WHILE (elpstm.LT.diftm)

     dt=diftm-elpstm
 
!    CALCUATE TIME RATE OF CHANGE OF ALL SPECIES.
!    CHECK D[NO]/DT FOR CONVERGENCE CRITERIA
 
     CALL ABCALC(NOSPEC,A,B,K,R,CON,active1,active2,cisop)

!    IF(abs(A(2)/CON(2)-B(2)).GT.0.005)THEN
        LMBDA=0.001
!    ELSE
!       LMBDA=0.01
!    ENDIF
 
!    WORK OUT MAXIMUM TIME STEP. ONLY CONSIDER SPECIES
!    WITH 'APPRECIABLE' CONCENTRATIONS (NO+NO2+O3)
!    INTEGRATE FORWARD.
 
     CONSUM=CON(1)+CON(2)+CON(3)
     CONSUM=CONSUM*0.001

     DO I=1,NOSPEC
        IF(B(I).GT.0)THEN
           AOB(I)=A(I)/B(I)
           IF(CON(I).GT.CONSUM)THEN
              if(con(i)-aob(i).ne.0.0)then
                 F=LMBDA*AMAX1(AOB(I),CON(I))/ABS(CON(I)-AOB(I))
              end if
              IF(F.LT.1.)DT=AMIN1(DT,-ALOG(1.0-F)/B(I))
           ENDIF
        ENDIF
     END DO           

     if(dt.le.0)THEN
        write(KF21,*)'*ERROR*: GRS chems subroutine dt<0'
        dt=0.001
     END IF

     CALL INTGRT(NOSPEC,A,B,AOB,DT,CON,CON1,ZERO)
 
!    CORRECTED SOLUTION
 
     CALL ABCALC(NOSPEC,AP,B,K,R,CON1,active1,active2,cisop)
     DO I=1,NOSPEC
        A(I)=(A(I)+AP(I))*0.5
        IF(B(I).GT.0)AOB(I)=A(I)/B(I)
     END DO
 
     CALL INTGRT(NOSPEC,A,B,AOB,DT,CON,CON1,ZERO)
     DO I=1,NOSPEC
        CON(I)=CON1(I)
     END DO
     elpstm=elpstm+DT

  END DO
 
  END SUBROUTINE gchems

!=======================================================

SUBROUTINE ABCALC(NOSPEC,A,B,K,R,CON,active1,active2,cisop)
 
! ROUTINE CALCULATES FORMATION AND REACTION
! RATES FOR ALL SPECIES
 
  REAL NO,NO2
  REAL A(7),B(7),K(7),R(7)
  REAL CONSP(200),con(7)
 
! WORKING ARRAY EQUIVALENCED TO SPECIES MENOMICS FOR CONVIENCE
 
  EQUIVALENCE (consp(1),no2),(consp(2),no),(consp(3),o3),                  &
              (consp(4),sgn),(consp(5),sngn),(consp(6),roc)

! ASSIGN CURRENT ARRAY TO WORKING ARRAY

  DO L=1,NOSPEC
     CONSP(L)=CON(L)
  END DO      

! EVALUATE RATE FUNCTIONS

  R(1)=K(1)*ROC+ active2*k(1)*cisop/active1
  R(3)=K(3)*NO2
  R(4)=K(4)*NO*O3
 
! calculate concentrations of RSP
 
  rx2=k(2)*no
  rx6=k(6)*no2
  rx7=k(7)*no2
  if(k(3).gt.0)then
     bquad=rx2+rx6+rx7
     root2=bquad*bquad+4.*k(5)*r(1)
     if(root2.ge.0.0)then
        rsp=-bquad+sqrt(root2)
        rsp=rsp*0.5/k(5)
     else
        rsp=0.
     end if
  else
     rsp=0.
  endif

  R(2)=rx2*RSP
  R(5)=k(5)*RSP*RSP
  R(6)=Rx6*RSP
  R(7)=Rx7*RSP

! FORMATION and LOSS TERMS

! NO2
      a(1)=R(2)+R(4)
      b(1)=R(3)+R(6)+R(7)
! NO
      a(2)=R(3)
      b(2)=R(2)+R(4)
! O3
      a(3)=R(3)
      b(3)=R(4)
! SGN
      a(4)=R(6)
      b(4)=0
! SNGN
      a(5)=R(7)
      b(5)=0.
! roc
      a(6)=0.
      b(6)=0.

  DO L=1,NOSPEC
     IF(CONSP(L).GT.0.0)B(L)=B(L)/CONSP(L)
  END DO      

END SUBROUTINE ABCALC

!========================================================

SUBROUTINE INTGRT(NOSPEC,A,B,AOB,DT,CON,CON1,ZERO)
 
! ROUTINE PERFORMS EXPONENTIAL PREDICTOR STEPCORRECTION STEP
! FOR CHEMISTRY.
 
  REAL A(7),B(7),AOB(7),CON(7),ZERO(7),CON1(7)
 
  DO I=1,NOSPEC
     T = B(I)*DT
 
!    CASE 1 B*DT<<1
     IF( T.LE.1.E-06)THEN
        CON1(I)=A(I)*DT+CON(I)
 
!    CASE 2 SMALL<B*DT<LARGE
     ELSEIF(T.LE.75.)THEN
        CON1(I)=AOB(I)+(CON(I)-AOB(I))*EXP(-T)
 
!    CASE 3 B*DT>>1
     ELSE
        CON1(I)=AOB(I)
     ENDIF

     CON1(I)=AMAX1(ZERO(I),CON1(I))

  END DO 

END SUBROUTINE INTGRT

!=======================================================

SUBROUTINE GRATES(ZENITH,TEMP,K,KP,ACTIVITY)
 
! ROUTINE CALCULATES ALL PHOTOLYTIC AND MOLAR REACTION
! RATES FOR THE GRS MECHANISM
 
  REAL K(7),KP
  PARAMETER (T1=0.00316)
 
! COMPUTE RATE CONSTANTS
! SET AND ADJUST PHOTOLYSIS RATES
 
  TEMPI=1.0/TEMP
! K(4)=9.24E+05*TEMPI*EXP(-1450.0*TEMPI)
  K(4)=2643.0*EXP(-1370.0*TEMPI)

  IF(ZENITH.LT.90.0)THEN
!    K(1)=KP*ACTIVITY*EXP(-1000.0*4.7*(TEMPI-T1))
     K(1)=KP*10000.0*EXP(-4710.0*TEMPI)

!    K(2)=(3.58E+06)*TEMPI
     K(2)=5482.0*EXP(242.0*TEMPI)

     K(3)=KP
     K(5)=10000.
     K(6)=125.
     K(7)=K(6)
 
! IF NIGHTIME THEN SET PHOTOLYSIS RATES TO ZERO
  ELSE
     K(1)=0.
     K(2)=0.
     K(3)=0.
     K(5)=0.
     K(6)=0.
     K(7)=0.
  ENDIF

END SUBROUTINE grates
