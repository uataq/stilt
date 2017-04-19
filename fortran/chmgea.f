!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CHMKIN           CHEMical Kinetics   
!   PRGMMR:    ARIEL STEIN      ORG:            DATE:97-07-17
!
! ABSTRACT:  THIS CODE MODIFIED AT THE DEPARTMENT OF METEOROLOGY
!            PENN STATE UNIVERSITY
!                                                                        
!*********************************************************************** 
!                                                                        
!     CHEMK IS A FORTRAN PROGRAM WHICH, WHEN GIVEN A PREDETERMINED       
!     KINETIC MECHANISM, COMPUTES THE CONCENTRATIONS OF THE VARIOUS      
!     REACTANTS IN TIME. IT WAS WRITTEN BY GARY Z. WHITTEN OF SYSTEMS    
!     APPLICATIONS INCORPORATED IN SAN RAFAEL, CA AND INCORPORATES       
!     THE GEAR INTEGRATION PACKAGE OF A. HINDMARSH OF LAWRENCE LIVERMORE 
!     LABORATORIES IN LIVERMORE, CA.  THE PRINTER-PLOT ROUTINE WAS       
!     PROVIDED BY DAVID C. WHITNEY OF SAI AND THE PROGRAM WAS CONVERTED  
!     ANSI FORTRAN BY JIM MEYER OF SAI. THE  COMPUTER CODE THAT FOLLOWS  
!     IS EXECUTABLE ON ANY IBM OR CDC MACHINE WITH A FORTRAN IV COMPILER 
!                                                                        
!     TWO SEPARATE SETS OF DATA ARE NEEDED TO EXECUTE THE PROGRAM. ONE   
!     SET CONTROLS THE COMPUTATIONAL PART OF THE PROGRAM AND IS          
!     NECESSARY TO ESTABLISH THE INTEGRATION ROUTINE.  THE SECOND SET    
!     CONTROLS THE PLOTTER OUTPUT AND PROVIDES THE PARAMETERS NECESSARY  
!     SPECIFY THE  OUTPUT FORMAT.  THE READER IS REFERED TO THE PROGRAM  
!     FOR FURTHER INFORMATION.                                           
!                                                                        
!                                    GARY Z. WHITTEN                     
!                                    SYSTEMS APPLICATIONS INCORPORATED   
!                                    950 NORTHGATE DRIVE                 
!                                    SAN RAFAEL, CALIFORNIA  94903       
!                                    (415) 472-4011                      
!									 
!******************************************************************* 
!
!
!
!     1/21/90
!     Converted to double precision and added "spread-sheet" type 
!       output file called "CHEMK.txt".
!                                            John R. Barker (Michigan)
!
!*********************************************************************** 
!  
! USAGE:  CALL CHMKIN(DIRT,DDT,CBSUM,SPCLTEMP,NUMTYP,CASUM,KLM)  
!   INPUT ARGUMENT LIST:
!     DDT	- integration time step (min)
!     CBSUM- initial concentration at the grid
!     SPCLTEMP	- average local temperature (K)
!     NUMTYP    - int	number of pollutant types
!     KLM- skip reading data file index
!   OUTPUT ARGUMENT LIST:
!     casum	- final concentration at the grid 
!   INPUT FILES:
!     CHEMK.DAT
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  DIOXIN
!
!$$$!                                                                        

      SUBROUTINE CHMGEA(DIRT,DDT,CBSUM,SPCLTEMP,NUMTYP,CASUM,KLM,RRK1,RKJ)                                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
  INCLUDE 'DEFCONC.INC'

  TYPE(pset),INTENT(IN)    :: dirt(:)           ! for each pollutant type 
  REAL*8,    INTENT(IN)    :: ddt               ! time step (min)
  REAL*8,    INTENT(IN)    :: cbsum(:)          ! concentration before chemistry
  REAL*8,    INTENT(IN)    :: spcltemp          ! temperature
  INTEGER,   INTENT(IN)    :: numtyp            ! number of pollutant types
  REAL*8,    INTENT(OUT)   :: casum(:)          ! concentration after chemistry
  INTEGER,   INTENT(INOUT) :: klm               ! only with klm=0 reads the chem mech
  REAL*8,    INTENT(IN)    :: rrk1              ! NO2 photolysis rate constant
  REAL*8,    INTENT(INOUT) :: rkj(:,:)          ! ratio of photolysis constants
 
      COMMON /DATA/ A(200),S(200),TEMP,ERR,START,STOPP,PC,SIG(91),     &
        R(200),BK,SG,DILUT,NP,NR,KR(200,7),IP(91),ITYPE(200),ITITLE(7)
      COMMON /NAMES/ SPECIS(91),REACT(91),NS                             
      COMMON /FRPLOT/ SAVCON(90,80),SAVTIM(80),JGRID(89,40),NIT(3),NT    
      COMMON /ALPHA/ IGO(4),IBLANK,MBLANK,JINTER                         
      COMMON /APLOT/ JVERT(52,2),JBLANK,JSTAR,JPLUS,JBAR                 
      COMMON /GEAR1/ T,GUESS,HMIN,HMAX,EPS1,UROUND,NC,MF1,KFLAG1,JSTART  
      COMMON /GEAR2/ YMAX(100)                                           
      COMMON /INOUT/ IN,IOUT,ITAPE, ITEXT                                
      COMMON /HEAT/ CV,Q(200),SC(200,7),ISC(200,3),ITEMP                 
      COMMON /SPARSE/ IA(91),JA(1000)                                    
      COMMON /STORE/ AST(35),IPL(7),TEMEND,NTEMP,TMI,NPHOT,PHI,IL,NFRST, &
      IPH(30),QM(100),PM(100),PSTOP                                      
      COMMON /PHOTR/ RDAT(80),RTIM(80),RR1,IN10,IPP 
      DIMENSION C(91),IRS(200,7),KRS(7),SC1(7) ,CAFTER(91),   &
                   DDEPTD(91),DDEPTW(91),AA(200),SS(200)
!      REAL*8 RRKJ(MAXTYP,2)	             
!      REAL*8 DDT,SPCLTEMP,CBSUM,CASUM
       REAL*8 CAFTER                          
      INTEGER SPECIS,PHI,TMI,AST,RID(91),KTIME
	CHARACTER*4 REACT,WREACT,PREACT
	CHARACTER*1 REATYPE(200)    
                                        
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,   &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,    &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
!                                                                        
!     ALPHAMERIC DATA ARE PREASSIGNED IN BLOCK DATA ALPHA1               
!                                                                        
      DATA NPLOT/4HSAVE/
      DATA IBLK/1H /
      DATA IN2,IN3,IN5,IN6,IN7,IN8,IN9/7*1/                              
      DATA XN4,XN6,XN9/3*1.0/                                            
      DATA KTIME/0/
!rrd  SAVE KTIME
	
!                                                                        
!     INITIAL PARAMETERS                                                 
!                                                                        
      IF(KLM.NE.0)GOTO 90
      
!     reads to get the character
       
      IF(KTIME.EQ.0)THEN
       DO I=1,NUMTYP
            READ(DIRT(I)%IDENT,'(A4)')RID(I)
         END DO
!         OPEN(15,FILE='IDENT.DAT')
!         DO I=1,NUMTYP
!            WRITE(15,221)DIRT(I)%IDENT
!         END DO
!         CLOSE(15)

!         OPEN(15,FILE='IDENT.DAT')
!         DO I=1,NUMTYP
!            READ(15,221)RID(I)
!         END DO
!         CLOSE(15)
         KTIME=1
      ENDIF
      
      IN=3                                                                   
      OPEN(IN,FILE='CHEMK.DAT',STATUS='OLD')
!      OPEN(IOUT,FILE='CHEMK.OUT',STATUS='UNKNOWN')
!      OPEN(ITEXT,FILE='CHEMK.txt',STATUS='UNKNOWN')
!                                                                        
!     CLEAR PATTERN MATRIX AND SET THE FIRST ELEMENTS                    
!     ALSO SET STOICHIOMETRIC COEFFICIENTS EQUAL TO 1                    
!                                                                        
      DO 5 J=1,200                                                 
      DO 5 K=1,7                                                   
      IRS(J,K)=IBLANK                                              

      SC(J,K)=1.                                                   
    5 KR(J,K)=0                                                    
      DO 10 J=1,200                                                
      A(J)=0.                                                      
      R(J)=0.                                                      
   10 S(J)=0.                                                       
      NR=0                                                          
      NS=0                                                               
!                                                                        
!     NX = LAST NUMBER OF THE REACTION SEQUENCE                          
!     NPLOT = PLOT OPTION                                                
!     DILUT= DILUTION FACTOR                                             
!                                                                        
      READ (IN,185) NX,NEPA,NPLIT,DILUT,NPHOT,PHI,IL,NTEMP,TMI,IPP       
      GO TO 20                                                            
   15 NS=NS-1                                                             
      READ (IN,185) NX,IN2,IN3,XN4,IN5,IN6,IN7,IN8,IN9,IN10               
      IF (IABS(IN2).NE.0) NEPA=IN2                                        
      NPLIT=IN3                                                          
      IF (ABS(XN4).NE.0) DILUT=XN4                                       
      IF (IABS(IN5).NE.0) NPHOT=IN5                                      
      IF (IABS(IN6).NE.0) PHI=IN6                                        
      IF (IABS(IN7).NE.0) IL=IN7                                         
      IF (IABS(IN8).NE.0) NTEMP=IN8                                      
      IF (IABS(IN9).NE.0) TMI=IN9                                        
      IF (IABS(IN10).NE.0) IPP=IN10                                      
   20 IF (IL.LE.0.OR.IN7.LE.0) GO TO 25                                  
      READ (IN,200) (IPH(I),I=1,IL)                                      
   25 IF (NPHOT.LE.0.OR.IN5.LE.0) GO TO 30                               
      NPHOT=NPHOT+1                                                      
      READ (IN,190) (PM(I),I=1,NPHOT)                                    
      PM(NPHOT+1)=2.0*PM(NPHOT)-PM(NPHOT-1)                              
      PM(NPHOT+2)=3.0*PM(NPHOT)-2.0*PM(NPHOT-1)                          
      PSTOP=DFLOAT((NPHOT-1)*PHI)                                        
   30 IF (NTEMP.LE.0.OR.IN8.LE.0) GO TO 35                               
      NTEMP=NTEMP+1                                                      
      READ (IN,190) (QM(I),I=1,NTEMP)                                    
      QM(NTEMP+1)=2.*QM(NTEMP)-QM(NTEMP-1)                               
      QM(NTEMP+2)=3.0*QM(NTEMP)-2.0*QM(NTEMP-1)                          
      TEMEND=DFLOAT((NTEMP-1)*TMI)                                       
!   35 IF (NX.GT.0) WRITE (IOUT,130)                                      
   35 IF (NX.LE.0) NS=NS+1                                               
      IF (NX.LE.0) WRITE(*,*)'ERROR IN CHEMK.DAT, # REACTIONS'                                             
!                                                                        
!     REACTION INPUT DATA                                                
!                                                                        
   40 IF (NEPA.LE.0) READ (IN,135) J,(IRS(J,I),I=1,7),A(J),S(J)          
      IF (NEPA.LE.0) GO TO 55                                            
      READ (IN,220) (KRS(I),I=1,3),J,SC1(1),JJ,KRS(4),(SC1(LL-3),KRS(LL)  &
      ,LL=5,6),AA(JJ*100+J),SS(JJ*100+J),REATYPE(JJ*100+J)                                                
      IF (JJ.GT.0) J=JJ*100+J                                            
      DO 45 II=1,6                                                       
   45 IRS(J,II)=KRS(II)





!	 determines the kind of reaction (ie, A=Arrhenius, P=Photolysis,
!                                         R=double phot,W=Aqueous,
!                                             F=photolysis in water
!                                             D=dry deposition
!                                              L=wet deposition)

      
      IF(REATYPE(J).EQ.'A')THEN                                               
      A(J)=AA(J)                                                           
      S(J)=SS(J)
	ELSEIF(REATYPE(J).EQ.'P')THEN
!##PREACT=KRS(1) 
         
	 DO KC=1,NUMTYP
	 
         IF(KRS(1).EQ.RID(KC))THEN
!	 WRITE(*,*)KRS(1), RID(KC),rkj(kc,1),rrk1
	A(J)=AA(J)*RRK1*rkj(KC,1)
	S(J)=SS(J)
	
	END IF
	 END DO
	ELSEIF(REATYPE(J).EQ.'R')THEN
!##PREACT=KRS(1) 
	 DO KC=1,NUMTYP
	IF(KRS(1).EQ.RID(KC))THEN
  	A(J)=AA(J)*RRK1*rkj(KC,2)
 	S(J)=SS(J)
	END IF
	 END DO

	ELSEIF(REATYPE(J).EQ.'W')THEN
!	A(J)=AA(J)*RRLWC*(SPCLTEMP/DBLE(298.0))
	S(J)=SS(J)
	ELSEIF(REATYPE(J).EQ.'F')THEN
	A(J)=AA(J)*RRK1
	S(J)=DBLE(0.0)
	ELSEIF(REATYPE(J).EQ.'D')THEN
	
!##WREACT=KRS(1) 
	DO KC=1,NUMTYP
	IF(KRS(1).EQ.RID(KC))THEN
!       A(J)=DDEPTD(KC)
        A(J)=0.0
	S(J)=DBLE(0.0)
	END IF
	END DO
	ELSEIF(REATYPE(J).EQ.'L')THEN
	
!##WREACT=KRS(1) 
	  DO KC=1,NUMTYP
	  IF(KRS(1).EQ.RID(KC))THEN
!         A(J)=DDEPTW(KC)
          A(J)=0.0
	  S(J)=DBLE(0.0)
	  END IF
	  END DO

      ELSEIF(REATYPE(J).EQ.'S')THEN
      COSA=DEXP(-4700.*((1/SPCLTEMP)-(1/316.)))
	A(J)=AA(J)*COSA*RRK1
	S(J)=DBLE(0.0) 
      END IF 
	                                                
      DO 50 II=4,6                                                       
   50 SC(J,II)=SC1(II-3)                                                 
   55 DO 60 I=1,7                                                        
   60 IF (SC(J,I).LE.0.) SC(J,I)=1.                                      
      DO 75 K=1,7                                                        
      NL=K*5                                                             
      NF=NL-4                                                            
      DO 65 LK=NF,NL                                                     
   65 AST(LK)=IBLK                                                       
      IPL(K)=IBLK                                                        
      IF (K.EQ.1.OR.K.EQ.4) GO TO 70                                     
      IF (IRS(J,K).NE.IBLANK) IPL(K-1)=JPLUS                             
   70 IF (IRS(J,K).NE.IBLANK) CALL VALU (SC(J,K),K,NF,NL)                
   75 CONTINUE                                                           
!      WRITE (IOUT,140) J,(AST(I),I=1,5),IRS(J,1),IPL(1),(AST(I),I=6,10), 
!     &IRS(J,2),IPL(2),(AST(I),I=11,15),IRS(J,3),IPL(3),(AST(I),I=16,20), 
!     &IRS(J,4),IPL(4),(AST(I),I=21,25),IRS(J,5),IPL(5),(AST(I),I=26,30), 
!     &IRS(J,6),IPL(6),(AST(I),I=31,35),IRS(J,7),IPL(7),A(J),S(J)         
      KR(J,1)=100                                                        
      IF (J.GT.NR) NR=J                                                  
      IF (J.NE.NX) GO TO 40  
      CLOSE(UNIT=IN)                                            
!                                                                        
!     ESTABLISH REACTION MATRIX AND SET UP SPARSE JACOBIAN VECTORS       
!                                                                        
      CALL MATRX (C,IRS)                                                 
      CALL SPARS (IA,JA,NS-2)                                            

   90   INX=1
      IF(KLM.NE.0)THEN 
   	DO J=1,NX
   	   IF(REATYPE(J).EQ.'P')THEN
!##PREACT=IRS(J,1)
 	     DO KC=1,NUMTYP
		IF(IRS(J,1).EQ.RID(KC))THEN
!	WRITE(*,*)'2nd round:',IRS(j,1), RID(KC),rkj(kc,1),rrk1
		A(J)=AA(J)*RRK1*rkj(KC,1)
		S(J)=SS(J)
	     END IF
		END DO
	ELSEIF(REATYPE(J).EQ.'R')THEN
!##PREACT=IRS(J,1) 
		DO KC=1,NUMTYP
		IF(IRS(J,1).EQ.RID(KC))THEN
  		A(J)=AA(J)*RRK1*rkj(KC,2)
 		S(J)=SS(J)
	   END IF
	 END DO

	   ELSEIF(REATYPE(J).EQ.'W')THEN
!	    A(J)=AA(J)*RRLWC*(SPCLTEMP/DBLE(298.0))
	    S(J)=SS(J)
	   ELSEIF(REATYPE(J).EQ.'F')THEN
	    A(J)=AA(J)*RRK1
	    S(J)=DBLE(0.0)
	   ELSEIF(REATYPE(J).EQ.'D')THEN
	     
!##WREACT=IRS(J,1) 
	     DO KC=1,NUMTYP
	      IF(IRS(J,1).EQ.RID(KC))THEN
!               A(J)=DDEPTD(KC)
                A(J)=0.0
	        S(J)=DBLE(0.0)
	      END IF
	     END DO
!       wet deposition
	ELSEIF(REATYPE(J).EQ.'L')THEN
	     
!##WREACT=IRS(J,1) 
	     DO KC=1,NUMTYP
	      IF(IRS(J,1).EQ.RID(KC))THEN
!               A(J)=DDEPTW(KC)
                A(J)=0.0
	        S(J)=DBLE(0.0)
	      END IF
	     END DO

	    ELSEIF(REATYPE(J).EQ.'S')THEN
          COSA=DEXP(-4700.*((1/SPCLTEMP)-(1/316.)))
	    A(J)=AA(J)*COSA*RRK1
  		S(J)=DBLE(0.0)
         END IF 
	END DO
	END IF
   
   	 
                                                             
!                                                                        
!     TITLE CARD AND OPTIONS                                             
!                                                                        
!   80 READ (IN,195) (ITITLE(I),I=1,7)                                    
!      NFL=0                                                              
!      IF (ITITLE(1).EQ.IGO(1)) GO TO 15                                  
!      IF (ITITLE(1).EQ.IGO(2)) GO TO 115                                 
!      IF (ITITLE(1).EQ.IGO(3)) GO TO 85                                  
!      IF (ITITLE(1).EQ.IBLANK) THEN                                      
!	 CLOSE (UNIT=IOUT)						 
!         CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=IN)
!	 RETURN
!	 END IF
!      GO TO 90                                                           
!                                                                        
!     CALL XPLOT                                                          
!                                                                        
!   85 READ (IN,180) (NIT(I),I=1,3),IDT,KALCMP                            
!      NT=NT-2                                                            
!      CALL XPLOT (NIT,NT,NS,SPECIS,SAVTIM,SAVCON,IDT,KALCMP,JGRID)        
!      GO TO 80                                                           
!                                                                        
!     OPTIONS CARD                                                       
!                                                                        
!   90 READ (IN,155) N,NPRNT,TPRNT,TSTEP,TFACT  
      N=NUMTYP + 2
      NPRNT=INT(DDT) + 1 
      TPRNT=DDT + 1.E-10
       TSTEP=DBLE(0.)
       TFACT=DBLE(0.)                             
!                                                                        
!     TIME STEP SKIP OPTION                                              
!                                                                        
      IF (NPRNT*NPRNT.EQ.0) NPRNT=100000                                 
!                                                                        
!     TIME STEP LENGTH OPTION                                            
!                                                                        
      IF (TPRNT*TPRNT.EQ.0.) TPRNT=1.E10                                 
!                                                                        
!     CONCENTRATION OF SPECIES INITIALLY PRESENT                         
!                                                                        
!      READ (IN,125) (REACT(I),I=1,N)                                     
!      READ (IN,190) (C(I),I=1,N)                                        
      DO I=1,NUMTYP
      REACT(I)=DIRT(I)%IDENT
!	IF(DIRT(I).IDENT.EQ.'CO2')THEN
!	C(I)=DBLE(356.0)
!	ELSE
      C(I)=CBSUM(I)
!	END IF
!      IF(C(I).NE.0.0)WRITE(*,*)'ENTRADA',   &
!     DIRT(I)%IDENT,RID(I),C(I)
!     :  C(I)
      END DO

      REACT(NUMTYP + 1)='N2'
      REACT(NUMTYP + 2)='O2'
      C(NUMTYP + 1)=DBLE(7.81E+05)
      C(NUMTYP + 2)=DBLE(2.09E+05)

!                                                                        
!     STARTING AND ENDING INTEGRATION TIMES                              
!                                                                        
!      READ (IN,175) START,STOPP,TEMP,ERR
       START=DBLE(0.)
      STOPP=DDT
      TEMP=SPCLTEMP
      ERR=DBLE(5.0E-04) 
!      READ(IN,*)ERR
                                 
      IF (START*START.EQ.0.) START=0.                                    
!                                                                        
!     SPECIFICATION OF THE TEMPERATURE IF UNSPECIFIED IN INPUT           
!                                                                        
      IF (TEMP.LE.0.) TEMP=298.                                          
!                                                                        
!     SPECIFICATION OF THE ERROR BOUND IF UNSPECIFIED IN INPUT           
!                                                                        
      IF (ERR*ERR.EQ.0.) ERR=1.E-2                                       
!                                                                        
!     OUTPUT OF THE INITIAL CONDITIONS                                   
!                                                                        
!      WRITE (IOUT,145) (ITITLE(I),I=1,7),(REACT(I),I=1,N)                
!   95 WRITE (IOUT,150) (C(I),I=1,N)
                                        
!                                                                          
!     COMPUTE THE NET RATES OF REACTION                                  
!                                                                        
      IF (NTEMP.GT.0) TEMP=QM(1)                                         
      IF (ITITLE(1).EQ.IGO(2).AND.NTEMP.GT.0) TEMP=TEMOLD                
      CALL RATES (C,N)                                                   
!                                                                        
!     OUTPUT OF THE TEMPERATURE AND ERROR BOUND                          
!                                                                        
!      WRITE (IOUT,160) TEMP,ERR                                          
!                                                                        
!     OUTPUT OF THE DILUTION FACTOR                                      
!                                                                        
      IF (DILUT*DILUT.EQ.0.) DILUT=0.                                    
!      IF (DILUT.NE.0.) WRITE (IOUT,170) DILUT                            
!                                                                        
!     SET LIMITS FOR TIMED OUTPUTS                                       
!                                                                        
      IF (TSTEP.NE.0.) YMAX(1)=TSTEP                                     
      IF (TSTEP.NE.0.) PC=TSTEP                                          
      IF (TSTEP*TSTEP.EQ.0.) YMAX(1)=1.E10                               
      IF (TSTEP*TSTEP.EQ.0.) PC=1.E10                                    
      IF (TFACT.NE.0.) YMAX(2)=TFACT                                     
      IF (TFACT*TFACT.EQ.0.) YMAX(2)=1.                                  
!                                                                        
!     RATE CONSTANTS                                                     
!                                                                        
!      WRITE (IOUT,165) (R(I),I=1,NR)                                     
!                                                                        
!     WRITE INITIAL CONDITIONS OF THE CELL                               
!                                                                        
      IF (NPHOT.LE.0) GO TO 105                                          
!      WRITE (IOUT,215) (IPH(I),I=1,IL)                                   
      DO 100 I=1,IL                                                      
      J=IPH(I)                                                           
      YMAX(10+I)=R(J)                                                    
  100 CONTINUE                                                           
!      WRITE (IOUT,225) (YMAX(10+I),I=1,IL)                               
      NPT=NPHOT                                                          
!      WRITE (IOUT,205) PHI,(PM(I),I=1,NPT)                               
  105 IF (NTEMP.LE.0) GO TO 110                                          
      NPT=NTEMP                                                          
!      WRITE (IOUT,210) TMI,(QM(I),I=1,NPT)                               
  110 GUESS=1.E-10                                                       
      T=START                                                            
      IF (NPLOT.EQ.IBLANK) YMAX(3)=0.                                    
      IF (NPLOT.NE.IBLANK) YMAX(3)=1.                                    
      IF (NPLIT.NE.IBLANK) YMAX(4)=1.                                    
      IF (NPLIT.EQ.IBLANK) YMAX(4)=0.                                    
!                                                                        
!     INITIALIZE PARAMETERS                                              
!                                                                        
      CALL YFIX (NS,TPRNT,C,NPRNT,INX,CAFTER,DDT,RID,numtyp)                                   
      INX=4                                                              
      TEMOLD=TEMP
      KLM=KLM+1 
	DO III=1,numtyp
	CASUM(III)=CAFTER(III)
	END DO                                                       
      RETURN
!      GO TO 80                                                           
!                                                                        
!     CONTINUATION OF DATA                                               
!                                                                        
!  115 READ (IN,175) START,STOPP,TEMP,ERR                                 
!      IF (TEMP.LE.0.) TEMP=298.                                          
!      IF (TEMOLD.NE.298..AND.TEMOLD.NE.0.) TEMP=TEMOLD                   
!      IF (ERR*ERR.EQ.0.) ERR=1.E-2                                       
!      DO 120 I=1,NS                                                      
!  120 REACT(I)=SPECIS(I)                                                 
!                                                                        
!     INITIAL CONDITIONS                                                 
!                                                                        
!      WRITE (IOUT,145) (ITITLE(I),I=1,7),(REACT(I),I=1,NS)               
!      N=NS                                                               
!      GO TO 95                                                           
!                                                                        
!                                                                        
  125 FORMAT (7(A4,6X))                                                  
  130 FORMAT (1H1,14H THE REACTIONS,86X,' RATE CONSTANT',2X,' ACT.  ENERGY (K)')
  135 FORMAT (I3,2X,7(A4,1X),2F10.0)                                     
  140 FORMAT (/,2X,I3,2X,3(5A1,1X,A4,2X,A1),1H=,2X,4(5A1,1X,A4,2X,A1),1P   &
      E11.3,2X,E13.3)                                                    
  145 FORMAT (1H1,30X,7A4,//,23H INITIAL CONCENTRATION ,//(10X,10(4X,A4,   &
      4X)))                                                              
  150 FORMAT (/(8X,1P10E12.3))                                           
  155 FORMAT (I3,17X,I6,4X,3F10.0)                                       
  160 FORMAT (/,34H THE TEMPERATURE OF THE CELL WAS =,1PE9.2,&
      '   AND THE ERROR TOLERANCE =',E9.2)
  165 FORMAT (/,29H THE RATE CONSTANTS USED WERE,/,(/,8X,1P10E12.3))
  170 FORMAT (/,30H THE OVERALL DILUTION RATE WAS,1PE9.2)
  175 FORMAT (8F10.0) 
  180 FORMAT (3A4,43X,2I5) 
  185 FORMAT (I3,2X,I5,A4,6X,F10.2,6I5)  
  190 FORMAT (7F10.3) 
  195 FORMAT (7A4)  
  200 FORMAT (7I10) 
  205 FORMAT (/44H0THE NO2 PHOTOLYSIS NUMBERS IN INTERVALS OF ,I3,'   TIME UNITS ARE',/,(1H0,1P10E13.3))
  210 FORMAT (/34H0THE TEMPERATURES IN INTERVALS OF ,I3,'   TIME UNITS ARE'/(1H0,1P10E13.3))
  215 FORMAT (/29H0THE PHOTOLYSIS REACTIONS ARE,/(1H0,9I13)) 
  220 FORMAT (2(A4,2X),A4,I2,F5.2,I1,A4,2X,2(F5.2,1X,A4,2X),F9.2,1X,F7.2,A1) 
  221 FORMAT (A4)                                                                   
  225 FORMAT ('   THE PHOTOLYSIS RATIOS OF EACH OF THE REACTIONS ABOVE ARE'/(1H0,4X,1P9E13.3))
      END                                                                


      SUBROUTINE YFIX (N0,TLAST,C,NQ,INDEX,CAFTER,DDT,RID,numtyp)                              
!                                                                        
!*********************************************************************** 
!                                                                        
!     THIS IS THE NORMAL OUTPUT ROUTINE                                  
!                                                                        
!*********************************************************************** 
!                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'DEFCONC.INC' 
      SAVE
      COMMON /SPARSE/ IA(91),JA(1000)                                    
      COMMON /DATA/ A(200),S(200),TEMP,ERR,START,STOPP,PC,SIG(91),    &
        R(200),BK,SG,DILUT,NP,NR,KR(200,7),IP(91),ITYPE(200),ITITLE(7)
      COMMON /NAMES/ SPECIS(91),REACT(91),NS                             
      COMMON /FRPLOT/ SAVCON(90,80),SAVTIM(80),JGRID(89,40),NIT(3),NT    
      COMMON /GEAR1/ D,H,DUM(4),IDUM(4)                                  
      COMMON /GEAR2/ YMAX(100)                                           
      COMMON /INOUT/ INPP,IOUT,ITAPE, ITEXT                              
      COMMON /HEAT/ CV,Q(200),SC(200,7),ISC(200,3),ITEMP                 
      COMMON /STORE/ AST(35),IPL(7),TEMEND,NTEMP,TMI,NPHOT,PHI,IL,NFRST, &
      IPH(30),QM(100),PM(100),PSTOP                                      
      COMMON /PHOTR/ RDAT(80),RTIM(80),RR1,IN10,IPP
!      COMMON / GBLCON / CONC, DIRT                        
      DIMENSION C(91), RT(200),RID(91),CAFTER(91) 
      REAL*8 DDT, CAFTER                                          
      INTEGER SPECIS,REACT,TMI,PHI,AST,RID,numtyp                                   
      CHARACTER*1 TAB							
      CHARACTER*4 ESPECIE(91)                      
                                                
      
      EXTERNAL RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,          &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,     &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DATA IGO/4HCONT/                                                   
      DATA NT1/1/                                                        
      DATA START1/0.0/                                                   
      TAB = '	'
      IPFLAG = 0             
      N=NS-1                                                             
      NS2=NS-2                                                           
      NFRST=1                                                            
!                                                                        
!     TIME INTERVALS                                                     
!                                                                        
      IF (ITITLE(1).NE.IGO) GO TO 30                                     
      IF (NT1.EQ.0) GO TO 30                                             
      TDC=(STOPP-START1)/80.                                             
      TD=START1+TDC                                                      
      J=0                                                                
      NTSV=NT                                                            
      DO 25 I=1,80                                                       
    5 J=J+1                                                              
      IF (J.GT.NTSV) GO TO 20                                            
      IF (SAVTIM(J).GE.TD) GO TO 10                                      
      GO TO 5                                                            
   10 SAVTIM(I)=SAVTIM(J)                                                
      DO 15 K=1,NS                                                       
   15 SAVCON(K,I)=SAVCON(K,J)                                            
      NT=I                                                               
      GO TO 25                                                           
   20 IF (I.EQ.1) NT=1                                                   
      SAVTIM(I)=TD                                                       
   25 TD=TD+TDC                                                          
      GO TO 40                                                           
   30 TDC=(STOPP-START)/80.                                              
      TD=START                                                           
      DO 35 I=1,80                                                       
      TD=TD+TDC                                                          
   35 SAVTIM(I)=TD                                                       
      START1=START                                                       
      NT=1                                                               
   40 NC=0                                                               
!                                                                        
!     INITIALIZE PARAMETERS                                              
!                                                                        
      NHS=0                                                              
      MS=IABS(NQ)                                                        
      TPRNT=TLAST+START                                                  
      TSTEP=PC                                                           
      TFACT=YMAX(2)                                                      
      IF (TFACT.GT.1.) TSTEP=0.                                          
      MT=50-((NS-1)/10)*3-(NR-1)/10                                      
      NN=NS+1-MIN0(NS/10,1)                                              
      IF (YMAX(3).EQ.0.) NPLOT=0
      IF (YMAX(3).NE.0.) NPLOT=1                                         
      IF (YMAX(4).EQ.0.) NOPT=0
      IF (YMAX(4).NE.0.) NOPT=1                                          
      IF (NOPT.EQ.0) MT=50-((NS-1)/10)*2                                 
      NTP=0                                                              
      NT1=NPLOT                                                          
      IF (NPLOT.EQ.0) NT=-1                                              
      IN=INDEX                                                           
      T=START                                                            
      TNEXT=START+1.0E-12      ! Originally: "..=START+1.0"              
      CALL DIFFUN (NS2,TNEXT,C,RT)                                       
      CALL SAVLIN (T,C,NS2)                                              
      NFRST=2                                                            
      IF (T-START) 60,60,50                                              
   45 TNEXT=MIN(TPRNT,STOPP)                                             
      IF (T.EQ.START) GO TO 55                                           
   50 IF (NQ.LT.0.AND.IABS(NQ).NE.0) IN=3                                
   55 CALL DRIVES (NS2,T,H,C,TNEXT,ERR,21,IN,IA,JA)                      
      T=TNEXT                                                            
      IF (T.GE.STOPP) TPRNT=STOPP                                        
   60 IF (T.EQ.START) TIMNW=START                                        
      IF (T.NE.START) TIMNW=TPRNT                                        
      IF (T.EQ.START) GO TO 65                                           
      IF (NQ.LT.0.AND.IABS(NQ).NE.0) TIMNW=T                             
      NHS=NHS+1                                                          
      IF (NHS.GE.MS.OR.T.GE.TPRNT) GO TO 65                              
      GO TO 45                                                           
   65 NHS=0                                                              
      CALL DIFFUN (NS2,TIMNW,C,RT)                                       
      IF (NOPT.LE.0) NTP=NTP+(N/10)+4                                    
      IF (NOPT.GT.0) NTP=NTP+(N/10)*3+(NR-1)/10+12                       
!      IF (NTP.GT.MT.OR.T.EQ.START) WRITE (IOUT,120)                      
      IF (T.EQ.START) GO TO 70                                           
      IF (NOPT.LE.0.AND.NTP.LE.MT) GO TO 75                              
   70 CONTINUE                                                           
      IF (NTP.GT.MT) NTP=0                                               
!      IF (NN.GE.11) WRITE (IOUT,115) (SPECIS(I),I=1,NN)                  
!      IF (NN.LE.10) WRITE (IOUT,125) (SPECIS(I),I=1,NN)                  
!      IF (IPFLAG.EQ.0) WRITE (ITEXT,126) (TAB,SPECIS(I),I=1,NN)              
      IPFLAG = 1                                                      
   75 CONTINUE                                                           
!      IF (NS.GE.11) WRITE (IOUT,130) TIMNW,(C(I),I=1,10),H,(C(I),I=11,NS 
!     &2),TEMP,SG                                                         
!      IF (NS.LE.10) WRITE (IOUT,130) TIMNW,(C(I),I=1,NS2),TEMP,SG,H      
!      WRITE (ITEXT,131) TIMNW,(TAB,C(I),I=1,NS2),TAB,TEMP,TAB,SG,TAB,H  
       
      IF (TIMNW.GE.DDT)  THEN
                DO J=1,ns2
!##ESPECIE(J)=SPECIS(J)
                CAFTER(J)=0.0
!          
                END DO
              DO  KK=1,numtyp
                  DO I=1,NS2
                     
                      IF(RID(KK).EQ.SPECIS(I)) THEN
      
                          CAFTER(KK)=C(I)
!             IF(C(I).NE.0.0)WRITE(*,*)'SALIDA',     &
!            SPECIS(I),RID(KK),CAFTER(KK)
                        

           

                      END IF                                     
                      
                  END DO
              END DO
                          
          ENDIF     
      
          
!      IF (NOPT.EQ.0) WRITE (IOUT,135)                                 
      IF (NOPT.EQ.0) GO TO 95                                            
      RT(NS)=0.                                                          
      RT(N)=0.                                                           
      DO 80 I=1,NS2                                                      
   80 RT(NS)=RT(NS)+RT(I)                                                
!      WRITE (IOUT,105) (RT(I),I=1,NS)                                    
      DO 90 I=1,NR                                                       
      J=KR(I,1)                                                          
      IF (J.EQ.0) RT(I)=0.                                               
      IF (J.EQ.0) GO TO 90                                               
      JT=ITYPE(I)                                                        
      XT=1.                                                              
      DO 85 L=1,JT                                                       
      J=KR(I,L)                                                          
      XJ=1.                                                              
      IF (J.GT.0.AND.J.LE.NS2) XJ=C(J)**ISC(I,L)                         
      IF (ISC(I,L).EQ.-1) XJ=C(J)**SC(I,L)                               
      IF (J.EQ.NS) XJ=SG                                                 
      XT=XT*XJ                                                           
   85 CONTINUE                                                           
      RT(I)=XT*R(I)                                                      
   90 CONTINUE                                                           
!      WRITE (IOUT,110) (RT(I),I=1,NR)                                    
   95 IF (TIMNW.EQ.START) GO TO 55                                       
      IF (T.GE.STOPP) RETURN                                             
      IF (IN.EQ.3.AND.T.LT.TPRNT) GO TO 45                               
      TPRNT=TPRNT*TFACT+TSTEP                                            
      GO TO 45                                                           
!                                                                        
!                                                                        
  105 FORMAT (/,10H NET RATES,2X,1P10E12.3,/,(12X,1P10E12.3))            
  110 FORMAT (//,2X,22HTHE REACTION RATES ARE,/,(1H ,1P10E13.2))         
  115 FORMAT (/,4X,5HTIME ,10(8X,A4),/,2X,8HINTERVAL,7X,A4,9(8X,A4),/,(9    &
      X,10(8X,A4)))                                                      
  120 FORMAT (1H1)                                                       
  125 FORMAT (/,4X,5HTIME ,10(8X,A4),/)                                  
  126 FORMAT ('* TIME',95(A1,A4))                                        
  130 FORMAT (/,1P11E12.3/,11E12.3/,(12X,10E12.3))                    
  131 FORMAT (96(1PE9.2,A1))                                            
  135 FORMAT (1H )                                                    
      END                                                                


      SUBROUTINE RATES (C,N)                                             
!                                                                        
!*********************************************************************** 
!                                                                        
!     THIS SUBROUTINE SETS THE INITIAL CONCENTRATIONS AND CALCULATES THE 
!     RATE CONSTANTS                                                     
!                                                                        
!*********************************************************************** 
!                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)				 
      INTEGER SPECIS,REACT                                               
      COMMON /DATA/ A(200),S(200),TEMP,ERR,START,STOPP,PC,SIG(91),       &
        R(200),BK,SG,DILUT,NP,NR,KR(200,7),IP(91),ITYPE(200),ITITLE(7)
      COMMON /NAMES/ SPECIS(91),REACT(91),NS                             
      COMMON /HEAT/ CV,Q(200),SC(200,7),ISC(200,3),ITEMP                 
      COMMON /ALPHA/ IGO(4),IBLANK,MBLANK,JINTER                         
      DIMENSION C(91)                                                    
      EXTERNAL YFIX, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,          &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,    &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      FCT=1./298.-1./TEMP                                                
      SIGG=0.                                                            
      DO 5 I=1,91                                                        
    5 SIG(I)=0.                                                          
      DO 15 I=1,N                                                        
      DO 10 J=1,NS                                                       
      IF (SPECIS(J).EQ.REACT(I)) SIG(J)=C(I)                             
      IF (SPECIS(J).EQ.REACT(I)) GO TO 15                                
   10 CONTINUE                                                           
      SIGG=SIGG+C(I)                                                     
   15 CONTINUE                                                           
      N=NS                                                               
      M=N-1                                                              
      C(N)=0.                                                            
      DO 20 I=1,M                                                        
      C(N)=C(N)+SIG(I)                                                   
   20 C(I)=SIG(I)                                                        
      BK=0.                                                              
!                                                                        
!     IF SIG(N) DOES NOT EQUAL ZERO, IT IMPLIES THAT THAT THE CONCENTRAT 
!     SPECIES N HAS BEEN READ. OTHERWISE M IS THE SUM OF THE INITIAL     
!     CONCENTRATIONS.                                                    
!                                                                        
      C(N)=C(N)+SIGG                                                     
      IF (SIG(N).NE.0.) BK=SIG(N)-C(N)                                   
      IF (SIG(N).NE.0.) C(N)=SIG(N)                                      
      NP=0                                                               
      DO 25 I=1,NR                                                       
      IF (KR(I,1).EQ.0) GO TO 25                                         
!                                                                        
!     CALCULATE THE RATE CONSTANTS                                       
!                                                                        
      IF (ABS(S(I)).EQ.0.) R(I)=A(I)                                     
      IF (ABS(S(I)).NE.0.) R(I)=A(I)*EXP(S(I)*FCT)                       
   25 CONTINUE                                                           
      IF (NS.LE.9) SPECIS(NS+1)=JINTER                                   
      RETURN                                                             
      END                                            
! 									 
      SUBROUTINE DIFFUN (N,T,X,XT)                                       
!                                                                        
!*********************************************************************** 
!                                                                        
!     THIS SUBROUTINE CALCULATES THE DERIVATIVE VECTOR OF THE ODE|S      
!                                                                        
!*********************************************************************** 
!                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DATA/ A(200),S(200),TEMP,ERR,START,STOPP,PC,SIG(91),          &
        R(200),BK,SG,DILUT,NP,NR,KR(200,7),IP(91),ITYPE(200),ITITLE(7)
      COMMON /NAMES/ SPECIS(91),REACT(91),NS                             
      COMMON /HEAT/ CV,Q(200),SC(200,7),ISC(200,3),ITEMP                 
      COMMON /STORE/ AST(35),IPL(7),TEMEND,NTEMP,TMI,NPHOT,PHI,IL,NFRST,     &
      IPH(30),QM(100),PM(100),PSTOP                                      
      COMMON /PHOTR/ RDAT(80),RTIM(80),RR1,IN10,IPP
      DIMENSION XT(N), X(N)                                              
      INTEGER PHI,TMI, AST                                               
      INTEGER SPECIS, REACT                                              
      EXTERNAL YFIX, RATES, MATRX, XPLOT, VALU, CONVT, SAVLIN,               &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,        &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      P=BK                                                               
      DO 5 I=1,N                                                         
      XT(I)=-DILUT*X(I)                                                  
    5 P=P+X(I)                                                           
      SG=P                                                               
      IF (NFRST.EQ.1) FCT=((3355.7046E-6)*TEMP-1.)/TEMP                  
      IF (NFRST.EQ.1) GO TO 10                                           
      IF (T.EQ.TOLD) GO TO 30                                            
      IF (T.LE.1.) GO TO 30                                              
   10 IF (NPHOT.LE.0) GO TO 30                                           
      IF (NFRST.NE.1) GO TO 15                                           
      PINT=DFLOAT(PHI)                                                   
      PINV=1./PINT                                                       
   15 CONTINUE                                                           
      IF (T.GT.PSTOP) GO TO 20                                           
      IZ=INT(T*PINV)+1                                                   
      Z=T*PINV-DFLOAT(IZ-1)                                              
      IF (T.LE.PINT) R1=PM(1)+(0.5*PM(3)*(Z-1.)+0.5*PM(1)*(Z-3.)-PM(2)*(    &
      Z-2.))*Z                                                           
      IF (T.LE.PINT) GO TO 20                                            
      R1=PM(IZ)+0.25*Z*(5.*PM(IZ+1)-3.0*PM(IZ)-PM(IZ-1)-PM(IZ+2)+(PM(IZ-    &
      1)-PM(IZ)-PM(IZ+1)+PM(IZ+2))*Z)                                    
   20 IF (R1.LT.0.) R1=0.                                                
      DO 25 IK=1,IL                                                      
      IR=IPH(IK)                                                         
      R(IR)=R1*A(IR)                                                     
   25 IF (R(IR).LT.0.) R(IR)=0.                                          
   30 TNOW=TEMP                                                          
      IF (NFRST.EQ.1) GO TO 35                                           
      IF (T.EQ.TOLD) GO TO 50                                            
      IF (T.LE.1.) GO TO 50                                              
   35 CONTINUE                                                           
      IF (T.GT.TEMEND) GO TO 50                                          
      IF (NFRST.NE.1) GO TO 40                                           
      TINT=DFLOAT(TMI)                                                   
      TINV=1./TINT                                                       
   40 CONTINUE                                                           
      IZ=INT(T*TINV)+1                                                   
      Z=T*TINV-DFLOAT(IZ-1)                                              
      IF (T.LE.TINT) TNOW=QM(1)+(0.5*QM(3)*(Z-1.)+0.5*QM(1)*(Z-3.)-QM(2)     &
      *(Z-2.))*Z                                                         
      IF (T.LE.TINT) GO TO 45                                            
      TNOW=QM(IZ)+0.25*Z*((5.*QM(IZ+1)-3.0*QM(IZ)-QM(IZ-1)-QM(IZ+2))+(QM     &
      (IZ-1)-QM(IZ)-QM(IZ+1)+QM(IZ+2))*Z)                                
   45 CONTINUE                                                           
      IF (TNOW.NE.TEMP.AND.TNOW.GT.0.) FCT=((3355.7046E-6)*TNOW-1.)/TNOW 
   50 DO 70 IR=1,NR                                                      
      I=KR(IR,1)                                                         
      IF (I.EQ.0.OR.R(IR).EQ.0.) GO TO 70                                
!                                                                        
!     CHECK FOR A ZEROTH ORDER REACTION                                  
!                                                                        
      RT=1.                                                              
      JT=ITYPE(IR)                                                       
      DO 55 L=1,JT                                                       
      I=KR(IR,L)                                                         
      IF (I.EQ.99.OR.I.EQ.0) GO TO 55                                    
      IF (I.NE.NS) RT=RT*X(I)                                            
      IF (I.EQ.NS) RT=RT*P                                               
   55 CONTINUE                                                           
      IF (ABS(S(IR)).NE.0..AND.T.GT.1.) R(IR)=A(IR)*EXP(S(IR)*FCT)       
      RT=RT*R(IR)                                                        
      DO 60 L=1,JT                                                       
      I=KR(IR,L)                                                         
      IF (I.EQ.0.OR.I.GT.N) GO TO 60                                     
      XT(I)=XT(I)-RT                                                     
   60 CONTINUE                                                           
      DO 65 K=4,7                                                        
      I=KR(IR,K)                                                         
      IF (I.EQ.0) GO TO 70                                               
!                                                                        
!     I WILL BE NEGATIVE IF CLEAN HAS BEEN CALLED                        
!                                                                        
      IF (I.LT.0) GO TO 65                                               
      XT(I)=XT(I)+SC(IR,K)*RT                                            
   65 CONTINUE                                                           
   70 CONTINUE                                                           
      TEMP=TNOW                                                          
      TOLD=T                                                             
      RETURN                                                             
      END                                                                
									 
      SUBROUTINE MATRX (C,KR)                                            
!                                                                        
!*********************************************************************** 
!                                                                        
!     MATRX CREATES A NR X 7 MATRIX OF THE REACTION SCHEME WHERE EACH    
!     ROW REPRESENTS A  REACTION.  THE FIRST THREE COLUMNS ARE REACTANTS 
!     AND THE LAST FOUR ARE THE PRODUCTS.  THE ELEMENTS CORRESPOND       
!     TO THE INDIVIDUAL SPECIES AND WILL BE USED AS SUBSCRIPTS.          
!                                                                        
!*********************************************************************** 
!                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER SPECIS,REACT                                               
      COMMON /DATA/ A(200),S(200),TEMP,ERR,START,STOPP,PC,SIG(91),        &
        R(200),BK,SG,DILUT,NP,NR,IR(200,7),IP(91),ITYPE(200),ITITLE(7)
      COMMON /NAMES/ SPECIS(91),REACT(91),NS                             
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPS1,UROUND,NC,MF1,KFLAG1,JSTART      
      COMMON /ALPHA/ IGO(4),IBLANK,MBLANK,JINTER                         
      COMMON /HEAT/ CV,Q(200),SC(200,7),ISC(200,3),ITEMP                 
      DIMENSION C(91), KR(200,7)                                         
      EXTERNAL YFIX, RATES, DIFFUN, XPLOT, VALU, CONVT, SAVLIN,		    &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,      &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      NOLD=NS+1                                                          
      DO 90 I=1,NR                                                       
!                                                                        
!     SKIP REACTIONS ALREADY PROCESSED                                   
!                                                                        
      IF (IR(I,2).EQ.NOLD.OR.IR(I,3).EQ.NOLD) IR(I,3)=99                 
      IF (IR(I,2).EQ.NOLD) IR(I,2)=0                                     
      IF (IR(I,1).EQ.NOLD) IR(I,1)=99                                    
      IF (IR(I,1).EQ.NOLD) IR(I,3)=99                                    
      IF (IABS(IR(I,1)).NE.100) GO TO 90                                 
!                                                                        
!     IF LESS THAN THREE REACTANTS, FILL FIRST SLOTS.                    
!                                                                        
      IF (KR(I,1).NE.IBLANK) GO TO 5                                     
      IF (KR(I,3).NE.IBLANK) KR(I,1)=KR(I,3)                             
      IF (KR(I,3).NE.IBLANK) SC(I,1)=SC(I,3)                             
      SC(I,3)=1.                                                         
      KR(I,3)=IBLANK                                                     
      IF (KR(I,1).NE.IBLANK) GO TO 5                                     
      KR(I,1)=KR(I,2)                                                    
      KR(I,2)=IBLANK                                                     
      SC(I,1)=SC(I,2)                                                    
      SC(I,2)=1.                                                         
    5 IF (KR(I,2).NE.IBLANK) GO TO 10                                    
      IF (KR(I,3).NE.IBLANK) SC(I,2)=SC(I,3)                             
      SC(I,3)=1.                                                         
      IF (KR(I,3).NE.IBLANK) KR(I,2)=KR(I,3)                             
      KR(I,3)=IBLANK                                                     
   10 DO 15 K=4,7                                                        
!                                                                        
!     GET RID OF M AS A PRODUCT                                          
!                                                                        
      IF (KR(I,K).EQ.MBLANK) KR(I,K)=IBLANK                              
   15 CONTINUE                                                           
!                                                                        
!     IF LESS THAN FOUR PRODUCTS, FILL FIRST SLOTS.                      
!                                                                        
      IF (KR(I,4).NE.IBLANK) GO TO 30                                    
      DO 20 K=1,3                                                        
      INDEX=8-K                                                          
      IF (KR(I,INDEX).NE.IBLANK) GO TO 25                                
   20 CONTINUE                                                           
      INDEX=5                                                            
   25 KR(I,4)=KR(I,INDEX)                                                
      KR(I,INDEX)=IBLANK                                                 
      SC(I,4)=SC(I,INDEX)                                                
      SC(I,INDEX)=1.                                                     
      IF (KR(I,1).NE.IBLANK.OR.KR(I,4).NE.IBLANK) GO TO 30               
      IR(I,1)=0                                                          
      GO TO 90                                                           
   30 CONTINUE                                                           
      IF (KR(I,5).NE.IBLANK) GO TO 35                                    
      IF (KR(I,7).NE.IBLANK) KR(I,5)=KR(I,7)                             
      IF (KR(I,7).NE.IBLANK) SC(I,5)=SC(I,7)                             
      SC(I,7)=1.                                                         
      KR(I,7)=IBLANK                                                     
      IF (KR(I,5).NE.IBLANK) GO TO 35                                    
      KR(I,5)=KR(I,6)                                                    
      KR(I,6)=IBLANK                                                     
      SC(I,5)=SC(I,6)                                                    
      SC(I,6)=1.                                                         
   35 K=KR(I,6)                                                          
      IF (K.NE.IBLANK) GO TO 40                                          
      KR(I,6)=KR(I,7)                                                    
      KR(I,7)=K                                                          
      SC(I,6)=SC(I,7)                                                    
      SC(I,7)=1.                                                         
   40 DO 85 J=1,7                                                        
      K=KR(I,J)                                                          
!                                                                        
!     PROCESS REACTANTS HERE                                             
!                                                                        
      IF (J.GT.3) GO TO 65                                               
!                                                                        
!     ALL M DEPENDENT REACTIONS ARE TO HAVE A 99 IN THE THIRD SLOT       
!                                                                        
      IF (K.NE.MBLANK) GO TO 60                                          
      GO TO (45,50,55),J                                                 
   45 KR(I,1)=KR(I,2)                                                    
      SC(I,1)=SC(I,2)                                                    
   50 KR(I,2)=KR(I,3)                                                    
      SC(I,2)=SC(I,3)                                                    
      SC(I,3)=1.                                                         
      KR(I,3)=MBLANK                                                     
   55 IR(I,3)=99                                                         
   60 K=KR(I,J)                                                          
!                                                                        
!     ZERO ORDER REACTIONS HAVE 99 IN FIRST SLOT                         
!                                                                        
      IF (KR(I,1).EQ.IBLANK) IR(I,1)=99                                  
      IF (KR(I,1).EQ.IBLANK) K=99                                        
!                                                                        
!     ALL BLANKS ARE SET EQUAL TO ZERO                                   
!                                                                        
      IF (J.EQ.3.AND.K.EQ.MBLANK) K=99                                   
   65 IF (K.EQ.MBLANK.OR.K.EQ.IBLANK) IR(I,J)=0                          
      IF (K.EQ.IBLANK.OR.K.EQ.99) GO TO 85                               
      IF (NS.NE.0) GO TO 70                                              
      NS=1                                                               
      GO TO 80                                                           
   70 DO 75 L=1,NS                                                       
      IF (K.NE.SPECIS(L)) GO TO 75                                       
!                                                                        
!     SLOT SET TO SPECIES NUMBER                                         
!                                                                        
      IR(I,J)=L                                                          
      GO TO 85                                                           
   75 CONTINUE                                                           
!                                                                        
!     IF NO SPECIES ARE FOUND, ADD ONE TO THE LIST                       
!                                                                        
      IF (SPECIS(NS).NE.ITEMP) NS=NS+1                                   
   80 SPECIS(NS)=K                                                       
      C(NS+2)=C(NS+1)                                                    
      C(NS+1)=C(NS)                                                      
      C(NS)=0.                                                           
      IR(I,J)=NS                                                         
   85 CONTINUE                                                           
   90 CONTINUE                                                           
      IF (SPECIS(NS).NE.ITEMP) NS=NS+1                                   
      SPECIS(NS)=ITEMP                                                   
      SPECIS(NS+1)=MBLANK                                                
      NS=NS+1                                                            
      DO 95 IK=1,NR                                                      
      DO 95 MT=1,3                                                       
      J=INT(SC(IK,MT)+UROUND)                                            
      ISC(IK,MT)=J                                                       
   95 IF (SC(IK,MT)-DFLOAT(J).GT.4.*UROUND) ISC(IK,MT)=-1                
      DO 100 I=1,NR                                                      
      IF (IR(I,1).EQ.0) GO TO 100                                        
      ITYPE(I)=2                                                         
!                                                                        
!     CHECK FOR 99 CODE AND SUBSTITUTE M WHICH IS THE LAST SPECIES       
!                                                                        
      IF (IR(I,1).EQ.99.AND.IR(I,3).EQ.99) IR(I,1)=NS                    
      IF (IR(I,1).EQ.NS.AND.IR(I,3).EQ.99) IR(I,3)=0                     
      IF (IR(I,2).EQ.0.AND.IR(I,3).EQ.99) IR(I,2)=NS                     
      IF (IR(I,2).EQ.NS.AND.IR(I,3).EQ.99) IR(I,3)=0                     
      IF (IR(I,2).EQ.0) ITYPE(I)=1                                       
      IF (IR(I,3).NE.0) ITYPE(I)=3                                       
      IF (IR(I,3).EQ.99) IR(I,3)=NS                                      
  100 CONTINUE                                                           
      RETURN                                                             
      END                                                                


      SUBROUTINE XPLOT (NTIT,NPNT,NTOT,NAME,SAVTIM,SAVCON,IDT,KCP,JGRID)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! SUBROUTINE ****** P L O T ******                                       
!                                                                        
! THIS SUBROUTINE READS THE PLOT CARDS AND PLOTS THE RESULTS AS PART     
! OF THE PRINTED OUTPUT -- IT DOES NOT DRIVE A PLOTTER.                  
!                                                                        
! SYMBOL DESCRIPTIONS --                                                 
!                                                                        
! CGRID    THE LENGTH OF THE VERTICAL AXIS, PPM                          
! CHIGH    HIGHEST CONCENTRATION VALUE, PPM                              
! CLOW     LOWEST CONCENTRATION VALUE, PPM                               
! CSPAN    CONCENTRATION NORMALIZATION FACTOR                            
! DATA     CONCENTRATION DATA POINTS, PPM, UP TO 80                      
! IDT      IF NOT ZERO THEN PRINT RAW DATA USED FOR PLOTS                
! J        DO-LOOP INDICES OR LOCAL POINTERS                             
! JBLANK   A HOLLERITH WORD OF FOUR BLANK CHARACTERS                     
! JCONC    CONCENTRATION LABELS                                          
! JGRID    THE PLOTTING GRID                                             
! JN       HOLLERITH SYMBOL FOR NO DATA POINTS                           
! JPLUS    THE CHARACTER >+>                                             
! JSTAR    THE CHARACTER >*>                                             
! JSYMB    SYMBOL TO BE USED FOR PLOTTING SAVED POINTS                   
! JVERT    VERTICAL LEGEND                                               
! JX       THE CHARACTER >X>                                             
! K        DO-LOOP INDICES OR LOCAL POINTERS                             
! KCP      IF NOT EQUAL TO ZERO SAVE DATA FOR OFFLINE PLOTTING           
! KCON     CONCENTRATION COORDINATE ON GRID                              
! KTIM     TIME COORDINATE ON GRID                                       
! L        DO-LOOP INDICES OR LOCAL POINTERS                             
! M        DO-LOOP INDICES OR LOCAL POINTERS                             
! MAXCON   LIMIT ON NUMBER OF VERTICAL POINTS                            
! MAXPNT   MAXIMUM NUMBER OF SAVED TIME AND CONCENTRATION POINTS         
! MAXTIM   LIMIT ON NUMBER OF HORIZONTAL POINTS                          
! N        DO-LOOP INDICES OR LOCAL POINTERS                             
! NAME     SPECIES NAMES, ONE PER SPECIES                                
! NDAT     NUMBER OF CONCENTRATION DATA POINTS                           
! NIN      THE FORTRAN INPUT UNIT (NORMALLY 5)                           
! NOUT     THE FORTRAN OUTPUT UNIT NUMBER (NORMALLY 6)                   
! NPLT     THE NUMBER OF SPECIES TO BE PLOTTED                           
! NPNT     NUMBER OF SAVED TIMES AND CONCENTRATIONS                      
! NTEST    SPECIES NAME FOR TESTING                                      
! NTIT     USER-INPUT TITLE FOR PRINTOUT, 3 FOUR-CHARACTER WORDS         
! NTOT     TOTAL NUMBER OF SPECIES                                       
! NX       FLAG FOR CORRECT SPECIES TEST                                 
! SAVCON   SPECIES CONCENTRATIONS, PPM, ONE PER SPECIES AT 80 TIMES      
! SAVTIM   TIMES THAT CONCENTRATIONS ARE SAVED, MIN, UP TO 80 VALUES     
! TGRID    THE LENGTH OF THE HORIZONTAL AXIS, MIN                        
! THIGH    HIGHEST TIME VALUE, MIN                                       
! TIME     TIMES AT WHICH CONCENTRATIONS ARE INPUT, MIN, UP TO 80        
! TLOW     LOWEST TIME VALUE, MIN                                        
! TPRINT   TIMES FOR PRINTOUT ON HORIZONTAL AXIS, MIN                    
! TSPAN    TIME NORMALIZATION FACTOR                                     
!                                                                        
! BEGINNING OF PROGRAM.                                                  
!                                                                        
! ENTRY POINT                                                            
!                                                                        
! SET DIMENSIONS OF INCOMING ARRAYS                                      
!                                                                        
      DIMENSION SAVCON(90,80), SAVTIM(80), NTIT(3), NAME(91)             
!                                                                        
! SET DIMENSIONS OF LOCAL ARRAYS                                         
!                                                                        
      DIMENSION JCONC(5), TIME(80), DATA(80), JN(3)                      
      DIMENSION JGRID(89,40), TPRINT(9), KVERT(52)                       
!                                                                        
! DEFINE THE VERTICAL LABEL AND ESTABLISH ALPHAMERIC DATA                
!                                                                        
      COMMON /APLOT/ JVERT(52,2),JBLANK,JSTAR,JPLUS,JBAR                 
      COMMON /PHOTR/ RDAT(80),RTIM(80),RR1,IN10,IPP
!                                                                        
! DEFINE MISCELLANEOUS DATA VALUES                                       
!                                                                        
      COMMON /INOUT/ NIN,NOUT,ITAPE, ITEXT                               
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, VALU, CONVT, SAVLIN,		  &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,    &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DATA MAXTIM/89/,MAXCON/40/,MAXPNT/80/                              
      DATA TGRID/88./,CGRID/40./,JX/1HX/,JNO/1H /                        
!      DATA KVERT/11*-' P  ',' H  ',' O  ',' T  ',' O  ',' L  ','  
!     &Y ',' S  ',' I  ',' S  ','    ',' C  ',' O  ',' N  ',' S  ',' T  '
!     &,' A  ',' N  ',' T  ',24*-      /                               
!                                                                        
! DEFINE VERTICAL AXIS VIA ASSIGNMENT STATEMENTS                         
!                                                                        
      IF (IDT*IDT.EQ.0) GO TO 10                                         
!     WRITE (NOUT,150)                                                   
!      WRITE (NOUT,145) (NAME(I),I=1,NTOT)                                
!      DO 5 J=1,NPNT                                                      
!    5 WRITE (NOUT,155) SAVTIM(J),(SAVCON(I,J),I=1,NTOT)                  
   10 DO 15 J=1,40                                                       
      JGRID(1,J)=JBAR                                                    
   15 CONTINUE                                                           
      JGRID(1,1)=JPLUS                                                   
      JGRID(1,11)=JPLUS                                                  
      JGRID(1,21)=JPLUS                                                  
      JGRID(1,31)=JPLUS                                                  
      JN(1)=JSTAR                                                        
      JN(2)=JPLUS                                                        
      JN(3)=JX                                                           
      TSPAN=-10.                                                         
!                                                                        
! READ PLOT CONTROL CARD                                                 
!                                                                        
   20 READ (NIN,130) NPLT,JCONC,CLOW,CHIGH,TLOW,THIGH                    
!                                                                        
! TEST FOR END PLOTTING                                                  
!                                                                        
      IF (NPLT.EQ.0) GO TO 105                                           
!      WRITE (NOUT,160) NTIT                                              
!                                                                        
! SET NORMALIZATION FACTORS AND VERTICAL LABELS                          
!                                                                        
      CSPAN=CGRID/(CHIGH-CLOW)                                           
      TSPAN=TGRID/(THIGH-TLOW)                                           
      JVERT(1,2)=JCONC(5)                                                
      JVERT(11,2)=JCONC(4)                                               
      JVERT(21,2)=JCONC(3)                                               
      JVERT(31,2)=JCONC(2)                                               
!                                                                        
! SET HORIZONTAL TIME LABELS                                             
!                                                                        
      DO 25 J=1,9                                                        
      TPRINT(J)=DFLOAT(J-1)/8.*(THIGH-TLOW)+TLOW                         
   25 CONTINUE                                                           
!                                                                        
! CLEAR GRID                                                             
!                                                                        
      DO 35 K=1,MAXCON                                                   
      DO 30 J=2,MAXTIM                                                   
      JGRID(J,K)=JBLANK                                                  
   30 CONTINUE                                                           
   35 CONTINUE                                                           
      NX=0                                                               
      DO 95 LK=1,NPLT                                                    
      READ (NIN,135) NTEST,NDAT,JSYMB                                    
!                                                                        
! TEST NUMBER OF DATA POINTS AND READ DATA                               
!                                                                        
      IF (NDAT.LE.0) GO TO 45                                            
      IF (NDAT.LE.MAXPNT) GO TO 40                                       
!      WRITE (NOUT,185) MAXPNT                                            
      GO TO 125                                                          
!                                                                        
! READ DATA POINTS                                                       
!                                                                        
   40 READ (NIN,140) (TIME(J),DATA(J),J=1,NDAT)                          
!                                                                        
! TEST FOR CORRECT SPECIES NAME                                          
!                                                                        
   45 DO 50 L=1,NTOT                                                     
      IF (NTEST.EQ.NAME(L)) GO TO 55                                     
   50 CONTINUE                                                           
!      WRITE (NOUT,190) NTEST                                             
      NX=NX+1                                                            
      IF (NPLT.EQ.1.OR.NX.EQ.3) GO TO 20                                 
      IF (NPLT.EQ.2.AND.NX.EQ.2) GO TO 20                                
      GO TO 95                                                           
!                                                                        
! IF THERE ARE CALCULATED POINTS, GET THEIR COORDINATES                  
!                                                                        
   55 IF (NPNT.LE.0) GO TO 65                                            
      DO 60 J=1,NPNT                                                     
      KTIM=INT((SAVTIM(J)-TLOW)*TSPAN+1.5)                               
      Y=(SAVCON(L,J)-CLOW)*CSPAN-0.5                                     
      IF (Y.LT.0.0) GO TO 60                                             
      KCON=INT(Y)                                                        
      KCON=MAXCON-KCON                                                   
!                                                                        
! CHECK FOR BEING WITHIN GRID, THEN PLACE ON GRID                        
!                                                                        
      IF (KTIM.LT.2) GO TO 60                                            
      IF (KCON.LT.1) GO TO 60                                            
      IF (KTIM.GT.MAXTIM) GO TO 60                                       
      IF (KCON.GT.MAXCON) GO TO 60                                       
      JGRID(KTIM,KCON)=JSYMB                                             
   60 CONTINUE                                                           
!                                                                        
! IF THERE ARE DATA POINTS, GET THEIR COORDINATES                        
!                                                                        
   65 IF (NDAT.LE.0) GO TO 75                                            
      DO 70 J=1,NDAT                                                     
      KTIM=INT((TIME(J)-TLOW)*TSPAN+1.5)                                 
      Y=(DATA(J)-CLOW)*CSPAN-0.5                                         
      IF (Y.LT.0.0) GO TO 70                                             
      KCON=INT(Y)                                                        
      KCON=MAXCON-KCON                                                   
!                                                                        
! CHECK FOR BEING WITHIN GRID, THEN PLACE ON GRID                        
!                                                                        
      IF (KTIM.LT.2) GO TO 70                                            
      IF (KCON.LT.1) GO TO 70                                            
      IF (KTIM.GT.MAXTIM) GO TO 70                                       
      IF (KCON.GT.MAXCON) GO TO 70                                       
      IF (LK.EQ.1) JGRID(KTIM,KCON)=JSTAR                                
      IF (LK.EQ.2) JGRID(KTIM,KCON)=JPLUS                                
      IF (LK.EQ.3) JGRID(KTIM,KCON)=JX                                   
   70 CONTINUE                                                           
   75 IF (NDAT.GT.0) GO TO 80                                            
!      WRITE (NOUT,195) NTEST,JNO,JSYMB                                   
      GO TO 90                                                           
   80 WRITE (*,195) NTEST,JN(LK),JSYMB                                
!                                                                        
! SAVE DATA FOR OFFLINE PLOTS IF DESIRED                                 
!                                                                        
   90 IF (KCP*KCP.EQ.0) GO TO 95                                         
!      WRITE (ITAPE) NTIT,NAME(L),CLOW,CHIGH,TLOW,THIGH,NDAT,NPNT,(TIME(J 
!     &),J=1,NDAT),(DATA(J),J=1,NDAT),(SAVTIM(J),J=1,NPNT),(SAVCON(L,J),J 
!     &=1,NPNT)                                                           
   95 CONTINUE                                                           
      DO 100 K=1,MAXCON                                                  
!      WRITE (NOUT,165) JVERT(K,1),JVERT(K,2),(JGRID(J,K),J=1,MAXTIM)     
  100 CONTINUE                                                           
!                                                                        
! PRINT THE HORIZONTAL AXIS AND LABELS                                   
!                                                                        
!      WRITE (NOUT,170) JCONC(1)                                          
!     IF (THIGH.GT.80.) WRITE (NOUT,175) TPRINT                          
!      IF (THIGH.LE.80.) WRITE (NOUT,180) TPRINT                          
      GO TO 20                                                           
!                                                                        
! END OF SUBROUTINE -- RETURN TO CALLER                                  
!                                                                        
  105 IF (IPP.LE.0) RETURN                                               
      IF (IN10.LE.0) RETURN                                              
      IF (TSPAN.LE.0.) RETURN                                            
      JVERT(1,2)= 0.80                                                 
      JVERT(11,2)= 0.60                                               
      JVERT(21,2)= 0.40                                                
      JVERT(31,2)= 0.20                                               
      JCONC(1)= 0.00                                                     
      CSPAN=CGRID/0.8                                                 
!      WRITE (NOUT,200) NTIT                                              
      DO 110 I=1,MAXCON                                                  
      DO 110 J=2,MAXTIM                                                  
  110 JGRID(J,I)=JBLANK                                                  
      DO 115 I=1,NPNT                                                    
      KTIM=INT((RTIM(I)-TLOW)*TSPAN+1.5)                                 
      Y=(RDAT(I)-CLOW)*CSPAN-0.5                                         
      IF (Y.LT.0.0) GO TO 115                                            
      KCON=INT(Y)                                                        
      KCON=MAXCON-KCON                                                   
      IF (KCON.LT.1) GO TO 115                                           
      IF (KTIM.LT.2) GO TO 115                                           
      IF (KTIM.GT.MAXTIM) GO TO 115                                      
      IF (KCON.GT.MAXCON) GO TO 115                                      
      JGRID(KTIM,KCON)=JSTAR                                             
  115 CONTINUE                                                           
!      DO 120 K=1,MAXCON                                                  
!      WRITE (NOUT,165) KVERT(K),JVERT(K,2),(JGRID(J,K),J=1,MAXTIM)       
!  120 CONTINUE                                                           
!      WRITE (NOUT,170) JCONC(1)                                          
!      IF (TPRINT(9).GT.80.) WRITE (NOUT,175) TPRINT                      
!      IF (TPRINT(9).LE.80.) WRITE (NOUT,180) TPRINT                      
      RETURN                                                             
  125 CONTINUE                                                           
!	 CLOSE (UNIT=NOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=NIN)
      STOP                                                               
!                                                                        
! LIST OF FORMAT STATEMENTS                                              
!                                                                        
!                                                                        
!                                                                        
  130 FORMAT (5X,I2,8X,5(A4,1X),4F10.0)                                  
  135 FORMAT (A4,1X,I2,1X,A1)                                            
  140 FORMAT (8F10.0)                                                    
  145 FORMAT (/,4X,5HTIME ,4X,10(4X,A4,4X),/,(13X,10(4X,A4,4X)))         
  150 FORMAT (1H1)                                                       
  155 FORMAT (1P11E12.4,/,(12X,10E12.4))                                 
  160 FORMAT (1H1,//////////,62X,3A4,/,32X,7HSPECIES,2X,5HEXPT.,2X,' SIM .')
  165 FORMAT (16X,2A4,89A1)                                              
  170 FORMAT (20X,A4,1H+,8(11H----------+))                              
  175 FORMAT (F25.1,2X,8F11.0,/,62X,14H     TIME     ,/)                 
  180 FORMAT (F25.3,2X,8F11.3,/,62X,14H     TIME     ,/)                 
  185 FORMAT (33H PROGRAM CANNOT HANDLE MORE THAN ,I4,'   PLOT POINTS -- JOB ABORTED.')
  190 FORMAT (33X,13HSPECIES NAME ,A4,26H NOT FOUND IN SPECIES LIST)     
  195 FORMAT (33X,A4,6X,A1,5X,A1)                                        
  200 FORMAT (1H1,//////////,62X,3A4,//)                                 
      END                                                                


      SUBROUTINE VALU (VAL,N,NF,NL)                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /STORE/ AST(35),IPL(7),TEMEND,NTEMP,TMI,NPHOT,PHI,IL,NFRST,    &
      IPH(30),QM(100),PM(100),PSTOP                                      
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPS1,UROUND,NC,MF1,KFLAG1,JSTART      
      DIMENSION AC(5), AD(5), BC(5)                                      
      INTEGER AC,AD,BC,AST,PHI,TMI                                       
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, CONVT, SAVLIN,             &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,        &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DATA IBLK/1H /,IPER/1H./,IZRO/1H0/                                 
      AVAL=LOG10(VAL)                                                    
      IF (AVAL.EQ.0.) RETURN                                             
      DO 5 I=1,5                                                         
      AC(I)=IBLK                                                         
      AD(I)=IZRO                                                         
    5 BC(I)=IZRO                                                         
      IAD=INT(AVAL+UROUND)                                               
      IF (IAD) 30,10,10                                                  
   10 IREM=3-IAD                                                         
      IAD=IAD+1                                                          
      JREM=INT(VAL+UROUND)                                               
      REM=VAL-DFLOAT(JREM)                                               
      IF (REM.GT.UROUND.AND.IREM.GT.0) GO TO 15                          
      CALL CONVT (JREM,AC,5)                                             
      GO TO 50                                                           
   15 CALL CONVT (JREM,AD,IAD)                                           
      JREM=INT(REM*(10.**IREM)+0.1)                                      
      CALL CONVT (JREM,BC,IREM)                                          
      DO 20 J=1,IAD                                                      
   20 AC(J)=AD(J)                                                        
      AC(IAD+1)=IPER                                                     
      DO 25 K=1,IREM                                                     
   25 AC(K+IAD+1)=BC(K)                                                  
      GO TO 40                                                           
   30 IAD=IABS(INT(AVAL-0.1))-1                                          
      IVAL=INT(VAL*10000.+0.1)                                           
      CALL CONVT (IVAL,AC,5)                                             
      AC(1)=IPER                                                         
      DO 35 I=1,IAD                                                      
   35 AC(I+1)=IZRO                                                       
   40 IF (AC(5).NE.IZRO) GO TO 50                                        
      DO 45 K=1,4                                                        
      L=6-K                                                              
   45 AC(L)=AC(L-1)                                                      
      AC(1)=IBLK                                                         
      GO TO 40                                                           
   50 DO 55 I=NF,NL                                                      
      K=I-NF+1                                                           
   55 AST(I)=AC(K)                                                       
      RETURN                                                             
      END                                                                


      SUBROUTINE CONVT (NUM,L,N)                                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!                                                                        
!  SUBROUTINE CONVT CONVERTS INTEGERS TO ALPHANUMERICS                   
!     FOR PRINTING                                                       
!                                                                        
!  ASSUMES VALUE OF INTEGER IS POSITIVE                                  
!                                                                        
      DIMENSION L(5), JDIGIT(10)                                         
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, SAVLIN,        &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,   &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
!                                                                        
      DATA JDIGIT/1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/               
      DATA JBLANK/1H /                                                   
!                                                                        
      NI=NUM                                                             
      DO 5 I=1,N                                                         
      L(I)=JBLANK                                                        
    5 CONTINUE                                                           
!                                                                        
      DO 10 K=1,N                                                        
      I=N-K+1                                                            
      NEXT=NI/10                                                         
      NDX=(NI-NEXT*10)+1                                                 
      L(I)=JDIGIT(NDX)                                                   
      IF (NEXT.LE.0) GO TO 15                                            
      NI=NEXT                                                            
   10 CONTINUE                                                           
!                                                                        
   15 RETURN                                                             
      END                                                                
!									 
      SUBROUTINE SAVLIN (T,C,N)                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)				 
      COMMON /FRPLOT/ SAVCON(90,80),SAVTIM(80),JGRID(89,40),NIT(3),NT    
      COMMON /PHOTR/ RDAT(80),RTIM(80),RR1,IN10,IPP
      DIMENSION C(91)                                                    
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT,           &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,    &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DATA NFRST/1/                                                      
      IF (NFRST.EQ.1) GO TO 5                                            
      IF (T.EQ.TOLD) RETURN                                              
    5 IF (NT.LT.0) RETURN                                                
      NFRST=2                                                            
      TOLD=T                                                             
      IF (NT.GT.80) RETURN                                               
      IF (T.LT.SAVTIM(NT)) RETURN                                        
      SAVTIM(NT)=T                                                       
      RTIM(NT)=T                                                         
      DO 10 I=1,N                                                        
   10 SAVCON(I,NT)=C(I)                                                  
      RDAT(NT)=RR1                                                       
      NT=NT+1                                                            
      RETURN                                                             
      END                                                                


      SUBROUTINE SPLNA (N,X,Y,J,D,C,W)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(10), Y(10), D(2), C(30), W(30)                         
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,    &
        DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,            &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
!     ------------------------------------------------------------------ 
!             OVER THE INTERVAL X(I) TO X(I+1), THE INTERPOLATING        
!             POLYNOMIAL                                                 
!                  Y=Y(I)+A(I)*Z+B(I)*Z**2+E(I)*Z**3                     
!             WHERE Z=(X-X(I))/(X(I+1)-X(I))                             
!             IS USED. THE COEFFICIENTS A(I),B(I) AND E(I) ARE COMPUTED  
!             BY SPLNA AND STORED IN LOCATIONS C(3*I-2),C(3*I-1) AND     
!             C(3*I) RESPECTIVELY.                                       
!             WHILE WORKING IN THE ITH INTERVAL,THE VARIABLE Q WILL      
!             REPRESENT Q=X(I+1) - X(I), AND Y(I) WILL REPRESENT         
!             Y(I+1)-Y(I)                                                
!     ------------------------------------------------------------------ 
!                                                                        
      Q=X(2)-X(1)                                                        
      YI=Y(2)-Y(1)                                                       
      IF (J.EQ.2) GO TO 5                                                
!     ------------------------------------------------------------------ 
!             IF THE FIRST DERIVATIVE AT THE END POINTS IS GIVEN,        
!             A(1) IS KNOWN, AND THE SECOND EQUATION BECOMES             
!             MERELY B(1)+E(1)=YI - Q*D(1).                              
!     ------------------------------------------------------------------ 
      C(1)=Q*D(1)                                                        
      C(2)=1.0                                                           
      W(2)=YI-C(1)                                                       
      GO TO 10                                                           
!     ------------------------------------------------------------------ 
!             IF THE SECOND DERIVATIVE AT THE END POINTS IS GIVEN        
!             B(1) IS KNOWN, THE SECOND EQUATION BECOMES                 
!             A(1)+E(1)=YI-0.5*Q*Q*D(1). DURING THE SOLUTION OF          
!             THE 3N-4 EQUATIONS,A1 WILL BE KEPT IN CELL C(2)            
!             INSTEAD OF C(1) TO RETAIN THE TRIDIAGONAL FORM OF THE      
!             COEFFICIENT MATRIX.                                        
!     ------------------------------------------------------------------ 
    5 C(2)=0.0                                                           
      W(2)=0.5*Q*Q*D(1)                                                  
   10 M=N-2                                                              
      IF (M.LE.0) GO TO 20                                               
!     ------------------------------------------------------------------ 
!             UPPER TRIANGULARIZATION OF THE TRIDIAGONAL SYSTEM OF       
!             EQUATIONS FOR THE COEFFICIENT MATRIX FOLLOWS--             
!     ------------------------------------------------------------------ 
      DO 15 I=1,M                                                        
      AI=Q                                                               
      Q=X(I+2)-X(I+1)                                                    
      H=AI/Q                                                             
      C(3*I)=-H/(2.0-C(3*I-1))                                           
      W(3*I)=(-YI-W(3*I-1))/(2.0-C(3*I-1))                               
      C(3*I+1)=-H*H/(H-C(3*I))                                           
      W(3*I+1)=(YI-W(3*I))/(H-C(3*I))                                    
      YI=Y(I+2)-Y(I+1)                                                   
      C(3*I+2)=1.0/(1.0-C(3*I+1))                                        
   15 W(3*I+2)=(YI-W(3*I+1))/(1.0-C(3*I+1))                              
!     ------------------------------------------------------------------ 
!             E(N-1) IS DETERMINED DIRECTLY FROM THE LAST EQUATION       
!             OBTAINED ABOVE, AND THE FIRST OR SECOND DERIVATIVE         
!             VALUE GIVEN AT THE END POINT.                              
!     ------------------------------------------------------------------ 
   20 IF (J.EQ.1) GO TO 25                                               
      C(3*N-3)=(Q*Q*D(2)/2.0-W(3*N-4))/(3.0-C(3*N-4))                    
      GO TO 30                                                           
   25 C(3*N-3)=(Q*D(2)-YI-W(3*N-4))/(2.0-C(3*N-4))                       
   30 M=3*N-6                                                            
      IF (M.LE.0) GO TO 40                                               
!     ------------------------------------------------------------------ 
!             BACK SOLUTION FOR ALL COEFFICENTS EXCEPT                   
!             A(1) AND B(1) FOLLOWS--                                    
!     ------------------------------------------------------------------ 
      DO 35 II=1,M                                                       
      I=M-II+3                                                           
   35 C(I)=W(I)-C(I)*C(I+1)                                              
   40 IF (J.EQ.1) GO TO 45                                               
!     ------------------------------------------------------------------ 
!             IF THE SECOND DERIVATIVE IS GIVEN AT THE END POINTS,       
!             A(1) CAN NOW BE COMPUTED FROM THE KNOWN VALUES OF          
!             B(1) AND E(1). THEN A(1) AND B(1) ARE PUT INTO THEIR       
!             PROPER PLACES IN THE C ARRAY.                              
!     ------------------------------------------------------------------ 
      C(1)=Y(2)-Y(1)-W(2)-C(3)                                           
      C(2)=W(2)                                                          
      RETURN                                                             
   45 C(2)=W(2)-C(3)                                                     
      RETURN                                                             
      END                                                                


      SUBROUTINE DRIVES (N,T0,H0,Y0,TOUT,EPS,MF,INDEX,IA,JA)             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION IA(1), JA(1), Y0(N)                                      
      DIMENSION Y(91,6)                                                  
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPSC,UROUND,NC,MFC,KFLAG,JSTART       
      COMMON /GEAR2/ YMAX(100)/GEAR3/ERROR(100)/GEAR4/W1(90,3)           
      COMMON /GEAR5/ IW1(91,9)/GEAR6/W2(2000)/GEAR7/IW2(2000)            
      COMMON /GEAR8/ EPSJ,IPTI2,IPTI3,IPTI4,IPTR2,IPTR3,NGRP             
      COMMON /GEAR9/ HUSED,NQUSED,NSTEP,NFE,NJE,NZA,NPL,NPU,NZL,NZU,NZRO 
      COMMON /INOUT/ INP,LOUT,ITAPE, ITEXT                               
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,        &
        SPLNA, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,                 &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DATA NMX/90/,LENW2/2000/,LENIW2/2000/                              
      NGP=0                                                              
      IF (INDEX.EQ.4) GO TO 15                                           
      IF (INDEX.EQ.0) GO TO 30                                           
      IF (INDEX.EQ.2) GO TO 35                                           
      IF (INDEX.EQ.-1) GO TO 40                                          
      IF (INDEX.EQ.3) GO TO 45                                           
      IF (INDEX.NE.1) GO TO 135                                          
      IF (EPS.LE.0.) GO TO 120                                           
      IF (N.LE.0) GO TO 125                                              
      IF ((T0-TOUT)*H0.GE.0.) GO TO 130                                  
      MITER=MF-10*(MF/10)                                                
      IF ((MITER.NE.1).AND.(MITER.NE.2)) GO TO 15                        
      NP1=N+1                                                            
      NZA=IA(NP1)-1                                                      
      NIMAX=LENIW2/2                                                     
      IPTI2=NIMAX+1                                                        
      CALL SORDER (N,IA,JA,IW1,IW1(1,5),NIMAX,IW2,IW2(IPTI2),IER)          
      IPTI2=NZA+1                                                        
      IF (IPTI2+NZA-1.GT.LENIW2) GO TO 145                               
      DO 5 I=1,NP1                                                       
    5 IW1(I,2)=IA(I)                                                     
      DO 10 I=1,NZA                                                      
   10 IW2(I)=JA(I)                                                       
      CALL NSCORD (N,IW1(1,2),IW2,IW1(1,3),IW2(IPTI2),IW1,IW1(1,5),IW1(1    &
      ,8))                                                               
      MAXPL=(LENIW2-NZA)/2                                               
      IPTI3=IPTI2+MAXPL                                                  
      MAXPU=LENIW2-IPTI3+1                                               
      CALL NSSFAC (N,IW1(1,2),IW2,MAXPL,IW1(1,3),IW2(IPTI2),IW1(1,4),MAXPU,    &
	IW1(1,5),IW2(IPTI3),IW1(1,6),Y(1,6),IW1(1,9),Y,Y(1,2),Y(1,3),IW1(1,7),   &
	IW1(1,8),Y(1,4),Y(1,5),IER)                                 
      NPL=IW1(N,4)                                                       
      NPU=IW1(N,6)                                                       
      NZL=IW1(N+1,3)                                                     
      NZU=IW1(N+1,5)                                                     
      IPTR2=NZA+1                                                        
      IPTR3=IPTR2+MAX0(NZA,NZL)                                          
      IF (IPTR3+MAX0(NZA,NZU)-1.GT.LENW2) GO TO 145                      
   15 DO 20 I=1,N                                                        
      YMAX(I)=ABS(Y0(I))                                                 
      IF (YMAX(I).EQ.0.) YMAX(I)=1.E-10                                  
   20 Y(I,1)=Y0(I)                                                       
      NC=N                                                               
      T=T0                                                               
      H=H0                                                               
      NZRO=0                                                             
      TST=EPS*1.E-10                                                     
      DO 25 I=1,N                                                        
   25 IF (Y(I,1).GT.TST) NZRO=NZRO+1                                     
      NZRO=MAX0(NZRO,1)                                                  
      NOLD=NZRO                                                          
      HMIN=ABS(H0)                                                       
      HMAX=ABS(T0-TOUT)*10.                                              
      HMAX=MIN(HMAX,20.D00)                                              
      EPSC=EPS                                                           
      MFC=MF                                                             
      JSTART=0                                                           
      N0=N                                                               
      NMX1=NMX+1                                                         
      EPSJ=SQRT(UROUND)                                                  
      NHCUT=0                                                            
      GO TO 50                                                           
   30 HMAX=ABS(TOUT-TOUTP)*10.                                           
      HMAX=MIN(HMAX,20.D00)                                              
      GO TO 80                                                           
   35 HMAX=ABS(TOUT-TOUTP)*10.                                           
      HMAX=MIN(HMAX,20.D00)                                              
      IF ((T-TOUT)*H.GE.0.) GO TO 150                                    
      GO TO 85                                                           
   40 IF ((T-TOUT)*H.GE.0.) GO TO 140                                    
      JSTART=-1                                                          
      NC=N                                                               
      EPSC=EPS                                                           
      MFC=MF                                                             
   45 CONTINUE                                                           
   50 CALL STIFFS (Y,N0,IA,JA,W1,NMX,IW1,NMX1)                           
      KGO=1-KFLAG                                                        
      GO TO (55,95,110,100),KGO                                          
   55 CONTINUE                                                           
      D=0.                                                               
      NZRO=0                                                             
      DO 70 I=1,NC                                                       
      IF (Y(I,1).GE.0.) GO TO 65                                         
      NGP=NGP+1                                                          
      DO 60 J=1,6                                                        
      K=(J-1)*N+I                                                        
   60 Y(K,1)=0.                                                          
   65 CONTINUE                                                           
      IF (Y(I,1).GT.TST) NZRO=NZRO+1                                     
      AYI=DABS(Y(I,1)) 
	 RATA=DBLE(1.0E-10)
      YMAX(I)=MAX(RATA,AYI)                                            
   70 D=D+(AYI/YMAX(I))**2                                               
      NZRO=MAX0(NZRO,1)                                                  
      IF (NZRO.NE.NOLD) JSTART=-1                                        
      D=D*(UROUND/EPS)**2                                                
      DO 75 II=1,NC                                                      
   75 Y0(II)=Y(II,1)                                                     
      CALL SAVLIN (T,Y0,N)                                               
      IF (D.GT.DFLOAT(N)) GO TO 115                                      
      IF (INDEX.EQ.3) GO TO 150                                          
      IF (INDEX.EQ.2) GO TO 85                                           
   80 IF ((T-TOUT)*H.LT.0.) GO TO 45                                     
      CALL INTERP (TOUT,Y,N0,Y0)                                         
      GO TO 160                                                          
   85 IF (T.GE.TOUT) GO TO 90                                            
      IF (((T+H)-TOUT).LE.0.) GO TO 45                                   
      H=(TOUT-T)*(1.+4.*UROUND)                                          
      JSTART=-1                                                          
      GO TO 45                                                           
   90 JSTART=-1                                                          
      H=MIN(H,1.D00)                                                     
      GO TO 150                                                          
   95 CONTINUE                                                           
  100 IF (NHCUT.EQ.20) GO TO 105    ! Originally: "..EQ.10).."           
      NHCUT=NHCUT+1                                                      
      HMIN=.1*HMIN                                                       
      H=.1*H                                                             
      JSTART=-1                                                          
      GO TO 45                                                           
  105 WRITE (*,165)                                                   
!      IF (KGO.EQ.4) WRITE (LOUT,180) T                                   
!	 CLOSE (UNIT=LOUT)						 
!         CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  110 WRITE (*,170) T,H                                               
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  115 WRITE (*,175) T                                                 
      KFLAG=-2                                                           
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  120 WRITE (*,185)                                                   
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  125 WRITE (*,190)                                                   
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  130 WRITE (*,195)                                                   
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  135 WRITE (*,200) INDEX                                             
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  140 WRITE (*,205) T,TOUT,H                                          
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  145 WRITE (*,210)                                                   
!	 CLOSE (UNIT=LOUT)
!	 CLOSE (UNIT=ITEXT)
!         CLOSE (UNIT=INP)
      STOP                                                               
  150 TOUT=T                                                             
      DO 155 I=1,N                                                       
  155 Y0(I)=Y(I,1)                                                       
      CALL SAVLIN (TOUT,Y0,N)                                            
  160 INDEX=KFLAG                                                        
      TOUTP=TOUT                                                         
      H0=HUSED                                                           
      IF (KFLAG.NE.0) H0=H                                               
      RETURN                                                             
!                                                                        
!                                                                        
  165 FORMAT (//44H PROBLEM APPEARS UNSOLVABLE WITH GIVEN INPUT//)       
  170 FORMAT (//35H KFLAG = -2 FROM INTEGRATOR AT T = ,E16.8,5H  H =,E16.8&
      /52H  THE REQUESTED ERROR IS SMALLER THAN CAN BE HANDLED//)      
  175 FORMAT (//37H INTEGRATION HALTED BY DRIVER AT T = ,E16.8,&
      '56H  EPS TOO SMALL TO BE ATTAINED FOR THE MACHINE PRECISION')
  180 FORMAT (//35H KFLAG = -3 FROM INTEGRATOR AT T = ,E16.8,&
      ' CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED')
  185 FORMAT (//28H ILLEGAL INPUT.. EPS .LE. 0.//)                       
  190 FORMAT (//25H ILLEGAL INPUT.. N .LE. 0//)                          
  195 FORMAT (//36H ILLEGAL INPUT.. (T0-TOUT)*H .GE. 0.//)               
  200 FORMAT (//24H ILLEGAL INPUT.. INDEX =,I5//)                        
  205 FORMAT (//44H INDEX = -1 ON INPUT WITH (T-TOUT)*H .GE. 0./4H T =,E    &
      16.8,9H   TOUT =,E16.8,6H   H =,E16.8)                             
  210 FORMAT ('     INSUFFICIENT WORKING STORAGE IN IW2 OR W2 ')         
      END                                                                


      SUBROUTINE STIFFS (Y,N0,IA,JA,W1,NMX,IW1,NMX1)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPS,UROUND,N,MF,KFLAG,JSTART          
      COMMON /GEAR2/ YMAX(100)/GEAR3/ERROR(100)                          
      COMMON /GEAR6/ W2(2000)/GEAR7/IW2(2000)                            
      COMMON /GEAR8/ EPSJ,IPTI2,IPTI3,IPTI4,IPTR2,IPTR3,NGRP             
      COMMON /GEAR9/ HUSED,NQUSED,NSTEP,NFE,NJE,IDUMMY(5),NZRO           
      COMMON /DATA/ A(200),S(200),TEMP,ERR,START,STOPP,PC,SIG(91),          &
        R(200),BK,SG,DILUT,NP,NR,KR(200,7),IP(91),ITYPE(200),ITITLE(7)
      COMMON /HEAT/ CV,Q(200),SC(200,7),ISC(200,3),ITEMP                 
      DIMENSION Y(N0,*), IA(1), JA(1), W1(NMX,*), IW1(NMX1,*)            
      DIMENSION EL(13), TQ(4), RT(3)                                     
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,       &
        SPLNA, DRIVES, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,                &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DATA EL(2)/1./,OLDL0/1./                                           
      KFLAG=0                                                            
      TOLD=T                                                             
      IF (JSTART.GT.0) GO TO 50                                          
      IF (JSTART.NE.0) GO TO 10                                          
      CALL DIFFUN (N,T,Y,W1)                                             
      DO 5 I=1,N                                                         
    5 Y(I,2)=H*W1(I,1)                                                   
      METH=MF/10                                                         
      MITER=MF-10*METH                                                   
      NQ=1                                                               
      L=2                                                                
      IDOUB=3                                                            
      RMAX=1.E4                                                          
      RC=0.                                                              
      CRATE=1.                                                           
      HOLD=H                                                             
      MFOLD=MF                                                           
      NSTEP=0                                                            
      NSTEPJ=0                                                           
      NFE=1                                                              
      NJE=0                                                              
      IRET=3                                                             
      GO TO 15                                                           
   10 IF (MF.EQ.MFOLD) GO TO 25                                          
      MEO=METH                                                           
      MIO=MITER                                                          
      METH=MF/10                                                         
      MITER=MF-10*METH                                                   
      MFOLD=MF                                                           
      IF (MITER.NE.MIO) IWEVAL=MITER                                     
      IF (METH.EQ.MEO) GO TO 25                                          
      IDOUB=L+1                                                          
      IRET=1                                                             
   15 CALL COSET (METH,NQ,EL,TQ,MAXDER)                                  
      LMAX=MAXDER+1                                                      
      RC=RC*EL(1)/OLDL0                                                  
      OLDL0=EL(1)                                                        
   20 FN=DFLOAT(NZRO)                                                    
      EDN=FN*(TQ(1)*EPS)**2                                              
      E=FN*(TQ(2)*EPS)**2                                                
      EUP=FN*(TQ(3)*EPS)**2                                              
      BND=FN*(TQ(4)*EPS)**2                                              
      EPSOLD=EPS                                                         
      NOLD=NZRO                                                          
      GO TO (30,35,50),IRET                                              
   25 IF (EPS.EQ.EPSOLD.AND.NZRO.EQ.NOLD) GO TO 30                       
      IRET=1                                                             
      GO TO 20                                                           
   30 IF (H.EQ.HOLD) GO TO 50                                            
      RH=H/HOLD                                                          
      H=HOLD                                                             
      IREDO=3                                                            
      GO TO 40                                                           
   35 RH=MAX(RH,HMIN/ABS(H))                                             
   40 RH=MIN(RH,HMAX/ABS(H),RMAX)                                        
      R1=1.                                                              
      DO 45 J=2,L                                                        
      R1=R1*RH                                                           
      DO 45 I=1,N                                                        
   45 Y(I,J)=Y(I,J)*R1                                                   
      H=H*RH                                                             
      RC=RC*RH                                                           
      IDOUB=L+1                                                          
      IF (IREDO.EQ.0) GO TO 290                                          
   50 IF (ABS(RC-1.).GT.0.3) IWEVAL=MITER                                
      IF (NSTEP.GE.NSTEPJ+20) IWEVAL=MITER                               
      T=T+H                                                              
      DO 55 J1=1,NQ                                                      
      DO 55 J2=J1,NQ                                                     
      J=(NQ+J1)-J2                                                       
      DO 55 I=1,N                                                        
   55 Y(I,J)=Y(I,J)+Y(I,J+1)                                             
   60 DO 65 I=1,N                                                        
   65 ERROR(I)=0.                                                        
      M=0                                                                
      CALL DIFFUN (N,T,Y,W1(1,2))                                        
      NFE=NFE+1                                                          
      IF (IWEVAL.LE.0) GO TO 140                                         
      IWEVAL=0                                                           
      RC=1.                                                              
      NJE=NJE+1                                                          
      NSTEPJ=NSTEP                                                       
      CON=-H*EL(1)                                                       
      ISV=M                                                              
      LSV=L                                                              
      NZ=IA(N+1)-1                                                       
      DO 70 I=1,NZ                                                       
   70 W2(I)=0.                                                           
      DO 125 IR=1,NR                                                     
      IF (KR(IR,1).EQ.0.OR.KR(IR,1).EQ.99) GO TO 125                     
      MT=ITYPE(IR)                                                       
      DO 85 I=1,MT                                                       
      IF (KR(IR,I).EQ.N+2) GO TO 85                                      
      JX=I+1-I/3*3                                                       
      KX=I+2-I/2*3                                                       
      J=KR(IR,JX)                                                        
      K=KR(IR,KX)                                                        
      XJ=1.                                                              
      IF (J.EQ.0) GO TO 75                                               
      XJ=Y(J,1)                                                          
      IF (J.EQ.N+2) XJ=SG                                                
   75 XK=1.                                                              
      IF (K.EQ.0) GO TO 80                                               
      XK=Y(K,1)                                                          
      IF (K.EQ.N+2) XK=SG                                                
   80 RT(I)=R(IR)*XJ*XK                                                  
   85 CONTINUE                                                           
      DO 120 K=1,MT                                                      
      I=KR(IR,K)                                                         
      IF (I.EQ.N+2) GO TO 120                                            
      DO 100 L=1,MT                                                      
      J=KR(IR,L)                                                         
      IF (J.EQ.N+2) GO TO 100                                            
      M=IA(J)-1                                                          
   90 M=M+1                                                              
      IF (I-JA(M)) 90,95,90                                              
   95 W2(M)=W2(M)-RT(L)                                                  
  100 CONTINUE                                                           
      DO 115 L=4,7                                                       
      J=KR(IR,L)                                                         
      M=IA(I)-1                                                          
      IF (J) 115,120,105                                                 
  105 M=M+1                                                              
      IF (J-JA(M)) 105,110,105                                           
  110 W2(M)=W2(M)+RT(K)*SC(IR,L)                                         
  115 CONTINUE                                                           
  120 CONTINUE                                                           
  125 CONTINUE                                                           
      DO 135 J=1,N                                                       
      KMIN=IA(J)                                                         
      KMAX=IA(J+1)-1                                                     
      DO 130 K=KMIN,KMAX                                                 
      W2(K)=W2(K)*CON                                                    
      IF (JA(K).EQ.J) W2(K)=W2(K)+1.-CON*DILUT                           
  130 CONTINUE                                                           
  135 CONTINUE                                                           
      CALL NSCORA (N,IA,JA,W2,IW1(1,2),W2(IPTR3),W2(IPTR2),IW1,IW1(1,7),     &
      IW1(1,8))                                                          
      CALL NSNFAC (N,IW1(1,2),IW2,W2,IW1(1,3),IW2(IPTI2),IW1(1,4),W2(IPTR2),   &
	W1(1,3),IW1(1,5),IW2(IPTI3),IW1(1,6),W2(IPTR3),W1,IW1(1,7),IW1(1,8),IER)                                                         
      M=ISV                                                              
      L=LSV                                                              
      IF (IER.NE.0) GO TO 160                                            
  140 DO 145 I=1,N                                                       
      IF (M.LE.0) GO TO 145                                              
      IF (-H*W1(I,2)*10..GT.Y(I,1)) GO TO 155                            
  145 W1(I,1)=H*W1(I,2)-(Y(I,2)+ERROR(I))                                
      CALL NSBSLV (N,IW1,IW1,IW1(1,3),IW2(IPTI2),IW1(1,4),W2(IPTR2),W1(1,3),    &
	IW1(1,5),IW2(IPTI3),IW1(1,6),W2(IPTR3),W1(1,2),W1,W2)          
      D=0.                                                               
      DO 150 I=1,N                                                       
      ERROR(I)=ERROR(I)+W1(I,2)                                          
      D=D+(W1(I,2)/YMAX(I))**2                                           
  150 W1(I,1)=Y(I,1)+EL(1)*ERROR(I)                                      
      IF (M.NE.0) CRATE=MAX(.9*CRATE,D/D1)                               
      IF ((D*MIN(1.D00,2.*CRATE)).LE.BND) GO TO 175                      
      D1=D                                                               
      M=M+1                                                              
      IF (M.EQ.3) GO TO 155                                              
      CALL DIFFUN (N,T,W1,W1(1,2))                                       
      GO TO 140                                                          
  155 NFE=NFE+2                                                          
      IF (IWEVAL.EQ.-1) GO TO 170                                        
  160 T=TOLD                                                             
      RMAX=2.                                                            
      DO 165 J1=1,NQ                                                     
      DO 165 J2=J1,NQ                                                    
      J=(NQ+J1)-J2                                                       
      DO 165 I=1,N                                                       
  165 Y(I,J)=Y(I,J)-Y(I,J+1)                                             
      IF (ABS(H).LE.HMIN*1.00001) GO TO 285                              
      RH=.25                                                             
      IREDO=1                                                            
      GO TO 35                                                           
  170 IWEVAL=MITER                                                       
      GO TO 60                                                           
  175 IF (MITER.NE.0) IWEVAL=-1                                          
      NFE=NFE+M                                                          
      D=0.                                                               
      DO 180 I=1,N                                                       
  180 D=D+(ERROR(I)/YMAX(I))**2                                          
      IF (D.GT.E) GO TO 195                                              
      KFLAG=0                                                            
      IREDO=0                                                            
      NSTEP=NSTEP+1                                                      
      HUSED=H                                                            
      NQUSED=NQ                                                          
      DO 185 J=1,L                                                       
      DO 185 I=1,N                                                       
  185 Y(I,J)=Y(I,J)+EL(J)*ERROR(I)                                       
      IF (IDOUB.EQ.1) GO TO 205                                          
      IDOUB=IDOUB-1                                                      
      IF (IDOUB.GT.1) GO TO 295                                          
      IF (L.EQ.LMAX) GO TO 295                                           
      DO 190 I=1,N                                                       
  190 Y(I,LMAX)=ERROR(I)                                                 
      GO TO 295                                                          
  195 KFLAG=KFLAG-1                                                      
      T=TOLD                                                             
      DO 200 J1=1,NQ                                                     
      DO 200 J2=J1,NQ                                                    
      J=(NQ+J1)-J2                                                       
      DO 200 I=1,N                                                       
  200 Y(I,J)=Y(I,J)-Y(I,J+1)                                             
      RMAX=2.                                                            
      IF (ABS(H).LE.HMIN*1.00001) GO TO 275                              
      IF (KFLAG.LE.-3) GO TO 265                                         
      IREDO=2                                                            
      PR3=1.E+20                                                         
      GO TO 215                                                          
  205 PR3=1.E+20                                                         
      IF (L.EQ.LMAX) GO TO 215                                           
      D1=0.                                                              
      DO 210 I=1,N                                                       
  210 D1=D1+((ERROR(I)-Y(I,LMAX))/YMAX(I))**2                            
      ENQ3=.5/DFLOAT(L+1)                                                
      PR3=((D1/EUP)**ENQ3)*1.4+1.4E-6                                    
  215 ENQ2=.5/DFLOAT(L)                                                  
      PR2=((D/E)**ENQ2)*1.2+1.2E-6                                       
      PR1=1.E+20                                                         
      IF (NQ.EQ.1) GO TO 225                                             
      D=0.                                                               
      DO 220 I=1,N                                                       
  220 D=D+(Y(I,L)/YMAX(I))**2                                            
      ENQ1=.5/DFLOAT(NQ)                                                 
      PR1=((D/EDN)**ENQ1)*1.3+1.3E-6                                     
  225 IF (PR2.LE.PR3) GO TO 230                                          
      IF (PR3.LT.PR1) GO TO 240                                          
      GO TO 235                                                          
  230 IF (PR2.GT.PR1) GO TO 235                                          
      NEWQ=NQ                                                            
      RH=1./PR2                                                          
      GO TO 255                                                          
  235 NEWQ=NQ-1                                                          
      RH=1./PR1                                                          
      GO TO 255                                                          
  240 NEWQ=L                                                             
      RH=1./PR3                                                          
      IF (RH.LT.1.1) GO TO 250                                           
      DO 245 I=1,N                                                       
  245 Y(I,NEWQ+1)=ERROR(I)*EL(L)/DFLOAT(L)                               
      GO TO 260                                                          
  250 IDOUB=10                                                           
      GO TO 295                                                          
  255 IF ((KFLAG.EQ.0).AND.(RH.LT.1.1)) GO TO 250                        
      IF (NEWQ.EQ.NQ) GO TO 35                                           
  260 NQ=NEWQ                                                            
      L=NQ+1                                                             
      IRET=2                                                             
      GO TO 15                                                           
  265 IF (KFLAG.EQ.-9) GO TO 280                                         
      RH=10.**KFLAG                                                      
      RH=MAX(HMIN/ABS(H),RH)                                             
      H=H*RH                                                             
      CALL DIFFUN (N,T,Y,W1)                                             
      NFE=NFE+1                                                          
      DO 270 I=1,N                                                       
  270 Y(I,2)=H*W1(I,1)                                                   
      IWEVAL=MITER                                                       
      IDOUB=10                                                           
      IF (NQ.EQ.1) GO TO 50                                              
      NQ=1                                                               
      L=2                                                                
      IRET=3                                                             
      GO TO 15                                                           
  275 KFLAG=-1                                                           
      GO TO 295                                                          
  280 KFLAG=-2                                                           
      GO TO 295                                                          
  285 KFLAG=-3                                                           
      GO TO 295                                                          
  290 RMAX=100.                                                          
  295 HOLD=H                                                             
      JSTART=NQ                                                          
      RETURN                                                             
      END                                                                


      SUBROUTINE COSET (METH,NQ,EL,TQ,MAXDER)                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PERTST(12,2,3), EL(13), TQ(4)                            
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,       &
        SPLNA, DRIVES, STIFFS, NSCORA, NSNFAC, NSBSLV, YSMER,               &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD 
      DATA PERTST/1.,1.,2.,1.,.3158,.07407,.01391,.002182,.0002945,.00003492,  &
      .000003692,.0000003524,1.,1.,.5,.1667,.04167,1.,1.,1.,1.,1.,1.,1.,       &
      2.,12.,24.,37.89,53.33,70.08,87.97,106.9,126.7,147.4,168.8,191.0,        &  
      2.0,4.5,7.333,10.42,13.7,1.,1.,1.,1.,1.,1.,1.,12.0,24.0,37.89,           &
      53.33,70.08,87.97,106.9,126.7,147.4,168.8,191.0,1.,3.0,6.0,9.167,12.5,   &
      1.,1.,1.,1.,1.,1.,1.,1./                                       
    5 MAXDER=5                                                           
      GO TO (10,15,20,25,30),NQ                                          
   10 EL(1)=1.0                                                          
      GO TO 35                                                           
   15 EL(1)=6.6666666666667E-01                                          
      EL(3)=3.3333333333333E-01                                          
      GO TO 35                                                           
   20 EL(1)=5.4545454545455E-01                                          
      EL(3)=EL(1)                                                        
      EL(4)=9.0909090909091E-02                                          
      GO TO 35                                                           
   25 EL(1)=0.48                                                         
      EL(3)=0.7                                                          
      EL(4)=0.2                                                          
      EL(5)=0.02                                                         
      GO TO 35                                                           
   30 EL(1)=4.3795620437956E-01                                          
      EL(3)=8.2116788321168E-01                                          
      EL(4)=3.1021897810219E-01                                          
      EL(5)=5.4744525547445E-02                                          
      EL(6)=3.6496350364964E-03                                          
   35 DO 40 K=1,3                                                        
   40 TQ(K)=PERTST(NQ,METH,K)                                            
      TQ(4)=.5*TQ(2)/DFLOAT(NQ+2)                                        
      RETURN                                                             
      END                                                                
! 									 
      SUBROUTINE NSCORA (N,IA,JA,A,IAP,JAWORK,AWORK,C,IR,ICT)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IA(1),JA(1),IAP(1),C(1),JAWORK(1),IR(1),ICT(1)             
      DOUBLE PRECISION A(1),AWORK(1)                                               
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,    &
        SPLNA, DRIVES, STIFFS, COSET, NSNFAC, NSBSLV, YSMER,             &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DO 5 K=1,N                                                         
      ICK=C(K)                                                           
    5 IR(ICK)=K                                                          
      JMIN=1                                                             
      DO 15 K=1,N                                                        
      ICK=C(K)                                                           
      JMAX=JMIN+IA(ICK+1)-IA(ICK)-1                                      
      IF (JMIN.GT.JMAX) GO TO 15                                         
      IAINK=IA(ICK)-1                                                    
      DO 10 J=JMIN,JMAX                                                  
      IAINK=IAINK+1                                                      
      JAOUTJ=JA(IAINK)                                                   
      JAOUTJ=IR(JAOUTJ)                                                  
      JAWORK(J)=JAOUTJ                                                   
   10 AWORK(J)=A(IAINK)                                                  
   15 JMIN=JMAX+1                                                        
      DO 20 I=1,N                                                        
   20 ICT(I)=IAP(I)                                                      
      JMIN=1                                                             
      DO 30 I=1,N                                                        
      ICK=C(I)                                                           
      JMAX=JMIN+IA(ICK+1)-IA(ICK)-1                                      
      IF (JMIN.GT.JMAX) GO TO 30                                         
      DO 25 J=JMIN,JMAX                                                  
      JAOUTJ=JAWORK(J)                                                   
      ICTJ=ICT(JAOUTJ)                                                   
      A(ICTJ)=AWORK(J)                                                   
      ICT(JAOUTJ)=ICTJ+1                                                 
   25 CONTINUE                                                           
   30 JMIN=JMAX+1                                                        
      RETURN                                                             
      END                                                                
!									 
      SUBROUTINE NSNFAC (N,IA,JA,A,IL,JL,ISL,L,D,IU,JU,ISU,U,X,IRL,JRL,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IA(1),JA(1),IL(1),JL(1),ISL(1)                             
      INTEGER IU(1),JU(1),ISU(1),IRL(1),JRL(1)                           
      DOUBLE PRECISION A(1),L(1),D(1),U(1),X(1)                                    
      DOUBLE PRECISION LKI                                                         
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,      &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSBSLV, YSMER,               &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      IER=0                                                              
      DO 5 K=1,N                                                         
      IRL(K)=IL(K)                                                       
    5 JRL(K)=0                                                           
      DO 90 K=1,N                                                        
      X(K)=0.                                                            
      I1=0                                                               
      IF (JRL(K).EQ.0) GO TO 15                                          
      I=JRL(K)                                                           
   10 I2=JRL(I)                                                          
      JRL(I)=I1                                                          
      I1=I                                                               
      X(I)=0.                                                            
      I=I2                                                               
      IF (I.NE.0) GO TO 10                                               
   15 JMIN=ISU(K)                                                        
      JMAX=JMIN+IU(K+1)-IU(K)-1                                          
      IF (JMIN.GT.JMAX) GO TO 25                                         
      DO 20 J=JMIN,JMAX                                                  
      JUJ=JU(J)                                                          
   20 X(JUJ)=0.                                                          
   25 JMIN=IA(K)                                                         
      JMAX=IA(K+1)-1                                                     
      DO 30 J=JMIN,JMAX                                                  
      JAJ=JA(J)                                                          
   30 X(JAJ)=A(J)                                                        
      I=I1                                                               
      IF (I.EQ.0) GO TO 50                                               
   35 IRLI=IRL(I)                                                        
      LKI=-X(I)                                                          
      L(IRLI)=-LKI                                                       
      JMIN=IU(I)                                                         
      JMAX=IU(I+1)-1                                                     
      IF (JMIN.GT.JMAX) GO TO 45                                         
      ISUB=ISU(I)-1                                                      
      DO 40 J=JMIN,JMAX                                                  
      ISUB=ISUB+1                                                        
      JUJ=JU(ISUB)                                                       
   40 X(JUJ)=X(JUJ)+LKI*U(J)                                             
   45 I=JRL(I)                                                           
      IF (I.NE.0) GO TO 35                                               
   50 IF (X(K).EQ.0.) GO TO 95                                           
      DK=1./X(K)                                                         
      D(K)=DK                                                            
      IF (K.EQ.N) GO TO 90                                               
      JMIN=IU(K)                                                         
      JMAX=IU(K+1)-1                                                     
      IF (JMIN.GT.JMAX) GO TO 60                                         
      ISUB=ISU(K)-1                                                      
      DO 55 J=JMIN,JMAX                                                  
      ISUB=ISUB+1                                                        
      JUJ=JU(ISUB)                                                       
   55 U(J)=X(JUJ)*DK                                                     
   60 CONTINUE                                                           
      I=I1                                                               
      IF (I.EQ.0) GO TO 85                                               
   65 IRL(I)=IRL(I)+1                                                    
      I1=JRL(I)                                                          
      IF (IRL(I).GE.IL(I+1)) GO TO 80                                    
      ISLB=IRL(I)-IL(I)+ISL(I)                                           
      J=JL(ISLB)                                                         
   70 IF (I.GT.JRL(J)) GO TO 75                                          
      J=JRL(J)                                                           
      GO TO 70                                                           
   75 JRL(I)=JRL(J)                                                      
      JRL(J)=I                                                           
   80 I=I1                                                               
      IF (I.NE.0) GO TO 65                                               
   85 ISLK=ISL(K)                                                        
      IF (IRL(K).GE.IL(K+1)) GO TO 90                                    
      J=JL(ISLK)                                                         
      JRL(K)=JRL(J)                                                      
      JRL(J)=K                                                           
   90 CONTINUE                                                           
      RETURN                                                             
   95 IER=K                                                              
      RETURN                                                             
      END                                                                
!									 
      SUBROUTINE NSBSLV (N,R,C,IL,JL,ISL,L,D,IU,JU,ISU,U,X,B,Y)          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(1), IL(1), JL(1), IU(1), JU(1), C(1), ISL(1), ISU(1)   
      DIMENSION L(1), X(1), B(1), U(1), Y(1), D(1)                       
      INTEGER R,RK,C,CK                                                  
      DOUBLE PRECISION L                                                           
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,        &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, YSMER,                 &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DO 5 K=1,N                                                         
      RK=R(K)                                                            
    5 Y(K)=B(RK)                                                         
      DO 15 K=1,N                                                        
      JMIN=IL(K)                                                         
      JMAX=IL(K+1)-1                                                     
      YK=-D(K)*Y(K)                                                      
      Y(K)=-YK                                                           
      IF (JMIN.GT.JMAX) GO TO 15                                         
      ISLB=ISL(K)-1                                                      
      DO 10 J=JMIN,JMAX                                                  
      ISLB=ISLB+1                                                        
      JLJ=JL(ISLB)                                                       
   10 Y(JLJ)=Y(JLJ)+YK*L(J)                                              
   15 CONTINUE                                                           
      K=N                                                                
      DO 30 I=1,N                                                        
      SUM=-Y(K)                                                          
      JMIN=IU(K)                                                         
      JMAX=IU(K+1)-1                                                     
      IF (JMIN.GT.JMAX) GO TO 25                                         
      ISUB=ISU(K)-1                                                      
      DO 20 J=JMIN,JMAX                                                  
      ISUB=ISUB+1                                                        
      JUJ=JU(ISUB)                                                       
   20 SUM=SUM+U(J)*Y(JUJ)                                                
   25 Y(K)=-SUM                                                          
      CK=C(K)                                                            
      X(CK)=-SUM                                                         
   30 K=K-1                                                              
      RETURN                                                             
      END                                                                


      SUBROUTINE YSMER (A,K,A1)                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,A1
	                                               
      COMMON /INOUT/ INP,LOUT,ITAPE, ITEXT                               
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,        &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV,                &
        INTERP, SPARS, SORDER, NSSFAC, NSCORD
      WRITE (*,5) A,K,A1                                     
      RETURN                                                             
!                                                                        
!                                                                        
!                                                                        
    5 FORMAT (1X,A10,I6,2A10)                                            
      END                                                                


      SUBROUTINE INTERP (TOUT,Y,N0,Y0)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /GEAR1/ T,H,DUMMY(4),N,IDUMMY(2),JSTART                     
      DIMENSION Y0(N0), Y(N0,1)                                          
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,   &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,    &
        SPARS, SORDER, NSSFAC, NSCORD
      DO 5 I=1,N                                                         
    5 Y0(I)=Y(I,1)                                                       
      L=JSTART+1                                                         
      S=(TOUT-T)/H                                                       
      S1=1.                                                              
      DO 15 J=2,L                                                        
      S1=S1*S                                                            
      DO 10 I=1,N                                                        
   10 Y0(I)=Y0(I)+S1*Y(I,J)                                              
   15 CONTINUE                                                           
      RETURN                                                             
      END                                                                
!									 
      SUBROUTINE SPARS (IA,JA,N)                                         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DATA/ A(200),S(200),TEMP,ERR,START,STOPP,PC,SIG(91),              &
        R(200),BK,SG,DILUT,NP,NR,KR(200,7),IP(91),ITYPE(200),ITITLE(7)
      DIMENSION IA(99), JA(1000)       ! size previously was set to N    
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,        &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,         &
        INTERP, SORDER, NSSFAC, NSCORD
      DO 5 I=1,N                                                         
    5 IA(I)=1                                                            
      KT=0                                                               
      IA(N+1)=1                                                          
      JA(1)=0                                                            
      DO 70 IR=1,NR                                                      
      IF (KR(IR,1).EQ.N+2) GO TO 70                                      
      IF (KR(IR,1).EQ.0.OR.KR(IR,1).EQ.99) GO TO 70                      
      MT=ITYPE(IR)                                                       
      DO 65 K=1,MT                                                       
      I=KR(IR,K)                                                         
      IF (KR(IR,K).EQ.N+2) GO TO 65                                      
      DO 30 L=1,MT                                                       
      J=KR(IR,L)                                                         
      IF (KR(IR,L).EQ.N+2) GO TO 30                                      
      K1=IA(J)                                                           
      K2=IA(J+1)-1                                                       
      IF (K1.GT.K2) GO TO 15                                             
      DO 10 M=K1,K2                                                      
   10 IF (I.EQ.JA(M)) GO TO 30                                           
   15 DO 20 M=J,N                                                        
   20 IA(M+1)=IA(M+1)+1                                                  
      KT=KT+1                                                            
      KD=KT-K2                                                           
      K2=K2+1                                                            
      DO 25 M=1,KD                                                       
   25 JA(KT+2-M)=JA(KT+1-M)                                              
      JA(K2)=I                                                           
   30 CONTINUE                                                           
      K1=IA(I)                                                           
      DO 60 L=4,7                                                        
      K2=IA(I+1)-1                                                       
      J=KR(IR,L)                                                         
      IF (KR(IR,L).EQ.N+2) GO TO 60                                      
      IF (J) 60,65,35                                                    
   35 IF (K1.GT.K2) GO TO 45                                             
      DO 40 M=K1,K2                                                      
      IF (J.EQ.JA(M)) GO TO 60                                           
   40 CONTINUE                                                           
   45 DO 50 M=I,N                                                        
   50 IA(M+1)=IA(M+1)+1                                                  
      KT=KT+1                                                            
      KD=KT-K2                                                           
      K2=K2+1                                                            
      DO 55 M=1,KD                                                       
   55 JA(KT+2-M)=JA(KT+1-M)                                              
      JA(K2)=J                                                           
   60 CONTINUE                                                           
   65 CONTINUE                                                           
   70 CONTINUE                                                           
      DO 80 I=1,N                                                        
      K1=IA(I)+1                                                         
      K2=IA(I+1)-1                                                       
      IF (K1.GT.K2) GO TO 80                                             
      MT=K2-K1+1                                                         
      DO 75 K=1,MT                                                       
      DO 75 M=K1,K2                                                      
      IF (JA(M).GT.JA(M-1)) GO TO 75                                     
      J=JA(M-1)                                                          
      JA(M-1)=JA(M)                                                      
      JA(M)=J                                                            
   75 CONTINUE                                                           
   80 CONTINUE                                                           
      M=N                                                                
      DO 95 I=1,M                                                        
      IF (IA(I+1).GT.IA(I)) GO TO 95                                     
      NM=I+1                                                             
      NN=N+1                                                             
      KMIN=IA(NM)                                                        
      KMAX=IA(NN)                                                        
      DO 85 J=KMIN,KMAX                                                  
      KM=KMAX+KMIN-J                                                     
   85 JA(KM)=JA(KM-1)                                                    
      KNOW=IA(I)                                                         
      JA(KNOW)=I                                                         
      DO 90 LL=NM,NN                                                     
   90 IA(LL)=IA(LL)+1                                                    
   95 CONTINUE                                                           
      RETURN                                                             
      END                                                                
!									 
      SUBROUTINE SORDER (N,IA,JA,P,Q,NIMAX,V,L,IER)                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IA(1),JA(1),P(1),Q(1),V(1),L(1)                            
      INTEGER S,SFS,PI,PJ,PK,VI,VJ,VK,QVK,DTHR,DMIN                      
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,        &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,         &
        INTERP, SPARS, NSSFAC, NSCORD
      IER=0                                                              
      DO 5 S=1,NIMAX                                                       
    5 L(S)=S+1                                                           
      SFS=1                                                              
      L(NIMAX)=0                                                           
      DO 10 K=1,N                                                        
      P(K)=K                                                             
      Q(K)=K                                                             
      V(K)=1                                                             
   10 L(K)=0                                                             
      SFS=SFS+N                                                          
      DO 50 K=1,N                                                        
      JMIN=IA(K)                                                         
      JMAX=IA(K+1)-1                                                     
      IF (JMIN.GT.JMAX+1) GO TO 145                                      
      KDIAG=0                                                            
      DO 45 J=JMIN,JMAX                                                  
      VJ=JA(J)                                                           
      IF (VJ.NE.K) GO TO 15                                              
      KDIAG=1                                                            
      GO TO 45                                                           
   15 LLK=K                                                              
   20 LK=LLK                                                             
      LLK=L(LK)                                                          
      IF (LLK.EQ.0) GO TO 25                                             
      IF (V(LLK)-VJ) 20,30,25                                            
   25 LLK=SFS                                                            
      IF (LLK.EQ.0) GO TO 150                                            
      SFS=L(SFS)                                                         
      V(K)=V(K)+1                                                        
      V(LLK)=VJ                                                          
      L(LLK)=L(LK)                                                       
      L(LK)=LLK                                                          
   30 LLK=VJ                                                             
   35 LK=LLK                                                             
      LLK=L(LK)                                                          
      IF (LLK.EQ.0) GO TO 40                                             
      IF (V(LLK)-K) 35,45,40                                             
   40 LLK=SFS                                                            
      IF (LLK.EQ.0) GO TO 150                                            
      SFS=L(SFS)                                                         
      V(VJ)=V(VJ)+1                                                      
      V(LLK)=K                                                           
      L(LLK)=L(LK)                                                       
      L(LK)=LLK                                                          
   45 CONTINUE                                                           
      IF (KDIAG.EQ.0) GO TO 160                                          
   50 CONTINUE                                                           
      J=0                                                                
      DTHR=0                                                             
      DMIN=N                                                             
      I=0                                                                
   55 I=I+1                                                              
      IF (I.GT.N) GO TO 140                                              
      JMIN=MAX0(J+1,I)                                                   
      IF (JMIN.GT.N) GO TO 70                                            
   60 DO 65 J=JMIN,N                                                     
      VI=P(J)                                                            
      IF (V(VI).LE.DTHR) GO TO 75                                        
      IF (V(VI).LT.DMIN) DMIN=V(VI)                                      
   65 CONTINUE                                                           
   70 DTHR=DMIN                                                          
      DMIN=N                                                             
      JMIN=I                                                             
      GO TO 60                                                           
   75 PJ=P(I)                                                            
      P(J)=PJ                                                            
      Q(PJ)=J                                                            
      PI=VI                                                              
      P(I)=PI                                                            
      Q(PI)=I                                                            
      LI=VI                                                              
   80 LI=L(LI)                                                           
      IF (LI.EQ.0) GO TO 105                                             
      VK=V(LI)                                                           
      LLK=VK                                                             
      LJ=VI                                                              
   85 LJ=L(LJ)                                                           
      IF (LJ.EQ.0) GO TO 100                                             
      VJ=V(LJ)                                                           
      IF (VJ.EQ.VK) GO TO 85                                             
   90 LK=LLK                                                             
      LLK=L(LK)                                                          
      IF (LLK.EQ.0) GO TO 95                                             
      IF (V(LLK)-VJ) 90,85,95                                            
   95 LLK=SFS                                                            
      IF (LLK.EQ.0) GO TO 155                                            
      SFS=L(SFS)                                                         
      V(VK)=V(VK)+1                                                      
      V(LLK)=VJ                                                          
      L(LLK)=L(LK)                                                       
      L(LK)=LLK                                                          
      GO TO 85                                                           
  100 IF (V(VK).GT.V(VI)) GO TO 80                                       
      I=I+1                                                              
      QVK=Q(VK)                                                          
      PI=P(I)                                                            
      P(QVK)=PI                                                          
      Q(PI)=QVK                                                          
      P(I)=VK                                                            
      Q(VK)=I                                                            
      GO TO 80                                                           
  105 LI=VI                                                              
  110 IF (L(LI).EQ.0) GO TO 135                                          
      LI=L(LI)                                                           
      VK=V(LI)                                                           
      LLK=VK                                                             
      QVK=MIN0(Q(VK),I)                                                  
  115 LK=LLK                                                             
      LLK=L(LK)                                                          
      IF (LLK.EQ.0) GO TO 120                                            
      VJ=V(LLK)                                                          
      IF (Q(VJ).GT.QVK) GO TO 115                                        
      V(VK)=V(VK)-1                                                      
      L(LK)=L(LLK)                                                       
      L(LLK)=SFS                                                         
      SFS=LLK                                                            
      LLK=LK                                                             
      GO TO 115                                                          
  120 IF (Q(VK).LE.I) GO TO 130                                          
      IF (V(VK).LE.DTHR) GO TO 125                                       
      IF ((DTHR.LT.V(VK)).AND.(V(VK).LT.DMIN)) DMIN=V(VK)                
      GO TO 110                                                          
  125 J=MIN0(Q(VK)-1,J)                                                  
      DTHR=V(VK)                                                         
      GO TO 110                                                          
  130 L(LK)=SFS                                                          
      SFS=L(VK)                                                          
      GO TO 110                                                          
  135 L(LI)=SFS                                                          
      SFS=L(VI)                                                          
      GO TO 55                                                           
  140 RETURN                                                             
  145 CALL YSMER (3HROW,K,13H OF A IS NULL)                              
      GO TO 165                                                          
  150 CALL YSMER (3HROW,K,16H EXCEEDS STORAGE)                           
      GO TO 165                                                          
  155 CALL YSMER (6HVERTEX,VI,16H EXCEEDS STORAGE)                       
      GO TO 165                                                          
  160 CALL YSMER (6HCOLUMN,K,19H.. DIAGONAL MISSING)                     
  165 IER=1                                                              
      RETURN                                                             
      END                                                                


      SUBROUTINE NSSFAC (N,IA,JA,MAXPL,IL,JL,ISL,MAXPU,             &
       IU,JU,ISU,P,V,IRA                                            &
      ,JRA,IRAC,IRL,JRL,IRU,JRU,IER)                                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!      INTEGER IA(1),JA(1),IL(1),JL(1),ISL(1)
!      INTEGER IU(1),JU(1),ISU(1),IRL(1),JRL(1),V(1)

!  NOTE THE FOLLOWING STATEMENT DECLARATION MUST BE CHANGED TO "REAL"
!   INSTEAD OF "INTEGER" ON  IBM, AMDAHL, AND HP COMPUTERS

!      REAL*8 P(1),IRAC(1),IRA(1),JRA(1),IRU(1),JRU(1)
!      REAL*8 VI,VJ,PK,PPK,PI,CEND,REND


      INTEGER IA(1),JA(1),IRA(1),JRA(1),IL(1),JL(1),ISL(1)               
      INTEGER IU(1),JU(1),ISU(1),IRL(1),JRL(1),IRU(1),JRU(1)             
      INTEGER P(1),V(1),IRAC(1)                                          
      INTEGER VI,VJ,VK,PK,PPK,PI,CEND,REND                               
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,        &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,         &
        INTERP, SPARS, SORDER, NSCORD
      IER=0                                                              
      DO 5 K=1,N                                                         
      IRAC(K)=0                                                          
      JRA(K)=0                                                           
      JRL(K)=0                                                           
    5 JRU(K)=0                                                           
      DO 10 K=1,N                                                        
      IAK=IA(K)                                                          
      IF (IAK.GT.IA(K+1)) GO TO 200                                      
      IF (JA(IAK).GT.K) GO TO 10                                         
      JAIAK=JA(IAK)                                                      
      JRA(K)=IRAC(JAIAK)                                                 
      IRAC(JAIAK)=K                                                      
   10 IRA(K)=IAK                                                         
      JLPTR=0                                                            
      IL(1)=1                                                            
      JUPTR=0                                                            
      IU(1)=1                                                            
      DO 195 K=1,N                                                       
      P(1)=1                                                             
      V(1)=N+1                                                           
      LSFS=2                                                             
      VJ=IRAC(K)                                                         
      IF (VJ.EQ.0) GO TO 30                                              
   15 PPK=1                                                              
   20 PK=PPK                                                             
      PPK=P(PK)                                                          
      IF (V(PPK)-VJ) 20,220,25                                           
   25 P(PK)=LSFS                                                         
      V(LSFS)=VJ                                                         
      P(LSFS)=PPK                                                        
      LSFS=LSFS+1                                                        
      VJ=JRA(VJ)                                                         
      IF (VJ.NE.0) GO TO 15                                              
   30 LASTI=0                                                            
      I=K                                                                
   35 I=JRU(I)                                                           
      IF (I.EQ.0) GO TO 60                                               
      PPK=1                                                              
      JMIN=IRL(I)                                                        
      JMAX=ISL(I)+IL(I+1)-IL(I)-1                                        
      IF (LASTI.GT.I) GO TO 40                                           
      LASTI=I                                                            
      LASTID=JMAX-JMIN                                                   
      IF (JL(JMIN).NE.K) LASTID=LASTID+1                                 
   40 IF (JMIN.GT.JMAX) GO TO 35                                         
      DO 55 J=JMIN,JMAX                                                  
      VJ=JL(J)                                                           
   45 PK=PPK                                                             
      PPK=P(PK)                                                          
      IF (V(PPK)-VJ) 45,55,50                                            
   50 P(PK)=LSFS                                                         
      V(LSFS)=VJ                                                         
      P(LSFS)=PPK                                                        
      PPK=LSFS                                                           
      LSFS=LSFS+1                                                        
   55 CONTINUE                                                           
      GO TO 35                                                           
   60 PI=P(1)                                                            
      IF (V(PI).NE.K) GO TO 225                                          
      IF (LASTI.EQ.0) GO TO 65                                           
      IF (LASTID.NE.LSFS-3) GO TO 65                                     
      IRLL=IRL(LASTI)                                                    
      ISL(K)=IRLL+1                                                      
      IF (JL(IRLL).NE.K) ISL(K)=ISL(K)-1                                 
      IL(K+1)=IL(K)+LASTID                                               
      IRL(K)=ISL(K)                                                      
      GO TO 80                                                           
   65 ISL(K)=JLPTR+1                                                     
      PI=P(1)                                                            
      PI=P(PI)                                                           
      VI=V(PI)                                                           
   70 IF (VI.GT.N) GO TO 75                                              
      JLPTR=JLPTR+1                                                      
      IF (JLPTR.GT.MAXPL) GO TO 230                                      
      JL(JLPTR)=VI                                                       
      PI=P(PI)                                                           
      VI=V(PI)                                                           
      GO TO 70                                                           
   75 IRL(K)=ISL(K)                                                      
      IL(K+1)=IL(K)+JLPTR-ISL(K)+1                                       
   80 P(1)=1                                                             
      V(1)=N+1                                                           
      LSFS=2                                                             
      JMIN=IRA(K)                                                        
      JMAX=IA(K+1)-1                                                     
      IF (JMIN.GT.JMAX) GO TO 100                                        
      DO 95 J=JMIN,JMAX                                                  
      VJ=JA(J)                                                           
      PPK=1                                                              
   85 PK=PPK                                                             
      PPK=P(PK)                                                          
      IF (V(PPK)-VJ) 85,205,90                                           
   90 P(PK)=LSFS                                                         
      V(LSFS)=VJ                                                         
      P(LSFS)=PPK                                                        
   95 LSFS=LSFS+1                                                        
  100 LASTI=0                                                            
      I=K                                                                
  105 I=JRL(I)                                                           
      IF (I.EQ.0) GO TO 130                                              
      PPK=1                                                              
      JMIN=IRU(I)                                                        
      JMAX=ISU(I)+IU(I+1)-IU(I)-1                                        
      IF (LASTI.GT.I) GO TO 110                                          
      LASTI=I                                                            
      LASTID=JMAX-JMIN                                                   
      IF (JU(JMIN).NE.K) LASTID=LASTID+1                                 
  110 IF (JMIN.GT.JMAX) GO TO 105                                        
      DO 125 J=JMIN,JMAX                                                 
      VJ=JU(J)                                                           
  115 PK=PPK                                                             
      PPK=P(PK)                                                          
      IF (V(PPK)-VJ) 115,125,120                                         
  120 P(PK)=LSFS                                                         
      V(LSFS)=VJ                                                         
      P(LSFS)=PPK                                                        
      PPK=LSFS                                                           
  125 LSFS=LSFS+1                                                        
      GO TO 105                                                          
  130 PI=P(1)                                                            
      IF (V(PI).NE.K) GO TO 210                                          
      IF (LASTI.EQ.0) GO TO 135                                          
      IF (LASTID.NE.LSFS-3) GO TO 135                                    
      IRUL=IRU(LASTI)                                                    
      ISU(K)=IRUL+1                                                      
      IF (JU(IRUL).NE.K) ISU(K)=ISU(K)-1                                 
      IU(K+1)=IU(K)+LASTID                                               
      IRU(K)=ISU(K)                                                      
      GO TO 150                                                          
  135 ISU(K)=JUPTR+1                                                     
      PI=P(1)                                                            
      PI=P(PI)                                                           
      VI=V(PI)                                                           
  140 IF (VI.GT.N) GO TO 145                                             
      JUPTR=JUPTR+1                                                      
      IF (JUPTR.GT.MAXPU) GO TO 215                                      
      JU(JUPTR)=VI                                                       
      PI=P(PI)                                                           
      VI=V(PI)                                                           
      GO TO 140                                                          
  145 IRU(K)=ISU(K)                                                      
      IU(K+1)=IU(K)+JUPTR-ISU(K)+1                                       
  150 I=K                                                                
  155 I1=JRL(I)                                                          
      CEND=ISL(I)+IL(I+1)-IL(I)                                          
      IF (IRL(I).GE.CEND) GO TO 160                                      
      IRLI=IRL(I)                                                        
      J=JL(IRLI)                                                         
      JRL(I)=JRL(J)                                                      
      JRL(J)=I                                                           
  160 I=I1                                                               
      IF (I.EQ.0) GO TO 165                                              
      IRL(I)=IRL(I)+1                                                    
      GO TO 155                                                          
  165 I=K                                                                
  170 I1=JRU(I)                                                          
      REND=ISU(I)+IU(I+1)-IU(I)                                          
      IF (IRU(I).GE.REND) GO TO 175                                      
      IRUI=IRU(I)                                                        
      J=JU(IRUI)                                                         
      JRU(I)=JRU(J)                                                      
      JRU(J)=I                                                           
  175 I=I1                                                               
      IF (I.EQ.0) GO TO 180                                              
      IRU(I)=IRU(I)+1                                                    
      GO TO 170                                                          
  180 I=IRAC(K)                                                          
      IF (I.EQ.0) GO TO 195                                              
  185 I1=JRA(I)                                                          
      IRA(I)=IRA(I)+1                                                    
      IF (IRA(I).GE.IA(I+1)) GO TO 190                                   
      IRAI=IRA(I)                                                        
      IF (JA(IRAI).GT.I) GO TO 190                                       
      JAIRAI=JA(IRAI)                                                    
      JRA(I)=IRAC(JAIRAI)                                                
      IRAC(JAIRAI)=I                                                     
  190 I=I1                                                               
      IF (I.NE.0) GO TO 185                                              
  195 CONTINUE                                                           
      ISL(N)=JLPTR                                                       
      ISU(N)=JUPTR                                                       
      RETURN                                                             
  200 CALL YSMER (3HROW,K,13H OF A IS NULL)                              
      GO TO 235                                                          
  205 CALL YSMER (3HROW,K,20H HAS DUPLICATE ENTRY)                       
      GO TO 235                                                          
  210 CALL YSMER (3HROW,K,17H HAS A NULL PIVOT)                          
      GO TO 235                                                          
  215 CALL YSMER (3HROW,K,19H EXCEEDS JU STORAGE)                        
      GO TO 235                                                          
  220 CALL YSMER (3HCOL,K,20H HAS DUPLICATE ENTRY)                       
      GO TO 235                                                          
  225 CALL YSMER (3HCOL,K,17H HAS A NULL PIVOT)                          
      GO TO 235                                                          
  230 CALL YSMER (3HCOL,K,19H EXCEEDS JL STORAGE)                        
  235 IER=1                                                              
      RETURN                                                             
      END                                                                


      SUBROUTINE NSCORD (N,IA,JA,IAWORK,JAWORK,C,IR,ICT)                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IA(1),JA(1),IAWORK(1),JAWORK(1),C(1),IR(1),ICT(1)          
      EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,        &
        SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,         &
        INTERP, SPARS, SORDER, NSSFAC
      DO 5 I=1,N                                                         
    5 ICT(I)=0                                                           
      IAWORK(1)=1                                                        
      DO 15 K=1,N                                                        
      ICK=C(K)                                                           
      JMIN=IAWORK(K)                                                     
      JMAX=JMIN+IA(ICK+1)-IA(ICK)-1                                      
      IAWORK(K+1)=JMAX+1                                                 
      IF (JMIN.GT.JMAX) GO TO 15                                         
      IAINK=IA(ICK)-1                                                    
      DO 10 J=JMIN,JMAX                                                  
      IAINK=IAINK+1                                                      
      JAOUTJ=JA(IAINK)                                                   
      JAOUTJ=IR(JAOUTJ)                                                  
      JAWORK(J)=JAOUTJ                                                   
   10 ICT(JAOUTJ)=ICT(JAOUTJ)+1                                          
   15 CONTINUE                                                           
      IA(1)=1                                                            
      DO 20 I=1,N                                                        
      IA(I+1)=IA(I)+ICT(I)                                               
   20 ICT(I)=IA(I)                                                       
      DO 30 I=1,N                                                        
      JMIN=IAWORK(I)                                                     
      JMAX=IAWORK(I+1)-1                                                 
      IF (JMIN.GT.JMAX) GO TO 30                                         
      DO 25 J=JMIN,JMAX                                                  
      JAOUTJ=JAWORK(J)                                                   
      ICTJ=ICT(JAOUTJ)                                                   
      JA(ICTJ)=I                                                         
   25 ICT(JAOUTJ)=ICTJ+1                                                 
   30 CONTINUE                                                           
      RETURN                                                             
      END                                                                

!	
      BLOCK DATA ALPHA1                                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /ALPHA/ IGO(4),IBLANK,MBLANK,JINTER                         
      COMMON /APLOT/ JVERT(52,2),JBLANK,JSTAR,JPLUS,JBAR                 
      COMMON /GEAR1/ T,GUESS,HMIN,HMAX,EPS1,UROUND,NC,MF1,KFLAG1,JSTART  
      COMMON /HEAT/ CV,Q(200),SC(200,7),ISC(200,3),ITEMP                 
      COMMON /STORE/ AST(35),IPL(7),TEMEND,NTEMP,TMI,NPHOT,PHI,IL,NFRST,   &
      IPH(30),QM(100),PM(100),PSTOP                                      
      COMMON /PHOTR/ RDAT(80),RTIM(80),RR1,IN10,IPP
      INTEGER PHI,TMI,AST                                              
!       EXTERNAL YFIX, RATES, DIFFUN, MATRX, XPLOT, VALU, CONVT, SAVLIN,
!     &  SPLNA, DRIVES, STIFFS, COSET, NSCORA, NSNFAC, NSBSLV, YSMER,
!     &  INTERP, SPARS, SORDER, NSSFAC, NSCORD
      DATA IGO(1)/4HMORE/,IGO(2)/4HCONT/,IGO(3)/4HPLOT/,IGO(4)/4H    /,IBLANK/4H    /,   &
      MBLANK/4HM   /,JINTER/4HINTV/                        
!      DATA JVERT/12*4H ,' C  ',' O  ',' N  ',' C  ',' E  ',' N  ','     
!     &T  ',' R  ',' A  ',' T  ',' I  ',' O  ',' N  ','    ',' P  ',' P
!     & ',' M  ', 75*4H   /
      DATA JBLANK/4H    /,JSTAR/1H*/,JPLUS/1H+/,JBAR/1HI/                
      DATA ITEMP/4HTEMP/,TEMEND/0.0/,PSTOP/0.0/,IN10/1/                  
      DATA UROUND/2.3E-16/                                               
      END                                                                
