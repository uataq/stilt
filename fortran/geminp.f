!###############################################################################
! GEMINP - Meteorology Input for one time period 
!-------------------------------------------------------------------------------
! LAST REVISED: 28 May 2008 (RRD) - initial version
!               01 Jul 2008 (RRD) - multiple time periods
!               17 Jul 2008 (RRD) - unit number reset with new file  
!-------------------------------------------------------------------------------

SUBROUTINE geminp (jet,back,ntim)

  USE gemkon   
  USE gemvar  
  USE gemcfg
  USE funits
  use module_defgrid

  IMPLICIT NONE


  INTEGER*4,     INTENT(IN) :: jet  
  LOGICAL,       INTENT(IN) :: back
  INTEGER*4,     INTENT(IN) :: ntim 

  CHARACTER(4)  :: VARB
  CHARACTER(50) :: LABEL
  INTEGER       :: MYR,MMO,MDA,MHR,MFH,LEV,KGRD
  INTEGER       :: K,L,NX,NY,NZ,NP,NEXP,KRET,JREC,KSUM
  REAL          :: VAR1,PREC

  INTEGER       :: KT  = 1     ! time period 
  INTEGER       :: MET = 0
  INTEGER       :: NUM = 6     ! default bucket rate

  CHARACTER(1), ALLOCATABLE :: BUF(:)
  REAL,         ALLOCATABLE :: VAR(:,:)

  SAVE buf,var,met,num,kt

  INTERFACE
  SUBROUTINE PAKINP(RVAR,CVAR,NX,NY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)
  REAL,          INTENT(OUT)   :: rvar (:,:)     ! real data unpacked
  CHARACTER(1),  INTENT(IN)    :: cvar (:)       ! packed input of NX*NY
  INTEGER,       INTENT(IN)    :: nx,ny          ! size of input array
  INTEGER,       INTENT(IN)    :: nx1,ny1        ! optional sub-grid left edge
  INTEGER,       INTENT(IN)    :: lx,ly          ! length of sub-grid
  REAL,          INTENT(IN)    :: prec           ! precision of packed data
  INTEGER,       INTENT(IN)    :: nexp           ! packing scaling exponent
  REAL,          INTENT(IN)    :: var1           ! value of array at (LX,LY)
  INTEGER,       INTENT(INOUT) :: ksum           ! rotating checksum
  END SUBROUTINE PAKINP
  END INTERFACE

  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  IF(KINIT.LT.0)RETURN

  ksum=1
  IF(.NOT.ALLOCATED(buf)) ALLOCATE (BUF(nx*ny))
  IF(.NOT.ALLOCATED(var)) ALLOCATE (VAR(nx,ny))

! check to determine if this is a data input time
  IF(MET.NE.0)THEN
     IF(BACK)THEN
        MET=JET-DREC(KGRD,KT)%DELTA/2
     ELSE
        MET=JET+DREC(KGRD,KT)%DELTA/2
     END IF
     IF(MOD(MET,DREC(KGRD,KT)%DELTA).NE.0)RETURN 
  ELSE
!    first time load data regardless of the time
     IF(BACK)THEN
        MET=FILE(KGRD,KT)%LAST%MACC-INT((FILE(KGRD,KT)%LAST%MACC-JET)/    &
            DREC(KGRD,KT)%DELTA)*DREC(KGRD,KT)%DELTA
     ELSE
        MET=FILE(KGRD,KT)%FIRST%MACC+INT((JET-FILE(KGRD,KT)%FIRST%MACC)/  &
            DREC(KGRD,KT)%DELTA)*DREC(KGRD,KT)%DELTA
     END IF
  END IF

! compute record numbers to index record for indicated time
  JREC=DREC(KGRD,KT)%REC_PER*(MET-FILE(KGRD,KT)%FIRST%MACC)/   &
       DREC(KGRD,KT)%DELTA+1

! hold read if model run starts before global data
  IF(JREC.LT.1)RETURN

! check for end of data
  IF(JREC.GT.FILE(KGRD,KT)%ENDREC)THEN
     IF(BACK)THEN
        KT=KT-1
        IF(KT.GE.1)THEN
           JREC=DREC(KGRD,KT)%REC_PER*(MET-FILE(KGRD,KT)%FIRST%MACC)/ &
                DREC(KGRD,KT)%DELTA+1
        ELSE
           WRITE(KF21,*)'*ERROR* geminp: reading before first meteo file'
           STOP 800
        END IF
     ELSE
        KT=KT+1
        IF(KT.LE.NTIM)THEN
           JREC=DREC(KGRD,KT)%REC_PER*(MET-FILE(KGRD,KT)%FIRST%MACC)/ &
                DREC(KGRD,KT)%DELTA+1
        ELSE
           WRITE(KF21,*)'*ERROR* geminp: reading after last meteo file'
           STOP 800
        END IF
     END IF
     KUMET=FILE(KGRD,KT)%KUNIT
     WRITE(KF21,*)' NOTICE geminp: starting new meteo file (unit) - ',KUMET 
  END IF

! read the index record for testing (krec always points to index)
  READ(KUMET,REC=JREC,IOSTAT=kret) LABEL 

  IF(kret.NE.0)THEN
     WRITE(KF21,*)'*ERROR* geminp: problem reading index record - ',jrec
     STOP 800
  ELSE
     READ(LABEL,'(6I2,2X,A4)') MYR,MMO,MDA,MHR,MFH,LEV,VARB
     WRITE(KF21,*)' NOTICE geminp: loaded new data (mo,da,hr) - ',MMO,MDA,MHR 
     IF(VARB.NE.'INDX')THEN
        WRITE(KF21,*)'*ERROR* geminp: index record not found - ',jrec
        WRITE(KF21,*)' ',LABEL
        STOP 800
     END IF
  END IF
  jrec=1+jrec ! skip past the index record

! zero out variables that may not be complete in the vertical
  WWW=0.0
  MMM=0.0
    
  DO L=0,nz
  DO K=1,nvar(L)

!    read one data record (variable at a level)
     READ(KUMET,REC=JREC)LABEL,BUF

!    decode index record for packing information
     READ(LABEL,'(6I2,2X,A4,I4,2E14.7)',IOSTAT=kret)       &
          MYR,MMO,MDA,MHR,MFH,LEV,VARB,NEXP,PREC,VAR1
     CALL PAKINP(VAR,BUF,NX,NY,1,1,NX,NY,PREC,NEXP,VAR1,KSUM)

     IF(kret.NE.0)THEN
!       error decoding label field
        WRITE(KF21,*)'*ERROR* geminp: decoding meteo label field'
        WRITE(KF21,*)' ',LABEL
        STOP 800
     END IF

     IF(mfh.LT.0.OR.varb.EQ.'NULL')THEN
!       missing data not permitted
        WRITE(KF21,*)'*ERROR* geminp: missing meteorological data'
        WRITE(KF21,*)' ',LABEL
        STOP 800
     END IF

     IF(lev.NE.L)THEN
!       level mismatch
        WRITE(KF21,*)'*ERROR* geminp: input data level mismatch'
        WRITE(KF21,*)' Prog: ',L,'  File: ',lev,'  Record: ',jrec 
        STOP 800
     END IF

     IF(L.GT.0)THEN
!       upper level variables
        IF(VARB.EQ.'HGTS') HHH(:,:,L)=VAR
        IF(VARB.EQ.'TEMP') TTT(:,:,L)=VAR
        IF(VARB.EQ.'UWND') UUU(:,:,L)=VAR
        IF(VARB.EQ.'VWND') VVV(:,:,L)=VAR
        IF(VARB.EQ.'WWND') WWW(:,:,L)=VAR
        IF(VARB.EQ.'RELH') MMM(:,:,L)=VAR

     ELSE
!       surface variables
        IF(VARB.EQ.'PRSS') SFC=VAR
        IF(VARB.EQ.'U10M') U10=VAR
        IF(VARB.EQ.'V10M') V10=VAR
        IF(VARB.EQ.'T02M') T02=VAR

        IF(VARB.EQ.'TPP6') THEN
!          convert accumulated totals into a rate (m/min)
!          assume global data all use a 6 hr bucket (360 min)
           IF(mhr.EQ.3.OR.mhr.EQ.9.OR.mhr.EQ.15.OR.mhr.EQ.21)THEN
              IF(num.EQ.6)num=3  
!             intermediate times (bucket reset at 0,6,12,18)
              P6H=VAR/180.0
           ELSE
              IF(num.EQ.6)THEN
!                simple rate conversion
                 P6H=VAR/360.0
              ELSE
!                subtract previous rate for new rate
                 P6H=VAR/180.0-P6H
              END IF
           END IF
        END IF
     END IF
     jrec=jrec+1

  END DO
  END DO

  CALL GEMDAT
  CALL GEMDIV
  CALL GEMSTB

END SUBROUTINE geminp
