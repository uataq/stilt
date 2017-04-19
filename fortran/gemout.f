!###############################################################################
! GEMOUT - Initialize packed concentration output file in standard HYSPLIT
! binary format. Horizontal grid is the same as the meteorological grid with the
! exception that the cells are centered over the meteo grid points and the pole 
! grid cell is circular about the pole. Routine supports kzavg=0 or kzavg=1 
!-------------------------------------------------------------------------------
! KZAVG  :  OUTPUT TYPE
! 0      =  one or more individual levels
! 1      =  single layer averaged concentration
! 2      =  single full-column integral
! 3      =  kzavg=2 + kzavg=0
! 4      =  kzavg=2 + kzavg=1
!-------------------------------------------------------------------------------
! LAST REVISED: 20 May 2008 (RRD) - initial version
!               28 Jul 2008 (RRD) - layer averaged concentrations
!-------------------------------------------------------------------------------

SUBROUTINE gemout (iy,im,id,ih,plat,plon,plvl)

  USE gemcfg
  USE gemkon   
  USE gemvar  
  USE funits

  IMPLICIT NONE

  INTEGER*4,    INTENT(IN) :: iy,im,id,ih
  REAL*4,       INTENT(IN) :: plat,plon,plvl

  INTEGER*4,   ALLOCATABLE :: conlev(:)

  LOGICAL       :: ftest
  INTEGER*4     :: i,j,k,n,kk,nx,ny,nz,np,kgrd,level
  REAL*4        :: clat1,clon1,dlat,dlon
  REAL*4        :: depth,zconc
  CHARACTER*80  :: label

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  SAVE conlev,level

  IF(KINIT.LT.0)RETURN
  IF(KZAVG.EQ.2)RETURN

!-------------------------------------------------
! output concentration to file if exists and is open

  INQUIRE(FILE=GEMCONC,EXIST=FTEST)
  IF(FTEST)THEN
     INQUIRE(FILE=GEMCONC,OPENED=FTEST)
     IF(FTEST)THEN

!       at appropriate interval
        IF(wfreq.EQ.0) RETURN 
        IF(ihour.GT.0) THEN  
           IF(ih.LT.ihour) RETURN
           ihour=0
        END IF
        IF(MOD(ih,wfreq).NE.0) RETURN  

!       exists and opened means write data to file
        WRITE(KF27)IY,IM,ID,IH,0,0
        WRITE(KF27)IY,IM,ID,IH,0,0

        DO N=1,NP
        DO K=1,LEVEL

!          convert to output units at STP
           DO J=1,NY
           DO I=1,NX

              IF(kzavg.EQ.1.OR.kzavg.EQ.4)THEN
                 ZCONC=0.0
                 CONC(I,J)=0.0
                 DO kk=kzbot,kztop
                    IF(kk.LT.nz)THEN
                       depth=MAX(DZMIN,HHH(i,j,kk+1)-HHH(i,j,kk))
                    ELSE
                       depth=MAX(DZMIN,HHH(i,j,kk)-HHH(i,j,kk-1))
                    END IF
                    ZCONC=ZCONC+DEPTH
                    CONC(I,J)=CONC(I,J)+REAL(XXX(i,j,kk,n))*CFACT*ROW*DEPTH 
                 END DO
                 CONC(I,J)=CONC(I,J)/ZCONC 

              ELSE
                 kk=k+kzbot-1
                 IF(kk.GT.0)THEN
                    CONC(I,J) = REAL(XXX(i,j,kk,n))*CFACT*ROW
                 ELSE
                    CONC(I,J) = REAL(DDD(i,j,n))*CFACT*ROW
                 END IF
              END IF

           END DO
           END DO

!          write all grid points 
           WRITE(LABEL,'(A1,I3.3)')'P',N
           WRITE(KF27)LABEL(1:4),CONLEV(K),((CONC(I,J),I=1,NX),J=1,NY)
        END DO
        END DO

!       diagnostic message
        WRITE(KF21,*)' NOTICE gemout: (da,hr,mass) - ',id,ih,REAL(qtot)

        RETURN
     END IF
  END IF
     
!--------------------------------------------------
! check initial output and frequency values for consistency

  ihour = MAX(0, MIN(23,ihour))
  IF (wfreq.GT.24) wfreq=MOD(wfreq-1,24)+1

!--------------------------------------------------
! file doesn't exist or is not opened, write headers

! two-dimensional concentration output file
  OPEN(KF27,FILE=GEMCONC,FORM='UNFORMATTED')

! Global Model,      Start date, Fcst Srcs Pack
  WRITE(KF27)'GEMx',IY,IM,ID,IH,    0,   1,   0

! calculation start ...   lat   long  level min
  WRITE(KF27)IY,IM,ID,IH, PLAT, PLON, PLVL,   0

! horizontal grid 
  WRITE(KF27)NY,NX,DLAT,DLON,CLAT1,CLON1

! vertical grid
  IF(kzavg.EQ.1.OR.kzavg.EQ.4)THEN
     level=1
  ELSE
     level=kztop-kzbot+1
  END IF

  ALLOCATE (conlev(level))

  IF(kzavg.EQ.1.OR.kzavg.EQ.4)THEN
     conlev(level)=SUM(hhh(:,:,kztop))/nx/ny
  ELSE
     DO k=1,level 
        kk=k+kzbot-1
        IF(kk.LE.0)THEN
           conlev(k)=0
        ELSE
           conlev(k)=SUM(hhh(:,:,kk))/nx/ny
        END IF
     END DO
  END IF
  WRITE(KF27) LEVEL, CONLEV  

! pollutant identification record 
  DO n=1,np
     WRITE(LABEL((N-1)*4+1:),'(A1,I3.3)')'P',N
  END DO
  WRITE(KF27)NP,LABEL(1:NP*4)

END SUBROUTINE gemout 
