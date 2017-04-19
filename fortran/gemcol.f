!###############################################################################
! GEMCOL - Initialize packed vertical column integrated concentration output 
! file in standard HYSPLIT binary format. Routine supports kzavg=2  
!-------------------------------------------------------------------------------
! KZAVG  :  OUTPUT TYPE
! 0      =  one or more individual levels
! 1      =  single layer averaged concentration
! 2      =  single full-column integral
! 3      =  kzavg=2 + kzavg=0
! 4      =  kzavg=2 + kzavg=1
!-------------------------------------------------------------------------------
! LAST REVISED: 20 May 2008 (RRD) - initial version
!               28 Jul 2008 (AS)  - vertically integrated output
!-------------------------------------------------------------------------------

SUBROUTINE gemcol (iy,im,id,ih,plat,plon,plvl)

  USE gemcfg
  USE gemkon   
  USE gemvar  
  USE funits

  IMPLICIT NONE

  INTEGER*4,    INTENT(IN) :: iy,im,id,ih
  REAL*4,       INTENT(IN) :: plat,plon,plvl

  LOGICAL       :: ftest
  INTEGER*4     :: i,j,k,n,kk,nx,ny,nz,np,kgrd,conlev
  REAL*4        :: clat1,clon1,dlat,dlon
  REAL*4        :: depth
  CHARACTER*80  :: label

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

  SAVE conlev

  IF(KINIT.LT.0)RETURN
  IF(KZAVG.LT.2)RETURN

!-------------------------------------------------
! output concentration to file if exists and is open

  INQUIRE(FILE=GEMZINT,EXIST=FTEST)
  IF(FTEST)THEN
     INQUIRE(FILE=GEMZINT,OPENED=FTEST)
     IF(FTEST)THEN

!       at appropriate interval
        IF(wfreq.EQ.0) RETURN 
        IF(ihour.GT.0) THEN  
           IF(ih.LT.ihour) RETURN
           ihour=0
        END IF
        IF(MOD(ih,wfreq).NE.0) RETURN  

!       exists and opened means write data to file
        WRITE(KF29)IY,IM,ID,IH,0,0
        WRITE(KF29)IY,IM,ID,IH,0,0

        DO N=1,NP

!          convert to output units at STP
           DO J=1,NY
           DO I=1,NX

              CONC(I,J)=0.0
              DO kk=1,nz
                 IF(kk.LT.nz)THEN
                    depth=MAX(DZMIN,HHH(i,j,kk+1)-HHH(i,j,kk))
                 ELSE
                    depth=MAX(DZMIN,HHH(i,j,kk)-HHH(i,j,kk-1))
                 END IF
                 CONC(I,J)=CONC(I,J)+REAL(XXX(i,j,kk,n))*CFACT*RRR(i,j,kk)*DEPTH 
              END DO

           END DO
           END DO

!          write all grid points 
           WRITE(LABEL,'(A1,I3.3)')'P',N
           WRITE(KF29)LABEL(1:4),CONLEV,((CONC(I,J),I=1,NX),J=1,NY)
        END DO

!       diagnostic message
        WRITE(KF21,*)' NOTICE gemcol: (da,hr) - ',id,ih

        RETURN
     END IF
  END IF
     
!--------------------------------------------------
! file doesn't exist or is not opened, write headers

  OPEN(KF29,FILE=GEMZINT,FORM='UNFORMATTED')
  WRITE(KF29)'GEMx',IY,IM,ID,IH,    0,   1,   0
  WRITE(KF29)IY,IM,ID,IH, PLAT, PLON, PLVL
  WRITE(KF29)NY,NX,DLAT,DLON,CLAT1,CLON1
  conlev=SUM(hhh(:,:,nz))/nx/ny
  WRITE(KF29) 1, CONLEV  
  DO n=1,np
     WRITE(LABEL((N-1)*4+1:),'(A1,I3.3)')'P',N
  END DO
  WRITE(KF29)NP,LABEL(1:NP*4)

END SUBROUTINE gemcol 
