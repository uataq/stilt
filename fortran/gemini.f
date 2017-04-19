!###############################################################################
! GEMINI - Set the intial concentration values according to several methods.
! =0 no initialization all zero     =1 by latitude bands from station.txt
! =2 pos interpolation station.txt  =3 gem 3D dump file gbldump.bin
! =4 gblmdl37 file (no deposition)  >4 single value CKON from namelist
!-------------------------------------------------------------------------------
! LAST REVISED: 20 May 2008 (RRD) - initial version
!               27 Jun 2008 (RRD) - 3D start without deposition field 
!               02 Jul 2008 (RRD) - eliminate normalization
!               09 Sep 2008 (RRD) - half grid point correction
!-------------------------------------------------------------------------------

SUBROUTINE gemini (IYR,IMO,IDA)

  USE gemcfg
  USE gemkon  
  USE gemvar
  USE funits

  IMPLICIT NONE

  INTEGER*4, INTENT(IN)  :: IYR,IMO,IDA      ! initialization time

  REAL*4,    ALLOCATABLE :: xcon(:)          ! concentration sampling data
  REAL*4,    ALLOCATABLE :: xvar(:), yvar(:) ! working array for regression
  INTEGER*4, ALLOCATABLE :: clat(:), clon(:) ! sampler positions

  LOGICAL   :: poslon = .true.               ! all positive longitudes
  LOGICAL   :: blend  = .false.              ! set the data blending flag

  REAL*8    :: height,conval                 ! binary input file values
  REAL*4    :: es,ea,dp                      ! dewpoint computation variables
  REAL*4    :: yintc_nh,slope_nh             ! linear regression
  REAL*4    :: yintc_sh,slope_sh
  REAL*4    :: cini,cini_nh,cini_sh
  REAL*4    :: xlat,xlon                     ! temporary position

  INTEGER*4 :: i,j,k,l,ii,jj,kk              ! grid indicies
  INTEGER*4 :: kret,numb                     ! return code and number of input records
  INTEGER*4 :: kyr,kmo,kda,khr,kmn           ! input data date

  INTEGER*4 :: nxg,nyg,nzg,npg               ! dimensions of 3D input grid
  REAL*4    :: clatg,clong,glat,glon         ! grid corners and spacing
  REAL*4,   ALLOCATABLE :: ggg(:)            ! vertical levels for input data
  REAL*8,   ALLOCATABLE :: gqc(:,:,:,:)      ! temp 3D mixing ratio array

  INTEGER*4 :: nx,ny,nz,np,kgrd              ! dimensions of current simulation grid
  REAL*4    :: clat1,clon1,dlat,dlon

  REAL*4            :: pres
  REAL*4, PARAMETER :: cbot = 0.87           ! lower stratosphere mass correction (0.87)
  REAL*4, PARAMETER :: ctop = 0.54           ! upper stratosphere mass correction (0.54)

                                             ! future versions may support multiple
  INTEGER*4 :: mp = 1                        ! initialization pollutant number

  COMMON /GEMGRID/ clat1,clon1,dlat,dlon
  COMMON /GEMDIMS/ nx,ny,nz,np,kgrd

!------------------------------------------------

  INTERFACE
    SUBROUTINE linreg (xvar,yvar,slope,yintc)
    REAL*4, INTENT(IN)  :: xvar(:), yvar(:)
    REAL*4, INTENT(OUT) :: slope, yintc
    END SUBROUTINE linreg
  END INTERFACE

  IF(KINIT.LT.0)RETURN
  XXX = 0.0
  DDD = 0.0

!-------------------------------------------------
! concentration dump file

  OPEN(KF28,FILE=GEMDUMP,FORM='UNFORMATTED')

!---------------------------------------------------
! read station file for initialization

  IF(KINIT.EQ.1.OR.KINIT.EQ.2)THEN

     OPEN (KF26,FILE='startup.txt')
     kret=0
     numb=0
     DO WHILE (kret.EQ.0)
        READ(KF26,*,IOSTAT=KRET) glat, glon
        IF(kret.EQ.0) numb=numb+1
     END DO
     REWIND(KF26)
     ALLOCATE (xcon(numb),clat(numb),clon(numb))

!    set the longitude system flag
     IF(clon1.LT.0.0)poslon=.false.

     jj=0
     kk=0
!    concentration always the dependent variable
     DO L=1,numb
        READ(KF26,*)CLAT(L),CLON(L),XCON(L)
        IF(clat(L).LE.0)jj=jj+1 
        IF(clat(L).GT.0)kk=kk+1

!       adjust the longitude system to reflect meteo data
        IF(poslon)THEN
           IF(clon(L).LT.0)clon(L)=clon(L)+360
        ELSE
           IF(clon(L).GT.180)clon(L)=360-clon(L)
        END IF

     END DO
     CLOSE(KF26)

     WRITE(KF21,*)' NOTICE gemini: station initialization from file - ',numb
     IF(jj.LT.2.OR.kk.LT.2)THEN
        WRITE(KF21,*)'Insufficient number for initialization!'
        WRITE(KF21,*)'Northern hemisphere: ',kk
        WRITE(KF21,*)'Southern hemisphere: ',jj
        STOP 800
     END IF
  END IF

!---------------------------------------------------
! no initialization

  IF(KINIT.EQ.0)THEN
     WRITE(KF21,*)' NOTICE gemini: all zero initialization for concentration'

!---------------------------------------------
! initialization from station file by latitude

  ELSEIF(KINIT.EQ.1)THEN

     ALLOCATE (xvar(kk),yvar(kk))
!    northern hemisphere regression
     kk=0
     DO L=1,numb
        IF(clat(L).GT.0)THEN
           kk=kk+1
           xvar(kk)=clat(L)
           yvar(kk)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_nh,yintc_nh)
     DEALLOCATE (xvar,yvar)

     ALLOCATE (xvar(jj),yvar(jj))
!    southern hemisphere regression
     jj=0
     DO L=1,numb
        IF(clat(L).LE.0)THEN
           jj=jj+1
           xvar(jj)=clat(L)
           yvar(jj)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_sh,yintc_sh)
     DEALLOCATE (xvar,yvar)

     DO j=1,ny
        xlat=(j-1)*dlat+clat1-dlat/2.0
        IF(xlat.GT.0.0)THEN
           CINI=MIN(CMAX,MAX(CMIN,yintc_nh+slope_nh*xlat))/CFACT/ROW
        ELSE
           CINI=MIN(CMAX,MAX(CMIN,yintc_sh+slope_sh*xlat))/CFACT/ROW
        END IF 
        DO i=1,nx
           DO k=1,nz
              PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP
              IF(pres.GT.150.0)THEN
                 XXX(i,j,k,mp)=CINI
              ELSEIF(pres.GT.70.0)THEN
                 XXX(i,j,k,mp)=CBOT*CINI
              ELSE
                 XXX(i,j,k,mp)=CTOP*CINI
              END IF
           END DO
        END DO
     END DO
     DEALLOCATE (clat,clon,xcon)
     WRITE(KF21,*)' NOTICE gemini: latitudinal initialization'

!----------------------------------------------
! dewpoint interpolation from station file

  ELSEIF(KINIT.EQ.2)THEN

     ALLOCATE (xvar(kk),yvar(kk))
!    northern hemisphere regression
     kk=0
     DO L=1,numb
        IF(clat(L).GT.0)THEN
           kk=kk+1
           j=1+(clat(L)-clat1+dlat/2.0)/dlat
           i=1+(clon(L)-clon1+dlon/2.0)/dlon
           IF(i.LT. 1) i=nx+i
           IF(i.GT.nx) i=i-nx
           ES=10.0**(9.4051-2353.0/TTT(I,J,1))
           EA=MAX(DBLE(1.0),MMM(I,J,1))*ES/100.0
           xvar(kk)=2353.0/(9.4051-ALOG10(EA))
           yvar(kk)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_nh,yintc_nh)
     DEALLOCATE (xvar,yvar)

     ALLOCATE (xvar(jj),yvar(jj))
!    southern hemisphere regression
     jj=0
     DO L=1,numb
        IF(clat(L).LE.0)THEN
           jj=jj+1
           j=1+(clat(L)-clat1+dlat/2.0)/dlat
           i=1+(clon(L)-clon1+dlon/2.0)/dlon
           IF(i.LT. 1) i=nx+i
           IF(i.GT.nx) i=i-nx
           ES=10.0**(9.4051-2353.0/TTT(I,J,1))
           EA=MAX(DBLE(1.0),MMM(I,J,1))*ES/100.0
           xvar(jj)=2353.0/(9.4051-ALOG10(EA))
           yvar(jj)=xcon(L)
        END IF
     END DO
     CALL linreg (xvar,yvar,slope_sh,yintc_sh)
     DEALLOCATE (xvar,yvar)

     DO J=1,ny
     DO I=1,nx

        ES=10.0**(9.4051-2353.0/TTT(I,J,1))
        EA=MAX(DBLE(1.0),MMM(I,J,1))*ES/100.0
        DP=2353.0/(9.4051-ALOG10(EA))

        xlat=(j-1)*dlat+clat1-dlat/2.0
        IF(xlat.GE.5.0)THEN
           CINI=MIN(CMAX,MAX(CMIN,yintc_nh+slope_nh*dp))/CFACT/ROW
        ELSEIF(xlat.LT.5.0.AND.xlat.GT.-5.0)THEN
           CINI_NH=MIN(CMAX,MAX(CMIN,yintc_nh+slope_nh*dp))/CFACT/ROW
           CINI_SH=MIN(CMAX,MAX(CMIN,yintc_sh+slope_sh*dp))/CFACT/ROW
           CINI=0.5*(CINI_NH+CINI_SH)
        ELSE
           CINI=MIN(CMAX,MAX(CMIN,yintc_sh+slope_sh*dp))/CFACT/ROW
        END IF 

        DO k=1,nz
           PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP
           IF(pres.GT.150.0)THEN
              XXX(i,j,k,mp)=CINI
           ELSEIF(pres.GT.70.0)THEN
              XXX(i,j,k,mp)=CBOT*CINI
           ELSE
              XXX(i,j,k,mp)=CTOP*CINI
           END IF
        END DO
		
     END DO
     END DO
     DEALLOCATE (clat,clon,xcon)
     WRITE(KF21,*)' NOTICE gemini: dew point initialization'

!----------------------------------------------
! initialization from previous simulation

  ELSEIF(KINIT.EQ.3.OR.KINIT.EQ.4)THEN

     READ(KF28) KYR,KMO,KDA,KHR,KMN
     IF(KINIT.EQ.3)THEN
        READ(KF28) NXG,NYG,NZG,NPG
     ELSE
        READ(KF28) NXG,NYG,NZG
        NPG=1
     END IF
     READ(KF28) CLATG,CLONG,GLAT,GLON
     ALLOCATE (GGG(nzg))
     READ(KF28) GGG

     IF(IYR.NE.KYR.OR.IMO.NE.KMO.OR.IDA.NE.KDA)THEN
        WRITE(KF21,*)'*ERROR* gemini: GBL3DM.BIN header inconsistent with start date'
        WRITE(KF21,*)' FILE date: ',KYR,KMO,KDA
        WRITE(KF21,*)' PROG date: ',IYR,IMO,IDA
        STOP 800
     END IF

!    if grids are not the same, set the data blending flag
     IF(nxg.NE.nx.OR.nyg.NE.ny.OR.nzg.NE.nz.OR.npg.NE.np) THEN
        BLEND=.true.
        WRITE(KF21,*)' NOTICE gemini: input data grid dimensions: ',nxg,nyg,nzg,npg
        WRITE(KF21,*)'                model data grid dimensions: ',nx ,ny ,nz ,np 
        WRITE(KF21,*)'  Horizontal domain limits: ',clatg,clong,glat,glon
        WRITE(KF21,*)'  Vertical levels (1 to 6): ',ggg(1:6)
     END IF

     IF(blend)THEN
        ALLOCATE (gqc(nxg,nyg,nzg,npg))
        IF(KINIT.EQ.3) READ(KF28) GQC(:,:,1,:) ! dummy for deposition
        READ(KF28) GQC                         ! concentration array

!       loop through the current grid and find position on input grid
        DO K=1,nz

!          find the corresponding height on the input grid (ggg)
           kk=nzg
           DO WHILE (kk.GT.1.AND.0.5*(ggg(kk)+ggg(kk-1)).LT.sig(k))
              kk=kk-1
           END DO
           kk=MAX(1,kk)
       
        DO J=1,ny
        DO I=1,nx

!          location on the current grid
           xlat=(j-1)*dlat+clat1-dlat/2.0
           xlon=(i-1)*dlon+clon1-dlon/2.0

!          index on the old grid
           jj=1+(xlat-clatg+glat/2.0)/glat
           ii=1+(xlon-clong+glon/2.0)/glon
           jj=MAX(1,MIN(jj,nyg))
           IF(ii.LT.  1) ii=nxg+ii
           IF(ii.GT.nxg) ii=ii-nxg

!          fill in values to compuational array
           XXX(i,j,k,:)=GQC(ii,jj,kk,:)

        END DO
        END DO
        END DO
        WRITE(KF21,*)' NOTICE gemini: blended initialization from GBL3DM.BIN'
        DEALLOCATE (gqc)

     ELSE
!       input grid is the same as the current simulation grid
        IF(KINIT.EQ.3) READ(KF28) DDD
        READ(KF28) XXX
        WRITE(KF21,*)' NOTICE gemini: direct initialization from GBL3DM.BIN'
     END IF

!    independent global model version normalized (gblmdl37)
     IF(KINIT.EQ.4) XXX=XXX/CFACT/ROW  

     REWIND(KF28)
     DEALLOCATE(ggg)

!---------------------------------------------
! single value initialization

  ELSE
     CINI = CKON/CFACT/ROW
     DO j=1,ny	 
        DO i=1,nx
           DO k=1,nz
              PRES=SIG(K)*(SFC(i,j)-PTOP)+PTOP
              IF(pres.GT.150.0)THEN
                 XXX(i,j,k,mp)=CINI
              ELSEIF(pres.GT.70.0)THEN
                 XXX(i,j,k,mp)=CBOT*CINI
              ELSE
                 XXX(i,j,k,mp)=CTOP*CINI
              END IF
           END DO
        END DO
     END DO
     WRITE(KF21,*)' NOTICE gemini: constant value initialization at - ',CKON 

  END IF

END SUBROUTINE gemini

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

SUBROUTINE linreg (xvar,yvar,slope,yintc)

  USE gemkon   

  IMPLICIT NONE

  REAL*4, INTENT(IN)  :: xvar(:), yvar(:)
  REAL*4, INTENT(OUT) :: slope, yintc

  INTEGER*4 :: numb
  REAL*4    :: sumx,sumy,sumxy,sumx2

! summations
  NUMB=SIZE(xvar,1)
  SUMX=SUM(xvar)
  SUMY=SUM(yvar)
  SUMXY=SUM(xvar*yvar)
  SUMX2=SUM(xvar*xvar)

! final coefficients
  IF(numb.GE.2)THEN
     SLOPE=(NUMB*SUMXY-SUMX*SUMY)/(NUMB*SUMX2-SUMX*SUMX)
     YINTC=(SUMY-SLOPE*SUMX)/NUMB
  ELSEIF(numb.EQ.1)THEN
     SLOPE=0.0
     YINTC=SUMY
  ELSE
     SLOPE=0.0
     YINTC=0.0
  END IF

END SUBROUTINE linreg
