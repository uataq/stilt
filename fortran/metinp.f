!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  METINP           METeorological INPut reads meteo data file
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   METEOROLOGICAL INPUT READS METEOROLOGICAL DATA FILE STARTING AT
!   AND FILLS DATA ARRAYS FOR ONLY ONE TIME PERIOD.  MISSING DATA ARE
!   SKIPPED AND OLD VALUES REMAIN IN THE ARRAY.  HOWEVER THIS MAY CAUSE
!   PROBLEMS IN OTHER ROUTINES WHEN DATA ARE REMAPPED TO NEW COORDINATES
!   OR UNIT CONVERSIONS. NO OTHER CHECKS ARE PERFORMED AND DATA FILL
!   ARRAY ACCORDING TO THE LEVEL NUMBER CONTAINED IN THE RECORD LABEL.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 10 Mar 1998 (RRD)
!                  16 Dec 1998 (RRD) - added end-of-file test KEND
!                  19 Apr 1999 (RRD) - added terrain height
!                  17 Jun 1999 (RRD) - correction to handle minutes
!                  02 Nov 1999 (RRD) - flag to test for data with sfc=0
!                  28 Jul 2000 (RRD) - generalized input of flux fields
!                  28 Sep 2000 (RRD) - fortran90 upgrade
!                  21 Jun 2001 (RRD) - variable change from t to a
!                  27 Sep 2001 (RRD) - simultaneous multiple meteorology
!                  18 Oct 2001 (RRD) - extended grid domains
!                  05 Dec 2001 (RRD) - additional sfc height name
!                  26 Feb 2002 (RRD) - downward shortwave flux 
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  23 Dec 2002 (RRD) - no missing time periods permitted
!                  16 Sep 2003 (RRD) - renamed sensible heat flux variable
!                  17 Oct 2003 (RRD) - added turbulent kinetic energy
!                  10 Nov 2003 (RRD) - check for velocity variances
!                  02 Apr 2004 (RRD) - generic file unit numbers
!                  03 Dec 2004 (BS)  - precipitation rate ('PRT')
!                  14 Apr 2006 (RRD) - handle various versions of WRF
!                  22 May 2006 (RRD) - mixed layer depth
!                  05 Mar 2008 (RRD) - initialize TKE at all levels
!                  28 Apr 2008 (RRD) - missing fields permitted 10->20  
!
! USAGE:  CALL METINP(BACK,KG,KT,KUNIT,KREC,LX1,LY1,NXS,NYS,NZS,MC,KEND,
!              IFHR,ZT,DS,P0,T0,U0,V0,UF,VF,HF,RT,ZI,U,V,W,A,Q,P,E,H,X)
!  
!   INPUT ARGUMENT LIST:       see below
!   OUTPUT ARGUMENT LIST:      see below
!   INPUT FILES:               units 10,11,12,etc - to # input files 
!   OUTPUT FILES:              unit KF21 for diagnostic MESSAGE file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

!dwen(20090730) SUBROUTINE METINP(BACK,KG,KT,KUNIT,KREC,LX1,LY1,NXS,NYS,NZS,MC,KEND,     &
!dwen(20090730)                  IFHR,ZT,DS,P0,T0,U0,V0,UF,VF,HF,RT,ZI,U,V,W,A,Q,P,E,H,X)

! CHG:(11/20/01) add RC conv precip, TC total cld, SW radiation
! CHG:(11/20/01) add LF latent heat flux
! CHG:(12/04/01) add LC low cloud cover
! CHG:(01/22/03) add SM soil moisture
! JCL:(03/27/03) add arrays that get assigned if reading RAMS fields
! CHG(09/23/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
! CHG(09/25/03) add RAMS convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
! CHG(09/25/03) add RAMS turb. kin. energy TKEN

subroutine metinp(BACK,KG,KT,KUNIT,KREC,LX1,LY1,NXS,NYS,NZS,MC,KEND,     &
     IFHR,ZT,DS,P0,T0,U0,V0,UF,VF,HF,RT,ZI,U,V,W,A,Q,P,E,H,X,&
     W0,fluxflg, deepflg, shallflg,              &
     muu,muv,mu,msfu,msfv,msft,RC,LF,TC,LC,SW,SM,  &
     TLRAMS,SIGWRAMS,cfu1,CFU2,CFD1,DFU1,       &
     DFU2,EFU1,EFU2,DFD1,EFD1,TKE)         


  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE


  !-------------------------------------------------------------------------------
  ! argument list variables
  !-------------------------------------------------------------------------------

  LOGICAL,      INTENT(IN)    :: back      ! integration direction
  INTEGER,      INTENT(IN)    :: kg        ! number of active grid
  INTEGER,      INTENT(IN)    :: kt        ! number of active time
  INTEGER,      INTENT(IN)    :: kunit     ! input device unit number
  INTEGER,      INTENT(IN)    :: krec      ! record number of index record
  INTEGER,      INTENT(IN)    :: lx1,ly1   ! lower left  of subgrid FG unit
  INTEGER,      INTENT(IN)    :: nxs,nys   ! dimensions of sub-grid
  INTEGER,      INTENT(IN)    :: nzs       ! number of data levels to read
  INTEGER,      INTENT(INOUT) :: mc        ! accumulated minutes of data read
  INTEGER,      INTENT(IN)    :: kend      ! last valid record number of file

  INTEGER,      INTENT(OUT)   :: ifhr      ! current forecast hour
  REAL,         INTENT(OUT)   :: zt (:,:)  ! terrain height
  REAL,         INTENT(OUT)   :: ds (:,:)  ! downward shortwave
  REAL,         INTENT(OUT)   :: p0 (:,:)  ! surface pressure 
  REAL,         INTENT(OUT)   :: u0 (:,:)  ! low level u wind
  REAL,         INTENT(OUT)   :: v0 (:,:)  ! low level v wind
  REAL,         INTENT(OUT)   :: t0 (:,:)  ! low level temperature
  REAL,         INTENT(OUT)   :: uf (:,:)  ! u momentum flux
  REAL,         INTENT(OUT)   :: vf (:,:)  ! v momentum flux
  REAL,         INTENT(OUT)   :: hf (:,:)  ! sensible heat flux
  REAL,         INTENT(OUT)   :: rt (:,:)  ! rainfall total
  REAL,         INTENT(OUT)   :: zi (:,:)  ! mixed layer depth
  REAL,         INTENT(OUT)   :: u (:,:,:) ! upper level u wind
  REAL,         INTENT(OUT)   :: v (:,:,:) ! upper level v wind
  REAL,         INTENT(OUT)   :: w (:,:,:) ! upper level w wind
  REAL,         INTENT(OUT)   :: a (:,:,:) ! upper level temperature
  REAL,         INTENT(OUT)   :: q (:,:,:) ! upper level moisture
  REAL,         INTENT(OUT)   :: p (:,:,:) ! upper level pressure
  REAL,         INTENT(OUT)   :: e (:,:,:) ! turbulent kinetic energy or ...
  ! velocity variance v'2
  REAL,         INTENT(OUT)   :: h (:,:,:) ! velocity variance u'2 
  REAL,         INTENT(OUT)   :: x (:,:,:) ! velocity variance w'2

  !dwen(20090730)  ***************************
  real,         intent(out)   :: w0(:,:)
  real,         intent(out)   :: muu(:,:)
  real,         intent(out)   :: muv(:,:)
  real,         intent(out)   :: mu(:,:)
  Logical,      intent(inout) :: fluxflg,deepflg,shallflg
  real,         intent(out)   :: msfu(:,:)
  real,         intent(out)   :: msfv(:,:)
  real,         intent(out)   :: msft(:,:)
  real,         intent(out)   :: rc(:,:)
  real,         intent(out)   :: lf(:,:)
  real,         intent(out)   :: tc(:,:)
  real,         intent(out)   :: lc(:,:)
  real,         intent(out)   :: sw(:,:)
  real,         intent(out)   :: sm(:,:)
  real,         intent(out)   :: tlrams(:,:,:)
  real,         intent(out)   :: sigwrams(:,:,:)
  real,         intent(out)   :: cfu1(:,:,:)
  real,         intent(out)   :: cfu2(:,:,:)
  real,         intent(out)   :: cfd1(:,:,:)
  real,         intent(out)   :: dfu1(:,:,:)
  real,         intent(out)   :: dfu2(:,:,:)
  real,         intent(out)   :: dfd1(:,:,:)
  real,         intent(out)   :: efu1(:,:,:)
  real,         intent(out)   :: efu2(:,:,:)
  real,         intent(out)   :: efd1(:,:,:)
  real,         intent(out)   :: tke(:,:,:)
  !dwen ***************************************
  !-------------------------------------------------------------------------------
  ! internal variables
  !-------------------------------------------------------------------------------

  LOGICAL                     :: sfc1
  INTEGER                     :: kret,misst,missv    
  CHARACTER(50)               :: label 
  CHARACTER(4)                :: varb  
  CHARACTER(108)              :: header
  INTEGER                     :: nxg,nyg,jrec,nndx,mnts,nvar,lenh
  INTEGER                     :: kk,nn,iy,im,id,ih,ic,ll,nexp
  REAL                        :: prec,var1
  CHARACTER(1), ALLOCATABLE   :: cdata(:)  ! packed data array of length NX*NY
  INTEGER                     :: ksum = -1 ! skip checksum calculation

  !dwen(20090730)*****************************************
  character(4)                :: fluxvars(2) = (/'MUU0','MUV0'/), &
       wrfvars(6) = (/'MUBA','MUPE','ALT0','MSFU','MSFV','MSFT'/), &
       deepvars(6) = (/'CFU1','CFD1','DFU1','DFD1','EFU1','EFD1'/), &
       shallvars(3) = (/'CFU2','DFU2','EFU2'/)
  logical                     :: haveflux(2),havewrf(6),havedeep(6),haveshall(3)
  logical                     :: ramsflg, ecmflg, sxtnbit, awrfflg, ppertflg

  real, allocatable           :: ppert(:,:,:), mupert(:,:)
  real                        :: pmin,pminp,ibad_pmin
  integer                     :: ig
  !dwen(20090730)*****************************************
  !-------------------------------------------------------------------------------
  ! external variables
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  INTERFACE
     !-------------------------------------------------------------------------------

     !dwen(20090819)  SUBROUTINE PAKINP(RVAR,CVAR,NX,NY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)
     SUBROUTINE PAKINP(sxtnbit,RVAR,CVAR,NX,NY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)
       IMPLICIT NONE
       REAL,          INTENT(OUT)   :: rvar (:,:)     ! real data unpacked
       CHARACTER(1),  INTENT(IN)    :: cvar (:)       ! packed input of NX*NY
       INTEGER,       INTENT(IN)    :: nx,ny          ! size of input array  
       INTEGER,       INTENT(IN)    :: nx1,ny1        ! optional sub-grid left edge 
       INTEGER,       INTENT(IN)    :: lx,ly          ! length of sub-grid
       REAL,          INTENT(IN)    :: prec           ! precision of packed data 
       INTEGER,       INTENT(IN)    :: nexp           ! packing scaling exponent
       REAL,          INTENT(IN)    :: var1           ! value of array at (1,1)
       INTEGER,       INTENT(INOUT) :: ksum           ! rotating checksum 
       !dwen(20090731) *****************
       logical,        intent(in)    :: sxtnbit        ! .TRUE.  for 16 bit per number in ARL
       !********************************
     END SUBROUTINE pakinp
     !-------------------------------------------------------------------------------
     SUBROUTINE TM2MIN(IY,IM,ID,IH,MN,MACC)
       IMPLICIT NONE
       INTEGER,  INTENT(IN)   :: iy,im,id,ih,mn       ! date and time
       INTEGER,  INTENT(OUT)  :: macc                 ! minutes since 1 Jan 1970
     END SUBROUTINE tm2min

     !-------------------------------------------------------------------------------
  END INTERFACE
  !-------------------------------------------------------------------------------

  ! check for end-of-file
  IF(KREC.LT.1)THEN        
     WRITE(*,*)'WARNING metinp: trying to read before file start - ',KREC
     RETURN
  END IF

  IF(KREC.GT.KEND)THEN     
     WRITE(*,*)'WARNING metinp: trying to read past end-of-file - ',KREC
     RETURN
  END IF

  ! set flag default to indicate that surface fields have the level
  ! index set to 0 (new style) rather than 1 (old style)
  SFC1=.FALSE.

  ppertflg = .FALSE. 

  ! see defgrid.inc for following definitions
  NXG=GRID(KG,KT)%NX
  NYG=GRID(KG,KT)%NY

  ! set starting record number at index record position
  JREC=KREC


  !**********************************
  ! CHG(10/02/03) add ramsflag, need for 16 bit...
  RAMSFLG = GRID(KG,kt)%MODEL_ID == 'RAMS'
  ECMFLG  = GRID(KG,kt)%MODEL_ID(1:2) == 'EC'
  AWRFFLG = GRID(KG,kt)%MODEL_ID(2:4) .EQ. 'WRF'
  SXTNBIT = RAMSFLG .OR. GRID(KG,kt)%MODEL_ID(1:3) == 'ECX' .OR. GRID(KG,kt)%MODEL_ID == 'DWRF'
  haveflux = .FALSE.
  havewrf = .FALSE.
  havedeep = .FALSE.
  haveshall = .FALSE.
  DREC(KG,kt)%TKEF = .FALSE.
  !************************************

  ! packed data buffer
  !dwen(20090826)  ALLOCATE (cdata(nxg*nyg),STAT=kret)
  if (SXTNBIT) then
     allocate(cdata(2*nxg*nyg), ppert(size(p,1),size(p,2),size(p,3)), &
          & mupert(size(p,1),size(p,2)),STAT=kret)
  else
     allocate(cdata(nxg*nyg), ppert(size(p,1),size(p,2),size(p,3)), &
          & mupert(size(p,1),size(p,2)),STAT=kret)
  endif
  IF(kret.NE.0)THEN
     WRITE(*,*)'*ERROR* metinp: memory allocation' 
     STOP 900
  END IF

  MISST=0

  ! loop to this position when missing data found
  inloop: DO 
     MISSV=0

     !-------------------------------------------------------------------------------
     ! read index record (or first sfc field) for data time
     !-------------------------------------------------------------------------------

     READ(KUNIT,REC=JREC,ERR=910)LABEL,HEADER
     READ(LABEL,'(7I2,A4,I4,2E14.7)')IY,IM,ID,IH,IFHR
     CALL TM2MIN(IY,IM,ID,IH,0,MC)

     IF((DREC(KG,KT)%TYPE).EQ.2)THEN
        READ(HEADER(1:108),'(4X,I3,97X,I4)')IFHR,LENH
        !    the number of index records per time period
        !dwen(20090730)****************************
        if (sxtnbit) then                      !copied from STILT
           nndx = lenh/(nxg*nyg*2)+1           !copied from STILT
        else                                   !copied from STILT
           !dwen(20090730)********************************
           NNDX=LENH/(NXG*NYG)+1
        endif                                  !copied from STILT
        !    new format data skip index record for subsequent reads
        JREC=JREC+NNDX
        !    decode extended portion for forecast hour, and minutes
        READ(HEADER,'(4X,I3,I2)')IFHR,MNTS
        !    add minutes field to accumulated time (6/17/99)
        MC=MC+MNTS
     END IF

     ! adjust read for data offset (old style SH)
     JREC=JREC+DREC(KG,KT)%OFFSET

     !-------------------------------------------------------------------------------
     ! loop through only number of levels in sub-grid
     !-------------------------------------------------------------------------------
     !dwen(20090730)****************************
     ppert = 0  !in case ppert is available at only some levels
     w0 = 0  !use zero w at sfc unless available from input file
     !dwen(20090730)****************************

     zloop : DO KK=0,NZS

        ! start with level 0 surface parameters
        NVAR=DREC(KG,KT)%NUM_VARB(KK+1)

        varloop : DO NN=1,NVAR

           READ(KUNIT,REC=JREC,ERR=920)LABEL,CDATA
           READ(LABEL,'(6I2,2X,A4,I4,2E14.7)') IY,IM,ID,IH,IC,LL,VARB,NEXP,PREC,VAR1

           !-------------------------------------------------------------------------------
           ! test for missing data
           !-------------------------------------------------------------------------------

           IF(IC.LT.0.OR.VARB.EQ.'NULL')THEN
              !    eta missing flux correction
              IF((.NOT.DREC(KG,KT)%UFLX).AND.(VARB.EQ.'UMOF'.OR.VARB.EQ.'VMOF'))THEN
                 JREC=JREC+1
                 CYCLE varloop
              ELSE
                 MISSV=MISSV+1
                 WRITE(KF21,*)'WARNING metinp: missing data - ',VARB
                 WRITE(KF21,*)'          Time: ',IY,IM,ID,IH,IC
                 WRITE(KF21,*)'         Level: ',LL,'   Variable: ',NN
                 JREC=JREC+1
                 CYCLE varloop
              END IF
           END IF

           ! All old format data files will always have the surface data
           ! at the beginning of each time period.  Prior to 1992 the surface
           ! index was set to #1 and after that the surface started at #0.

           IF(KK.EQ.0.AND.NN.EQ.1.AND.LL.EQ.1)SFC1=.TRUE.

           ! old style data that starts at level #1 adjust index by one
           IF(SFC1)LL=LL-1

           !------------------------------------------------------------------------------
           ! initialize certain variables
           !------------------------------------------------------------------------------

           IF(DREC(KG,KT)%TKEN) E=0.20

           !------------------------------------------------------------------------------
           ! upper levels
           !------------------------------------------------------------------------------


           IF(LL.GT.0)THEN
              SELECT CASE (VARB)

                 !dwen(20090730)  added SXTNBIT to the argument list of PAKINP

              CASE ('UWND') 
                 !dwen     CALL PAKINP(U(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 CALL PAKINP(SXTNBIT,U(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !IF(DREC(KG,KT)%WRF.EQ.3) U(:,:,LL)=U(:,:,LL)/P(:,:,LL)/100.0
              CASE ('VWND') 
                 CALL PAKINP(SXTNBIT,V(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !IF(DREC(KG,KT)%WRF.EQ.3) V(:,:,LL)=V(:,:,LL)/P(:,:,LL)/100.0
              CASE ('WWND')
                 CALL PAKINP(SXTNBIT,W(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !IF(DREC(KG,KT)%WRF.EQ.3) W(:,:,LL)=W(:,:,LL)/100.0
              CASE ('DZDT')
                 CALL PAKINP(SXTNBIT,W(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('TEMP')
                 CALL PAKINP(SXTNBIT,A(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('THET')
                 CALL PAKINP(SXTNBIT,A(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('SPHU','RELH') 
                 CALL PAKINP(SXTNBIT,Q(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('HGTS') 
                 !    pressure level input data requires heights of pressure surface
                 !    height will be replaced by pressure in prfprs subroutine
                 CALL PAKINP(SXTNBIT,P(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('PRES') 
                 !    wrf=3 base state pressure or wrf=1 sum base and perturbation
                 CALL PAKINP(SXTNBIT,P(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('PPRE') 
                 !    wrf=3 perturbation pressure (assume it comes after PRES and before U,V)
                 !     CALL PAKINP(SXTNBIT,E(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !    units appear to be hPa rather than  Pa as indicated in documentation
                 !     P(:,:,LL)=P(:,:,LL)+E(:,:,LL)
                 !     A(:,:,LL)=A(:,:,LL)*(P(:,:,LL)/1000.0)**0.286
                 !dwen ************************
                 ppertflg = .TRUE. 
                 call pakinp(SXTNBIT,ppert(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 ! JCL(03/27/03):new RAMS variables
                 ! JCL(03/27/03):turbulence variables are directly calculated by RAMS, based on TKE
              case ('TLGR')
                 if (GRID(KG,kt)%MODEL_ID == 'RAMS')    &
                      call pakinp(SXTNBIT,tlrams(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('SIGW')
                 if (GRID(KG,kt)%MODEL_ID == 'RAMS')    &
                      call pakinp(SXTNBIT,sigwrams(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 ! convective fluxes CFU1 CFU2 CFD1 DFU1 DFU2 EFU1 EFU2 DFD1 EFD1
              case ('CFU1')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        & 
                      call pakinp(SXTNBIT,cfu1(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('CFU2')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,cfu2(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('CFD1')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,cfd1(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('DFU1')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,dfu1(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('DFU2')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,dfu2(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('EFU1')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,efu1(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('EFU2')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,efu2(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('DFD1')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,dfd1(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              case ('EFD1')
                 if ( RAMSFLG .OR. ECMFLG .OR. awrfflg )        &
                      call pakinp(SXTNBIT,efd1(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 ! CHG(09/25/03) RAMS turb. kin. energy TKEN
              case ('TKEN')
                 CALL PAKINP(SXTNBIT,E(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 if ( RAMSFLG .OR. awrfflg ) then 
                    DREC(KG,kt)%TKEF = .TRUE.
                    call pakinp(SXTNBIT,tke(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 end if

              case ('ALT0')
                 IF (awrfflg) THEN
                    call pakinp(SXTNBIT,tlrams(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 endif

                 !dwen ************************

              CASE ('UVAR') 
                 CALL PAKINP(SXTNBIT,H(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('VVAR') 
                 CALL PAKINP(SXTNBIT,E(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              CASE ('WVAR') 
                 CALL PAKINP(SXTNBIT,X(:,:,LL),CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)

              END SELECT
              IF (awrfflg) THEN
                 ! check for additional 3d variables needed for momentum and convective fluxes
                 where (VARB .EQ. wrfvars) havewrf = .TRUE.
                 where (VARB .EQ. deepvars) havedeep = .TRUE.
                 where (VARB .EQ. shallvars) haveshall = .TRUE.
              endif

              !-------------------------------------------------------------------------------
              ! surface level fields
              !-------------------------------------------------------------------------------

           ELSE
              !dwen(20090730)  SELECT CASE (VARB(1:3))
              !dwen(20090730)  used if statement instead of case statement
              !dwen(20090731)    so can VARB(1:3) and VARB(1:4) to discriminate 
              !dwen(20090730)      MSFU, MSFV adn MSFT

              !dwen  CASE ('PRSS')
              !dwen(20090730)  add SXTNBIT to argument list of PAKINP
              !dwen     CALL PAKINP(P0,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              if (varb == 'PRSS')   & 
                   CALL PAKINP(SXTNBIT,P0,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen  CASE ('TMPS','T02M') 
              if (varb == 'TMPS' .OR. varb == 'T02M')        &
                   CALL PAKINP(SXTNBIT,T0,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen  CASE ('U10M')
              if (varb == 'U10M')                            &
                   CALL PAKINP(SXTNBIT,U0,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen  CASE ('V10M')
              if (varb == 'V10M')                            &
                   CALL PAKINP(SXTNBIT,V0,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen(20090730)  CASE ('SHG','HGT')
              !dwen(20090730)  IF (VARB == 'SHGT' .OR. VARB == 'TERR')  in STILT
              !dwen  CASE ('SHG','HGT')
              if (varb == 'SHGT' .OR. varb == 'TERR')        &
                   CALL PAKINP(SXTNBIT,ZT,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen  CASE ('MXH','HPB','PBL')
              if (varb(1:3)=='MXH' .OR. varb(1:3)=='HPB' .OR. varb(1:3)=='PBL') then
                 CALL PAKINP(SXTNBIT,ZI,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              end if
              !dwen  case ('WWND','DZDT')
              if (varb == 'WWND' .OR. varb == 'DZDT')        &
                   call pakinp(sxtnbit,w0,cdata,nxg,nyg,lx1,ly1,nxs,nys,prec,nexp,var1,ksum)
              !*************************
              !dwen  CASE ('TPP')
              if (varb(1:3) == 'TPP')                        &
                   CALL PAKINP(SXTNBIT,RT,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen  CASE ('PRT')
              if (varb(1:3) == 'PRT')                        &
                   CALL PAKINP(SXTNBIT,RT,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              IF(DREC(KG,KT)%USTR)  then       
                 !dwen  CASE ('UST')
                 if (varb(1:3) == 'UST')                        &
                      CALL PAKINP(SXTNBIT,UF,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              else
                 !dwen  CASE ('UMO','EXC')
                 if (varb(1:3) == 'UMO' .OR. varb(1:3) == 'EXC')    &
                      CALL PAKINP(SXTNBIT,UF,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !dwen  CASE ('VMO')
                 if (varb(1:3) == 'VMO')                            &
                      CALL PAKINP(SXTNBIT,VF,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              endif
              IF(DREC(KG,KT)%TSTR) then         
                 !dwen  CASE ('TST')
                 if (varb(1:3) == 'TST')                            &
                      CALL PAKINP(SXTNBIT,HF,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              else
                 !dwen  CASE ('SHT','HFL')
                 if (varb(1:3) == 'SHT' .OR. varb(1:3) == 'HFL') then
                    CALL PAKINP(SXTNBIT,HF,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                    !  correction for sign reversal in ECMWF/ETA/EDAS sensible heat flux
                    IF (GRID(KG,kt)%MODEL_ID(1:2) == 'EC'  .OR.                 &
                         GRID(KG,kt)%MODEL_ID      == ' ETA' .OR. GRID(KG,kt)%MODEL_ID == 'EDAS') HF = -HF
                 ENDIF
              endif

              !dwen  CASE ('LHT','WFL')
              if (varb(1:3) == 'LHT' .OR. varb(1:3) == 'WFL' ) then
                 CALL PAKINP(SXTNBIT,LF,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !  correction for sign reversal in ECMWF/ETA/EDAS sensible heat flux
                 IF (GRID(KG,kt)%MODEL_ID(1:2) == 'EC'  .OR.                 &
                      GRID(KG,kt)%MODEL_ID      == ' ETA' .OR. GRID(KG,kt)%MODEL_ID == 'EDAS') LF = -LF
              ENDIF

              !dwen   case ('CPP','CPR')
              if (varb(1:3) == 'CPP' .OR. varb(1:3) == 'CPR')                &
                   CALL PAKINP(SXTNBIT,RC,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen   case ('TCL')
              if (varb(1:3) == 'TCL')                                        &
                   CALL PAKINP(SXTNBIT,TC,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen   case ('LCL')
              if (varb(1:3) == 'LCL')                                         &
                   CALL PAKINP(SXTNBIT,LC,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              !dwen   case ('DSW')
              if (varb(1:3) == 'DSW') then
                 CALL PAKINP(SXTNBIT,SW,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 sw = max(sw,0.0)
              endif
              !dwen   case ('SOL')
              if (varb(1:3) == 'SOL')                                         &
                   CALL PAKINP(SXTNBIT,SM,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              IF (awrfflg) THEN
                 ! check for additional 2d variables needed for momentum fluxes
                 where (VARB .EQ. fluxvars) haveflux = .TRUE.
                 where (VARB .EQ. wrfvars) havewrf = .TRUE.  !wrfvars has both 2d and 3d
                 !dwen     case ('MUU0')
                 if ( varb == 'MUU0')                                           &
                      CALL PAKINP(SXTNBIT,muu,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !dwen     case ('MUV0')
                 if ( varb == 'MUV0')                                           &
                      CALL PAKINP(SXTNBIT,muv,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !dwen     case ('MUBA')
                 if ( varb == 'MUBA')                                           &
                      CALL PAKINP(SXTNBIT,mu,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !dwen     case ('MUPE')
                 if ( varb == 'MUPE')                                           &
                      CALL PAKINP(SXTNBIT,mupert,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !dwen     case ('MSFU')
                 if ( varb == 'MSFU')                                           &
                      CALL PAKINP(SXTNBIT,msfu,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !dwen     case ('MSFV')
                 if ( varb == 'MSFV')                                           &
                      CALL PAKINP(SXTNBIT,msfv,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
                 !dwen     case ('MSFT')
                 if ( varb == 'MSFT')                                           &
                      CALL PAKINP(SXTNBIT,msft,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              endif

              !dwen  CASE ('DSW')
              if (varb(1:3) == 'DSW') then
                 IF(DREC(KG,KT)%DSWF)         &
                      CALL PAKINP(SXTNBIT,DS,CDATA,NXG,NYG,LX1,LY1,NXS,NYS,PREC,NEXP,VAR1,KSUM)
              endif

              !dwen  END SELECT
           END IF

           JREC=JREC+1

        END DO varloop
     END DO zloop

     !dwen(20090731) ***************************
     !dwen            copied from STILT
     IF (ppertflg) THEN
        p = p + ppert
        pminp = 1e10
        DO ll=1,nzs
           ibad_pmin = 0
           pmin = minval(p(:,:,ll))
           IF (pmin < 0 .OR. pmin > pminp) THEN


              ibad_pmin = ll
              WRITE (KFJCLMSG,*) 'metinp: bad pressure at time,ll,pminp,pmin : ',IY,IM,ID,IH,IC,ll,pminp,pmin
              WRITE (*,*) 'metinp: bad pressure at time,ll,pminp,pmin : ',IY,IM,ID,IH,IC,ll,pminp,pmin
           ENDIF
           pminp = pmin
        ENDDO
        IF (ibad_pmin .NE. 0) THEN
           DO ll=1,nzs
              write (KFJCLMSG,'(a,i4,a,2g15.6)') 'll=',ll,' Min,max ptot=',minval(p(:,:,ll)),maxval(p(:,:,ll))
              IF (ppertflg) write (KFJCLMSG,'(a,i4,a,2g15.6) ') &
                   'll=',ll,' Min,max ppert=',minval(ppert(:,:,ll)),maxval(ppert(:,:,ll))
           ENDDO
           stop 'metinp: bad pressure'
        ENDIF
     ENDIF
     IF (awrfflg) THEN
        fluxflg = any(haveflux)
        IF (fluxflg .AND. .NOT. all(haveflux)) THEN
           WRITE (KFJCLMSG,*) 'metinp: missing one or more fields for mom flux input:', &
                ' haveflux= ',haveflux,' fluxvars= ',fluxvars
           WRITE (*,*) 'metinp: missing one or more fields for mom flux input:', &
                ' haveflux= ',haveflux,' fluxvars= ',fluxvars
           stop 'metinp: missing one or more fields for mom flux input'
        ENDIF
        IF (.NOT. all(havewrf)) THEN
           WRITE (KFJCLMSG,*) 'metinp: missing one or more fields for wrf input:', &
                ' havewrf= ',havewrf,' wrfvars= ',wrfvars
           WRITE (*,*) 'metinp: missing one or more fields for wrf input:', &
                ' havewrf= ',havewrf,' wrfvars= ',wrfvars
           stop 'metinp: missing one or more fields for wrf input'
        ENDIF
        deepflg = all(havedeep)
        shallflg = all(haveshall)
        IF (.NOT. DREC(KG,kt)%TKEF) THEN
           tke(:,:,:) = 0.
        END IF
        IF (.NOT. deepflg) THEN
           cfu1(:,:,:) = 0.
           cfd1(:,:,:) = 0.
           dfu1(:,:,:) = 0.
           dfd1(:,:,:) = 0.
           efu1(:,:,:) = 0.
           efd1(:,:,:) = 0.
        ENDIF
        IF (.NOT. shallflg) THEN
           cfu2(:,:,:) = 0.
           dfu2(:,:,:) = 0.
           efu2(:,:,:) = 0.
        ENDIF
        mu = mu + mupert
     ENDIF !awrfflg
     !dwen ******************************

     !-------------------------------------------------------------------------------
     ! up to five missing variables per time period permitted
     !-------------------------------------------------------------------------------

     IF(MISSV.GT.20)THEN
        MISST=MISST+1
        !    only one missing time period permitted
        !    IF(MISST.LE.1)THEN  (rrd - 12/23/2002)
        IF(MISST.LE.0)THEN
           !       skip back to previous time period, otherwise continue
           IF(BACK)JREC=JREC-2*DREC(KG,KT)%REC_PER
           CYCLE inloop
        ELSE
           WRITE(KF21,*)'*ERROR* metinp: too many missing times - ',MISST
           WRITE(KF21,*)' Current version no missing time permitted '         
           STOP 900
        END IF

     ELSE
        EXIT inloop
     END IF

  END DO inloop

  ! optional diagnostic information
  WRITE(KF21,*)' NOTICE metinp: ',GRID(KG,KT)%MODEL_ID,KG,KREC,NXS,NYS,  &
       LX1,LY1,MC,IY,IM,ID,IH

  DEALLOCATE (cdata,ppert,mupert)
  RETURN

  !-------------------------------------------------------------------------------
  ! error termination messages
  !-------------------------------------------------------------------------------

910 WRITE(*,*)'*ERROR* metinp: reading index record in meteo file'
  write (*,*) 'KG,KT,KUNIT,KREC,jrec=',KG,KT,KUNIT,KREC,jrec
  STOP 900
920 WRITE(*,*)'*ERROR* metinp: reading data record in meteo file'
  write (*,*) 'KG,KT,KUNIT,KREC,kk,nn,JREC=',KG,KT,KUNIT,KREC,kk,nn,jrec
  STOP 900

END SUBROUTINE metinp
