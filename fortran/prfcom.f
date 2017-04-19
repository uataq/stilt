!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PRFCOM           PRoFile COMmon driver for input meteo
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PROFILE COMMON IS THE COMMON DRIVER FOR THE FOUR METEOROLGICAL
!   PROFILE ANALYSIS PROGRAMS.  ROUTINE EXAMINES THE PROFILE AT EACH
!   SUB-GRID NODE AND CONVERTS THE INPUT DATA TO COMMON UNITS AND
!   INTERPOLATES THE SOUNDING TO THE INTERNAL MODEL SIGMA SYSTEM.
!   THE DIFFERENT INPUT DATA COORDINATE SYSTEMS ARE SUPPORTED FOR
!   CONVERSION: ABSOLUTE PRESSURE (PRFPRS), SIGMA PRESSURE (PRFSIG),
!   TERRAIN FOLLOWING SIGMA (PRFTER), AND ECMWF HYBRID.  AFTER UNITS
!   CONVERSION THE STABILITY ANALYSIS ROUTINES ARE CALLED AT EACH POINT.
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 01 Apr 1998 (RRD)
!                  18 Aug 1998 (RRD) - added isotropic turbulence variable
!                  04 Mar 1999 (RRD) - test for proper RH limits
!                  19 Apr 1999 (RRD) - added terrain array
!                  01 May 2000 (RRD) - change to prfter argument list
!                  28 Jul 2000 (RRD) - eta/edas sensible heat flux sign
!                  14 Aug 2000 (RRD) - added zflg to prfsig call
!                  29 Sep 2000 (RRD) - fortran90 upgrade
!                  16 Mar 2001 (RRD) - optional global lat lon grid
!                  02 Jul 2001 (RRD) - saving ambient temperature
!                  29 Aug 2001 (RRD) - simultaneous multiple meteorology
!                  17 Oct 2001 (RRD) - added mixing depth variable
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  25 Oct 2002 (RRD) - support for NMM pressure offset sigma
!                  07 Feb 2003 (RRD) - additional diagnostic arguments stbanl
!                  16 Sep 2003 (RRD) - stability function independent variable
!                  17 Oct 2003 (RRD) - added turbulent kinetic energy
!                  10 Nov 2003 (RRD) - incorporated U',V',W' mixing
!                  02 Dec 2003 (RRD) - new variance equations for short-range
!                  21 May 2004 (RRD) - test for negative variances
!                  02 Jul 2004 (RRD) - added heat flux logical test
!                  11 May 2005 (RRD) - argument list change to stbsnd
!                  08 Mar 2006 (RRD) - static stability parameter (ss)
!                  14 Apr 2006 (RRD) - wrf vertical motion options
!                  21 Nov 2006 (RRD) - day night tke partition
!                  12 Jun 2008 (RRD) - additional mixing depth options
!                  15 Aug 2008 (RRD) - refined stability analysis tests
!
! USAGE:  CALL PRFCOM(TKERD,TKERN,KG,KT,KSFC,GX,GY,Z0,ZT,NXS,NYS,NZS,ZMDL,
!              ZSG,NLVL,VMIX,KMIXD,KMIX0,ZI,P0,T0,U0,V0,UF,VF,SF,SS,U,V,W,
!              A,T,Q,P,E,H,X)
!
!   INPUT ARGUMENT LIST:      see below
!   OUTPUT ARGUMENT LIST:     see below
!   INPUT FILES:              none
!   OUTPUT FILES:             none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

! JCL:add the array of Lagrangian timescale 'TL' &
!     stddev of vertical velocity 'SIGW' to argument list
! JCL:(4/28/00)add array of mixed-layer heights 'ZML' as input & output
! CHG:(12/04/01)add array of lim. of conv. heights 'ZLOC' as input & output
! JCL:(4/3/02)add mass violation array as argument
! JCL:(9/16/02) ZISCALE is scaling factor for mixed-layer ht--used to prescribe mixed-layer ht
! CHG:(9/17/02) add 'ICONVECT' as convection flag
! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)

!dwen(20090803)  SUBROUTINE PRFCOM(TKERD,TKERN,KG,KT,KSFC,GX,GY,Z0,ZT,NXS,NYS,NZS,ZMDL,ZSG,    &
!dwen(20090803)                  NLVL,VMIX,KMIXD,KMIX0,ZI,P0,T0,U0,V0,UF,VF,HF,SF,SS,U,V,W,  &
!dwen(20090803)                  A,T,Q,P,E,H,X)
!dwen(20090806)  add hm to output horizontal mixing coefficient

subroutine PRFCOM(TKERD,TKERN,KG,KT,KSFC,GX,GY,Z0,              &
     ZT,NXS,NYS,NZS,ZMDL,ZSG,NLVL,VMIX,KMIXD,KMIX0,         &
     ZI,P0,T0,U0,V0,UF,VF,HF,SF,SS,U,V,W,A,T,Q,P,E,H,X,     &
     iconvect,w0,tl,d,sigw,zloc,delmass,           &
     ziscale,alt0,mu,muu,muv,msfu,msfv,msft,              &
     fluxflg, deepflg, shallflg, cfxup1,cfxup2,cfxdn1,      &
     dfxup1,dfxup2,efxup1,efxup2,dfxdn1,efxdn1,tke,xm,hm)     


  USE funits
  use module_defgrid ! meteorological array 
  IMPLICIT NONE


  !-------------------------------------------------------------------------------
  ! argument list variables
  !-------------------------------------------------------------------------------

  REAL,    INTENT(IN)    :: tkerd      ! day turbulent kinetic eneregy ratio
  REAL,    INTENT(IN)    :: tkern      ! night turbulent kinetic eneregy ratio
  INTEGER, INTENT(IN)    :: kg         ! grid indicator
  INTEGER, INTENT(IN)    :: kt         ! time indicator
  INTEGER, INTENT(IN)    :: ksfc       ! index top of the sfc layer
  !dwen(20090823) 
  !  REAL,    INTENT(IN)    :: gx (:,:)   ! grid size (m)
  !  REAL,    INTENT(IN)    :: gy (:,:)   ! grid size (m)
  REAL,    INTENT(INout)    :: gx (:,:)   ! grid size (m)
  REAL,    INTENT(INout)    :: gy (:,:)   ! grid size (m)
  REAL,    INTENT(IN)    :: z0 (:,:)   ! aerodynamic roughness length (m)
  REAL,    INTENT(INOUT) :: zt (:,:)   ! terrain height elevations (m)
  INTEGER, INTENT(IN)    :: nxs,nys    ! subgrid dimensions
  INTEGER, INTENT(IN)    :: nzs        ! number of input data levels
  REAL,    INTENT(IN)    :: zmdl       ! internal model top (meters)
  REAL,    INTENT(IN)    :: zsg (:)    ! array internal sigma levels
  INTEGER, INTENT(IN)    :: nlvl       ! number of output levels
  LOGICAL, INTENT(IN)    :: vmix       ! indicator for mixing computation
  INTEGER, INTENT(IN)    :: kmixd      ! mixed layer depth options
  INTEGER, INTENT(IN)    :: kmix0      ! minimum mixing depth

  REAL,    INTENT(INOUT) :: zi (:,:)   ! mixing depth  
  REAL,    INTENT(INOUT) :: p0 (:,:)   ! pressue surface variable 
  REAL,    INTENT(IN)    :: t0 (:,:)   ! temperature 
  REAL,    INTENT(IN)    :: u0 (:,:)   ! u velocity 
  REAL,    INTENT(IN)    :: v0 (:,:)   ! v velocity 
  REAL,    INTENT(INOUT) :: uf (:,:)   ! u momentum flux and u*
  REAL,    INTENT(INOUT) :: vf (:,:)   ! v momentum flux and v*
  REAL,    INTENT(IN)    :: hf (:,:)   ! sensible heat flux
  REAL,    INTENT(OUT)   :: sf (:,:)   ! stability function
  REAL,    INTENT(OUT)   :: ss (:,:)   ! static stability  

  REAL,    INTENT(INOUT) :: u  (:,:,:) ! wind 3d variables       
  REAL,    INTENT(INOUT) :: v  (:,:,:) ! v component       
  REAL,    INTENT(INOUT) :: w  (:,:,:) ! w component        
  REAL,    INTENT(INOUT) :: a  (:,:,:) ! ambient temperature        
  REAL,    INTENT(OUT)   :: t  (:,:,:) ! potential temperature        
  REAL,    INTENT(INOUT) :: q  (:,:,:) ! moisture           
  REAL,    INTENT(INOUT) :: p  (:,:,:) ! pressure           
  REAL,    INTENT(INOUT) :: e (:,:,:)  ! turbulent kinetic energy or ...
  ! velocity variance v'2
  REAL,    INTENT(INOUT) :: h (:,:,:)  ! velocity variance u'2
  REAL,    INTENT(INOUT) :: x (:,:,:)  ! velocity variance w'2

  !dwen(20090803) ******************* 
  real,      intent(in)       :: w0(:,:)
  real,      intent(inout)    :: d(:,:,:)
  real,      intent(out)      :: tl(:,:,:)
  real,      intent(out)      :: sigw(:,:,:)
  real,      intent(out)      :: hm(:,:,:) !horizontal mixing coefficient, in order to discriminate from h

  ! CHG:(12/04/2001)array to store lim. of conv. heights
  real,      intent(out)         :: zloc(:,:)
  ! JCL:(4/3/02)mass violation grid [fraction of gridcell/min]
  real,      intent(out)   :: delmass(:,:,:)
  real,      intent(in)         :: ziscale
  real,      intent(in)         :: alt0(:,:,:)
  real,      intent(in)         :: mu(:,:)
  real,      intent(in)         :: muu(:,:)
  real,      intent(in)         :: muv(:,:)
  real,      intent(in)         :: msfu(:,:)
  real,      intent(in)         :: msfv(:,:)
  real,      intent(in)         :: msft(:,:)
  logical,   intent(inout) :: fluxflg  ! flag for WRF momentum flux input
  logical,   intent(inout) :: deepflg, shallflg
  real,      intent(inout)         :: cfxup1(:,:,:)
  real,      intent(inout)         :: cfxup2(:,:,:)
  real,      intent(inout)         :: cfxdn1(:,:,:)
  real,      intent(inout)         :: dfxup1(:,:,:)
  real,      intent(inout)         :: dfxup2(:,:,:)
  real,      intent(inout)         :: efxup1(:,:,:)
  real,      intent(inout)         :: efxup2(:,:,:)
  real,      intent(inout)         :: dfxdn1(:,:,:)
  real,      intent(inout)         :: efxdn1(:,:,:)
  real,      intent(inout)         :: tke(:,:,:)  !turbulence kinetic energy
  real,      intent(out)           :: xm(:,:,:)
  integer,   intent(in)            :: iconvect

  !************************************


  !-------------------------------------------------------------------------------
  ! internal variables
  !-------------------------------------------------------------------------------

  ! temporary holding arrays (note output arrays can be larger
  ! than input data due to finer interpolation

  REAL, ALLOCATABLE :: psg (:) ! input level data array
  REAL, ALLOCATABLE :: esg (:) ! working input level array

  REAL, ALLOCATABLE :: U1(:),V1(:),W1(:),T1(:),Q1(:),P1(:)
  REAL, ALLOCATABLE :: U2(:),V2(:),W2(:),T2(:),Q2(:),P2(:)

  REAL, ALLOCATABLE ::                   X1(:),H1(:),E1(:)
  REAL, ALLOCATABLE :: A2(:),Z2(:),D2(:),X2(:),H2(:),E2(:)

  !dwen(20090803) *******************************
  real, allocatable :: xv(:)
  real, allocatable :: alt1(:)
  real, allocatable :: cfxup1_1(:),cfxup2_1(:),cfxdn1_1(:),    &
       dfxup1_1(:),dfxup2_1(:),dfxdn1_1(:),    &
       efxup1_1(:),efxup2_1(:),efxdn1_1(:),    &
       tke_1(:)
  logical           :: ramsflg,awrfflg

  real, allocatable :: cfxup1_2(:),cfxup2_2(:),cfxdn1_2(:),    &
       dfxup1_2(:),dfxup2_2(:),dfxdn1_2(:),    &
       efxup1_2(:),efxup2_2(:),efxdn1_2(:),    &
       tke_2(:)
  ! JCL:temporary holding arrays for TL & SIGW
  real, allocatable :: tl2(:),sigw2(:)

  ! CHG(09/18/03) add mass flux fields for DMASS calculation
  real, allocatable :: umf(:,:,:),vmf(:,:,:),wmf(:,:,:)

  ! JCL:(4/5/02) parameters to convert P=>Density for mass violation claculation
  !     gravity (m/s2)       dry air (J/Kg-K)   mb to j/m3
  !dwen(20090822) ***********
  real             :: grav,rdry,p2jm
  DATA GRAV/9.80616/,  RDRY/287.04/,      P2JM/100.0/

  !dwen ****************************************

  INTEGER           :: i,j,k,kret,kzmix,kbls,kk

  REAL              :: offset      ! pressure offset for NMM
  REAL              :: zmdlt,gspdx,gspdy,hbot,xbot  
  REAL              :: pblh,ebot,shtf,tvmix,zlocn,zmix,zmixnew,        &
       gspd_x,gspd_y,ztt1,ztt2,ztt3,ztt4,ztt,          &
       fmap_x,fmap_y,hdivu,hdivv,hdiv,wmfs1,wmfs2,     &
       dz,vdiv,zsigw_dn 
  real              :: zsigw_up,dzsig,dudxup,dvdyup,rhobar,wup,wdn,    &
       dudxdn,dvdydn,rhodn,rhoup
  LOGICAL           :: mixd,hflx,eflx,uflx,ustx,tstx,tken,velv
  LOGICAL           :: zflg,qflg,uflg,tflg,prss,shgt,dzdt

  !-------------------------------------------------------------------------------
  ! external variables
  !-------------------------------------------------------------------------------


  INCLUDE 'DEFARG4.INC' ! interface statements for local subroutines

  !-------------------------------------------------------------------------------

  ! place structure in temporary variables   

  MIXD = DREC(KG,KT)%MIXD
  ZFLG = DREC(KG,KT)%ZFLG
  QFLG = DREC(KG,KT)%QFLG
  UFLG = DREC(KG,KT)%UFLG
  TFLG = DREC(KG,KT)%TFLG
  PRSS = DREC(KG,KT)%PRSS
  SHGT = DREC(KG,KT)%SHGT
  EFLX = DREC(KG,KT)%EFLX
  UFLX = DREC(KG,KT)%UFLX
  HFLX = DREC(KG,KT)%HFLX
  USTX = DREC(KG,KT)%USTR
  TSTX = DREC(KG,KT)%TSTR
  TKEN = DREC(KG,KT)%TKEN
  VELV = DREC(KG,KT)%VELV
  KBLS = DREC(KG,KT)%KBLS
  DZDT = DREC(KG,KT)%DZDT
  ZMDLT= DREC(KG,KT)%ZMDLT
  KZMIX= DREC(KG,KT)%KZMIX
  TVMIX= DREC(KG,KT)%TVMIX

  !-------------------------------------------------------------------------------
  ! allocate working arrays
  !-------------------------------------------------------------------------------

  IF(.NOT.ALLOCATED(PSG))THEN
     ALLOCATE (ESG(NZS),PSG(NZS),STAT=kret)
     IF(kret.NE.0)THEN
        WRITE(*,*)'*ERROR* prfcom: memory allocation - sigma' 
        STOP 900
     END IF

     ALLOCATE (U1(NZS),V1(NZS),W1(NZS),T1(NZS),Q1(NZS),P1(NZS),STAT=kret)
     IF(kret.NE.0)THEN
        WRITE(*,*)'*ERROR* prfcom: memory allocation - input' 
        STOP 900
     END IF

     ALLOCATE (E1(NZS),H1(NZS),X1(NZS),STAT=kret)
     IF(kret.NE.0)THEN
        WRITE(*,*)'*ERROR* prfcom: memory allocation - turbulence' 
        STOP 900
     END IF

     ALLOCATE (U2(NLVL),V2(NLVL),W2(NLVL),T2(NLVL),Q2(NLVL),P2(NLVL),STAT=kret)
     IF(kret.NE.0)THEN
        WRITE(*,*)'*ERROR* prfcom: memory allocation - output' 
        STOP 900
     END IF

     ALLOCATE (E2(NLVL),H2(NLVL),X2(NLVL),STAT=kret)
     IF(kret.NE.0)THEN
        WRITE(*,*)'*ERROR* prfcom: memory allocation - mixing' 
        STOP 900
     END IF
     E2=0.0
     H2=0.0
     X2=0.0

     ALLOCATE (A2(NLVL),Z2(NLVL),D2(NLVL),STAT=kret)
     IF(kret.NE.0)THEN
        WRITE(*,*)'*ERROR* prfcom: memory allocation - diagnostic' 
        STOP 900
     END IF

     !dwen(20090803) *************************
     allocate (alt1(nzs),cfxup1_1(nzs),cfxup2_1(nzs),cfxdn1_1(nzs),    &
          dfxup1_1(nzs),dfxup2_1(nzs),dfxdn1_1(nzs),efxup1_1(nzs), &
          efxup2_1(nzs),efxdn1_1(nzs),tke_1(nzs),stat=kret)
     if(kret.NE.0) then
        write(*,*) '*ERROR* prfcom: memory allocation - diagnostic'
        stop 900
     end if


     allocate (xv(nlvl),cfxup1_2(nlvl),cfxup2_2(nlvl),cfxdn1_2(nlvl),    &
          dfxup1_2(nlvl),dfxup2_2(nlvl),dfxdn1_2(nlvl),efxup1_2(nlvl), &
          efxup2_2(nlvl),efxdn1_2(nlvl),tke_2(nlvl),stat=kret)
     if(kret.NE.0) then
        write(*,*) '*ERROR* prfcom: memory allocation - diagnostic'
        stop 900
     end if


     allocate (tl2(nlvl),sigw2(nlvl),    &
          umf(nxs,nys,nlvl),vmf(nxs,nys,nlvl),wmf(nxs,nys,nlvl),stat=kret) 
     if(kret.NE.0) then
        write(*,*) '*ERROR* prfcom: memory allocation - diagnostic'
        stop 900
     end if

  END IF


  !-------------------------------------------------------------------------------
  ! remap vertical input data heights into input array
  !-------------------------------------------------------------------------------
  !dwen(20090803) ***********************
  ! JCL(03/27/03): flag specifying whether data from RAMS or not
  RAMSFLG=GRID(KG,kt)%MODEL_ID.EQ.'RAMS'

  ! JCL(050725): flag specifying whether data from WRFFLG or not
  AWRFFLG=GRID(KG,kt)%MODEL_ID(2:4).EQ.'WRF'
  !******************************************

  DO K=1,NZS
     PSG(K)=DREC(KG,KT)%HEIGHT(K+1)
  END DO

  OFFSET=GRID(KG,KT)%DUMMY

  !-------------------------------------------------------------------------------
  ! process each node on subgrid
  !-------------------------------------------------------------------------------

  DO J=1,NYS
     DO I=1,NXS

        !    wind speed conversion (m/s) to (grid/min)
        GSPDX=60.0/GX(I,J)
        GSPDY=60.0/GY(I,J)

        !-------------------------------------------------------------------------------
        !    place subgrid profile into temporary variables
        !-------------------------------------------------------------------------------

        !    special preprocessing to fill in below-ground zero TKE values (ETA model)
        IF(TKEN)THEN
           K=1
           DO WHILE (K.LE.NZS)
              IF(E(I,J,K).GT.0.0)THEN
                 EBOT=E(I,J,K)
                 K=NZS
              END IF
              K=K+1
           END DO

           !    special preprocessing to fill in below-ground zero variance values
        ELSEIF(VELV)THEN
           K=1
           DO WHILE (K.LE.NZS)
              IF(E(I,J,K).GT.0.0)THEN
                 EBOT=E(I,J,K)
                 HBOT=H(I,J,K)
                 XBOT=X(I,J,K)
                 K=NZS
              END IF
              K=K+1
           END DO
        END IF

        datain : DO K=1,NZS

           !       each level should have u,v,t(ambient)
           U1(K) = U(I,J,K)
           V1(K) = V(I,J,K)
           T1(K) = A(I,J,K)

           !       if zflg is set then field contains valid pressure data
           !       sigma coordinates    - pressure is computed from hyspometric eq
           !       pressure coordinates - contains height the replaced by pressure
           P1(K) = P(I,J,K)

           !       vertical motion field may not be present at all levels
           IF(DREC(KG,KT)%WFLG(K))THEN
              W1(K) = W(I,J,K)
           ELSE
              W1(K) = 0.0
           END IF

           !       humidity field may not be present at all levels
           IF(DREC(KG,KT)%RFLG(K))THEN
              Q1(K) = Q(I,J,K)
           ELSE
              Q1(K) = 0.0
           END IF

           !(dwen20090803) **************************
           !dwen                copied from STILT codes
           if (awrfflg .and. (deepflg .or. shallflg)) tke_1(k) = tke(i,j,k)
           if (awrfflg .and. deepflg) then
              CFXUP1_1(k) = CFXUP1(i,j,k)
              CFXDN1_1(k) = CFXDN1(i,j,k)
              DFXUP1_1(k) = DFXUP1(i,j,k)
              EFXUP1_1(k) = EFXUP1(i,j,k)
              DFXDN1_1(k) = DFXDN1(i,j,k)
              EFXDN1_1(k) = EFXDN1(i,j,k)
           endif
           if (awrfflg .and. shallflg) then
              CFXUP2_1(k) = CFXUP2(i,j,k)
              DFXUP2_1(k) = DFXUP2(i,j,k)
              EFXUP2_1(k) = EFXUP2(i,j,k)
           endif
           !           only needed for WRF:
           if (awrfflg) then
              alt1(k)=alt0(i,j,k)
           else
              alt1(k)=0.0
           endif

           !         END DO

           !dwen **************************************



           !       skip remaining variables for non-mixing applications
           IF(.NOT.VMIX) CYCLE datain

           !       turbulent kinetic energy field may not be present at all levels
           IF(TKEN)THEN
              E1(K) = E(I,J,K)
              IF(E1(K).EQ.0.0) E1(K)=EBOT
           ELSEIF(VELV)THEN
              E1(K) = E(I,J,K)
              H1(K) = H(I,J,K)
              X1(K) = X(I,J,K)
              IF(E1(K).EQ.0.0) THEN
                 E1(K)=EBOT
                 H1(K)=HBOT
                 X1(K)=XBOT
              END IF
           END IF

        END DO datain

        !-------------------------------------------------------------------------------
        !    remap data according to coordinate system
        !-------------------------------------------------------------------------------

        !dwen(20090804)     SELECT CASE (DREC(KG,KT)%Z_FLAG)

        !dwen(20090804)     CASE (1)
        !dwen(20090804) ************************
        ! TN (20050824): add wrf vertical coordinate system
        IF( AWRFFLG )THEN
           !dwen(20090804) use if statement instead of select due to awrfflg
           !              input data on WRF pressure-sigma system
           !dwen(20090804)  rename TKEN used in PRFERF to TKE because there is a variable named TKEN in PRFCOM 
           !dwen               CALL PRFwrf(VMIX,DREC(KG)%QFLG,DREC(KG)%UFLG,            &
           !dwen                 DREC(KG)%TFLG,DREC(KG)%PRSS,DREC(KG)%SHGT,ZT(I,J),    &
           !dwen                 P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),w0(i,j),alt1, &
           !dwen                 mu(i,j),muu(i,j),muv(i,j),msfu(i,j),msfv(i,j),msft(i,j),fluxflg,&
           !dwen                 NZS,PSG,P1,U1,V1,W1,T1,Q1,NLVL,ZMDL,ZSG,              &
           !dwen                 P2,U2,V2,W2,T2,Z2,Q2,D2, &
           !dwen                 deepflg, shallflg, &
           !dwen           CFXUP1_1,       &
           !dwen           CFXUP2_1,CFXDN1_1,DFXUP1_1,         &
           !dwen           DFXUP2_1,EFXUP1_1,EFXUP2_1,         &
           !dwen           DFXDN1_1,EFXDN1_1,TKE_1,           &
           !dwen           CFXUP1_2,       &
           !dwen           cFXUP2_2,CFXDN1_2,DFXUP1_2,         &
           !dwen           DFXUP2_2,EFXUP1_2,EFXUP2_2,         &
           !dwen           DFXDN1_2,EFXDN1_2,TKE_2 )

           CALL PRFwrf(VMIX,QFLG,UFLG,            &
                TFLG,PRSS,SHGT,ZT(I,J),    &
                P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),w0(i,j),alt1, &
                mu(i,j),muu(i,j),muv(i,j),msfu(i,j),msfv(i,j),msft(i,j),fluxflg,&
                NZS,PSG,P1,U1,V1,W1,T1,Q1,NLVL,ZMDL,ZSG,              &
                P2,U2,V2,W2,T2,Z2,Q2,D2, &
                deepflg, shallflg, &
                CFXUP1_1,       &
                CFXUP2_1,CFXDN1_1,DFXUP1_1,         &
                DFXUP2_1,EFXUP1_1,EFXUP2_1,         &
                DFXDN1_1,EFXDN1_1,TKE_1,           &
                CFXUP1_2,       &
                cFXUP2_2,CFXDN1_2,DFXUP1_2,         &
                DFXUP2_2,EFXUP1_2,EFXUP2_2,         &
                DFXDN1_2,EFXDN1_2,TKE_2 )
           !dwen ****************************************

           !dwen(20090804)     CASE (1)

        ELSEIF( (DREC(KG,KT)%Z_FLAG).EQ.1)THEN
           !       input data on pressure-sigma system
           CALL PRFSIG(VMIX,TKEN,VELV,ZFLG,QFLG,UFLG,TFLG,PRSS,DZDT,OFFSET,       &
                ZT(I,J),P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),                  &
                NZS,PSG,P1,U1,V1,W1,T1,Q1,E1,H1,X1,NLVL,ZMDL,ZSG,                 &
                P2,U2,V2,W2,T2,Z2,Q2,D2,A2,E2,H2,X2)

           !dwen(20090804)     CASE (2)
        ELSEIF( (DREC(KG,KT)%Z_FLAG).EQ.2)THEN
           !       input data on pressure-absolute system
           CALL PRFPRS(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PRSS,SHGT,                   &
                ZT(I,J),P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),                  &
                NZS,PSG,P1,U1,V1,W1,T1,Q1,E1,H1,X1,NLVL,ZMDL,ZSG,                 &
                P2,U2,V2,W2,T2,Z2,Q2,D2,A2,E2,H2,X2)

           !dwen(20090804)     CASE (3) 

        ELSEIF( (DREC(KG,KT)%Z_FLAG).EQ.3)THEN 
           !       input data on terrain-sigma system
           !dwen(20090804) *********************************
           !               add ramsflg in PRFTER according to STILT codes
           !               add zmdlt according to HYSPLIT4.9 codes, zmdlt is assigned in runset.f
           CALL PRFTER(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PRSS,SHGT,                   &
                ZT(I,J),ZMDLT,P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),            &
                NZS,PSG,P1,U1,V1,W1,T1,Q1,E1,H1,X1,NLVL,ZMDL,ZSG,                 &
                P2,U2,V2,W2,T2,Z2,Q2,D2,A2,E2,H2,X2,ramsflg)

           !dwen        CALL PRFTER(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PRSS,SHGT,                   &
           !dwen             ZT(I,J),ZMDLT,P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),            &
           !dwen             NZS,PSG,P1,U1,V1,W1,T1,Q1,E1,H1,X1,NLVL,ZMDL,ZSG,                 &
           !dwen             P2,U2,V2,W2,T2,Z2,Q2,D2,A2,E2,H2,X2)
           !dwen **************************************************



           !dwen(20090804)     CASE (4)
        ELSEIF( (DREC(KG,KT)%Z_FLAG).EQ.4)THEN 
           !       input data on ecmwf hybrid sigma system
           CALL PRFECM(VMIX,TKEN,VELV,QFLG,UFLG,TFLG,PRSS,SHGT,                   &
                ZT(I,J),P0(I,J),U0(I,J),V0(I,J),T0(I,J),Z0(I,J),                  &
                NZS,PSG,ESG,U1,V1,W1,T1,Q1,E1,H1,X1,NLVL,ZMDL,ZSG,                &
                P2,U2,V2,W2,T2,Z2,Q2,D2,A2,E2,H2,X2)

           !dwen(20090804)     CASE DEFAULT
        ELSE
           WRITE(*,*)'*ERROR* prfcom: no vertical data remapping'
           WRITE(*,*)' Unknown vertical type:', DREC(KG,KT)%Z_FLAG
           STOP 900

           !dwen(20090804)     END SELECT
        END IF

        !-------------------------------------------------------------------------------
        !    stability analysis for concentration simulations
        !-------------------------------------------------------------------------------

        IF(VMIX)THEN
           !       correction for sign reversal in ETA/EDAS sensible heat flux
           SHTF=HF(I,J)
           IF(GRID(KG,KT)%MODEL_ID.EQ.' ETA'.OR.                                  &
                GRID(KG,KT)%MODEL_ID.EQ.'EDAS')SHTF=-HF(I,J)

           !       analyze surface stability (z/L) from fluxes or wind/temp profile
           !dwen(20090804)        CALL STBANL(KBLS,I,J,KSFC,KMIX0,KMIXD,MIXD,TKEN,HFLX,EFLX,UFLX,USTX,   &
           !dwen(20090804)                    TSTX,Z0(I,J),UF(I,J),VF(I,J),SHTF,NLVL,U2,V2,T2,Z2,E2,D2,  &
           !dwen(20090804)                    ZI(I,J),SF(I,J),SS(I,J))

           !dwen(20090804) change STBANL according STILT code

           ! CHG:(11/19/01) pass on humidity (RH fraction) to subroutine STBANL as well
           ! CHG:(12/04/01) get lim. of convection height (ZLOCN) from stbanl
           !           analyze surface stability (z/L) from fluxes  if available
           ! CHG:(9/17/02) add 'ICONVECT' as convection flag

           !            CALL STBANL(KSFC,DREC(KG)%EFLX,DREC(KG)%HFLX,DREC(KG)%UFLX, &
           !     &         DREC(KG)%USTR,DREC(KG)%TSTR,                             &
           !     &         Z0(I,J),UF(I,J),VF(I,J),SF(I,J),NLVL,U2,V2,T2,Z2,D2,Q2,  &
           !     &         USTR,TSTR,WSTR,ZMIX,SLEN,PSI,ZLOCN,ICONVECT)
           !dwen(20090804)
           !                                     STILT      HYSPLIT4.9
           ! sensible heat flux                   SF            SHTF
           ! turbulent kinetic Energy                            E
           ! relative humidity                    Q
           ! mixed layer depth                   ZMIX           ZI
           ! integrated stability function heat  PSI            SF
           ! static stability parameter                         SS
           ! friction veloctiy                   USTR
           ! friction temperature                TSTR
           ! convective velocity scale           WSTR
           ! Obukhov stability length            SLEN
           ! note that SF in STBANL here represents stability function, not sensible heat flux
           !     ustr,tstr,wstr,slen,vstr,zmix,psi were put into COMMON block in HYSPLIT4.9

           CALL STBANL(KBLS,I,J,KSFC,KMIX0,KMIXD,MIXD,TKEN,HFLX,EFLX,UFLX,USTX,   &
                TSTX,Z0(I,J),UF(I,J),VF(I,J),SHTF,NLVL,U2,V2,T2,Z2,E2,D2,  &
                ZI(I,J),SF(I,J),SS(I,J), &
                q2,zlocn,iconvect,ziscale)


           IF(DREC(KG,KT)%KBLT.EQ.1)THEN
              !          compute mixing Beljaars-Holtslag 

              !dwen(20090804) modify STBSND according to current version of STILT codes
              ! JCL:      add temporary holding arrays TL2&SIGW2 to hold results of
              !             calculations for Lagrangian timescale & stddev of W
              ! JCL:      add the roughness length as argument as well
              !dwen(20090804)            CALL STBSND(ISOT,KSFC,TSTR,USTR,WSTR,SLEN,ZMIX,             &
              !dwen(20090804)     &         NLVL,U2,V2,T2,Z2,XV,TL2,SIGW2,Z0(I,J))
              !dwen(20090806)    add xv as vertical mixing coefficient

              !           
              !dwen(20090804)  ustr,vstr,tstr,wstr,slen,zmix were put into common block in HYSPLIT4.9

              CALL STBSND(TKERD,TKERN,KZMIX,TVMIX,KSFC,NLVL,U2,V2,T2,Z2,E2,H2,X2, &
                   xv,tl2,sigw2,z0(i,j))

              !dwen(20090804)           CALL STBSND(TKERD,TKERN,KZMIX,TVMIX,KSFC,NLVL,U2,V2,T2,Z2,E2,H2,X2)

           ELSEIF(DREC(KG,KT)%KBLT.EQ.2)THEN
              !          compute mixing Kantha-Clayson 
              CALL STBVAR(TKERD,TKERN,KZMIX,TVMIX,KSFC,NLVL,U2,V2,T2,Z2,E2,H2,X2)

           ELSEIF(DREC(KG,KT)%KBLT.EQ.3)THEN
              !          velocity variances from input TKE
              CALL STBTKE(TKERD,TKERN,KZMIX,TVMIX,KSFC,NLVL,U2,V2,   Z2,E2,H2,X2)

           ELSEIF(DREC(KG,KT)%KBLT.EQ.4)THEN
              !          velocity variances from input meteorological data file
              !          no processing is required
              CONTINUE

           ELSE
              WRITE(*,*)'*ERROR* prfcom: turbulence method undefined'
              WRITE(*,*)' Unknown selection: ',DREC(KG,KT)%KBLT
              STOP 900
           END IF

           ! CHG:(12/04/01)store lim. of conv. height in array
           ZLOC(I,J)=ZLOCN
           !dwen *********************************************

           !       component friction velocities (reverse of 10m wind vectors)
           UF(I,J)=-UF(I,J)*GSPDX
           VF(I,J)=-VF(I,J)*GSPDY

        END IF

        !dwen(20090804) **********************
        !               copied from STILT codes
        ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
        ! JCL:   need following line to prevent GD(I,J)=0, causing model to crash
        !        because would be dividing by 0!
        IF(GX(I,J).EQ.0)GX(I,J)=GX(I,J-1)
        IF(GY(I,J).EQ.0)GY(I,J)=GY(I,J-1)

        !        wind speed conversion (m/s) to (grid/min)
        ! CHG(09/18/03) use interpolated grid size for RAMS
        IF(.NOT.RAMSFLG)THEN
           ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
           GSPD_X=60.0/GX(I,J)
           GSPD_Y=60.0/GY(I,J)
           if (awrfflg .and. fluxflg) then
              ! GX is valid at the mass point; GX*msfu/msft is the gx at the staggered u-point
              gspd_x = gspd_x * msft(i,j)/msfu(i,j)
              gspd_y = gspd_y * msft(i,j)/msfv(i,j)
           end if
        ELSE
           ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
           !at u flux position i
           IF(I.LT.NXS)GSPDX=60.0*2.0/(GX(I+1,J)+GX(I,J))
           !at v flux position j
           IF(J.LT.NYS)GSPDY=60.0*2.0/(GY(I,J+1)+GY(I,J))
           IF(I.EQ.NXS)GSPDX=60.0/GX(I,J)
           IF(J.EQ.NYS)GSPDY=60.0/GY(I,J)
        END IF
        !dwen  **********************************************

        !-------------------------------------------------------------------------------
        !    put adjusted data back into 3d array now interpolated to internal grid
        !-------------------------------------------------------------------------------

        dataout : DO K=1,NLVL
           !dwen(20090804) ****************
           !               copied from STILT codes
           if(.not.ramsflg) then
              U(I,J,K)=U2(K)*GSPD_X
              V(I,J,K)=V2(K)*GSPD_Y
           ELSE
              UMF(I,J,K)=U2(K)
              VMF(I,J,K)=V2(K)
              WMF(I,J,K)=W2(K)
              !dwen ***************************
              U(I,J,K)=U2(K)*GSPDX                 ! wind speed to grid / min
              V(I,J,K)=V2(K)*GSPDY
              !dwen(20090804) *****************
           end if
           !dwen***************************
           W(I,J,K)=W2(K)*60.0                  ! convert dp/dt to per minute
           T(I,J,K)=T2(K)                       ! converted to potential
           A(I,J,K)=A2(K)                       ! local ambient temperature
           Q(I,J,K)=MAX(0.0,MIN(1.0,Q2(K)))     ! humidity is relative fraction
           P(I,J,K)=P2(K)                       ! local pressure (mb)
           !dwen(20090804) *******************
           !           local air density (kg/m3)
           D(I,J,K)=D2(K)
           !           vertical mixing coefficient (m2/s)
           X(I,J,K)=XV(K)
           ! JCL:      Lagrangian timescale (s)
           TL(I,J,K)=TL2(K)
           ! JCL:      stddev of vertical velocity (m/s)
           SIGW(I,J,K)=SIGW2(K)
           if (awrfflg .and. (deepflg .or. shallflg)) tke(i,j,k) = tke_2(k)
           if (awrfflg .and. deepflg) then
              CFXUP1(i,j,k) = CFXUP1_2(k)
              CFXDN1(i,j,k) = CFXDN1_2(k)
              DFXUP1(i,j,k) = DFXUP1_2(k)
              EFXUP1(i,j,k) = EFXUP1_2(k)
              DFXDN1(i,j,k) = DFXDN1_2(k)
              EFXDN1(i,j,k) = EFXDN1_2(k)
           endif
           if (awrfflg .and. shallflg) then
              CFXUP2(i,j,k) = CFXUP2_2(k)
              DFXUP2(i,j,k) = DFXUP2_2(k)
              EFXUP2(i,j,k) = EFXUP2_2(k)
           endif
           !dwen **************************************

           IF(.NOT.VMIX) CYCLE dataout
           H(I,J,K)=MAX(0.0,H2(K))              ! U component turbulence (m2/s2)
           E(I,J,K)=MAX(0.0,E2(K))              ! V component turbulence (m2/s2) 
           X(I,J,K)=MAX(0.0,X2(K))              ! W component turbulence (m2/s2)
        END DO dataout

     END DO
  END DO

  !-------------------------------------------------------------------------------
  ! Deformation method can be used to compute the horizontal mixing 
  !-------------------------------------------------------------------------------

  !dwen(2009084) *******************************
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! JCL(4/3/02): following lines quantify MASS VIOLATION
  !     NOTE:  still need to add the (Drho/Dt)s term--will be done in ADVPNT
  !     loop over each grid position
  !     mass violation array has one less element in each dimension,
  !         except for the VERTICAL (b/c impose no-slip b.c. at the surface)
  ! CHG(09/18/03) use mass fluxes
  Model: IF (RAMSFLG) THEN

     DO K=0,NLVL-1
        DO J=1,NYS-1
           DO I=1,NXS-1
              ! to get from mass flux density to mass flux (kg/s), need proper area for hor. fluxes
              ! need to interpolate terrain height to locations
              !at u-flux location I+1
              ZTT1=ZMDL-(ZT(I+1,J+1)+ZT(MIN0(I+2,NXS),J+1))/2
              !at u-flux location I
              ZTT2=ZMDL-(ZT(I,J+1)+ZT(I+1,J+1))/2
              !at v-flux location J+1
              ZTT3=ZMDL-(ZT(I+1,J+1)+ZT(I+1,MIN0(J+2,NYS)))/2
              !at v-flux location J
              ZTT4=ZMDL-(ZT(I+1,J)+ZT(I+1,J+1))/2
              ZTT=ZMDL-ZT(I+1,J+1)

              ! horizontal divergence
              ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
              ! CHG need to apply correct GD(i,j) ("delta-x" etc.)
              !       FMAP=45000.0/GD(I+1,J+1)
              ! WJS (4/1/05): Make grid size dynamic; multiply by 1000 to convert to meters
              FMAP_X=(GRID(KG,kt)%SIZE * 1000.0)/GX(I+1,J+1)
              FMAP_Y=(GRID(KG,kt)%SIZE * 1000.0)/GY(I+1,J+1)
              HDIVU=(UMF(I+1,J+1,K+1)/ZTT-UMF(I,J+1,K+1)/ZTT)*ZMDL
              HDIVV=(VMF(I+1,J+1,K+1)/ZTT-VMF(I+1,J,K+1)/ZTT)*ZMDL

              ! JCL:(07/14/2004) split up GD => GX & GY for global grids (non-conformal)
              ! WJS (4/1/05): Make grid size dynamic; multiply by 1000 to convert to meters
              HDIVU=HDIVU*FMAP_X*FMAP_X/(GRID(KG,kt)%SIZE * 1000.0)
              HDIVV=HDIVV*FMAP_Y*FMAP_Y/(GRID(KG,kt)%SIZE * 1000.0)

              HDIV=HDIVU+HDIVV
              ! CHG use different map factors, specific for each flux location
              ! ->only minimal changes, not significant...
              !w* fluxes
              WMFS1=-WMF(I+1,J+1,K+1)*ZMDL
              IF(K.GT.0)WMFS2=-WMF(I+1,J+1,K)*ZMDL
              IF(K.EQ.0)WMFS2=0.0
              !z* altitude difference
              IF(K.GT.0)DZ=-1.0*(ZSG(K+1)-ZSG(K))*ZMDL
              IF(K.EQ.0)DZ=-1.0*(ZSG(K+1)-1.0)*ZMDL
              VDIV=(WMFS1-WMFS2)/DZ
              !also convert to 1/minute
              DELMASS(I+1,J+1,K+1)=(HDIV+VDIV)*60.0
              !of DO I=1,NXS-1
           END DO
           !of DO J=1,NYS-1
        END DO
        !of DO K=0,NLVL-1
     END DO
     ! for awrf:
  ELSE IF (awrfflg) THEN Model
     ! C-grid setup in WRF: for x-staggered u and unstaggered m grid points:
     !                        i-1  i   i  i+1 i+1
     !                         m   u   m   u   m
     !
     !                      for z-staggered w and unstaggered m (and u,v) grid points:
     !                        k-1  k   k  k+1 k+1 k+2
     !                         m   w   m   w   m   w
     DO K=1,NLVL
        DO J=1,NYS
           DO I=1,NXS
              ! couple decoupled velocities at staggered u,v-grid points with density
              UMF(i,j,k)=0.5*(D(max(i-1,1),j,k)+D(i,j,k))*u(i,j,k)
              VMF(i,j,k)=0.5*(D(i,max(j-1,1),k)+D(i,j,k))*v(i,j,k)
              ! couple decoupled w at staggered w-grid points with density
              WMF(i,j,k)=0.5*(D(i,j,max(1,k-1))+D(i,j,k))*w(i,j,k)
           enddo
        enddo
     enddo
     DO K=0,NLVL-1
        ! Compute depth of layer bounded by w-levels:
        if (k .le. 0) then
           ! kk=k+1=1: layer bounded by zsigw(1) and zsigw(2)
           zsigw_dn=1.
           zsigw_up=2*zsg(k+1) - 1. !zs of second w-level
        else
           ! kk=k+1: layer bounded by zsigw(kk) and zsigw(kk+1)
           zsigw_up=2.*zsg(k+1)-zsigw_dn  !zs of w-level kk+1
        endif
        DZSIG=zsigw_up-zsigw_dn
        zsigw_dn=zsigw_up
        DO J=1,NYS-1
           DO I=1,NXS-1
              ! compute zmdl-terrain height
              !at u-flux location I+1
              ZTT1=ZMDL-(ZT(I,J)+ZT(MIN(I+1,NXS),J))/2
              !at u-flux location I
              ZTT2=ZMDL-(ZT(max(I-1,1),J)+ZT(I,J))/2
              !at v-flux location J+1
              ZTT3=ZMDL-(ZT(I,J)+ZT(I,MIN(J+1,NYS)))/2
              !at v-flux location J
              ZTT4=ZMDL-(ZT(I,max(J-1,1))+ZT(I,J))/2
              ZTT=ZMDL-ZT(I,J)
              !        hor divergence term (not have to divide by GD b/c already done so in wind speed conversion above)
              !        (computed at mass level,
              DUDXUP=UMF(I+1,J,K+1)*ztt1 - UMF(I,J,K+1)*ztt2
              DVDYUP=VMF(I,J+1,K+1)*ztt3 - VMF(I,J,K+1)*ztt4
              rhobar=d(i,j,k+1)
              if (k+1 .ge. NLVL) then
                 wup = 0         !upper b.c. for w
              else
                 wup=wmf(i,j,k+2)
              endif
              wdn=wmf(i,j,k+1)
              ! same as for other non-RAMS models (see below)         
              HDIV=(DUDXUP+DVDYUP)/ZTT
              VDIV=(WUP-WDN)/DZSIG
              DELMASS(I,J,K+1)=(HDIV+VDIV)/RHOBAR
           enddo
        enddo
     enddo
     ! for other than RAMS or awrf
  ELSE Model
     UMF(1:NXS,1:NYS,1:NLVL)=D(1:NXS,1:NYS,1:NLVL)*                  &
          &                          U(1:NXS,1:NYS,1:NLVL)
     VMF(1:NXS,1:NYS,1:NLVL)=D(1:NXS,1:NYS,1:NLVL)*                  &
          &                          V(1:NXS,1:NYS,1:NLVL)
     WMF(1:NXS,1:NYS,1:NLVL)=D(1:NXS,1:NYS,1:NLVL)*                  &
          &                          W(1:NXS,1:NYS,1:NLVL)
     DO K=0,NLVL-1
        DO J=1,NYS-1
           DO I=1,NXS-1
              ZTT1=ZMDL-ZT(I,J)
              ZTT2=ZMDL-ZT(I+1,J)
              ZTT3=ZMDL-ZT(I,J+1)
              ZTT4=ZMDL-ZT(I+1,J+1)
              ZTT=(ZTT1+ZTT2+ZTT3+ZTT4)/4.0
              !        hor divergence term(not have to divide by GD b/c already done so in wind speed conversion above
              DUDXUP=(UMF(I+1,J+1,K+1)*ZTT4+                                 &
                   &           UMF(I+1,J,K+1)*ZTT2)/2.0-                              &
                   &          (UMF(I,J+1,K+1)*ZTT3+                                   &
                   &           UMF(I,J,K+1)*ZTT1)/2.0
              DVDYUP=(VMF(I,J+1,K+1)*ZTT3+                                   &
                   &           VMF(I+1,J+1,K+1)*ZTT4)/2.0-                            &
                   &          (VMF(I+1,J,K+1)*ZTT2+                                   &
                   &           VMF(I,J,K+1)*ZTT1)/2.0
              WUP=(WMF(I,J,K+1)+WMF(I+1,J,K+1)                               &
                   &       +WMF(I,J+1,K+1)+WMF(I+1,J+1,K+1))/4.0

              !        NOT near ground
              IF(K.GT.0)THEN
                 !           horizontal divergence in lower model layer
                 DUDXDN=(UMF(I+1,J+1,K)*ZTT4+                                &
                      &           UMF(I+1,J,K)*ZTT2)/2.0-                                &
                      &          (UMF(I,J+1,K)*ZTT3+                                     &
                      &           UMF(I,J,K)*ZTT1)/2.0
                 DVDYDN=(VMF(I,J+1,K)*ZTT3+                                  &
                      &           VMF(I+1,J+1,K)*ZTT4)/2.0-                              &
                      &          (VMF(I+1,J,K)*ZTT2+                                     &
                      &           VMF(I,J,K)*ZTT1)/2.0
                 !           vertical velocity at lower model layer
                 WDN=(WMF(I,J,K)+WMF(I+1,J,K)                                &
                      &       +WMF(I,J+1,K)+WMF(I+1,J+1,K))/4.0
                 RHODN=(D(I,J,K)+D(I+1,J,K)+D(I,J+1,K)+D(I+1,J+1,K))/4.0
                 RHOUP=(D(I,J,K+1)+D(I+1,J,K+1)+                             &
                      &             D(I,J+1,K+1)+D(I+1,J+1,K+1))/4.0
                 RHOBAR=(RHODN+RHOUP)/2.0
                 DZSIG=ZSG(K+1)-ZSG(K)
                 !        near ground
              ELSE
                 !           no-slip b.c. at ground surface
                 DUDXDN=0.0
                 DVDYDN=0.0
                 WDN=0.0
                 RHOUP=(D(I,J,K+1)+D(I+1,J,K+1)+                             &
                      &             D(I,J+1,K+1)+D(I+1,J+1,K+1))/4.0
                 ! JCL:(5/2/02) NGM data doesn't have surface temp data, so simply assign RHOUP to RHODN
                 IF(GRID(KG,kt)%MODEL_ID.EQ.'NGM')THEN
                    RHODN=RHOUP
                 ELSE
                    !              density at the surface--from Ideal Gas Law
                    RHODN=(P2JM/RDRY)*(P0(I,J)/T0(I,J)+P0(I+1,J)/T0(I+1,J)+  &
                         &            P0(I,J+1)/T0(I,J+1)+P0(I+1,J+1)/T0(I+1,J+1))/4.0
                 END IF
                 RHOBAR=(RHODN+RHOUP)/2.0
                 DZSIG=ZSG(K+1)-1.0
              END IF

              !        horizontal divergence term--average over upper and lower model layers
              ! fixed bug: was (DVDYDN+DVDYDN)
              HDIV=((DUDXDN+DUDXUP)/2.0+(DVDYDN+DVDYUP)/2.0)/ZTT
              !        vertical divergence term
              VDIV=(WUP-WDN)/DZSIG

              !        assign element to mass violation array; still need to add the (Drho/Dt)s term
              !        note convention: mass viol betw levels 1 & 2, e.g., gets stored in level 2 (K+1) of mass viol grid
              !        mass violation is in units of [fraction of mass in gridcell/min]
              DELMASS(I,J,K+1)=(HDIV+VDIV)/RHOBAR

              !of DO I=1,NXS-1
           END DO
           !of DO J=1,NYS-1
        END DO
        !of DO K=0,NLVL-1
     END DO
     !of IF(RAMSFLAG) ... ELSE ...
  END IF Model


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !dwen ****************************************


  IF(VMIX.AND.(DREC(KG,KT)%KDEF.EQ.1)) CALL STBHOR(NXS,NYS,NLVL,GX,GY,U,V,H,hm,E)
  xm(:,:,:) = 0. !not used other than to assign metz%wmxc in advmet
  DEALLOCATE (ESG,PSG)
  DEALLOCATE (U1,V1,W1,T1,Q1,P1)
  DEALLOCATE (E1,H1,X1)
  DEALLOCATE (U2,V2,W2,T2,Q2,P2)
  DEALLOCATE (E2,H2,X2)
  DEALLOCATE (A2,Z2,D2)
  deallocate (alt1,cfxup1_1,cfxup2_1,cfxdn1_1,    &
       dfxup1_1,dfxup2_1,dfxdn1_1,efxup1_1, &
       efxup2_1,efxdn1_1,tke_1)
  deallocate (xv,cfxup1_2,cfxup2_2,cfxdn1_2,    &
       dfxup1_2,dfxup2_2,dfxdn1_2,efxup1_2, &
       efxup2_2,efxdn1_2,tke_2)
  deallocate (tl2,sigw2,umf,vmf,wmf)

END SUBROUTINE prfcom
