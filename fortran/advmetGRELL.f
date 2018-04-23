!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  CONCALL           calls Convective redistribution using GRELL scheme
!   PRGMMR:    CHRISTOPH GERBIG   DATE:09-23-03
!
! ABSTRACT:
!   extracts convective fluxes for gridcell and calls CGRELL convective
!   redistribution of particles using updraft & downdraft mass
!   Comparable to ADVMET
! USAGE:  CALL ADVMETGRELL(X,Y,CFXUP1, CFXUP2, CFXDN1, DFXUP1, DFXUP2, EFXUP1,
!                          EFXUP2, DFXDN1, EFXDN1, RAUP1, RAUP2, RADN1)
!   INPUT ARGUMENT LIST:
!     x,y       - real  horizontal particle position (grid units)
!     ?FX*      - real  Cloud/Entrainment/Detrainment Mass Fluxes
!                       (Updraft, downdraft; 1:deep; 2: shallow)
!                       UNITS: kg/m^2s; convention: upward=positive
!                       TO GET MASSFLUX (kg/s):
!                       MASSFLUX=?FX* * DX*DY (up/dndraft & hor. fluxes)
!                       numbers represent ro*w*a w/ a=fractional coverage
!                       to get mass flux[kg/s]: ?FX* * area (gridcell)
!     RA*       - real  updraft area coverage (fract.),
!                       (Updraft, downdraft; 1:deep; 2: shallow)
!     NXS,NYS   - int   horizontal grid dimensions
!     NLVL      - int   number of vertical levels
!   OUTPUT ARGUMENT LIST:
!     GBLMET COMMON WITH VARIABLES DEFINED IN DEFMETz.INC
!
! $Id: advmetGRELL.f,v 1.2 2009/11/19 19:55:01 jel Exp $
!
!$$$
! CHG (09/23/03) Radius up/downdraft (RAUP1,RAUP2,RADN1) not yet from RAMS, so assign small value
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
! JCL:(07/12/2004) added cyclic boundary condition flag, global boundaries
! dwen(20090812)   use ADV3NT instead
! dwen(20090812)   remove NXS,NYS
! dwen(20090812)   fortran90 upgrade
! dwen(20090812)   add METz

!dwen(20090812)      SUBROUTINE ADVMETGRELL(X,Y,CFXUP1,CFXUP2,CFXDN1,DFXUP1,DFXUP2,    &
!dwen(20090812)     &  EFXUP1,EFXUP2,DFXDN1,EFXDN1,TKEN,NXS,NYS,NLVL,GLOBAL,NXP,NYP)
      SUBROUTINE ADVMETGRELL(METz,X,Y,CFLXUP1,CFLXUP2,CFLXDN1,DFLXUP1,DFLXUP2,    &
                       EFLXUP1,EFLXUP2,DFLXDN1,EFLXDN1,TKEG,NLVL,GLOBAL,NXP,NYP)

!dwen(20090812)      use module_defmeto
      IMPLICIT NONE 

      include 'DEFMETO.INC'
!------------------------------------------
!dwen(20090824)  TYPE(aset),INTENT(OUT)   :: meto       ! surface advection variables
  TYPE(bset),INTENT(OUT)   :: metz(:)       ! surface advection variables
  real,      intent(in)    :: x,y        
  real,      intent(in)    :: cflxup1(:,:,:),cflxup2(:,:,:)
  real,      intent(in)    :: cflxdn1(:,:,:),dflxup1(:,:,:)
  real,      intent(in)    :: dflxup2(:,:,:),dflxdn1(:,:,:)
  real,      intent(in)    :: eflxup1(:,:,:),eflxup2(:,:,:)
  real,      intent(in)    :: eflxdn1(:,:,:),tkeg(:,:,:)
  integer,   intent(in)    :: nlvl,nxp,nyp 
  character(2),   intent(in)    :: global
!--------------------------------------------
!   internal variables
!-------------------------------------------
  real       ::xx,yy,var1,zk
  integer    ::kl

 interface

  SUBROUTINE ADV3NT(S,XP,YP,ZX,SS,GLOBAL,NXP,NYP)
  IMPLICIT NONE
  REAL,      INTENT(IN)    :: s(:,:,:)      ! field for interpolation
  REAL,      INTENT(IN)    :: xp,yp         ! position of interpolated value
  REAL,      INTENT(IN)    :: zx            ! vertical interpolation fraction
  REAL,      INTENT(OUT)   :: ss            ! value of S at x1,y1,z1
  CHARACTER(2),   INTENT(IN)    :: global        ! global cyclic boundary conditions
  INTEGER,   INTENT(IN)    :: nxp,nyp       ! global boundary values
  END SUBROUTINE adv3nt

end interface

! CHG (09/23/03) Radius up/downdraft not yet from RAMS, so assign small value
! CHG(09/25/03) add RAMS turb. kin. energy TKEN


!dwen(20090825)      XX=DNINT(X)
!dwen(20090825)      YY=DNINT(Y)
      XX=aNINT(X)
      YY=aNINT(Y)
      DO KL=1,NLVL
!     set vertical interpolation point to index position
        ZK=KL
!dwen        CALL ADVINT(CFXUP1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen     &              VAR1)
        CALL ADV3NT(CFLXUP1,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%CFXUP1(KL)=VAR1
        METz(kl)%CFXUP1=VAR1
!dwen        CALL ADVINT(CFXUP2,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen    &              VAR1)
        CALL ADV3NT(CFLXUP2,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%CFXUP2(KL)=VAR1
        METz(kl)%CFXUP2=VAR1
!dwen        CALL ADVINT(CFXDN1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen      &              VAR1)
        CALL ADV3NT(CFLXDN1,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%CFXDN1(KL)=VAR1
        METz(kl)%CFXDN1=VAR1
!dwen        CALL ADVINT(DFXUP1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen     &              VAR1)
        CALL ADV3NT(DFLXUP1,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%DFXUP1(KL)=VAR1
        METz(kl)%DFXUP1=VAR1
!dwen        CALL ADVINT(DFXUP2,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen     &              VAR1)
        CALL ADV3NT(DFLXUP2,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%DFXUP2(KL)=VAR1
        METz(kl)%DFXUP2=VAR1
!dwen        CALL ADVINT(EFXUP1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen     &              VAR1)
        CALL ADV3NT(EFLXUP1,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
! CHG(09/24/03) need to shift?
                             !no shift
!dwen(20090826)        METz%EFXUP1(KL)=VAR1
        METz(kl)%EFXUP1=VAR1
!        IF(KL.LT.NLVL)METz%EFXUP1(KL+1)=VAR1 !leaves EFXUP1(1) zero...
!dwen        CALL ADVINT(EFXUP2,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen     &              VAR1)
        CALL ADV3NT(EFLXUP2,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%EFXUP2(KL)=VAR1
        METz(kl)%EFXUP2=VAR1
!dwen        CALL ADVINT(DFXDN1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen     &              VAR1)
        CALL ADV3NT(DFLXDN1,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%DFXDN1(KL)=VAR1
        METz(kl)%DFXDN1=VAR1
!dwen        CALL ADVINT(EFXDN1,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,  &
!dwen     &              VAR1)
        CALL ADV3NT(EFLXDN1,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%EFXDN1(KL)=VAR1
        METz(kl)%EFXDN1=VAR1
! CHG(09/25/03) add RAMS turb. kin. energy TKEN
!dwen        CALL ADVINT(TKEN,NXS,NYS,NZM,XX,YY,ZK,GLOBAL,NXP,NYP,    &
!dwen     &              VAR1)
        CALL ADV3NT(TKEG,XX,YY,ZK,var1,GLOBAL,NXP,NYP)
!dwen(20090826)        METz%TKEN(KL)=VAR1
        METz(kl)%TKEN=VAR1
      END DO

      RETURN
      END
