!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKINP           PAcK INPut converts char*1 to real
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACK INPUT DOES THE CONVERSION OF CHAR*1 PACKED ARRAY TO
!   A REAL*8*4 DATA ARRAY
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 14 Feb 1997 - RRD
!
! USAGE:  CALL PAKINP(RVAR,CVAR,NX,NY,NXY,NX1,NY1,LX,LY,
!              PREC,NEXP,VAR1,KSUM)
!   INPUT ARGUMENT LIST:
!     CVAR  - char      packed char*1 input array of length NXY=NX*NY
!     NX1   - int optional sub-grid left edge in NX units
!     NY1   - int optional sub-grid lower edge in NY units
!     LX    - int first dimension element length of sub-grid
!     LY    - int second dimension element lenght of sub-grid
!     PREC  - real      precision of packed data array
!     NEXP  - int packing scaling exponent
!     VAR1  - real      value of real array at position (1,1)
!   OUTPUT ARGUMENT LIST:
!     RVAR  - real      work array of dimensions LX,LY
!     KSUM  - int rotating checksum of packed data array
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
! $Id: pakinp.f,v 1.1 2009/10/26 15:36:54 jel Exp $
!
!$$$
!dwen(20090830)********************************
!          remove nxy
!    SUBROUTINE PAKINP (SXTNBIT,RVAR,CVAR,NX,NY,NXY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)
    SUBROUTINE PAKINP (SXTNBIT,RVAR,CVAR,NX,NY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)


!dwen(20090831)      IMPLICIT REAL*8 (A-H,O-Z)

      IMPLICIT none

      LOGICAL  , INTENT(IN)    :: SXTNBIT                   ! .TRUE. for 16 bit per number in ARL
      INTEGER,   INTENT(IN)    :: NX,NY                     ! size of input array
      INTEGER,   INTENT(IN)    :: NX1,NY1                   ! optional sub-grid left edge
      INTEGER,   INTENT(IN)    :: LX,LY                     ! length of sub-grid
!dwen(20090830)      CHARACTER, INTENT(IN)    :: CVAR(NXY)                 ! packed input of NX*NY
      CHARACTER, INTENT(IN)    :: CVAR(:)                 ! packed input of NX*NY
!dwen(20090830)      REAL*8,    INTENT(IN)    :: PREC                      ! precision of packed data
      REAL,      INTENT(IN)    :: PREC                      ! precision of packed data
      INTEGER,   INTENT(IN)    :: NEXP                      ! packing scaling exponent
!dwen(20090830)      REAL*8,    INTENT(IN)    :: VAR1                      ! value of array at (1,1)
      REAL,      INTENT(IN)    :: VAR1                      ! value of array at (1,1)

      INTEGER,   INTENT(INOUT) :: KSUM                      ! rotating checksum
!dwen(20090830)      REAL*8,    INTENT(OUT)   :: RVAR(LX,LY)               ! real data unpacked
      REAL,      INTENT(OUT)   :: RVAR(:,:)               ! real data unpacked

!dwen(20090831) *******************
      real                     :: scale,rold,rnew
      integer                  :: j,k,itmp,itmp1,itmp2,jj,i,ii
!*********************************

!---------------------------------------------------------------------------------------------------
!     scaling exponent
! CHG (10/01/03) change from 8 to 16 bit
      IF (.NOT.SXTNBIT) THEN
         SCALE=2.0**(7-NEXP)
      ELSE
         SCALE=2.0**(15-NEXP)
      ENDIF

!     unpack initial value for each row of column 1
      ROLD=VAR1
      DO J=1,NY
!        compute position in 1D packed array at column 1
! CHG (10/01/03) change from 8 to 16 bit, need 2 bytes
         IF(.NOT.SXTNBIT)K=(J-1)*NX+1
         IF(SXTNBIT)K=(J*2-2)*NX+2
!        new value from old value of previous row
! JCL:(4/24/02) to deal with cases in which ICHAR returns negative integers
!         RNEW=((ICHAR(CVAR(K))-127.0)/SCALE)+ROLD
! CHG (10/01/03) change from 8 to 16 bit
         IF(.NOT.SXTNBIT)THEN
           ITMP=ICHAR(CVAR(K))
           IF(ITMP.LT.0)ITMP=ITMP+256
           RNEW=((ITMP-127.0)/SCALE)+ROLD
         ELSE
           ITMP1=ICHAR(CVAR(K-1))
           ITMP2=ICHAR(CVAR(K))
           ITMP=ITMP2*256+ITMP1
           IF(ITMP.LT.0)ITMP=ITMP+65536
           RNEW=((ITMP-32767.0)/SCALE)+ROLD
         END IF
         ROLD=RNEW
!        compute index in output sub-grid
         JJ=J-NY1+1
!        special case of sub-grid starting at left edge
         IF(JJ.GE.1.AND.JJ.LE.LY)RVAR(1,JJ)=RNEW
      END DO

!     only unpack within J-subgrid
      DO J=NY1,(NY1+LY-1)
!        sub-grid array element (1 to LY)
         JJ=J-NY1+1

         ROLD=RVAR(1,JJ)
!        unpack I from element 1 to I sub-grid max
         DO I=2,(NX1+LX-1)
! CHG (10/01/03) change from 8 to 16 bit, need 2 bytes
            IF(.NOT.SXTNBIT)K=(J-1)*NX+I
            IF(SXTNBIT)K=(J*2-2)*NX+I*2
! JCL:(4/24/02) to deal with cases in which ICHAR returns negative integers
!           RNEW=((ICHAR(CVAR(K))-127.0)/SCALE)+ROLD
! CHG (10/01/03) change from 8 to 16 bit
            IF(.NOT.SXTNBIT)THEN
              ITMP=ICHAR(CVAR(K))
              IF(ITMP.LT.0)ITMP=ITMP+256
              RNEW=((ITMP-127.0)/SCALE)+ROLD
            ELSE
              ITMP1=ICHAR(CVAR(K-1))
              ITMP2=ICHAR(CVAR(K))
              ITMP=ITMP2*256+ITMP1
              IF(ITMP.LT.0)ITMP=ITMP+65536
              RNEW=((ITMP-32767.0)/SCALE)+ROLD
            END IF
            ROLD=RNEW
!           sub-grid array element (1 to LX)
            II=I-NX1+1
!           check against precision for true zero
            IF(ABS(RNEW).LT.PREC)RNEW=0.0
            IF(II.GE.1.AND.II.LE.LX)RVAR(II,JJ)=RNEW
         END DO
      END DO

!     only do full-grid checksum when KSUM=0
      IF(KSUM.EQ.0)THEN
!dwen(20090830)         DO K=1,NXY
         DO K=1,nx*ny
! JCL:(4/24/02) to deal with cases in which ICHAR returns negative integers
!            KSUM=KSUM+ICHAR(CVAR(K))
! CHG (10/01/03) change from 8 to 16 bit
            IF(.NOT.SXTNBIT)THEN
              ITMP=ICHAR(CVAR(K))
              IF(ITMP.LT.0)ITMP=ITMP+256
            ELSEIF(MOD(K,2).EQ.0)THEN
              ITMP=ICHAR(CVAR(K))*256+ICHAR(CVAR(K-1))
              IF(ITMP.LT.0)ITMP=ITMP+65536
            END IF
            KSUM=KSUM+ITMP
!           if sum carries over the eighth bit add one
! CHG (10/01/03) change from 8 to 16 bit
            IF(.NOT.SXTNBIT.AND.KSUM.GE.256)KSUM=KSUM-255
            IF(SXTNBIT.AND.KSUM.GE.65536)KSUM=KSUM-65535
         END DO
      END IF

      END SUBROUTINE PAKINP
