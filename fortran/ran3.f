!#######################################################################
      REAL FUNCTION RAN3(IDUMI)
!
!  RETURNS A UNIFORM RANDOM DEVIATE BETWEEN 0.0 AND 1.0.  SET IDUM TO
!  ANY NEGATIVE VALUE TO INITIALIZE OR REINITIALIZE THE SEQUENCE.
!  THIS FUNCTION IS TAKEN FROM W.H. PRESS', "NUMERICAL RECIPES" P. 199.
!
!
!  $Id: ran3.f,v 1.1 2009/10/26 15:36:55 jel Exp $
!
      IMPLICIT REAL (A-H,O-Z)

      INTEGER, INTENT(IN) :: IDUMI
!      SAVE
!      IMPLICIT REAL*4 (M)
      SAVE INEXT,INEXTP
      SAVE MA
!      PARAMETER (MBIG = 4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      PARAMETER (MBIG = 1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!
!  ACCORDING TO KNUTH, ANY LARGE MBIG, AND ANY SMALLER (BUT STILL LARGE)
!  MSEED CAN BE SUBSTITUTED FOR THE ABOVE VALUES.
      DIMENSION MA(55)
      DATA IFF /0/

!---------------------------------------------------------------------------------------------------
       IDUM = IDUMI
   10  IF (IDUM < 0 .OR. IFF == 0) THEN
         IFF = 1
         MJ = MSEED-IABS(IDUM)
!         MJ = MSEED-DBLE(IABS(IDUM))
         MJ = MOD(MJ,MBIG)
!         MJ = DMOD(MJ,MBIG)
         MA(55) = MJ
         MK = 1
         DO 11 I=1,54
            II = MOD(21*I,55)
            MA(II) = MK
            MK = MJ-MK
            IF (MK < MZ) MK = MK+MBIG
            MJ = MA(II)
   11    CONTINUE
         DO 13 K=1,4
            DO 12 I = 1,55
               MA(I) = MA(I)-MA(1+MOD(I+30,55))
               IF (MA(I) < MZ) MA(I) = MA(I)+MBIG
   12       CONTINUE
   13    CONTINUE
         INEXT = 0
         INEXTP = 31
         IDUM = 1
      ENDIF
      INEXT = INEXT+1
      IF (INEXT == 56) INEXT = 1
      INEXTP = INEXTP+1
      IF (INEXTP == 56) INEXTP = 1
      MJ = MA(INEXT)-MA(INEXTP)
      IF (MJ < MZ) MJ = MJ+MBIG
      MA(INEXT) = MJ
      RAN3 = MJ*FAC
! JCL:FOR SOME REASON, THE ALGORITHM RETURNS NEGATIVE NUMBERS A LOT OF
!     TIMES--NEED TO REDO IF THIS HAPPENS
      IF (RAN3 < 0) GOTO 10

      END FUNCTION RAN3
