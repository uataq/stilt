!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKSET           PACking SETup to initialize routines
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACKING SETUP INITIALIZE PACKING ROUTINES PAKREC AND PAKNDX
!   ONE ENTRY PER DATA FILE TO PACK.  MULTIPLE FILES CAN BE PACKED A
!   SAME TIME (UP TO MGRD - SEE STRUCTURE). A UNIQUE UNIT NUMBER MUST
!   BE DEFINED FOR EACH OUTPUT FILE.  THE FILE NAME CONTAINS THE METEO
!   STRUCTURE INFORMATION. AFTER THIS ROUTINE THE OUTPUT FILE SHOULD
!   OPENED TO THE SPECIFIED UNIT.
!
! PROGRAM HISTORY LOG:
!   Last Revision: 14 Feb 1997 (RRD)
!                  13 Jul 1999 (RRD)
!                  21 Jan 2000 (RRD) - set krec=krec1
!                  13 Feb 2001 (RRD) - fortran90 upgrade
!                  18 Oct 2001 (RRD) - support for large grids
!                  09 Sep 2002 (RRD) - fortran coding standards
!
! USAGE:  CALL PAKSET(LUNIT,FNAME,KREC1,NXP,NYP,NZP)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            METDATA.CFG unless otherwise defined
!   OUTPUT FILES:           none 
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PAKSET(LUNIT,FNAME,KREC1,NXP,NYP,NZP)

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER,       INTENT(IN)    :: lunit     ! output unit number
  CHARACTER(*),  INTENT(INOUT) :: fname     ! file name of METDATA.CFG
  INTEGER,       INTENT(IN)    :: krec1     ! position of index record at time-1
  INTEGER,       INTENT(OUT)   :: nxp, nyp  ! horizontal grid dimensions
  INTEGER,       INTENT(OUT)   :: nzp       ! vertical grid dimension (incl sfc)

!-------------------------------------------------------------------------------

  LOGICAL  :: FTEST ! file existence test
  INTEGER  :: knx,kny,ntot,krec,k,i,nlvl,kg,ng,nv,l,nl,nvar 

! define and save structure between routines
  INCLUDE 'DEFPACK.INC'
  COMMON / PAKCOM / GV, KG

!==>save grid number counter

  SAVE NG
  DATA NG /0/

  NG=NG+1
  IF(NG.GT.MGRD)THEN
     WRITE(*,*)'*ERROR* pakset: '
     WRITE(*,*)' Trying to intialize grid numbr: ', NG
     WRITE(*,*)' Exceeding max numbers of grids: ', MGRD
     STOP 900
  END IF
  KG=NG

!==>check for output format configuration file

  INQUIRE(FILE=FNAME,EXIST=FTEST)
  IF(.NOT.FTEST)THEN
     FNAME='METDATA.CFG'
     INQUIRE(FILE=FNAME,EXIST=FTEST)
     IF(.NOT.FTEST)THEN
        WRITE(*,*)'Unable to find data configuration file: ',FNAME
        WRITE(*,*)'Re-enter the file name ...'
        READ(*,'(A)') FNAME
     END IF
  END IF

  OPEN(LUNIT,FILE=FNAME)
! character field meteo model identification
  READ(LUNIT,'(20X,A4)') GV(NG)%MODEL

! integer grid number and vertical coordinate system
! ksys - 1:sigma  2: pressure  3: z-terrain  4: hybrid
  READ(LUNIT,'(20X,I4)') GV(NG)%IG, GV(NG)%KSYS

! read 12 grid mapping variables
  READ(LUNIT,'(20X,F10.0)')(GV(NG)%GRIDS(I),I=1,12)

! grid size in x,y,z
  READ(LUNIT,'(20X,I4)') GV(NG)%NXG, GV(NG)%NYG, GV(NG)%NLVL

! level heights, variables per level, variable id
  NLVL=GV(NG)%NLVL
  DO NL=1,NLVL
     READ(LUNIT,'(20X,F6.0,I3,99(1X,A4))')                                 &
          GV(NG)%HEIGHT(NL),NVAR,(GV(NG)%VARB(NV,NL),NV=1,NVAR)
          GV(NG)%NVAR(NL)=NVAR
  END DO
  CLOSE (LUNIT)

!==>adjust grid parameters to deal with large grid sizes

  IF(GV(NG)%NXG.GE.1000.OR.GV(NG)%NYG.GE.1000)THEN
!    header record does not support grids of more than 999 
!    therefore in those situations the grid number is
!    converted to character to represent the 1000s digit
!    e.g. @(64)=<1000, A(65)=1000, B(66)=2000, etc

     KNX=GV(NG)%NXG/1000
     KNY=GV(NG)%NYG/1000
     WRITE(GV(NG)%IGC,'(A2)') CHAR(KNX+64)//CHAR(KNY+64)
     GV(NG)%XGPT=.TRUE.
  ELSE 
     WRITE(GV(NG)%IGC,'(I2)') GV(NG)%IG
     GV(NG)%XGPT=.FALSE. 
  END IF

!==>compute extended header length

  GV(NG)%LENH=108
  NLVL=GV(NG)%NLVL
  DO L=1,NLVL
     GV(NG)%LENH=GV(NG)%LENH+8
     NVAR=GV(NG)%NVAR(L)
     DO K=1,NVAR
        GV(NG)%LENH=GV(NG)%LENH+8
     END DO
  END DO

!==>check for multiple extended header records

  GV(NG)%LREC=(GV(NG)%NXG*GV(NG)%NYG)
! check limits
  IF((GV(NG)%LENH).GT.MLEN.OR.(GV(NG)%LREC).LT.108)THEN
     WRITE(*,*)'*ERROR* pakset: Extended header exceeds limits...'
     WRITE(*,*)' Maximum header length : ',MLEN
     WRITE(*,*)' Required HEADER length: ',GV(NG)%LENH
     WRITE(*,*)' Available (nx*xy)     : ',GV(NG)%LREC
     STOP 900
  END IF

! number of header records
  GV(NG)%NHREC=GV(NG)%LENH/GV(NG)%LREC+1
! bytes in last header record
  GV(NG)%NHBYT=GV(NG)%LENH-(GV(NG)%NHREC-1)*GV(NG)%LREC

!==>compute record count offset from index record for each level

  NTOT=GV(NG)%NHREC
  GV(NG)%NREC(1)=NTOT
  NLVL=GV(NG)%NLVL
  DO K=2,NLVL
     NTOT=NTOT+GV(NG)%NVAR(K-1)
     GV(NG)%NREC(K)=NTOT
  END DO

! records per time period
  NTOT=NTOT+GV(NG)%NVAR(NLVL)
  GV(NG)%NRPT=NTOT

! check validity of index record
  KREC=MAX(KREC1,1)
  IF(MOD(KREC-1,NTOT).NE.0)THEN
     WRITE(*,*)'*ERROR* pakset: position record not even multiple'
     WRITE(*,*)'  Record numb in argument:',KREC
     WRITE(*,*)'  Records per time period:',NTOT
     STOP 900
  ELSE
     WRITE(*,*)' NOTICE pakset:'
     WRITE(*,*)' Number of index records = ',GV(NG)%NHREC
     WRITE(*,*)' Number of records /time = ',NTOT
  END IF

!==>set remaining variables

  GV(NG)%MREC=KREC-NTOT     ! internal record counter points to index record
  GV(NG)%NEWT=.TRUE.        ! set flag for first time group for initialization
  GV(NG)%KUNIT=LUNIT        ! unit number

!==>return sizes for variable dimensions

  NXP=GV(NG)%NXG
  NYP=GV(NG)%NYG
  NZP=GV(NG)%NLVL

END SUBROUTINE pakset
