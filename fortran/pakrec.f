!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKREC           PAcK a RECord writes one meteo record
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACK A RECORD WRITES A SPECIFIC RECORD OF DATA FIELD TO UNIT KUN
!   PREVIOUSLY OPENED FOR DIRECT UNFORMATTED ACCESS
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 16 Oct 1997 (RRD) 
!                 02 Feb 2001 (RRD) - fortran90 upgrade
!                 18 Oct 2001 (RRD) - support large grid domains
!                 11 Apr 2002 (RRD) - intent cvar to inout
!
! USAGE:  CALL PAKREC(LUNIT,RVAR,CVAR,NX,NY,NXY,KVAR,IY,IM,ID,IH,MN,IC,LL,KINI)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           defined by LUNIT
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PAKREC(LUNIT,RVAR,CVAR,NX,NY,NXY,KVAR,IY,IM,ID,IH,MN,IC,LL,KINI)

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER,      INTENT(IN)    :: LUNIT       ! output unit number
  INTEGER,      INTENT(IN)    :: NX,NY       ! dimensions of RVAR
  INTEGER,      INTENT(IN)    :: NXY         ! dimensions of CVAR
  REAL,         INTENT(IN)    :: RVAR(NX,NY) ! input data to be packed
  CHARACTER(1), INTENT(INOUT) :: CVAR(NXY)   ! packed data array
  CHARACTER(4), INTENT(IN)    :: KVAR        ! descriptor of variable written
  INTEGER,      INTENT(IN)    :: IY,IM,ID    ! date identification
  INTEGER,      INTENT(IN)    :: IH,MN       ! time identification (MN-minutes)
  INTEGER,      INTENT(IN)    :: IC          ! forecast hour, ICX hour for >99
  INTEGER,      INTENT(IN)    :: LL          ! level indicator 
  INTEGER,      INTENT(IN)    :: KINI        ! initialization (0-no 1-yes)

!-------------------------------------------------------------------------------

  CHARACTER(50) :: LABEL                  ! standard record label
  REAL          :: prec,var1 
  INTEGER       :: nvar,nxyg,nexp,jrec,ksum 
  INTEGER       :: k,kk,ng,kg,nv,il,icw 

! pass structure to other routines
  INCLUDE 'DEFPACK.INC'
  COMMON / PAKCOM / GV, NG

!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE PAKOUT(RVAR,CVAR,NX,NY,NXY,PREC,NEXP,VAR1,KSUM)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)  :: nx,ny,nxy   ! dimension limits
  REAL,      INTENT(IN)  :: rvar(nx,ny) ! data array to be packed  
  CHARACTER, INTENT(OUT) :: cvar(nxy)   ! packed char*1 output array
  REAL,      INTENT(OUT) :: prec        ! precision of packed data array
  INTEGER,   INTENT(OUT) :: nexp        ! packing scaling exponent
  REAL,      INTENT(OUT) :: var1        ! value of real array at position (1,1)
  INTEGER,   INTENT(OUT) :: ksum        ! rotating checksum of packed data 
  END SUBROUTINE pakout
  END INTERFACE

!-------------------------------------------------------------------------------

!==>check if requested unit number properly opened

  KG=0
  DO KK=1,NG
     IF(LUNIT.EQ.GV(KK)%KUNIT)KG=KK
  END DO
  IF(KG.EQ.0)THEN
     WRITE(*,*)'*ERROR* pakrec: Requesting uninitialized unit'
     WRITE(*,*)' Require initial call to PAKSET for unit: ',LUNIT
     STOP 900
  END IF

!==>test grid dimensions for consistency with initialization

  NXYG=GV(KG)%NXG*GV(KG)%NYG
  IF(GV(KG)%NXG.NE.NX.OR.GV(KG)%NYG.NE.NY.OR.NXYG.NE.NXY)THEN
     WRITE(*,*)'*ERROR* pakrec: file dimensions do not match'
     WRITE(*,*)'Initialization: ',GV(KG)%NXG,GV(KG)%NYG, NXYG
     WRITE(*,*)'Argument list : ',NX,NY,NXY
     STOP 900
  END IF

!==>standard forecast hour to write cannot exceed two digits

  ICW=MIN(IC,99)

!==>set all base variables with first entry at each new time

  IF(GV(KG)%NEWT)THEN

!    increment internal counter to next index record
     GV(KG)%MREC=GV(KG)%MREC+GV(KG)%NRPT

!    extended forecast hour (3 digit)
     GV(KG)%ICX=IC

!    save initial times for headers
     GV(KG)%IY0=IY
     GV(KG)%IM0=IM
     GV(KG)%ID0=ID
     GV(KG)%IH0=IH
     GV(KG)%MN0=MN
     GV(KG)%IC0=ICW

!    set switch, reset by pakndx
     GV(KG)%NEWT=.FALSE.

!    initialize all records in this time group to NULL
     IF(KINI.EQ.1)CALL PAKINI(KG,CVAR,NXY)

  ELSE

!==>check current if current record consistent with first

     IF(IY.NE.GV(KG)%IY0.OR.IM.NE.GV(KG)%IM0.OR.                           &
        ID.NE.GV(KG)%ID0.OR.IH.NE.GV(KG)%IH0)THEN

        WRITE(*,*)'*ERROR* pakrec - at index: ',GV(KG)%MREC
        WRITE(*,*)'  Argument list times    : ',IY,IM,ID,IH
        WRITE(*,*)'  Do not match initial   : ',GV(KG)%IY0,                  &
                   GV(KG)%IM0, GV(KG)%ID0, GV(KG)%IH0
        STOP 900

     END IF
  END IF

!==>when no data is supplied just do a normal return
!   normally used in conjunction with the initialization flag

  IF(KVAR.EQ.'NULL')RETURN

!==>check vertical index

  IF(LL.LT.1.OR.LL.GT.GV(KG)%NLVL)THEN
     WRITE(*,*)'*ERROR* pakrec  : Level indicator out of range'
     WRITE(*,*)'  Argument level: ',LL
     WRITE(*,*)'  Valid Range   : 1 to ',GV(KG)%NLVL
     STOP 900
  ELSE
!    level indicator should =0 at the surface
     IL=LL-1
  END IF

!==>compute the record offset based upon variable match

  NV=0
  NVAR=GV(KG)%NVAR(LL)
  DO K=1,NVAR
     IF(KVAR.EQ.GV(KG)%VARB(K,LL))NV=K
  END DO

  IF(NV.EQ.0)THEN
     WRITE(*,*)'*ERROR* pakrec: Variable not in CFG file'
     WRITE(*,*)' Argument list variable: ',KVAR
     WRITE(*,*)' At level: ',LL  
     WRITE(*,*)' File list: ',GV(KG)%VARB(:,LL)
     STOP 900
  END IF

!==>pack data and write

! convert real to packed character
  CALL PAKOUT(RVAR,CVAR,NX,NY,NXY,PREC,NEXP,VAR1,KSUM)

! save checksum in table
  GV(KG)%CHKS(NV,LL)=KSUM

! write index portion of record
  IF(GV(KG)%XGPT)THEN
     WRITE(LABEL,'(6I2,A2,A4,I4,2E14.7)')            &
     IY,IM,ID,IH,ICW,IL,GV(KG)%IGC,KVAR,NEXP,PREC,VAR1
  ELSE
     WRITE(LABEL,'(7I2,A4,I4,2E14.7)')              &
     IY,IM,ID,IH,ICW,IL,GV(KG)%IG,KVAR,NEXP,PREC,VAR1
  END IF

! compute record based upon variable offset
  JREC=GV(KG)%MREC+GV(KG)%NREC(LL)+NV-1
  IF(JREC.LE.1)THEN
     WRITE(*,*)'*ERROR* pakrec: output record <=1'
     WRITE(*,*)'  Index record: ',GV(KG)%MREC
     WRITE(*,*)'  Level offset: ',GV(KG)%NREC(LL)
     WRITE(*,*)'  Varbl offset: ',NV
     STOP 900
  ELSE
     WRITE(GV(KG)%KUNIT,REC=JREC)LABEL,CVAR
  END IF

END SUBROUTINE pakrec
