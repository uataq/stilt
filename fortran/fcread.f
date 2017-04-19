SUBROUTINE FCREAD(HANDLE,BUFF,KLEN,KNUM,*)

!-------------------------------------------------------------------------------
! Reading of sequential records using direct access methods
! Last Revised: 06 Mar 2001
!               05 Mar 2003 (RRD) - cleaner EOF exit when buffer exceeds file
!               14 Aug 2008 (RRD) - klen<0 for little-big endian swaps
!               09 Oct 2008 (RRD) - insure klen a number for gfortran
!-------------------------------------------------------------------------------

  LOGICAL            :: FEND
  INTEGER            :: KRET, HANDLE
  INTEGER, PARAMETER :: MBUFLEN = 2048

! master data buffer for up to 80 files (units 10-90)
  CHARACTER BUFF(ABS(KLEN*1)*KNUM)*1, JBUFF(MBUFLEN,90)*1

! buffer pointers for next IO request
  COMMON /KINDEX/ KBYTE(90), KREC(90), KBUFLEN

! buffer pointers for data that have been read
  INTEGER JREC(90)
  DATA JREC/90*0/
  SAVE JREC,JBUFF

! check buffer length
  IF(MBUFLEN.NE.KBUFLEN)THEN
     WRITE(*,*)'ERROR fcread: recompile with buffer - ',KBUFLEN
     STOP
  END IF

! start with file end parameter set to false
  FEND=.FALSE.

! load data for initial request if request record doesn't match the 
! record number of the last record read into the data buffer
  IF(JREC(HANDLE).NE.KREC(HANDLE))THEN
     READ(HANDLE,REC=KREC(HANDLE),IOSTAT=KRET)JBUFF(:,HANDLE)
!    end of file or end of record condition
     IF(KRET.NE.0) FEND=.TRUE.
!    update internal record number to reflect data in buffer 
     JREC(HANDLE)=KREC(HANDLE)
  END IF

! output buffer byte count
  KK=0 

! reset error flag
  KRET=0

! number of variables
  DO I=1,KNUM

!    output buffer for reverse transfer
     LL=KK+ABS(KLEN)+1
   
!    number of bytes per variable
     DO J=1,ABS(KLEN) 

!       transfer data to output buffer
        KK=KK+1
        LL=LL-1
        IF(KLEN.GE.0)THEN
           BUFF(KK)=JBUFF(KBYTE(HANDLE),HANDLE)
        ELSE
           BUFF(LL)=JBUFF(KBYTE(HANDLE),HANDLE)
        END IF

!       check if next byte exceeds input buffer, read records if required
        KBYTE(HANDLE)=KBYTE(HANDLE)+1
        IF(KBYTE(HANDLE).GT.KBUFLEN)THEN
           KREC(HANDLE)=KREC(HANDLE)+1

           IF(.NOT.FEND)                                               &
              READ(HANDLE,REC=KREC(HANDLE),IOSTAT=KRET)JBUFF(:,HANDLE)

!          end of file condition when buffer exceeds file length
           IF(KRET.NE.0)THEN
!             WRITE(*,*)'ERROR: ',KRET,' from fcread on unit - ',HANDLE
!             WRITE(*,*)'Record number ',KREC(HANDLE)
!             WRITE(*,*)'Buffer byte   ',KBYTE(HANDLE) 
!             WRITE(*,*)'Variable      ',I,KNUM
!             WRITE(*,*)'Length        ',J,KLEN
              FEND=.TRUE.
              KRET=0
           END IF

           KBYTE(HANDLE)=1
           JREC=KREC(HANDLE)
        END IF

     END DO
  END DO

  IF(FEND)RETURN 1
  RETURN

END SUBROUTINE fcread
