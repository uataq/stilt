SUBROUTINE FCPTPS(HANDLE,NBYTE,*)

!-------------------------------------------------------------------------------
! Sets byte and record pointers within the data buffer.  The value of NBYTE is 
! assumed to be the last byte read, therefore the next read will be at NBYTE+1.
! Last Revised: 06 Mar 2001
!-------------------------------------------------------------------------------

  INTEGER HANDLE
  COMMON /KINDEX/ KBYTE(90), KREC(90), KBUFLEN

! compute record number assuming buffer length
  KREC(HANDLE)=INT(NBYTE/KBUFLEN)+1
  IF(KREC(HANDLE).LT.1)THEN
     WRITE(*,*)'ERROR fcptps: position to record - ',KREC(HANDLE)
     RETURN 1
  END IF

! compute next byte position to read within buffer
  KBYTE(HANDLE)=MOD(NBYTE,KBUFLEN)+1

END SUBROUTINE fcptps
