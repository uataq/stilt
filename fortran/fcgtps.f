INTEGER FUNCTION FCGTPS(HANDLE)

!-------------------------------------------------------------------------------
! Returns the current position within the file specified by HANDLE, i.e. how   
! many bytes have been read since the beginning of the file. 
! Last Revised: 06 Mar 2001
!-------------------------------------------------------------------------------

  INTEGER HANDLE
  COMMON /KINDEX/ KBYTE(90), KREC(90), KBUFLEN

  FCGTPS = (KREC(HANDLE)-1)*KBUFLEN+KBYTE(HANDLE)-1

END FUNCTION fcgtps
