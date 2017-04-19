SUBROUTINE FCCLOS(HANDLE,*)

!-------------------------------------------------------------------------------
! Clean exit for all open units
! Last Revised: 06 Mar 2001
!-------------------------------------------------------------------------------

  INTEGER HANDLE

  CLOSE(HANDLE,ERR=900)
  RETURN

  900 RETURN 1

END SUBROUTINE fcclos
