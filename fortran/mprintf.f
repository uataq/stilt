module module_debug
!-------------------------------------------------------------------------------
!
!  Purpose
!
!
!-------------------------------------------------------------------------------
!  History
!  -------
!
!  3/26/2008  Initial version. Thomas Nehrkorn
!
!  $Id: mprintf.f,v 1.1 2009/11/23 15:46:36 jel Exp $
!-------------------------------------------------------------------------------

   integer, parameter :: QUIET=-100, LOGFILE=-2, DEBUG=0, INFORM=1, WARN=2, ERROR=3, STDOUT=100

   integer :: the_debug_level = DEBUG

   logical :: have_set_logname = .false.

   contains

   subroutine set_debug_level(ilev)

      implicit none

      ! Arguments
      integer, intent(in) :: ilev

      the_debug_level = ilev

   end subroutine set_debug_level


   subroutine mprintf(assertion, level, message)

      implicit none

      ! Arguments
      integer, intent(in) :: level
      logical, intent(in) :: assertion
      character (len=*), intent(in) :: message

      character (len=1024) :: ctemp

      if (assertion) then
         ! If this is a debug message give up if level is not high enough
         if (level == DEBUG .and. the_debug_level > DEBUG) return

         if (level == DEBUG) then
            write(ctemp,'(a)') 'DEBUG: '
         else if (level == INFORM) then
            write(ctemp,'(a)') 'INFORM: '
         else if (level == WARN) then
            write(ctemp,'(a)') 'WARNING: '
         else if (level == ERROR) then
            write(ctemp,'(a)') 'ERROR: '
         end if

         write (*,*) trim(ctemp) // ' ' // trim(message)

         if (level == ERROR) stop 'abort'
      end if
   end subroutine mprintf

end module module_debug
