!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DALLOC           DeALLOCate meteorological array 
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:06-05-12
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DEALLOCATES ALL THE METEOROLOGOCAL VARIABLES DEFINED IN ADVPNT  
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 May 2006 (RRD) - initial version
!
! USAGE:  CALL DALLOC 
!     
!   INPUT ARGUMENT LIST:    none     
!   OUTPUT ARGUMENT LIST:   none      
!   INPUT FILES:            none 
!   OUTPUT FILES:           none 
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE DALLOC 

  USE metval

  CHARACTER(80)            :: ecode
  integer :: kret
  ECODE=''
  DEALLOCATE( u,v,w,a,t,q,e,p,x,h,stat=kret)
  if (kret .ne. 0) ECODE=ECODE//' 3d-a '
  DEALLOCATE( p0,rt,u0,v0,t0,zi,stat=kret)
  if (kret .ne. 0) ECODE=ECODE//' 2d-a '
  DEALLOCATE( ds,ss,uf,vf,hf,sf,stat=kret)
  if (kret .ne. 0) ECODE=ECODE//' flux '
  DEALLOCATE( gx,gy,z0,zt,lu,stat=kret)
  if (kret .ne. 0) ECODE=ECODE//' fixed '

  DEALLOCATE(d,tlrams,sigwrams,cfxup1,cfxup2,cfxdn1,&
             dfxup1,dfxup2,dfxdn1,efxup1,efxup2,efxdn1,&
             tke,tl,sigw,dmass,xm,hm,stat=kret)
  if (kret .ne. 0) ECODE=ECODE//' 3d-b '
  DEALLOCATE(rc,tc,sw,lc,sm,w0,muu,muv, &
              mu,msfu,msfv,msft,lf,zloc,stat=kret)
  if (kret .ne. 0) ECODE=ECODE//' 2d-b '

  if (ECODE .ne. '') write (*,'(2a)') &
       & 'dalloc warning: error in final deallocation of the following variable(s): ', &
       & ecode

END SUBROUTINE dalloc 
