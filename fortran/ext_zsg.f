! This file contains subroutines for dealing with external specification of the
! internal model levels
! $Id: ext_zsg.f,v 1.4 2009/12/04 20:23:50 trn Exp $
! PROGRAM HISTORY LOG:
!   $Log: ext_zsg.f,v $
!   Revision 1.4  2009/12/04 20:23:50  trn
!   Activate diagnostic message
!
!   Revision 1.3  2009/12/01 18:51:26  trn
!   Added bounds check flag for humboldt, some minor cleanups
!
!   Revision 1.2  2009/11/23 15:17:59  jel
!   single precision for consistency with adviec
!
!   Revision 1.2  2006/05/29 19:19:19  skoerner
!   cosmetics only
!
!   Revision 1.1  2005/12/14 17:05:58  tnehrkor
!   Added support for WRF model output; includes a bug fix for dmass computation in prfcom
!
!   Revision 1.3  2005/10/11 15:54:08  trn
!   Corrected ind_zsg for h < h1
!
!   Revision 1.2  2005/10/07 19:39:35  trn
!   Fixed ind_zsg bug
!
!   Revision 1.1  2005/09/30 18:20:19  trn
!   Intermediate changes for wrf winds/fluxes.  Problems: dmass, stblanl for fluxflg=T
!

subroutine read_zsg(fname,iu,zmdl,zsg,nzm,nlvl,ksfc,sfcl,aa,bb,cc,ret_code)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
! SUBPROGRAM:  read_zsg           read in internal model level heights (m AGL)
!   PRGMMR:    Thomas Nehrkorn  AER, Inc.  September 2005
!
! ABSTRACT:  
!   Provides optional input of internal model level heights (as m AGL) from file
!   In case of I/O errors, the input values are not changed

  USE funits
  implicit none

! Interface:
  CHARACTER(LEN=*), INTENT(IN)    :: fname               !input filename
  INTEGER         , INTENT(IN)    :: iu, nzm             !input unit number, max array dimension for zsg
  INTEGER         , INTENT(INOUT) :: nlvl, ksfc          !number of levels, sfc layer index
  REAL          , INTENT(INOUT) :: zmdl,zsg(nzm),sfcl  !model top height, model z-sigma levels, top of sfc layer
  REAL          , INTENT(INOUT) :: aa,bb,cc            !parameters for default quadratic relationship
  INTEGER, INTENT(OUT) :: ret_code !Return code; =0: ok; <0: nlvl>nzm; >0: other errors

! Automatic arrays:
  REAL             :: hgt(nzm)
! Other local variables
  REAL             :: zmdl_new
  INTEGER            :: ierr, iline, k, nlvl_new, ksfc_new
  CHARACTER(LEN=256) :: reason, line
  LOGICAL            :: ftest

!---------------------------------------------------------------------------------------------------
! Open and read in header info from file
  ret_code=1
  reason = 'File does not exist: '//trim(fname)
  INQUIRE(FILE=fname,EXIST=FTEST)
  if (.not. ftest) goto 900
  ret_code=2
  reason = 'Error opening file: '//trim(fname)
  OPEN(iu,FILE=fname,STATUS='OLD',ACTION='READ',iostat=ierr)
  if (ierr .ne. 0) goto 900

! First line must contain (as list-directed I/O): 
!    the number of levels, the sfc layer index, and the model top height
  ret_code=3
  iline = 1
  write(reason,'(a,i4,2a)') 'Error reading line',iline,' from file: ',trim(fname)
  read(iu,*,iostat=ierr) nlvl_new, ksfc_new, zmdl_new
  if (ierr .ne. 0) goto 900
  if (nlvl_new .gt. nzm) then
     nlvl=nlvl_new
     write(reason,'(a,i4,3a,i4)') 'Specified number of levels:',nlvl_new, &
          & ' in file: ',trim(fname), &
          & ' exceeds max allowable nzm=',nzm
     ret_code=-1
     goto 900
  end if

! Header must be terminated by 'ENDHEADER' line  
  ret_code=4
  line = ' '
  do while (line(1:9) .ne. 'ENDHEADER')
     iline=iline+1
     write(reason,'(a,i4,2a)') 'Error reading header: line',iline,' from file: ',trim(fname)
     read(iu,'(a)',iostat=ierr) line
     if (ierr .ne. 0) goto 900
  end do

! Read in data (One value per line, list-directed I/O)
  ret_code=5
  do k=1,nlvl_new
     iline=iline+1
     write(reason,'(a,i4,2a)') 'Error reading data: line',iline,' from file: ',trim(fname)
     read(iu,*,iostat=ierr) hgt(k)
     if (ierr .ne. 0) goto 900
  end do

  ret_code=6
  if (hgt(1) .le. 0) ierr=1
  do k=2,nlvl_new
     if (ierr .ne. 0) exit
     if (hgt(k) .le. hgt(k-1)) ierr=k
  end do
  if (ierr .ne. 0) then
     write (reason,'(a,i4,a,g15.6)') 'Bad value for height specified for level ',ierr,' : ',hgt(ierr)
     goto 900
  endif
  ret_code=0

! Assign new values to output
! disable default quadratic relationship
  aa = -HUGE(aa)
  bb = aa
  cc = aa

! Set sfc layer values:
  ksfc = ksfc_new
  sfcl = hgt(ksfc)
! First set new model top:
  nlvl = nlvl_new
  zmdl = zmdl_new
! Compute sigma from zmdl:
  zsg(1:nlvl) = 1.0 - hgt(1:nlvl)/zmdl

! Print out results:
  write (*,'(2a)') 'Read_zsg: Model levels as determined from file:',trim(fname)
  goto 950

! error exit:
900 continue
  write (*,'(a/1x,a)') 'Read_zsg: zsg model levels not changed because optional reading from file failed:',&
       & trim(reason)
  if (ret_code .ge. 0 .and. nlvl .le. nzm) hgt(1:nlvl) = zmdl*(1.0-zsg(1:nlvl))

! generate this output in either case:
950 continue
  close(unit=iu,iostat=ierr)
  if (ierr .ne. 0) then
     write (*,*) 'read_zsg warning: error closing file: '//trim(fname)
  end if

  if (ret_code .ge. 0 .and. nlvl .le. nzm) then
     write (*,'(1x,a,i4)') 'Number of levels: ',nlvl,'Sfc layer index ksfc: ',ksfc
     write (*,'(1x,a,g15.6)') 'Model top zmdl: ',zmdl,'Sfc layer top: ',sfcl
     write (*,'(1x,a8,2a15)') ' Level ',' Height (m AGL) ',' z-sigma (ZSG) '
     write (*,'(1x,i8,2g15.6)') (k,hgt(k),zsg(k),k=1,nlvl)
  end if

  return
end subroutine read_zsg

subroutine ind_zsg(zmdl,zsg,nlvl,zt,zx,aa,bb,cc)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
! SUBPROGRAM:  ind_zsg           compute vertical level index (as in adviec)
!   PRGMMR:    Thomas Nehrkorn  AER, Inc.  September 2005
!
! ABSTRACT:  
!   From the input point sigma value, and the model sigma levels, compute
!   the vertical level index.  This routine uses a piecewise linear
!   function, unless the values for aa,bb,cc indicate that
!   the default quadratic relationship is to be used (in this case,
!   the quadratic equation (see runset) is inverted)
!   

  implicit none

  integer, intent(in) :: nlvl !number of levels
  real, intent(in) :: zmdl, zsg(nlvl), zt !model top, model level and point sigma
  real, intent(in) :: aa,bb,cc !parameters for default quadratic relationship
  real, intent(out) :: zx !vertical level index (can be < 1, or > nlvl)

! AGL heights (as in adviec)
  real :: hgt(nlvl), zz
  integer :: k,level

  zz = ZMDL*(1.0-MIN(1.0,ZT))
  IF (AA < -9999e0) THEN
     IF (BB >= -9999e0 .OR. CC >= -9999e0) then
        write (*,*) 'In ind_zsg, either all of AA, BB, CC need to be set, or none. Stop.'
        STOP 'In ind_zsg, either all of AA, BB, CC need to be set, or none. Stop.'
     end IF
     
     hgt(1:nlvl) = zmdl*(1.0-zsg(1:nlvl))
     level = nlvl !treat heights above the top level the same as those in the top layer
     do k=1,nlvl
        if (zz .lt. hgt(k)) then
           level = k
           exit
        end if
     end do

     if (level .eq. 1) then
        zx = 0.5 + 0.5*zz/hgt(1) !below lowest level
     else
        zx = float(level-1) + (zz - hgt(level-1)) / (hgt(level)-hgt(level-1))
     end if

  else !aa,bb,cc are set, use the quadratic relationship:
     ZX=(-BB+SQRT(BB*BB-4.0*AA*(CC-ZZ)))/(2.0*AA)
  endif

  return
end subroutine ind_zsg
