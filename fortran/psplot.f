!--------------------------------------------------------------------
! PSPLOT - Postscript graphics library
! Contains the source code of the PSPLOT graphics library.  The 
! plotting routines use the subroutines in this library.  The 
! source code is available from www.nova.edu/ocean/psplot.html
! Note that some of the routines have been changed from the orignal 
! source code available from www.nova.edu to support GIS 
! conversions and Fortran 90 compilation.
!********************************************************************
! ARL Addition
! 26 Sep 2002 (RRD) - changed format on contour value for GIS
! 26 Oct 2004 (RJP) - added dbf capability
! 20 Dec 2004 (RRD) - g95 compatibility
! 13 Oct 2005 (RRD) - expanded GIS capability (IGIS = 1,2)
! 03 Jan 2008 (GDR) - end of loop exits to fix complicated contours
!                     causing text dumps (denoted gdr20080102)
! 29 Apr 2008 (RRD) - increased ntotpt=20000
! 09 Oct 2008 (RRD) - recoded obsolete fortran statements
!********************************************************************
!*****GISDUMP
      subroutine gisdump(x,y,npts,cv)

      dimension x(npts),y(npts)

!rrd  write(16,'(E7.1,I6)')cv,npts
      write(16,'(E14.7,I6)')cv,npts

      do 10 nn=1,npts
        write(16,'(F10.7,F10.7)')x(nn),y(nn)
   10 continue
      return
      end
!*****ARC
      subroutine arc(xc,yc,rad,ang1,ang2)
      character*132 cmdstr
      common/plt1/cmdstr
      common/cnvcom/conver
      radi=rad*conver
      xci=xc*conver
      yci=yc*conver
      cmdstr=' '
      write(cmdstr,                                                            &
      '(f8.2,'' '',f8.2,'' '',f8.2,'' '',2f8.2,'' arcit'')') xci,yci,          &
      radi,ang1,ang2
      call filler
      return
      end
!*****ARCTO
      subroutine arcto(x1,y1,x2,y2,rad)

      character*132 cmdstr
      character*80 ifrmt
      common/plt1/cmdstr
      common/cnvcom/conver

      radi=rad*conver
      x1i=x1*conver
      y1i=y1*conver
      x2i=x2*conver
      y2i=y2*conver
      cmdstr=' '
      ifrmt= '(f8.2,'' '',f8.2,'' '',f8.2,'' '',f8.2,'' '',f8.2,'//            &
       ''' arcto 4 {pop} repeat'')'
      write(cmdstr,ifrmt)x1i,y1i,x2i,y2i,radi
      call filler
      return
      end
!*****AROHED
      subroutine arohed(xpp,ypp,dir,arolnp,sprang,locxy)
!  xpp,ypp are coordinates of tip of arrowhead
!  dir is angle, in degrees, of arrowhead, measured east from north
!  arolnp is the length, in inches, of the sides of the arrowhead
!  sprang is the angle, in degrees, from one side of the arrowhead
!  to the arrow, that is half the angular spread of the arrowhead
!  locxy=0 x,y at arrow point
!  locxy=1 x,y at center of arrow
!  locxy=2 x,y,at tail of arrow
      parameter (rdpdeg=3.14159265358979/180.)
      dimension ymul(3),x(3),y(3)
      equivalence(x(1),x1),(x(2),x2),(x(3),x3)
      equivalence(y(1),y1),(y(2),y2),(y(3),y3)
      character*132 cmdstr
      common/plt1/cmdstr
      common/cnvcom/conver
      data ymul/0.,.5,1./
!      data rdpdeg/.01745329/

      xp=xpp
      yp=ypp
      ya=dir*rdpdeg
      cosdir=cos(ya)
      sindir=sin(ya)
      arolen=arolnp
      sprrad=sprang*rdpdeg
      xa=sin(sprrad)*arolen
      ya=cos(sprrad)*arolen
      i=locxy
      x1=0.
      y1=ya*ymul(i+1)
      x2=xa
      y2=y1-ya
      x3=-xa
      y3=y2

!  Now rotate
      do i=1,3
        xa=x(i)
        ya=y(i)
        x(i)=((xa*cosdir)+(ya*sindir))+xp
        y(i)=((ya*cosdir)-(xa*sindir))+yp
      enddo

      do n=1,3
        x(n)=x(n)*conver
        y(n)=y(n)*conver
      enddo

      write(cmdstr,'(6f8.2,'' Ah'')')x3,y3,x1,y1,x2,y2
      call filler
      return
      end
!*****ARROW
      subroutine arrow(xss,yss,xpp,ypp,arolnp,sprang,locxy)
!  xpp,ypp are coordinates of tip of arrowhead
!  dir is angle, in degrees, of arrowhead, measured east from north
!  arolnp is the length, in inches, of the sides of the arrowhead
!  sprang is the angle, in degrees, from one side of the arrowhead
!  to the arrow, that is half the angular spread of the arrowhead
!  locxy=0 x,y at arrow point
!  locxy=1 x,y at center of arrow
!  locxy=2 x,y,at tail of arrow
      parameter (rdpdeg=3.14159265358979/180.)
      dimension ymul(3),x(3),y(3)
      equivalence(x(1),x1),(x(2),x2),(x(3),x3)
      equivalence(y(1),y1),(y(2),y2),(y(3),y3)
      character*132 cmdstr
      common/plt1/cmdstr
      common/cnvcom/conver
      data ymul/0.,.5,1./
!      data rdpdeg/.01745329/

      pi=4.*abs(atan(1.))
      dv=ypp-yss
      du=xpp-xss
      dir=90.-atan2(dv,du)*180./pi
      xp=xpp
      yp=ypp
      ya=dir*rdpdeg
      cosdir=cos(ya)
      sindir=sin(ya)
      arolen=arolnp
      sprrad=sprang*rdpdeg
      xa=sin(sprrad)*arolen
      ya=cos(sprrad)*arolen
      i=locxy
      x1=0.
      y1=ya*ymul(i+1)
      x2=xa
      y2=y1-ya
      x3=-xa
      y3=y2

!  Now rotate
      do i=1,3
        xa=x(i)
        ya=y(i)
        x(i)=((xa*cosdir)+(ya*sindir))+xp
        y(i)=((ya*cosdir)-(xa*sindir))+yp
      enddo

      do n=1,3
        x(n)=x(n)*conver
        y(n)=y(n)*conver
      enddo

      xssi=xss*conver
      yssi=yss*conver
      write(cmdstr,'(8f8.2,'' Ar'')')x3,y3,x2,y2,x1,y1,xssi,yssi
      call filler
      return
      end
!*****AXIS
      subroutine axis (x,y,iscr,nc,size,theta,ymin,dy)
!  postscript version
!
!      this routine draws a labled axis and produces a tic mark at
!      every inch with the value of the coordinate plotted above
!      (or below) it.
!  inputs:
!      x            x axis starting position
!      y            y axis starting position
!      iscr            axis title (h)
!      nc            number of characters in title
!      size      length of axis in in/cm (fp)
!      theta      angle to project axis (fp)
!      ymin      annotation starting value (fp)
!      dy      scaling increment between tic marks (fp)
!
      character*80 scr
      character * (*) iscr
      parameter (rdpdeg=3.14159265358979/180.)
!      dimension iscr(20),nscr(20)
!      data rdpdeg/.01745329/

      siztic=.1
      siznum=.13
      sizttl=.15
      rnc=nc+.1
      sig=sign(1.,rnc)
      nac=iabs(nc)
      th=theta*rdpdeg
      n = size + 0.50
      cth = cos  (th)
      sth = sin  (th)
      tn = n
      xb = x
      yb = y
      xa = x - siztic * sig * sth
      ya = y + siztic * sig * cth
      call plot (xa,ya,3)
      do i =1,n
        call plot (xb,yb,2)
        xc = xb + cth
        yc = yb + sth
        call plot (xc,yc,2)
        xa = xa + cth
        ya = ya + sth
        call plot (xa,ya,2)
        xb = xc
        yb = yc
      enddo
      ady = dy
      absv = ymin + ady * tn
      exp = 0.0
!# modified: 10/9/2008
      if ( ady.NE.0.0 ) then      
         do while (ady.GE.100.0)
            ady = ady / 10.0
            absv = absv / 10.0
            exp = exp + 1.0
         end do
         do while (ady.LT.0.01)
            ady = ady * 10.0
            absv = absv * 10.0
            exp = exp - 1.0
         end do
      end if
!# modified: 10/9/2008
      if(sig.ge.0) then
        dist = .15
      else
        dist=-(.15+siznum)
      endif
      xa = xb - dist * sth
      ya = yb + dist * cth
      n = n + 1
      do i = 1,n
        if(i.ne.n) then
          call kekflt (xa,ya,siznum,absv,theta,2,1)
        else
          call kekflt (xa,ya,siznum,absv,theta,2,0)
        endif
        absv = absv - ady
        xa = xa - cth
        ya = ya - sth
      enddo
      if(abs(exp).lt.1.e-5) then
        tnc=nac
      else
        tnc = nac + 7
      endif
      if(sig.ge.0)then
        dist = .3+siznum
      else
        dist=-(.3+siznum+sizttl)
      endif

      xa = x + (size / 2.0)*cth - dist* sth
      ya = y + (size / 2.0)*sth + dist* cth
      call keksym (xa,ya,sizttl,iscr,theta,nac,1)

      if(abs(exp).gt.1.e-5) then
      call keksym (999.,999.,.8*sizttl,' (X10',theta,5,0)
      write(scr,'(i10)')int(exp)
      call blkstp(scr,80,scr,lsc)
!      read(scr,'(20a4)')nscr
      call super(scr,lsc,.8*sizttl,theta)
      call keksym(999.,999.,.8*sizttl,')',theta,1,0)
      endif
   45 return
      end
!*****BLKSTP
      subroutine blkstp(ch,ndim,a,leng)
!rrd  character*1 ch(ndim),a(ndim)
      character*80 ch,a
!  Strip out blanks only (leave in esc, etc.)
      i=1
      leng=0
   10 continue
!rrd  if(ichar(ch(i)).ne.32)then
      if(ichar(ch(i:i)).ne.32)then
        leng=leng+1
!rrd    a(leng)=ch(i)
        a(leng:leng)=ch(i:i)
      endif

      if(i.eq.ndim) then
!       Blankfill remainder of output array
        do l=leng+1,ndim
!rrd      a(l)=' '
          a(l:l)=' '
        enddo
        return
      endif

      i=i+1
      goto 10
      end
!*****BORDER
      subroutine border(xlen,ylen,itic,ibord,majrx,minrx,majry,minry)
!  draws a rectangular border with tic marks
!  xlen,ylen are x,y axes lengths in inches
!  itic is a 4-digit integer indicating which borders will have tic marks(1=tic)
!  if itic<0, tic marks are drawn to outside of axes
!  ibord is 4-digit integer indicating which borders are to be drawn (digit=1)
!  or omitted (digit=0)
!  the following refers to both itic and ibord:
!    first  digit refers to left vertical axis
!    second digit refers to bottom horizontal axis
!    third  digit refers to right vertical axis
!    fourth digit refers to top horizontal axis
!  majrx and majry are number of major divisions for x and y axes
!  minrx and minry are number of minor divisions per major division
!  program does not alter current plot origin
      logical ticit,bordit
      parameter (pi=3.14159265358979,pi2=pi/2.)
!      pi2=2.*abs(atan(1.)) ! no need to use transcendental function to find pi
      jtic=iabs(itic)
      jbord=iabs(ibord)
      itsign=isign(1,itic)
      tmaj=.1
      tmin=tmaj/2.

      do k=1,4
        it=jtic/(10**(4-k))
        ib=jbord/(10**(4-k))
        if(mod(it,2).ne.0) then
          ticit=.true.
        else
          ticit=.false.
        endif
        if(mod(ib,2).ne.0) then
          bordit=.true.
        else
          bordit=.false.
        endif

        if((ticit.or.bordit)) then
          if(mod(k,2).eq.1) then
            xs=(xlen/2.)*(k-1.)
            ys=(ylen/2.)*(3.-k)
            numdiv=majry*minry
            minor=minry
            xdiv=0.
            ydiv=(k-2.)*ylen/float(numdiv)
          else
            xs=(xlen/2)*(k-2)
            ys=(ylen/2)*(k-2)
            numdiv=majrx*minrx
            minor=minrx
            xdiv=(3.-k)*xlen/float(numdiv)
            ydiv=0.
          endif

          if(bordit) then
            call plot(xs,ys,3)
            call plot(xs+numdiv*xdiv,ys+numdiv*ydiv,2)
          endif

          if(ticit) then
            ang=(k-1)*pi2
            do n=0,numdiv
              if(mod(n,minor).eq.0) then
              tlen=itsign*tmaj
            else
              tlen=itsign*tmin
            endif
            x1=xs+n*xdiv
            y1=ys+n*ydiv
            x2=x1+tlen*cos(ang)
            y2=y1+tlen*sin(ang)
            call sldlin(x1,y1,x2,y2,0.)
          enddo
        endif
      endif

      enddo

      return
      end
!*****CHOPIT
      subroutine chopit(xorg,yorg)
!  This routine logically closes the present page and opens another one.
!  xorg,yorg are the coordinates of the origin for the first plot
      logical prtrt
      character*80 fileout,scr
      character*132 cmdstr

      common/plt1/cmdstr
      common/plt2/fac
      common/chpcom/ientry,prtrt
      common/pagcom/npage
      common/io/fileout,inew

      if(.not.prtrt) then
        cmdstr='0 8.5 inch translate -90 rotate'
        call filler
      endif

      cmdstr='%----------------chopit'
      call filler
      cmdstr='stroke showpage'
      call filler
      npage=npage+1
      write(scr,'(i4)')npage
      call blkstp(scr,80,scr,nch)
      cmdstr='%%Page: '//scr(1:nch)//' '//scr(1:nch)
      call filler
!rrd-start
      cmdstr='save'
      call filler
!rrd-stop
      cmdstr='newpath 0 0 moveto'
      call filler

      inew=inew+1
      ientry=999        ! Set so psinich will use factor
      facsav=fac
!  Assume user wants initial origin coordinates independent of factor
      call psinich(prtrt)
      call plot(xorg,yorg,-3)
      call factor(facsav)
      return
      end
!*****CIRCLE
      subroutine circle(xc,yc,rad,fill)
      character*132 cmdstr,scrc
      common/plt1/cmdstr
      common/cnvcom/conver
      logical fill
      xci=xc*conver
      yci=yc*conver
      radi=rad*conver
      scrc=' '
      write(scrc,'(f8.2,'' '',f8.2,'' '',f8.2,'' C'')') xci,yci,               &
      radi
      if(fill) then
        cmdstr='Np '//scrc(1:lenstr(scrc,132))//' fill'
      else
        cmdstr='Np '//scrc(1:lenstr(scrc,132))//' stroke'
      endif
      call filler
      return
      end
!*****CLGEN
      subroutine clgen (z,mx,nx,nny,cclo,chi,cinc,nla,nlm,cl,ncl,icnst)
      dimension cl(*),z(mx,nny)
      common/conre1/ioffp,spval
!
! clgen puts the values of the contour levels in cl.
! variable names match those in conrec, with the following additions.
!         ncl     -number of contour levels put in cl.
!         icnst   -flag to tell conrec if a constant field was detected.
!                 .icnst=0 means non-constant field.
!                 .icnst non-zero means constant field.
!
! to produce non-uniform contour level spacing, replace the code in this
! subroutine with code to produce whatever spacing is desired.
!

      icnst = 0
      ny = nny
      clo = cclo
      glo = clo
      ha = chi
      fanc = cinc
      crat = nla
!# modified: 10/9/2008
      if (ha-glo.LT.0.0) then
         glo = ha
         ha = clo
         go to 45
      elseif (ha-glo.GT.0.0) then
         go to 45
      end if
!# modified: 10/9/2008
      if (glo .ne. 0.) go to 95
      glo = 1.e32
      ha = -glo
      if (ioffp .eq. 0) go to 30
      do 25 j=1,ny
        do 20 i=1,nx
          zz = z(i,j)
          if (zz .eq. spval) go to 20
          glo = amin1(zz,glo)
          ha = amax1(zz,ha)
   20   continue
   25 continue
      go to 45
   30 do 40 j=1,ny
        do 35 i=1,nx
          glo = amin1(z(i,j),glo)
          ha = amax1(z(i,j),ha)
   35   continue
   40 continue
!# modified: 10/9/2008
   45 if (fanc.LE.0.0) then 
         if (fanc.LT.0.0) crat = amax1(1.,-fanc)
         fanc = (ha-glo)/crat
         if (fanc.LE.0.0) go to 90
         p = 10.**(ifix(alog10(fanc)+500.)-500)
         fanc = aint(fanc/p)*p
      end if
      if (chi.EQ.clo) then     
         glo = aint(glo/fanc)*fanc
         ha = aint(ha/fanc)*fanc*(1.+sign(1.e-6,ha))
      end if    
!# modified: 10/9/2008
      do 80 k=1,nlm
        cc = glo+float(k-1)*fanc
        if (cc .gt. ha) go to 85
        kk = k
        cl(k) = cc
   80 continue
   85 ncl = kk
      cclo = cl(1)
      chi = cl(ncl)
      cinc = fanc
      return
   90 icnst = 1
      ncl = 1
      cclo = glo
      return
   95 cl(1) = glo
      ncl = 1
      return
      end
!*****CLIP
      subroutine clip
      character*132 cmdstr
      common/plt1/cmdstr

      cmdstr='clip'
      call filler
      return
      end
!*****CLIPBOX
      subroutine clipbox(xpts,ypts,npts)
! This routine creates a clip region bounded by arrays xpts and ypts
! xpts,ypts are in inches.  Uses first npts points of xpts and ypts.
      dimension xpts(npts),ypts(npts)
      character*132 cmdstr,scr
      common/plt1/cmdstr
      common/cnvcom/conver

      cmdstr=' '
      do nn=1,npts
        xx=xpts(nn)*conver
        yy=ypts(nn)*conver
        scr=' '
        write(scr,'('' '',f8.2,'' '',f8.2)')xx,yy
        lens=lenstr(scr,132)
        if(nn.eq.1) then
          cmdstr=scr(2:lens)    !Don't need initial space
        else
          lc=lenstr(cmdstr,132)
          if((lc+lens).gt.132) then
            call filler
            cmdstr=scr(1:lens)
          else
            cmdstr=cmdstr(1:lc)//scr(1:lens)
          endif
        endif
      enddo
      call filler
      nm1=npts-1
      write(cmdstr,'(i6,'' Cln'')')nm1
      cmdstr=cmdstr(1:lenstr(cmdstr,132))
      call filler
      return
      end
!*****COLBOX
      subroutine colbox(xpts,ypts,npts,red,green,blue)
! This routine fills region bounded by arrays xpts and ypts
!  with colors red,green,blue.  xpts,ypts are in inches.  Uses first npts
!  points of xpts and ypts
      dimension xpts(npts),ypts(npts)
      character*132 cmdstr,scr
      common/plt1/cmdstr
      common/cnvcom/conver
      common/colrcom/cred,cgreen,cblue,cgry

      if(red.ne.cred.or.green.ne.cgreen.or.blue.ne.cblue)                      &
         call setcolr(red,green,blue)

      cmdstr=' '
      do nn=1,npts
        xx=xpts(nn)*conver
        yy=ypts(nn)*conver
        scr=' '
        write(scr,'('' '',f8.2,'' '',f8.2)')xx,yy
        lens=lenstr(scr,132)
        if(nn.eq.1) then
          cmdstr=scr(2:lens)    !Don't need initial space
        else
          lc=lenstr(cmdstr,132)
          if((lc+lens).gt.132) then
            call filler
            cmdstr=scr(1:lens)
          else
            cmdstr=cmdstr(1:lc)//scr(1:lens)
          endif
        endif
      enddo

      nm1=npts-1
      scr=' '
      write(scr,'(i6,'' Fbn'')')nm1
      lc=lenstr(cmdstr,132)
      ls=lenstr(scr,132)
      if((lc+ls).gt.132) then
        call filler
        cmdstr=scr(1:ls)
      else
        cmdstr=cmdstr(1:lc)//scr(1:ls)
      endif
      call filler
      return
      end
!*****COLBOXC
      subroutine colboxc(xpts,ypts,npts,ioff,joff,red,green,blue)
! this routine fills region bounded by arrays xpts and ypts
!  with colors red,green,blue.  xpts,ypts are in inches.  Uses first npts
!  points of xpts and ypts
!  Called by concolr
      dimension xpts(npts),ypts(npts)
      character*132 cmdstr,scr
      common/plt1/cmdstr
      common/cnvcom/conver
      common/colrcom/cred,cgreen,cblue,cgry
      common/boxcom/dx,dy
      if(red.ne.cred.or.green.ne.cgreen.or.blue.ne.cblue)                      &
          call setcolr(red,green,blue)

      cmdstr=' '
      do nn=1,npts
        scr=' '
        xp=(ioff+xpts(nn))*dx*conver
        yp=(joff+ypts(nn))*dy*conver
        write(scr,'('' '',f8.2,'' '',f8.2)')xp,yp
        lens=lenstr(scr,132)
        if(nn.eq.1) then
          cmdstr=scr(2:lens)    !Don't need initial space
        else
          lc=lenstr(cmdstr,132)
          if((lc+16).gt.132) then
            call filler
            cmdstr=scr(1:lens)
          else
            cmdstr=cmdstr(1:lc)//scr(1:lens)
          endif
        endif
      enddo

      nm1=npts-1
      scr=' '
      write(scr,'(i6,'' Fbnc'')')nm1
      lc=lenstr(cmdstr,132)
      ls=lenstr(scr,132)
      if((lc+ls).gt.132) then
        call filler
        cmdstr=scr(1:ls)
      else
        cmdstr=cmdstr(1:lc)//scr(1:ls)
      endif
      call filler
      return
      end
!*****CONCOLR
      subroutine concolr(arr,idim,imax,jmax,xlen,ylen,cvalo,coloro,            &
      nvall,ioffp,spval)

!  This routine placed certain variables and arrays in double precision
!  to aid in making decisions regarding box corners, contour points, etc.

!  This routine produces color contour maps of array arr.
!  idim is leading dimension of arr in calling program
!  imax,jmax are i,j extent of arr to contour
!  xlen,ylen are lengths of contour box in x,y directions
!  cval is array containing color demarcations (similar to contour values
!  in conrec)
!  color is the array containing the RGB values for each of the contour values
!  nval is number of color demarcations

      parameter(intmax=100)   !Max number of contour intersections/grid box
      parameter(maxpt=100,maxcrv=50)
      parameter(maxfil=20)
      double precision xfill(maxfil),yfill(maxfil)
      double precision disttot,distptmn,distpt,difx,dify,slown,slnew,          &
      xp1,xp2,yp1,yp2,xp,yp,dl,dr,db,dt
      dimension xfills(maxfil),yfills(maxfil)
      character*132 cmdstr,curfnt
      character*132 ifrmt
      common/plt1/cmdstr
      common/cnvcom/conver
      common/fntcom/curfnt,ifntsz,nfont
      common/colrcom/cred,cgreen,cblue,cgry
      common/boxcom/dx,dy
      dimension arr(idim,jmax),cvalo(nvall),coloro(3,nvall)
      dimension cval(100),color(3,100),colormax(3)
      double precision xx,yy
      common/crvcomdp/xx(maxpt,maxcrv),yy(maxpt,maxcrv),rcval(maxcrv),         &
      nump(maxcrv),isvncv,ncrv
      dimension xa(5),ya(5)
      double precision smarr(2,2)
      dimension cnval(100)
      integer*2 isider(maxpt,maxcrv)

      dimension val(4)
      dimension direct(4)
      dimension xtt(4),ytt(4)

      common/conpar/ispec,ioffpp,spvall,ilegg,ilabb,nhii,ndeccn,nlbll,         &
      lscal,ldash,hgtlab

      data direct/1.,1.,-1.,-1./

!  Save current colors
      redsv=cred
      greensv=cgreen
      bluesv=cblue

!  Save current value of ilabb
      ilabsv=ilabb

      isvncv=1

!  Fill cval and color
      do nn=1,nvall
        cval(nn)=cvalo(nn)
        do kk=1,3
          color(kk,nn)=coloro(kk,nn)
        enddo
      enddo

!  Ensure that contours enclose total array range; if not, add extreme
!  contour and color values to either end.
      nval=nvall
      armin= 1.e20
      armax=-1.e20
      do jj=1,jmax
        do ii=1,imax
          if(arr(ii,jj).ne.spval) then
            armin=amin1(armin,arr(ii,jj))
            armax=amax1(armax,arr(ii,jj))
          endif
        enddo
      enddo

!     print *,'***CONCOLR*** concolr array range ',armin,armax
      if(cval(1).ge.armin) then
        do n=nval,1,-1
          cval(n+1)=cval(n)
          do kkk=1,3
            color(kkk,n+1)=color(kkk,n)
          enddo
        enddo
        cval(1)=-1.e20
        do kkk=1,3
          color(kkk,1)=color(kkk,2)
        enddo
        nval=nval+1
!       print *,'***CONCOLR*** Adding low-end cval,color',
!    +       cval(1),color(1,1),color(2,1),color(3,1)
!       print *,'Incrementing nval, nval now equals ',nval
      endif
      if(cval(nval).le.armax) then
        cval(nval+1)=1.e20
        do kkk=1,3
          color(kkk,nval+1)=color(kkk,nval)
        enddo
        nval=nval+1
!       print *,'***CONCOLR*** Adding high-end cval,color',
!    +         cval(nval),color(1,nval),color(2,nval),color(3,nval)
!       print *,'Incrementing nval, nval now equals ',nval
      endif

      colormax(1)=color(1,nval)
      colormax(2)=color(2,nval)
      colormax(3)=color(3,nval)

! Fill entire region with highest color level
      call setcolr(colormax(1),colormax(2),colormax(3))

      xa(1)=0.
      xa(2)=xlen
      xa(3)=xa(2)
      xa(4)=xa(1)
      ya(1)=0.
      ya(2)=ya(1)
      ya(3)=ylen
      ya(4)=ya(3)
      call filrgnc(xa,ya,4)

      imaxm=imax-1
      jmaxm=jmax-1
      dx=xlen/float(imaxm)
      dy=ylen/float(jmaxm)
!     ifrmt='(''/dx '',f11.7,i5,'' div 72 mul def /dy '',f11.7,i5,'//
!    +''' div 72 mul def'')'
      ifrmt='(''/dx '',f11.7,i5,'' div '',f7.4,'' mul def /dy '',f11.7,'       &
             //'i5,'' div '',f7.4,'' mul def'')'
      cmdstr=' '
      write(cmdstr,ifrmt)xlen,imaxm,conver,ylen,jmaxm,conver
      call filler
      cmdstr='/Dx dx def /Dy dy def'
      call filler

      xtt(1)=0.
      xtt(2)=1.
      xtt(3)=1.
      xtt(4)=0.
      ytt(1)=0.
      ytt(2)=0.
      ytt(3)=1.
      ytt(4)=1.

!  Outer box loop
      do 55 jj=1,jmax-1
        do 50 ii=1,imax-1

          val(1)=arr(ii,jj)
          val(2)=arr(ii+1,jj)
          val(3)=arr(ii+1,jj+1)
          val(4)=arr(ii,jj+1)

!  Don't bother with points completely surrounded by spval
          if(ioffp.ne.0.and.(val(1).eq.spval.and.val(2).eq.spval .and.         &
          val(3).eq.spval.and.val(4).eq.spval)) goto 50

!  Deal with partial "special" boxes. Fill box with gray level of average value.
          if(ioffp.ne.0.and.(val(1).eq.spval.or.val(2).eq.spval                &
           .or.val(3).eq.spval.or.val(4).eq.spval)) then
            sum=0.
            nc=0
            do nn=1,4
              if(val(nn).ne.spval) then
                sum=sum+val(nn)
                nc=nc+1
              endif
            enddo
            sum=sum/float(nc)
            do nn=nval,1,-1
              if(sum.le.cval(nn))nsav=nn
            enddo
            if(color(1,nsav).ne.cred.or.color(2,nsav).ne.cgreen.or.            &
            color(3,nsav).ne.cblue)                                            &
            call setcolr(color(1,nsav),color(2,nsav),color(3,nsav))
            write(cmdstr,'(f8.2,'' '',f8.2,'' Fb'')')dx*(ii-1)*conver,         &
                        dy*(jj-1)*conver
            call filler
            goto 50  !Done with this box
          endif

!  See if this box can be filled immediately (no contour intersections)

! Loop on gray levels so we don't have to constantly fill the same box with
!  increasing gray levels

          nsav=0
          do n=1,nval-1
            rcrit1=cval(n)
            rcrit2=cval(n+1)

!      Deal with certain scenarios immediately
!      if(val(1).ge.rcrit1.and.val(2).ge.rcrit1.and.val(3).ge.rcrit1.and.
!    1    val(4).ge.rcrit1.and.
!      Change below says that boxes at or below the contour value will have
!      the corresponding color value. This only matters for wide areas
!      having the exact contour value.
            if(val(1).gt.rcrit1.and.val(2).gt.rcrit1.and.                      &
               val(3).gt.rcrit1.and.val(4).gt.rcrit1.and.                      &
               val(1).le.rcrit2.and.val(2).le.rcrit2.and.                      &
               val(3).le.rcrit2 .and.val(4).le.rcrit2) nsav=n+1
          enddo

          if(nsav.ne.0) then
            if( colormax(1).ne.color(1,nsav).or.                               &
                colormax(2).ne.color(2,nsav).or.                               &
                colormax(3).ne.color(3,nsav))then

              if(color(1,nsav).ne.cred.or.color(2,nsav).ne.cgreen.or.          &
                 color(3,nsav).ne.cblue)                                       &
                 call setcolr(color(1,nsav),color(2,nsav),color(3,nsav))
              write(cmdstr,'(f8.2,'' '',f8.2,'' Fb'')')dx*(ii-1)*conver,       &
                                                       dy*(jj-1)*conver
              call filler
            endif
            goto 50    !Done with this box
          endif

          armin= 1.e20
          armax=-1.e20
          do n=1,4
            armin=amin1(armin,val(n))
            armax=amax1(armax,val(n))
          enddo

!  fill smarr array
          smarr(1,1)=arr(ii,jj)
          smarr(1,2)=arr(ii,jj+1)
          smarr(2,1)=arr(ii+1,jj)
          smarr(2,2)=arr(ii+1,jj+1)

!  Get all contours intersecting this box
          ncon=0
          do n=1,nval-1
            ctar=cval(n)
            if(ctar.ge.armin.and.ctar.le.armax) then
              ncon=ncon+1
              cnval(ncon)=ctar
              nbig=n+1
            endif
          enddo

!  First, fill entire box with color(nbig) to fill in next higher level
          if(color(1,nbig).ne.cred.or.color(2,nbig).ne.cgreen.or. color        &
          (3,nbig).ne.cblue)                                                   &
          call setcolr(color(1,nbig),color(2,nbig),color(3,nbig))
          write(cmdstr,'(f8.2,'' '',f8.2,'' Fb'')')dx*(ii-1)*conver,           &
                                                   dy*(jj-1)*conver
          call filler

!  Have ncon contours, now get intersecting points
          ilabb=-999
          call conqck(smarr,2,2,2,1.,1.,cnval,ncon)

!  Scan curves to determine contour orientation
!     thresh=1.e-7
          thresh=1.e-13 !dp
          xln=0.
          xrn=1.
          ybn=0.
          ytn=1.

!  Fill isider array
          do 15 nc=1,ncrv
            if(nump(nc).eq.0) goto 15
            do 10 np=1,nump(nc)
              xp=xx(np,nc)
              yp=yy(np,nc)
              dl=dabs(xp-xln)
              dr=dabs(xp-xrn)
              db=dabs(yp-ybn)
              dt=dabs(yp-ytn)
              if(dl.le.thresh) then
                if(dt.le.thresh) then
                  isider(np,nc)=3
                else
                  isider(np,nc)=4
                endif
              else if(dr.le.thresh) then
                if(db.le.thresh)then
                  isider(np,nc)=1
                else
                  isider(np,nc)=2
                endif
              else if(db.le.thresh) then
                if(dl.le.thresh)then
                  isider(np,nc)=4
                else
                  isider(np,nc)=1
                endif
              else if(dt.le.thresh) then
                if(dr.le.thresh)then
                  isider(np,nc)=2
                else
                  isider(np,nc)=3
                endif
              endif
   10       continue
   15     continue

          do 45 nc=1,ncrv
            if(nump(nc).eq.0) goto 45
!      Find what contour value we're on
            do nnn=1,nval
              if(rcval(nc).eq.cval(nnn)) then
                nlev=nnn
                goto 20
              endif
            enddo
   20       continue

            np=1
            npp=min0(np+1,nump(nc))
            xp1=xx(np,nc)
            yp1=yy(np,nc)
            xp2=xx(npp,nc)
            yp2=yy(npp,nc)
            ncstart=nc
            npstart=np

            xfill(1)=xp1
            yfill(1)=yp1
            if(npp.ne.np+1) then !Only one point in curve
              nplast=np
              nclast=nc
              ksider=isider(np,nc)
              isidenow=mod(ksider,4)+1 !Corner point, so start search on
                                        !next side
              nfill=1
            else
              xfill(2)=xp2
              yfill(2)=yp2
              nfill=2
              nplast=npp
              nclast=nc
              isidenow=isider(npp,nc)
            endif
   25       continue
            isnext=mod(isidenow,4)+1
            difx=xfill(nfill)-xtt(isnext)      !Check distance xtt,ytt(isnext)
            dify=yfill(nfill)-ytt(isnext)
            disttot=sqrt(difx*difx+dify*dify)

!  See if we encounter any points between here and xtt,ytt(isnext)
            ifind=0              !Flag set if we bump into a contour point on
            distptmn=1.e20       !the way to the next box corner
            do 35 ncc=1,ncrv
              do 30 npt=1,nump(ncc)
                if(ncc.eq.nclast.and.npt.eq.nplast) goto 35 !Don't think about y
                if(isider(npt,ncc).eq.isidenow) then !Found contour pt on same b
                  difx=xx(npt,ncc)-xfill(nfill)
                  dify=yy(npt,ncc)-yfill(nfill)
                  difxs=difx
                  difys=dify
!       Make sure we check in proper direction only
                  if(mod(isidenow,2).eq.1.and.sign(1.,difxs).ne.direct         &
                  (isidenow)) goto 35
                  if(mod(isidenow,2).eq.0.and.sign(1.,difys).ne.direct         &
                  (isidenow)) goto 35
                  distpt=sqrt(difx*difx+dify*dify)
                  if(distpt.le.disttot.and.distpt.le.distptmn)then !Make le to g
!                                                         advantage to
!                                                         contour pt rather
!                                                         than box corner.
!         We now have a point of another contour. Don't use this point if it's
!         where we are now, as this could lead us down another contour rather
!         than lead us back to our beginning point.
!         Changes made 1/4/95 kek
!          c-out next 4 lines and replace with subsequent do loop
!           if(distpt.eq.0.) then ! Don't use this point
!           print *,'possible other contour with distpt= 0.',npt,ncc
!           goto 7600
!          endif

                    if(distpt.eq.0.) then ! Don't use this point
!gdr20080102          print *,'possible other contour with distpt= 0.',        &
!gdr20080102          npt,ncc
!           Compare slopes and don't use if slope is greater than your own
                      slown=(yfill(nfill)-yfill(nfill-1)) /(xfill(nfill)       &
                      -xfill(nfill-1))
                      nother=3-npt
                      slnew=(yy(nother,ncc)-yy(npt,ncc)) /(xx                  &
                      (nother,ncc)-xx(npt,ncc))
!gdr20080102          print *,'slown,slnew= ',slown,slnew
                      if(slnew.gt.slown) goto 35
                    endif
                    ifind=1 !Bumped into a point, but keep checking to
                    ncuse=ncc !ensure it's the closest
                    npuse=npt
                    distptmn=distpt
                  endif
                endif
   30         continue
   35       continue
            if(ifind.eq.0) then !All clear to box corner
              nfill=nfill+1

!gdr20080102
              if(nfill.gt.maxfil) goto 45
!             if(nfill.gt.maxfil) then
!               print *,'***Error(1) in concolr***'
!               print *,'nc,rcval= ',nc,rcval(nc)
!               print 100,(xfill(kkk),yfill(kkk),kkk=1,nfill-1)
  100 format(1x,2f20.10)
!             endif
!gdr20080102

              xfill(nfill)=xtt(isnext)
              yfill(nfill)=ytt(isnext)
!       print *,'nfill1 ',nfill,xfill(nfill),yfill(nfill)
              nplast=-999
              nclast=-999
              isidenow=isnext
              goto 25
            else
              if(ncuse.eq.ncstart.and.npuse.eq.npstart) then
                goto 40 !We're back to beginning, so quit
              endif
!       Also check to see if found contour point is really the starting point
!       of another contour. This can occur if two contours intersect at a box
!       corner.
              arg1=xx(npuse,ncuse)-xx(npstart,ncstart)
              arg2=yy(npuse,ncuse)-yy(npstart,ncstart)
              distchk=sqrt(arg1*arg1+arg2*arg2)
              if(distchk.eq.0.) then
!          It's the same point, so we're done
                goto 40
              endif
              nfill=nfill+1


!gdr20080102
              if(nfill.gt.maxfil) goto 45
!             if(nfill.gt.maxfil) then
!               print *,'***Error(2) in concolr***'
!               print *,'nc,rcval= ',nc,rcval(nc)
!               print 100,(xfill(kkk),yfill(kkk),kkk=1,nfill-1)
!             endif
!gdr20080102

              xfill(nfill)=xx(npuse,ncuse)
              yfill(nfill)=yy(npuse,ncuse)
!       print *,'nfill2 ',nfill,xfill(nfill),yfill(nfill)
!  Follow contour to other side
!       nother=3-npuse           !changed 1/19/95 kek
              nother=min0(3-npuse,nump(ncuse))
              nfill=nfill+1

!gdr20080102
              if(nfill.gt.maxfil) then
                print *,'***Error(3) in concolr***'
!gdr20080102    print *,'nc,rcval= ',nc,rcval(nc)
                print 100,(xfill(kkk),yfill(kkk),kkk=1,nfill-1)
              endif
!gdr20080102

              xfill(nfill)=xx(nother,ncuse)
              yfill(nfill)=yy(nother,ncuse)
!       print *,'nfill3 ',nfill,xfill(nfill),yfill(nfill)
              isidenow=isider(nother,ncuse)
              nplast=nother
              nclast=ncuse
              goto 25
            endif

   40       continue

            do nnnn=1,nfill
              xfills(nnnn)=xfill(nnnn)
              yfills(nnnn)=yfill(nnnn)
            enddo
            call colboxc(xfills,yfills,nfill,ii-1,jj-1,color(1,nlev),          &
            color(2,nlev),color(3,nlev))
   45     continue

   50   continue
   55 continue

      ilabb=ilabsv
!  Reset colors
      if(redsv.ne.cred.or.greensv.ne.cgreen.or. bluesv.ne.cblue)               &
      call setcolr(redsv,greensv,bluesv)

      return
      end
!*****CONFILL
      subroutine confill(arr,idim,imax,jmax,xlen,ylen,cvalo,grylevo,           &
      nvall,ioffp,spval)

!  This routines places certain variables and arrays in double precision
!  to aid in making decisions regarding box corners, contour points, etc.

!  This routine produces grayscale contour maps of array arr.
!  idim is leading dimension of arr in calling program
!  imax,jmax are i,j extent of arr to contour
!  xlen,ylen are lengths of contour box in x,y directions
!  cval is array containing color demarcations (similar to contour values
!  in conrec)
!  grylev is array containing grayscale values for each of the contour values
!  nval is number of grayscale demarcations

      parameter(intmax=100)   !Max number of contour intersections/grid box
      parameter(maxpt=100,maxcrv=50)
      parameter(maxfil=20)
      double precision xfill(maxfil),yfill(maxfil)
      double precision disttot,distptmn,distpt,difx,dify,slown,slnew,          &
      xp1,xp2,yp1,yp2,xp,yp,dl,dr,db,dt
      dimension xfills(maxfil),yfills(maxfil)
      character*132 cmdstr
      character*132 ifrmt
      common/plt1/cmdstr
      common/cnvcom/conver
      common/colrcom/cred,cgreen,cblue,cgry
      common/boxcom/dx,dy
      dimension arr(idim,jmax),cvalo(nvall),grylevo(nvall)
      dimension cval(100),grylev(100)
      double precision xx,yy
      common/crvcomdp/xx(maxpt,maxcrv),yy(maxpt,maxcrv),rcval(maxcrv),         &
      nump(maxcrv),isvncv,ncrv
      dimension xa(5),ya(5)
      double precision smarr(2,2)
      dimension cnval(100)
      integer*2 isider(maxpt,maxcrv)

      dimension val(4)
      dimension direct(4)
      dimension xtt(4),ytt(4)

      common/conpar/ispec,ioffpp,spvall,ilegg,ilabb,nhii,ndeccn,nlbll,         &
      lscal,ldash,hgtlab

      data direct/1.,1.,-1.,-1./

      isvncv=1

!  Fill cval,grylev
      do nn=1,nvall
        cval(nn)=cvalo(nn)
        grylev(nn)=grylevo(nn)
      enddo

!  Ensure that contours enclose total array range; if not, add extreme
!  contour and grayscale values to either end
      nval=nvall
      armin= 1.e20
      armax=-1.e20
      do jj=1,jmax
        do ii=1,imax
          if(arr(ii,jj).ne.spval) then
            armin=amin1(armin,arr(ii,jj))
            armax=amax1(armax,arr(ii,jj))
          endif
        enddo
      enddo

!     print *,'***CONFILL*** confill array range ',armin,armax
      if(cval(1).ge.armin) then
        do n=nval,1,-1
          cval(n+1)=cval(n)
          grylev(n+1)=grylev(n)
        enddo
        cval(1)=-1.e20
        grylev(1)=grylev(2)
        nval=nval+1
!       print *,'***CONFILL*** Adding low-end cval,grylev', cval(1),
!    +  grylev(1)
!       print *,'Incrementing nval, nval now equals ',nval
      endif
      if(cval(nval).le.armax) then
        cval(nval+1)=1.e20
        grylev(nval+1)=grylev(nval)
        nval=nval+1
!       print *,'***CONFILL*** Adding high-end cval,grylev', cval(nval),
!    +  grylev(nval)
!       print *,'Incrementing nval, nval now equals ',nval
      endif

! Fill entire region with highest gray level
      grymax=grylev(nval)

      xa(1)=0.
      xa(2)=xlen
      xa(3)=xa(2)
      xa(4)=xa(1)
      ya(1)=0.
      ya(2)=ya(1)
      ya(3)=ylen
      ya(4)=ya(3)
      call filrgn(xa,ya,4,grymax)

!  Save current gray level
      grysav=cgry
!  Save current value of ilabb
      ilabsv=ilabb

      imaxm=imax-1
      jmaxm=jmax-1
      dx=xlen/float(imaxm)
      dy=ylen/float(jmaxm)

!     ifrmt='(''/dx '',f11.7,i5,'' div 72 mul def /dy '',f11.7,i5,'//
!    +''' div 72 mul def'')'
      ifrmt='(''/dx '',f11.7,i5,'' div '',f7.4,'' mul def /dy '',f11.7,'       &
             //'i5,'' div '',f7.4,'' mul def'')'
      cmdstr=' '
      write(cmdstr,ifrmt)xlen,imaxm,conver,ylen,jmaxm,conver
      call filler
      cmdstr='/Dx dx def /Dy dy def'
      call filler

      xtt(1)=0.
      xtt(2)=1.
      xtt(3)=1.
      xtt(4)=0.
      ytt(1)=0.
      ytt(2)=0.
      ytt(3)=1.
      ytt(4)=1.

!  Outer box loop
      do 65 jj=1,jmax-1
        do 60 ii=1,imax-1

          val(1)=arr(ii,jj)
          val(2)=arr(ii+1,jj)
          val(3)=arr(ii+1,jj+1)
          val(4)=arr(ii,jj+1)

!  Don't bother with points completely surrounded by spval
          if(ioffp.ne.0.and.(val(1).eq.spval.and.val(2).eq.spval .and.         &
          val(3).eq.spval.and.val(4).eq.spval)) goto 60

!  Deal with partial "special" boxes
          if(ioffp.ne.0.and.(val(1).eq.spval.or.val(2).eq.spval .or.val        &
          (3).eq.spval.or.val(4).eq.spval)) then

!       Fill box with gray level of average value
            sum=0.
            nc=0
            do nn=1,4
              if(val(nn).ne.spval) then
                sum=sum+val(nn)
                nc=nc+1
              endif
            enddo
            sum=sum/float(nc)
            do nn=NVAL,1,-1
              if(sum.LE.cval(nn))nsav=nn
            enddo
            if(grylev(nsav).ne.cgry) call setgry(grylev(nsav))
            write(cmdstr,'(f8.2,'' '',f8.2,'' Fb'')')dx*(ii-1)*conver,         &
                                                     dy*(jj-1)*conver
            call filler
            goto 60    !Done with this box
          endif

!  See if this box can be filled immediately (no contour intersections)

! Loop on gray levels so we don't have to constantly fill the same box with
!  increasing gray levels
          nsav=0
          do 10 n=1,nval-1
            rcrit1=cval(n)
            rcrit2=cval(n+1)

!      Deal with certain scenarios immediately
!      if(val(1).ge.rcrit1.and.val(2).ge.rcrit1.and.val(3).ge.rcrit1.and.
!    1    val(4).ge.rcrit1.and.
!      Change below says that boxes at or below the contour value will have
!      the corresponding grayscale value. This only matters for wide areas
!      having the exact contour value.
            if(val(1).gt.rcrit1.and.val(2).gt.rcrit1.and.val(3).gt.            &
            rcrit1 .and.val(4).gt.rcrit1.and. val(1).le.rcrit2.and.val         &
            (2).le.rcrit2.and.val(3).le.rcrit2 .and.val(4).le.rcrit2)          &
            nsav=n+1
   10     continue

          if(nsav.ne.0) then
            if(grylev(nsav).ne.grymax) then
              if(grylev(nsav).ne.cgry) call setgry(grylev(nsav))
              write(cmdstr,'(f8.2,'' '',f8.2,'' Fb'')')dx*(ii-1)*conver,       &
                                                       dy*(jj-1)*conver
              call filler
            endif
            goto 60    !Done with this box
          endif

          armin= 1.e20
          armax=-1.e20
          do n=1,4
            armin=amin1(armin,val(n))
            armax=amax1(armax,val(n))
          enddo

!  fill smarr array
          smarr(1,1)=arr(ii,jj)
          smarr(1,2)=arr(ii,jj+1)
          smarr(2,1)=arr(ii+1,jj)
          smarr(2,2)=arr(ii+1,jj+1)

!  Get all contours intersecting this box
          ncon=0
          do n=1,nval-1
            ctar=cval(n)
            if(ctar.ge.armin.and.ctar.le.armax) then
              ncon=ncon+1
              cnval(ncon)=ctar
              nbig=n+1
            endif
          enddo

!  First, fill entire box with grylev(nbig) to fill in next higher level
          if(grylev(nbig).ne.cgry) call setgry(grylev(nbig))
          write(cmdstr,'(f8.2,'' '',f8.2,'' Fb'')')dx*(ii-1)*conver,           &
                                                   dy*(jj-1)*conver
          call filler

!  Have ncon contours, now get intersecting points
          ilabb=-999
          call conqck(smarr,2,2,2,1.,1.,cnval,ncon)

!  Scan curves to determine contour orientation
!     thresh=1.e-7
          thresh=1.e-13 !dp
          xln=0.
          xrn=1.
          ybn=0.
          ytn=1.

!  Fill isider array
          do 20 nc=1,ncrv
            if(nump(nc).eq.0) goto 20
            do 15 np=1,nump(nc)
              xp=xx(np,nc)
              yp=yy(np,nc)
              dl=dabs(xp-xln)
              dr=dabs(xp-xrn)
              db=dabs(yp-ybn)
              dt=dabs(yp-ytn)
              if(dl.le.thresh) then
                if(dt.le.thresh) then
                  isider(np,nc)=3
                else
                  isider(np,nc)=4
                endif
              else if(dr.le.thresh) then
                if(db.le.thresh)then
                  isider(np,nc)=1
                else
                  isider(np,nc)=2
                endif
              else if(db.le.thresh) then
                if(dl.le.thresh)then
                  isider(np,nc)=4
                else
                  isider(np,nc)=1
                endif
              else if(dt.le.thresh) then
                if(dr.le.thresh)then
                  isider(np,nc)=2
                else
                  isider(np,nc)=3
                endif
              endif
   15       continue
   20     continue

!  Curve loop
          do 55 nc=1,ncrv
            if(nump(nc).eq.0) goto 55
!      Find what contour value we`re on
            do 25 nnn=1,nval
              if(rcval(nc).eq.cval(nnn)) then
                nlev=nnn
                goto 30
              endif
   25       continue
   30       continue

            np=1
            npp=min0(np+1,nump(nc))
            xp1=xx(np,nc)
            yp1=yy(np,nc)
            xp2=xx(npp,nc)
            yp2=yy(npp,nc)
            ncstart=nc
            npstart=np

            xfill(1)=xp1
            yfill(1)=yp1
            if(npp.ne.np+1) then !Only one point in curve
              nplast=np
              nclast=nc
              ksider=isider(np,nc)
              isidenow=mod(ksider,4)+1 !Corner point, so start Search on
                                        !next side
              nfill=1
            else
              xfill(2)=xp2
              yfill(2)=yp2
              nfill=2
              nplast=npp
              nclast=nc
              isidenow=isider(npp,nc)
            endif
   35       continue
            isnext=mod(isidenow,4)+1
            difx=xfill(nfill)-xtt(isnext)      !Check distance xtt,ytt(isnext)
            dify=yfill(nfill)-ytt(isnext)
            disttot=sqrt(difx*difx+dify*dify)
!  See if we encounter any points between here and xtt,ytt(isnext)
            ifind=0                 !Flag set if we bump into a contour point
            distptmn=1.e20          !on the way to the next box corner
            do 40 ncc=1,ncrv
              do 40 npt=1,nump(ncc)
                if(ncc.eq.nclast.and.npt.eq.nplast) goto 40 !Don't think about y
                if(isider(npt,ncc).eq.isidenow) then !Found contour pt on same b
                  difx=xx(npt,ncc)-xfill(nfill)
                  dify=yy(npt,ncc)-yfill(nfill)
!       Make sure we check in proper direction only
                  difxs=difx
                  difys=dify
                  if(mod(isidenow,2).eq.1.and.sign(1.,difxs).ne.direct         &
                  (isidenow)) goto 40
                  if(mod(isidenow,2).eq.0.and.sign(1.,difys).ne.direct         &
                  (isidenow)) goto 40
                  distpt=sqrt(difx*difx+dify*dify)

                  if(distpt.le.disttot.and.distpt.le.distptmn)then !Make le to g
!               advantage to contour pt rather than box corner.

!  We now have a point of another contour. Don't use this point if it's
!  where we are now, as this could lead us down another contour rather than
!  lead us back to our beginning point.
!  Changes made 1/4/95 kek
!  c-out next 4 lines and replace with subsequent do loop
!         if(distpt.eq.0.) then ! Don't use this point
!          print *,'possible other contour with distpt= 0.',npt,ncc
!          goto 7600
!         endif

                    if(distpt.eq.0.) then ! Don't use this point
!gdr20080102          print *,'possible other contour with distpt= 0.',        &
!gdr20080102          npt,ncc, rcval(nc)
!           Compare slopes and don't use if slope is greater than your own
                      slown=(yfill(nfill)-yfill(nfill-1)) /(xfill(nfill)       &
                      -xfill(nfill-1))
                      nother=3-npt
                      slnew=(yy(nother,ncc)-yy(npt,ncc)) /(xx                  &
                      (nother,ncc)-xx(npt,ncc))
!gdr20080102          print *,'slown,slnew= ',slown,slnew
                      if(slnew.gt.slown) goto 40
                    endif
                    ifind=1 !Bumped into a point, but keep checking to
                    ncuse=ncc !ensure it's the closest
                    npuse=npt
                    distptmn=distpt
                  endif
                endif
   40       continue
            if(ifind.eq.0) then !All clear to box corner
              nfill=nfill+1

!gdr20080102
              if(nfill.gt.maxfil) goto 45
!             if(nfill.gt.maxfil) then
!               print *,'***Error(1) in confill***'
!               print *,'nc,rcval= ',nc,rcval(nc)
!               print *,ii,jj,smarr(1,1),smarr(2,1),smarr(2,2),smarr(1,2)
!               do n=1,ncrv
!                 do nn=1,nump(n)
!                   print *,n,nn,rcval(n),xx(nn,n),yy(nn,n)
!                 enddo
!               enddo
!               print 100,(xfill(kkk),yfill(kkk),kkk=1,nfill-1)
! 100 format(1x,2f20.10)
!               stop
!             endif
!gdr20080102

              xfill(nfill)=xtt(isnext)
              yfill(nfill)=ytt(isnext)
              nplast=-999
              nclast=-999
              isidenow=isnext
              goto 35
            else
              if(ncuse.eq.ncstart.and.npuse.eq.npstart) then
                goto 45 !We're back to beginning, so quit
              endif
!       Also check to see if found contour point is really the starting point
!       of another contour. This can occur if two contours intersect at a box
!       corner.
              arg1=xx(npuse,ncuse)-xx(npstart,ncstart)
              arg2=yy(npuse,ncuse)-yy(npstart,ncstart)
              distchk=sqrt(arg1*arg1+arg2*arg2)
              if(distchk.eq.0.) then
!          It's the same point, so we're done
                goto 45
              endif
              nfill=nfill+1

!gdr20080102
              if(nfill.gt.maxfil) goto 45
!             if(nfill.gt.maxfil) then
!               print *,'***Error(2) in confill***'
!               print *,'nc,rcval= ',nc,rcval(nc)
!               print *,ii,jj,smarr(1,1),smarr(2,1),smarr(2,2),smarr(1,2)
!               do n=1,ncrv
!                 do nn=1,nump(n)
!                   print *,n,nn,rcval(n),xx(nn,n),yy(nn,n)
!                 enddo
!               enddo
!               print 100,(xfill(kkk),yfill(kkk),kkk=1,nfill-1)
!               stop
!             endif
!gdr20080102

              xfill(nfill)=xx(npuse,ncuse)
              yfill(nfill)=yy(npuse,ncuse)
!  Follow contour to other side
!       nother=3-npuse    !Changed 1/19/95 kek
              nother=min0(3-npuse,nump(ncuse))
              nfill=nfill+1

!gdr20080102
              if(nfill.gt.maxfil) goto 45
!             if(nfill.gt.maxfil) then
!               print *,'***Error(3) in confill***'
!               print *,'nc,rcval= ',nc,rcval(nc)
!               print *,ii,jj,smarr(1,1),smarr(2,1),smarr(2,2),smarr(1,2)
!               do n=1,ncrv
!                 do nn=1,nump(n)
!                   print *,n,nn,rcval(n),xx(nn,n),yy(nn,n)
!                 enddo
!               enddo
!               print 100,(xfill(kkk),yfill(kkk),kkk=1,nfill-1)
!               stop
!             endif
!gdr20080102

              xfill(nfill)=xx(nother,ncuse)
              yfill(nfill)=yy(nother,ncuse)
              isidenow=isider(nother,ncuse)
              nplast=nother
              nclast=ncuse
              goto 35
            endif

   45       continue
            do nnnn=1,nfill
              xfills(nnnn)=xfill(nnnn)
              yfills(nnnn)=yfill(nnnn)
            enddo
            call fillboxc(xfills,yfills,nfill,ii-1,jj-1,grylev(nlev))
   50       continue
   55     continue

   60   continue
   65 continue

      ilabb=ilabsv
!  Reset gray level
      if(grysav.ne.cgry)call setgry(grysav)
      return
      end
!*****CONQCK
      subroutine conqck(z,mx,nx,ny,xlen,ylen,cval,nval)
!  "quick" version of conrec called by confill and concolr
!  contour plots of array z

      parameter(maxpt=100,maxcrv=50)
      double precision z(mx,ny)
      dimension cl(100)
      dimension cval(*)
      double precision xor,xsc,yor,ysc
      common/scalesdp/xor,xsc,yor,ysc
      common/contyp/ispcon,idsh
      common/konval/conval(100),numcon
      double precision xx,yy
      common/crvcomdp/xx(maxpt,maxcrv),yy(maxpt,maxcrv),rcval(maxcrv),         &
      nump(maxcrv),isvncv,ncrv

      data mscal/0/
      data ioffd/0/
      data flo,hi,finc/0., 0., 0./

      ncrv=0
      idsh=0

      gl = flo
      ha = hi
      gp = finc
      l=mx
      m=nx
      n=ny

!  Set scales
      xor=1.
      yor=1.
      xsc=float(nx-1)/xlen
      ysc=float(ny-1)/ylen

!  Set contour levels.

      ncl=nval
      do i=1,ncl
        cl(i)=cval(i)
      enddo
      icnst=0
      gl=cl(1)
      ha=cl(ncl)
      numcon=ncl
      do nn=1,numcon
        conval(nn)=cl(nn)
      enddo

      cmx=amax1(abs(cl(1)),abs(cl(ncl)))
      iscl=0
      if(mscal.eq.0.and.cmx.ne.0.)iscl=498+ifix(alog10(cmx)-500.)
   10 scal=10.**(-iscl)
      c1prt=cl(1)*scal
      c2prt=cl(ncl)*scal
      cldif=(cl(2)-cl(1))*scal

! find major and minor lines

      nml=0
      nmlp1=nml+1

! set up label scaling
      ioffdt = ioffd
      if(ioffdt.eq.0) go to 15
      if(gl.ne.0.0.and.(abs(gl).lt.0.1.or.abs(gl).ge.1.e5))ioffdt=1
      if(ha.ne.0.0.and.(abs(ha).lt.0.1.or.abs(ha).ge.1.e5))ioffdt=1
      ash=10.**(3-ifix(alog10(amax1(abs(gl),abs(ha),abs(gp)))-500.)-500)
   15 if (ioffdt.eq.0) ash=1.
      do i=1,ncl
        contr = cl(i)
        ilabl=0
        call stlinedp (z,mx,nx,ny,contr,ilabl,scal)
      enddo
      return
      end
!*****CONREC
      subroutine conrec(z,mx,nx,ny,XXLEN,YYLEN,cval,nval)
!  contour plots of array z
!  z is dimensioned (mx,my) in calling program where ny.le.my
!  nx and ny are the number of values to be contoured in the x and y directions
!  xxlen and yylen are the lengths of the x and y axes in inches. if one or
!  both are less than zero, the contour values will be drawn along the contour
!  lines, rather than at the ends, using a simple minded approach.
!  xlen and ylen are the lengths of the x and y axes in inches
!  program assumes x=1., y=1. at origin of plot
!                  x=nx at xlen inches
!                  y=ny at ylen inches
!  program does not alter current plot origin
!  cval (1 to nval) are function values to contour
!  if nval=0, program chooses values to contour
!
      parameter(maxpt=3000,maxcrv=350)
      character*132 cmdstr
      character*80 scr
      character*80 iscr
      dimension z(mx,ny),cl(100),cwork(100)
      dimension cval(*)
      common/plt1/cmdstr
      common/conre1/ioffp,spval
      common/scales/xor,xsc,yor,ysc
      common/contyp/ispcon,idsh
      common/conpar/ispec,ioffpp,spvall,ilegg,ilabb,nhii,ndeccn,nlbll,         &
      mscall,ldsh,hgtlab
      common/konval/conval(100),numcon
      common/crvcom/xx(maxpt,maxcrv),yy(maxpt,maxcrv),rcval(maxcrv),           &
      nump(maxcrv),isvncv,ncrv

      data mscal/0/
      data ileg/1/
      data sizem,sizep,nla,nlm,ilab/ 2.5, 1., 16, 40, 1/

      data ioffd,nulbll/0,3/
      data flo,hi,finc,nhi/ 0., 0., 0.,  0/

      data ldash/0/

      isvncv=0
      ncrv=0
      if(ispec.ne.0) then   !user specified labelling and min/max plotting
        ioffp=ioffpp
        spval=spvall
        ilab=ilabb
        nhi=nhii
        ileg=ilegg
        nulbll=abs(nlbll)   !nlbll could be <0 (draw drlin2 labels with pen 2)
        if(nulbll.eq.999) nulbll=0  !case of nulbll=0, but use pen 2
        mscal=mscall
        ldash=ldsh
      endif

      idsh=ldash

      xlen=abs(xxlen)
      ylen=abs(yylen)

!  Set ispcon=1 always, regardless of xlen,ylen sign. This will label
!  contours along the contour lines.
      ispcon=1

      gl = flo
      ha = hi
      gp = finc
      l=mx
      m=nx
      n=ny

!  Set scales
      xor=1.
      yor=1.
      xsc=float(nx-1)/xlen
      ysc=float(ny-1)/ylen
!
!  Set contour levels.
!
      if (nval.eq.0) then
        call clgen(z,mx,nx,ny,gl,ha,gp,nla,nlm,cl,ncl,icnst)
      else
        ncl=nval
        do i=1,ncl
          cl(i)=cval(i)
        enddo
        icnst=0
        gl=cl(1)
        ha=cl(ncl)
      endif
      numcon=ncl
      do nn=1,numcon
        conval(nn)=cl(nn)
      enddo

!      if(ilabb.ne.-999)write (6,110) (cl(i),i=1,ncl)
!  110 FORMAT('0','FUNCTION CONTOURED FOR FOLLOWING VALUES'/(1X,10E12.5))
      cmx=amax1(abs(cl(1)),abs(cl(ncl)))
      iscl=0
      if(mscal.eq.0.and.cmx.ne.0.)iscl=498+ifix(alog10(cmx)-500.)
   10 scal=10.**(-iscl)
      c1prt=cl(1)*scal
      c2prt=cl(ncl)*scal
      cldif=(cl(2)-cl(1))*scal
      sizen=.13
      sizen=.11    !decrease from lxy value
!  get min and max values of z for legend print
      zmin=1.e20
      zmax = -zmin
      do 20 j=1,ny
        do 15 i=1,nx
          zz = z(i,j)
          if (ioffp .ne. 0.and. zz .eq. spval) go to 15  !Skip special values
          zmin = amin1(zz,zmin)
          zmax = amax1(zz,zmax)
   15   continue
   20 continue
!      if(ilabb.ne.-999)print 120, zmin,zmax
!  120 FORMAT('0','MIN AND MAX CONTOURED ARRAY VALUES ARE ',2E15.5)

      if(ilabb.ne.-999.and.ileg.ne.0) then
        if (icnst .ne. 0) then
          call keksymc(0.00,-.25,sizen,'CONSTANT FIELD DETECTED', 0.,23,       &
          0)
          call keksymc(0.00,-.45,sizen,'CONSTANT VALUE= ',0.,16,0)
          call kekexp(999.,999.,sizen,gl,0.,4,0)
          go to 30
        else
          call keksymc(0.00,-.25,sizen,'CMIN: ',0.,6,0)
          call kekexp(1.5,-.25,sizen,c1prt,0.,2,2)
          call keksymc(0.,-.45,sizen,'CMAX: ',0.,6,0)
          call kekexp(1.5,-.45,sizen,c2prt,0.,2,2)
          call keksymc(1.7,-.25,sizen,'STEP :',0.,6,0)
          call kekexp(3.3,-.25,sizen,cldif,0.,2,2)
          if (iscl.ne.0) then
            write(scr,'(''10**'',i10)')iscl
            call blkstp(scr,80,scr,nch)
            read(scr,'(20a4)')iscr
            call keksymc(1.7,-.45,sizen,'SCALE:',0.,6,0)
            call keksymc(3.3,-.45,sizen,iscr,0.,nch,2)
          endif
          call keksymc(3.5,-.25,sizen,'ZMIN: ',0.,6,0)
          call kekexp(5.0,-.25,sizen,zmin,0.,2,2)
          call keksymc(3.5,-.45,sizen,'ZMAX: ',0.,6,0)
          call kekexp(5.0,-.45,sizen,zmax,0.,2,2)
        endif
      endif

! Find major and minor lines

      if(ilab.ne.0.and.ilabb.ne.-999) then
        call reord (cl,ncl,cwork,nml,nulbll+1)
      else
        nml=0
      endif
      nmlp1=nml+1

! Set up label scaling

      ioffdt = ioffd
      if (ioffdt.eq.0) go to 25
      if (gl.ne.0.0.and.(abs(gl).lt.0.1.or.abs(gl).ge.1.e5))ioffdt=1
      if (ha.ne.0.0.and.(abs(ha).lt.0.1.or.abs(ha).ge.1.e5))ioffdt=1
      ash=10.**(3-ifix(alog10(amax1(abs(gl),abs(ha),abs(gp)))-500.)-500)
   25 if (ioffdt.eq.0) ash=1.
      do i=1,ncl
        contr = cl(i)
        ilabl=1
        if (i.gt.nml) ilabl=0
        cmdstr='S'  !Stroke previous path
        call filler
        call stline (z,mx,nx,ny,contr,ilabl,scal)
      enddo

!  Stroke final path just to be sure
      cmdstr='S'
      call filler

! Find relative minimums and maximums if wanted, and mark values if wanted.
      if (nhi .eq. 0..and.ilabb.ne.-999) call minmax (z,mx,nx,ny,sizem,        &
      ash,ioffdt,scal)
      if (nhi .gt. 0.and.ilabb.ne.-999) call pntval (z,mx,nx,ny,sizep,         &
      ash,ioffdt,scal)
      return
   30 write (6,130) gl
  130 format('0','CONSTANT FIELD DETECTED, VALUE=',E15.7)
      return
      end
!*****CURVE
      subroutine curve(x1,y1,x2,y2,x3,y3,x4,y4,contin)
!  This routines uses PostScript operator curveto
!  The curve is not "stroked" implicitly.
!  8/15/95  kek
!      Add logical contin: If contin=.true.,do not set initial point
!       explicitly.

      character*132 cmdstr
      character*132 ifrmt
      common/plt1/cmdstr
      common/cnvcom/conver
      logical contin

      if(.not.contin)call plot(x1,y1,3)
      ifrmt='(f8.2,'' '',f8.2,'' '',f8.2,'' '',f8.2,'' '',f8.2,'//             &
      ''' '',f8.2,   '' curveto'')'
      cmdstr=' '
      write(cmdstr,ifrmt)x2*conver,y2*conver,x3*conver,y3*conver,              &
                         x4*conver,y4*conver
      call filler
      return
      end
!*****DRLIN2
      subroutine drlin2 (z,l,mm,nn,ilabl,scal)
!  This routine differs from old drlin2 in that points comprising contours
!  are saved for a call by confill with ilabb=-999.
! this routine traces a contour line when given the beginning by stline.
! transformations can be added by deleting the statement functions for
! fx and fy in drline and minmax and adding external functions.
! x=1. at z(1,j), x=float(m) at z(m,j). x takes on non-integer values.
! y=1. at z(i,1), y=float(n) at z(i,n). y takes on non-integer values.

! **************************************************************
! ARL Addition
! **************************************************************
! igis>=1 means to dump lat/lon pairs of contours to file for GIS
      common/gis/gisout,gisatt,gistmp,igis
      character*80 gisout,gisatt,gistmp
! **************************************************************

! rrd parameter (ntotpt=10000)
      parameter (ntotpt=20000)
      parameter(maxpt=3000,maxcrv=350)
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      common/plt1/cmdstr
      dimension z(l,nn)

      dimension xsav(ntotpt),ysav(ntotpt),s(ntotpt),ipuord(ntotpt)
      dimension isav(ntotpt),jsav(ntotpt)
      dimension nsava(ntotpt),smaxa(ntotpt)


      common/pltparam/curlin
      common/crvcom/xx(maxpt,maxcrv),yy(maxpt,maxcrv),rcval(maxcrv),           &
      nump(maxcrv),isvncv,ncrv

      common /conre2/ ix ,iy ,idx ,idy , is ,iss ,np ,cv , inx(8) ,iny         &
      (8) ,ir(10000) ,nr

      common/conre1/ioffp,spval
      common/scales/xor,xsc,yor,ysc

      common/conpar/ispec,ioffpp,spvall,ilegg,ilabb,nhii,ndeccn,nlbll,         &
      lscal,ldash,hgtlab

      common/konval/conval(100),numcon
      common/kkplot/szrat
      logical         ipen       ,ipeno
      data ipen,ipeno/.true.,.true./
!  hgt is height of label
!  slab1f is fraction of total contour length to skip before first label on
!  contour
!  dslab is distance between labels on same contour

      data hgt/.1/,dslab/3.5/,slab1f/.1/,wide/1./,big/1.e30/  !.857=wide

!     data ncrv/0/

      fx(x,y) = x
      fy(x,y) = y
      lc16(k) = k*65536
!     lc16(k) = k*'200000'O   !VMS
!     lc16(k) = k*8#200000    !MS Fortran
      cfcn(p1,p2) = (p1-cv)/(p1-p2)

      if(isvncv.eq.0)ncrv=0    !Don't save curve numbers

!  Allow user to specify height of label,; if 0, use .1
      if(hgtlab.ne.0.)hgt=hgtlab

      ncrvs=ncrv

!      pi=4.*abs(atan(1.)) ! no need to use transcendental function to find pi
      m = mm
      n = nn
      if(ispec.ne.0) then
        ndec=ndeccn
      else
        ndec=1
      endif
      cvs=cv*scal
      xmin=1.e20
      ymin=1.e20
      nsav=0

      if(cvs.eq.0.) then
        nchar=2+ndec
      else
        nchar=max1(1.,alog10(abs(cvs))+1.)+1+ndec
      endif
      if(cvs.lt.0.)nchar=nchar+1

! ! Add .05" for space on either side of label
      width=hgt*szrat*nchar+.13
      hhgt=hgt*wide/2.
!# modified: 10/9/2008
      if (ioffp .ne. 0) then
        jump1=35
        jump2=60
      else
        jump1=45
        jump2=65
      endif
!# modified: 10/9/2008
   10 ix0 = ix
      iy0 = iy
      is0 = is
      if (ioffp .eq. 0) go to 15
      ix2 = ix+inx(is)
      iy2 = iy+iny(is)
      ipen = z(ix,iy).ne.spval .and. z(ix2,iy2).ne.spval
      ipeno = ipen
   15 if (idx .ne. 0) then
        y = iy
        isub = ix+idx
        x = cfcn(z(ix,iy),z(isub,iy))*float(idx)+float(ix)
      else
        x = ix
        isub = iy+idy
        y = cfcn(z(ix,iy),z(ix,isub))*float(idy)+float(iy)
      endif
   20 xx1=(fx(x,y)-xor)/xsc
      yy1=(fy(x,y)-yor)/ysc
      if(ipen) then
        if(nsav.eq.0) then
          nsav=nsav+1
          if(nsav.gt.ntotpt) then
            print *,'nsav too big in drlin2, program abandoned'
            stop
          endif
          xsav(nsav)=xx1
          ysav(nsav)=yy1
          isav(nsav)=ix
          jsav(nsav)=iy
          ipuord(nsav)=3
          xmin=amin1(xmin,xx1)
          ymin=amin1(ymin,yy1)
        else if(abs(xx1-xsav(nsav)).gt.1.e-6.or. abs(yy1-ysav(nsav)).gt.       &
        1.e-6) then
          nsav=nsav+1
          if(nsav.gt.ntotpt) then
            print *,'nsav too big in drlin2, program abandoned'
            stop
          endif
          xsav(nsav)=xx1
          ysav(nsav)=yy1
          isav(nsav)=ix
          jsav(nsav)=iy
          ipuord(nsav)=3
          xmin=amin1(xmin,xx1)
          ymin=amin1(ymin,yy1)
        endif
        if(ilabl.eq.1) then
          cvs=cv*scal
        endif
      endif
   25 is = is+1
      if (is .gt. 8) is = is-8
      idx = inx(is)
      idy = iny(is)
      ix2 = ix+idx
      iy2 = iy+idy
      if (iss .ne. 0) go to 30
      if (ix2.gt.m .or. iy2.gt.n .or. ix2.lt. 1.or. iy2.lt.1) go to 80
   30 if ((cv-z(ix2,iy2)).le.0) then
        is = is+4
        ix = ix2
        iy = iy2
        go to 25
      else
        if(mod(is,2).eq.0) goto 25
      endif
!# modified: 10/9/2008
      if (jump1.EQ.35) then
         isbig = is+(8-is)/6*8
         ix3 = ix+inx(isbig-1)
         iy3 = iy+iny(isbig-1)
         ix4 = ix+inx(isbig-2)
         iy4 = iy+iny(isbig-2)
         ipeno = ipen
         if (iss .EQ. 0) then      
            if (ix3.gt.m .or. iy3.gt.n .or. ix3.lt. 1.or. iy3.lt.1) go to 80
            if (ix4.gt.m .or. iy4.gt.n .or. ix4.lt. 1.or. iy4.lt.1) go to 80
         end if
         ipen = z(ix,iy).ne.spval .and. z(ix2,iy2).ne.spval .and. z        &
               (ix3,iy3).ne.spval .and. z(ix4,iy4).ne.spval
      end if
!# modified: 10/9/2008
      if (idx .eq. 0) go to 50
      y = iy
      isub = ix+idx
      x = cfcn(z(ix,iy),z(isub,iy))*float(idx)+float(ix)
      go to 55
   50 x = ix
      isub = iy+idy
      y = cfcn(z(ix,iy),z(ix,isub))*float(idy)+float(iy)
   55 continue 
!# modified: 10/9/2008
      if (jump2.EQ.65) go to 65
!# modified: 10/9/2008
      if (.not.ipen) go to 70
      if (ipeno) go to 65

! end of line segment

      xx1=(fx(xold,yold)-xor)/xsc
      yy1=(fy(xold,yold)-yor)/ysc
      nsav=nsav+1
      if(nsav.gt.ntotpt) then
        print *,'nsav too big in drlin2, program abandoned'
        stop
      endif
      xsav(nsav)=xx1
      ysav(nsav)=yy1
      isav(nsav)=ix
      jsav(nsav)=iy
      ipuord(nsav)=3
      xmin=amin1(xmin,xx1)
      ymin=amin1(ymin,yy1)
      if(ilabl.eq.1) then
        cvs=cv*scal
      endif

! continue line segment

   65 xx2=(fx(x,y)-xor)/xsc
      yy2=(fy(x,y)-yor)/ysc
      if(nsav.eq.0) then
        nsav=nsav+1
        if(nsav.gt.ntotpt) then
          print *,'nsav too big in drlin2, program abandoned'
          stop
        endif
        xsav(nsav)=xx2
        ysav(nsav)=yy2
        isav(nsav)=ix
        jsav(nsav)=iy
        ipuord(nsav)=2
        xmin=amin1(xmin,xx2)
        ymin=amin1(ymin,yy2)
      elseif(abs(xx2-xsav(nsav)).gt.1.e-6.or.abs(yy2-ysav(nsav)).gt.1.e-6) then
        nsav=nsav+1
        if(nsav.gt.ntotpt) then
          print *,'nsav too big in drlin2, program abandoned'
          stop
        endif
        xsav(nsav)=xx2
        ysav(nsav)=yy2
        isav(nsav)=ix
        jsav(nsav)=iy
        ipuord(nsav)=2
        xmin=amin1(xmin,xx2)
        ymin=amin1(ymin,yy2)
      endif

   70 xold = x
      yold = y
      if (is .ne. 1) go to 75
      np = np+1
      if (np .gt. nr) go to 80
      ir(np) = lc16(ix)+iy
   75 if (iss .eq. 0) go to 25
      if (ix.ne.ix0 .or. iy.ne.iy0 .or. is.ne.is0) go to 25

! end of line

   80 continue

!  ---------------------------------------------------------------------
!  next section is modfied jimar code
!
!                 The arrays X, Y are now complete. calculate
!                 the distance S along the contour. Start plotting.
      if(ilabb.eq.-999.or.ilabl.ne.1) then
        if(nsav.eq.0) return
        ncs=0
   85   ncs=ncs+1
        if(ncs.gt.nsav) goto 90
        if(ipuord(ncs).eq.3) then
          if(ncrv.ne.ncrvs.and.ilabb.ne.-999) call sldcrv(xx(1,ncrv),yy        &
          (1,ncrv),npts,0.)
          npts=1
          ncrv=ncrv+1
          rcval(ncrv)=cv
          xx(npts,ncrv)=xsav(ncs)
          yy(npts,ncrv)=ysav(ncs)
          nump(ncrv)=npts
        elseif(ipuord(ncs).eq.2) then
          npts=npts+1
          xx(npts,ncrv)=xsav(ncs)
          yy(npts,ncrv)=ysav(ncs)
          nump(ncrv)=npts
        endif
        goto 85
   90   continue
        if(npts.ne.1) then
          nump(ncrv)=npts
!       **************************************************************
!       ARL Addition

          if(igis.ge.1)call gisdump(xx(1,ncrv),yy(1,ncrv),npts,cv)
!       **************************************************************
          if(ilabb.ne.-999)call sldcrv(xx(1,ncrv),yy(1,ncrv),npts,0.)
        endif
        return
      endif
!-----------------------------------------------------------------------------
!  find out what contour level we're on
      do nc=1,numcon
        if(abs(cv-conval(nc)).lt.1.e-6) then
          lev=nc
          goto 95
        endif
      enddo
      print *,'Cannot find current contour level, program abandoned.'
      stop
   95 continue

      nsav0=nsav
      s(1)=0.
      if(nsav .le. 1) return

!  calculate s(k)'s normally
      smaxx=-1.e20
      do k=2,nsav0
        km1 = k - 1
        kp1 = k + 1
        dxx = xsav(k)-xsav(km1)
        dyy = ysav(k)-ysav(km1)
        s(k) = s(km1) + sqrt(dxx*dxx+dyy*dyy)
      enddo

      do k=2,nsav0
        nsava(k)=nsav0
        smaxa(k)=s(nsav0)
      enddo

!  Check for non-continuous contours and reset nsava,smaxa if necessary
      k1=1
      do k=2,nsav0
        km1 = k - 1
        kp1 = k + 1
        if(ipuord(k).eq.2.and.ipuord(kp1).eq.3) then
          do kk=k1,k
!         reset k's in this segment
            nsava(kk)=k
            smaxa(kk)=s(k)
          enddo
          k1=kp1
        endif
      enddo

      slab1 = smaxa(1)*slab1f !fraction of total length
      stest = dslab - slab1 !set so first label is at slab1

      k = 1
      if(ipuord(k).eq.3) then
        cmdstr='S'
        call filler
        call movet(xsav(k),ysav(k))
      else
        call linet(xsav(k),ysav(k))
      endif

! Check conditions for labelling.

!-----------------------------------------------------------------------------

  100 continue                    ! k loop
      nsav=nsava(k)
      smax = smaxa(k)

      if(ndec.le.-2)go to 125

!  Don't label if not enough contour left
      if(smax-s(k).le.width)go to 125
      km1= max0(k-1,1)
      stest = stest + s(k)-s(km1)
      if(stest.lt.dslab)go to 125
      kp1=k+1
      if(lev.eq.1) go to 105
      dlev=abs(conval(lev)-conval(lev-1))

!  Is there enough space between adjacent contours?

      i=(xsav(k)-xmin)*xsc+1.01
      i=min0(i,mm-1)
      j=(ysav(k)-ymin)*ysc+1.01
      j=min0(j,nn-1)

      i=isav(k)
      j=jsav(k)
      i=min0(i,mm-1)
      j=min0(j,nn-1)

      dzdx=(z(i+1,j)-z(i,j))*xsc
      if( dzdx .ge. big ) go to 125
      dzdy=(z(i,j+1)-z(i,j))*ysc
      if( dzdy .ge. big ) go to 125
      dzdg=sqrt(dzdx*dzdx+dzdy*dzdy)
      if(dzdg.eq.0.) go to 105
      cspace=dlev/dzdg
      if(cspace.lt.hgt/2.) go to 125 !label drawn at mid-contour

  105 continue

!     Check if enough space for contour.
!     If so, check that curvature is not too great.
!     Determine # of points to drop to fit label.

      do 110 kk = kp1,nsav
        dxx = xsav(kk)-xsav(k)
        dyy = ysav(kk)-ysav(k)
        space = sqrt(dxx*dxx+dyy*dyy)
        if(space.lt.width)go to 110
        ark=s(kk)-s(k)
        if(ark.eq.0.0) goto 125
        if(space/ark.lt.0.80)goto 125
        go to 115
  110 continue
      go to 125

!  Draw the label
  115 stest=0.
      kkm1 = kk - 1

!  Determine how much of the available space is actually needed.

      do kkk=1,5
        t=0.2*kkk
        dxx=xsav(kkm1)+(xsav(kk)-xsav(kkm1))*t-xsav(k)
        dyy=ysav(kkm1)+(ysav(kk)-ysav(kkm1))*t-ysav(k)
        space=sqrt(dxx*dxx+dyy*dyy)
        if(space.gt.width) goto 120
      enddo

  120 xendl = xsav(k) + width*dxx/space
      yendl = ysav(k) + width*dyy/space

      if(dxx.eq.0) then
        angle = 90.
        if(dyy.lt.0.)angle=270.
      else
        angle = atan(dyy/dxx)/rdpdeg
      endif

      if(dxx.ge.0.0)then
        xlab = xsav(k)+ hhgt*(dxx+dyy)/space
        ylab = ysav(k)+ hhgt*(dyy-dxx)/space
      else
        xlab = xendl-hhgt*(dxx+dyy)/space
        ylab = yendl-hhgt*(dyy-dxx)/space
      endif

      call keknum(xlab,ylab,hgt,cvs,angle,ndec,0)

      call movet(XENDL,YENDL)
      call linet(XSAV(KK),YSAV(KK))
      k=kk


!  Plot the segment from xk,yk to  xsav(k+1),ysav(k+1).

  125 continue
      if((k+1).le.nsav0) then
        if(ipuord(k+1).eq.3) then
          cmdstr='S'
          call filler
          call movet(xsav(k+1),ysav(k+1))
        else
          call linet(xsav(k+1),ysav(k+1))
        endif
      endif

      k=k+1
      if( k .lt. nsav0 ) go to 100

      return
      end
!*****DRLIN2DP
      subroutine drlin2dp (z,l,mm,nn,ilabl,scal)
!  Called by stlinedp (which is called by conqck)

      parameter(maxpt=100,maxcrv=50)

      double precision x,y,xx1,yy1,xx2,yy2,xmin,ymin,xold,yold, xx,yy,         &
      xsav,ysav
      double precision fx,fy,cfcn,p1,p2
      double precision z(l,nn)
      double precision cv
      common/conre3/cv
      dimension xsav(maxpt),ysav(maxpt),ipuord(maxpt)
      dimension isav(maxpt),jsav(maxpt)
      common/crvcomdp/xx(maxpt,maxcrv),yy(maxpt,maxcrv),rcval(maxcrv),         &
      nump(maxcrv),isvncv,ncrv

      common /conre2/ ix ,iy ,idx ,idy , is ,iss ,np ,cvdummy , inx(8) ,       &
      iny(8) ,ir(10000) ,nr

      common/conre1/ioffp,spval
      double precision xor,xsc,yor,ysc
      common/scalesdp/xor,xsc,yor,ysc

      common/conpar/ispec,ioffpp,spvall,ilegg,ilabb,nhii,ndeccn,nlbll,         &
      lscal,ldash,hgtlab

      common/kkplot/szrat
      logical ipen,ipeno

!     **************************************************************
!     ARL Addition
!     **************************************************************
      common/gis/gisout,gisatt,gistmp,igis
      character*80 gisout,gisatt,gistmp
!     **************************************************************

      data ipen,ipeno/.true.,.true./

!  hgt is height of label
      data hgt/.1/wide/1./  !.857=wide

!     data ncrv/0/

      fx(x,y) = x
      fy(x,y) = y

      lc16(k) = k*65536
!     lc16(k) = k*'200000'O   !VMS
!     lc16(k) = k*8#200000    !MS
      cfcn(p1,p2) = (p1-cv)/(p1-p2)

      if(isvncv.eq.0)ncrv=0    !Don't save curve numbers

!  Allow user to specify height of label,; if 0, use .1
      if(hgtlab.ne.0.)hgt=hgtlab

      ncrvs=ncrv

      m = mm
      n = nn
      if(ispec.ne.0) then
        ndec=ndeccn
      else
        ndec=1
      endif
      cvs=cv*scal
      xmin=1.e20
      ymin=1.e20
      nsav=0

      if(cvs.eq.0.) then
        nchar=2+ndec
      else
        nchar=max1(1.,alog10(abs(cvs))+1.)+1+ndec
      endif
      if(cvs.lt.0.)nchar=nchar+1

!  Add .05" for space on either side of label
      width=hgt*szrat*nchar+.13

      hhgt=hgt*wide/2.
!# modified: 10/9/2008
      if (ioffp .ne. 0) then
        jump1=35
        jump2=60
      else
        jump1=45
        jump2=65
      endif
!# modified: 10/9/2008
   10 ix0 = ix
      iy0 = iy
      is0 = is
      if (ioffp .eq. 0) go to 15
      ix2 = ix+inx(is)
      iy2 = iy+iny(is)
      ipen = z(ix,iy).ne.spval .and. z(ix2,iy2).ne.spval
      ipeno = ipen
   15 if (idx .ne. 0) then
        y = iy
        isub = ix+idx
        x = cfcn(z(ix,iy),z(isub,iy))*float(idx)+float(ix)
      else
        x = ix
        isub = iy+idy
        y = cfcn(z(ix,iy),z(ix,isub))*float(idy)+float(iy)
      endif
   20 xx1=(fx(x,y)-xor)/xsc
      yy1=(fy(x,y)-yor)/ysc
      if(ipen) then
        if(nsav.eq.0) then
          nsav=nsav+1
          if(nsav.gt.maxpt) then
            print *,'nsav too big in drlin2, program abandoned'
            stop
          endif
          xsav(nsav)=xx1
          ysav(nsav)=yy1
          isav(nsav)=ix
          jsav(nsav)=iy
          ipuord(nsav)=3
          xmin=dmin1(xmin,xx1)
          ymin=dmin1(ymin,yy1)
        else if(dabs(xx1-xsav(nsav)).gt.1.e-14.or. dabs(yy1-ysav(nsav))        &
        .gt.1.e-14) then
          nsav=nsav+1
          if(nsav.gt.maxpt) then
            print *,'nsav too big in drlin2, program abandoned'
            stop
          endif
          xsav(nsav)=xx1
          ysav(nsav)=yy1
          isav(nsav)=ix
          jsav(nsav)=iy
          ipuord(nsav)=3
          xmin=dmin1(xmin,xx1)
          ymin=dmin1(ymin,yy1)
        endif
        if(ilabl.eq.1) cvs=cv*scal
      endif
   25 is = is+1
      if (is .gt. 8) is = is-8
      idx = inx(is)
      idy = iny(is)
      ix2 = ix+idx
      iy2 = iy+idy
      if (iss .ne. 0) go to 30
      if (ix2.gt.m .or. iy2.gt.n .or. ix2.lt. 1.or. iy2.lt.1) go to 80
   30 if ((cv-z(ix2,iy2)).le.0) then
        is = is+4
        ix = ix2
        iy = iy2
        go to 25
      else
        if(mod(is,2).eq.0) goto 25
      endif
!# modified: 10/9/2008
      if (jump1.EQ.35) then     
         isbig = is+(8-is)/6*8
         ix3 = ix+inx(isbig-1)
         iy3 = iy+iny(isbig-1)
         ix4 = ix+inx(isbig-2)
         iy4 = iy+iny(isbig-2)
         ipeno = ipen
         if (iss .EQ. 0) then      
            if (ix3.gt.m .or. iy3.gt.n .or. ix3.lt. 1.or. iy3.lt.1) go to 80
            if (ix4.gt.m .or. iy4.gt.n .or. ix4.lt. 1.or. iy4.lt.1) go to 80
         end if
         ipen = z(ix,iy).ne.spval .and. z(ix2,iy2).ne.spval .and. z        &
               (ix3,iy3).ne.spval .and. z(ix4,iy4).ne.spval
      end if
!# modified: 10/9/2008
      if (idx .eq. 0) go to 50
      y = iy
      isub = ix+idx
      x = cfcn(z(ix,iy),z(isub,iy))*float(idx)+float(ix)
      go to 55
   50 x = ix
      isub = iy+idy
      y = cfcn(z(ix,iy),z(ix,isub))*float(idy)+float(iy)
   55 continue 
!# modified: 10/9/2008
      if (jump2.EQ.65) goto 65
!# modified: 10/9/2008
      if (.not.ipen) go to 70
      if (ipeno) go to 65

! End of line segment

      xx1=(fx(xold,yold)-xor)/xsc
      yy1=(fy(xold,yold)-yor)/ysc
      nsav=nsav+1
      if(nsav.gt.maxpt) then
        print *,'nsav too big in drlin2, program abandoned'
        stop
      endif
      xsav(nsav)=xx1
      ysav(nsav)=yy1
      isav(nsav)=ix
      jsav(nsav)=iy
      ipuord(nsav)=3
      xmin=dmin1(xmin,xx1)
      ymin=dmin1(ymin,yy1)
      if(ilabl.eq.1) cvs=cv*scal

! Continue line segment

   65 xx2=(fx(x,y)-xor)/xsc
      yy2=(fy(x,y)-yor)/ysc
      if(nsav.eq.0) then
        nsav=nsav+1
        if(nsav.gt.maxpt) then
          print *,'nsav too big in drlin2, program abandoned'
          stop
        endif
        xsav(nsav)=xx2
        ysav(nsav)=yy2
        isav(nsav)=ix
        jsav(nsav)=iy
        ipuord(nsav)=2
        xmin=dmin1(xmin,xx2)
        ymin=dmin1(ymin,yy2)
      else if(dabs(xx2-xsav(nsav)).gt.1.e-14.or. dabs(yy2-ysav(nsav))          &
      .gt.1.e-14) then
        nsav=nsav+1
        if(nsav.gt.maxpt) then
          print *,'nsav too big in drlin2, program abandoned'
          stop
        endif
        xsav(nsav)=xx2
        ysav(nsav)=yy2
        isav(nsav)=ix
        jsav(nsav)=iy
        ipuord(nsav)=2
        xmin=dmin1(xmin,xx2)
        ymin=dmin1(ymin,yy2)
      endif
   70 xold = x
      yold = y
      if (is .ne. 1) go to 75
      np = np+1
      if (np .gt. nr) go to 80
      ir(np) = lc16(ix)+iy
   75 if (iss .eq. 0) go to 25
      if (ix.ne.ix0 .or. iy.ne.iy0 .or. is.ne.is0) go to 25

! End of line

   80 continue
!
!  -----------------------------------------------------------------------
!  NEXT SECTION IS MODFIED JIMAR CODE
!
!                 The arrays X, Y are now complete. calculate
!                 the distance S along the contour. Start plotting.
!
      if(ilabb.eq.-999.or.ilabl.ne.1) then
        if(nsav.eq.0) return
        ncs=0
   85   ncs=ncs+1
        if(ncs.gt.nsav) goto 90
        if(ipuord(ncs).eq.3) then
          npts=1
          ncrv=ncrv+1
          rcval(ncrv)=cv
          xx(npts,ncrv)=xsav(ncs)
          yy(npts,ncrv)=ysav(ncs)
          nump(ncrv)=npts
        elseif(ipuord(ncs).eq.2) then
          npts=npts+1
          xx(npts,ncrv)=xsav(ncs)
          yy(npts,ncrv)=ysav(ncs)
          nump(ncrv)=npts
        endif
        goto 85
   90   continue
        if(npts.ne.1) then
          nump(ncrv)=npts
!     **************************************************************
!     ARL Addition
!     **************************************************************
!         if(igis.ge.1)call gisdump(xx(1,ncrv),yy(1,ncrv),npts,cv)
!     **************************************************************
        endif

        RETURN
      ENDIF
      END
!*****DRLINDSH
      subroutine drlindsh (z,l,mm,nn,ilabl,scal)
      parameter(maxpt=10000)
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      common/plt1/cmdstr
      dimension       z(l,nn)
      dimension xsav(maxpt),ysav(maxpt),s(maxpt),ipuord(maxpt)
      dimension isav(maxpt),jsav(maxpt)
      dimension nsava(maxpt),smaxa(maxpt)
      common/contyp/ispcon,idsh
      common/pltparam/curlin
      dimension XX(maxpt),YY(maxpt)
      dimension xdash(maxpt),ydash(maxpt)
!
! this routine traces a contour line when given the beginning by stline.
! transformations can be added by deleting the statement functions for
! fx and fy in drline and minmax and adding external functions.
! x=1. at z(1,j), x=float(m) at z(m,j). x takes on non-integer values.
! y=1. at z(i,1), y=float(n) at z(i,n). y takes on non-integer values.
!  this routine is same as drline, but labels contours along the contour lines
!  rather than at the edges, using a modified version of code obtained from
!  jimar(contur.for).
!
      common /conre2/ ix ,iy ,idx ,idy , is ,iss ,np ,cv , inx(8) ,iny         &
      (8) ,ir(10000) ,nr

      common/conre1/ioffp,spval
      common/scales/xor,xsc,yor,ysc

      common/conpar/ispec,ioffpp,spvall,ilegg,ilabb,nhii,ndeccn,nlbll,         &
      lscal,ldash,hgtlab

      common/konval/conval(100),numcon
      common/kkplot/szrat

!     **************************************************************
!     ARL Addition
!     **************************************************************
      common/gis/gisout,gisatt,gistmp,igis
      character*80 gisout,gisatt,gistmp
!     **************************************************************

      logical         ipen       ,ipeno
      data ipen,ipeno/.true.,.true./

!  hgt is height of label
!  slab1f is fraction of total contour length to skip before first label on
!  contour
!  dslab is distance between labels on same contour
!
      data hgt/.1/,dslab/3.5/,slab1f/.1/,wide/1./,big/1.e30/  !.857=wide
!
      fx(x,y) = x
      fy(x,y) = y
      lc16(k) = k*65536
!     lc16(k) = k*'200000'O   !VMS
!     lc16(k) = k*8#200000    !MS
      cfcn(p1,p2) = (p1-cv)/(p1-p2)
      if(hgtlab.ne.0.)hgt=hgtlab

!      pi=4.*abs(atan(1.)) ! no need for transcendental function to find pi
      m = mm
      n = nn
      if(ispec.ne.0) then
        ndec=ndeccn
      else
        ndec=1
      endif
      cvs=cv*scal
      xmin=1.e20
      ymin=1.e20
      nsav=0

      if(cvs.eq.0.) then
        nchar=2+ndec
      else
        nchar=max1(1.,alog10(abs(cvs))+1.)+1+ndec
      endif
      if(cvs.lt.0.)nchar=nchar+1

!  Add .05" for space on either side of label
      width=hgt*szrat*nchar+.13

      hhgt=hgt*wide/2.
!# modified: 10/9/2008
      if (ioffp .ne. 0) then
        jump1=35
        jump2=60
      else
        jump1=45
        jump2=65
      endif
!# modified: 10/9/2008
   10 ix0 = ix
      iy0 = iy
      is0 = is
      if (ioffp .eq. 0) go to 15
      ix2 = ix+inx(is)
      iy2 = iy+iny(is)
      ipen = z(ix,iy).ne.spval .and. z(ix2,iy2).ne.spval
      ipeno = ipen
   15 if (idx .ne. 0) then
        y = iy
        isub = ix+idx
        x = cfcn(z(ix,iy),z(isub,iy))*float(idx)+float(ix)
      else
        x = ix
        isub = iy+idy
        y = cfcn(z(ix,iy),z(ix,isub))*float(idy)+float(iy)
      endif
   20 xx1=(fx(x,y)-xor)/xsc
      yy1=(fy(x,y)-yor)/ysc
      if(ipen) then
        if(nsav.eq.0) then
          nsav=nsav+1
          if(nsav.gt.maxpt) then
            print *,'nsav too big in drlindsh, program abandoned'
            stop
          endif
          xsav(nsav)=xx1
          ysav(nsav)=yy1
          isav(nsav)=ix
          jsav(nsav)=iy
          ipuord(nsav)=3
          xmin=amin1(xmin,xx1)
          ymin=amin1(ymin,yy1)
        else if(abs(xx1-xsav(nsav)).gt.1.e-4.or. abs(yy1-ysav(nsav)).gt.       &
        1.e-4) then
          nsav=nsav+1
          if(nsav.gt.maxpt) then
            print *,'nsav too big in drlindsh, program abandoned'
            stop
          endif
          xsav(nsav)=xx1
          ysav(nsav)=yy1
          isav(nsav)=ix
          jsav(nsav)=iy
          ipuord(nsav)=3
          xmin=amin1(xmin,xx1)
          ymin=amin1(ymin,yy1)
        endif
        if(ilabl.eq.1) cvs=cv*scal
      endif
   25 is = is+1
      if (is .gt. 8) is = is-8
      idx = inx(is)
      idy = iny(is)
      ix2 = ix+idx
      iy2 = iy+idy
      if (iss .ne. 0) go to 30
      if (ix2.gt.m .or. iy2.gt.n .or. ix2.lt. 1.or. iy2.lt.1) go to 80
   30 if ((cv-z(ix2,iy2)).le.0) then
        is = is+4
        ix = ix2
        iy = iy2
        go to 25
      else
        if(mod(is,2).eq.0) goto 25
      endif
!# modified: 10/9/2008
      if (jump1.EQ.35) then 
         isbig = is+(8-is)/6*8
         ix3 = ix+inx(isbig-1)
         iy3 = iy+iny(isbig-1)
         ix4 = ix+inx(isbig-2)
         iy4 = iy+iny(isbig-2)
         ipeno = ipen
         if (iss .EQ. 0) then     
            if (ix3.gt.m .or. iy3.gt.n .or. ix3.lt. 1.or. iy3.lt.1) go to 80
            if (ix4.gt.m .or. iy4.gt.n .or. ix4.lt. 1.or. iy4.lt.1) go to 80
         end if
         ipen = z(ix,iy).ne.spval .and. z(ix2,iy2).ne.spval .and. z        &
               (ix3,iy3).ne.spval .and. z(ix4,iy4).ne.spval
      end if
!# modified: 10/9/2008
      if (idx .eq. 0) go to 50
      y = iy
      isub = ix+idx
      x = cfcn(z(ix,iy),z(isub,iy))*float(idx)+float(ix)
      go to 55
   50 x = ix
      isub = iy+idy
      y = cfcn(z(ix,iy),z(ix,isub))*float(idy)+float(iy)
   55 continue
!# modified: 10/9/2008
      if (jump2.EQ.65) go to 65
!# modified: 10/9/2008
      if (.not.ipen) go to 70
      if (ipeno) go to 65

! End of line segment

      xx1=(fx(xold,yold)-xor)/xsc
      yy1=(fy(xold,yold)-yor)/ysc
      nsav=nsav+1
      if(nsav.gt.maxpt) then
        print *,'nsav too big in drlindsh, program abandoned'
        stop
      endif
      xsav(nsav)=xx1
      ysav(nsav)=yy1
      isav(nsav)=ix
      jsav(nsav)=iy
      ipuord(nsav)=3
      xmin=amin1(xmin,xx1)
      ymin=amin1(ymin,yy1)
      if(ilabl.eq.1) cvs=cv*scal

! Continue line segment

   65 xx2=(fx(x,y)-xor)/xsc
      yy2=(fy(x,y)-yor)/ysc
      if(nsav.eq.0) then
        nsav=nsav+1
        if(nsav.gt.maxpt) then
          print *,'nsav too big in drlindsh, program abandoned'
          stop
        endif
        xsav(nsav)=xx2
        ysav(nsav)=yy2
        isav(nsav)=ix
        jsav(nsav)=iy
        ipuord(nsav)=2
        xmin=amin1(xmin,xx2)
        ymin=amin1(ymin,yy2)
      elseif(abs(xx2-xsav(nsav)).gt.1.e-4.or.abs(yy2-ysav(nsav)).gt.1.e-4)then
        nsav=nsav+1
        if(nsav.gt.maxpt) then
          print *,'nsav too big in drlindsh, program abandoned'
          stop
        endif
        xsav(nsav)=xx2
        ysav(nsav)=yy2
        isav(nsav)=ix
        jsav(nsav)=iy
        ipuord(nsav)=2
        xmin=amin1(xmin,xx2)
        ymin=amin1(ymin,yy2)
      endif
   70 xold = x
      yold = y
      if (is .ne. 1) go to 75
      np = np+1
      if (np .gt. nr) go to 80
      ir(np) = lc16(ix)+iy
   75 if (iss .eq. 0) go to 25
      if (ix.ne.ix0 .or. iy.ne.iy0 .or. is.ne.is0) go to 25

! End of line

   80 continue

!  --------------------------------------------------------------------
!  NEXT SECTION IS MODFIED JIMAR CODE
!
!                 The arrays X, Y are now complete. calculate
!                 the distance S along the contour. Start plotting.
!
      if(ilabl.ne.1) then
        if(nsav.eq.0) return
        ncs=0
        icrv=0
   85   ncs=ncs+1
        if(ncs.gt.nsav) goto 90
        if(ipuord(ncs).eq.3) then
          icrv=icrv+1
          if(icrv.ne.1) call dshcrv(xx,yy,npts,IDSH,0.)
          npts=1
          xx(npts)=xsav(ncs)
          yy(npts)=ysav(ncs)
        elseif(ipuord(ncs).eq.2) then
          npts=npts+1
          xx(npts)=xsav(ncs)
          yy(npts)=ysav(ncs)
        endif
        goto 85
   90   continue
!       **************************************************************
!       ARL Addition
!       **************************************************************
        if(npts.ne.1.and.igis.ge.1)call gisdump(xx,yy,npts,cv)
!       **************************************************************
        if(npts.ne.1) call dshcrv(xx,yy,npts,IDSH,0.)

        return
      endif
!------------------------------------------------------------------------------
!  Find out what contour level we're on
      do nc=1,numcon
        if(abs(cv-conval(nc)).lt.1.e-4) then
          lev=nc
          goto 95
        endif
      enddo
      print *,'cannot find current contour level, program abandoned'
      stop
   95 continue

      nsav0=nsav
      s(1)=0.
      if(nsav .le. 1) return

!  Calculate s(k)'s normally
      smaxx=-1.e20
      do k=2,nsav0
        km1 = k - 1
        kp1 = k + 1
        dxx = xsav(k)-xsav(km1)
        dyy = ysav(k)-ysav(km1)
        s(k) = s(km1) + sqrt(dxx*dxx+dyy*dyy)
      enddo

      do k=2,nsav0
        nsava(k)=nsav0
        smaxa(k)=s(nsav0)
      enddo

!  Check for non-continuous contours and reset nsava,smaxa if necessary
      k1=1
      do k=2,nsav0
        km1 = k - 1
        kp1 = k + 1
        if(ipuord(k).eq.2.and.ipuord(kp1).eq.3) then
          do kk=k1,k
            nsava(kk)=k
            smaxa(kk)=s(k)
          enddo

          k1=kp1
        endif
      enddo

      slab1 = smaxa(1)*slab1f !fraction of total length
      stest = dslab - slab1 !set so first label is at slab1
      k = 1
      ndsh=1
      xdash(ndsh)=xsav(k)
      ydash(ndsh)=ysav(k)
      call plot(xsav(k),ysav(k),ipuord(k))

!  Check conditions for labelling.

  100 continue                    ! k loop
      nsav=nsava(k)
      smax = smaxa(k)

      if(ndec.le.-2)go to 125
!  Don't label if not enough contour left
      if(smax-s(k).le.width)go to 125
      km1= max0(k-1,1)
      stest = stest + s(k)-s(km1)
      if(stest.lt.dslab)go to 125
      kp1=k+1
      if(lev.eq.1) go to 105
      dlev=abs(conval(lev)-conval(lev-1))

!  Is there enough space between adjacent contours?

      i=(xsav(k)-xmin)*xsc+1.01
      i=min0(i,mm-1)
      j=(ysav(k)-ymin)*ysc+1.01
      j=min0(j,nn-1)

      i=isav(k)
      j=jsav(k)
      i=min0(i,mm-1)
      j=min0(j,nn-1)

      dzdx=(z(i+1,j)-z(i,j))*xsc
      if( dzdx .ge. big ) go to 125
      dzdy=(z(i,j+1)-z(i,j))*ysc
      if( dzdy .ge. big ) go to 125
      dzdg=sqrt(dzdx*dzdx+dzdy*dzdy)
      if(dzdg.eq.0.) go to 105
      cspace=dlev/dzdg
      if(cspace.lt.hgt/2.) go to 125 !label drawn at mid-contour

  105 continue

!   Check if enough space for contour.
!   If so, check that curvature is not too great.
!   Determine # of points to drop to fit label.

      do 110 kk = kp1,nsav
        dxx = xsav(kk)-xsav(k)
        dyy = ysav(kk)-ysav(k)
        space = sqrt(dxx*dxx+dyy*dyy)
        if(space.lt.width)go to 110
        ark=s(kk)-s(k)
        if(ark.eq.0.0) goto 125
        if(space/ark.lt.0.80)go to 125
        go to 115
  110 continue
      go to 125

!  Draw the label
  115 stest=0.
      kkm1 = kk - 1

!  Determine how much of the available space is actually needed.

      do kkk=1,5
        t=0.2*kkk
        dxx=xsav(kkm1)+(xsav(kk)-xsav(kkm1))*t-xsav(k)
        dyy=ysav(kkm1)+(ysav(kk)-ysav(kkm1))*t-ysav(k)
        space=sqrt(dxx*dxx+dyy*dyy)
        if(space.gt.width) goto 120
      enddo

  120 xendl = xsav(k) + width*dxx/space
      yendl = ysav(k) + width*dyy/space

      if(dxx.eq.0) then
        angle = 90.
        if(dyy.lt.0.)angle=270.
      else
        angle = atan(dyy/dxx)/rdpdeg
      endif

      if(dxx.ge.0.0)then
        xlab = xsav(k)+ hhgt*(dxx+dyy)/space
        ylab = ysav(k)+ hhgt*(dyy-dxx)/space
      else
        xlab = xendl-hhgt*(dxx+dyy)/space
        ylab = yendl-hhgt*(dyy-dxx)/space
      endif

      call keknum(xlab,ylab,hgt,cvs,angle,ndec,0)
      call dshcrv(xdash,ydash,ndsh,idsh,0.)
      ndsh=1
      xdash(ndsh)=xendl
      ydash(ndsh)=yendl
      ndsh=ndsh+1
      xdash(ndsh)=xsav(kk)
      ydash(ndsh)=ysav(kk)
      k=kk

!  Plot the segment from xk,yk to  xsav(k+1),ysav(k+1).

  125 continue
      if((k+1).le.nsav0) then
        if(ipuord(k+1).eq.3) then
          call dshcrv(xdash,ydash,ndsh,idsh,0.)
          ndsh=1
          xdash(ndsh)=xsav(k+1)
          ydash(ndsh)=ysav(k+1)
        else
          ndsh=ndsh+1
          xdash(ndsh)=xsav(k+1)
          ydash(ndsh)=ysav(k+1)
        endif
      endif

      k=k+1
      if( k .lt. nsav0 ) go to 100
      if(ndsh.ne.1) call dshcrv(xdash,ydash,ndsh,idsh,0.)

      return
      end
!*****DRWCRV
      subroutine drwcrv(xarr,yarr,n,thkk,closer)
      logical closer
      character*132 cmdstr
      common/plt1/cmdstr
      dimension xarr(n),yarr(n)

! Curwid is set in routine setlw
      common/lcom/curwid

      cmdstr='Np'
      call filler

      if(thkk.eq.0.) then      !Use the current linewidth
      else
        curwids=curwid
        call setlw(thkk)
      endif
      call movet(xarr(1),yarr(1))
      do 10 nn=2,n
        call linet(xarr(nn),yarr(nn))
        if(mod(nn,100).eq.0) then !Stroke now to avoid possible limitcheck
          call stroke
          call movet(xarr(nn),yarr(nn))
        endif
   10 continue

      if(closer)then
        cmdstr='Cs'
      else
        cmdstr='S'
      endif
      call filler
      if(thkk.ne.0.) call setlw(curwids)   !Reset linewidth
      return
      end
!*****DRWTRI
      subroutine drwtri(xc,yc,side,thkk)
      dimension xarr(3),yarr(3)
      common/lcom/curwid

      sqrt2=sqrt(2.)
      sqrt3=sqrt(3.)
      y1=yc-side/(2.*sqrt3)
      x1=xc-side/2.
      x2=xc
      y2=yc+side/sqrt3
      x3=xc+side/2.
      y3=y1
      if(thkk.eq.0.) then    !Use the current linewidth
        thk=curwid
      else
        thk=thkk
      endif
      if(thk.eq.0.) then
        call plot(x1,y1,3)
        call plot(x2,y2,2)
        call plot(x3,y3,2)
        call plot(x1,y1,2)
      else
        xarr(1)=x1
        xarr(2)=x2
        xarr(3)=x3
        yarr(1)=y1
        yarr(2)=y2
        yarr(3)=y3
        call drwcrv(xarr,yarr,3,thk,.true.)
      endif
      return
      end
!*****DSHCRV
      subroutine dshcrv(x,y,npts,idshpn,thk)
!  this routine draws a dashed curve connecting points x(n),yn(n),n=1,npts
!  x,y arrays are in real inches
!  idshpn is the number of solid/blank patterns per inch
!  thk is the thickness of the dashes(inches).  if zero, a single line is drawn
!  for the dashes
      parameter (pi=3.14159265358979)
      double precision dlen,dleft,xsd,ysd,xed,yed
      dimension x(npts),y(npts)
      common/lcom/curwid

!      pi=4.*abs(atan(1.)) ! no need to use transcendental function to find pi

      iset=0
      if(thk.eq.0.) then       ! Use current linewidth
      elseif(thk.ne.curwid) then
        iset=1
        curwids=curwid
        call setlw(thk)
      endif

      thko2=thk/2.
      dlen=1./(2.*idshpn)
      ipen=3
      call plot(x(1),y(1),ipen)
      xsd=x(1)
      ysd=y(1)
      n=2

   10 continue
      if(n.gt.npts) goto 15         ! we're done
      ipen=5-ipen !(3-->2, 2-->3)

!  Find angle to next data point
      absx=abs(x(n)-xsd)
      if(absx.lt.1.e-5) then
        ang=pi/2.*dsign(dble(1.),y(n)-ysd)
      else
        ang=atan2( (y(n)-ysd),(x(n)-xsd) )
      endif

!  Calculate end point of dash or space
      xed=xsd+dlen*cos(ang)
      yed=ysd+dlen*sin(ang)
      xe=xed
      ye=yed

!  See if we can fit this in before next data point
      dist=sqrt((x(n)-xsd)**2+(y(n)-ysd)**2)
      if(dist.le.dlen)  then  !can't fit
        dleft=dlen
        do nn=n,npts
!        Find angle to next data point
        absx=abs(x(nn)-xsd)
        if(absx.lt.1.e-5) then
          ang=pi/2.*dsign(dble(1.),y(nn)-ysd)
        else
          ang=atan2( (y(nn)-ysd),(x(nn)-xsd) )
        endif
        dist=sqrt((x(nn)-xsd)**2+(y(nn)-ysd)**2)
        if(dist.le.dleft) then    !Draw line to next data pt & calc rest
          call plot(x(nn),y(nn),ipen)
          xsd=x(nn)
          ysd=y(nn)
          dleft=dleft-dist               !Remaining portion
          if(dleft.eq.0.) then           !Exact fit, so just do next point
            n=nn+1
            goto 10
          endif
        else                             !Draw remaining portion and move toward
          xnew=xsd+dleft*cos(ang)
          ynew=ysd+dleft*sin(ang)
          call plot(xnew,ynew,ipen)
          xsd=xnew
          ysd=ynew
          n=nn
          goto 10
        endif
      enddo
      else                               !Can fit, so just draw normally
!       Check to see if we should just move to data point itself since roundoff
!       could cause weird angles
      dist=sqrt((x(n)-xed)**2+(y(n)-yed)**2)
      if(dist.le.1.e-5) then
        xed=x(n)
        yed=y(n)
        xe=xed
        ye=yed
        n=n+1
      endif
      call plot(xe,ye,ipen)
      xsd=xed
      ysd=yed
      goto 10
      endif

   15 continue
      if(iset.eq.1) call setlw(curwids) !reset linewidth to original
      return
      end
!*****DSHLIN
      subroutine dshlin(x1,y1,x2,y2,idshpn,thk)
!  This routine draws a dashed line between points (x1,y1) and (x2,y2).
!  idshpn is the number of solid/blank patterns per inch
!  thk is the thickness of the dashes(inches).  if zero, a single line is drawn
!  for the dashes
      parameter (pi=3.14159265358979)
      common/lcom/curwid

!      pi=4.*abs(atan(1.)) ! no need for transcendental function to find pi
      iset=0
      if(thk.eq.0.) then             !Use current width
      elseif(thk.ne.curwid) then
        iset=1
        curwids=curwid
        call setlw(thk)
      endif

      dlen=1./(2.*idshpn)
      tlen=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
      np=tlen*float(idshpn)
      if(abs(x2-x1).lt.1.e-5) then
        ang=pi/2.*sign(1.,y2-y1)
      else
        ang=atan2( (y2-y1),(x2-x1) )
      endif
!  do np patterns
      tdist=0.
      xs=x1
      ys=y1
      do n=1,np
        xe=xs+dlen*cos(ang)
        ye=ys+dlen*sin(ang)
        call plot(xs,ys,3)
        call plot(xe,ye,2)
        tdist=tdist+2.*dlen
        xs=xe+dlen*cos(ang)
        ys=ye+dlen*sin(ang)
      enddo

!  Finish partial pattern
      call plot(xs,ys,3)
      call plot(x2,y2,2)
      if(iset.eq.1) call setlw(curwids) !Reset linewidth to original

      return
      end
!*****FACTOR
      subroutine factor(facc)
      common/plt2/fac
      character*132 cmdstr
      common/plt1/cmdstr

!  Unscale previous scaling
      recipx=1./fac
      recipy=1./fac
      write(cmdstr,'(2f7.3,'' scale'')')recipx,recipy
      call filler
      fac=facc
      write(cmdstr,'(2F7.3,'' scale'')')fac,fac
      call filler
      return
      end
!*****FAROHED
      subroutine farohed(xpp,ypp,dir,arolnp,sprang,locxy,fill)
!  This routine is similar to arohed, but style of arrowhead is fancier
!  xpp,ypp are coordinates of tip of arrowhead
!  dir is angle, in degrees, of arrowhead, measured east from north
!  arolnp is the length, in inches, of the sides of the arrowhead
!  sprang is the angle, in degrees, from one side of the arrowhead
!  to the arrow, that is half the angular spread of the arrowhead
!  locxy=0 x,y at arrow point
!  locxy=1 x,y at center of arrow
!  locxy=2 x,y,at tail of arrow
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      logical fill
      dimension ymul(3)
      dimension x(4),y(4)
      equivalence(x(1),x1),(x(2),x2),(x(3),x3),(x(4),x4)
      equivalence(y(1),y1),(y(2),y2),(y(3),y3),(y(4),y4)
      character*132 cmdstr
      character*80 ifrmt
      common/plt1/cmdstr
      common/cnvcom/conver
      data ymul/0.,.5,1./

!      pi=4*abs(atan(1.)) ! no need to use transcendental function to find pi
!      rdpdeg=pi/180.

      xp=xpp
      yp=ypp
      ya=dir*rdpdeg
      cosdir=cos(ya)
      sindir=sin(ya)
      arolen=arolnp
      sprrad=sprang*rdpdeg
      xa=sin(sprrad)*arolen
      ya=cos(sprrad)*arolen
      i=locxy
      x1=0.
      y1=ya*ymul(i+1)
      x2=xa
      y2=y1-ya
      x3=-xa
      y3=y2
      x4=x1         !x center of arc
      rfac=1.1
      y4=-rfac*ya   !y center of arc

!  Now rotate
      do i=1,4
        xa=x(i)
        ya=y(i)
        x(i)=((xa*cosdir)+(ya*sindir))+xp
        y(i)=((ya*cosdir)-(xa*sindir))+yp
      enddo

      call movet(x(3),y(3))
      call linet(x(1),y(1))
      call linet(x(2),y(2))

!  Define some arc stuff
      rad=sqrt((x2-x4)**2+(y2-y4)**2)
      ang1=atan2((y2-y4),(x2-x4))/rdpdeg
      ang2=atan2((y3-y4),(x3-x4))/rdpdeg

      ifrmt='(f8.2,'' '',f8.2,'' '',f8.2,'' '',2f8.2,'//                       &
      ''' arc closepath'')'
      cmdstr=' '
      write(cmdstr,ifrmt)x4*conver,y4*conver,rad*conver,ang1,ang2
      call filler

      if(fill) then
        cmdstr='fill'
      else
        cmdstr='stroke'
      endif
      call filler

      return
      end
!*****FILL
      subroutine fill
      character*132 cmdstr
      common/plt1/cmdstr

      cmdstr='fill'
      call filler
      return
      end
!*****FILLBOX
      subroutine fillbox(xpts,ypts,npts,grylev)
!  This routine fills region bounded by arrays xpts and ypts
!  with gray level grylev.  xpts,ypts are in inches.  Uses first npts
!  points of xpts and ypts
      dimension xpts(npts),ypts(npts)
      character*132 cmdstr,scr
      common/plt1/cmdstr
      common/cnvcom/conver
      common/colrcom/cred,cgreen,cblue,cgry

      if(grylev.ne.cgry) call setgry(grylev)

      cmdstr=' '
      do 10 nn=1,npts
        xx=xpts(nn)*conver
        yy=ypts(nn)*conver
        scr=' '
        write(scr,'('' '',f8.2,'' '',f8.2)')xx,yy
        lens=lenstr(scr,132)
        if(nn.eq.1) then
          cmdstr=scr(2:lens)   !Don't need initial space
        else
          lc=lenstr(cmdstr,132)
          if((lc+lens).gt.132) then
            call filler
            cmdstr=scr(1:lens)
          else
            cmdstr=cmdstr(1:lc)//scr(1:lens)
          endif
        endif
   10 continue

      nm1=npts-1
      scr=' '
      write(scr,'(i6,'' Fbn'')')nm1
      lc=lenstr(cmdstr,132)
      ls=lenstr(scr,132)
      if((lc+ls).gt.132) then
        call filler
        cmdstr=scr(1:ls)
      else
        cmdstr=cmdstr(1:lc)//scr(1:ls)
      endif
      call filler
      return
      end
!*****FILLBOXC
      subroutine fillboxc(xpts,ypts,npts,IOFF,JOFF,GRYLEV)
!  This routine fills region bounded by arrays xpts and ypts
!  with gray level grylev.
!  xpts,ypts are in inches.  Uses first npts points of xpts and ypts
!  Called by confill.

      dimension xpts(npts),ypts(npts)
      character*132 cmdstr,scr
      common/plt1/cmdstr
      common/cnvcom/conver
      common/colrcom/cred,cgreen,cblue,cgry
      common/boxcom/dx,dy

      IF(grylev.ne.cgry)call setgry(grylev)

      cmdstr=' '
      do nn=1,npts
        scr=' '
        xp=(ioff+xpts(nn))*dx*conver
        yp=(joff+ypts(nn))*dy*conver
        write(scr,'('' '',f8.2,'' '',f8.2)')xp,yp
        lens=lenstr(scr,132)
        if(nn.eq.1) then
          cmdstr=scr(2:lens)    !Don't need initial space
        else
          lc=lenstr(cmdstr,132)
          if((lc+lens).gt.132) then
            call filler
            cmdstr=scr(1:lens)
          else
            cmdstr=cmdstr(1:lc)//scr(1:lens)
          endif
        endif
      enddo

      nm1=npts-1
      scr=' '
      write(scr,'(i6,'' Fbnc'')')nm1
      lc=lenstr(cmdstr,132)
      ls=lenstr(scr,132)
      if((lc+ls).gt.132) then
        call filler
        cmdstr=scr(1:ls)
      else
        cmdstr=cmdstr(1:lc)//scr(1:ls)
      endif
      call filler
      return
      end
!*****FILLER
      subroutine filler
!  nfild is the last position filled in the compressed aaa buffer array
!  work is a work array used to load array aaa
      character*132 cmdstr,cmdc(132)*1
      common/plt1/cmdstr
      common/outcom/iunit
      logical ispace
!     equivalence (cmdstr,cmdc(1))
      ibslash=92
      lc=0
      lcc=lenstr(cmdstr,132)
      ispace=.false.

!  itot is running total of left/right parentheses in text string
!  if itot=0 then we are not in text mode, i.e. left=right

      itot=0
      icclst=-999
      do 10 l=1,lcc
        icc=ichar(cmdstr(l:l))
        if(icc.eq.32.and.ispace.and.itot.eq.0) goto 10 !Don't place 2 or more
!                                                       spaces together if
!                                                       not in text mode
        if(icc.ge.32.and.icc.le.127) then
          lc=lc+1
          cmdc(lc)=cmdstr(l:l)
        endif
        if(icc.eq.32) then
          ispace=.true.
        else
          ispace=.false.
        endif
        if(icc.eq.40.and.icclst.ne.ibslash)itot=itot+1
        if(icc.eq.41.and.icclst.ne.ibslash)itot=itot-1
        icclst=icc
   10 continue

!  Write cmdstr to output file
      write(iunit,'(132a1)')(cmdc(ii),ii=1,lc)
      return
      end
!*****FILRGN
      subroutine filrgn(xpts,ypts,npts,grylev)
! this routine fills region bounded by arrays xpts and ypts
!  with gray level grylev.  xpts,ypts are in inches
      dimension xpts(npts),ypts(npts)
      common/colrcom/cred,cgreen,cblue,cgry

      cgrysv=cgry   !Save initial gray value
      call fillbox(xpts,ypts,npts,grylev)
      if(grylev.ne.cgrysv) call setgry(cgrysv) !reset to original
      return
      end
!*****FILRGNC
      subroutine filrgnc(xpts,ypts,npts)
! this routine fills region bounded by arrays xpts and ypts
!  with current gray level or colors.  xpts,ypts are in inches

      dimension xpts(npts),ypts(npts)
      character*132 cmdstr
      common/plt1/cmdstr

      do i=1,npts
        if(i.eq.1) then
          call movet(xpts(i),ypts(i))
        else
          call linet(xpts(i),ypts(i))
        endif
      enddo
      cmdstr='closepath fill'
      call filler
      return
      end
!*****GREST
      subroutine grest
!  this routine calls the postscript operator grestore
!  it is useful after clipping
      character*132 cmdstr,curfnt
      common/plt1/cmdstr
      common/fntcom/curfnt,ifntsz,nfont
      common/savcom/ifnt,nfnt

      cmdstr='grestore'
      call filler
      ifntsz=ifnt
      call setfnt(nfnt)
      return
      end
!*****GRKSYM
      subroutine grksym(xp,yp,size,ititle,ang,nchar,mjus)
!  this routine uses postscript symbol font to generate greek symbols
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*80 grkchar,charout*1
      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

      data grkchar/'ABGDEZHQIKLMNXOPRSTUFCYWabgdezhJiklmnxoprstujcywqf'/

      bslash=char(92)
!      pi=4.*abs(atan(1.)) ! no need for transcendental function to find pi

!  Save current font
      nfsav=nfont

!  choose proper character height
      mchar=iabs(nchar)

!  set character size
      iht=size*conver/.6     !.6 FACTOR IS EMPIRICAL

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')')       &
      IHT
      call filler

      charout=grkchar(ititle:ititle)
      if(charout.eq.'?') then
        print *,'Greek character ',ititle,' is not available ',                &
        'in postscript'
        print *,'a blank is being substituted'
        charout=' '
      endif

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size

!  character space height is 2.0 x char height
!  character space width is 1.5 x char width
!  actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) THEN
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        NJUS=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'Incorrect justification code ',i5,'found in ',                &
      'grksym, zero used')
        njus=0
      endif
      shift=njus*strlen/2.

      if(charout.ne.'U') then
        cmdstr='('//charout//') Lend'
      else
        cmdstr='('//bslash//'241) Lend'
      endif

      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      if(xarg.ne.0.) then
        scrc=' '
        write(scrc,'(f5.3,'' Xposjd'')')xarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr               &
        (scrc,132))
      endif

      if(yarg.ne.0.) then
        scrc=' '
        write(scrc,'(f5.3,'' Yposjd'')')yarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr               &
        (scrc,132))
      endif

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        if(xarg.ne.0.) then
          scrc=' '
          write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr             &
          (scrc,132))
        endif

        if(yarg.ne.0.) then
          scrc=' '
          write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr             &
          (scrc,132))
        endif
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr               &
        (scrc,132))
      endif

      if(charout.ne.'U') then
        scrc='('//charout//') show'
      else
        scrc='('//bslash//'241) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' ' //scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg=arg / rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

!  Reset current font
      call setfnt(nfsav)

      return
      end
!*****GSAV
      subroutine gsav
!  this routine calls the postscript operator gsave
!  it is useful before clipping
      character*132 cmdstr,curfnt
      common/plt1/cmdstr
      common/fntcom/curfnt,ifntsz,nfont
      common/savcom/ifnt,nfnt
      ifnt=ifntsz
      nfnt=nfont
      cmdstr='gsave'
      call filler
      return
      end
!*****HILITEC
      subroutine hilitec(xp,yp,size,titlein,ang,edg,jusx,jusy,fred,            &
      fgreen,fblue,bred,bgreen,bblue)
!  xp,yp specifies the coordinates of the box, subject to justification
!  given by jusx,jusy
!  this routine highlights the input text with a rectangle of color
!  bred,bgreen,bblue
!  the text is drawn with color fred,fgreen,fblue
!  edg is the fraction of text height to use as an edge border
!  jusx,y is the justification in the x and y directions, respectively.
!  it accepts character strings instead of integer (hollerith) strings
!  this routine does not support octal codes for special characters
!  (use keksym instead)

!  Explanation of coordinates:
!  If the text were drawn at (0.,0.), then the lower left corner of
!  the text bounding box would be at (llx,lly), and the lower left
!  corner of the shaded box would be at (llx-edg,lly-edg)
!  This programs shifts the local coordinate system such that the
!  lower left corner of the shaded box is at (0.,0.), taking into
!  account jusx and jusy
!  Therefore, in the new coordinate system, the lower left corner
!  of the text bounding box is at (edg,edg) and the lower left
!  corner of the text is at (edg-llx,edg-lly).

      character*132 cmdstr
      character*(*) titlein
      character*80 scr
      common/plt1/cmdstr
      common/cnvcom/conver
      common/colrcom/cred,cgreen,cblue,cgry
      data itime/0/
      itime=itime+1

!  Save initial colors
      redsv=cred
      greensv=cgreen
      bluesv=cblue

      lt=len(titlein)

      if(itime.eq.1) then
        cmdstr='/Fwidth { newpath 0 0 moveto text false charpath '             &
             //'flattenpath pathbbox'
        call filler
        cmdstr='/ury exch def /urx exch def /lly exch def /llx exch def'
        call filler
        cmdstr='/swidth urx llx sub def /sheight ury lly sub def} def'
        call filler
      endif

!  Get bounding box before rotation, since it is inaccurate if rotation
!  is not a multiple of 90 degrees (PS ref, p.461)

!  This call sets the font size
      ihgt=size*conver/.6    !Must match keksym, keksymc, etc.
      write(cmdstr,'(i3,'' Setf'')')ihgt
      call filler
      cmdstr='/text ('//titlein(1:lt)//') def'
      call filler

!  This defines llx, lly, swidth, etc.
      cmdstr='gsave Fwidth grestore'    !This eliminates charpath from
      call filler                       !stroking the character outline

!  Define edge width (a fraction of bounding box height)
      write(scr,'(f5.2)')edg
      call blkstp(scr,80,scr,lsc)
      cmdstr='/edg '//scr(1:lsc)//' sheight mul def'
      call filler
      cmdstr='/tothgt sheight 2 edg mul add def' //                            &
      ' /totwid swidth 2 edg mul add def'
      call filler

!  Shift origin. Remember, xp,yp represent the coords of the box,
!  not the text
      call plot(xp,yp,-3)

!  Set angle.  Then everything is in relation to text coordinate frame
      if(ang.ne.0.) call rotate(ang)

!  Figure out where to put the lower left corner of the box ......
      write(cmdstr,'(''/xsh '',i1,'' totwid mul 2 div def'')')jusx
      call filler
      write(cmdstr,'(''/ysh '',i1,'' tothgt mul 2 div def'')')jusy
      call filler

!  ......and move the new origin there
      cmdstr='xsh neg ysh neg translate'
      call filler

!  Draw and fill in the background box (lower left corner is at (0,0))

!  Set background (fill) color
      call setcolr(bred,bgreen,bblue)

      cmdstr= 'Np 0 0 moveto totwid 0 rlineto'
      call filler
      cmdstr='0 tothgt rlineto totwid neg 0 rlineto cf'
      call filler

!  Set foreground (text) color
      call setcolr(fred,fgreen,fblue)

!  Text now starts at (edg-llx,edg-lly)
      cmdstr='edg llx sub edg lly sub moveto text show'
      call filler

!  Reset color to original
      call setcolr(redsv,greensv,bluesv)

!  Unshift origin
      cmdstr='xsh ysh translate'
      call filler

!  Reset angle
      if(ang.ne.0.) call rotate(-ang)

!  Reset origin
      call plot(-xp,-yp,-3)
      return
      end
!*****HILITEG
      subroutine hiliteg(xp,yp,size,titlein,ang,edg,jusx,jusy,fgry,            &
      bgry)
!  xp,yp specifies the coordinates of the box, subject to justification
!  given by jusx,jusy
!  this routine highlights the input text with a rectangle of graylevel bgry.
!  the text is drawn with graylevel fgry.
!  edg is the fraction of text height to use as an edge border
!  jusx,y is the justification in the x and y directions, respectively.
!  it accepts character strings instead of integer (hollerith) strings
!  this routine does not support octal codes for special characters
!  (use keksym instead)

!  Explanation of coordinates:
!  If the text were drawn at (0.,0.), then the lower left corner of
!  the text bounding box would be at (llx,lly), and the lower left
!  corner of the shaded box would be at (llx-edg,lly-edg)
!  This programs shifts the local coordinate system such that the
!  lower left corner of the shaded box is at (0.,0.), taking into
!  account jusx and jusy
!  Therefore, in the new coordinate system, the lower left corner
!  of the text bounding box is at (edg,edg) and the lower left
!  corner of the text is at (edg-llx,edg-lly).

      character*132 cmdstr
      character*(*) titlein
      character*80 scr
      common/plt1/cmdstr
      common/cnvcom/conver
      common/colrcom/cred,cgreen,cblue,cgry
      data itime/0/
      itime=itime+1

!  Save initial graylevel
      cgrysv=cgry

      lt=len(titlein)

      if(itime.eq.1) then
        cmdstr='/Fwidth { newpath 0 0 moveto text false charpath '             &
             //'flattenpath pathbbox'
        call filler
        cmdstr='/ury exch def /urx exch def /lly exch def /llx exch def'
        call filler
        cmdstr='/swidth urx llx sub def /sheight ury lly sub def} def'
        call filler
      endif

!  Get bounding box before rotation, since it is inaccurate if rotation
!  is not a multiple of 90 degrees (PS ref, p.461)

!  This call sets the font size
      ihgt=size*conver/.6    !Must match keksym, keksymc, etc.
      write(cmdstr,'(i3,'' Setf'')')ihgt
      call filler
      cmdstr='/text ('//titlein(1:lt)//') def'
      call filler

!  This defines llx, lly, swidth, etc.
      cmdstr='gsave Fwidth grestore'     !This eliminates charpath from
      call filler                        !stroking the character outline

!  Define edge width (a fraction of bounding box height)
      write(scr,'(f5.2)')edg
      call blkstp(scr,80,scr,lsc)
      cmdstr='/edg '//scr(1:lsc)//' sheight mul def'
      call filler
      cmdstr='/tothgt sheight 2 edg mul add def' //                            &
      ' /totwid swidth 2 edg mul add def'
      call filler

!  Shift origin. Remember, xp,yp represent the coords of the box,
!  not the text
      call plot(xp,yp,-3)

!  Set angle.  Then everything is in relation to text coordinate frame
      if(ang.ne.0.) call rotate(ang)

!  Figure out where to put the lower left corner of the box ......
      write(cmdstr,'(''/xsh '',i1,'' totwid mul 2 div def'')')jusx
      call filler
      write(cmdstr,'(''/ysh '',i1,'' tothgt mul 2 div def'')')jusy
      call filler

!  ......and move the new origin there
      cmdstr='xsh neg ysh neg translate'
      call filler

!  Draw and fill in the background box (lower left corner is at (0,0))

!  Set background (fill) graylevel
      call setgry(bgry)        !Background fill

      cmdstr= 'Np 0 0 moveto totwid 0 rlineto'
      call filler
      cmdstr='0 tothgt rlineto totwid neg 0 rlineto cf'
      call filler

!  Set foreground (text) graylevel
      call setgry(fgry)        !Foreground letters

!  Draw characters
!  Text now starts at (edg-llx,edg-lly)
      cmdstr='edg llx sub edg lly sub moveto text show'
      call filler

!  Reset graylevel to original
      call setgry(cgrysv)

!  Unshift origin
      cmdstr='xsh ysh translate'
      call filler

!  Reset angle
      if(ang.ne.0.) call rotate(-ang)

!  Reset origin
      call plot(-xp,-yp,-3)
      return
      end
!*****INTEGRAL
      subroutine integral(xp,yp,size,ang,lower,nl,lupper,nu)
!  this routines plots an integral at position xp,yp.
!  size is the height of the (subsequent) integrand.
!  the integral has been empirically enlarged.
!  lower,nl are the hollerith lower limit and number of chars
!  lupper,nu are the hollerith upper limit and number of chars
!  stringlength of integral has been increased to slightly shift the
!  integrand (call by keksym with 999., etc.)

      character*132 cmdstr
      character*132 curfnt,scrc
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titlec
      character*1 bslash
      character * (*) lower,lupper
!      dimension lower(nl),lupper(nu)
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

      mjus=0    !Force left justified

!  Stroke previous paths before this write
      cmdstr='S'
      call filler

      bslash=char(92)

!      pi=4*abs(atan(1.)) ! no need to call transcendental function to get pi
      nchar=1
      ititle=362
   10 continue

!  Choose proper font height, using current font
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      iht=iht*1.5         !Enlarge integral

      offset=iht*.13      !Define amount of offset to shift integral "down"

      arg=ang*rdpdeg
      xoff=offset*sin(arg)
      yoff=-offset*cos(arg)

!  Insert limits
      call plot(xp+xoff/conver,yp+yoff/conver,-3)

!  Set angle
      if(ang.ne.0.) then
        cmdstr=' '
        write(cmdstr,'(F7.1,'' rotate'')') ang
        call filler
      endif

      limiht=iht*.2      !Limits height
      limiht=max0(limiht,1)
      ss=limiht/conver
      nfontsv=nfont
      if(nu.ne.0)then
        if(nu.eq.-999) call setfnt(29)
        call keksym(.3*iht/conver,.82*iht/conver,ss,lupper,0.,nu,0)
        if(nu.eq.-999) call setfnt(nfontsv)
      endif

      if(nl.ne.0)then
        if(nl.eq.-999) call setfnt(29)
        call keksym(.2*iht/conver,-.13*iht/conver,ss,lower,0.,nl,0)
        if(nl.eq.-999) call setfnt(nfontsv)
      endif
      call setfnt(nfontsv)

!  Reset angle
      if(ang.ne.0.) then
        cmdstr=' '
        write(cmdstr,'(F7.1,'' rotate'')') -ang
        call filler
      endif

      call plot(-(xp+xoff/conver),-(yp+yoff/conver),-3)

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')')       &
      IHT
      call filler

      write(titlec,'(A1,I10)')bslash,ititle
      call blkstp(titlec,80,titlec,numc)

   15 continue
      mchar=numc
      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif

      rsize=size
!  character space height is 2.0 x char height
!  character space width is 1.5 x char width
!  actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*rdpdeg
      xoff=offset*sin(arg)
      if(xpos.eq.999.) then
        write(cmdstr,'(''/Xpos Xpres '',f8.2,''add def'')')xoff
      else
        cmdstr=' '
        xnew=xp*conver+xoff
        write(cmdstr,'(f8.2,'' Xposd'')')xnew
      endif
      call filler

      yoff=-offset*cos(arg)
      if(ypos.eq.999.) then
        write(cmdstr,'(''/Ypos Ypres '',f8.2,''add def'')')yoff
      else
        cmdstr=' '
        ynew=yp*conver+yoff
        write(cmdstr,'(f8.2,'' Yposd'')')ynew
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'integral, zero used')
        njus=0
      endif

      shift=njus*strlen/2.
!     xpos=xpos-cos(arg)*shift
!     ypos=ypos-sin(arg)*shift
!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lendi'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      if(xarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Xposjd'')')xarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      if(yarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Yposjd'')')yarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      call filler

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Save current integral coordinates & size, even if not used again
      write(cmdstr,                                                            &
      '(''/Xpint Xp def /Ypint Yp def /Intsz '',i3,'' def'')')iht
      call filler

!  Reset font
      call setfnt(nfont)

!  Start next char at .5 char width away
      argdeg=arg/rdpdeg

!  Reset Xpos &Ypos from integral offset
      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos '',f8.2,'' sub def'')')xoff

      call filler
      write(cmdstr,'(''  /Ypos Ypos '',f8.2,'' sub def'')')yoff
      call filler

      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

      return
      end
!*****KEKEXP
      subroutine kekexp(xp,yp,size,fpn,ang,ndec,mjus)
      character*80 ichrnum

      fnum=fpn
!  Get number in character form
      call numsym(fnum,ndec,ichrnum,ndigit,.true.)
      call keksymc(xp,yp,size,ichrnum,ang,ndigit,mjus)
      return
      end
!*****KEKFLT
      subroutine kekflt(xp,yp,size,fpn,ang,ndec,mjus)
      character*80 ichrnum

      fnum=fpn
!  Get number in character form
      call numsym(fnum,ndec,ichrnum,ndigit,.false.)
      call keksymc(xp,yp,size,ichrnum,ang,ndigit,mjus)
      return
      end
!*****KEKNUM
      subroutine keknum(xp,yp,size,fpn,ang,ndec,mjus)
!  Just assume that user really wants kekflt.
      call kekflt(xp,yp,size,fpn,ang,ndec,mjus)
      return
      end
!*****KEKSYM
      subroutine keksym(xp,yp,size,ltitle1,ang,nchar1,mjus)
      parameter (pi=3.14159265358979, rdpdeg=pi/180.)
      character * (*) ltitle1
      character * 80 ltitle
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titlec,titleb
      character*1 bslash
      equivalence(ititle,ltitle)
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

!  Stroke previous paths before this write
      cmdstr='S'
      call filler

      bslash=char(92)

!      pi=4*abs(atan(1.)) !


      if(nchar1.eq.-999) then !octal code
        ltitle = ltitle1(1:4)  ! copy first 4 bytes to equivalence area
        nchar = 1
      else
        nchar = nchar1
        ltitle = ltitle1(1:iabs(nchar))
        titlec = ltitle
      endif

!  Choose proper font height, using current font
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 FACTOR IS EMPIRICAL

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')iht
        call filler
        ifntsz=iht
      endif

      if(nchar1.eq.-999) then        !Octal code
        write(titlec,'(A1,I10)')bslash,ititle
        call blkstp(titlec,80,titlec,numc)
      else if (nchar1.gt.0) then
!        Check if titlec contains ( or ) or \.  These characters must be treated
!        specially by preceding them with a "\".  Do this to ( and ) even though
!        they might be balanced, i.e. () within a string, which can be treated
!        normally.

          titleb=titlec
          numc=0
          do m=1,mchar
            if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m)        &
              .eq.bslash) then
               numc=numc+1
               titlec(numc:numc)=bslash
             endif
             numc=numc+1
             titlec(numc:numc)=titleb(m:m)
           enddo
        else
!         Unedited character string.  User must be sure of correctness.
          numc = - nchar1
        endif

      mchar=numc

      xpos=xp
      ypos=yp
      if(nchar.eq.-1) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat

      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(ypos.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'KEKSYM, zero used')
        njus=0
      endif
!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'
      arg=ang*4.*abs(atan(1.))/180.
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.

      if(xarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Xposjd'')')xarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if(yarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Yposjd'')')yarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        if(xarg.ne.0.) then
          scrc=' '
          write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr             &
          (scrc,132))
        endif

        if(yarg.ne.0.) then
          scrc=' '
          write(scrc,                                                          &
          '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr             &
          (scrc,132))
        endif
      endif

      call filler

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg = arg / rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

      return
      end
!*****KEKSYMC
!  this routine is the same as keksym, but it accepts character strings instead
!  of integer (hollerith) strings
!  this routine does not support octal codes for special characters
!  (use keksym instead)
      subroutine keksymc(xp,yp,size,titlein,ang,nchar1,mjus)
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*(*) titlein
      character*80 titleb,titlec
      character*1 bslash
      dimension ltitle(20)
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

!  Stroke previous paths before this write
      cmdstr='S'
      call filler

      bslash=char(92)

!      pi=4.*abs(atan(1.)) ! no need for transcendental function to find pi

      nchar=nchar1
      titlec=titlein
      if(iabs(nchar).lt.80)titlec(nchar+1:80)=' '
! rrd change for g95 compatibility (line below not needed)
!     read(titlec,'(20a4)')ltitle

   10 continue

!  Choose proper font height, using current font
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')iht
        call filler
        ifntsz=iht
      endif

!  Check if titlec contains ( or ) or \.  These characters must be treated
!  specially by preceding them with a "\".  Do this to ( and ) even though
!  they might be balanced, i.e. () within a string, which can be treated
!  normally.

      titleb=titlec
      numc=0
      do m=1,mchar
        if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m).eq.        &
        bslash) then
          numc=numc+1
          titlec(numc:numc)=bslash
        endif
        numc=numc+1
        titlec(numc:numc)=titleb(m:m)
      enddo

      mchar=numc
      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.

      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(ypos.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'Incorrect justification code ',i5,'found in ',                &
      'KEKSYMC, zero used')
        njus=0
      endif

!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      if(xarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Xposjd'')')xarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      if(yarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Yposjd'')')yarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        if(xarg.ne.0.) then
          scrc=' '
          write(scrc,                                                          &
          '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr             &
          (scrc,132))
        endif
        if(yarg.ne.0.) then
          scrc=' '
          write(scrc,                                                          &
          '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr             &
          (scrc,132))
        endif
      endif
      call filler

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg=arg/rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

      return
      end
!*****KEKSYMO
      subroutine keksymo(xp,yp,size,ltitle1,ang,nchar1,mjus)
!  this routine is same as keksym, but outlines text rather than
!  filling
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titlec,titleb
      character*1 bslash
      dimension ltitle(20),ltitle1(20)
      equivalence(ititle,ltitle(1))
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

      bslash=char(92)

!      pi=4.*abs(atan(1.)) ! no need for transcendental function to find pi

      if(nchar1.eq.-999) then !octal code
        do n=1,20
          ltitle(n)=ltitle1(n)
        enddo
        nchar=1
      else
        nchar=nchar1
        write(titlec,'(20a4)')ltitle1
        if(iabs(nchar).lt.80)titlec(nchar+1:80)=' '
        read(titlec,'(20a4)')ltitle
      endif

!  Choose proper font height, using current font
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 FACTOR IS EMPIRICAL

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')IHT
        call filler
        ifntsz=iht
      endif

      if(nchar1.eq.-999) then        !octal code
        write(titlec,'(a1,i10)')bslash,ititle
        call blkstp(titlec,80,titlec,numc)
      else
        write(titlec,'(20a4)')ltitle
!       Check if titlec contains ( or ) or \.  These characters must be treated
!       specially by preceding them with a "\".  Do this to ( and ) even though
!       they might be balanced, i.e. () within a string, which can be treated
!       normally.

        titleb=titlec
        numc=0
        do m=1,mchar
          if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m)          &
          .eq.bslash) then
            numc=numc+1
            titlec(numc:numc)=bslash
          endif
          numc=numc+1
          titlec(numc:numc)=titleb(m:m)
        enddo
      endif

      mchar=numc
      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.

      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(ypos.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'Incorrect justification code ',i5,'found in ',                &
      'keksymo, zero used')
        njus=0
      endif

!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      if(xarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Xposjd'')')xarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if(yarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Yposjd'')')yarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        if(xarg.ne.0.) then
          scrc=' '
          write(scrc,                                                          &
          '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr             &
          (scrc,132))
        endif

        if(yarg.ne.0.) then
          scrc=' '
          write(scrc,                                                          &
          '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr             &
          (scrc,132))
        endif
      endif

      call filler

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') false charpath stroke'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg=arg/rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

      return
      end
!*****LENSTR
      function lenstr(string,ls)

!  This routine finds actual length of string by eliminating trailing blanks
      character*(*) string

      do i=ls,1,-1
        is=i
        if(string(i:i).ne.char(32)) goto 10
      enddo
      is=0
   10 lenstr=is
      return
      end
!*****LINET
      subroutine linet(xcall,ycall)
!  This routine is same as call to plot with ip=2, but does not
!  perform the "stroke".
      character*132 cmdstr
      common/plt1/cmdstr
      common/cnvcom/conver

      write(cmdstr,'(2f8.2,'' L'')')xcall*conver,ycall*conver
      call filler

      return
      end
!*****MINMAX
      subroutine minmax (z,l,mm,nn,ssizem,aash,joffdt,scal)

! This routine finds relative minimums and maximums.  a relative minimum
! (or maximum) is defined to be the lowest (or highest) point within
! a certain neighborhood of the point.  the neighborhood used here
! is + or - mn in the x direction and + or - nm in the y direction.

! Originator       david kennison

      dimension       z(l,nn)

      common/conre1/ioffp,spval
      common/scales/xor,xsc,yor,ysc

      fx(x,y) = x
      fy(x,y) = y

      m = mm
      n = nn
      sizem = ssizem
      ash = aash
      ioffdt = joffdt
      mn=min0(15,max0(2,ifix(float(m)/10.)))
      nm=min0(15,max0(2,ifix(float(n)/10.)))
      nm1 = n-1
      mm1 = m-1
!
! line loop follows - the complete two-dimensional test for a minimum or
! maximum of the field is only performed for points which are minima or
! maxima along some line - finding these candidates is made efficient by
! using a count of consecutive increases or decreases of the function
! along the line
!
      do 110 jp=2,nm1
        im = mn-1
        ip = -1
        go to 105
!
! control returns to statement 10 as long as the function is increasing
! along the line - we seek a possible maximum
!
   10   ip = ip+1
        aa = an
        if (ip .eq. mm1) go to 25
        an = z(ip+1,jp)
        if (ioffp.ne. 0.and. an.eq.spval) go to 100
!# modified: 10/9/2008
        if (aa.LT.an) then
           im = im+1
        elseif (aa.EQ.an) then
           im = 0
        end if
        go to 10
!# modified: 10/9/2008

! Function decreased - test for maximum on line

   25   if (im .ge. mn) go to 30
        is = max0(1,ip-mn)
        it = ip-im-1
        if (is .gt. it) go to 30
        do ii=is,it
        if (aa .le. z(ii,jp)) go to 50
      enddo
   30 is = ip+2
      it = min0(m,ip+mn)
      if (is .gt. it) go to 40
      do ii=is,it
      if (ioffp.eq. 0.or. z(ii,jp).ne.spval) go to 35
      ip = ii-1
      go to 100
   35 if (aa .le. z(ii,jp)) go to 50
      enddo

! We have maximum on line - do two-dimensional test for maximum of field

   40 js = max0(1,jp-nm)
      jt = min0(n,jp+nm)
      is = max0(1,ip-mn)
      it = min0(m,ip+mn)
      do 45 jk=js,jt
      if (jk .eq. jp) go to 45
      do ik=is,it
      if (z(ik,jk).ge.aa .or. (ioffp.ne. 0.and. z(ik,jk).eq.spval))            &
      go to 50
      enddo
   45 continue

      x = float(ip)
      y = float(jp)
      xx1=(fx(x,y)-xor)/xsc
      yy1=(fy(x,y)-yor)/ysc
      call keksymc(xx1-.06,yy1-.06,.13,'H',0.,1,0)
      as=aa*scal
      call kekflt(xx1,yy1-.25,.1,as,0.,1,1)

   50 im = 1
!# modified: 10/9/2008
      if (ip.GE.mm1) go to 110
!# modified: 10/9/2008

! Control returns to statement 20 as long as the function is decreasing
! along the line - we seek a possible minimum

   55 ip = ip+1
      aa = an
      if (ip .eq. mm1) go to 70
      an = z(ip+1,jp)
      if (ioffp.ne. 0.and. an.eq.spval) go to 100
!# modified: 10/9/2008
      if (aa.GT.an) then
         im = im+1
      elseif (aa.EQ.an) then
         im = 0
      end if
      go to 55
!# modified: 10/9/2008

! Function increased - test for minimum on line

   70 if (im .ge. mn) go to 75
      is = max0(1,ip-mn)
      it = ip-im-1
      if (is .gt. it) go to 75
      do ii=is,it
      if (aa .ge. z(ii,jp)) go to 95
      enddo
   75 is = ip+2
      it = min0(m,ip+mn)
      if (is .gt. it) go to 85
      do ii=is,it
      if (ioffp.eq. 0.or. z(ii,jp).ne.spval) go to 80
      ip = ii-1
      go to 100
   80 if (aa .ge. z(ii,jp)) go to 95
      enddo

! We have minimum on line - do two-dimensional test for minimum of field

   85 js = max0(1,jp-nm)
      jt = min0(n,jp+nm)
      is = max0(1,ip-mn)
      it = min0(m,ip+mn)
      do 90 jk=js,jt
      if (jk .eq. jp) go to 90
      do ik=is,it
      if (z(ik,jk).le.aa .or. (ioffp.ne. 0.and. z(ik,jk).eq.spval))            &
      go to 95
      enddo
   90 continue

      x = float(ip)
      y = float(jp)
      xx1=(fx(x,y)-xor)/xsc
      yy1=(fy(x,y)-yor)/ysc
      call keksymc(xx1-.06,yy1-.06,.13,'L',0.,1,0)
      as=aa*scal
      call kekflt(xx1,yy1-.25,.1,as,0.,1,1)
   95 im = 1
!# modified: 10/9/2008
      if (ip.LT.mm1) go to 10
      if (ip.GE.mm1) go to 110
!# modified: 10/9/2008

! Skip special values on line

  100 im = 0
  105 ip = ip+1
      if (ip .ge. mm1) go to 110
      if (ioffp.ne. 0.and. z(ip+1,jp).eq.spval) go to 100
      im = im+1
      if (im .le. mn) go to 105
      im = 1
      an = z(ip+1,jp)
!# modified: 10/9/2008
      if (z(ip,jp).LT.an) then
         go to 10
      elseif (z(ip,jp).EQ.an) then
         im = 0
         go to 10
      else   
         go to 55
      end if 
!# modified: 10/9/2008
  110 continue

      return

! ****************************** entry pntval **************************
      entry pntval (z,l,mm,nn,ssizem,aash,joffdt,scal)

      m = mm
      n = nn
      sizem = ssizem
      ash = aash
      ioffdt = joffdt
      ii = (m-1+24)/24
      jj = (n-1+48)/48
      niq = 1
      njq = 1
      do j=njq,n,jj
      y = j
      do i=niq,m,ii
      x = i
      zz = z(i,j)
      enddo
      enddo
      return
      end
!*****MOVET
      subroutine movet(xcall,ycall)
!  This routine is same as call to plot with ip=3, but does not
!  perform the "stroke".
      character*132 cmdstr
      common/plt1/cmdstr
      common/cnvcom/conver

      write(cmdstr,'(2f8.2,'' M'')')xcall*conver, ycall*conver
      call filler

      return
      end
!*****NEWDEV
      subroutine newdev (iusnb,ich)
!
!  This routine changes the output device/filename to the user
!  specifications. it must be called before psinit.
!  iusnb: new file name ascii string

!  Make ich non-essential.   Just use whole string

      character*80 fileout
      character*(*) iusnb
      common/io/fileout,inew
      lus=lenstr(iusnb,len(iusnb))
      fileout=' '
      fileout=iusnb(1:lus)
!  Set inew to 999 so that init knows to use this name
      inew=999
      iich=ich      !to avoid unused condition
      return
      end
!*****NUMBER
      subroutine number(xp,yp,size,fpn,ang,ndec)
      mj=0
      call kekflt(xp,yp,size,fpn,ang,ndec,mj)
      return
      end
!*****NUMSYM
      subroutine numsym(fpn,ndec,itext,nchar,eform)

!  eform:  true for exponential format
!          false for floating pt format

      logical eform
      character*4 itext(20)
      character*10 ifrmt
      character*80 a
      a=' '

!  Check if ndec is valid
      if(ndec.gt.15) then
        print 100, ndec
  100 format(1x,'In call to numsym, ndec gt 15 ',i20,' program',               &
      ' abandoned')
        stop
      else if(eform.and.ndec.lt.0) then
        print 110, ndec
  110 format(1x,'in numsym, exponential format specified with ','ndec= '       &
      ,i3,' program abandoned')
        stop
      endif

      if(eform) then
        write(ifrmt,'(''(1pe16.'',i2,'')'')')ndec
      else if(ndec.lt.0) then
        ifrmt='(f16.1)'
      else
        write(ifrmt,'(''(f16.'',i2,'')'')')ndec
      endif

      write(a,ifrmt)fpn

!  Strip off all blanks in a
      call blkstp(a,80,a,nchar)
      ipos=index(a,'.')
      if(eform) then
        if(ndec.eq.0) then
!         Delete characters between '.' and 'E'
          ie=index(a,'E')
          do n=ie,nchar
          nind=ipos+n-ie+1
          a(nind:nind)=a(n:n)
        enddo
        nchar=nchar-(ie-ipos-1)
      endif
      else if(ndec.lt.0) then
      nchar=ipos-1
      else if(ndec.eq.0) then
      nchar=ipos
      endif
      read(a,'(20a4)')itext

      return
      end
!*****ONEHLF
      subroutine onehlf(xp,yp,size,ang,mjus)
!  this routine draws the 1/2 character
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat
      common/fntcom/curfnt,ifntsz,nfont

      bslash=char(92)

!RRD unitialized variable
      xpos=xp
      ypos=yp

!  Save current font
      nfsav=nfont

!  Make new font with iso latin 1 encoding
!  This code taken from PostScript Lang. Ref. Man. Section 5.6.1

      cmdstr='Curfnt dup length dict begin'
      call filler
      cmdstr='  {1 index /FID ne {def} {pop pop } ifelse} forall'
      call filler
      cmdstr='  /Encoding ISOLatin1Encoding def'
      call filler
      cmdstr='  currentdict'
      call filler
      cmdstr='end'
      call filler
      cmdstr='/Latino exch definefont pop'
      call filler

!      pi=4.*abs(atan(1.)) ! don't use transcendental function to find pi
      nchar=1
!  Choose proper character height
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 FACTOR IS EMPIRICAL

      cmdstr=' '
      write(cmdstr,'(''/Latino findfont '',I3,'' scalefont setfont'')')        &
      iht
      call filler

      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width

      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'ONEHLF, zero used')
        njus=0
      endif

      cmdstr='('//bslash//'275) Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif


!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      scrc='('//bslash//'275) show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg=arg/rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')') argdeg
      call filler

!  Reset font and size
      call setfnt(nfsav)

      cmdstr=' '
      write(cmdstr,'(I3,'' Setf'')')IFNTSZ
      call filler

      return
      end
!*****OVERBAR
      subroutine overbar(xp,yp,size,ltitle,ang,nchar,mjus)
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titleb,titlec
      character*1 bslash
      dimension ltitle(20)
      common/lcom/curwid
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

!      pi=4.*abs(atan(1.)) ! get pi more readily by parm
      bslash=char(92)

!  Choose proper font height, using current font
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')iht
        call filler
        ifntsz=iht
      endif

      write(titlec,'(20a4)')ltitle

      titleb=titlec
      numc=0
      do m=1,mchar
        if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m)            &
          .eq.bslash) then
            numc=numc+1
            titlec(numc:numc)=bslash
          endif
          numc=numc+1
          titlec(numc:numc)=titleb(m:m)
      enddo
      mchar=numc

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*rdpdeg
      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(ypos.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'OVERBAR, zero used')
        njus=0
      endif

!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))

      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

!  Determine how high to draw overbar. if any capitals exist in text,
!  raise overbar up slightly.
      do n=1,nchar
        if(ichar(titlec(n:n)).ge.65.and.ichar(titlec(n:n)).le.90) then
          dist=size+.08
        elseif(ichar(titlec(n:n)).eq.98.or.ichar(titlec(n:n)).eq.100.or.  &
               ichar(titlec(n:n)).eq.102.or.ichar(titlec(n:n)).eq.104.or. &
               ichar(titlec(n:n)).eq.105.or.ichar(titlec(n:n)).eq.106.or.  &
               ichar(titlec(n:n)).eq.107.or.ichar(titlec(n:n)).eq.108.or.  &
               ichar(titlec(n:n)).eq.116) then

          dist=size+.08
        else
          dist=size+.01
        endif
      enddo

!  Determine start and end points for overbar
      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

      cmdstr=' '
      write(cmdstr, '(''/Xpos Strlen '',F6.1,'' cos mul Xpos add def'')'       &
      ) ang

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Strlen '',F6.1,'' sin mul Ypos add def'')'       &
      ) ang

      call filler

!  Get current linewidth
      curwids=curwid

      call setlw(.01)       !Thickness of overbar
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')ang
      call filler

      return
      end
!*****OVERSBSP
      subroutine oversbsp(xp,yp,size,ltitle,ang,nchar,mjus, isub,nsub,         &
      isup,nsup)
!  This routines draws an overbar over a subscripted variable, i.e.
!  the subscript is also overbarred.
!  isub is the subscript character(s); nsub is the number of subscript chars.
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titleb,titlec
      character*1 bslash
!rrd  dimension ltitle(20),isub(20),isup(20)
      character*80 ltitle,isub,isup
      common/lcom/curwid
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

!      pi=4.*abs(atan(1.)) ! don't use atan to get pi
      bslash=char(92)

!  Choose proper font height, using current font
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(i3,'' Setf'')')iht
        call filler
        ifntsz=iht
      endif

!rrd  write(titlec,'(20a4)')ltitle
      titlec=ltitle

      titleb=titlec
      numc=0
      do m=1,mchar
        if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m)            &
          .eq.bslash) then
            numc=numc+1
            titlec(numc:numc)=bslash
          endif
          numc=numc+1
          titlec(numc:numc)=titleb(m:m)
      enddo
      mchar=numc

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size

!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*rdpdeg
      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        NJUS=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Xposd'')')xp*conver
      endif
      call filler

      IF(YPOS.EQ.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'OVERSBSP, zero used')
        njus=0
      endif

!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc,                                                            &
        '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
        scrc=' '
        write(scrc,                                                            &
        '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

      dist=size+.12

!  Determine start and end points for overbar

      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

!  Define length of subscript (remember font size is 3/4 of size)
!rrd  write(titlec,'(20A4)')isub
      titlec=isub
      cmdstr='('//titlec(1:nsub)//') Lenssd'
      call filler

!  Rename subscript length to another variable
      cmdstr=' '
      cmdstr='/Strlensb Strlenss def'
      call filler

!  Define length of superscript (remember font size is 3/4 of size)
!rrd  write(titlec,'(20A4)')isup
      titlec=isup
      cmdstr='('//titlec(1:nsup)//') Lenssd'
      call filler

!  Redefine Strlenss to be longer of Strlensb, Strlenss; i.e set
!  Strlenss to Strlensb if Strlensb is longer

      cmdstr=' '
      cmdstr='Strlensb Strlenss ge {/Strlenss Strlensb def} if'
      call filler

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Strlen Strlenss add '',F6.1,'' cos mul Xpos add def'')')       &
      ang
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Strlen Strlenss add '',F6.1,'' sin mul Ypos add def'')')       &
      ang
      call filler

!  Get current linewidth
      curwids=curwid
      call setlw(.02)
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

!  draw subscript and superscript
!  even though we have drawn title already, draw again to get in position
!  for call to subsup

      call keksym(xp,yp,size,ltitle,ang,nchar,mjus)
      call subsup(isub,nsub,isup,nsup,size,ang)

      return
      end
!*****OVERSBSPG
      subroutine oversbspg(xp,yp,size,ititle,ang,nchar,mjus,isub,nsub,         &
      isup,nsup)
!  this routine draws an overbarred subscripted and superscripted Greek symbol
      character*80 grkchar,titlec,charout*1
!rrd  dimension isup(20),isub(20)
      character*80 isup,isub
      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/lcom/curwid
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat
      data grkchar/'ABGDEZHQIKLMNXOPRSTUFCYWabgdezhJiklmnxoprstujcywqf'/

      bslash=char(92)

      if(xp.eq.999..or.yp.eq.999.) then
!       Save beginning character postion
        cmdstr='/Xpbeg Xpres def /Ypbeg Ypos def'
        call filler
      endif

!  Save current font
      nfsav=nfont

!  Choose proper character height
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')')       &
      IHT
      call filler

      charout=grkchar(ititle:ititle)
      if(charout.eq.'?') then
        print *,'Greek character ',ititle,' is not available ',                &
        'in postscript'
        print *,'a blank is being substituted'
        charout=' '
      endif

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'keksym, zero used')
        njus=0
      endif

      if(charout.ne.'U') then
        cmdstr='('//charout//') Lend'
      else
        cmdstr='('//bslash//'241) Lend'
      endif
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif


!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if(charout.ne.'U') then
        scrc='('//charout//') show'
      else
        scrc='('//bslash//'241) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )


!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

!  Determine how high to draw overbar. if any capitals exist in text,
!  raise overbar up slightly.
      dist=size+.07
      icap=0
      do n=1,nchar
        if(ititle.le.24.or.ititle.eq.26.or.ititle.eq.28.or.ititle.eq.30.or. &
           ititle.eq.32.or.ititle.eq.35.or.ititle.eq.38.or.ititle.eq.49.or. &
           ititle.eq.50) then
          dist=size+.1
          icap=1
          goto 10
        endif
      enddo
   10 continue

!  Determine start and end points for overbar
      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

!  Define length of subscript (remember font size is 3/4 of size)
!rrd  write(titlec,'(20A4)')isub
      titlec=isub
      cmdstr='('//titlec(1:nsub)//') Lenssd'
      call filler

!  Rename subscript length to another variable
      cmdstr=' '
      cmdstr='/Strlensb Strlenss def'
      call filler

!  Define length of superscript (remember font size if 3/4 of size)
!rrd  write(titlec,'(20A4)')isup
      titlec=isup
      cmdstr='('//titlec(1:nsup)//') Lenssd'
      call filler

!  Redefine Strlenss to be longer of Strlensb, Strlenss; i.e set
!  Strlenss to Strlensb if Strlensb is longer

      cmdstr=' '
      cmdstr='Strlensb Strlenss ge {/Strlenss Strlensb def} if'
      call filler

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Strlen Strlenss add '',F6.1,'' cos mul Xpos add def'')')       &
      ang
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Strlen Strlenss add '',F6.1,'' sin mul Ypos add def'')')       &
      ang
      call filler

!  Get current linewidth
      curwids=curwid
      call setlw(.02)
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

!  Reset current font
      call setfnt(nfsav)

!  draw subscript and superscript
!  even though we have drawn title already, draw again to get in position
!  for call to subsup.

      call grksym(xp,yp,size,ititle,ang,nchar,mjus)

      if(icap.eq.0) then
        call super(isup,-nsup,size,ang)
      else
        call super(isup,nsup,size,ang)
      endif

!  We must reset position to the original if continuation.
      if(xp.eq.999..or.yp.eq.999.)then
        cmdstr='/Xpres Xpbeg def /Ypres Ypbeg def'
        call filler
      endif

      call grksym(xp,yp,size,ititle,ang,nchar,mjus)
      call subber(isub,nsub,size,ang)

      return
      end
!*****OVERSUB
      subroutine oversub(xp,yp,size,ltitle,ang,nchar,mjuS,isub,nsub)
!  This routines draws an overbar over a subscripted variable, i.e.
!  the subscript is also overbarred.
!  isub is the subscript character(s); nsub is the number of subscript chars.
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titleb,titlec
      character*1 bslash
!rrd  dimension ltitle(20),ISUB(20)
      character*80 ltitle,isub
      common/lcom/curwid
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

!      pi=4.*abs(atan(1.)) don't use atan to get pi
      bslash=char(92)

!  Choose proper font height, using current font
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')iht
        call filler
        ifntsz=iht
      endif

      write(titlec,'(20a4)')ltitle

      titleb=titlec
      numc=0
      do m=1,mchar
        if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m)            &
          .eq.bslash) then
            numc=numc+1
            titlec(numc:numc)=bslash
          endif
          numc=numc+1
          titlec(numc:numc)=titleb(m:m)
      enddo
      mchar=numc

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg = ang * rdpdeg
      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Xposd'')')xp*conver
      endif
      call filler

      if(ypos.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'keksym, zero used')
        njus=0
      endif

!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

!  Determine how high to draw overbar. if any capitals exist in text,
!  raise overbar up slightly.
      dist=size+.01
      do n=1,nchar
        if(ichar(titlec(n:n)).ge.65.and.ichar(titlec(n:n)).le.90) then
          dist=size+.08
          goto 10
        elseif(ichar(titlec(n:n)).eq.98..or.ichar(titlec(n:n)).eq.100.or.     &
            ichar(titlec(n:n)).eq.102.or.ichar(titlec(n:n)).eq.104.or.        &
            ichar(titlec(n:n)).eq.105.or.ichar(titlec(n:n)).eq.106.or.        &
            ichar(titlec(n:n)).eq.107.or.ichar(titlec(n:n)).eq.108.or.        &
            ichar(titlec(n:n)).eq.116) then
          dist=size+.08
          goto 10
        endif
      enddo
   10 continue

!  Determine start and end points for overbar
      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

!  Define length of subscript (remember font size is 3/4 of size)
!rrd  write(titlec,'(20A4)')isub
      titlec=isub
      cmdstr='('//titlec(1:nsub)//') Lenssd'
      call filler

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Strlen Strlenss add '',F6.1,'' cos mul Xpos add def'')') ang
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Strlen Strlenss add '',F6.1,'' sin mul Ypos add def'')') ang
      call filler

!  Get current linewidth
      curwids=curwid
      call setlw(.02)
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

!  Start next char at .5 char width away
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')ang
      call filler

!  Draw subscript
      call subber(isub,nsub,size,ang)

      return
      end
!*****OVERSUBG
      subroutine oversubg(xp,yp,size,ititle,ang,nchar,mjus,isub,nsub)
!  this routine draws an overbarred subscripted greek symbol
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*80 grkchar,titlec,charout*1
!rrd  dimension isub(20)
      character*80 isub
      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/lcom/curwid
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

      data grkchar/'ABGDEZHQIKLMNXOPRSTUFCYWabgdezhJiklmnxoprstujcywqf'/

      bslash=char(92)
!      pi=4.*abs(atan(1.)) ! don't use atan to get pi

!  Save current font
      nfsav=nfont

!  Choose proper character height
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')') IHT
      call filler

      charout=grkchar(ititle:ititle)
      if(charout.eq.'?') then
        print *,'Greek character ',ititle,' is not available ', 'in postscript'
        print *,'A blank is being substituted'
        charout=' '
      endif

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'OVERSUBG, zero used')
        njus=0
      endif

      if(charout.ne.'U') then
        cmdstr='('//charout//') Lend'
      else
        cmdstr='('//bslash//'241) Lend'
      endif
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      if(charout.ne.'U') then
        scrc='('//charout//') show'
      else
        scrc='('//bslash//'241) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

!  Determine how high to draw overbar. if any capitals exist in text,
!  raise overbar up slightly.
      dist=size+.05
      icap=0
      do n=1,nchar
        if(ititle.le.24.or.ititle.eq.26.or.ititle.eq.28.or.ititle.eq.30.or. &
           ititle.eq.32.or.ititle.eq.35.or.ititle.eq.38.or.ititle.eq.49.or. &
           ititle.eq.50) then
          dist=size+.08
        else
          dist=size+.01
        endif
      enddo

!  Determine start and end points for overbar

      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

!  Define length of subscript (remember font size if 3/4 of size)
!  Use defined function LENSSD even though asthetically should be LENSBD
!rrd  write(titlec,'(20A4)')isub
      titlec=isub
      cmdstr='('//titlec(1:nsub)//') Lenssd'
      call filler

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Strlen Strlenss add '',F6.1,'' cos mul Xpos add def'')') ang
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Strlen Strlenss add '',F6.1,'' sin mul Ypos add def'')') ang
      call filler

!  Get current linewidth
      curwids=curwid
      call setlw(.02)
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

!  Start next char at .5 char width away
      argdeg=arg / rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

!  Reset current font
      call setfnt(nfsav)

      call subber(isub,nsub,size,ang)

      return
      end
!*****OVERSUP
      subroutine oversup(xp,yp,size,ltitle,ang,nchar,mjus,isup,nsup)
!  This routines draws an overbar over a superscripted variable, i.e.
!  the superscript is also overbarred.
!  isup is the superscript character(s)
!  nsup is the number of superscript chars.
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titleb,titlec
      character*1 bslash
!rrd  dimension ltitle(20),ISUP(20)
      character*80 ltitle,isup
      common/lcom/curwid
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

!      pi=4.*abs(atan(1.)) ! don't use atan to find pi
      bslash=char(92)

!  Choose proper font height, using current font
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(i3,'' Setf'')')iht
        call filler
        ifntsz=iht
      endif

      write(titlec,'(20a4)')ltitle
      titleb=titlec
      numc=0
      do m=1,mchar
       if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')'.or.titleb(m:m).eq.bslash)then
            numc=numc+1
            titlec(numc:numc)=bslash
          endif
          numc=numc+1
          titlec(numc:numc)=titleb(m:m)
      enddo
      mchar=numc

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg = ang * rdpdeg
      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(ypos.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'OVERSUP, zero used')
        njus=0
      endif

!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

!  Determine how high to draw overbar. if any capitals exist in text,
!  raise overbar up slightly.
      dist=size+.08
      icap=0
      do n=1,nchar
        if(ichar(titlec(n:n)).ge.65.and.ichar(titlec(n:n)).le.90) then
          dist=size+.12    !Caps
          icap=1
          goto 10
        elseif(ichar(titlec(n:n)).eq.98 .or.ichar(titlec(n:n)).eq.100.or.  &
               ichar(titlec(n:n)).eq.102.or.ichar(titlec(n:n)).eq.104.or.  &
               ichar(titlec(n:n)).eq.105.or.ichar(titlec(n:n)).eq.106.or.  &
               ichar(titlec(n:n)).eq.107.or.ichar(titlec(n:n)).eq.108.or.  &
               ichar(titlec(n:n)).eq.116) then

          dist=size+.12
          icap=1
          goto 10
        endif
      enddo
   10 continue

!  Determine start and end points for overbar
      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

!  Define length of superscript (remember font size if 3/4 of size)
!rrd  write(titlec,'(20A4)')isup
      titlec=isup
      cmdstr='('//titlec(1:nsup)//') Lenssd'
      call filler

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Strlen Strlenss add '',F6.1,'' cos mul Xpos add def'')') ang
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Strlen Strlenss add '',F6.1,'' sin mul Ypos add def'')') ang
      call filler

!  Get current linewidth
      curwids=curwid
      call setlw(.02)
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

!  Start next char at .5 char width away
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')ang
      call filler

!  Draw superscript
!  Make nsup negative if we are superscripting a lower case character.
!  Super then knows we're calling from oversup and lowers superscript slightly
      if(icap.eq.0) then
        call super(isup,-nsup,size,ang)
      else
        call super(isup,nsup,size,ang)
      endif
      return
      end
!*****OVERSUPG
      subroutine oversupg(xp,yp,size,ititle,ang,nchar,mjus,isup,nsup)
!  this routine draws an overbarred superscripted greek symbol
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*80 grkchar,titlec,charout*1
!rrd  dimension isup(20)
      character*80 isup
      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/lcom/curwid
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

      data grkchar/'ABGDEZHQIKLMNXOPRSTUFCYWabgdezhJiklmnxoprstujcywqf'/

      bslash=char(92)
!      pi=4.*abs(atan(1.)) !

!  Save current font
      nfsav=nfont

!  Choose proper character height
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')') iht
      call filler

      charout=grkchar(ititle:ititle)
      if(charout.eq.'?') then
        print *,'Greek character ',ititle,' is not available ', 'in postscript'
        print *,'a blank is being substituted'
        charout=' '
      endif

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2, '' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'Incorrect justification code ',i5,'found in ',                &
      'OVERSUPG, zero used')
        njus=0
      endif

      if(charout.ne.'U') then
        cmdstr='('//charout//') Lend'
      else
        cmdstr='('//bslash//'241) Lend'
      endif
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      endif


!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      if(charout.ne.'U') then
        scrc='('//charout//') show'
      else
        scrc='('//bslash//'241) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))


!  set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

!  Determine how high to draw overbar. if any capitals exist in text,
!  raise overbar up slightly.
      dist=size+.05
      icap=0
      do  n=1,nchar
        if(ititle.le.24.or.ititle.eq.26.or.ititle.eq.28.or.ititle.eq.30.or.  &
           ititle.eq.32.or.ititle.eq.35.or.ititle.eq.38.or.ititle.eq.49.or.  &
           ititle.eq.50) then
          dist=size+.08
          icap=1
          goto 10
        endif
      enddo
   10 continue

!  Determine start and end points for overbar
      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

!  Define length of superscript (remember font size if 3/4 of size)
!rrd  write(titlec,'(20A4)')isup
      titlec=isup 
      cmdstr='('//titlec(1:nsup)//') Lenssd'
      call filler

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Strlen Strlenss add '',F6.1,'' cos mul Xpos add def'')') ang
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Strlen Strlenss add '',F6.1,'' sin mul Ypos add def'')') ang
      call filler

!  Get current linewidth
      curwids=curwid
      call setlw(.02)
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

!  Start next char at .5 char width away
      argdeg = arg / rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

!  Reset current font
      call setfnt(nfsav)

      if(icap.eq.0) then
        call super(isup,-nsup,size,ang)
      else
        call super(isup,nsup,size,ang)
      endif

      return
      end
!*****OVRGRK
      subroutine ovrgrk(xp,yp,size,ititle,ang,nchar,mjus)
!  this routine uses postscript symbol to generate greek symbols.
!  An overbar is drawn over the greek character.
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*80 grkchar,charout*1

      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/lcom/curwid
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

      data grkchar/'ABGDEZHQIKLMNXOPRSTUFCYWabgdezhJiklmnxoprstujcywqf'/

      bslash=char(92)
!      pi=4.*abs(atan(1.)) !

!  Save current font
      nfsav=nfont

!  Choose proper character height
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')') IHT
      call filler

      charout=grkchar(ititle:ititle)
      if(charout.eq.'?') then
        print *,'Greek character ',ititle,' is not available ', 'in postscript'
        print *,'a blank is being substituted'
        charout=' '
      endif

      xpos=xp
      ypos=yp
      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  character space height is 2.0 x char height
!  character space width is 1.5 x char width
!  actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'incorrect justification code ',i5,'found in ',                &
      'OVRGRK, zero used')
        njus=0
      endif

      if(charout.ne.'U') then
        cmdstr='('//charout//') Lend'
      else
        cmdstr='('//bslash//'241) Lend'
      endif
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

!
!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      if(charout.ne.'U') then
        scrc='('//charout//') show'
      else
        scrc='('//bslash//'241) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Save current character postion
      cmdstr='/Xpsav Xpos def /Ypsav Ypos def'
      call filler

!  Determine how high to draw overbar. if any capitals exist in text,
!  raise overbar up slightly.
      do n=1,nchar
        if(ititle.le.24.or.ititle.eq.26.or.ititle.eq.28.or.ititle.eq.30.or.  &
           ititle.eq.32.or.ititle.eq.35.or.ititle.eq.38.or.ititle.eq.49.or.  &
           ititle.eq.50) then
          dist=size+.08
        else
          dist=size+.01
        endif
      enddo

!  Determine start and end points for overbar
      xarg=dist*sin(arg)
      yarg=dist*cos(arg)

      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos'',f8.2,'' sub def'')')xarg*conver

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

      call filler

      cmdstr='xydef mover'
      call filler

      cmdstr=' '
      write(cmdstr, '(''/Xpos Strlen '',F6.1,'' cos mul Xpos add def'')') ang

      call filler
      cmdstr=' '
      write(cmdstr, '(''/Ypos Strlen '',F6.1,'' sin mul Ypos add def'')') ang

      call filler

!  Get current linewidth
      curwids=curwid
      call setlw(.02)
!  Draw overbar
      cmdstr='xydef lsm /Xpos Xpsav def /Ypos Ypsav def'
      call filler
!  Reset linewidth
      call setlw(curwids)

!  Start next char at .5 char width away
      argdeg = arg / rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

!  Reset current font
      call setfnt(nfsav)

      return
      end
!*****PLOT
      subroutine plot(xcall,ycall,ip)
      character*132 cmdstr
      character*80 scr

      common/pagcom/npage
      common/plt1/cmdstr
      common/cnvcom/conver
      common/outcom/iunit

      ipp=iabs(ip)

      if(ip.eq.999) then   !    Terminate plot session.
!rrd-start
        cmdstr='restore'
        call filler
!rrd-stop
        cmdstr='stroke showpage'
        call filler
        cmdstr='%%Trailer'
        call filler
        write(scr,'(i4)')npage
        call blkstp(scr,80,scr,nch)
        cmdstr='%%Pages: '//scr(1:nch)
        call filler
        cmdstr='%%EOF'
        call filler
        close(iunit)
        return
      endif

!rrd-start
      if(xcall.gt.1.0E+25)xcall=0.0
!rrd-stop

!  Moving pen
      if(ipp.eq.3) then    !Stroke to paint previous path, then moveto
        write(cmdstr,'(2f8.2,'' SM'')')xcall*conver,ycall*conver
      else                 !Lineto
        write(cmdstr,'(2f8.2,'' LSM'')')xcall*conver,ycall*conver
      endif
      call filler

!  Reset origin if ip.lt.0
      if(ip.lt.0) then
        write(cmdstr,'(2f8.2,'' translate'')')xcall*conver,ycall*conver
        call filler
        ipen=ipp
      endif

      return
      end
!*****PLOTND
      subroutine plotnd
      call plot(0.,0.,999)
      return
      end
!*****PLSMIN
      subroutine plsmin(xp,yp,size,ang,mjus)
!  this routine draws the plus/minus character
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat
      common/fntcom/curfnt,ifntsz,nfont

!RRD unitialized variable
      xpos=xp
      ypos=yp

      bslash=char(92)

!  Save current font
      nfsav=nfont

!      pi=4.*abs(atan(1.)) !
      nchar=1

!  Choose proper character height
      mchar=iabs(nchar)

!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')') IHT
      call filler

      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'Incorrect justification code ',i5,'found in ',                &
      'PLSMIN, zero used')
        njus=0
      endif

      cmdstr='('//bslash//'261) Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif


!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      scrc='('//bslash//'261) show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg = arg / rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')') argdeg
      call filler

!  Reset current font
      call setfnt(nfsav)

      return
      end
!*****PRIME
      subroutine prime(xp,yp,size,ang,mjus)
!  this routine draws the prime character (not an apostrophe)
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr,curfnt,scrc
      character*1 bslash
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat
      common/fntcom/curfnt,ifntsz,nfont

!RRD unitialized variable
      xpos=xp
      ypos=yp

      bslash=char(92)
!      pi=4.*abs(atan(1.)) !
      nchar=1
!  Choose proper character height
      mchar=iabs(nchar)

!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')') IHT
      call filler

      if(nchar.lt.0) then
        njus=0
      else
        njus=mjus
      endif
      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat
      arg=ang*4.*abs(atan(1.))/180.
      if(xp.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Xposd'')')xp*conver
      endif
      call filler

      if(yp.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(f8.2,'' Yposd'')')yp*conver
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'Incorrect justification code ',i5,'found in ',                &
      'PRIME, zero used')
        njus=0
      endif

      cmdstr='('//bslash//'242) Lend'
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      scrc=' '
      write(scrc,'(f7.3,'' Xposjd'')')xarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      scrc=' '
      write(scrc,'(f7.3,'' Yposjd'')')yarg
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        scrc=' '
        write(scrc, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        scrc=' '
        write(scrc, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

!
!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      scrc='('//bslash//'242) show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg = arg / rdpdeg
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

!  Reset font
      call setfnt(nfont)

      return
      end
!*****PSINICH
      subroutine psinich(portrait)
!  same as psinit, but is called only by chopit with a single argument
      logical first,portrait,prtrt
      character*132 cmdstr,curfnt
      character tim*8,dat*9
      character*1 timer(8),dater(9)
      equivalence(timer(1),tim),(dater(1),dat)
      common/conre1/ioffp,spval
      common/plt1/cmdstr
      common/cnvcom/conver
      common/plt2/fac
      common/kkplot/szrat
      common/chpcom/ientry,prtrt
      common/fntcom/curfnt,ifntsz,nfont
      common/pagcom/npage
!     data ioffp,spval/0,0.0/

      first=.false.

!  Szrat is the ratio of width to height of characters. Determined empirically
      szrat=.6

!  Set initial font to helvetica, 12 point
      ifntsz=12
      if(ientry.eq.999) then !Chopit was called, so use last font
        call setfnt(nfont)
      else
        call setfnt(20)
      endif

!  set factor to 1 for initialization, reset later if chopit called
      fac=1.
      call factor(fac)

!  Set initial lineweight to 0
      call setlw(0.)
!  Set initial grayscale to 0
      call setgry(0.)
!  Set initial rgb colors to black(0)
      call setcolr(0.,0.,0.)

      if(.not.portrait) then
        cmdstr='90 rotate 0 -8.5 inch translate '
        call filler
      endif

      call plot(.25*72./conver,.25*72./conver,-3)

      if(portrait) then
        xsh=.25*72./conver
        ysh=0.
      else
        xsh=0.
        ysh=.25*72./conver
      endif
      call plot(xsh,ysh,-3)

      return
      end
!*****PSINIT
      subroutine psinit(portrait)

!  RRD new date time function
      integer date_time (8)
      character (len=12) real_clock(3)
      character (len=3) month(12)
      
!  initializes plot for hp plotter
      logical first,portrait,prtrt,iopen
      character*132 cmdstr,curfnt
      character*80 fileout
      character tim*8,dat*9
      character*1 timer(8),dater(9)
      equivalence(timer(1),tim),(dater(1),dat)
      common/conre1/ioffp,spval
      common/plt1/cmdstr
      common/cnvcom/conver
      common/plt2/fac
      common/io/fileout,inew
      common/kkplot/szrat
      common/chpcom/ientry,prtrt
      common/fntcom/curfnt,ifntsz,nfont
      common/outcom/iunit
      common/pagcom/npage
!     data ioffp,spval/0,0.0/

      data month /'JAN','FEB','MAR','APR','MAY','JUN', &
                  'JUL','AUG','SEP','OCT','NOV','DEC'/
      ioffp=0
      spval=0.0

!  Set conversion factor (conver=72. for inches, conver=72./25.4 for mm, etc.)
!  conver
      conver=72.

      npage=1

      prtrt=portrait

      first=.true.

!  Use default name unless newdev has already been called (inew=999).
      if(inew.eq.0) then
        fileout='psplot.ps'
        inew=1
      else if(inew.eq.999)then
        inew=1
      endif

!  Open output file
      iunit=60
   10 inquire(iunit,opened=iopen)    ! get a non-existent unit
      if(iopen) then
        iunit=iunit+1
        goto 10
      endif
!  VMS
!     open(iunit,file=fileout,status='new',form='formatted', recl=132,
!    +carriagecontrol='list')
!  Everyone else
!     open(iunit,file=fileout,status='new',form='formatted', recl=132)
      open(iunit,file=fileout)

      cmdstr='%!PS-Adobe-2.0'              !Disable EPS RRD 6/27/01
!     cmdstr='%!PS-Adobe-2.0 EPSF-2.0'     #Use EPS format although it isn't.
                                           !Then just have to add BoundingBox
      call filler

      cmdstr= '%%Title: '//fileout(1:lenstr(fileout,80))

      call filler

! not all systems have date function
!RRD  dat='01-JUN-99'
!RRD  tim='hh:mm:ss'
!RRD  call time(tim)
!RRD  call date(dat)

      call date_and_time (real_clock(1),real_clock(2),real_clock(3),date_time)
      write(dat,'(i2.2,a1,a3,a1,i2.2)')date_time(3),'-',month(date_time(2)), &
            '-',(date_time(1)-2000)
      write(tim,'(i2.2,a1,i2.2,a1,i2.2)')date_time(5),':',date_time(6),      &
            ':',date_time(7)

      if(timer(1).eq.' ')timer(1)='0'
      if(dater(1).eq.' ')dater(1)='0'
      cmdstr= '%%CreationDate: '//DAT//' '//TIM

      call filler

      cmdstr= '%%Created with: PSPLOT PostScript Plotting Package'

      call filler
      cmdstr='%%Library Creator: Kevin E. Kohler <kevin@ocean.nova.edu>'
      call filler

!RRD - add string to indicate that page count comes at end
      cmdstr='%%Pages: (atend)'
      call filler

      cmdstr='%%EndComments'
      call filler
        cmdstr='/inch {72 mul} bind def'
      call filler

        cmdstr='/Ah {moveto lineto lineto stroke} def'
      call filler

        cmdstr='/Ar {moveto 2 copy lineto 4 -2 roll'
      call filler
        cmdstr='     moveto lineto lineto stroke } def'
      call filler

        cmdstr='/arcit {S /A2 exch def /A1 exch def /Rad exch def'
      call filler
      cmdstr='        /Yc exch def /Xc exch def'
      call filler
      cmdstr='        Xc Rad A1 cos mul add Yc Rad A1 sin mul add'
      call filler
      cmdstr='        moveto newpath'
      call filler
        cmdstr='        Xc Yc Rad A1 A2 arc stroke} def'
      call filler

        cmdstr='/C {/Rad exch def /Yc exch def /Xc exch def'
      call filler
      cmdstr='         Xc Yc Rad 0 360 arc closepath'    !closepath needed to
                                                         !avoid notch
      call filler                                        !with fat line width
        cmdstr='        } def'
      call filler

        cmdstr='/c0sf {closepath 0 setgray fill} def'
      call filler

        cmdstr='/cf {closepath fill} def'
      call filler

        cmdstr='/Cs {closepath stroke} def'
      call filler

        cmdstr='/Cln {newpath 3 1 roll'
      call filler
        cmdstr='      moveto {lineto} repeat clip newpath'
      call filler
        cmdstr='     } def'
      call filler

        cmdstr='/Cs {closepath stroke} def'
      call filler

        cmdstr='/Fb {newpath moveto '
      call filler
      cmdstr=' Dx 0 rlineto 0 Dy rlineto Dx neg 0 rlineto closepath'
      call filler
        cmdstr=' fill } def'
      call filler

        cmdstr='/Fbn { newpath 3 1 roll moveto {lineto} repeat'
      call filler
        cmdstr='       closepath fill } def'
      call filler

        cmdstr='/Fbnc { newpath 3 1 roll moveto'
      call filler
        cmdstr='       {lineto} repeat closepath fill } def'
      call filler

      cmdstr='/L /lineto load def'
      call filler

        cmdstr='/Lend {/Strlen exch stringwidth pop def} def'
      call filler

!  Define stringlength slightly increased for integrand placement
        cmdstr='/Lendi {/Strlen exch stringwidth pop 1.5 mul def} def'
      call filler

!  Define stringlength slightly increased for summation placement
        cmdstr='/Lends {/Strlen exch stringwidth pop 1.1 mul def} def'
      call filler

       cmdstr='/Lenssd '//'{/Strlenss exch stringwidth pop 3 mul 4 div def} def'
      call filler

        cmdstr='/LSM {2 copy lineto stroke moveto} def'
      call filler

        cmdstr='/lsm {Xp Yp lineto stroke mover} def'
      call filler

      cmdstr='/M /moveto load def'
      call filler

        cmdstr='/mover {Xp Yp moveto} def'
      call filler

        cmdstr='/Np {newpath} def'
      call filler

      cmdstr='/S /stroke load def'
      call filler

        cmdstr='/Sc {setrgbcolor} def'
      call filler

        cmdstr='/Sg {setgray} def'
      call filler

        cmdstr='/Setf {Curfnt exch scalefont setfont} def'
      call filler

        cmdstr='/SM {stroke moveto} def'
      call filler

        cmdstr='/sm {stroke mover} def'
      call filler

!       cmdstr='/Slw {inch setlinewidth} def'
      write(cmdstr,'(''/Slw {'',f7.4,'' mul setlinewidth} def'')') conver
      call filler

        cmdstr='/Slw0 {.24 setlinewidth} bind def'  !Minimum line width 300 dpi
      call filler

!  Add this for fun
      cmdstr= '%Line Breaking Procedure'
      call filler

      cmdstr='/TurnLineFL'
      call filler
        cmdstr='   { /T exch def /spacewidth space stringwidth pop def'
      call filler
      cmdstr='     /currentw 0 def /wordspace_count 0 def'
      call filler
      cmdstr='     /restart 0 def  /remainder T def'
      call filler
        cmdstr='     {remainder space search'
      call filler
        cmdstr='       {/nextword exch def pop'
      call filler
      cmdstr='        /remainder exch def'
      call filler
      cmdstr='        /nextwordwidth nextword stringwidth pop def'
      call filler
      cmdstr='        currentw nextwordwidth add lw gt'
      call filler
        cmdstr='        {T restart wordspace_count restart sub'
      call filler
      cmdstr='         getinterval showline'
      call filler
      cmdstr='         /restart wordspace_count def'
      call filler
      cmdstr='         /currentw nextwordwidth spacewidth add def'
      call filler
        cmdstr='        }'
      call filler
        cmdstr='        {/currentw currentw nextwordwidth add'
      call filler
      cmdstr='         spacewidth add def'
      call filler
        cmdstr='        } '
      call filler
      cmdstr='        ifelse'
      call filler
      cmdstr='        /wordspace_count wordspace_count'
      call filler
      cmdstr='        nextword length add 1 add def'
      call filler
        cmdstr='       }'
      call filler
        cmdstr='       {pop exit}'
      call filler
      cmdstr='       ifelse'
      call filler
        cmdstr='     } loop'
      call filler
      cmdstr='     /lrem remainder stringwidth pop def'
      call filler
      cmdstr='     currentw lrem add lw gt'
      call filler
        cmdstr='     {T restart wordspace_count restart sub '
      call filler
        cmdstr='      getinterval showline remainder showline}'
      call filler
        cmdstr='     {/lastchar T length def'
      call filler
      cmdstr='      T restart lastchar restart sub getinterval '
      call filler
        cmdstr='      lm y moveto show}'
      call filler
      cmdstr='     ifelse'
      call filler
        cmdstr='   } def'
      call filler

        cmdstr=' /parms {/y exch def /lm exch def /rm exch def'
      call filler
      cmdstr='         /leading exch def /pointsize exch def'
      call filler
      cmdstr='         /lw rm lm sub def'
      call filler
      cmdstr='         findfont pointsize scalefont setfont '
      call filler
        cmdstr='         /showline {lm y moveto show'
      call filler
        cmdstr='         /y y leading sub def} def'
      call filler
        cmdstr='         lm y moveto } def'
      call filler

        cmdstr='/Xposd {/Xpos exch def} def'
      call filler

        cmdstr='/Xposjd '//                                                    &
         ' {/Xpos exch Xpos exch Strlen mul sub def} def'
      call filler

        cmdstr=                                                                &
       '/xydef {/Xp Xpos def /Yp Ypos def} def'
      call filler

        cmdstr='%/Xypd {/Yp exch def /Xp exch def} def'
      call filler

        cmdstr='/Xypos0d {/Xpos0 Xpres def /Ypos0 Ypres def} def'
      call filler

        cmdstr='/Xyprset {dup '//                                              &
        '/Xpres exch cos Strlen mul Xpos add def'
      call filler
        cmdstr='              '//                                              &
       '/Ypres exch sin Strlen mul Ypos add def} def'
      call filler

        cmdstr='/Xyprset0 {dup '//                                             &
         '/Xpres exch cos Strlen mul Xpos0 add def'
      call filler
        cmdstr='               '//                                             &
         '/Ypres exch sin Strlen mul Ypos0 add def} def'
      call filler

        cmdstr='/Yposd {/Ypos exch def} def'
      call filler

        cmdstr='/Yposjd '//                                                    &
         ' {/Ypos exch Ypos exch Strlen mul sub def} def'
      call filler

      cmdstr='%%EndProlog'
      call filler
      cmdstr='/space ( ) def'
      call filler
      cmdstr='%%Page: 1 1'
      call filler
!rrd-start
      cmdstr='save'
      call filler
!rrd-stop

!  Szrat is the ratio of width to height of characters. Determined empirically.
      szrat=.6
!  Set initial font to helvetica, 12 point
      ifntsz=12
      call setfnt(20)
!  Set factor to 1 for initialization, reset later if chopit called
      fac=1.
      call factor(fac)

!  Set initial lineweight to 0
      call setlw(0.)
!  Set initial grayscale to 0
      call setgry(0.)
!  Set initial rgb colors to black(0)
      call setcolr(0.,0.,0.)

      if(.not.portrait) then
        cmdstr='90 rotate 0 -8.5 inch translate '
        call filler
      endif

      call plot(.25*72./conver,.25*72./conver,-3)

      if(portrait) then
        xsh=.25*72./conver
        ysh=0.
      else
        xsh=0.
        ysh=.25*72./conver
      endif
      call plot(xsh,ysh,-3)
      return
      end
!*****RECT
      subroutine rect(xx1,yy1,xx2,yy2,height)
!  this routines draws a rectangle with lower left
!  coordinate (x1,y1), lower right coordinate (x2,y2) and
!  height height
      parameter (pi=3.14159265358979)
      dimension xa(4),ya(4)

!      pi=4.*abs(atan(1.)) !
      dx=(xx2-xx1)
      dy=(yy2-yy1)
      plmi1=sign(1.,yy1-yy2)
      if(abs(dx).lt.1.e-5) then
        xinc=height*plmi1
        yinc=0.
      else
        xinc=height*cos(atan2(dy,dx)+pi/2.)
        yinc=height*sin(atan2(dy,dx)+pi/2.)
      endif
      xa(1)=xx1
      xa(2)=xx2
      ya(1)=yy1
      ya(2)=yy2
      xa(4)=xa(1)+xinc
      ya(4)=ya(1)+yinc
      xa(3)=xa(2)+xinc
      ya(3)=ya(2)+yinc
      call drwcrv(xa,ya,4,0.,.true.)
      return
      end
!*****RECTFILC
      subroutine rectfilc(xx1,yy1,xx2,yy2,height,red,green,blue)
!  this routines draws a rectangle with lower left
!  coordinate (x1,y1), lower right coordinate (x2,y2) and
!  height height. it then fills the rectangle with color(red,green,blue)
      parameter (pi=3.14159265358979)
      dimension xa(4),ya(4)
      common/colrcom/cred,cgreen,cblue,cgry

!  save current colors
      redsv=cred
      greensv=cgreen
      bluesv=cblue

!      pi=4.*abs(atan(1.)) !
      dx=(xx2-xx1)
      dy=(yy2-yy1)
      plmi1=sign(1.,yy1-yy2)
      if(abs(dx).lt.1.e-5) then
        xinc=height*plmi1
        yinc=0.
      else
        xinc=height*cos(atan2(dy,dx)+pi/2.)
        yinc=height*sin(atan2(dy,dx)+pi/2.)
      endif
      xa(1)=xx1
      xa(2)=xx2
      ya(1)=yy1
      ya(2)=yy2
      xa(4)=xa(1)+xinc
      ya(4)=ya(1)+yinc
      xa(3)=xa(2)+xinc
      ya(3)=ya(2)+yinc
      call setcolr(red,green,blue)
      call filrgnc(xa,ya,4)
      call setcolr(redsv,greensv,bluesv)
      return
      end
!*****RECTFILG
      subroutine rectfilg(xx1,yy1,xx2,yy2,height,gry)
!  this routines draws a rectangle with lower left
!  coordinate (x1,y1), lower right coordinate (x2,y2) and
!  height height. it then fills the rectangle with graylevel gry
      parameter (pi=3.14159265358979)
      dimension xa(4),ya(4)
      common/colrcom/cred,cgreen,cblue,cgry

!  save current graylevel
      cgrysv=cgry

!      pi=4.*abs(atan(1.)) !
      dx=(xx2-xx1)
      dy=(yy2-yy1)
      plmi1=sign(1.,yy1-yy2)
      if(abs(dx).lt.1.e-5) then
        xinc=height*plmi1
        yinc=0.
      else
        xinc=height*cos(atan2(dy,dx)+pi/2.)
        yinc=height*sin(atan2(dy,dx)+pi/2.)
      endif
      xa(1)=xx1
      xa(2)=xx2
      ya(1)=yy1
      ya(2)=yy2
      xa(4)=xa(1)+xinc
      ya(4)=ya(1)+yinc
      xa(3)=xa(2)+xinc
      ya(3)=ya(2)+yinc
      call setgry(gry)
      call filrgnc(xa,ya,4)
      call setgry(cgrysv)
      return
      end
!*****REORD
      subroutine reord (cl,ncl,c1,mark,nmg)
! This routine puts the major (labeled) levels in the beginning of cl
! and the minor (unlabeled) levels in end of cl.  The number of major
! levels is returned in mark.  c1 is used as a work space.  nmg is the
! number of minor gaps (one more than the number of minor levels between
! major levels).

      dimension cl(*),c1(*)

      nl = ncl
      if (nl.le. 4.or. nmg.le.1) go to 25
      nml = nmg-1
      if (nl .le. 10) nml = 1
!
! check for zero or other nice number for a major line
!
      nmlp1 = nml+1
      do i=1,nl
        isave = i
        if (cl(i) .eq. 0.) go to 10
      enddo
      l = nl/2
      l = alog10(abs(cl(l)))+1.
      q = 10.**l
      do j=1,3
        q = q/10.
        do i=1,nl
          isave = i
          if (amod(abs(cl(i)+1.e-9*cl(i))/q,float(nmlp1)).le.0.0001) goto 10
        enddo
      enddo
      isave = nl/2

! Put major levels in c1

   10 istart = mod(isave,nmlp1)
      if (istart .eq. 0) istart = nmlp1
      nmajl = 0
      do i=istart,nl,nmlp1
        nmajl = nmajl+1
        c1(nmajl) = cl(i)
      enddo
      mark = nmajl
      l = nmajl

! Put minor levels in c1

      if (istart .eq. 1) go to 15
      do i=2,istart
        isub = l+i-1
        c1(isub) = cl(i-1)
      enddo
   15 l = nmajl+istart-1
      do i=2,nmajl
        do j=1,nml
          l = l+1
          isub = istart+(i-2)*nmlp1+j
          c1(l) = cl(isub)
        enddo
      enddo
      nlml = nl-l
      if (l .eq. nl) go to 20
      do i=1,nlml
        l = l+1
        c1(l) = cl(l)
      enddo

! Put reordered array back in original place

   20 do i=1,nl
        cl(i) = c1(i)
      enddo
      return
   25 mark = nl
      return
      end
!*****ROTATE
      subroutine rotate(ang)

      character*132 cmdstr
      common/plt1/cmdstr

      cmdstr=' '
      write(cmdstr,'(f8.2,'' rotate'')')ang
      call filler
      return
      end
!*****RRECT
      subroutine rrect(x1,y1,width,height,rad,ang,fill)
!  This routine draws a rounded rectangle with lower left corner (x1,y1)
!  and corner radii rad
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)

      character*132 cmdstr
      common/plt1/cmdstr
      logical fill

!      pi=4.*abs(atan(1.)) !
!      rd=pi/180.           !
      cosa=cos(ang*rdpdeg)
      sina=sin(ang*rdpdeg)

      call plot(x1,y1,-3)
      xx1=0.
      xx2=xx1+width
      xx3=xx2
      xx4=xx1

      yy1=0.
      yy2=yy1
      yy3=yy2+height
      yy4=yy3

      if(ang.ne.0.) call rotate(ang)
      cmdstr='newpath'
      call filler
      call movet(.5*(xx1+xx2),yy1)
      call arcto(xx2,yy2,xx3,yy3,rad)
      call linet(xx2,.5*(yy2+yy3))
      call arcto(xx3,yy3,xx4,yy4,rad)
      call linet(.5*(xx3+xx4),yy3)
      call arcto(xx4,yy4,xx1,yy1,rad)
      call linet(xx1,.5*(yy4+yy1))
      call arcto(xx1,yy1,xx2,yy2,rad)

      if(fill) then
        cmdstr='closepath fill'
      else
        cmdstr='closepath stroke'
      endif
      call filler

      if(ang.ne.0.) call rotate(-ang)
      call plot(-x1,-y1,-3)

      return
      end
!*****SETCOLR
      subroutine setcolr(red,green,blue)
!  this routines sets the current color
!  red, green blue are the saturation ratios between 0 and 1
      character*132 cmdstr
      common/plt1/cmdstr
      common/colrcom/cred,cgreen,cblue,cgry

      r=red
      r=amin1(1.,r)
      r=amax1(0.,r)
      g=green
      g=amin1(1.,g)
      g=amax1(0.,g)
      b=blue
      b=amin1(1.,b)
      b=amax1(0.,b)

      cmdstr=' '
      write(cmdstr,'(3F7.3,'' Sc'')')r,g,b
      call filler
      cred=r
      cgreen=g
      cblue=b
      return
      end
!*****SETFNT
      subroutine setfnt(numfnt)
!  This routines changes the typeface of the current font
      character*132 cmdstr,scrc
      character*132 curfnt
      common/fntcom/curfnt,ifntsz,nfont
      character*40 fntnam(35)
      common/plt1/cmdstr
      data fntnam/'AvantGarde-Book','AvantGarde-BookOblique',                  &
      'AvantGarde-Demi','AvantGarde-DemiOblique','Bookman-Demi',               &
      'Bookman-DemiItalic','Bookman-Light','Bookman-LightItalic',              &
      'Courier-Bold','Courier-BoldOblique','Courier-Oblique', 'Courier',       &
      'Helvetica-Bold','Helvetica-BoldOblique', 'Helvetica-Narrow-Bold',       &
      'Helvetica-Narrow-BoldOblique', 'Helvetica-Narrow-Oblique',              &
      'Helvetica-Narrow', 'Helvetica-Oblique','Helvetica',                     &
      'NewCenturySchlbk-Bold','NewCenturySchlbk-BoldItalic',                   &
      'NewCenturySchlbk-Italic','NewCenturySchlbk-Roman',                      &
      'Palatino-Bold','Palatino-BoldItalic','Palatino-Italic',                 &
      'Palatino-Roman','Symbol','Times-Bold','Times-BoldItalic',               &
      'Times-Italic','Times-Roman','ZapfChancery-MediumItalic',                &
      'ZapfDingbats'/


      nfont=numfnt
      if(numfnt.lt.1.or.numfnt.gt.35) then
        print *,'Invalid font number encountered in **setfnt**'
        print *,'Using Helvetica default'
        nfont=20
      endif
      scrc=fntnam(nfont)
      cmdstr='/Curfnt /'//scrc(1:lenstr(scrc,132))//' findfont def'
      call filler
      write(cmdstr,'(i3,'' Setf'')')ifntsz
      call filler
      return
      end
!*****SETGRY
      subroutine setgry(gry)
!  This routines sets the current gray level
!  Gry is set to be between 0 and 1
      character*132 cmdstr
      common/plt1/cmdstr
      common/colrcom/cred,cgreen,cblue,cgry

      g=gry
      g=amin1(1.,g)
      g=amax1(0.,g)

      cmdstr=' '
      write(cmdstr,'(F7.3,'' Sg'')')g
      call filler
      cgry=g
      return
      end
!*****SETLW
      subroutine setlw(rlwi)
!  this routines sets the current linewidth
!  rlwi is linewidth in inches
      character*132 cmdstr
      common/plt1/cmdstr
      common/lcom/curwid

      if(abs(rlwi).lt.1.e-5) then  !0
        cmdstr='Slw0'
      else
        cmdstr=' '
        write(cmdstr,'(F7.3,'' Slw'')')rlwi
      endif
      call filler
      curwid=rlwi
      return
      end
!*****SIGMA
      subroutine sigma(xp,yp,size,ang,lower,nl,lupper,nu)
!  this routines plots the summation symbol sigma at position xp,yp.
!  size is the height of the (subsequent) integrand.
!  sigma character has been empirically enlarged.
!  lower,nl are the hollerith lower limit and number of chars
!  lupper,nu are the hollerith upper limit and number of chars
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      character*132 cmdstr
      character*132 curfnt,scrc
      common/fntcom/curfnt,ifntsz,nfont
      character*80 titlec
      character*1 bslash
      character * (*) lower,lupper
!      dimension lower(20),lupper(20)
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

!  Stroke previous paths before this write
      cmdstr='S'
      call filler

      bslash=char(92)

!      pi=4*abs(atan(1.)) !
      nchar=1
      ititle=345

!  Choose proper font height, using current font
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      iht=iht*1.5         !Enlarge sigma

      offset=iht*.13       !Define amount of offset to shift sigma "down"

      arg = ang * rdpdeg
      xoff=offset*sin(arg)
      yoff=-offset*cos(arg)

      call plot(xp+xoff/conver,yp+yoff/conver,-3)

!  Set angle
      if(ang.ne.0.) then
        cmdstr=' '
        write(cmdstr,'(F7.1,'' rotate'')') ang
        call filler
      endif

      limiht=iht*.2      !Limits height
      limiht=max0(limiht,1)
      ss=limiht/conver
      nfontsv=nfont
      if(nu.ne.0)then
        if(nu.eq.-999) then
          call setfnt(29)
          call keksym(0.,.85*iht/conver,1.4*ss,lupper,0.,nu,1)
          call setfnt(nfontsv)
        else
          call keksym(0.,.85*iht/conver,ss,lupper,0.,nu,1)
        endif
      endif

      if(nl.ne.0)then
        if(nl.eq.-999) then
          call setfnt(29)
          call keksym(0.,-1.7*ss,1.5*ss,lower,0.,nl,1)
          call setfnt(nfontsv)
        else
          call keksym(0.,-1.7*ss,ss,lower,0.,nl,1)
        endif
      endif
      call setfnt(nfontsv)

!  Reset angle
      if(ang.ne.0.) then
        cmdstr=' '
        write(cmdstr,'(F7.1,'' rotate'')') -ang
        call filler
      endif

      call plot(-(xp+xoff/conver),-(yp+yoff/conver),-3)

      cmdstr=' '
      write(cmdstr, '(''/Symbol findfont '',I3,'' scalefont setfont'')') IHT
      call filler
      write(titlec,'(a1,i10)')bslash,ititle
      call blkstp(titlec,80,titlec,numc)

      mchar=numc
      xpos=xp
      ypos=yp

      njus=1      !Force to be centered at xp,yp

      rsize=size
!  Character space height is 2.0 x char height
!  Character space width is 1.5 x char width
!  Actual string length is (nc-1)*1.5*char width + char width
      strlen=(rsize*szrat)*1.5*(mchar-1.)+rsize*szrat

      if(xpos.eq.999.) then
        write(cmdstr,'(''/Xpos Xpres '',f8.2,''add def'')')xoff
      else
        cmdstr=' '
        xnew=xp*conver+xoff
        write(cmdstr,'(f8.2,'' Xposd'')')xnew
      endif
      call filler

      if(ypos.eq.999.) then
        write(cmdstr,'(''/Ypos Ypres '',f8.2,''add def'')')yoff
      else
        cmdstr=' '
        ynew=yp*conver+yoff
        write(cmdstr,'(f8.2,'' Yposd'')')ynew
      endif
      call filler

      if(njus.ne.0.and.njus.ne.1.and.njus.ne.2) then
        print 110, njus
  110 format(1x,'Incorrect justification code ',i5,'found in ',                &
      'SIGMA, zero used')
        njus=0
      endif

!  Strlen has already been "factored" by the choice of font height
!  Since it will eventually be factored again, we must divide by
!  factor now.
      cmdstr='('//titlec(1:mchar)//') Lend'  !Use actual length for centering
      xarg=cos(arg)*njus/2.
      yarg=sin(arg)*njus/2.
      if(xarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Xposjd'')')xarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      if(yarg.ne.0.) then
        scrc=' '
        write(scrc,'(f7.3,'' Yposjd'')')yarg
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      call filler

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      scrc='('//Titlec(1:mchar)//') show'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Start next char at .5 char width away
      argdeg = arg / rdpdeg

!  Reset Xpos &Ypos from sigma offset
      cmdstr=' '
      write(cmdstr, '(''/Xpos Xpos '',f8.2,'' sub def'')')xoff

      call filler
      write(cmdstr,'(''  /Ypos Ypos '',f8.2,'' sub def'')')yoff
      call filler

      cmdstr='('//titlec(1:mchar)//') Lends'  !Shift next position slightly
      call filler

      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset'')')argdeg
      call filler

!  Reset font
      call setfnt(nfont)

      return
      end
!*****SLDCRV
      subroutine sldcrv(x,y,npts,thk)
! this routine draws a solid curve connected data points in arrays x andy.
! x and y arrays are in real inches.
! npts is the number of points to connect
! thk is the thickness of the line(inches).  if zero, a single line is drawn.
! the line is centered on the line joining the two points
      dimension x(npts),y(npts)

      call drwcrv(x,y,npts,thk,.false.)

      return
      end
!*****SLDLIN
      subroutine sldlin(x1,y1,x2,y2,thk)
!  this routine draws a solid line between points (x1,y1) and (x2,y2).
! thk is the thickness of the line(inches).  if zero, a single line is drawn.
! the line is centered on the line joining the two points
      common/lcom/curwid

      iset=0
      if(thk.eq.0.) then    !Use current line width
      else if(thk.ne.curwid) then
        iset=1
        curwids=curwid
        call setlw(thk)
      endif
      call plot(x1,y1,3)
      call plot(x2,y2,2)

      if(iset.eq.1) call setlw(curwids)  !Reset linewidth to original

      return
      end
!*****SQRSGN
      subroutine sqrsgn(xx,yy,hgt,rlen)
!  xx is the x-coordinate of the radical sign overbar break
!  yy is the y-coordinate of the bottom of the radical sign
!  hgt is height of radical sign
!  rlen is length of overbar

      dimension xa(8),ya(8)
      xa(1)=xx
      xa(2)=xx+hgt*.24
      xa(3)=xx+hgt*.5
      xa(4)=xx+hgt*1.
      xa(5)=xx+hgt*1.
      xa(6)=xx+hgt*.46
      xa(7)=xx+hgt*.12
      xa(8)=xx

      do 10 i=1,8
   10 xa(i)=xa(i)-hgt

      ya(1)=yy+hgt*.48
      ya(2)=yy+hgt*.6
      ya(3)=yy+hgt*.2
      ya(4)=yy+hgt*1.
      ya(5)=yy+hgt*.9
      ya(6)=yy
      ya(7)=yy+hgt*.44
      ya(8)=yy+hgt*.38

      call filrgn(xa,ya,8,0.)
!  Draw overbar
      xa(1)=xa(4)
      xa(2)=xa(1)+rlen
      xa(3)=xa(2)
      xa(4)=xa(5)
      ya(1)=ya(4)
      ya(2)=ya(1)
      ya(3)=ya(5)
      ya(4)=ya(3)
      call filrgn(xa,ya,4,0.)
      return
      end
!*****SQUARE
      subroutine square(xc,yc,side)
      dimension xa(4), ya(4)
      xa(1)=xc-side/2.
      ya(1)=yc-side/2.
      xa(2)=xa(1)+side
      ya(2)=ya(1)
      xa(3)=xa(2)
      ya(3)=ya(2)+side
      xa(4)=xa(1)
      ya(4)=ya(3)
      call drwcrv(xa,ya,4,0.,.true.)
      return
      end
!*****STLINE
      subroutine stline (z,ll,mm,nn,conv,ilabl,scal)
      dimension       z(ll,nn)
!
! This routine finds the beginnings of all contour lines at level conv.
! first the edges are searched for lines intersecting the edge (open
! lines) then the interior is searched for lines which do not intersect
! the edge (closed lines).  beginnings are stored in ir to prevent re-
! tracing of lines.  if ir is filled, the search is stopped for this
! conv.
!
      common /conre2/ ix ,iy ,idx ,idy , is ,iss ,np ,cv , inx(8) ,iny         &
      (8) ,ir(10000) ,nr

      common/contyp/ispcon,idsh
!  Changed 9/24/96 kek
!     data inx(1),inx(2),inx(3),inx(4),inx(5),inx(6),inx(7),inx(8)/ - 1,
!    +- 1,  0,  1,  1,  1,  0, - 1/
!     data iny(1),iny(2),iny(3),iny(4),iny(5),iny(6),iny(7),iny(8)/ 0,
!    +1, 1, 1, 0, - 1, - 1, - 1/
!     data nr/10000/

      lc16(k) = k*65536
!     lc16(k) = k*'200000'O    !VMS
!     lc16(k) = k*8#200000    !MS

      inx(1)=-1
      inx(2)=-1
      inx(3)=0
      inx(4)=1
      inx(5)=1
      inx(6)=1
      inx(7)=0
      inx(8)=-1

      iny(1)=0
      iny(2)=1
      iny(3)=1
      iny(4)=1
      iny(5)=0
      iny(6)=-1
      iny(7)=-1
      iny(8)=-1

      nr=10000

      l = ll
      m = mm
      n = nn
      cv = conv
      np = 0
      iss = 0
      do 15 ip1=2,m
        i = ip1-1
        if (z(i,1).ge.cv .or. z(ip1,1).lt.cv) go to 10
        ix = ip1
        iy = 1
        idx = -1
        idy = 0
        is = 1
        if(idsh.ne.0) then
          call drlindsh(z,l,m,n,ilabl,scal)
        else
          call drlin2 (z,l,m,n,ilabl,scal)
        endif
   10   if (z(ip1,n).ge.cv .or. z(i,n).lt.cv) go to 15
        ix = i
        iy = n
        idx = 1
        idy = 0
        is = 5
        if(idsh.ne.0) then
          call drlindsh(z,l,m,n,ilabl,scal)
        else
          call drlin2(z,l,m,n,ilabl,scal)
        endif
   15 continue
      do 25 jp1=2,n
        j = jp1-1
        if (z(m,j).ge.cv .or. z(m,jp1).lt.cv) go to 20
        ix = m
        iy = jp1
        idx = 0
        idy = -1
        is = 7
        if(idsh.ne.0) then
          call drlindsh(z,l,m,n,ilabl,scal)
        else
          call drlin2 (z,l,m,n,ilabl,scal)
        endif
   20   if (z(1,jp1).ge.cv .or. z(1,j).lt.cv) go to 25
        ix = 1
        iy = j
        idx = 0
        idy = 1
        is = 3
        if(idsh.ne.0) then
          call drlindsh(z,l,m,n,ilabl,scal)
        else
          call drlin2 (z,l,m,n,ilabl,scal)
        endif
   25 continue
      iss = 1
      do 45 jp1=3,n
        j = jp1-1
        do 40 ip1=2,m
          i = ip1-1
          if (z(i,j).ge.cv .or. z(ip1,j).lt.cv) go to 40
          ixy = lc16(ip1)+j
          if (np .eq. 0) go to 35
          do 30 k=1,np
            if (ir(k) .eq. ixy) go to 40
   30     continue
   35     np = np+1
          if (np .gt. nr) return
          ir(np) = ixy
          ix = ip1
          iy = j
          idx = -1
          idy = 0
          is = 1
          if(idsh.ne.0) then
            call drlindsh(z,l,m,n,ilabl,scal)
          else
            call drlin2 (z,l,m,n,ilabl,scal)
          endif
   40   continue
   45 continue
      return
      end

!*****STLINEDP
      subroutine stlinedp (z,ll,mm,nn,conv,ilabl,scal)
      double precision z(ll,nn)
!
! this routine finds the beginnings of all contour lines at level conv.
! first the edges are searched for lines intersecting the edge (open
! lines) then the interior is searched for lines which do not intersect
! the edge (closed lines).  beginnings are stored in ir to prevent re-
! tracing of lines.  if ir is filled, the search is stopped for this
! conv.
!
      common /conre2/ix,iy,idx,idy,is,iss,np,cv,inx(8),iny(8),ir(10000),nr

      common/contyp/ispcon,idsh
      double precision cvdp
      common/conre3/cvdp
!  Changed 9/24/96 kek
!     data inx(1),inx(2),inx(3),inx(4),inx(5),inx(6),inx(7),inx(8)/ - 1,
!    +- 1,  0,  1,  1,  1,  0, - 1/
!     data iny(1),iny(2),iny(3),iny(4),iny(5),iny(6),iny(7),iny(8)/ 0,
!    +1, 1, 1, 0, - 1, - 1, - 1/
!     data nr/10000/

      lc16(k) = k*65536
!     lc16(k) = k*'200000'O    !VMS
!     lc16(k) = k*8#200000    !MS

      inx(1)=-1
      inx(2)=-1
      inx(3)=0
      inx(4)=1
      inx(5)=1
      inx(6)=1
      inx(7)=0
      inx(8)=-1

      iny(1)=0
      iny(2)=1
      iny(3)=1
      iny(4)=1
      iny(5)=0
      iny(6)=-1
      iny(7)=-1
      iny(8)=-1

      nr=10000

      l = ll
      m = mm
      n = nn
      cvdp = conv
      np = 0
      iss = 0
      do 15 ip1=2,m
        i = ip1-1
        if (z(i,1).ge.cvdp .or. z(ip1,1).lt.cvdp) go to 10
        ix = ip1
        iy = 1
        idx = -1
        idy = 0
        is = 1
!  Changed 9/20/96 kek
!       if(idsh.ne.0) then
!         call drlindsh(z,l,m,n,ilabl,scal)
!       else
!         call drlin2dp (z,l,m,n,ilabl,scal)
!       endif
        call drlin2dp (z,l,m,n,ilabl,scal)
   10   if (z(ip1,n).ge.cvdp .or. z(i,n).lt.cvdp) go to 15
        ix = i
        iy = n
        idx = 1
        idy = 0
        is = 5
!  Changed 9/20/96 kek
!       if(idsh.ne.0) then
!         call drlindsh(z,l,m,n,ilabl,scal)
!       else
!         call drlin2dp (z,l,m,n,ilabl,scal)
!       endif
        call drlin2dp (z,l,m,n,ilabl,scal)
   15 continue
      do 25 jp1=2,n
        j = jp1-1
        if (z(m,j).ge.cvdp .or. z(m,jp1).lt.cvdp) go to 20
        ix = m
        iy = jp1
        idx = 0
        idy = -1
        is = 7
!  Changed 9/20/96 kek
!       if(idsh.ne.0) then
!         call drlindsh(z,l,m,n,ilabl,scal)
!       else
!         call drlin2dp (z,l,m,n,ilabl,scal)
!       endif
        call drlin2dp (z,l,m,n,ilabl,scal)
   20   if (z(1,jp1).ge.cvdp .or. z(1,j).lt.cvdp) go to 25
        ix = 1
        iy = j
        idx = 0
        idy = 1
        is = 3
!  Changed 9/20/96 kek
!       if(idsh.ne.0) then
!         call drlindsh(z,l,m,n,ilabl,scal)
!       else
!         call drlin2dp (z,l,m,n,ilabl,scal)
!       endif
        call drlin2dp (z,l,m,n,ilabl,scal)
   25 continue
      iss = 1
      do 45 jp1=3,n
        j = jp1-1
        do 40 ip1=2,m
          i = ip1-1
          if (z(i,j).ge.cvdp .or. z(ip1,j).lt.cvdp) go to 40
          ixy = lc16(ip1)+j
          if (np .eq. 0) go to 35
          do 30 k=1,np
            if (ir(k) .eq. ixy) go to 40
   30     continue
   35     np = np+1
          if (np .gt. nr) return
          ir(np) = ixy
          ix = ip1
          iy = j
          idx = -1
          idy = 0
          is = 1
!  Changed 9/20/96 kek
!         if(idsh.ne.0) then
!           call drlindsh(z,l,m,n,ilabl,scal)
!         else
!           call drlin2dp (z,l,m,n,ilabl,scal)
!         endif
          call drlin2dp (z,l,m,n,ilabl,scal)
   40   continue
   45 continue
      return
      end
!*****STROKE
      subroutine stroke
      character*132 cmdstr
      common/plt1/cmdstr

      cmdstr='stroke'
      call filler
      return
      end
!*****SUBBER
      subroutine subber(ititle,nchar,size,ang)
!  ititle is subscript character(s)
!  size is height of subscripted variable
      parameter (pi=3.14159265358979, rdpdeg=pi/180.)
      character * (*) ititle
!      dimension ititle(20)    !hollerith
      character*132 cmdstr,curfnt,scrc,titlec*80
      character*1 bslash
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver

      bslash=char(92)
!      pi=4.*abs(atan(1.)) !
      angrad = ang * rdpdeg
      ss=3.*size/4.
      down=ss/2.
      xshift=sin(angrad)*down
      yshift=cos(angrad)*down

      cmdstr='Xypos0d'
      call filler

!      write(titlec,'(20a4)')ititle
      titlec=ititle

      cmdstr=' '
      write(cmdstr,'(''/Xpos Xpres '',f8.2,'' add def '')')xshift*conver
      call filler
      cmdstr=' '
      write(cmdstr,'(''/Ypos Ypres '',F8.2,'' sub def '')')yshift*conver
      call filler

!  Choose font proper height using currently set font
      mchar=iabs(nchar)

!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

      ireset=0
      if(titlec(1:mchar).ne.char(39)) then
        if(iht.ne.ifntsz) then
          cmdstr=' '
          write(cmdstr,'(I3,'' Setf'')')iht
          call filler
          ifntsz=iht
        endif
      else        !Use prime, not apostrophe (do not change ifntsz)
        ireset=1
        cmdstr=' '
        write(cmdstr,'(''/Symbol findfont '',I3,'' scalefont setfont'')')iht
        call filler
      endif

      rsize=size

      if(titlec(1:mchar).ne.char(39)) then
        cmdstr='('//titlec(1:mchar)//') Lend'
      else
        cmdstr='('//bslash//'242) Lend'
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif

      if(titlec(1:mchar).ne.char(39)) then
        scrc='('//Titlec(1:mchar)//') show'
      else
        scrc='('//bslash//'242) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
      endif
      call filler

!  Compute next character position
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset0'')')ang
      call filler

!  Reset font if we drew a prime
      if(ireset.eq.1) call setfnt(nfont)

      return
      end
!*****SUBBERSP
      subroutine subbersp(nset,nfnt,ititle,nchr,size,ang)
!  this routines allows subscripts of "special" characters
!  ititle is superscript character(s)
!  size is height of superscripted variable

!  nset is the number of different font sets needed
!  nfnt are the font numbers
!  ititle holds both the octal codes for fonts 29 and 35, and the
!  hollerith characters for other fonts
!  nchr holds the number of characters in each set (usually one, but can
!  be gt 1 if font is not symbol or dingbats)
      parameter (pi=3.14159265358979,rdpdeg=pi/180.)
      dimension nfnt(nset),nchr(nset)
      dimension ititle(20,nset)    !hollerith
      character*132 cmdstr,curfnt,scrc
      character*80 titleb,titlec
      character*1 bslash
      logical spchar
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver

!      pi=4.*abs(atan(1.)) !
      bslash=char(92)

!  Save current font
      nfsav=nfont

      angrad = ang * rdpdeg
      ss=3.*size/4.
      down=ss/2.
      xshift=sin(angrad)*down
      yshift=cos(angrad)*down

!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

!  Loop on character sets
      do 10 nc=1,nset

        cmdstr='Xypos0d'
        call filler
        cmdstr=' '
        write(cmdstr,'(''/Xpos Xpres '',f8.2,'' add def '')') xshift*conver

        call filler
        cmdstr=' '
        write(cmdstr,'(''/Ypos Ypres '',F8.2,'' sub def '')') yshift*conver

        call filler

!  Set the proper font
        call setfnt(nfnt(nc))

!  Set font height
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')iht
        call filler

        if(nfnt(nc).eq.29.or.nfnt(nc).eq.35) then
          spchar=.true.
        else
          spchar=.false.
        endif

        mchar=iabs(nchr(nc))

        if(spchar) then   !octal code
          write(titlec,'(a1,i4)')bslash,ititle(1,nc)
          call blkstp(titlec,80,titlec,mchar)
          cmdstr='('//titlec(1:mchar)//') Lend'
        ELSE
          mchar=nchr(nc)
          write(titlec,'(20a4)')(ititle(ii,nc),ii=1,20)

!       Check if titlec contains ( or ) or \.  These characters must be treated
!       specially by preceding them with a "\".  Do this to ( and ) even though
!       they might be balanced, i.e. () within a string, which can be treated
!       normally.

          titleb=titlec
          numc=0
          do m=1,mchar
          if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')'.or.titleb(m:m).eq.bslash)&
          then
              numc=numc+1
              titlec(numc:numc)=bslash
            endif
            numc=numc+1
            titlec(numc:numc)=titleb(m:m)
          enddo
          mchar=numc
          cmdstr='('//TITLEC(1:MCHAR)//') Lend'
        endif

!  Move pen to proper coordinates
        scrc='xydef mover'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif

        scrc='('//Titlec(1:mchar)//') show'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') -ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif
        call filler

!  Compute next character position
        cmdstr=' '
        write(cmdstr,'(f6.1,'' Xyprset0'')')ang
        call filler

   10 continue

!  Reset font
      call setfnt(nfsav)

      return
      end
!*****SUBSUP
      subroutine subsup(isub,nsub,isup,nsup,size,ang)
!  isub is subscript character(s); nsub is no. of subscript characters
!  isup is subscript character(s); nsup is no. of subscript characters
!  size is height of subscripted variable
      parameter (pi=3.14159265358979, rdpdeg=pi/180.)
!rrd  dimension isub(20),isup(20)    !hollerith
      character*80 isub,isup
      character*132 cmdstr,curfnt,scrc,titlec*80
      character*1 bslash
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver

      bslash=char(92)
!      pi=4.*abs(atan(1.)) !
      angrad = ang * rdpdeg
      ss=3.*size/4.

!  Subscript part
      down=ss/2.
      xshift=sin(angrad)*down
      yshift=cos(angrad)*down
      cmdstr='Xypos0d'
      call filler

!rrd  write(titlec,'(20a4)')isub
	  titlec=isub

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Xpres '',f8.2,'' add def '')')xshift*conver
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Ypres '',F8.2,'' sub def '')')yshift*conver
      call filler

!  Choose proper height
      mchar=iabs(nsub)
!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

      ireset=0  !apostrophe
      if(titlec(1:mchar).ne.char(39)) then
        if(iht.ne.ifntsz) then
          cmdstr=' '
          write(cmdstr,'(I3,'' Setf'')')iht
          call filler
          ifntsz=iht
        endif
      else        !Use prime, not apostrophe (do not change ifntsz)
        ireset=1
        cmdstr=' '
        write(cmdstr,                                                          &
        '(''/Symbol findfont '',I3,'' scalefont setfont'')')iht
        call filler
      endif

      rsize=size

      if(titlec(1:mchar).ne.char(39)) then
        cmdstr='('//titlec(1:mchar)//') Lend'
      else
        cmdstr='('//bslash//'242) Lend'
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if(titlec(1:mchar).ne.char(39)) then
        scrc='('//Titlec(1:mchar)//') show'
      else
        scrc='('//bslash//'242) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Superscript part
      write(titlec,'(20a4)')isup

      up=4.*ss/5.
      xshift=sin(angrad)*up
      yshift=cos(angrad)*up

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Xpres '',f8.2,'' sub def '')')xshift*conver
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Ypres '',F8.2,'' add def '')')yshift*conver
      call filler

!  Choose proper height
      mchar=iabs(nsup)

!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

      if(titlec(1:mchar).ne.char(39)) then
        if(iht.ne.ifntsz) then
          cmdstr=' '
          write(cmdstr,'(I3,'' Setf'')')iht
          call filler
          ifntsz=iht
        endif
      else        !Use prime, not apostrophe (do not change ifntsz)
        ireset=1
        cmdstr=' '
        write(cmdstr,                                                          &
        '(''/Symbol findfont '',I3,'' scalefont setfont'')')iht
        call filler
      endif

      rsize=size

      if(titlec(1:mchar).ne.char(39)) then
        cmdstr='('//titlec(1:mchar)//') Lend'
      else
        cmdstr='('//bslash//'242) Lend'
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if(titlec(1:mchar).ne.char(39)) then
        scrc='('//Titlec(1:mchar)//') show'
      else
        scrc='('//bslash//'242) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Reset angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Compute next character position
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset0'')')ang
      call filler

!  Reset font if we drew a prime
      if(ireset.eq.1)call setfnt(nfont)

      return
      end
!*****SUBSUPSP
      subroutine subsupsp(nsub,nfntsb,ititlesb,nchrsb,nsup,nfntsp,             &
      ititlesp,nchrsp, size,ang)

!  this routines allows subscripts and superscripts of "special" characters
!  ititle is superscript character(s)
!  size is height of superscripted variable

!  nsub,nsup are the number of different font sets needed for subs and supers
!  nfntsb,nfntsp are the font numbers
!  ititlesb,sp holds both the octal codes for fonts 29 and 35, and the
!  hollerith characters for other fonts
!  nchrsb,sp holds the number of characters in each set (usually one, but can
!  be gt 1 if font is not symbol or dingbats)
      parameter (pi=3.14159265358979, rdpdeg=pi/180.)
      dimension nfntsb(nsub),nchrsb(nsub)
      dimension ititlesb(20,nsub)    !hollerith
      dimension nfntsp(nsup),nchrsp(nsup)
      dimension ititlesp(20,nsup)    !hollerith
      character*132 cmdstr,curfnt,scrc
      character*80 titleb,titlec
      character*1 bslash
      logical spchar
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver

!      pi=4.*abs(atan(1.)) !
      bslash=char(92)

!  Save current font
      nfsav=nfont

!  Subscript part
      angrad = ang * rdpdeg
      ss=3.*size/4.
      down=ss/2.
      xshift=sin(angrad)*down
      yshift=cos(angrad)*down

!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

!  Loop on character sets
      do 10 nc=1,nsub

        cmdstr='Xypos0d'
        call filler

        if(nc.eq.1) then  !save this coordinate for the superscript part
          cmdstr='/Xposv Xpos0 def  /Yposv Ypos0 def'
          call filler
        endif

        cmdstr=' '
        write(cmdstr,'(''/Xpos Xpres '',f8.2,'' add def '')')                  &
              xshift*conver

        call filler
        cmdstr=' '
        write(cmdstr,'(''/Ypos Ypres '',F8.2,'' sub def '')')                  &
              yshift*conver

        call filler

!  Set the proper font
        call setfnt(nfntsb(nc))
!  Set font height
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')IHT
        call filler

        if(nfntsb(nc).eq.29.or.nfntsb(nc).eq.35) then
          spchar=.true.
        else
          spchar=.false.
        endif

        mchar=iabs(nchrsb(nc))

        if(spchar) then   !octal code
          write(titlec,'(a1,i4)')bslash,ititlesb(1,nc)
          call blkstp(titlec,80,titlec,mchar)
          cmdstr='('//titlec(1:mchar)//') Lend'
        else
          mchar=nchrsb(nc)
          write(TITLEC,'(20A4)')(ITITLESB(II,NC),II=1,20)

!       Check if titlec contains ( or ) or \.  These characters must be treated
!       specially by preceding them with a "\".  Do this to ( and ) even though
!       they might be balanced, i.e. () within a string, which can be treated
!       normally.

          titleb=titlec
          numc=0
          do m=1,mchar
            if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')'.or.titleb(m:m)        &
            .eq.bslash) then
              numc=numc+1
              titlec(numc:numc)=bslash
            endif
            numc=numc+1
            titlec(numc:numc)=titleb(m:m)
          enddo
          mchar=numc
          cmdstr='('//TITLEC(1:MCHAR)//') Lend'
        endif

!  Move pen to proper coordinates
        scrc='xydef mover'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif

        scrc='('//Titlec(1:mchar)//') show'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') -ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif
        call filler

!  Compute next character position
        cmdstr=' '
        write(cmdstr,'(f6.1,'' Xyprset0'')')ang
        call filler

   10 continue

!  Superscript part
      up=4.*ss/5.
      xshift=sin(angrad)*up
      yshift=cos(angrad)*up

!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

!  Reset xpres,ypres from to initial value
      cmdstr='/Xpres Xposv def /Ypres Yposv def'
      call filler

!  Loop on character sets
      do 15 nc=1,nsup

        cmdstr='Xypos0d'
        call filler

        cmdstr=' '
        write(cmdstr,'(''/Xpos Xpres '',f8.2,'' sub def '')') xshift*conver

        call filler
        cmdstr=' '
        write(cmdstr,'(''/Ypos Ypres '',F8.2,'' add def '')') yshift*conver

        call filler

!  Set the proper font
        call setfnt(nfntsp(nc))
!  Set font height
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')iht
        call filler

        if(nfntsp(nc).eq.29.or.nfntsp(nc).eq.35) then
          spchar=.true.
        else
          spchar=.false.
        endif

        mchar=iabs(nchrsp(nc))

        if(spchar) then   !octal code
          write(titlec,'(a1,i4)')bslash,ititlesp(1,nc)
          call blkstp(titlec,80,titlec,mchar)
          cmdstr='('//titlec(1:mchar)//') Lend'
        else
          mchar=nchrsp(nc)
          write(titlec,'(20a4)')(ititlesp(ii,nc),ii=1,20)

!       Check if titlec contains ( or ) or \.  These characters must be treated
!       specially by preceding them with a "\".  Do this to ( and ) even though
!       they might be balanced, i.e. () within a string, which can be treated
!       normally.

          titleb=titlec
          numc=0
          do m=1,mchar
            if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m)        &
            .eq.bslash) then
              numc=numc+1
              titlec(numc:numc)=bslash
            endif
            numc=numc+1
            titlec(numc:numc)=titleb(m:m)
          enddo
          mchar=numc
          cmdstr='('//titlec(1:mchar)//') Lend'
        endif

!  Move pen to proper coordinates
        scrc='xydef mover'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif

        scrc='('//Titlec(1:mchar)//') show'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Set angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') -ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif
        call filler

!  Compute next character position
        cmdstr=' '
        write(cmdstr,'(f6.1,'' Xyprset0'')')ang
        call filler

   15 continue

!  Reset font
      call setfnt(nfsav)

      return
      end
!*****SUPER
      subroutine super(ititle,ncharr,size,ang)
!  ititle is superscript character(s)
!  size is height of superscripted variable
      parameter (pi=3.14159265358979, rdpdeg=pi/180.)
!rrd  dimension ititle(20)    !hollerith
      character*80 ititle
      character*132 cmdstr,curfnt,scrc,titlec*80
      character*1 bslash
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver
      bslash=char(92)

!     write(titlec,'(20a4)')ititle
      write(titlec,'(a)')ititle

!  If ncharr is lt 0, we're calling from oversup and we should
!  lower superscript slightly

      nchar=abs(ncharr)

!      pi=4.*abs(atan(1.)) !
      angrad = ang * rdpdeg
      ss=3.*size/4.
      if(ncharr.ge.0.) then
        up=4.*ss/5.
      else
        up=1.*ss/2.
      endif
      xshift=sin(angrad)*up
      yshift=cos(angrad)*up
      cmdstr='Xypos0d'
      call filler

      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Xpos Xpres '',f8.2,'' sub def '')')XSHIFT*conver
      call filler
      cmdstr=' '
      write(cmdstr,                                                            &
      '(''/Ypos Ypres '',F8.2,'' add def '')')YSHIFT*conver
      call filler

!  Choose proper height
      mchar=iabs(nchar)
!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

      ireset=0
      if(titlec(1:mchar).ne.char(39)) then
        if(iht.ne.ifntsz) then
          cmdstr=' '
          write(cmdstr,'(I3,'' Setf'')')IHT
          call filler
          ifntsz=iht
        endif
      else        !Use prime, not apostrophe (do not change ifntsz)
        ireset=1
        cmdstr=' '
        write(cmdstr,                                                          &
        '(''/Symbol findfont '',I3,'' scalefont setfont'')')iht
        call filler
      endif

      rsize=size

      if(titlec(1:mchar).ne.char(39)) then   !DRAW PRIME, NOT APOSTROPHE
        cmdstr='('//titlec(1:mchar)//') Lend'
      else
        cmdstr='('//bslash//'242) Lend'
      endif

!  Move pen to proper coordinates
      scrc='xydef mover'
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif

      if(titlec(1:mchar).ne.char(39)) then
        scrc='('//Titlec(1:mchar)//') show'
      else
        scrc='('//bslash//'242) show'
      endif
      cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132)       &
      )

!  Set angle
      if(ang.ne.0.) then
        scrc=' '
        write(scrc,'(F7.1,'' rotate'')') -ang
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))
      endif
      call filler

!  Compute next character position
      cmdstr=' '
      write(cmdstr,'(f6.1,'' Xyprset0'')')ang
      call filler

!  Reset font if we drew a prime
      if(ireset.eq.1) call setfnt(nfont)

      return
      end
!*****SUPERSP
      subroutine supersp(nset,nfnt,ititle,nchr,size,ang)
!  this routines allows superscripts of "special" characters
!  ititle is superscript character(s)
!  size is height of superscripted variable

!  nset is the number of different font sets needed
!  nfnt are the font numbers
!  ititle holds both the octal codes for fonts 29 and 35, and the
!  hollerith characters for other fonts
!  nchr holds the number of characters in each set (usually one, but can
!  be gt 1 if font is not symbol or dingbats)
      parameter (pi=3.14159265358979, rdpdeg=pi/180.)
      dimension nfnt(nset),nchr(nset)
      dimension ititle(20,nset)    !hollerith
      character*132 cmdstr,curfnt,scrc
      character*80 titleb,titlec
      character*1 bslash
      logical spchar
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver

!      pi=4.*abs(atan(1.)) !
      bslash=char(92)

!  Save current font
      nfsav=nfont

      angrad = ang * rdpdeg
      ss=3.*size/4.
      UP=4.*ss/5.
      xshift=sin(angrad)*up
      yshift=cos(angrad)*up


!  Set character size
      iht=ss*conver/.6     !.6 factor is empirical

!  Loop on character sets
      do 10 nc=1,nset

        cmdstr='Xypos0d'
        call filler
        cmdstr=' '
        write(cmdstr,'(''/Xpos Xpres '',f8.2,'' sub def '')')                  &
              XSHIFT*conver

        call filler
        cmdstr=' '
        write(cmdstr,'(''/Ypos Ypres '',F8.2,'' add def '')')                  &
              YSHIFT*conver

        call filler

!  Set the proper font
        call setfnt(nfnt(nc))
!  Set font height
        cmdstr=' '
        write(cmdstr,'(I3,'' Setf'')')iht
        call filler

        if(nfnt(nc).eq.29.or.nfnt(nc).eq.35) then
          spchar=.true.
        else
          spchar=.false.
        endif

        mchar=iabs(nchr(nc))

        if(spchar) then   !octal code
          write(titlec,'(a1,i4)')bslash,ititle(1,nc)
          call blkstp(titlec,80,titlec,mchar)
          cmdstr='('//titlec(1:mchar)//') Lend'
        else
          mchar=nchr(nc)
          write(titlec,'(20a4)')(ititle(ii,nc),ii=1,20)

!       Check if titlec contains ( or ) or \.  These characters must be treated
!       specially by preceding them with a "\".  Do this to ( and ) even though
!       they might be balanced, i.e. () within a string, which can be treated
!       normally.

          titleb=titlec
          numc=0
          do m=1,mchar
            if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m)        &
            .eq.bslash) then
              numc=numc+1
              titlec(numc:numc)=bslash
            endif
            numc=numc+1
            titlec(numc:numc)=titleb(m:m)
          enddo
          mchar=numc
          cmdstr='('//TITLEC(1:MCHAR)//') Lend'
        endif

!  Move pen to proper coordinates
        scrc='xydef mover'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr               &
        (scrc,132))

!  Set angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif

        scrc='('//Titlec(1:mchar)//') show'
        cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))

!  Reset angle
        if(ang.ne.0.) then
          scrc=' '
          write(scrc,'(F7.1,'' rotate'')') -ang
          cmdstr=cmdstr(1:lenstr(cmdstr,132))//' '// scrc(1:lenstr(scrc,132))
        endif
        call filler

!  Compute next character position
        cmdstr=' '
        write(cmdstr,'(f6.1,'' Xyprset0'')')ang
        call filler

   10 continue

!  Reset font
      call setfnt(nfsav)

      return
      end
!*****SYMBOL
      subroutine symbol(xp,yp,size,ltitle,ang,nchar)
      character*132 cmdstr,curfnt
      character*80 titlec,titleb
      character*1 bslash
      dimension ltitle(20)
      common/fntcom/curfnt,ifntsz,nfont
      common/plt1/cmdstr
      common/cnvcom/conver
      common/kkplot/szrat

      bslash=char(92)

!  Choose proper height
      mchar=iabs(nchar)
!  Set character size
      iht=size*conver/.6     !.6 factor is empirical

      if(iht.ne.ifntsz) then
        cmdstr=' '
        write(cmdstr,'(i3,'' Setf'')')iht
        call filler
        ifntsz=iht
      endif

      write(titlec,'(20a4)')ltitle

!  Check if titlec contains ( or ) or \.  These characters must be treated
!  specially by preceding them with a "\".  Do this to ( and ) even though
!  they might be balanced, i.e. () within a string, which can be treated
!  normally.

      titleb=titlec
      numc=0
      do m=1,mchar
        if(titleb(m:m).eq.'('.or.titleb(m:m).eq.')' .or.titleb(m:m).eq.        &
        bslash) then
          numc=numc+1
          titlec(numc:numc)=bslash
        endif
        numc=numc+1
        titlec(numc:numc)=titleb(m:m)
      enddo

      mchar=numc
      xpos=xp
      ypos=yp

      rsize=size

      arg=ang*4.*abs(atan(1.))/180.
      if(xpos.eq.999.) then
        cmdstr='/Xpos Xpres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(''/Xpos '',e10.4,'' def'')')xp*conver
      endif
      call filler

      if(ypos.eq.999.) then
        cmdstr='/Ypos Ypres def'
        njus=0
      else
        cmdstr=' '
        write(cmdstr,'(''/Ypos '',e10.4,'' def'')')yp*conver
      endif
      call filler

      if (nchar.eq.-1) then    !centered symbol
        high=rsize
        wide=high*szrat

        xpos=xpos+high/2.*sin(arg)-wide/2.*cos(arg)
        ypos=ypos-high/2.*cos(arg)-wide/2.*sin(arg)

        xarg=high/2.*sin(arg)-wide/2.*cos(arg)
        yarg=-high/2.*cos(arg)-wide/2.*sin(arg)

        cmdstr=' '
        write(cmdstr, '(''/Xpos Xpos'',f8.2,'' add def'')')xarg*conver

        call filler
        cmdstr=' '
        write(cmdstr, '(''/Ypos Ypos'',f8.2,'' add def'')')yarg*conver

        call filler
      endif

!  Move pen to proper coordinates
      cmdstr='xydef mover'
      call filler

!  Set angle
      if(ang.ne.0.) then
        cmdstr=' '
        write(cmdstr,'(F7.1,'' rotate'')') ang
        call filler
      endif

      cmdstr='('//Titlec(1:mchar)//') show'
      call filler

!  Reset angle
      if(ang.ne.0.) then
        cmdstr=' '
        write(cmdstr,'(F7.1,'' rotate'')') -ang
        call filler
      endif

      return
      end
