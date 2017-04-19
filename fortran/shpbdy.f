!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  SHPBDY           SHaPe BounDarY creates a map background
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:98-12-07
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   CREATES A MAP BACKGROUND ACCORDING TO THE PROJECTION
!   SPECIFIED BY THE PARMAP ARGUMENT USING ESRI FORMATTED SHAPE FILES
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 12 Aug 2008 (RRD) - initial version from MAPBDY
!
! USAGE:  CALL SHPBDY(KPROJ,LLGRID,PARMAP,AINTR,XTHICK,FNAME)
!   INPUT ARGUMENT LIST:
!     LLGRID - logical flag
!     PARMAP - real 9 element array defining the map projection
!     AINTR - real lat/lon line intervals (0 for none)
!     XTHICK - map background thickness (0.5 = half  1.0 = normal  2.0=twice)
!     FNAME - background file name
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     UNIT 15 - ascii map background vector file
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE SHPBDY(KPROJ,LLGRID,PARMAP,AINTR,XTHICK,AFILE)

  IMPLICIT NONE 

  INTEGER           :: KPROJ
  LOGICAL           :: LLGRID, FTEST, PTEST1, PTEST2, CLIP
  CHARACTER(80)     :: AFILE
  CHARACTER(20)     :: POINTS
  CHARACTER(6)      :: LATLON
  REAL              :: PARMAP(9)
  REAL, ALLOCATABLE :: PLAT(:),PLON(:),XARR(:),YARR(:)

  TYPE plot 
     SEQUENCE
     INTEGER(4)    :: DL   ! dash line = 0 (repeat/inch)
     REAL(4)       :: TH   ! line thickness = 0.01
     REAL(4)       :: RR   ! red = 0.0 
     REAL(4)       :: GG   ! green = 0.0 
     REAL(4)       :: BB   ! blue = 0.0
     CHARACTER(80) :: FN   ! shape file name in quotes 
  END TYPE
  TYPE (plot),     ALLOCATABLE :: shp(:) 

  CHARACTER(1),    ALLOCATABLE :: BUFF(:)
  INTEGER(4),      ALLOCATABLE :: PARTS(:)

  CHARACTER(80) :: fname = 'country.shp' 
  INTEGER(4)    :: handle,fcopen  ! file i/o handle 
  INTEGER(4)    :: file_length    ! length in bytes
  INTEGER(4)    :: kbyte = 0      ! file byte pointer
  INTEGER(4)    :: klen = 100     ! bytes to read      
  INTEGER(4)    :: mlen = 100     ! maximum buffer length
  INTEGER(4)    :: kret           ! decoder return code
  INTEGER(4)    :: rec_numb       ! record counter
  INTEGER(4)    :: num_parts      ! number of parts 
  INTEGER(4)    :: max_parts = 1  ! number of parts 
  INTEGER(4)    :: num_points     ! total number of points
  INTEGER(4)    :: num_pairs      ! x,y pairs in a part
  REAL(8)       :: Xpnt,Ypnt 

  REAL(4)       :: aintr,thick,xthick  
  REAL(4)       :: xx,yy,xp,yp,xc1,yc1,xc2,yc2,xu1,yu1,xu2,yu2  
  REAL(4)       :: clon,clat,xlon,xlat,xnew,ynew,xold,yold  

  INTEGER(4)    :: j,k,n,jbyte,nfiles 
  INTEGER(4)    :: lonx,latx,npts,kpts,mpts,marr,kline,intra  

  LOGICAL(4)    :: big_endian, diag
  COMMON /shpdesc/ big_endian, diag

  COMMON /MAPBND/ XU1,XU2,YU1,YU2
  COMMON /MAPCLP/ XC1,XC2,YC1,YC2

!-----------------------------------------------

  INTERFACE
  SUBROUTINE shphead(buff,file_length)
  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)  :: buff(*)
  INTEGER(4),   INTENT(OUT) :: file_length
  END SUBROUTINE shphead

  SUBROUTINE shpdat0(buff,rec_numb,rec_length)
  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)  :: buff(*)
  INTEGER(4),   INTENT(OUT) :: rec_numb
  INTEGER(4),   INTENT(OUT) :: rec_length
  END SUBROUTINE shpdat0

  SUBROUTINE shprec5(buff,num_parts,num_points)
  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)  :: buff(*)
  INTEGER,      INTENT(OUT) :: num_parts
  INTEGER,      INTENT(OUT) :: num_points
  END SUBROUTINE shprec5

  SUBROUTINE shploc5(buff,jbyte,num_parts,parts)
  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)    :: buff(*)
  INTEGER(4),   INTENT(INOUT) :: jbyte
  INTEGER(4),   INTENT(IN)    :: num_parts
  INTEGER(4),   INTENT(OUT)   :: parts(*)
  END SUBROUTINE shploc5
 
  SUBROUTINE shpvec5(buff,jbyte,Xpnt,Ypnt)
  IMPLICIT NONE
  CHARACTER(1), INTENT(IN)    :: buff(*)
  INTEGER,      INTENT(INOUT) :: jbyte
  REAL(8),      INTENT(OUT)   :: Xpnt,Ypnt
  END SUBROUTINE shpvec5
  END INTERFACE

!-----------------------------------------------

!==>limits check

      PTEST1(XP,YP)=(XP.GE.XU1.AND.XP.LE.XU2.AND.YP.GE.YU1.AND.YP.LE.YU2) 
      PTEST2(XP,YP)=(XP.GE.XU2.AND.XP.LE.XU1.AND.YP.GE.YU2.AND.YP.LE.YU1)

!     check for boundary clip region
      CLIP=.TRUE.

!     shape file diagnostics
      DIAG=.FALSE.

!=>initial  allocations

      MARR=5000 ! temporary array space for drawing
      ALLOCATE (XARR(MARR),YARR(MARR))
      MPTS=500  ! map background vector array
      ALLOCATE (PLAT(MPTS),PLON(MPTS))

!==>check range on line thickness multiplier

      THICK=AMAX1(0.01, AMIN1(100.0, XTHICK))

!==>set map paper units clip region

      CALL MAP2XY(XU1,YU1,XC1,YC1)
      CALL MAP2XY(XU2,YU2,XC2,YC2)

!------------------------------------------------------------------------------
! mark lat/lon line labels
!------------------------------------------------------------------------------

!     determine position of map center
      IF(LLGRID)THEN
         CALL GBL2LL(1,1,10.5*(XU1+XU2),0.5*(YU1+YU2),CLAT,CLON)
      ELSEIF(KPROJ.EQ.4)THEN 
         CALL CYL2LL(0.5*(XU1+XU2),0.5*(YU1+YU2),CLAT,CLON)
      ELSE
         CALL CXY2LL(PARMAP,0.5*(XU1+XU2),0.5*(YU1+YU2),CLAT,CLON)
      END IF

      INTRA=AINTR*10.0
!     print out lat/lon labels
      IF(AINTR.GT.0.0)THEN 

         CLON=NINT(CLON/AINTR)*AINTR
         CLAT=NINT(CLAT/AINTR)*AINTR

!        longitude labels
         IF(CLAT.GT.80.0.OR.CLAT.LT.-80.0)THEN
           XLAT=CLAT/2.0+AINTR/2.0
         ELSE
           XLAT=CLAT+AINTR/2.0
         END IF

         DO LONX=-(1800-INTRA),(1800-INTRA),INTRA
            XLON=LONX/10.0
            IF(LLGRID)THEN
               CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
            ELSEIF(KPROJ.EQ.4)THEN
               CALL CYL2XY(XLAT,XLON,XP,YP)
            ELSE
               CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
            END IF

            IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
               CALL MAP2XY(XP,YP,XX,YY)
               IF(AINTR.LT.1.0)THEN
                  WRITE(LATLON,'(F6.1)')XLON
                  CALL KEKSYMC(XX,YY,0.10,LATLON,0.0,6,1)
               ELSE
                  WRITE(LATLON,'(I4)')INT(XLON)
                  CALL KEKSYMC(XX,YY,0.10,LATLON,0.0,4,1)
               END IF
            END IF
         END DO

!        latitude labels
         XLON=CLON+AINTR/2.0

         DO LATX=-(900-INTRA),(900-INTRA),INTRA
            XLAT=LATX/10.0 
            IF(LLGRID)THEN
               CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
            ELSEIF(KPROJ.EQ.4)THEN
               CALL CYL2XY(XLAT,XLON,XP,YP)
            ELSE
               CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
            END IF

            IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
               CALL MAP2XY(XP,YP,XX,YY)
               IF(AINTR.LT.1.0)THEN
                  WRITE(LATLON,'(F5.1)')XLAT
                  CALL KEKSYMC(XX,(YY-0.05),0.10,LATLON,0.0,5,1)
               ELSE
                  WRITE(LATLON,'(I3)')INT(XLAT)
                  CALL KEKSYMC(XX,(YY-0.05),0.10,LATLON,0.0,3,1)
               END IF
            END IF
         END DO
      END IF

      IF(AINTR.NE.0.0)THEN 

!------------------------------------------------------------------------------
! draw lines of longitude  
!------------------------------------------------------------------------------

         DO LONX=-1800,(1800-INTRA),INTRA
            XLON=LONX/10.0
            N=0
            NPTS=0
            FTEST=.FALSE.

            DO LATX=-890,890,MAX(MIN(1,INTRA/10),1)
               XLAT=LATX/10.0
               N=N+1
  
               IF(LLGRID)THEN
                  CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
               ELSEIF(KPROJ.EQ.4)THEN
                  CALL CYL2XY(XLAT,XLON,XP,YP)
               ELSE
                  CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
               END IF

               IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
                  IF(.NOT.FTEST.AND.N.GT.1.AND.CLIP)THEN
!                    add previous position just outside domain
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL MAP2XY(XP,YP,XX,YY)
                  NPTS=MIN(NPTS+1,MARR)
                  XARR(NPTS)=XX
                  YARR(NPTS)=YY
                  FTEST=.TRUE.

               ELSEIF(FTEST)THEN
                  IF(CLIP)THEN
!                    add current off domain position to line
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
                  NPTS=0
                  FTEST=.FALSE.
               END IF

!              save previous value
               XOLD=XP
               YOLD=YP
            END DO
            IF(FTEST)CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
            IF(NPTS.EQ.MARR)WRITE(*,*)'Exceeding (XARR,YARR) array limits'
         END DO

!------------------------------------------------------------------------------
! draw lines of latitude  
!------------------------------------------------------------------------------

         DO LATX=-(900-INTRA),(900-INTRA),INTRA
            XLAT=LATX/10.0
            N=0
            NPTS=0
            FTEST=.FALSE.

            DO LONX=-1800,1800,MAX(MIN(1,INTRA/10),1)
               XLON=LONX/10.0
               N=N+1
  
               IF(LLGRID)THEN
                  CALL GBL2XY(1,1,XLAT,XLON,XP,YP)
               ELSEIF(KPROJ.EQ.4)THEN
                  CALL CYL2XY(XLAT,XLON,XP,YP)
               ELSE
                  CALL CLL2XY(PARMAP,XLAT,XLON,XP,YP)
               END IF

               IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
                  IF(.NOT.FTEST.AND.N.GT.1.AND.CLIP)THEN
!                    add previous position just outside domain
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL MAP2XY(XP,YP,XX,YY)

                  NPTS=MIN(NPTS+1,MARR)
                  XARR(NPTS)=XX
                  YARR(NPTS)=YY
                  FTEST=.TRUE.

               ELSEIF(FTEST)THEN
                  IF(CLIP)THEN
!                    add current off domain position to line
                     CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
                     CALL MAP2XY(XNEW,YNEW,XX,YY)
                     NPTS=MIN(NPTS+1,MARR)
                     XARR(NPTS)=XX
                     YARR(NPTS)=YY
                  END IF
                  CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
                  NPTS=0
                  FTEST=.FALSE.
               END IF

!              save previous value
               XOLD=XP
               YOLD=YP
            END DO
            IF(FTEST)CALL DSHCRV(XARR,YARR,NPTS,20,0.005*THICK)
            IF(NPTS.EQ.MARR)WRITE(*,*)'Exceeding (XARR,YARR) array limits'
         END DO
      END IF

!------------------------------------------------------------------------------
! open map background vector data file
!------------------------------------------------------------------------------

  INQUIRE(FILE=AFILE,EXIST=FTEST)
  IF(FTEST)THEN
     OPEN(15,FILE=AFILE,STATUS='OLD')
     KRET=0
     NFILES=0
     DO WHILE (kret.EQ.0)
        READ(15,*,IOSTAT=kret)
        NFILES=NFILES+1
     END DO
     REWIND(15)
     NFILES=NFILES-1
     ALLOCATE (shp(nfiles),STAT=kret)
     DO j=1,nfiles
        READ(15,*)shp(j)%fn,shp(j)%dl,shp(j)%th,shp(j)%rr,shp(j)%gg,shp(j)%bb
     END DO
     CLOSE(15)
  ELSE
     WRITE(*,*)'ERROR: shapefiles.txt file not found'
     WRITE(*,*)'Record format: ''file.shp'' dash thick red green blue'
     WRITE(*,*)'file.shp = /dir/name of input shapefile in quotes'
     WRITE(*,*)'dash = {0} for solid; {dashes}/in; <0 color fill'
     WRITE(*,*)'thick = line thickness in inches (default = 0.0)'
     WRITE(*,*)'Red Green Blue = RGB values (0.0 0.0 0.0 is black)'
     WRITE(*,*)'Sample file has been created in the working directory!'
     OPEN(15,FILE='shapefiles.txt',STATUS='NEW')
     WRITE(15,'(A)')'''arlmap.shp'' 0 0.005 0.4 0.6 0.8'
     CLOSE(15) 
     STOP
  END IF

  DO J=1,NFILES

! set color for this file (use defaults if thickness set to zero)
  IF(shp(j)%th.NE.0.0)THEN
     CALL setcolr(shp(j)%rr,shp(j)%gg,shp(j)%bb)
     shp(j)%th=0.01
  END IF

  HANDLE=FCOPEN(shp(j)%fn,'r')
  IF(.NOT.ALLOCATED(buff))  ALLOCATE (BUFF(mlen),STAT=kret)
  IF(.NOT.ALLOCATED(parts)) ALLOCATE (PARTS(max_parts+1),STAT=kret)

! read the file header record
  IF(diag) write(*,*)'Reading 1 start,length: ',kbyte,klen
  CALL FCPTPS(handle,kbyte,*910)
  CALL FCREAD(handle,buff,1,klen,*920)
  kbyte=kbyte+klen
  CALL shphead(buff,file_length)

!==>read vectors and draw lines

100 CONTINUE  
  IF(kbyte.GE.file_length)GOTO 900

! read the data record header
  klen=8
  IF(diag) write(*,*)'Reading 2 start,length: ',kbyte,klen
  CALL FCPTPS(handle,kbyte,*910)
  CALL FCREAD(handle,buff,1,klen,*920)
  kbyte=kbyte+klen
  CALL shpdat0(buff,rec_numb,klen)
  IF(klen.GT.mlen)THEN
     mlen=klen
     DEALLOCATE(buff)
     ALLOCATE(buff(mlen),STAT=kret)
     IF(diag) write(*,*)'Re-allocated buffer: ',mlen
  END IF
  IF(klen.EQ.0)GOTO 100 ! no contents, just rec numb

! read the data record contents
  IF(diag) write(*,*)'Reading 3 start,length: ',kbyte,klen
  CALL FCPTPS(handle,kbyte,*910)
  CALL FCREAD(handle,buff,1,klen,*820)
  820 CONTINUE ! last buffer read can go beyond EOF
  kbyte=kbyte+klen
  CALL shprec5(buff,num_parts,num_points)

  IF(num_points.GT.MPTS)THEN
     DEALLOCATE (PLAT,PLON)
     MPTS=num_points
     ALLOCATE (PLAT(MPTS),PLON(MPTS))
  END IF

! read the location points to each vector
  jbyte=45
  IF(num_parts.GT.max_parts)THEN
     DEALLOCATE(parts)
     ALLOCATE(parts(num_parts+1),STAT=kret)
     IF(diag) write(*,*)'Re-allocated num_parts: ',num_parts
     max_parts=num_parts
  END IF
  CALL shploc5(buff,jbyte,num_parts,parts)
  parts(num_parts+1)=num_points

! return the vectors for each part
  DO k=1,num_points 
     CALL shpvec5(buff,jbyte,Xpnt,Ypnt)
     IF(diag) write(*,*)'Vector point (byte,x,y): ',jbyte,Xpnt,Ypnt
     PLON(k)=Xpnt
     PLAT(k)=Ypnt 
  END DO
  IF(num_points.LE.1)GO TO 100

   DO KLINE=1,num_parts

      KPTS=0
      DO k=parts(kline)+1,parts(kline+1)
         KPTS=KPTS+1
         IF(KLINE.GT.1)THEN
            PLON(kpts)=PLON(k)
            PLAT(kpts)=PLAT(k)
         END IF 
      END DO

      NPTS=0
      FTEST=.FALSE.
      dline : DO N=1,KPTS

         IF(LLGRID)THEN
            CALL GBL2XY(1,1,PLAT(N),PLON(N),XP,YP)
         ELSEIF(KPROJ.EQ.4)THEN
            CALL CYL2XY(PLAT(N),PLON(N),XP,YP)
         ELSE
            CALL CLL2XY(PARMAP,PLAT(N),PLON(N),XP,YP)
         END IF

!        test for too long vectors due to coordinate change (prime meridian)
!        modification introduced 10/18/2007
         IF(N.GT.1)THEN
            IF(ABS(XP-XOLD).GT.0.5*(XU2-XU1).OR.    &
               ABS(YP-YOLD).GT.0.5*(YU2-YU1))THEN
               IF(NPTS.GT.1) THEN
                  IF(shp(j)%dl.EQ.0)THEN
                     CALL SLDCRV(XARR,YARR,NPTS,shp(j)%th*THICK)
                  ELSEIF(shp(j)%dl.GT.0)THEN
                     CALL DSHCRV(XARR,YARR,NPTS,shp(j)%dl,shp(j)%th*THICK)
                  ELSE
                     CALL FILRGNC(XARR,YARR,(NPTS-1))
                  END IF
               END IF
               NPTS=0
               FTEST=.FALSE.
               CYCLE dline 
            END IF
         ELSE
            XOLD=XP
            YOLD=YP
         END IF

         IF(PTEST1(XP,YP).OR.PTEST2(XP,YP))THEN
            IF(.NOT.FTEST.AND.N.GT.1.AND.CLIP)THEN
!              add previous position just outside domain
               CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
               CALL MAP2XY(XNEW,YNEW,XX,YY)
               NPTS=MIN(NPTS+1,MARR)
               XARR(NPTS)=XX
               YARR(NPTS)=YY
            END IF
            CALL MAP2XY(XP,YP,XX,YY)
            NPTS=MIN(NPTS+1,MARR)
            XARR(NPTS)=XX
            YARR(NPTS)=YY
            FTEST=.TRUE.

         ELSEIF(FTEST)THEN
            IF(CLIP)THEN
!              add current off domain position to line
               CALL BOUNDS(XOLD,YOLD,XP,YP,XNEW,YNEW)
               CALL MAP2XY(XNEW,YNEW,XX,YY)
               NPTS=MIN(NPTS+1,MARR)
               XARR(NPTS)=XX
               YARR(NPTS)=YY
            END IF
!           draw curve if anything in array
            IF(shp(j)%dl.EQ.0)THEN
               CALL SLDCRV(XARR,YARR,NPTS,shp(j)%th*THICK)
            ELSEIF(shp(j)%dl.GT.0)THEN
               CALL DSHCRV(XARR,YARR,NPTS,shp(j)%dl,shp(j)%th*THICK)
            ELSE
               CALL FILRGNC(XARR,YARR,(NPTS-1))
            END IF
!           start new line if out of window
            NPTS=0
            FTEST=.FALSE.
         END IF

!        save previous value
         XOLD=XP
         YOLD=YP
      END DO dline

!     complete curve if anything in array
      IF(FTEST)THEN
         IF(shp(j)%dl.EQ.0)THEN
            CALL SLDCRV(XARR,YARR,NPTS,shp(j)%th*THICK)
         ELSEIF(shp(j)%dl.GT.0)THEN
            CALL DSHCRV(XARR,YARR,NPTS,shp(j)%dl,shp(j)%th*THICK)
         ELSE
            CALL FILRGNC(XARR,YARR,(NPTS-1))
         END IF
      END IF 
      IF(NPTS.EQ.MARR)WRITE(*,*)'Exceeding (XARR,YARR) array limits'

    END DO 
    GOTO 100

910 KRET=1  
    write(*,*)'Shapefile decoder FCPTPS error'
    write(*,*) kbyte, file_length
    CALL FCCLOS(handle,*900)
    STOP 

920 KRET=1 
    write(*,*)'Shapefile decoder FCREAD error'
    write(*,*)'at byte: ', kbyte, file_length
    CALL FCCLOS(handle,*900)
    STOP 

900 KRET=0    !   normal end of file
    CALL FCCLOS(handle,*900)

    kbyte=0
    klen=100
! shape file loop
  END DO

! draw map boundary box in black before returning

  CALL SETCOLR(0.,0.,0.)
  CALL MAP2XY(XU1,YU1,XARR(1),YARR(1))
  CALL MAP2XY(XU1,YU2,XARR(2),YARR(2))
  CALL MAP2XY(XU2,YU2,XARR(3),YARR(3))
  CALL MAP2XY(XU2,YU1,XARR(4),YARR(4))
  CALL MAP2XY(XU1,YU1,XARR(5),YARR(5))
  CALL SLDCRV(XARR,YARR,5,0.015*THICK)

  DEALLOCATE (PLAT,PLON,XARR,YARR)

  END SUBROUTINE shpbdy
