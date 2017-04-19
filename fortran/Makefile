#-----------------------------------------------------------------
# Master Makefile for the ../source directory
# Last Revised: 10 Aug 2004
#               10 Jan 2005 - g95 fortran 
#               25 May 2005 - added massum 
#               13 Oct 2005 - lagrangian sampling
#               12 May 2006 - deallocation
#               12 Jan 2007 - cyl2xy
#               09 Aug 2007 - support for gfortran
#               30 May 2008 - added GEM routines
#               12 Aug 2008 - shapefile reader subroutines
#-----------------------------------------------------------------

SHELL = /bin/sh

SRC=./
LIB=./libhysplit.a

SUN5 = -O -free -I. 
AIX8 = -O3 -qarch=auto -qmaxmem=32768 -I. -qstrict
AIX6 = -O -bmaxdata:2000000000 -bmaxstack:256000000 -qarch=com -qmaxmem=8192 -I.
AIX5 = -O -qarch=com -qmaxmem=8192 -I.
DEC3 = -O -assume byterecl -I. 
SGI5 = -O -bytereclen -I. -freeform 
ABSF = -O1 -f free -I. 
#PGF9 = -O -Mfree -byteswapio -Mlfs -I.
#PGF9 = -O -g -Mfree -Mbounds -byteswapio -Mlfs -I.
PGF9 = -O3 -g -Mfree -Mbounds -byteswapio -Mlfs -I.
INTL = -O -FR -assume byterecl -convert big_endian -I.
GF95 = -O2 -g -fendian=big -ffree-form -ftrace=full -I.
GFOR = -O2 -fconvert=big-endian -frecord-marker=4 -fbounds-check -ffree-form -I.
#for Pleiades:
#FC = ifort
IFOR = -check bounds -free -O -w -I. -assume byterecl -mcmodel=medium #bounds checking will slow down the run
IFOR = -free -O -w -I. -assume byterecl -mcmodel=medium
#LDFLAGS = -shared-intel #needed for old stilt on Pleiades, not needed for multi (?)

HUM = -qsource -O2 -qmaxmem=-1 -q64 -C #AER-humboldt for xlf compilers; -C does bounds checking
#CFLAGS = $(GFOR)
CFLAGS = $(PGF9)
# LFLAGS = -vr -X64   (NCEP implementation)           
LFLAGS = -vr           
#FC = f90
#FC = gfortran
FC = pgf90
#FC = pgf95
CODE =          \
       funits.f \
       metval.f \
       constants_module.f \
       misc_definitions_module.f \
       mprintf.f \
       module_map_utils.f \
       gasdev.f \
       stbcon.f \
       mapbox.f \
       iercfg.f \
       adviec.f \
       advrkm.f \
       advrnt.f \
       adv2nt.f \
       adv3nt.f \
       adviso.f \
       advmet.f \
       advpnt.f \
       advsfc.f \
       advrng.f \
       basg2m.f \
       basm2g.f \
       chem02.f \
       cgrnll.f \
       cgszll_stilt.f \
       condsk.f \
       conini.f \
       conset.f \
       consum.f \
       conzro.f \
       cpolll_stilt.f \
       cyl2xy.f \
       datset.f \
       dalloc.f \
       decodi.f \
       decodr.f \
       depdry.f \
       depelm.f \
       deprad.f \
       depset.f \
       depsus.f \
       emsdat.f \
       emsgrd.f \
       emsini.f \
       emspnt.f \
       emsset.f \
       emsmat.f \
       emstmp.f \
       emrise.f \
       geplot.f \
       gbl2ll.f \
       gbl2xy.f \
       gbldxy.f \
       gbldll.f \
       gblset.f \
       geo_ll.f \
       grfbdy.f \
       lagems.f \
       lagout.f \
       lagsum.f \
       limmath.f \
       ll_geo.f \
       map_xe.f \
       map_xy.f \
       massum.f \
       metdiv.f \
       metgrd.f \
       metini.f \
       metinp.f \
       metlvl.f \
       metold.f \
       metpos.f \
       metset.f \
       metslp.f \
       metsub.f \
       metsum.f \
       metter.f \
       metwnd.f \
       mpstrt.f \
       p10adj.f \
       pakinp.f \
       pakout.f \
       pakset.f \
       pakrec.f \
       pakini.f \
       pakndx.f \
       pardsp.f \
       parinp.f \
       parpuf.f \
       parout.f \
       parvar.f \
       prfcom.f \
       prfecm.f \
       prfprs.f \
       prfsig.f \
       prfter.f \
       proj_3d.f \
       pufdel.f \
       pufdsp.f \
       pufmin.f \
       pufmrg.f \
       pufpar.f \
       pufrnd.f \
       pufsph.f \
       pufspv.f \
       pufsrt.f \
       runset.f \
       sfcinp.f \
       sobstr.f \
       stbanl.f \
       stbhor.f \
       stbsnd.f \
       stbvar.f \
       stbtke.f \
       stcm1p_stilt.f \
       stlmbr_stilt.f \
       sunang.f \
       sunflx.f \
       tm2day.f \
       tm2jul.f \
       tm2min.f \
       tminit.f \
       tmplus.f \
       trjdsk.f \
       trjset.f \
       x_prod.f \
       xe_xy.f \
       xy_map.f \
       con2xy.f \
       conavg.f \
       coninp.f \
       congrd.f \
       conmap.f \
       conndx.f \
       consmt.f \
       contur.f \
       contxt.f \
       map2xy.f \
       mapbdy.f \
       mapini.f \
       monset.f \
       module_defgrid.f \
       shpbdy.f \
       trjinp.f \
       trjlmt.f \
       trjmap.f \
       trjndx.f \
       cc2gll.f \
       cc2gxy.f \
       ccrvll.f \
       ccrvxy.f \
       cg2cll.f \
       cg2cxy.f \
       cg2wll.f \
       cg2wxy.f \
       cgszll.f \
       cgszxy.f \
       cgszxy_wps.f \
       cll2xy.f \
       cll2xy_wps.f \
       cnllxy.f \
       cnxyll.f \
       cpolll.f \
       cpolxy.f \
       cspanf.f \
       cw2gll.f \
       cw2gxy.f \
       cxy2ll.f \
       cxy2ll_wps.f \
       eqvlat.f \
       stcm1p.f \
       stcm2p.f \
       stlmbr.f \
       fcclos.f \
       fcopen.f \
       fcptps.f \
       fcgtps.f \
       fcread.f \
       w3code.f \
       iercon.f \
       ierno2.f \
       ierset.f \
       iersum.f \
       iervoc.f \
       grseqn.f \
       grscon.f \
       grsdst.f \
       prchem.f \
       chmphg.f \
       viscrt.f \
       gemcfg.f \
       gemkon.f \
       gemvar.f \
       gemchk.f \
       gemcol.f \
       gemdat.f \
       gemday.f \
       gemdep.f \
       gemdiv.f \
       gemeqn.f \
       gemgrd.f \
       gemini.f \
       geminp.f \
       gemout.f \
       gempar.f \
       gemset.f \
       gemstb.f \
       gemswf.f \
       gemsum.f \
       gasdev.f \
       ext_zsg.f \
       cgrell.f \
       output.f \
       ran3.f \
       adv3ntZML.f \
       adv3ntZLOC.f \
       advmetGRELL.f \
       adv3ntWIND.f \
       prfwrf.f \
       rsat.f \
       rtsafe.f \
       thes.f \
       sunclr.f \
       sunave.f \
       gemtst.f

all: $(LIB) hymodelc

$(LIB) : $(addsuffix .o, $(basename $(CODE)))
	$(AR) $(LFLAGS) $(LIB) $^
	ranlib $(LIB)

hymodelc : hymodelc.f90 $(LIB)
	$(FC) -o $@ $(LDFLAGS) $(CFLAGS) hymodelc.f90 -L./ -lhysplit 

.f.o:
	$(FC) -c $(CFLAGS) $<

clean :
	rm -f hymodelc
	rm -f *.o
	rm -f *.a
	rm -f *.mod
	rm -f *.lst

advmet.o: DEFMETO.INC
advmetGRELL.o: DEFMETO.INC
advpnt.o: funits.o metval.o module_defgrid.o DEFARG2.INC DEFARG3.INC \
	DEFCONC.INC DEFMETO.INC
advrng.o: funits.o module_defgrid.o
cc2gxy.o: module_defgrid.o
cgrell.o: funits.o
cgszxy_wps.o: module_defgrid.o module_map_utils.o
chem02.o: funits.o DEFCONC.INC
chmphg.o: DEFCONC.INC
cll2xy_wps.o: module_defgrid.o module_map_utils.o
condsk.o: funits.o module_defgrid.o DEFCONC.INC
conini.o: funits.o module_defgrid.o DEFCONC.INC DEFSPOT.INC
conndx.o: mapbox.o
conset.o: funits.o DEFCONC.INC
consum.o: funits.o DEFCONC.INC
conzro.o: DEFCONC.INC
cxy2ll_wps.o: module_defgrid.o module_map_utils.o
dalloc.o: metval.o
datset.o: funits.o
depdry.o: DEFCONC.INC
depelm.o: DEFCONC.INC
deprad.o: DEFCONC.INC
depset.o: funits.o DEFCONC.INC
depsus.o: funits.o metval.o module_defgrid.o DEFCONC.INC
emrise.o: funits.o
emsdat.o: funits.o
emsgrd.o: funits.o module_defgrid.o DEFCONC.INC DEFSPOT.INC
emsini.o: funits.o
emsmat.o: funits.o module_defgrid.o DEFSPRT.INC
emspnt.o: funits.o DEFCONC.INC DEFSPOT.INC
emsset.o: funits.o DEFCONC.INC
emstmp.o: funits.o DEFSPRT.INC
ext_zsg.o: funits.o
gbl2ll.o: module_defgrid.o
gbl2xy.o: module_defgrid.o
gbldll.o: module_defgrid.o
gbldxy.o: module_defgrid.o
gblset.o: module_defgrid.o
gemchk.o: funits.o gemcfg.o gemkon.o gemvar.o module_defgrid.o
gemcol.o: funits.o gemcfg.o gemkon.o gemvar.o
gemdat.o: funits.o gemkon.o gemvar.o
gemday.o: funits.o gemcfg.o gemkon.o gemvar.o
gemdep.o: gemcfg.o gemkon.o gemvar.o DEFCONC.INC
gemdiv.o: gemcfg.o gemkon.o gemvar.o
gemeqn.o: gemcfg.o gemkon.o gemvar.o
gemgrd.o: funits.o gemcfg.o gemkon.o gemvar.o
gemini.o: funits.o gemcfg.o gemkon.o gemvar.o
geminp.o: funits.o gemcfg.o gemkon.o gemvar.o module_defgrid.o
gemout.o: funits.o gemcfg.o gemkon.o gemvar.o
gempar.o: gemcfg.o gemkon.o gemvar.o module_defgrid.o
gemset.o: funits.o gemkon.o gemvar.o
gemstb.o: funits.o gemcfg.o gemkon.o gemvar.o
gemsum.o: gemcfg.o gemvar.o DEFCONC.INC
gemswf.o: funits.o gemcfg.o gemvar.o
gemtst.o: gemcfg.o gemvar.o DEFMETO.INC
grscon.o: funits.o metval.o module_defgrid.o DEFCONC.INC
grsdst.o: metval.o module_defgrid.o DEFCONC.INC
grseqn.o: funits.o
iercon.o: funits.o iercfg.o metval.o module_defgrid.o DEFCONC.INC
ierset.o: funits.o iercfg.o
iersum.o: funits.o iercfg.o DEFCONC.INC
iervoc.o: funits.o metval.o module_defgrid.o DEFSPOT.INC
lagems.o: funits.o DEFLAGS.INC
lagout.o: funits.o module_defgrid.o DEFLAGS.INC
lagsum.o: funits.o DEFCONC.INC DEFLAGS.INC
mapini.o: mapbox.o
massum.o: funits.o DEFCONC.INC
metgrd.o: funits.o module_defgrid.o
metini.o: funits.o module_defgrid.o module_map_utils.o
metinp.o: funits.o module_defgrid.o
metlvl.o: module_defgrid.o
metold.o: module_defgrid.o
metpos.o: funits.o module_defgrid.o
metset.o: funits.o module_defgrid.o
metsub.o: funits.o module_defgrid.o
metsum.o: funits.o DEFCONC.INC
metwnd.o: funits.o
module_defgrid.o: DEFGRID.INC
module_map_utils.o: constants_module.o misc_definitions_module.o mprintf.o
pakini.o: DEFPACK.INC
pakndx.o: DEFPACK.INC
pakrec.o: DEFPACK.INC
pakset.o: DEFPACK.INC
pardsp.o: DEFARG2.INC
parinp.o: funits.o module_defgrid.o
parout.o: funits.o module_defgrid.o
parpuf.o: funits.o
prchem.o: funits.o
prfcom.o: funits.o module_defgrid.o DEFARG4.INC
prfecm.o: funits.o
prfprs.o: funits.o
prfsig.o: funits.o
prfter.o: funits.o
prfwrf.o: funits.o
pufpar.o: funits.o
pufsph.o: funits.o
runset.o: funits.o module_defgrid.o DEFSPOT.INC
sfcinp.o: funits.o
stbanl.o: funits.o stbcon.o
stbsnd.o: stbcon.o
stbtke.o: stbcon.o
stbvar.o: stbcon.o
trjdsk.o: funits.o
trjlmt.o: mapbox.o
trjset.o: funits.o module_defgrid.o DEFSPOT.INC
hymodelc.o: funits.o module_defgrid.o DEFARG1.INC DEFARG2.INC DEFARG3.INC \
	DEFCHEM.INC DEFCONC.INC DEFLAGS.INC DEFMETO.INC DEFSPOT.INC \
	DEFSPRT.INC
