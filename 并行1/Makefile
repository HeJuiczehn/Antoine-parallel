FC= mpif90
FFLAGS=-O2 -w   -malign-double -ffixed-line-length-none -fno-range-check
LIBS= -L../../lib -lantoine 
MASTER= ../../bin/antoine.out
MASTERSRC= aaa9.f
OBJECTS= acas.o achoi.o adime.o basis.o bfonc.o blok.o bmath0.o		\
bmath.o cmot.o cpiv.o ejlanc.o elanc.o emini.o endia.o epiv.o epivn.o	\
fbb.o hamil.o  hdiag.o imat.o lancj.o lanczos.o matbas.o	\
matgt.o matid.o matjz.o matop.o matpi.o matsp.o matt2.o mnnh1.o		\
mnnh2.o mpnh.o pivlec.o pnmat.o ppmat.o proj2.o prort.o prot2.o		\
qbase.o qgt.o qjz.o qprn.o qro.o qsp.o roop.o srand.o tconf.o tran.o tranpi.o	\
trbig.o trunq.o tsum.o xadd.o  xbibli.o xfile.o xjj.o xmosh.o xtest.o	\
xvecto.o ycopy.o  ylin.o ylist.o yorder.o yover.o yprec.o yspec.o	\
zcdm.o zcomp.o zcoul.o zdens.o zforma.o zidem.o zmadd.o zmodif.o  \
zmono.o zmsub.o

lib: $(OBJECTS)
	ar crv ../../lib/libantoine.a $(OBJECTS)
	$(FC)  -o $(MASTER) $(FFLAGS) $(MASTERSRC) $(LIBS)  
 
# clean out the directory
clean:
	rm -fr *~ *.o 
	rm -fr ../../lib/libantoine.a
	rm -fr ../../bin/antoine.out
