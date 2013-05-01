CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config --cppflags`
LDFLAGS=`"${R_HOME}/bin/R" CMD config --ldflags`

all: phymlexe phytimeexe raxml-SSE3-PTHREADSexe mrbayesexe fasttreeexe lapexe argtablelib clustaloexe muscleexe

phymlexe:
	(cd phyml; ./configure; make)

phytimeexe:
	(cd phytime; ./configure --enable-phytime; make)

raxml-SSE3-PTHREADSexe:
	(cd raxml; make -f Makefile.SSE3.PTHREADS.gcc)

mrbayesexe:
	(cd mrbayes; autoconf --force; ./configure --with-beagle=no --enable-sse=yes; make)

fasttreeexe:
	(cd fasttree; $(CC) -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm)

lapexe:
	(cd lap; $(CC) $(CFLAGS) -o myLap myLap.c)

argtablelib:
	(cd argtable; ./configure --prefix=$(PWD)/build; make; make install)

clustaloexe:
	(cd clustalo; CPPFLAGS=-I$(PWD)/build/include LDFLAGS=-L$(PWD)/build/lib ./configure; make)

muscleexe:
	(cd muscle; make)

clean:
	(cd phyml; make clean)
	(cd phytime; make clean)
	(cd raxml; rm *.o raxml-SSE3-HPC)
	(cd mrbayes; make clean)
	(cd fasttree; rm *.o)
	(cd clustalo; make clean)
	(cd muscle; rm *.o)

