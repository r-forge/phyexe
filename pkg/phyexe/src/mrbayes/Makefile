srcdir    = .
RM        = rm
CC        = gcc
CFLAGS    = -O3 -ffast-math -Wall -DNDEBUG -DUSECONFIG_H
LDFLAGS   =  -lreadline -lm
INSTALL   = /usr/bin/install -c

SRC       =  bayes.c command.c mbmath.c mcmc.c model.c plot.c sump.c sumt.c tree.c utils.c mbbeagle.c best.c
OBJECTS   =  bayes.o command.o mbmath.o mcmc.o model.o plot.o sump.o sumt.o tree.o utils.o mbbeagle.o best.o
PROGS     =  mb

all: $(PROGS)

mb: $(OBJECTS)

clean:
	$(RM) -f *.o *~ mb

distclean: clean
	$(RM) -f Makefile config.log config.status config.cache config.h 

showdep:
	@$(CC) -MM $(SRC)
	
bayes.o: bayes.c mb.h globals.h bayes.h command.h mcmc.h model.h utils.h
command.o: command.c mb.h globals.h command.h bayes.h model.h mcmc.h \
  plot.h sump.h sumt.h tree.h utils.h
mbmath.o: mbmath.c mb.h globals.h mbmath.h bayes.h model.h utils.h
mcmc.o: mcmc.c mb.h globals.h bayes.h mcmc.h model.h command.h mbmath.h \
  sump.h sumt.h plot.h tree.h utils.h mbbeagle.h
model.o: model.c mb.h globals.h bayes.h model.h command.h mcmc.h mbmath.h \
  tree.h utils.h
plot.o: plot.c mb.h globals.h command.h bayes.h plot.h sump.h utils.h
  sump.o: sump.c mb.h globals.h command.h bayes.h sump.h mcmc.h utils.h
sumt.o: sumt.c mb.h globals.h command.h bayes.h mbmath.h sumt.h sump.h \
  mcmc.h model.h utils.h
tree.o: tree.c bayes.h mb.h mcmc.h mbmath.h model.h globals.h tree.h \
  utils.h
utils.o: utils.c mb.h globals.h utils.h
mbbeagle.o: mbbeagle.c mbbeagle.h utils.h
best.o: best.c best.h
