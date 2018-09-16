# GRR20176143 Cláudio Torres Júnior
# GRR20171607 Gabriela Stein
# -----------------------------------------------------------------------------

    CC     = gcc -std=c11 -g
    CFLAGS = 
    LFLAGS = -lm -Wall

      PROG = cgSolver
      OBJS = utils.o \
             sistemarandom.o \
             gradienteconjugado.o \
             $(PROG).o

.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -f *.o core a.out
	@rm -f $(PROG)
