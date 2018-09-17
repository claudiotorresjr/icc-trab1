# GRR20176143 Cláudio Torres Júnior
# GRR20171607 Gabriela Stein
# -----------------------------------------------------------------------------

    CC     = gcc -std=c11 -g
    CFLAGS = -Wall
    LFLAGS = -lm 

      PROG = cgSolver
      OBJS = utils.o \
             sistemarandom.o \
             gradienteconjugado.o \
             $(PROG).o

.PHONY: doc purge clean all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS) $(INCLUDES) 


%.o: %.c %.h 
	$(CC) -c $(CFLAGS) -o $@ $<

all: $(PROG) doc

doc: Doc
	doxygen $<

Doc:
	doxygen -g

clean:
	rm -rf *~ *.bak 

purge: clean
	rm -rf html
	rm -f *.o $(PROG)
