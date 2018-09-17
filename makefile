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

.PHONY: doc purge clean

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS) $(INCLUDES) 


%.o: %.c %.h 
	$(CC) -c $(CFLAGS) -o $@ $<

all: $(PROG) doc

doc: Doxyfile
	doxygen $<

Doxyfile:
	doxygen -g
	sed -e "s;OPTIMIZE_OUTPUT_FOR_C *= *.*;OPTIMIZE_OUTPUT_FOR_C = YES;g" $@ >$@.new
	sed -e "s;EXTRACT_ALL *= *.*;EXTRACT_ALL = YES;g" $@.new >$@
	sed -e "s;EXTRACT_PRIVATE *= *.*;EXTRACT_PRIVATE = YES;g" $@ >$@.new
	sed -e "s;EXTRACT_STATIC *= *.*;EXTRACT_STATIC = YES;g" $@.new >$@
	sed -e "s;EXTRACT_LOCAL_METHODS *= *.*;EXTRACT_LOCAL_METHODS = YES;g" $@ >$@.new
	sed -e "s;GENERATE_LATEX *= *.*;GENERATE_LATEX = NO;g" $@.new >$@
	sed -e "s;SOURCE_BROWSE *= *.*;SOURCE_BROWSE = YES;g" $@ >$@.new
	sed -e "s;VERBATIM_HEADER *= *.*;VERBATIM_HEADER = YES;g" $@.new >$@
	rm -f $@.new

clean:
	rm -rf *~ *.bak

purge: clean
	rm -rf Doxyfile html latex
	rm -f *.o $(PROG)