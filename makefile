include ./makefile.pub

ALLPGS = solver_strong_D solver_strong_MxP
# ALLPGS = SPMV_test
# ALLPGS = solver_strong_MxP_v1

all: $(ALLPGS)

default: all

solver_src.o: ./src/solver_src.c
	$(CC) -o $@ ${CFLAGS} -c $<

solver_strong_F: ./src/solver_strong_F.o
	$(CC) -o $@ ${CFLAGS} $^ $(ALLLIBS)

solver_strong_D: ./src/solver_strong_D.o
	$(CC) -o $@ ${CFLAGS} $^ $(ALLLIBS)

solver_strong_MxP: ./src/solver_strong_MxP.o
	$(CC) -o $@ ${CFLAGS} $^ $(ALLLIBS)

solver_strong_MxP_v1: ./src/solver_strong_MxP_v1.o
	$(CC) -o $@ ${CFLAGS} $^ $(ALLLIBS)

solver_strong_MxP_v2: ./src/solver_strong_MxP_v2.o solver_src.o
	$(CC) -o $@ ${CFLAGS} $^ $(ALLLIBS)

SPMV_test: ./src/SPMV_test.o
	$(CC) -o $@ ${CFLAGS} $^ $(ALLLIBS)

clean:
	-rm -f *.o *~
	-rm -f $(ALLPGS)

