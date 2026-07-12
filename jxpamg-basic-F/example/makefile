include ../makefile.pub

ALLPGS =  solver_strong

all: $(ALLPGS)

default: all

solver_laplace: par_laplace.o solver_laplace.o
	$(CC) -o $@ $^ $(ALLLIBS)

solver_asetup: solver_asetup.o
	$(CC) -o $@ $^ $(ALLLIBS)

solver_strong: solver_strong.o
	$(CC) $(C_COMPILE_FLAGS) -o $@ $^ $(ALLLIBS)

solver_weak: solver_weak.o
	$(CC) -o $@ $^ $(ALLLIBS)

solver_strong_ij: solver_strong_ij.o
	$(CC) -o $@ $^ $(ALLLIBS)

solver_blockprec: solver_blockprec.o
	$(CC) -o $@ $^ $(ALLLIBS)

hx_driver: hx_driver.o
	$(CC) -o $@ $^ $(ALLLIBS)
clean:
	-rm -f *.o *~
	-rm -f $(ALLPGS)
