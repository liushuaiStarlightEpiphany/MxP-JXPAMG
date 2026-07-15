#========================================================================#
#  MxP-JXPAMG BSR Example Makefile                                       #
#========================================================================#

JXPAMGDIR   = /vol8/home/xtu_pcy/l_s/jxpamg_all/JXPAMG-bsr/jxpamg-0220
JXFPAMGDIR  = /vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/jxfpamg-0220
MPIDIR      = /vol8/appsoftware/mpi-x

JXPAMGINCLUDES = -I $(JXPAMGDIR)/include
JXFPAMGINCLUDES = -I $(JXFPAMGDIR)/include
JXPAMGLIB  = -L $(JXPAMGDIR)/lib -lJXPAMG
JXFPAMGLIB = -L $(JXFPAMGDIR)/lib -lJXFPAMG
MPILIB     = -L $(MPIDIR)/lib -lmpi

CFLAGS = -O3 -fopenmp -Wall \
         -DJX_USING_MPI=1 -DJX_USING_OPENMP=1 \
         -DJX_USING_BIG_INT=0 -DJX_USING_BIG_DOUBLE=0 \
         $(JXPAMGINCLUDES) $(JXFPAMGINCLUDES) -I $(MPIDIR)/include

LDLIBS = $(JXFPAMGLIB) $(JXPAMGLIB) $(MPILIB) -lm -lstdc++ -Wl,--allow-multiple-definition

ALL = test_bsr_cprbicgstab test_bsr_cprbicgstab_hybrid test_bsr_cprgbicgstab_mixed

all: $(ALL)

test_bsr_cprbicgstab: test_bsr_cprbicgstab.o
	mpicc -fopenmp -o $@ $^ $(LDLIBS)

test_bsr_cprbicgstab_hybrid: test_bsr_cprbicgstab_hybrid.o
	mpicc -fopenmp -o $@ $^ $(LDLIBS)

test_bsr_cprgbicgstab_mixed: test_bsr_cprgbicgstab_mixed.o
	mpicc -fopenmp -o $@ $^ $(LDLIBS)

.c.o:
	mpicc -fopenmp -o $@ -c $(CFLAGS) $<

data:
	ln -sf /vol8/home/xtu_pcy/l_s/jxpamg_all/spe10_bsr .

clean:
	rm -f *.o $(ALL)
