#========================================================================#
#  MxP-JXPAMG BSR Example Makefile                                       #
#========================================================================#

JXPAMGDIR   = /vol8/home/xtu_pcy/l_s/jxpamg_all/JXPAMG-bsr/jxpamg-0220
JXFPAMGDIR  = /vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/jxfpamg-0220
JXPAMG_ROOT = /vol8/home/xtu_pcy/l_s/jxpamg_all/JXPAMG-bsr
JXFPANG_ROOT = /vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr
MPIDIR      = /vol8/appsoftware/mpi-x

JXPAMGINCLUDES = -I $(JXPAMG_ROOT) -I $(JXPAMGDIR)/include
JXFPAMGINCLUDES = -I $(JXFPANG_ROOT) -I $(JXFPAMGDIR)/include
JXPAMGLIB  = -L $(JXPAMGDIR)/lib -lJXPAMG
JXFPAMGLIB = -L $(JXFPAMGDIR)/lib -lJXFPAMG
MPILIB     = -L $(MPIDIR)/lib -lmpi

BASE_CFLAGS = -O3 -fopenmp -Wall
BASE_CFLAGS += -Wno-unused-result -Wno-unused-variable
BASE_CFLAGS += -Wno-unused-but-set-variable -Wno-maybe-uninitialized
BASE_CFLAGS += -DJX_USING_MPI=1 -DJX_USING_OPENMP=1
BASE_CFLAGS += -DJX_USING_BIG_INT=0 -DJX_USING_BIG_DOUBLE=0
BASE_CFLAGS += $(JXPAMGINCLUDES) $(JXFPAMGINCLUDES) -I $(MPIDIR)/include

JXF_CFLAGS = -DJXF_USING_BIG_INT=0 -DJXF_USING_BIG_DOUBLE=0
JXF_CFLAGS += -DJXF_SINGLE=0 -DJXF_USING_OPENMP=1

LDLIBS = $(JXFPAMGLIB) $(JXPAMGLIB) $(MPILIB) -lm -lstdc++
LDLIBS += -Wl,--allow-multiple-definition

ALL = test_bsr_cprbicgstab test_bsr_cprbicgstab_hybrid test_bsr_cprgbicgstab_mixed

all: $(ALL)

test_bsr_cprbicgstab: test_bsr_cprbicgstab.o
	mpicc -fopenmp -o $@ $< $(LDLIBS)

test_bsr_cprbicgstab_hybrid: test_bsr_cprbicgstab_hybrid.o
	mpicc -fopenmp -o $@ $< $(LDLIBS)

test_bsr_cprgbicgstab_mixed: test_bsr_cprgbicgstab_mixed.o
	mpicc -fopenmp -o $@ $< $(LDLIBS)

.c.o:
	mpicc -fopenmp -o $@ -c $(BASE_CFLAGS) $(CFLAGS_EXTRA) $<

test_bsr_cprbicgstab_hybrid.o: CFLAGS_EXTRA = $(JXF_CFLAGS)
test_bsr_cprgbicgstab_mixed.o: CFLAGS_EXTRA = $(JXF_CFLAGS)

data:
	ln -sf /vol8/home/xtu_pcy/l_s/jxpamg_all/spe10_bsr .

clean:
	rm -f *.o $(ALL)
