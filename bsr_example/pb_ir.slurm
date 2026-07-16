#!/bin/bash
#SBATCH -N 1 -n 1 --partition=mt_test
export PATH=/vol8/appsoftware/mpi-x/bin:$PATH
export LD_LIBRARY_PATH=/vol8/appsoftware/mpi-x/lib:$LD_LIBRARY_PATH
cd /vol8/home/xtu_pcy/l_s/MxP-JXPAMG/bsr_example
yhrun --mpi=pmix --partition=mt_test --mpi=pmix -n 1 ./test_bsr_cprgbicgstab_mixed spe10_bsr/A_bsr_spe10.bin 1 1 spe10_bsr/b_spe10.dat 4
