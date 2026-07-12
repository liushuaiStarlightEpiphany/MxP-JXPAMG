#!/bin/bash
#SBATCH --job-name=compare-all
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --partition=mt_test
#SBATCH --exclusive

cd /vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/example

OUTDIR=/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/example/log
N=1

echo "================================================================"
echo "  1/5  Baseline: ipt=60, ai_ipt=1 (original AHI)"
echo "================================================================"
yhrun --mpi=pmix -n $N ./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 1 -ipt 60 2>&1 | tee $OUTDIR/compare_baseline.log

echo ""
echo "================================================================"
echo "  2/5  ipt=61, ai_ipt=2 (classical+mix, th=0.8)"
echo "================================================================"
yhrun --mpi=pmix -n $N ./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 2 -ipt 61 2>&1 | tee $OUTDIR/compare_ipt61.log

echo ""
echo "================================================================"
echo "  3/5  ipt=62, ai_ipt=2 (classical+fix, th=0.8)"
echo "================================================================"
yhrun --mpi=pmix -n $N ./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 2 -ipt 62 2>&1 | tee $OUTDIR/compare_ipt62.log

echo ""
echo "Done. Results:"
echo "  $OUTDIR/compare_baseline.log"
echo "  $OUTDIR/compare_ipt61.log"
echo "  $OUTDIR/compare_ipt62.log"
