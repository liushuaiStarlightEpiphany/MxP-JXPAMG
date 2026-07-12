#!/bin/bash
#SBATCH --job-name=K-spmtv          # 作业名
#SBATCH -N 1                        # 节点数：1 到 4 个节点
#SBATCH -n 16                       # MPI任务数
#SBATCH --partition=mt_test         # 计算分区
#SBATCH --exclusive                 # 不共享节点

export UCX_TLS=^cma

PID="1"
N_TASKS="4"

# for pid in $PID
# do
#     for n in $N_TASKS
#     do
#         yhrun --mpi=pmix -n $n ./solver_strong_D -nts 1 -pid $pid -sid 22 -ipt 0
#     done
# done

for pid in $PID
do
    for n in $N_TASKS
    do
        yhrun --mpi=pmix -n $n ./solver_strong_MxP -nts 1 -pid $pid -sid 22 -ipt 0
    done
done