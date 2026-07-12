#!/bin/bash
#SBATCH --job-name=K-spmtv          # 作业名
#SBATCH -N 1                        # 节点数：1 到 4 个节点
#SBATCH -n 4                        # MPI任务数
#SBATCH --partition=mt_test         # 计算分区
#SBATCH --exclusive                 # 不共享节点

# # ================= 配置区域 =================
# PID="3 6 25"
# PID="31"
PID="3" 
N_TASKS="4" 
# # ===========================================
# # ===========================================

for pid in $PID 
do 
    for n in $N_TASKS
    do 
        mpirun -n $n ./solver_strong -nts 1 -pid $pid -sid 22 -ai_mt 0 -ait 0 -ipt 0 -rlx 3
    done
done

# for pid in $PID 
# do 
#     for n in $N_TASKS
#     do 
#         mpirun -n $n ./solver_strong -nts 1 -pid $pid -sid 22 -ai_mt 1 -ait 0 -ipt 60 -rlx 3
#     done
# done

# for pid in $PID 
# do 
#     for n in $N_TASKS
#     do 
#         mpirun -n $n ./solver_strong -nts 1 -pid $pid -sid 22 -ai_mt 1 -ait 2 -ipt 60 -rlx 3
#     done
# done

# for pid in $PID 
# do 
#     for n in $N_TASKS
#     do 
#         mpirun -n $n ./solver_strong -nts 1 -pid $pid -sid 22 -ai_mt 0 -ait 0 -ipt 6 -rlx 3
#     done
# done


