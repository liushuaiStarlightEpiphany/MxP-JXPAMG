#!/bin/bash

source /opt/intel/oneapi/setvars.sh --force
export MPIDIR=/opt/intel/oneapi/mpi/2021.6.0/bin
export MKL_DIR=/opt/intel/oneapi/mkl/2022.1.0

cd /home/spring/xingxin/jxpamg-basic-mxp/example

# # ================= 配置区域 =================
# PID="3 6 25"
# PID="31"
PID="6" 
N_TASKS="1" 
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


