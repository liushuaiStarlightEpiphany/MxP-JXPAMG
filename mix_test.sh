#!/bin/bash
#SBATCH --job-name=K-spmtv          # 作业名
#SBATCH -N 1                        # 节点数：1 到 4 个节点
#SBATCH -n 16                       # MPI任务数
#SBATCH --partition=mt_test         # 计算分区
#SBATCH --exclusive                 # 不共享节点

export UCX_TLS=^cma

# ================= 配置区域 =================
PID="30"
SID=42
IPT="0 60 6"    
N_TASKS="4" 
# ===========================================

# echo ""
# echo "========================================================"
# echo "Starting Benchmark: PID=$PID, SID=$SID, IPT=$IPT"
# echo "========================================================"

# # =========================================================
# # 1. 跑基准对照组 (Baseline: FP64)
# # =========================================================
echo ""
echo "-------------------------------------------------------"
echo " [Baseline] Solver: Strong (FP64)"
echo "-------------------------------------------------------"
    
for i in {1..1}
do
    echo "  -> Run #$i (FP64 Baseline)"
    for IPT in $IPT
    do
        yhrun --mpi=pmix -n $N_TASKS ./solver_strong -nts 1 -pid $PID -sid $SID -ipt $IPT
        sleep 2
    done
done

# # =========================================================
# # 2. 跑实验组 (MxP / FP32)
# # =========================================================
# echo ""
# echo "-------------------------------------------------------"
# echo " [Experiment] Solver: MxP (FP32 Mode)"   
# echo "-------------------------------------------------------"

# for i in {1..1}
# do
#     echo "  -> Run #$i (MxP, FP16_Lvl=0)" 
#     for PID in $PID
#     do
#         yhrun --mpi=pmix -n $N_TASKS ./solver_MxP -nts 1 -pid $PID -sid $SID -ipt $IPT
#         sleep 2
#     done
# done

# # # =========================================================
# # # 3. 跑实验组 (MxP / FP16)
# # # =========================================================
# echo ""
# echo "-------------------------------------------------------"
# echo " [Experiment] Solver: MxP (FP16 + FP32)"  
# echo "-------------------------------------------------------"

# for i in {1..1}
# do
#     echo "  -> Run #$i (MxP, FP16_Lvl=1)"      
#     for PID in $PID
#     do
#         yhrun --mpi=pmix -n $N_TASKS ./solver_MxP -nts 1 -pid $PID -sid $SID -ipt $IPT -fp16lvl 1
#         sleep 2
#     done
# done
    
echo ""
echo "========================================================"
echo "All Tests Completed."
echo "========================================================"


# # ================= 配置区域 =================
# SID=42
# IPT=60  
# N_TASKS_LIST="1 4 8 16"         # 弱扩展测试点
# # ===========================================

# echo ""
# echo "-------------------------------------------------------"
# echo " [Benchmark] Solver: Weak Scaling (FP64)"
# echo "-------------------------------------------------------"
    
# for i in {1..1}
# do
#     echo "  -> Run #$i"
    
#     for n in $N_TASKS_LIST
#     do
#         # ==========================================
#         # 核心逻辑：根据进程数 n 选择对应的 PID
#         # ==========================================
#         case $n in
#             1)
#                 CURRENT_PID="2"   # 64X64X64
#                 ;;
#             4)
#                 CURRENT_PID="27"  # 64X128X128
#                 ;;
#             8)
#                 CURRENT_PID="8"   # 128X128X128
#                 ;;
#             16)
#                 CURRENT_PID="28"  # 128X128X256
#                 ;;
#             *)
#                 echo "Error: Undefined PID for N_TASKS=$n"
#                 exit 1
#                 ;;
#         esac

#         echo "     [Running] Tasks: $n | Problem ID: $CURRENT_PID"

#         yhrun --mpi=pmix -n $n ./solver_strong -nts 1 -pid $CURRENT_PID -sid $SID -ipt $IPT
        
#         sleep 2
#     done
# done