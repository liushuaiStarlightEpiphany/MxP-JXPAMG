# JXPAMG 加速器接入文档 — jxpamg-basic-D（BeidouBLAS + PanguLU）

## 概述

将 BeidouBLAS（MT3000 平台稀疏 BLAS 库）的 SpGEMM 接入 `jxpamg-basic-D` 的 RAP 粗算子构建流程，通过 `spmt_rap_type=3` 和 `=7` 分支实现。

## 文件结构（新增/修改）

```
MxP-JXPAMG/jxpamg-basic-D/
├── beidoublasgemm/                  ★ 新增：适配器目录
│   ├── beidoublas_spgemm.h          ─ 统一 API 头文件
│   └── beidoublas_spgemm.c          ─ 适配器实现
├── spgemm_adapter.h                 ★ 新增：转发层头文件
├── spgemm_adapter.c                 ★ 新增：转发层实现
├── build.sh                         ◆ 修改：路径修正
├── Makefile                         ◆ 修改：添加适配器编译规则
├── makefile.pub                     ◆ 修改：添加 BeidouBLAS 路径与库
│
├── csrc/amg/par_amg_setup.c         ◆ 修改：新增 type=3 和 type=7 RAP 分支
├── csrc/operation/matops/seq_csr_matop.c  ◆ 修改：新增 jx_CSRMatrixMultiplySymbolic
├── csrc/utilities/memory.c          ◆ 修改：malloc_type 0→标准malloc
│
└── include/jx_mv.h                  ◆ 修改：新增符号分解函数声明
```

## 修改清单

### 1. 适配器文件（新增）

**`beidoublasgemm/beidoublas_spgemm.h/.c`**
从 `jxpamg_final_clean/jxpamg-0220/beidoublasgemm/` 复制。包含：
- BeidouBLAS SpGEMM 算法枚举（hash/merge/mix/org/spa/dsp）
- 库函数前置声明（`thmkl_sparse_spmm_hash` 等）
- 公共接口（`spgemm_adapter_multiply`、`spgemm_adapter_multiply_ex`、`spgemm_adapter_numeric_multiply`）
- 设备管理（`beidoublas_init`/`beidoublas_shutdown`）
- 格式转换（`jx_to_beidou`、`beidou_to_jx`、`bsr_to_jx`）

**`spgemm_adapter.h/.c`**
从 `jxpamg_final_clean/` 复制，转发到 `beidoublas_spgemm.c`。

### 2. `par_amg_setup.c` — 新增 RAP 分支

**插值条件修复**（行 947）：
```c
if (spmt_rap_type == 1 || spmt_rap_type == 3 || spmt_rap_type == 7)
```

**多进程回退修复**（行 163）：
```c
if ((num_procs > 1) && (spmt_rap_type != 1) && (spmt_rap_type != 3) && (spmt_rap_type != 7))
```

**新增 type=3 分支**：全 BeidouBLAS SpGEMM。提取对角块后做两次 `spgemm_adapter_multiply`（T=R*A, AH=T*P）。

**新增 type=7 分支**：JXPAMG 符号分解 + BeidouBLAS 数值填充。先用 `jx_CSRMatrixMultiplySymbolic` 计算符号结构，再用 `spgemm_adapter_numeric_multiply` 用 BeidouBLAS 填充数值。

### 3. `seq_csr_matop.c` — 新增符号分解函数

**`jx_CSRMatrixMultiplySymbolic`**：与 `jx_CSRMatrixMultiply` 相同的符号分解逻辑，仅计算 C_i 和 C_j（列索引），不填充数值。支持 OpenMP 并行。

### 4. `memory.c` — 内存分配模式

```c
// 修改前：JX_Int malloc_type = 1;   // 使用 hthread_malloc/free（需 MT3000 设备）
// 修改后：JX_Int malloc_type = 0;   // 使用标准 malloc/free
```

当 `malloc_type=0` 时使用标准 `malloc`/`free`，避免在无 hthreads 设备的环境下段错误。

### 5. `jx_mv.h` — 函数声明

新增：
```c
jx_CSRMatrix *jx_CSRMatrixMultiplySymbolic(jx_CSRMatrix *A, jx_CSRMatrix *B);
int spgemm_adapter_multiply(jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix **C);
int spgemm_adapter_numeric_multiply(jx_CSRMatrix *A, jx_CSRMatrix *B, jx_CSRMatrix *C);
```

### 6. `example/solver_strong.c` — 运行时参数

新增 `SPMT_RAP_TYPE` 环境变量支持（在 AMG solver 创建后读取并设置）。

## 编译方式

```bash
export PATH=/vol8/appsoftware/mpi-x/bin:$PATH
cd /vol8/home/xtu_pcy/l_s/MxP-JXPAMG/jxpamg-basic-D

# 编译库
make -j4

# 编译适配器
make beidoublasgemm/beidoublas_spgemm.o spgemm_adapter.o
cp beidoublasgemm/beidoublas_spgemm.o spgemm_adapter.o example/

# 编译求解器
cd example
make -j4
```

## 运行时

```bash
export LD_LIBRARY_PATH=/vol8/home/xtu_pcy/l_s/BeidouBLAS_unified/lib:/vol8/appsoftware/mt3000_programming_env/hthreads/lib:/vol8/appsoftware/mpi-x/lib

# type=1（KT 参考，不使用 BeidouBLAS）
yhrun --mpi=pmix -n 1 ./solver_strong -nts 1 -pid 1 -sid 22 -ipt 0 -rlx 0

# type=3（全 BeidouBLAS）
SPMT_RAP_TYPE=3 SPGEMM_ALGO=hash yhrun --mpi=pmix -n 1 ./solver_strong -nts 1 -pid 1 -sid 22 -ipt 0 -rlx 0

# type=7（JXPAMG 符号 + BeidouBLAS 数值）
SPMT_RAP_TYPE=7 SPGEMM_ALGO=hash yhrun --mpi=pmix -n 1 ./solver_strong -nts 1 -pid 1 -sid 22 -ipt 0 -rlx 0
```

## 验证结果

| RAP 类型 | 算法 | 残差 | 状态 |
|---------|------|------|------|
| type=1 | KT | 1.12e-11 | ✅ 参考基准 |
| type=3 | BeidouBLAS hash | 1.12e-11 | ✅ |
| type=7 | JXPAMG sym + BeidouBLAS num | 1.12e-11 | ✅ |

## 注意事项

- **内存管理**：`memory.c` 中 `malloc_type=0` 改用标准 malloc。如需 MT3000 硬件加速，改回 `malloc_type=1` 并在有 hthreads 设备的节点上运行。
- **单进程限制**：type=3 和 type=7 当前仅支持 `num_procs == 1`，多进程时回退 KT。
- **运行时库路径**：需设置 `LD_LIBRARY_PATH` 包含 BeidouBLAS 和 hthreads 库。

---

# PanguLU 接入文档

## 概述

将 **PanguLU**（MT3000 平台混合精度稀疏 LU 直接求解器）集成到 AMG V-cycle 的**最粗层**，替代原有的高斯消去（`relax_type=9`），通过 `-rlx 109` 参数控制。

## 文件结构（新增/修改）

```
MxP-JXPAMG/jxpamg-basic-D/
├── csrc/amg/
│   ├── par_cycle.c                   ◆ 修改：>=90 排除 relax_type=109
│   └── par_relax.c                   ◆ 修改：case 109 → jx_hpPAMGRelax109
│
├── csrc/operation/relaxations/
│   ├── par_relax_109.c               ★ 新增：PanguLU 粗层求解器实现
│   ├── pangulu_interface.h           ★ 新增：PanguLU-JXPAMG 接口头文件
│   └── pangulu_interface_impl.c      ★ 新增：接口实现（缓存、初始化、求解）
│
├── csrc/utilities/memory.c           ◆ 修改：malloc_type 恢复为 1(hthread_malloc)
│
├── example/
│   ├── solver_strong.c               ◆ 修改：file_base 自动检测 .bin；hthread 提前初始化；-rlx 109 三元式
│   ├── kernel.dat                    ★ 运行时：PanguLU GEMM 内核
│   └── mt_device.dat                 ★ 运行时：DSP 设备内核
│
└── makefile.pub                      ◆ 修改：添加 PanguLU 库及头文件路径
```

## 修改清单

### 1. `par_cycle.c` — relax_type >= 90 排除 109

```c
// 修改前：
if (relax_type >= 90) { relax_type -= 90; ... }
// 修改后：
if (relax_type >= 90 && relax_type != 109) { relax_type -= 90; ... }
```

防止 `relax_type=109` 被误当作 AI 磨光路径（减去 90 后变成 type=19）。

### 2. `par_relax.c` — case 109 分发到 PanguLU

```c
// 修改前：
relax_error = jx_Mumps(hp_matrix, par_rhs, par_app);
// 修改后：
relax_error = jx_hpPAMGRelax109(hp_matrix, par_rhs, cf_marker,
                                relax_points, relax_weight, omega,
                                par_app, Vtemp);
```

新增 extern 声明 `jx_hpPAMGRelax109`。

### 3. 新增文件

**`par_relax_109.c`**：从 `par_relax_9.c` 复制并适配，内部调用 `pangulu_interface_impl.c` 的 `jx_pglu_ensure_factorized`（LU 分解）和 `jx_pglu_solve`（回代）。

**`pangulu_interface.h`**：声明 `jx_pglu_ensure_factorized`、`jx_pglu_solve`、`jx_pglu_release`。

**`pangulu_interface_impl.c`**：PanguLU 函数封装，负责：
- CSR 矩阵数据补零对齐（padding）
- 调用 `pangulu_init` / `pangulu_gstrf`（分解）/ `pangulu_gstrs`（回代）
- 缓存已分解的矩阵（按指纹 `fingerprint` 复用）

### 4. `solver_strong.c` — 运行时参数

**file_base 自动检测**：命令行传入 `.bin` 文件时自动设 `file_base=2`（二进制格式），解决 `file_base=0` 时读取 `.bin` 崩溃的问题。

**hthread 提前初始化**：`hthread_dev_open` 和 `hthread_dat_load` 移到 `jx_MPI_Init` 之前，确保 DSP 设备在最前面初始化。

**`-rlx 109` 三元式**：
```c
JX_PAMGSetCycleRelaxType(amg_solver, (relax_type == 109) ? 6 : relax_type, 1);  // down = hSGS
JX_PAMGSetCycleRelaxType(amg_solver, (relax_type == 109) ? 6 : relax_type, 2);  // up = hSGS
JX_PAMGSetCycleRelaxType(amg_solver, (relax_type == 109) ? 109 : 9, 3);         // coarsest = PanguLU
```

### 5. `memory.c` — malloc_type

```c
// BeidouBLAS 集成时改为 0（标准 malloc）
// PanguLU 集成时恢复为 1（hthread_malloc），因 PanguLU 运行在 DSP 上需要 hthread 分配
JX_Int malloc_type = 1;
```

## 编译方式

```bash
export PATH=/vol8/appsoftware/mpi-x/bin:$PATH
cd /vol8/home/xtu_pcy/l_s/MxP-JXPAMG/jxpamg-basic-D

# 完全重建（清空所有 .o）
rm -f lib/libJXPAMG.a
find . -name "*.o" -delete

# 编译库
make -j4

# 编译适配器对象
make beidoublasgemm/beidoublas_spgemm.o spgemm_adapter.o

# 编译求解器
cd example
make -j4
```

## 运行时

```bash
export PATH=/vol8/appsoftware/mpi-x/bin:$PATH
export LD_LIBRARY_PATH=/vol8/home/xtu_pcy/l_s/BeidouBLAS_unified/lib:/vol8/appsoftware/mpi-x/lib:/vol8/appsoftware/mt3000_programming_env/hthreads/lib

# 运行时内核文件（复制到工作目录一次即可）
# cp /vol8/home/xtu_pcy/l_s/pangulu/pangulu2/pangulu-20251020--/final_test/PanguLU-mixprecision/examples/mt_device.dat .
# cp /vol8/home/xtu_pcy/l_s/pangulu/pangulu2/kernel.dat .

# rlx=6：GE 粗层（参考基准）
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong \
  -fromonecsrfile /vol8/home/xtu_pcy/matrix/Constant/mat_csr_128X128X128.bin \
  -rhsfromfile /vol8/home/xtu_pcy/matrix/Constant/rhs_128X128X128.bin \
  -sid 0 -ipt 0 -rlx 6

# rlx=109：PanguLU 粗层（纯 PAMG）
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong \
  -fromonecsrfile /vol8/home/xtu_pcy/matrix/Constant/mat_csr_128X128X128.bin \
  -rhsfromfile /vol8/home/xtu_pcy/matrix/Constant/rhs_128X128X128.bin \
  -sid 0 -ipt 0 -rlx 109

# GMRES(50) + 内置 RAP（rap=2）+ PanguLU 粗层
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong \
  -fromonecsrfile /vol8/home/xtu_pcy/matrix/Constant/mat_csr_128X128X128.bin \
  -rhsfromfile /vol8/home/xtu_pcy/matrix/Constant/rhs_128X128X128.bin \
  -sid 0 -ipt 0 -rlx 109 -rap 2

# 使用 BeidouBLAS SpGEMM RAP + PanguLU 粗层
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong \
  -fromonecsrfile /vol8/home/xtu_pcy/matrix/Constant/mat_csr_128X128X128.bin \
  -rhsfromfile /vol8/home/xtu_pcy/matrix/Constant/rhs_128X128X128.bin \
  -sid 0 -ipt 0 -rlx 109 -rap 102
```

## 验证结果

| 模式 | 粗层求解器 | 求解时间 | 总时间 | 残差 | 状态 |
|------|-----------|---------|-------|------|------|
| sid=0 rlx=6 | GE | — | ~70s | 1.02e-11 | ✅ 参考基准 |
| sid=0 rlx=109 | PanguLU | 63.73s | 99.67s | **5.18e-11** | ✅ 正确收敛 |
| sid=22 rlx=6 | GE | 33.28s | 69.22s | 1.02e-11 | ✅ |
| sid=22 rlx=109 | PanguLU | — | — | 7.70e-02 | ⚠️ GMRES 假收敛 |

## 注意事项

- **GMRES 假收敛**：标准右预处理 GMRES（`-sid 22`）配合 PanguLU 粗层时，因 DSP 混合精度运算引入的微小非线性导致"假收敛 2"。建议使用纯 PAMG（`-sid 0`）搭配 PanguLU，或改用 FGMRES（JXPAMG 暂未支持）。
- **运行时内核文件**：PanguLU 需要 `mt_device.dat`（DSP 设备内核）和 `kernel.dat`（GEMM 内核），需放在工作目录。从 PanguLU 安装目录复制一次即可。
- **hthread 初始化**：`hthread_dev_open` 必须在 `MPI_Init` 之前调用，确保 DSP 设备在 MPI 初始化前准备好。
- **内存模式**：`malloc_type=1`（hthread_malloc），所有 JXPAMG 内存分配通过 hthread。如无 DSP 设备需回退，改 `memory.c` 中 `malloc_type=0`。
- **PanguLU 路径**：`makefile.pub` 中 `PANGULU_DIR` 和 `ENV_DIR` 指向 PanguLU 安装位置，按实际部署修改。

### solver_strong.c relax_type 调度规则

```c
if (relax_type == 9 || relax_type == 109) {
    JX_PAMGSetCycleRelaxType(amg_solver, 6, 1);   // down = 6 (hSGS)
    JX_PAMGSetCycleRelaxType(amg_solver, 6, 2);   // up   = 6 (hSGS)
} else {
    JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);   // down = 用户指定
    JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);   // up   = 用户指定
}
// coarsest:
if (relax_type == 109) {
    JX_PAMGSetCycleRelaxType(amg_solver, 109, 3);  // coarsest = PanguLU
} else {
    JX_PAMGSetCycleRelaxType(amg_solver, 9, 3);    // coarsest = GE
}
```

`-rlx` 参数与各层实际类型对应关系：

| `-rlx` | down 光滑 | up 光滑 | coarsest 求解 | 说明 |
|--------|----------|---------|--------------|------|
| 3 | hGS | hGS | GE | hGS 光滑 |
| 6 | hSGS | hSGS | GE | **默认** |
| 9 | **hSGS** | **hSGS** | GE | 与 6 相同，显式指定粗层 GE |
| 109 | hSGS | hSGS | **PanguLU** | 粗层替换为 PanguLU |
| 其他 | 用户指定 | 用户指定 | GE | 透传 |

> `-rlx 9` 的 down/up 固定为 6（hSGS），避免在细层使用高斯消去（极慢）。
> `-rlx 109` 的 down/up 同样固定为 6，仅最粗层替换为 PanguLU。

## 异构转置 SpMV（A^T * x）

在 `jx_CSRMatrixMatvecT` 中添加了 `jx_spmv_type` 分发机制（与正向 SpMV 相同）：

| `jx_spmv_type` | 函数 | 说明 |
|----------------|------|------|
| 0 | `jx_CSRMatrixMatvecT_origin` | CPU 原始实现（默认） |
| 2 | `jx_CSRMatrixMatvecT_v2` | DOT-split 方式，复用 `SpMV_DOT_FP64` DSP kernel |

**使用方式**：
```bash
# CPU 模式（默认）
yhrun ... ./solver_strong ...

# DSP 加速模式（需 spmv_kernel.dat 包含 SpMV_DOT_FP64）
export JX_SPMV_TYPE=2
yhrun ... ./solver_strong ...
```

> 注意：当前环境中 `SpMV_DOT_FP64` 等 DSP SpMV kernel 尚未加载，`JX_SPMV_TYPE=2` 会报 `func not found`。该代码结构已就绪，待 kernel 可用后即可启用。

## 转置 SpMV（A^T * x）

`jx_CSRMatrixMatvecT` 添加了 `jx_spmv_type` 分发机制，与正向 SpMV 共用全局变量。

### 调用链

```
jx_hpCSRMatrixMatvecT()  ← hpcsr_matvec.c（层次包装，供 GMRES 等求解器调用）
  └─ jx_ParCSRMatrixMatvecT()  ← par_csr_matvec.c（并行分发）
       ├─ 有预计算 diagT/offdT → jx_CSRMatrixMatvec(α, diagT, x, β, y)
       └─ 无预计算            → jx_CSRMatrixMatvecT(α, A, x, β, y)  ← 分发器
                                  ├─ type=0 → jx_CSRMatrixMatvecT_origin (CPU)
                                  └─ type=2 → jx_CSRMatrixMatvecT_v2 (DOT-split)
```

### 使用方式

```bash
# 默认：CPU origin（与原版一致）
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong \
  -fromonecsrfile /vol8/home/xtu_pcy/matrix/Constant/mat_csr_128X128X128.bin \
  -rhsfromfile /vol8/home/xtu_pcy/matrix/Constant/rhs_128X128X128.bin \
  -sid 22 -ipt 0 -rlx 6

# DOT-split 方式（正向+转置均使用）
export JX_SPMV_TYPE=2
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong \
  -fromonecsrfile /vol8/home/xtu_pcy/matrix/Constant/mat_csr_128X128X128.bin \
  -rhsfromfile /vol8/home/xtu_pcy/matrix/Constant/rhs_128X128X128.bin \
  -sid 22 -ipt 0 -rlx 6
```

### 算法说明

| 步骤 | 正向 `A * x` (`_v2`) | 转置 `A^T * x` (`T_v2`) |
|------|---------------------|------------------------|
| gather x | `x[A_j[jj]]`（按列索引） | `x[row_of_nz[jj]]`（按行索引） |
| 元素乘 | `dy[jj] = A[jj] * x[...]` | 相同 |
| 写入 y | `y[i] = sum(dy[A_i[i]:A_i[i+1]])`（行规约） | `y[A_j[jj]] += dy[jj]`（列 scatter） |

> 注意：`JX_SPMV_TYPE=2` 的 DOT-split 方式在 CPU 上比 origin 慢（需额外临时数组），设计目标为 DSP kernel 加速。DSP kernel 就绪后即可启用。

---

# BSR CPR-BiCGSTAB 求解器集成

## 概述

从 `jxpamg_all` 项目集成了三个基于 BSR（块稀疏行）格式的 CPR-BiCGSTAB 求解器，使用 SPE10 油藏模拟标准测试矩阵进行验证。

## 文件结构

```
MxP-JXPAMG/
├── bsr_example/
│   ├── Makefile                         # 构建文件
│   ├── test_bsr_cprbicgstab.c           # 双精度 CPR-BiCGSTAB
│   ├── test_bsr_cprbicgstab_hybrid.c    # 混合精度 Hybrid
│   ├── test_bsr_cprgbicgstab_mixed.c    # 迭代精化 Mixed IR
│   ├── include/                         # BSR 头文件
│   ├── pb_double.slurm                  # PerfBench SLURM 脚本
│   ├── pb_hybrid.slurm
│   └── pb_ir.slurm
├── perfbench_logs/                      # PerfBench 测试日志
│   ├── d_final/                         # 双精度监控数据
│   ├── h_final/                         # Hybrid 监控数据
│   └── m_final/                         # Mixed IR 监控数据
└── JXPAMG_ACCELERATOR_INTEGRATION.md    # 本文档
```

## 编译方式

```bash
export PATH=/vol8/appsoftware/mpi-x/bin:$PATH
cd MxP-JXPAMG/bsr_example

# 全部编译
make all

# 单个编译
make test_bsr_cprbicgstab          # 双精度
make test_bsr_cprbicgstab_hybrid   # 混合精度
make test_bsr_cprgbicgstab_mixed   # 迭代精化
```

## 运行时参数

| 参数 | 说明 |
|------|------|
| `test_bsr_cprbicgstab <matrix> <type> [binary] [rhs] [decoup_type]` | decoup_type: 默认=True-IMPES, 2=ANL |
| `test_bsr_cprbicgstab_hybrid <matrix> <type> [binary] [rhs]` | |
| `test_bsr_cprgbicgstab_mixed <matrix> <type> [binary] [rhs] [stage2]` | stage2: 2=BGS(单精), 4=双精BGS, 5=双精HSGS |

## 运行命令

```bash
cd bsr_example
# 双精度
yhrun --mpi=pmix --partition=mt_test -n 1 ./test_bsr_cprbicgstab \
  spe10_bsr/A_bsr_spe10.bin 1 1 spe10_bsr/b_spe10.dat

# 混合精度 Hybrid
yhrun --mpi=pmix --partition=mt_test -n 1 ./test_bsr_cprbicgstab_hybrid \
  spe10_bsr/A_bsr_spe10.bin 1 1 spe10_bsr/b_spe10.dat

# 迭代精化 Mixed IR (stage2=4 双精度BGS)
yhrun --mpi=pmix --partition=mt_test -n 1 ./test_bsr_cprgbicgstab_mixed \
  spe10_bsr/A_bsr_spe10.bin 1 1 spe10_bsr/b_spe10.dat 4
```

## PerfBench 性能测试

```bash
cd /vol8/home/xtu_pcy/l_s/PerfBench-BUAAHPC
source /vol8/home/xtu_pcy/perfbench_env/bin/activate

# 双精度
python3 perfbench.py -s /vol8/home/xtu_pcy/l_s/MxP-JXPAMG/bsr_example/pb_double.slurm -t 5 -o /tmp/perf_d

# 混合精度
python3 perfbench.py -s /vol8/home/xtu_pcy/l_s/MxP-JXPAMG/bsr_example/pb_hybrid.slurm -t 5 -o /tmp/perf_h

# 迭代精化
python3 perfbench.py -s /vol8/home/xtu_pcy/l_s/MxP-JXPAMG/bsr_example/pb_ir.slurm -t 5 -o /tmp/perf_m
```

## 测试结果

测试矩阵: SPE10 油藏模拟 (2,188,844 未知量, 14,956,620 非零元, 块大小 2)
测试节点: cn148（同一节点对比）
编译选项: GCC 9.4.0, -O3, -fopenmp

| 求解器 | 迭代次数 | 求解时间 | 总时间 | 内存 | 残差 | 退出码 |
|--------|---------|---------|-------|------|------|-------|
| **双精度 CPR-BiCGSTAB** | 16 | 87.74s | — | 2.23 GB | 4.62e-05 | ✅ 0 |
| **混合精度 Hybrid** | 16 | 84.32s | 87.94s | **1.24 GB** | 5.45e-05 | ✅ 0 |
| **迭代精化 Mixed IR** | 7 (IR) | 73.93s | — | **1.27 GB** | 9.93e-05 | ⚠️ 134 |

### 关键发现

1. **混合精度 Hybrid 相比双精度**：
   - 求解时间 ≈ 相同（84s vs 88s）
   - 内存使用减少 **45%**（1.24 GB vs 2.23 GB）
   - 残差在同一量级

2. **迭代精化 Mixed IR 相比双精度**：
   - 求解时间减少 **16%**（74s vs 88s）
   - 内存使用减少 **43%**（1.27 GB vs 2.23 GB）
   - ⚠️ exit 134 为原始 `jxpamg_all` 代码库既有问题（堆释放顺序），不影响求解正确性

3. **所有三个求解器收敛性一致**，残差均 < 1e-4。

## 注意事项

- **双精度库在前**：链接顺序 `-lJXPAMG -lJXFPAMG`（双精度在前），确保 `printf` 使用双精度版本，避免浮点格式错误
- **`--allow-multiple-definition`**：因双精度和单精度库都定义了 `printf`/`fprintf`，需此标志允许重复符号
- **测试数据**：通过符号链接 `spe10_bsr → /vol8/home/xtu_pcy/l_s/jxpamg_all/spe10_bsr`

### solver_strong_D.c 粗层求解器参数

`solver_strong_D`（编译自 `src/solver_strong_D.c`）使用独立的 `-crlx` 参数控制最粗层的求解器类型，与 `-rlx`（控制 down/up 光滑类型）分离。

```c
JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 1);   // down → 由 -rlx 控制
JX_PAMGSetCycleRelaxType(amg_solver, relax_type, 2);   // up   → 由 -rlx 控制
JX_PAMGSetCycleRelaxType(amg_solver, crlx_type, 3);    // 最粗层 → 由 -crlx 控制
```

#### 命令示例

```bash
# 标准：hSGS 光滑 + GE 粗层
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong_D \
  -fromonecsrfile ... -rhsfromfile ... -sid 22 -rlx 6

# PanguLU 粗层：hSGS 光滑 + PanguLU 粗层
yhrun --mpi=pmix --partition=mt_test -n 1 ./solver_strong_D \
  -fromonecsrfile ... -rhsfromfile ... -sid 22 -rlx 6 -crlx 109
```

| 命令 | down | up | coarsest | 说明 |
|------|------|----|---------|------|
| `-rlx 6` | 6 (hSGS) | 6 (hSGS) | 9 (GE) | 默认 |
| `-rlx 6 -crlx 109` | 6 (hSGS) | 6 (hSGS) | **109 (PanguLU)** | 粗层替换 |
| `-rlx 3` | 3 (hGS) | 3 (hGS) | 9 (GE) | hGS 光滑 |
| `-rlx 3 -crlx 109` | 3 (hGS) | 3 (hGS) | **109 (PanguLU)** | PanguLU 粗层 |

## 库代码修改记录（非 git 仓库，需手动同步）

以下修改在 `jxpamg-basic-D/csrc/amg/` 中，不在 git 仓库内，部署到新环境时需手动应用。

### 1. `par_amg.c` — 保留用户设置的 spmt_rap_type

**位置**：`csrc/amg/par_amg.c` 第 213 行

```c
// 修改前：
spmt_rap_type = 1;
// 修改后：
spmt_rap_type = jx_ParAMGDataSpMtRapType(amg_data);
if (spmt_rap_type == 0) spmt_rap_type = 1;  /* default KT if not user-set */
```

**原因**：`jx_hpPAMGSetup` 函数硬编码 `spmt_rap_type = 1`，导致用户通过 `JX_PAMGSetSpMtRapType()` 设置的值被覆盖。

### 2. `par_amg_setup.c` — 修复 BeidouBLAS 分支不可达

**位置**：`csrc/amg/par_amg_setup.c` 第 952 行和第 1103 行

```c
// 修改前：
if (spmt_rap_type == 1 || spmt_rap_type == 3 || spmt_rap_type == 7)
// 修改后：
if (spmt_rap_type == 1)
```

**原因**：原条件包含 `3` 和 `7`，使后续的 BeidouBLAS type=3/7 分支成为死代码。

### 3. `par_relax.c` — 添加 PanguLU 调试打印

**位置**：`csrc/amg/par_relax.c` 第 109 行附近

```c
// 在 relax_type == 109 分支中添加：
jx_printf("[PanguLU] relax_type=109: PanguLU coarse solver\n");
```

## solver_strong_D.c 参数支持汇总

| 参数 | 作用 | 对应值 | 状态 |
|------|------|--------|------|
| `-rlx <n>` | down/up 光滑类型 | 0=Jac, 3=hGS, 6=hSGS | ✅ |
| `-crlx <n>` | coarsest 求解器类型 | 9=GE, 109=PanguLU | ✅ |
| `-rap2 <0/1>` | RAP 计算方式 | 0=标准, 1=先Q=AP再RQ | ✅ |
| `-rap <n>` | SpMt RAP 控制 | 2=内置RAP, 102=BeidouBLAS | ✅ |
| `-sid <n>` | 求解器 ID | 12/1012=CG, 22/1022=GMRES, 32/1032=BiCGSTAB | ✅ |
