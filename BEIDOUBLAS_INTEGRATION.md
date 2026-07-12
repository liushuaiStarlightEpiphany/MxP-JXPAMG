# BeidouBLAS SpGEMM 接入文档 — jxpamg-basic-D

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
