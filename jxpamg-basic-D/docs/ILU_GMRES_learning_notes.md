# ILU-GMRES 模块学习笔记

本文面向 `solver_id == 25` 的 ILU-GMRES 路径，目的是弄清当前 ILU 预条件器如何接入 GMRES、关键参数如何影响算法，并为后续“基于局部残差校正的组合型 ILU 预条件方法”设计做准备。

## 1. 入口与调用链

主入口在 `solver_ILU.c` 的 `case 25: /* ILU-GMRES */`。

核心调用顺序如下：

1. 创建 ILU 预条件器：
   - `JX_ILUCreate(&ilu_precond)`
   - 一组 `JX_ILUSet...` 设置 ILU 参数
2. 创建 GMRES 求解器：
   - `JX_GMRESCreate(comm, &solver)`
   - 设置 `k_dim`、`max_iter`、`tol`、`print_level`
3. 把 ILU 注册为 GMRES 预条件器：
   - `JX_GMRESSetPrecond(solver, JX_ILUSolve, JX_ILUSetup, ilu_precond)`
4. 显式执行 ILU setup：
   - `JX_ILUSetup(ilu_precond, hp_matrix, par_rhs, par_sol)`
5. GMRES setup 与 solve：
   - `JX_GMRESSetup(solver, hp_matrix, par_rhs, par_sol)`
   - `JX_GMRESSolve(solver, hp_matrix, hp_matrix, par_rhs, par_sol)`

需要注意：`case 25` 里先手动调用了一次 `JX_ILUSetup`，然后又把 ILU setup 函数注册进 GMRES。当前 `jx_GMRESSetup` 本身主要分配 Krylov 工作空间，并没有从这段已读代码中看到它主动调用预条件器 setup；因此手动 setup 是实际生效的关键步骤。

## 2. 当前 ILU-GMRES 是什么算法

外层是重启 GMRES(m)，其中 `m = k_dim`。当前实现采用右预条件形式：

```text
Krylov 基向量 v_i
z_i = M^{-1} v_i        // 调用 JX_ILUSolve
w_i = A z_i
对 w_i 做 Arnoldi 正交化
更新 Hessenberg 最小二乘问题
```

也就是说，GMRES 实际在近似求解：

```text
A M^{-1} y = b,  x = M^{-1} y
```

这里 `M^{-1}` 由 ILU 预条件器给出。由于是右预条件，打印和收敛判断仍围绕真实残差 `r = b - A x` 展开，比较适合后续加入“局部残差校正”后观察真实残差是否改善。

## 3. ILU 预条件器的两阶段

### 3.1 Setup 阶段

实现位置：

- `csrc/ilu/par_ilu.c`: 创建数据结构、参数 setter
- `csrc/ilu/par_ilu_setup.c`: 实际分解

`jx_ILUSetup` 做的事情：

1. 释放旧的 L/U/D/临时向量。
2. 创建工作向量 `Utemp`、`Ftemp`。
3. 保存矩阵、右端、解向量指针。
4. 若未给定 `perm`，调用 `jx_ILUGetLocalPerm` 生成局部重排序。
5. 根据 `ilu_type` 做分解：
   - `0`: Block Jacobi + ILU(0)
   - `1`: Block Jacobi + ILU(k)
   - `2`: Block Jacobi + ILUT
6. 统计 `operator_complexity = (n + nnz(L) + nnz(U)) / nnz(A)`。

当前简化版只支持 `ilu_type = 0, 1, 2`。源码里还保留了更复杂的类型说明，例如 Schur-GMRES、NSH、RAS、ddPQ、RAP-Modified-ILU，但在当前 `par_ilu_setup.c` 的主 switch 中不再启用。

### 3.2 Solve 阶段

实现位置：`csrc/ilu/par_ilu_solve.c`。

`jx_ILUSolve` 不是单纯覆盖式求解，而是做一次或多次预条件校正：

```text
r = f - A u
解 M e = r
u = u + e
```

在 GMRES 右预条件中，调用时输入通常是一个 Krylov 基向量，输出为预条件后的向量。`solver_ILU.c` 中把 ILU 的 `max_iter` 固定设为 `1`，所以每次预条件应用只做一次 ILU 校正。

`tri_solve = 1` 时使用直接前代/回代：

```text
L y = r
D/U 相关回代
u += e
```

`tri_solve = 0` 时用 Jacobi 迭代近似三角解，迭代次数由 `lower_jacobi_iters`、`upper_jacobi_iters` 控制。这通常更并行友好，但预条件质量可能弱于精确三角解。

## 4. 参数含义与影响

### 4.1 GMRES 参数

| 参数 | 默认值 | 位置 | 含义 | 主要影响 |
|---|---:|---|---|---|
| `solver_id` | `22`，运行 ILU-GMRES 需设为 `25` | `solver_ILU.c` | 选择求解器 | `25` 进入 ILU-GMRES |
| `tol` | `1e-6` | `solver_ILU.c` | GMRES 相对残差收敛阈值 | 越小精度越高，迭代可能更多 |
| `max_iter` | `500` | `solver_ILU.c` | GMRES 最大迭代次数 | 防止不收敛无限运行 |
| `k_dim` | `50` | `solver_ILU.c` | GMRES 重启长度 m | 越大 Krylov 子空间越丰富，但内存和正交化成本更高 |
| `is_check_restarted` | `1` | `solver_ILU.c` | 重启/假收敛检查 | 为 1 时遇到假收敛会重算真实残差并继续检查 |
| `print_level` | `3` | `solver_ILU.c` | GMRES 打印等级 | 影响日志详细程度 |

当前命令行已解析：

- `-sid 25`
- `-mxit <val>`
- `-kdim <val>`
- `-ptlv <val>`

没有看到 `tol` 的命令行解析分支；如果要批量实验不同精度，需要补一个参数，例如 `-tol <val>`。

### 4.2 ILU 类型与填充参数

| 参数 | 默认值 | 含义 | 主要影响 |
|---|---:|---|---|
| `ilu_type` | `0` | ILU 分解类型 | 当前支持 `0/1/2` |
| `ilu_lfil` | `0` | ILU(k) 的填充层级 | 只对 `ilu_type=1` 有效；越大预条件更强但更耗内存 |
| `ilu_droptol` | `1e-2` | ILUT 丢弃阈值 | 只对 `ilu_type=2` 有效；越小保留更多小元素 |
| `ilu_max_row_nnz` | `1000` | ILUT 每行最多保留非零数 | 控制 ILUT 内存上限和稀疏度 |

类型解释：

- `ilu_type = 0`: ILU(0)，只在原矩阵已有非零模式内做不完全分解，不允许额外填充。便宜、稳定可控，但预条件可能偏弱。
- `ilu_type = 1`: ILU(k)，允许 k 层 fill-in。`ilu_lfil=0` 时会退化到 ILU(0)。
- `ilu_type = 2`: ILUT，通过阈值和每行非零上限筛选填充。更灵活，但对 `droptol` 和 `max_row_nnz` 敏感。

重要实现细节：ILUT 中每行先计算输入行的平均绝对值 `inorm`，再设置行尺度阈值：

```text
itolb  = droptol[0] * inorm
itolef = droptol[1] * inorm
```

然后在消元时丢弃小于阈值的新增项，并在 L/U 中最多保留绝对值较大的 `lfil` 个条目。当前 `solver_ILU.c` 将 `ilu_max_row_nnz` 传给 ILUT 的 `lfil` 形参。

### 4.3 局部重排序参数

| 参数 | 默认值 | 含义 | 主要影响 |
|---|---:|---|---|
| `ilu_reordering` | `1` | 局部矩阵重排序 | `0` 不重排序；`1` 使用局部 RCM |

当前 `jx_ILUGetLocalPerm` 支持：

- `0`: identity permutation
- `1`: local RCM
- 其他值：默认也走 local RCM

RCM 主要目标是降低带宽，可能减少三角求解中的数据跨度，也可能影响 ILU 稳定性和填充结构。它是局部重排序，不跨 MPI 分区重排。

### 4.4 三角求解参数

| 参数 | 默认值 | 含义 | 主要影响 |
|---|---:|---|---|
| `ilu_tri_solve` | `1` | 三角解方式 | `1` 直接前代/回代；`0` Jacobi 近似三角解 |
| `ilu_ljac_iters` | `5` | L 部分 Jacobi 次数 | 只在 `ilu_tri_solve=0` 时有效 |
| `ilu_ujac_iters` | `5` | U 部分 Jacobi 次数 | 只在 `ilu_tri_solve=0` 时有效 |

直接三角解预条件效果更强，但串行依赖更重。Jacobi 三角解并行性更好，适合 GPU 或线程并行方向，但会引入“预条件器近似求解误差”。后续若做局部残差校正，这里是一个天然切入点：用局部残差判断哪些行/块需要更多 Jacobi 修正或更强局部解。

### 4.5 ILU 作为内层校正器的参数

| 参数 | 默认值 | 当前设置 | 含义 |
|---|---:|---:|---|
| `pc_tol` | `0.0` | `JX_ILUSetTol(ilu_precond, pc_tol)` | ILU solve 自身残差阈值 |
| `ilu_sm_max_iter` | `1` | 未直接使用；实际写死 `JX_ILUSetMaxIter(..., 1)` | 每次 ILU solve 的最大校正次数 |

由于 `pc_tol=0` 且 `max_iter=1`，ILU 预条件应用固定做一次校正，不追求内层收敛。这符合预条件器的常见用法，也意味着局部残差校正可以设计为“每次预条件应用内部的增强步骤”，而不是把 ILU solve 改成昂贵的全局迭代。

### 4.6 当前定义但未生效或未暴露的参数

以下变量在 `solver_ILU.c` 中定义或设置，但当前路径里没有实际影响，或没有命令行入口：

- `ilu_iter_setup_type`
- `ilu_iter_setup_option`
- `ilu_iter_setup_max_iter`
- `ilu_iter_setup_tolerance`
- `ilu_schur_max_iter`
- `ilu_nsh_droptol`
- `ilu_sm_max_iter`

已补充 `-ilu_dtol` 命令行解析，用于改变 ILUT 的 `ilu_droptol`。旧的 `-dtol` help 文案已替换。

## 5. 当前实现与历史设计之间的差异

`jx_ILUWriteSolverParams` 里还列出了多种 ILU 类型：

- `0`: Block Jacobi with ILU(k)
- `1`: Block Jacobi with ILUT
- `10/11`: ILU-GMRES with ILU(k)/ILUT
- `20/21`: Newton-Schulz-Hotelling with ILU(k)/ILUT
- `30/31`: RAS with ILU(k)/ILUT
- `40/41`: ddPQ-ILU-GMRES with ILU(k)/ILUT
- `50`: RAP-Modified-ILU

但当前 `par_ilu_setup.c` 的实际 switch 只启用：

- `0`: BJ + ILU(0)
- `1`: BJ + ILU(k)
- `2`: BJ + ILUT

这说明代码经历过裁剪或简化。后续做组合型 ILU 时，建议优先基于当前真实可运行的 `0/1/2` 路径扩展，而不是直接复活旧的 Schur/RAS/ddPQ 路径。

## 6. 对后续“局部残差校正组合型 ILU”的启发

当前 ILU solve 已经是校正形式：

```text
r = b - A x
e = M^{-1} r
x = x + e
```

这和局部残差校正非常契合。可以考虑在 `jx_ILUSolve` 或 `jx_ILUSolveLU` 后增加局部增强逻辑：

1. 计算局部残差指标：
   - 行残差 `r_i`
   - 块残差范数
   - 局部相对残差 `|r_i| / (|b_i| + eps)`
2. 根据阈值选择“坏区域”：
   - 残差大的行
   - 残差下降慢的块
   - 界面行或强耦合行
3. 对坏区域应用附加校正：
   - 增加局部 Jacobi/Gauss-Seidel 步数
   - 局部小块 ILU/ILUT
   - 局部稠密小解
   - 与 AMG 或 block correction 组合
4. 把校正作为预条件器的一部分返回给 GMRES。

一个保守的第一版可以设计为：

```text
标准 ILU solve 得到 e0
r1 = r - A e0
找出局部残差大的行/块
在这些行/块上求局部校正 ec
返回 e = e0 + ec
```

这样不改变 GMRES 外层，只增强 `M^{-1}`。风险和调试成本较低。

## 7. 建议下一步工作

1. 先补齐实验参数入口：
   - `-ilu_type`
   - `-ilu_lfil`
   - `-ilu_dtol`
   - `-ilu_max_row_nnz`
   - `-ilu_reordering`
   - `-ilu_tri_solve`
   - `-ilu_ljac`
   - `-ilu_ujac`
   - `-tol`
2. 建立基线实验表：
   - ILU(0) vs ILU(k) vs ILUT
   - RCM on/off
   - 直接三角解 vs Jacobi 三角解
   - 记录 GMRES 迭代数、setup 时间、solve 时间、operator complexity、最终残差
3. 在 `jx_ILUSolve` 中加入可选残差观测开关：
   - 每次预条件应用后的局部残差统计
   - 最大残差行、分位数、局部块范数
4. 再实现局部残差校正原型：
   - 先从纯局部 Jacobi 加强版开始
   - 然后再考虑块 ILU/ILUT 或界面区域校正

## 8. 一个需要优先修正的小问题

原始 `solver_ILU.c` 的 help 里曾写了：

```text
-dtol <val> : Drop-tolerance for ILU(0) factorization
```

当时没有对应解析逻辑，而且描述也不准确：drop tolerance 实际主要用于 ILUT，而不是 ILU(0)。现在已改成：

```text
-ilu_dtol <val> : drop tolerance for ILUT
```

并真正赋值给 `ilu_droptol`。
