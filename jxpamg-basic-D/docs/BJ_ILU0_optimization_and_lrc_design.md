# BJ+ILU(0) 当前实现与局部残差校正原型

## 1. 当前 BJ+ILU(0) 路径

入口：

```text
solver_ILU.c
  case 25
    JX_ILUSetType(ilu_precond, ilu_type)
    JX_ILUSetup(...)
    JX_GMRESSetPrecond(..., JX_ILUSolve, JX_ILUSetup, ilu_precond)
    JX_GMRESSolve(...)
```

当 `ilu_type = 0` 时：

```text
jx_ILUSetup
  jx_ILUGetLocalPerm              // 默认 local RCM
  jx_ILUSetupILU0
    jx_ILUSetupMILU0(..., modified=0)
```

Solve 阶段：

```text
jx_ILUSolve
  jx_ILUSolveLU                   // 直接三角解，默认
    r = f - A u
    L y = r
    U/D 回代得到 e
    u = u + e
```

也就是说，BJ+ILU(0) 是局部分块 Jacobi 预条件：每个 MPI rank 只对本地 diag block 做 ILU(0)，外部耦合只通过 `r = f - A u` 的残差计算进入，三角解本身没有跨 rank 耦合。

## 2. 待优化点

1. 参数不可实验

   原来 `solver_ILU.c` 里定义了 `ilu_type`、`ilu_lfil`、`ilu_reordering`、`ilu_tri_solve` 等变量，但命令行没有解析。这样很难比较 ILU(0)、ILU(k)、ILUT、RCM on/off、直接三角解/迭代三角解。

2. ILU(0) 只做局部 block Jacobi

   这是最主要的算法瓶颈。跨分区耦合只体现在 residual matvec 中，预条件校正阶段没有处理界面误差，因此强耦合问题或分区边界误差可能下降慢。

3. 三角解串行依赖较强

   默认 `tri_solve=1` 是直接前代/回代，预条件质量较好，但并行性有限。`tri_solve=0` 的 Jacobi 三角解更并行，但质量可能弱。后续可以让局部残差选择决定是否追加修正。

4. 没有局部残差诊断

   当前只记录整体相对残差。对组合型预条件方法来说，需要知道哪些行/块残差大、是否集中在界面、校正后局部残差是否下降。

5. ILU(0) pivot 保护较粗

   setup 中当对角元过小时直接设为 `1e-6` 再取倒数。这能避免除零，但可能引入数值扰动。后续可考虑行尺度相关的 pivot perturbation。

6. 重排序策略单一

   当前局部排序只有 identity 和 local RCM。RCM 对带宽有帮助，但未必总能改善 ILU 稳定性。后续可比较自然序、RCM、残差热点优先排序等。

## 3. 本次加入的局部残差校正原型

目标是不改变外层 GMRES，只增强 ILU 预条件器内部的 `M^{-1}` 应用：

```text
标准 ILU 校正:
  r0 = f - A u
  e0 = ILU(0)^{-1} r0
  u  = u + e0

局部残差选择:
  r1 = f - A u
  选择 |r1_i| >= theta * max_j |r1_j| 的局部行

局部增强校正:
  u_i = u_i + omega * r1_i / A_ii
```

这是一个对角 Jacobi 型的第一版增强。它的好处是：

- 默认关闭，不影响原始 BJ+ILU(0)。
- 插入点很小，在 `jx_ILUSolve` 内部完成。
- 使用真实残差选择热点行。
- 不需要改动 ILU(0) factor 数据结构。

局限也很明确：

- 目前只做点级对角校正，还不是块校正。
- 每一轮局部校正多一次 `A*u` matvec，会增加预条件器成本。
- 选择阈值基于全局最大残差，对离群值敏感；后续可以改为分位数或块范数。

## 4. 新增实验参数

GMRES/ILU 参数：

```text
-tol <val>
-ilu_type <0|1|2>
-ilu_lfil <val>
-ilu_dtol <val>
-ilu_max_row_nnz <val>
-ilu_reordering <0|1>
-ilu_tri_solve <0|1>
-ilu_ljac <val>
-ilu_ujac <val>
-ilu_ptlv <val>
```

局部残差校正参数：

```text
-ilu_lrc_type <0|1>        // 0 off, 1 diagonal residual correction
-ilu_lrc_iters <val>       // 局部校正 sweep 数
-ilu_lrc_threshold <val>   // 选择 |r_i| >= val * max |r_i| 的行
-ilu_lrc_omega <val>       // 对角校正松弛因子
```

建议第一组实验：

```text
原始基线:
  -sid 25 -ilu_type 0 -ilu_lrc_type 0

轻量局部校正:
  -sid 25 -ilu_type 0 -ilu_lrc_type 1 -ilu_lrc_iters 1 -ilu_lrc_threshold 0.5 -ilu_lrc_omega 0.8

更激进选择:
  -sid 25 -ilu_type 0 -ilu_lrc_type 1 -ilu_lrc_iters 1 -ilu_lrc_threshold 0.2 -ilu_lrc_omega 0.8

多 sweep:
  -sid 25 -ilu_type 0 -ilu_lrc_type 1 -ilu_lrc_iters 2 -ilu_lrc_threshold 0.5 -ilu_lrc_omega 0.6
```

## 5. 下一步更强版本

如果这个第一版能减少 GMRES 迭代数，后续可以沿三条线增强：

1. 从点校正升级到小块校正：
   - 对残差热点行扩展一层邻接行。
   - 抽取局部小块 `A_BB`。
   - 做小规模 ILU(0)/ILUT 或稠密解。

2. 从最大残差阈值升级到块范数选择：
   - 按 MPI 本地块、界面块、网格 patch 统计残差。
   - 对高残差块做组合型校正。

3. 面向界面误差：
   - 优先选择有 offd 耦合的行。
   - 将 BJ+ILU(0) 与一个轻量 Schur/interface correction 组合。
