content = """# AHI 自适应插值方法：基于理想插值后验误差的 EI 点选择

## 背景

现有 AHI（Adaptive Hybrid Interpolation）方法在 `par_interp_mix.c` 中实现，
通过 `jx_PAMGMeasureAI_ipt` 计算基于 S^T S 图度量的 `AI_measure`，
当 `AI_measure[i] >= mai_threshold` 时将该 F 点标记为 EI 点，
采用长距离插值；否则使用经典修正插值。

## 新方法原理

新增函数 `jx_PAMGMeasureEI_ideal`（`csrc/amg/par_measure_ai.c`），
基于理想插值后验误差思想选择 EI 点。

对每个 F 点 i，计算误差指示子：

```
eta_i = (|a_ii| - sum_{j in C_i} |a_ij|) / |a_ii|
```

其中 C_i 是 i 的强连接 C 点邻居集合（含跨处理器连接 S_offd 部分）。

物理含义：
- eta_i ~ 0：C 点覆盖了行的大部分量，经典插值已足够
- eta_i ~ 1：C 点覆盖不足，需要长距离插值（EI）

## 实现策略

### 替换模式（当前实现）

粗化阶段使用 `jx_PAMGMeasureAI`（基本 AI 度量），粗化后（此时 CF_marker 可用）调用
`jx_PAMGMeasureEI_ideal` 计算 eta_i，**替换** `AI_measure` 为 eta_i 值，
供后续插值阶段使用。

优点：不依赖 AI 粗化路径，保持粗化结果不变，仅影响插值选择。

## 使用方式

运行时参数 `-ai_ipt 2` 启用，例如：

```
./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 2 -ipt 60
```

阈值在 `par_amg_setup.c` 中 `mai_threshold` 控制（ipt_type=2 默认 0.9）。

## 阈值调参结果

### 测试环境
- 矩阵：Jump 128x128x128，2,097,152 行，14,581,760 非零元
- 求解器：PAMG-GMRES(30)
- 进程数：1
- 粗化：PMIS，ai_measure_type=1

### 不同阈值对比

| 阈值 | Level 0 EI | Setup | Solve | 总时间 | 算子复杂度 |
|---|---|---|---|---|---|
| Baseline(S^T S) | 17,293 (1.2%) | 16.55s | 26.21s | 42.78s | 2.557 |
| 0.3 | 1,283,262 (88.6%) | 35.77s | 19.83s | 55.69s | 2.483 |
| 0.5 | 1,283,262 (88.6%) | 35.77s | 19.83s | 55.69s | 2.483 |
| 0.7 | 349,057 (24.1%) | 32.46s | 22.16s | 54.71s | 1.371 |
| **0.9** | **0 (0%)** | **13.54s** | 28.41s | **41.99s** | 1.387 |

### 分析

1. **阈值 0.3-0.5**：几乎所有 F 点 eta_i 都超过 0.3，EI 点爆炸，Setup 翻倍
2. **阈值 0.7**：EI 降至 24%，仍有大量长距离插值
3. **阈值 0.9**：Level 0 无 EI 点（全经典插值），深层少量 EI。总时间最优

### 最佳结果（threshold=0.9）

| 指标 | Baseline | 新方法 (threshold=0.9) | 变化 |
|---|---|---|---|
| Setup | 16.55s | **13.54s** | -18.2% |
| Solve | 26.21s | 28.41s | +8.4% |
| **总时间** | **42.78s** | **41.99s** | **-1.8%** |
| EI 总数 | 459,611 | 60,239 | -86.9% |
| 算子复杂度 | 2.557 | **1.387** | -45.8% |
| 网格复杂度 | 1.385 | 1.387 | 持平 |

## 代码位置

| 文件 | 内容 |
|---|---|
| `csrc/amg/par_measure_ai.c` | `jx_PAMGMeasureEI_ideal` 函数实现 |
| `csrc/amg/par_amg_setup.c` | ipt_type==2 分支：粗化后替换 AI_measure |
| `include/jx_pamg.h` | 函数声明 |

## 分支

功能在 `feature/ideal-interp-adaptive` 分支上开发。

## 测试命令

```
# Baseline (原有 AHI)
yhrun --mpi=pmix -n 1 ./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 1 -ipt 60

# 新方法 (阈值 0.9)
yhrun --mpi=pmix -n 1 ./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 2 -ipt 60
```
"""

import os
fpath = "/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/AHI_EI_IDEAL_README.md"
with open(fpath, "w") as f:
    f.write(content)
print("Written", os.path.getsize(fpath), "bytes")
