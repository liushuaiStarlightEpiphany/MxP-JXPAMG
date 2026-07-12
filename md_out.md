# 基于经典修正插值 + 理想插值后验误差的自适应插值方法

## 背景

以 MSI 经典修正插值（interp_type=0）为基础，通过理想插值后验误差检测精度不足的节点，
在这些节点上局部替换为长距离插值（扩展插值），从而得到比纯经典插值更好、比 AHI 更高效的算子。

## 核心思想

对每个 F 点 i，计算误差指示子：

```
eta_i = (|a_ii| - sum_{j in C_i} |a_ij|) / |a_ii|
```

- eta_i ~ 0：C 点覆盖了行的大部分量，经典插值已足够
- eta_i ~ 1：C 点覆盖不足，需要长距离插值

当 eta_i >= threshold 时，对该点使用扩展插值（包含 F-F 耦合的三点路径），
否则使用经典修正插值。

## 实现方式

新增 interp_type=61，独立于 AHI（interp_type=60）：

1. 粗化阶段使用 jx_PAMGMeasureAI（基本 AI 度量）
2. 粗化后（CF_marker 可用时）调用 jx_PAMGMeasureEI_ideal 计算 eta_i
3. 替换 AI_measure 为 eta_i 值
4. 插值阶段根据 AI_measure[i] >= mai_threshold 决定 EI 点

## 使用方式

```
./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 2 -ipt 61
```

阈值由 mai_threshold 控制（ipt_type=2 默认 0.9）。

## 阈值调参结果

### 测试环境
- 矩阵：Jump 128x128x128，2,097,152 行，14,581,760 非零元
- 求解器：PAMG-GMRES(30)
- 进程数：1
- 粗化：PMIS，ai_measure_type=1

### 不同阈值对比

| 阈值 | Level 0 EI | Setup | Solve | 总时间 | 算子复杂度 |
|---|---|---|---|---|---|
| Baseline (interp=60) | 17,293 (1.2%) | 16.55s | 26.21s | 42.78s | 2.557 |
| 0.3 | 1,283,262 (88.6%) | 35.77s | 19.83s | 55.69s | 2.483 |
| 0.5 | 1,283,262 (88.6%) | 35.77s | 19.83s | 55.69s | 2.483 |
| 0.7 | 349,057 (24.1%) | 32.46s | 22.16s | 54.71s | 1.371 |
| **0.9** | **0 (0%)** | **13.61s** | 28.51s | **42.12s** | 1.387 |

### 最佳结果（threshold=0.9）

| 指标 | Baseline (interp=60) | 新方法 (interp=61, th=0.9) | 变化 |
|---|---|---|---|
| Setup | 16.55s | **13.61s** | -17.8% |
| Solve | 26.21s | 28.51s | +8.8% |
| **总时间** | **42.78s** | **42.12s** | **-1.5%** |
| EI 总数 | 459,611 | 60,239 | -86.9% |
| 算子复杂度 | 2.557 | **1.387** | -45.8% |
| 网格复杂度 | 1.385 | 1.387 | 持平 |

Level 0 无 EI 点（全经典插值），深层少量 EI。Setup 加速显著。

## 代码位置

| 文件 | 内容 |
|---|---|
| csrc/amg/par_measure_ai.c | jx_PAMGMeasureEI_ideal 函数实现 |
| csrc/amg/par_amg_setup.c | interp_type=61 分支 + ipt_type=2 替换逻辑 |
| include/jx_pamg.h | 函数声明 |

## 分支

功能在 feature/ideal-interp-adaptive 分支上开发。

## 测试命令

```
# Baseline (原有 AHI)
yhrun --mpi=pmix -n 1 ./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 1 -ipt 60

# 新方法 (interp=61, 阈值 0.9)
yhrun --mpi=pmix -n 1 ./solver_strong -nts 1 -pid 3 -sid 22 -ai_mt 1 -ai_ipt 2 -ipt 61
```
