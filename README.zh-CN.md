# mrioTensor：基于张量表示的多区域投入产出分析

[![R-CMD-check](https://github.com/MaEconStat/mrioTensor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MaEconStat/mrioTensor/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

`mrioTensor` 是一个基于 R 原生（base-R）实现的多区域投入产出（MRIO）张量计算与环境核算工具包。它将传统的“国家—部门”复合二维矩阵转化为显式的四阶张量，提供张量模式乘积、Leontief 总产出求解、隐含环境流归因核算以及高阶奇异值分解（HOSVD）。

---

## 选题背景

在国民经济核算与环境扩展投入产出分析（EEIOA）中，全球多区域投入产出表通常以二维分块大矩阵 $Z \in \mathbb{R}^{N \times N}$ 进行存储与计算。每一个复合索引 $i = (r - 1)S + s$ 将空间地理维度（国家/地区 $r$）与产业生产技术维度（部门 $s$）合并。

二维拉平虽然便于传统的矩阵求逆，但抹平了全球产业链的多维网络特征：

1. **结构语义混淆**：跨国贸易流动（国家 $\to$ 国家）与产业间供需依赖（部门 $\to$ 部门）被揉合在单一维度中，难以正交解耦；
2. **多线性操作繁琐**：若需剥离特定维度的影响（如汇总所有产业得到纯国家间宏观贸易转移矩阵，或汇总所有国家得到全球通用生产技术矩阵），二维切片操作复杂；
3. **高阶拓展受限**：当进一步纳入多重环境要素（水、碳、能、土地）或跨年份面板数据时，数据天然构成 5 阶或 6 阶张量，亟需多线性代数（Multilinear Algebra）工具的支持。

`mrioTensor` 旨在打通投入产出分析与高阶张量代数之间的桥梁，提供规范、轻量且无外部第三方依赖的 R 语言基础工具。

---

## 核心接口

| 函数 | 说明 | 数学定义 |
| :--- | :--- | :--- |
| `as_mrio_tensor()` | 二维复合矩阵转换为四阶张量 | $\mathcal{T} \leftarrow Z$ via $i = (r-1)S + s$ |
| `as_mrio_matrix()` | 四阶张量无损还原为二维复合矩阵 | $Z \leftarrow \mathcal{T}$ via $i = (r-1)S + s$ |
| `aggregate_mrio()` | 沿空间或产业维度进行边际汇总 | $\sum_{s, s'} \mathcal{T}$（区域）或 $\sum_{r, r'} \mathcal{T}$（部门） |
| `mode_product()` | 张量与矩阵的模式乘积（Mode-$n$ 缩并） | $\mathcal{Y} = \mathcal{X} \times_n M$ |
| `leontief_output()` | 投入产出线性方程总产出求解 | $(I - A) x = y$ |
| `embodied_transfer()` | 最终需求归因的环境溢出矩阵 | $F = \hat{f} (I - A)^{-1} \hat{y}$ |
| `path_order_contributions()` | 生产层级（Neumann 级数幂次）贡献分解 | $c_k = f^\top A^k y \quad (k = 0, 1, \dots)$ |
| `mrio_hosvd()` | 截断高阶奇异值分解（Tucker 分解） | $\mathcal{T} \approx \mathcal{G} \times_1 U_1 \times_2 U_2 \times_3 U_3 \times_4 U_4$ |
| `example_mrio()` | 两国两部门基准算例 | 用于单元测试与算法验证 |

---

## 快速上手

```r
library(mrioTensor)

# 加载示例数据（2 个地区 x 2 个部门）
e <- example_mrio()

# 1. 2D 矩阵与 4D 张量无损转换
T_arr <- as_mrio_tensor(e$A, regions = e$regions, sectors = e$sectors)
dim(T_arr)
# [1] 2 2 2 2
# 维度顺序：(流出地区, 流出部门, 流入地区, 流入部门)

# 2. 边际聚合：分别提取纯区域转移矩阵与纯部门技术矩阵
reg_flows  <- aggregate_mrio(T_arr, by = "region")  # 2 x 2 国家间流动
sec_flows  <- aggregate_mrio(T_arr, by = "sector")  # 2 x 2 部门间流动

# 3. Leontief 总产出求解：(I - A) x = y
x <- leontief_output(e$A, e$y)

# 4. 隐含环境流最终需求归因核算
account <- embodied_transfer(
  A = e$A, y = e$y, intensity = e$intensity,
  regions = e$regions, sectors = e$sectors
)
# 验证全要素足迹恒等式：sum(F) == f' x
stopifnot(isTRUE(all.equal(account$total, sum(e$intensity * x))))

# 5. 生产层级 Neumann 级数幂次分解
path_order_contributions(e$A, e$y, e$intensity, max_order = 4)

# 6. 张量截断高阶 SVD (HOSVD)
fit <- mrio_hosvd(T_arr, ranks = c(2, 2, 2, 2))
cat("相对拟合误差:", fit$relative_error, "\n")
```

---

## 数学表达与核算恒等式

### 1. 张量索引与排布约定

设网络包含 $R$ 个地区，每个地区包含 $S$ 个行业部门，复合节点总数 $N = RS$。在标准的**地区优先、部门次之**（region-major, sector-minor）排布下：

$$
i = (r - 1)S + s, \quad j = (r' - 1)S + s'
$$

四阶张量 $\mathcal{T} \in \mathbb{R}^{R \times S \times R \times S}$ 与拉平后的中间使用/直接消耗矩阵 $Z \in \mathbb{R}^{N \times N}$ 满足严格的一一对应关系：

$$
\mathcal{T}_{r, s, r', s'} = Z_{(r-1)S+s, \, (r'-1)S+s'}
$$

### 2. Leontief 供需平衡方程

给定直接消耗系数矩阵 $A = Z \hat{x}^{-1}$ 及最终需求向量 $y$，经典投入产出平衡方程：

$$
x = A x + y \iff (I - A)x = y
$$

在包内通过 LU 分解（`solve(I - A, y)`）直接求解线性方程组，避免显式计算耗时的逆矩阵。

当谱半径 $\rho(A) < 1$ 时，Leontief 逆矩阵满足收敛的 Neumann 级数展开：

$$
L = (I - A)^{-1} = \sum_{k=0}^{\infty} A^k
$$

生产层级贡献 $c_k = f^\top A^k y$ 采用序列矩阵-向量递推乘法 $v_k = A v_{k-1}$（$v_0 = y$）进行评估，时间复杂度控制在 $O(k N^2)$。

### 3. 最终需求隐含环境溢出归因

设行业直接环境强度向量为 $f \in \mathbb{R}^N$（如单位总产出的碳排放量），最终需求驱动的双边环境溢出归因矩阵 $F \in \mathbb{R}^{N \times N}$ 为：

$$
F = \hat{f} (I - A)^{-1} \hat{y}
$$

分量形式即为：

$$
F_{ij} = f_i L_{ij} y_j
$$

行和 $\sum_j F_{ij}$ 刻画生产侧直接排放责任，列和 $\sum_i F_{ij}$ 刻画消费侧完全环境足迹，系统总足迹满足恒等式：

$$
\sum_{i=1}^N \sum_{j=1}^N F_{ij} = f^\top x
$$

---

## 计算规模说明

* **内存消耗**：四阶数组包含 $R^2 S^2 = N^2$ 个双精度浮点数。以典型中型投入产出表（$R=189, S=26 \implies N=4,914$）为例，单个密集数组占用内存约 $192\text{ MB}$。
* **大规模求逆策略**：在超大规模网络下，应优先使用 `leontief_output()` 直接求解线性方程组，避免显式计算密集的 $L = (I - A)^{-1}$。`embodied_transfer()` 内置了规模安全保护阈值（默认 `max_nodes = 1000L`）。

---

## 关联项目

* **[LN-MMSBM](https://github.com/MaEconStat/ln-mmsbm-python)**：基于正则化张量分解（TWIST）的有向加权多层网络社区结构识别 Python 原型。

---

## 许可证

MIT License. Copyright (c) 2026 Zhen Ma (MaEconStat).
