# Homepage project entry · 中文稿

## mrioTensor

**把多区域投入产出模型中的复合节点矩阵，还原为具有经济含义的多维数据结构。**

传统 MRIO 计算通常以“国家—部门”为复合节点组织矩阵。`mrioTensor` 将其转换成
国家 × 部门 × 国家 × 部门四阶数组，在保留既有投入产出恒等式的同时，
使国家与部门维度能够独立索引、汇总和参与多线性运算。

这个 R 包目前公开的是方法基础设施：标准 Leontief 线性方程、
最终需求归因、张量模态乘积，以及小规模 HOSVD 示例。
所有演示使用合成数据，不涉及未发表的实证结果或专有算法。

**关键词**：R package · Input–Output Economics · Tensor Algebra · Reproducible Software

**关联项目**：[LN-MMSBM](https://github.com/MaEconStat/ln-mmsbm)（Python，多层网络统计模型）

## Technical snapshot · English

*An inspectable R toolkit connecting the familiar region–sector matrix
representation of MRIO systems with explicit fourth-order tensors. It includes
standard Leontief accounting, exact index transformations, mode-wise algebra,
and synthetic, reproducible examples. This public release deliberately
excludes unpublished empirical results.*

## Suggested homepage presentation

- Card title: `mrioTensor`
- Category: `Methods / Research software`
- Brief value: `From composite MRIO matrices to semantically explicit tensors`
- Primary link: GitHub repository (after publication)
- Secondary link: Related `ln-mmsbm` project
- Technical detail page: show a **verified** 2D→4D indexing diagram,
  public function interfaces, and one synthetic `2 regions × 2 sectors`
  worked example. Avoid benchmark claims without measurements.
