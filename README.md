# 基于 scVI 的单细胞转录组批次效应整合 —— 乳腺癌肺转移机制挖掘

> 第十一届全国大学生统计建模大赛 ｜ 省级二等奖 ｜ 第一负责人
> Python（scVI / Scanpy）+ R（Harmony / ComBat）

多批次单细胞转录组数据之间存在**批次效应**，直接合并会把「样本差异」误当成「生物学差异」。
本项目用**变分自编码器 scVI** 做批次整合与潜空间降维，并与多种经典方法做定量对比，
再在整合后的潜空间上做聚类、注释与通路富集，挖掘与乳腺癌**肺转移**相关的细胞亚群。

---

## 一、数据

| 项目 | 内容 |
|---|---|
| 数据集 | **GSE202695**（GEO 公开数据） |
| 规模 | **1523 个**肿瘤细胞（原发瘤 Primary + 肺转移瘤 Lung metastasis） |
| 批次 | 5 个样本批次（HBRX2353 / MDAMB231 / HBRX3078 / HBRX1921 / HBRX2344） |
| 元数据 | 细胞类型、细胞周期 Phase（G1/S/G2M）、PDG 分级（L1~L4） |

> 原始表达矩阵体积较大（>100 MB），**未随仓库分发**。
> 复现时请从 GEO 下载 `GSE202695`，或按下文「运行」一节的说明准备目录结构。

---

## 二、方法

### 1. 批次效应整合

- **scVI**（single-cell Variational Inference）：以负二项分布建模计数、用变分自编码器
  学习低维潜空间，把批次作为条件变量，从而在潜空间中消除批次效应。
  配置：`n_layers=5`、`n_latent=20`、`gene_likelihood=zinb`、`use_layer_norm=both`、
  `max_epochs=10000`、`batch_size=512`、`early_stopping_patience=30`、`lr=1e-4`
- **对照方法**：ComBat、Harmony、CCA、fastMNN（R 实现）

### 2. 批次整合效果的定量评价

用 `kBET`（批次混合检验，越低越好）、`LISI`（局部逆辛普森指数，越高越好）、
`Silhouette`、`Calinski-Harabasz` 四项指标综合评价：

| Method | kBET | LISI_model | LISI_Phase | Silhouette | Calinski |
|---|---|---|---|---|---|
| ComBat | 0.972 | 1.054 | 1.622 | 0.193 | 217.292 |
| Harmony | 0.967 | 1.055 | 1.632 | 0.173 | 211.034 |
| CCA | 0.676 | 2.222 | 1.896 | 0.044 | 59.384 |
| fastMNN | 0.693 | 2.034 | 1.995 | 0.076 | 132.266 |
| **scVI** | **0.549** | **4.121** | **2.988** | **0.195** | 108.783 |

**结论：scVI 在 kBET 与 LISI 两项批次混合指标上均为最优**，说明其整合后不同批次的细胞
混合得最充分；同时保留了较好的聚类可分性（Silhouette 最高）。

### 3. 下游分析

```
潜空间 → Leiden 聚类（resolution=0.5）→ 12 类细胞
       → GSEA 超聚类归并 → A / B / C 三个超级细胞簇
       → 差异表达分析 + GO / GSEA 通路富集 → 肺转移相关亚群与特征通路
       → Monocle3 拟时序分析（K=20 近邻，A 类为起点）→ 关键基因 MALAT1 / ZEB1 梯度表达
```

---

## 三、文件说明

```
code/
  main.py                        论文复现主流程：QC → scVI 训练 → 聚类 → 指标评价
  scvi_analysis.py               scVI 建模与潜空间提取
  batch_correction_baselines.R   ComBat / Harmony / CCA / fastMNN 对照实验
  gsea_go_enrichment.py          差异表达 + GO / GSEA 富集分析
  pseudotime_analysis.py         Monocle3 拟时序分析
results/
  scvi-400epochs/                scVI 400 epoch 的 t-SNE 可视化（按 Model / Type / Phase / PDG 着色）
  scvi-400epochs/metrics.txt     对应的四项聚类/混合指标
  reproduction/                  论文复现实验的聚类结果图与指标（含 Louvain 聚类对比）
```

## 四、运行

```bash
pip install scvi-tools scanpy anndata matplotlib pandas scikit-learn
# 将 GSE202695 的表达矩阵与 metadata 放到 code/ 同级的数据目录后：
python code/main.py
```

R 脚本需要 `harmony`、`sva`（ComBat）、`batchelor`（fastMNN）等包。

---

## 说明

本项目为竞赛参赛作品，代码由参赛时的分析脚本整理而成（部分原始文件为 `.txt`，
在此按用途重命名为 `.py` / `.R`，内容未改动）。论文与证书等材料未随仓库分发。
