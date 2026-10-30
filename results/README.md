# results/ 说明：每个数字从哪来

本目录只保存**能直接跑出来、体积又小**的产物（指标文本 + t-SNE 图）。
原始表达矩阵超过 100 MB 未随仓库分发，因此部分指标无法在当前仓库内直接复算，
这里逐项标明来源，避免「表里的数字对不上文件」的尴尬。

## 目录

| 路径 | 内容 | 怎么产生的 |
|---|---|---|
| `scvi-400epochs/metrics.txt` | scVI 400 epoch 版本的 Silhouette / Calinski-Harabasz / Davies-Bouldin（按 MODEL / TYPE / PDG / PHASE 四个标签分别计算） | `code/main.py` 中的聚类指标函数，读完潜空间后直接打印并写盘 |
| `scvi-400epochs/tsne_*.png` | 同一版本的 t-SNE，按 Model / Type / Phase / PDG 着色 | `code/scvi_analysis.py` |
| `reproduction/metrics.txt` | 论文复现流程（完整训练 + Louvain 聚类）的同三项指标 | `code/main.py` |
| `reproduction/*.png` | 复现流程的 t-SNE 与 Louvain 聚类图 | 同上 |

## 主 README 里那张对照表的出处

主 README「批次整合效果的定量评价」表格包含 kBET / LISI / Silhouette / Calinski 四项，
方法为 ComBat、Harmony、CCA、fastMNN、scVI：

- **kBET / LISI**：由 `code/batch_correction_baselines.R` 计算
  （R 包 `kBET` 的 `kBET()` 与 `lisi` 的 `compute_lisi()`），
  Python 侧 `code/scvi_analysis.py` 里也有等价实现 `calculate_kbet` / `calculate_lisi`；
- **Silhouette / Calinski**：由 `code/main.py` 在潜空间上计算。

⚠️ 这两个指标（kBET/LISI）的**原始输出文本没有随仓库保存**，表里的数值来自当时
完整跑完对照实验后的结果。要复核需要：
1. 按主 README 第四节准备 `GSE202695` 数据；
2. 先跑 `code/main.py` 得到 scVI 潜空间，再跑 `code/batch_correction_baselines.R`
   （R 环境需装 `kBET`、`lisi`、`batchelor`、`harmony` 等包），最后汇总成同一张表。

## 为什么两个 metrics.txt 的数不一样

`scvi-400epochs/` 是**训练 400 epoch 的中间版本**，`reproduction/` 是完整训练流程。
两者标签体系与聚类算法（KMeans vs Louvain）不同，指标不可直接横向比较，
更不代表「400 epoch 更好或更差」。对外只引用主 README 那张完整对照表。
