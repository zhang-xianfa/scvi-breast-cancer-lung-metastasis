import pandas as pd
import scanpy as sc
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import leidenalg
import igraph as ig

# 用这个函数计算出来PDG并写入metadata
# def calculate_pdg(adata):
#     sc.pp.filter_cells(adata, min_genes=200)
#     sc.pp.filter_genes(adata, min_cells=3)
#
#     detected_genes = (adata.X > 0).sum(axis=1) # 对每行求和
#     total_genes = adata.n_vars
#     adata.obs['PDG'] = (detected_genes / total_genes * 100)
#     return adata

metadata = pd.read_csv('DATA/metadata_with_CCandPDG.csv', index_col=0)
df = pd.read_csv('DATA/AfterQC_count.csv', index_col=0)
adata = sc.AnnData(X=df.T)

sc.pp.normalize_total(adata, target_sum=1e7)
sc.pp.log1p(adata)

sc.tl.tsne(
    adata,
    use_rep='X',
    perplexity=30
)


'''
以下部分画出了数据处理之前的model, type, PDG, Phase的tsne图，并给出聚类指标值
'''
metadata['tSNE_1'] = adata.obsm['X_tsne'][:, 0]
metadata['tSNE_2'] = adata.obsm['X_tsne'][:, 1]
metadata['PDG_category'] = pd.qcut(metadata['PDG'], q=4, labels=['L1', 'L2', 'L3', 'L4'])

def plot_tsne(metadata, color_by, title):
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='tSNE_1', y='tSNE_2',
        hue=color_by, data=metadata,
        palette='tab20', s=50, alpha=0.8
    )
    plt.title(title)
    plt.xlabel('tSNE 1')
    plt.ylabel('tSNE 2')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()
print(metadata['PDG_category'].value_counts())
# 分别绘制 model, type, Phase pdg 的 t-SNE 图
plot_tsne(metadata, color_by='model', title='t-SNE by Model (log-normalized)')
plot_tsne(metadata, color_by='type', title='t-SNE by Type (log-normalized)')
plot_tsne(metadata, color_by='Phase', title='t-SNE by Phase (log-normalized)')
plot_tsne(metadata, color_by='PDG_category', title='t-SNE by PDG Category (L1, L2, L3, L4)')


from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score


def calculate_silhouette_score(metadata, label_column):
    labels = metadata[label_column].astype('category').cat.codes  # 将类别转换为数值
    X = metadata[['tSNE_1', 'tSNE_2']]  # 使用 t-SNE 坐标
    return silhouette_score(X, labels)


def calculate_ch_score(metadata, label_column):
    labels = metadata[label_column].astype('category').cat.codes  # 将类别转换为数值
    X = metadata[['tSNE_1', 'tSNE_2']]  # 使用 t-SNE 坐标
    return calinski_harabasz_score(X, labels)


def calculate_db_score(metadata, label_column):
    labels = metadata[label_column].astype('category').cat.codes  # 将类别转换为数值
    X = metadata[['tSNE_1', 'tSNE_2']]  # 使用 t-SNE 坐标
    return davies_bouldin_score(X, labels)


metrics = {
    'MODEL': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'model'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'model'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'model')
    },
    'TYPE': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'type'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'type'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'type')
    },
    'PDG': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'PDG_category'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'PDG_category'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'PDG_category')
    },
    'PHASE': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'Phase'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'Phase'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'Phase')
    }
}


for variable, scores in metrics.items():
    print(f"Metrics for {variable}:")
    for metric, value in scores.items():
        print(f"  {metric}: {value:.4f}")



'''
以下部分用louvain算法对数据处理前的细胞进行聚类，总共分了11个类别(我一点都看不懂这个代码在写什么:( )
'''

# 构建 KNN 图
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)  # 使用 PCA 结果构建 KNN 图
knn_graph = adata.obsp['connectivities']

# 将 KNN 图转换为 igraph 对象
sources, targets = knn_graph.nonzero()
weights = knn_graph.data.astype(np.float64)  # 显式转换为浮点数

# 创建边列表（去重处理）
seen_edges = set()
edges = []
weights_list = []
for s, t, w in zip(sources, targets, weights):
    if (s, t) not in seen_edges and (t, s) not in seen_edges:
        seen_edges.add((s, t))
        edges.append((str(s), str(t)))  # 节点ID字符串化
        weights_list.append(float(w))   # 权重列表

# 构建 igraph 图并添加权重
g = ig.Graph.TupleList(
    edges,
    directed=False,
    edge_attrs={'weight': weights_list}  # 显式传递权重
)

# Louvain 聚类执行
try:
    partition = g.community_multilevel(weights=weights_list)
except ValueError as e:
    print("Louvain 聚类失败，错误信息：", e)
    partition = g.community_multilevel()

# 结果映射
clusters = np.array(partition.membership)
print("聚类结果分布:", np.unique(clusters, return_counts=True))

# 将聚类结果添加到 metadata 中
metadata['cluster'] = pd.Categorical(clusters)
plot_tsne(metadata, color_by='cluster', title='t-SNE by Cluster (Louvain Algorithm)')

'''
以下内容，把model, type, Phase当作批次信息，利用SCVI处理批次效应
'''
import scvi

# 前面进行了对数转换，但是SCVI需要原数据，这里重新导入数据
df = pd.read_csv('DATA/AfterQC_count.csv', index_col=0)
adata = sc.AnnData(X=df.T)
adata.obs['PDG'] = metadata['PDG']
adata.obs['model'] = metadata['model']
adata.obs['type'] = metadata['type']
adata.obs['Phase'] = metadata['Phase']
adata.obs['PDG_category'] = pd.qcut(adata.obs['PDG'], q=4, labels=['L1', 'L2', 'L3', 'L4'])

adata.obs['batch'] = (
    adata.obs['model'].astype(str) + "_" +
    adata.obs['type'].astype(str) + "_" +
    adata.obs['Phase'].astype(str) + "_" +
    adata.obs['PDG_category'].astype(str)
)

# 使用 SCVI 进行批次效应校正
scvi.model.SCVI.setup_anndata(adata, batch_key='batch')  # 设置批次信息
model = scvi.model.SCVI(adata)  # 初始化 SCVI 模型
model.train(max_epochs=400)  # 训练模型

# 获取校正后的 latent representation
latent = model.get_latent_representation()
adata.obsm['X_scvi'] = latent  # 将校正后的低维表示添加到 adata

# 使用校正后的数据进行 t-SNE 降维
sc.tl.tsne(adata, use_rep='X_scvi', perplexity=30)

# 将 t-SNE 结果写入 metadata
metadata['tSNE_1'] = adata.obsm['X_tsne'][:, 0]
metadata['tSNE_2'] = adata.obsm['X_tsne'][:, 1]
adata.obs['PDG_category'] = pd.qcut(adata.obs['PDG'], q=4, labels=['L1', 'L2', 'L3', 'L4'])

plot_tsne(metadata, color_by='model', title='t-SNE by Model (SCVI corrected)')
plot_tsne(metadata, color_by='type', title='t-SNE by Type (SCVI corrected)')
plot_tsne(metadata, color_by='Phase', title='t-SNE by Phase (SCVI corrected)')
plot_tsne(metadata, color_by='PDG_category', title='t-SNE by PDG Category (L1, L2, L3, L4)')

metrics = {
    'MODEL': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'model'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'model'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'model')
    },
    'TYPE': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'type'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'type'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'type')
    },
    'PDG': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'PDG_category'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'PDG_category'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'PDG_category')
    },
    'PHASE': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'Phase'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'Phase'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'Phase')
    }
}


for variable, scores in metrics.items():
    print(f"Metrics for {variable}:")
    for metric, value in scores.items():
        print(f"  {metric}: {value:.4f}")



'''
以下部分我将用SENET改进后的SCVI进行训练
'''

from sklearn.metrics import silhouette_score
import scanpy as sc
import scvi
import torch
import torch.nn as nn
from scvi.module import VAE
from scvi.nn import FCLayers
import pandas as pd

# ---------------------------
# 自定义 Self-Attention 模块
# ---------------------------
class SelfAttention(nn.Module):
    def __init__(self, embed_size, heads):
        super(SelfAttention, self).__init__()
        self.embed_size = embed_size
        self.heads = heads
        self.head_dim = embed_size // heads

        assert (
            self.head_dim * heads == embed_size
        ), "Embedding size needs to be divisible by heads"

        self.values = nn.Linear(self.head_dim, self.head_dim, bias=False)
        self.keys = nn.Linear(self.head_dim, self.head_dim, bias=False)
        self.queries = nn.Linear(self.head_dim, self.head_dim, bias=False)
        self.fc_out = nn.Linear(heads * self.head_dim, embed_size)

    def forward(self, values, keys, query):
        N = query.shape[0]
        value_len, key_len, query_len = values.shape[1], keys.shape[1], query.shape[1]

        # Split the embedding into self.heads different pieces
        values = values.reshape(N, value_len, self.heads, self.head_dim)
        keys = keys.reshape(N, key_len, self.heads, self.head_dim)
        queries = query.reshape(N, query_len, self.heads, self.head_dim)

        values = self.values(values)
        keys = self.keys(keys)
        queries = self.queries(queries)

        # Scaled dot-product attention
        energy = torch.einsum("nqhd,nkhd->nhqk", [queries, keys])
        attention = torch.softmax(energy / (self.embed_size ** (1 / 2)), dim=3)

        out = torch.einsum("nhql,nlhd->nqhd", [attention, values]).reshape(
            N, query_len, self.heads * self.head_dim
        )

        out = self.fc_out(out)
        return out


# ---------------------------
# 自定义带 Self-Attention 的编码器
# ---------------------------
class SelfAttentionEncoder(nn.Module):
    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_hidden: int = 128,
        n_latent: int = 30,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        **kwargs
    ):
        super().__init__()
        # 原始全连接层
        self.fc_layers = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_layer_norm=True,
            **kwargs
        )

        # 添加 Self-Attention 模块
        self.attention = SelfAttention(embed_size=n_hidden, heads=8)  # 8 个头

        # 输出层
        self.mean_encoder = nn.Linear(n_hidden, n_latent)
        self.var_encoder = nn.Linear(n_hidden, n_latent)

    def forward(self, x: torch.Tensor):
        q = self.fc_layers(x)
        q = self.attention(q, q, q)  # 应用 Self-Attention
        q_m = self.mean_encoder(q)
        q_v = self.var_encoder(q)
        return q_m, q_v


# ---------------------------
# 自定义 VAE 模型
# ---------------------------
class SelfAttentionVAE(VAE):
    def __init__(
        self,
        n_input: int,
        n_hidden: int = 128,
        n_latent: int = 30,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        **kwargs
    ):
        super().__init__(
            n_input=n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **kwargs
        )
        # 替换原始编码器
        self.encoder = SelfAttentionEncoder(
            n_input=n_input,
            n_output=2 * n_latent,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate
        )


# ---------------------------
# 修改后的 SCVI 模型
# ---------------------------
class SCVIWithSelfAttention(scvi.model.SCVI):
    def __init__(
        self,
        adata: sc.AnnData,
        n_hidden: int = 128,
        n_latent: int = 30,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        **model_kwargs
    ):
        # 获取输入维度
        n_input = adata.n_vars
        super().__init__(adata, **model_kwargs)
        # 替换原始 VAE 模块
        self.module = SelfAttentionVAE(
            n_input=n_input,  # 显式传递 n_input
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **model_kwargs
        )
        self._model_summary_string = "SCVIWithSelfAttention"
        self.init_params_ = self._get_init_params(locals())


# ---------------------------
# 数据处理和模型训练
# ---------------------------
# 加载数据
metadata = pd.read_csv('DATA/metadata_with_CCandPDG.csv', index_col=0)
df = pd.read_csv('DATA/AfterQC_count.csv', index_col=0)
adata = sc.AnnData(X=df.T)

# 添加 metadata 到 adata
adata.obs['PDG'] = metadata['PDG']
adata.obs['model'] = metadata['model']
adata.obs['type'] = metadata['type']
adata.obs['Phase'] = metadata['Phase']
adata.obs['PDG_category'] = pd.qcut(adata.obs['PDG'], q=4, labels=['L1', 'L2', 'L3', 'L4'])

# 设置批次信息
adata.obs['batch'] = (
    adata.obs['model'].astype(str) + "_" +
    adata.obs['type'].astype(str) + "_" +
    adata.obs['Phase'].astype(str) + "_" +
    adata.obs['PDG_category'].astype(str)
)

# 使用 Self-Attention SCVI 进行批次效应校正
SCVIWithSelfAttention.setup_anndata(adata, batch_key='batch')
model = SCVIWithSelfAttention(
    adata,
    n_hidden=128,
    n_latent=30,
    n_layers=2,
    dropout_rate=0.2
)
model.train(max_epochs=400)

# 获取校正后的 latent representation
latent = model.get_latent_representation()
adata.obsm['X_scvi'] = latent  # 将校正后的低维表示添加到 adata

# 使用校正后的数据进行 t-SNE 降维
sc.tl.tsne(adata, use_rep='X_scvi', perplexity=30)

# 将 t-SNE 结果写入 metadata
metadata['tSNE_1'] = adata.obsm['X_tsne'][:, 0]
metadata['tSNE_2'] = adata.obsm['X_tsne'][:, 1]

# 可视化 t-SNE 结果
def plot_tsne(metadata, color_by, title):
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        x='tSNE_1', y='tSNE_2',
        hue=color_by, data=metadata,
        palette='tab20', s=50, alpha=0.8
    )
    plt.title(title)
    plt.xlabel('tSNE 1')
    plt.ylabel('tSNE 2')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()

plot_tsne(metadata, color_by='model', title='t-SNE by Model (Self-Attention SCVI corrected)')
plot_tsne(metadata, color_by='type', title='t-SNE by Type (Self-Attention SCVI corrected)')
plot_tsne(metadata, color_by='Phase', title='t-SNE by Phase (Self-Attention SCVI corrected)')
plot_tsne(metadata, color_by='PDG_category', title='t-SNE by PDG Category (L1, L2, L3, L4) Self-Attention SCVI')

# 计算聚类指标
def calculate_silhouette_score(metadata, label_column):
    labels = metadata[label_column].astype('category').cat.codes
    X = metadata[['tSNE_1', 'tSNE_2']]
    return silhouette_score(X, labels)

def calculate_ch_score(metadata, label_column):
    labels = metadata[label_column].astype('category').cat.codes
    X = metadata[['tSNE_1', 'tSNE_2']]
    return calinski_harabasz_score(X, labels)

def calculate_db_score(metadata, label_column):
    labels = metadata[label_column].astype('category').cat.codes
    X = metadata[['tSNE_1', 'tSNE_2']]
    return davies_bouldin_score(X, labels)

metrics = {
    'MODEL': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'model'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'model'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'model')
    },
    'TYPE': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'type'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'type'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'type')
    },
    'PDG': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'PDG_category'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'PDG_category'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'PDG_category')
    },
    'PHASE': {
        'Silhouette Score': calculate_silhouette_score(metadata, 'Phase'),
        'Calinski-Harabasz Index': calculate_ch_score(metadata, 'Phase'),
        'Davies-Bouldin Index': calculate_db_score(metadata, 'Phase')
    }
}

for variable, scores in metrics.items():
    print(f"Metrics for {variable}:")
    for metric, value in scores.items():
        print(f"  {metric}: {value:.4f}")













