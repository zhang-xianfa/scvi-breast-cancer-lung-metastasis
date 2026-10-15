import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.neighbors import NearestNeighbors
import scvi

np.random.seed(123)

metadata = pd.read_csv('metadata_with_superclusters.csv', index_col=0)
df = pd.read_csv('DATA/AfterQC_count.csv', index_col=0)
adata = sc.AnnData(X=df.T)

# 数据预处理
sc.pp.filter_genes(adata, max_counts=1e5)
sc.pp.filter_cells(adata, min_genes=200)
adata.layers["counts"] = adata.X.copy()

adata.obs['model'] = metadata['model']
adata.obs['Phase'] = metadata['Phase']
adata.obs['batch'] = adata.obs['model'].astype(str) + "_" + adata.obs['Phase'].astype(str)

# SCVI模型设置
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key='batch',
    layer="counts"
)

model = scvi.model.SCVI(
    adata,
    n_layers=5,  # Start simpler (2-3 layers usually sufficient)
    n_latent=20,
    gene_likelihood="zinb",  # Zero-inflated negative binomial often better

    use_layer_norm="both"  # Helps with stability
)
model.train(
    max_epochs=10000,  # Start with fewer epochs
    batch_size=512,  # Larger batch size often helps
    early_stopping=True,
    early_stopping_monitor="elbo_train",  # Monitor training ELBO
    early_stopping_patience=30,
    plan_kwargs={
        "lr": 1e-4,  # More standard learning rate
        "weight_decay": 0.1  # Stronger regularization
    },
    check_val_every_n_epoch=1,
    validation_size=0.1
)
latent = model.get_latent_representation()
adata.obsm['X_scvi'] = latent


sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added='seurat_clusters')


def calculate_kbet(adata, batch_key='model', n_neighbors=15):
    X = adata.obsm['X_scvi']
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(X)
    _, indices = nbrs.kneighbors(X)
    batch = adata.obs[batch_key].values
    same_batch = (batch[:, None] == batch[indices]).mean(axis=1)
    return same_batch.mean()

def calculate_lisi(adata, vars, n_neighbors=30):
    X = adata.obsm['X_scvi']
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(X)
    _, indices = nbrs.kneighbors(X)
    scores = {}
    for var in vars:
        categories = adata.obs[var].astype('category')
        neighbor_cats = categories.values[indices]
        diversity = [len(np.unique(row)) for row in neighbor_cats]
        scores[var] = np.mean(diversity)
    return scores


results = {
    'Method': 'SCVI',
    'kBET': calculate_kbet(adata),
    'LISI_model': calculate_lisi(adata, ['model'])['model'],
    'LISI_Phase': calculate_lisi(adata, ['Phase'])['Phase']
}

clusters = adata.obs['seurat_clusters'].astype('category').cat.codes
X = adata.obsm['X_scvi'][:, :30]
results['Silhouette'] = silhouette_score(X, clusters)
results['Calinski'] = calinski_harabasz_score(X, clusters)
results['N_Clusters'] = len(adata.obs['seurat_clusters'].unique())

pd.DataFrame([results]).to_csv('SCVI_correction_metrics.csv', index=False)

fig, axes = plt.subplots(1, 3, figsize=(24, 6))

sc.pl.umap(adata, color='model', ax=axes[0], show=False, palette='tab20')
axes[0].set_title('SCVI Corrected by Model')

sc.pl.umap(adata, color='Phase', ax=axes[1], show=False, palette='tab20')
axes[1].set_title('SCVI Corrected by Phase')

sc.pl.umap(adata, color='seurat_clusters', ax=axes[2], show=False, palette='tab20')
axes[2].set_title('SCVI Clusters (Leiden)')

plt.tight_layout()
plt.savefig('SCVI_correction_umap.png')
plt.close()

plt.figure(figsize=(10, 5))
plt.plot(model.history['elbo_train'], label='Training ELBO')
plt.xlabel('Epochs')
plt.ylabel('ELBO')
plt.title('SCVI Training Convergence')
plt.savefig('SCVI_training_curve.png')
plt.close()

pd.DataFrame(latent, index=adata.obs_names).to_csv("scvi_latent.csv")
adata.obs.to_csv('metadata_with_SCVI_clusters.csv')
