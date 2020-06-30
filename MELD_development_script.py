# -*- coding: utf-8 -*-
"""
MELD test script

This is composed using the development dataset.
"""

import pandas as pd
import numpy as np
import graphtools as gt
import phate
import magic
import scprep
import meld
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(42)

corrected_SCT = scprep.io.load_csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_SCT.csv")
corrected_logNorm = scprep.io.load_csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_LogNorm.csv")
colData = scprep.io.load_csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_colData.csv")
rowData = scprep.io.load_csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_rowData.csv")
pca = scprep.io.load_csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_pca_embeddings.csv")
umap = scprep.io.load_csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_umap_embeddings.csv")

metadata = colData[['Tissue', 'whole_experiment_PHATE_clusters', 'Sample', 'Replicate']]
metadata['RES'] = np.array([-1 if tiss.startswith('Primary') else 1 for tiss in metadata['Tissue']])
metadata = metadata.rename(columns={'whole_experiment_PHATE_clusters': 'Cluster'})

G = gt.Graph(corrected_SCT, knn=10, decay=40, n_pca=100, use_pygsp=True, n_jobs=-2, verbose=True)

phate_op = phate.PHATE(knn_dist='precomputed', gamma=1, n_jobs=-2, n_components=2)
data_phate = phate_op.fit_transform(G.kernel)

samples_cdict = {'Primary': '#E69F00',
                 'LN': '#56B4E9',
                 'Liver': '#009E73',
                 'Lung': '#F0E442'}

scprep.plot.scatter2d(data_phate, c=metadata['Tissue'], legend_anchor=(1, 1), figsize=(5, 5), s=10, label_prefix='PHATE', ticks=False, cmap= samples_cdict)


fig, ax = plt.subplots(1, figsize=(4,4))
sns.kdeplot(data_phate[:,0], data_phate[:,1], n_levels=100, shade=True, cmap='inferno', zorder=0, ax=ax)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('PHATE 1', fontsize=18)
ax.set_ylabel('PHATE 2', fontsize=18)
ax.set_title('Whole experiment Primary and Mets', fontsize=20)
fig.tight_layout()

clusters = phate.cluster.kmeans(phate_op, n_clusters=3)
scprep.plot.scatter2d(data_phate, c=clusters, cmap=sns.husl_palette(5), s=1, figsize=(4.3,4), ticks=None, label_prefix='PHATE', legend_anchor=(1,1), fontsize=12, title='PHATE clusters')
clusters = phate.cluster.kmeans(phate_op, n_clusters=4)
scprep.plot.scatter2d(data_phate, c=clusters, cmap=sns.husl_palette(5), s=1, figsize=(4.3,4), ticks=None, label_prefix='PHATE', legend_anchor=(1,1), fontsize=12, title='PHATE clusters')
clusters = phate.cluster.kmeans(phate_op, n_clusters=5)
scprep.plot.scatter2d(data_phate, c=clusters, cmap=sns.husl_palette(5), s=1, figsize=(4.3,4), ticks=None, label_prefix='PHATE', legend_anchor=(1,1), fontsize=12, title='PHATE clusters')
clusters = phate.cluster.kmeans(phate_op, n_clusters=6)
scprep.plot.scatter2d(data_phate, c=clusters, cmap=sns.husl_palette(5), s=1, figsize=(4.3,4), ticks=None, label_prefix='PHATE', legend_anchor=(1,1), fontsize=12, title='PHATE clusters')
clusters = phate.cluster.kmeans(phate_op, n_clusters=7)
scprep.plot.scatter2d(data_phate, c=clusters, cmap=sns.husl_palette(5), s=1, figsize=(4.3,4), ticks=None, label_prefix='PHATE', legend_anchor=(1,1), fontsize=12, title='PHATE clusters')


magic_op = magic.MAGIC(knn=G.knn, decay=G.decay)
data_magic = magic_op.fit_transform(corrected_logNorm, graph=G)

meld_op = meld.MELD()
genotype_ees = meld_op.fit_transform(G, metadata['RES'])
genotype_ees = genotype_ees - np.mean(genotype_ees)
metadata['EES'] = genotype_ees

scprep.plot.scatter2d(data_phate, c=metadata['EES'], cmap='viridis', figsize=(6, 5), s=5, label_prefix='PHATE', ticks=False)

samples_cvec = np.array([samples_cdict[s] for s in metadata['Tissue']])


def plot_EES_in_clusters(ax, clusters, EES, c=None):

    clusters_idx = np.arange(len(set(clusters)))
    n_clusts = len(set(clusters_idx))

    # Calculate means
    means = np.zeros(n_clusts)
    for i, cl in enumerate(np.unique(clusters)):
        means[i] = np.mean(EES[clusters == cl])

    # Plotting cells
    x = clusters + np.random.normal(0, .1, len(clusters))
    y = EES
    r = np.random.choice(len(y), len(y), replace=False)
    ax.scatter(x[r], y[r], c=c[r], s=2)

    # Plotting means
    ax.scatter(np.arange(n_clusts), means, c='#cbc9ff', edgecolors='black', lw=1.5, marker='o', zorder=3, s=100)

    # Plotting vetical lines
    for i in np.unique(clusters):
        ax.axvline(i, c='k', lw=.1, zorder=0)


fig, ax = plt.subplots(1, figsize=(7, 5))
plot_EES_in_clusters(ax, clusters=metadata['Cluster'], EES=metadata['EES'], c=samples_cvec)
ax.set_xticklabels(metadata['Cluster'], rotation=45, ha='right', fontsize=14)
scprep.plot.utils.shift_ticklabels(ax.xaxis, dx=0.1)
ax.set_ylabel('Enhanced Experimental Signal', fontsize=18)
ax.tick_params(labelsize=14)
fig.tight_layout()


metadata.groupby('Cluster')






