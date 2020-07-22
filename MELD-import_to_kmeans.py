#!/usr/bin/env python
# coding: utf-8

# #### Load libraries, data and set input/output paths for local and cluster environments

# In[1]:


import os
import sklearn
import pickle
import pandas as pd
import numpy as np
import graphtools as gt
import phate
import magic
import scprep
import meld
print(meld.version)
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

np.random.seed(42)
font = {'size'   : 14}
mpl.rc('font', **font)

if os.path.isdir('/Users/mac/cloudstor/') == True:
    place = 'local'
    indir = '/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/'
    outdir = '/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data/'
    dataset = 'practice'
else:
    place = 'wolfpack'
    indir = '/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/'
    outdir = '/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data/'
    dataset = 'all'

print('Place: ', place,
       'indir: ', indir,
       'outdir: ', outdir)


# Load matrix of scaled, batch-corrected, SCT-normalised pearson-residuals from Seurat (for Graph construction)

# In[2]:


corrected_SCT = scprep.io.load_csv(os.path.join(indir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_SCT.csv'))
print(corrected_SCT.head())


# Load matrix of scaled, log-normalised, batch-corrected counts from Seurat (for imputation)

# In[3]:


corrected_logNorm = scprep.io.load_csv(os.path.join(indir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_LogNorm.csv'))
print(corrected_logNorm.head())


# Load cell-data (to become metadata object)

# In[4]:


colData = scprep.io.load_csv(os.path.join(indir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_colData.csv'))


# Create metadata object from cell data

# In[5]:


metadata = colData[['Sample', 'Tissue', 'Replicate']].astype(str)
str_reps = metadata['Replicate'].astype(str)
metadata.loc[:, ('Replicate')] = str_reps
print(metadata.head())


# #### Calculate Graph object for PHATE embedding, MAGIC imputation and later, MELD signals

# Create Graph object based on scaled, normalised, batch corrected expression values (SCTranform or vst normalised pearson residuals)

# In[6]:


G = gt.Graph(corrected_SCT, knn=10, decay=40, n_pca=100, use_pygsp=True, n_jobs=-2, verbose=True, random_state=42)


# In[7]:


pickle.dump(G, open(os.path.join(outdir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_GRAPH.pkl'), "wb"), protocol=4)


# Use Graph to calculate PHATE embeddings in 2 dimensions

# In[8]:


phate_op = phate.PHATE(knn_dist='precomputed', t = 'auto', gamma=1, n_jobs=-2, n_components=2, random_state=42)
data_phate = phate_op.fit_transform(G.kernel)


# In[9]:


pickle.dump(data_phate, open(os.path.join(outdir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_PHATE.pkl'), "wb"), protocol=4)


# Add phate embeddings to metadata object

# In[10]:


df_data_phate = pd.DataFrame(data = data_phate[0:,0:], columns = np.array(['PHATE1', 'PHATE2']), index = corrected_SCT.index)
metadata = pd.concat([metadata, df_data_phate], axis=1)
print(metadata.head())


# Impute log-normalised counts based on Graph

# In[11]:


magic_op = magic.MAGIC(knn=G.knn, decay=G.decay, t = 'auto', n_jobs=-2, random_state=42)
data_magic = magic_op.fit_transform(corrected_logNorm, graph=G)


# In[12]:


pickle.dump(data_magic, open(os.path.join(outdir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_MAGIC.pkl'), "wb" ), protocol=4)

data_magic.to_csv(os.path.join(outdir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_MAGIC_counts.csv'), index = True, header = True)


# #### Visualise PHATE embedding and perform kmeans clustering to establishing starting clusters for MELD analysis

# Create color maps for each experimental grouping (by replicate/tissue/sample)

# In[13]:


rep_cmap = {'1': '#42a5f0',
            '2': '#a5f042',
            '3': '#e442f0',
            '4': '#f0424e'}

tissue_cmap = {'Primary': '#E69F00',
               'LN': '#56B4E9',
               'Liver': '#009E73',
               'Lung': '#F0E442'}

sample_cmap = {'Primary_1': '#E69F00',
               'Primary_2': '#ffb000',
               'Primary_3': '#ffb81a',
               'Primary_4': '#ffc034',
               'LN_1': '#56B4E9',
               'LN_2': '#6dbeec',
               'LN_3': '#83c8ef',
               'LN_4': '#9ad2f2',
               'Liver_1': '#009E73',
               'Liver_2': '#00b886',
               'Liver_3': '#00d198',
               'Liver_4': '#00ebab',
               'Lung_1': '#F0E442',
               'Lung_2': '#f2e75a',
               'Lung_3': '#f4eb71',
               'Lung_4': '#f6ee89'}

cmaps = {'rep_cmap': rep_cmap, 'tissue_cmap': tissue_cmap, 'sample_cmap': sample_cmap}
pickle.dump(cmaps, open(os.path.join(outdir, 'MELD_cmaps.pkl'), "wb" ))


# In[14]:


fig, ax = plt.subplots(1)
fig.set_size_inches(10,10)
scprep.plot.scatter2d(data_phate, c=metadata['Sample'], cmap=cmaps['sample_cmap'], ticks=False, ax=ax, title='Sample')
fig.tight_layout()
fig.savefig(os.path.join(outdir, 'Sample_labeled_PHATE.png'), dpi=300)


# In[15]:


fig, ax = plt.subplots(1)
fig.set_size_inches(10,10)
scprep.plot.scatter2d(data_phate, c=metadata['Tissue'], cmap=cmaps['tissue_cmap'], ticks=False, ax=ax, title='Tissue')
fig.tight_layout()
fig.savefig(os.path.join(outdir, 'Tissue_labeled_PHATE.png'), dpi=300)


# In[16]:


fig, ax = plt.subplots(1)

groups, counts = np.unique(metadata['Sample'], return_counts=True)
for i, c in enumerate(counts):
    ax.bar(i, c, color=cmaps['sample_cmap'][groups[i]])
    
ax.set_xticks(np.arange(i+1))
ax.set_xticklabels(groups)
ax.set_ylabel('# cells')

fig.tight_layout()

fig.savefig(os.path.join(outdir, 'Cells_per_sample.png'), dpi=300)


# Perform k-means clustering on the PHATE operator. At this stage its better to under-cluster than over-cluster as VCF will partition clusters with more information.

# In[17]:


for k in range(3, 20):
    clusters = phate.cluster.kmeans(phate_op, n_clusters=k, random_state=42)
    metadata['Kmeans_' + str(k) + '_cluster'] = clusters

    ax = scprep.plot.scatter2d(data_phate, c=clusters, cmap=sns.husl_palette(5), s=1,
                      figsize=(4.3,4), ticks=None, label_prefix='PHATE',
                     legend_anchor=(1,1), fontsize=12, title='PHATE ' + str(k) + ' clusters')
    fig = ax.figure
    fig.savefig(os.path.join(outdir, 'PHATE_' + str(k) + '_kmeans.png'), dpi=300)


# #### Save metadata object, inspect K-means plots to determine optimum starting clusters for MELD

# In[18]:


pickle.dump(metadata, open(os.path.join(outdir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_MELD_metadata_Kmeans.pkl'), "wb"), protocol=4)

metadata.to_csv(os.path.join(outdir, 'Prefiltered_experiment_' + str(dataset) + '_seurat_integrated_MELD_metadata_Kmeans.csv'), index = True, header = True)

