#!/usr/bin/env python
# coding: utf-8

# # This is the notebook for QCing the data
# 
# ## Use sc_gpu/rapids_sc/scv environment for this part

# In[1]:


import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import os
import rapids_singlecell as rsc
import scvi


# In[2]:


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")


# In[3]:


import logging
#reduce logging amount
# Configure logging
logging.basicConfig(level=logging.WARNING)


# In[18]:


input_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/scrna/data/velocyto_cellbender_processed/before_qc/combined_raw/"


# In[19]:


# Process all screens
samples = [f'screen-{i}' for i in range(1, 5)]


# In[20]:


adata_list = []
for sample in samples:
    #check which sample we're using
    print(f'merging {sample}' )
    # read in object
    adata = sc.read_h5ad(input_dir + sample + '_combined.h5ad')
    adata_list.append(adata)


# In[21]:


adata_list


# In[22]:


aav = ['AAV1', 'AAV7', 'AAVDJ8', 'AAVRH10']
for adata in adata_list:
    print(sum(adata.var_names.isin(aav)))


# # Preprocessing: Initial Filtration and Doublet Removal

# ## Basic Filtering

# In[23]:


# filter out useless genes
for adata in adata_list:
    adata.var_names_make_unique
    sc.pp.filter_genes(adata, min_cells= 1)
    sc.pp.filter_cells(adata,min_genes = 200)


# In[24]:


for adata in adata_list:
    print(sum(adata.var_names.isin(aav)))


# In[25]:


#check for dupes
[adata.var_names for adata in adata_list]


# ### We will use the scvi SOLO package. For this package, we will accelerate the doublet removal process using a GPU and we will shorten down our matrix by finding the highly variable genes. To clarify, we will only be using the RAW counts

# In[26]:


def find_doublets(adata):
    sc.pp.highly_variable_genes(adata,
                                n_top_genes = 3000,
                                subset = True,
                                flavor = 'seurat_v3'
                                )
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train(accelerator = 'gpu',
              log_every_n_steps = 250)
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train(accelerator = 'gpu')
    df = solo.predict()
    df['prediction'] = solo.predict(soft = False)
    print(df)
    print(df.groupby('prediction').count())
    return df


# In[27]:


adata_list


# In[28]:


adata_list_predoubs = adata_list.copy()


# In[29]:


adata_list_predoubs


# In[30]:


df_list = []
for adata in adata_list_predoubs:
    df = find_doublets(adata)
    df_list.append(df)


# In[31]:


[print(df.groupby('prediction').count()) for df in df_list]


# In[32]:


def calc_diff(df):
    df['dif'] = df.doublet - df.singlet
    
[calc_diff(df) for df in df_list]
import seaborn as sns
[sns.displot(df[df.prediction == 'doublet'], x = 'dif') for df in df_list]


# In[33]:


doublet_list = []
def doublet_cutoff(df, cutoff = 1):
    doublet_score_cutoff = cutoff
    doublets = df[(df.prediction == 'doublet') & (df.dif > doublet_score_cutoff)]
    return doublets

[doublet_list.append(doublet_cutoff(df)) for df in df_list]


# In[34]:


[print(doublet.shape) for doublet in doublet_list]
[print(df.shape) for df in df_list]


# In[35]:


[doublet.index[1:10] for doublet in doublet_list]


# In[40]:


del adata_list


# In[42]:


adata_list = []
for sample in samples:
    #check which sample we're using
    print(f'merging {sample}' )
    # read in object
    adata = sc.read_h5ad(input_dir + sample + '_combined.h5ad')
    adata_list.append(adata)


# In[43]:


# filter out useless genes
for adata in adata_list:
    adata.var_names_make_unique
    sc.pp.filter_genes(adata, min_cells= 1)
    sc.pp.filter_cells(adata,min_genes = 200)


# In[44]:


adata_list


# In[50]:


adata_list_nodoub = []
for i,adata in enumerate(adata_list):
    adata.obs['doublet'] = adata.obs.index.isin(doublet_list[i].index)
    print(adata.obs['doublet'].value_counts())
    adata = adata[~adata.obs.doublet]
    adata_list_nodoub.append(adata)


# In[51]:


adata_list_nodoub


# In[52]:


aav


# In[53]:


[sum(adata.var_names.isin(aav)) for adata in adata_list_nodoub]


# In[60]:


output_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/scrna/data/velocyto_cellbender_processed/before_qc/doublet_removal"
for i,adata in enumerate(adata_list_nodoub):
    adata.var_names_make_unique()
    adata.write(os.path.join(output_dir, samples[i]+'-nodoublets_combined.h5ad'))


# In[61]:


def process_data(adata):
    #remove duplicate names first
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # ribosomal genes    
    adata.var["ribo"] = adata.var_names.str.startswith("Rpl") | adata.var_names.str.startswith("Rps")
    # calculate qc metrics
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo"],percent_top=[20], inplace=True, log1p=True)
    #print the violin plot of metrics
    print(sc.pl.violin(adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4, multi_panel=True,))
    print(sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt"))


# In[62]:


[process_data(adata) for adata in adata_list_nodoub]


# In[63]:


adata_list_nodoub


# In[67]:


df2 = pd.concat([x.obs for x in adata_list_nodoub])
df2 = df2.sort_values('batch')


# In[72]:


qc_output_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/scrna/results/qc/qc_metrics"


# In[73]:


adata.obs_keys


# In[74]:


qc_metrics = ['total_counts', 'log1p_total_counts',
              'n_genes', 'n_genes_by_counts', 
              'pct_counts_mt',
              'pct_counts_in_top_20_genes']


# In[78]:


def generate_qc_figs(df, value):
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    g = sns.FacetGrid(df, row="batch", hue="batch", aspect=15, height=0.5, palette="tab20")

    g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)

    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, value)

    g.figure.subplots_adjust(hspace=-.6)

    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    for ax in g.axes.flat:
        ax.axvline(x=df[value].median(), color='r', linestyle='-')
    
    plt.show()
    return g
    #g.savefig(os.path.join(qc_output_dir, value+'_beforeqc.png'), dpi = 300, bbox_inches = 'tight')


# In[77]:


#value = "pct_counts_mt"
#value = "n_genes"
#value = 'pct_counts_in_top_20_genes'
#value = "log1p_total_counts"
for value in qc_metrics:
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    g = sns.FacetGrid(df2, row="batch", hue="batch", aspect=15, height=0.5, palette="tab20")

    g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)

    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, value)

    g.figure.subplots_adjust(hspace=-.6)

    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    for ax in g.axes.flat:
        ax.axvline(x=df2[value].median(), color='r', linestyle='-')


    plt.show()
    g.savefig(os.path.join(qc_output_dir, value+'_beforeqc.png'), dpi = 300, bbox_inches = 'tight')


# In[84]:


# set initial cutoffs for pct_counts mt at 5, total_counts at 5k, and n_genes by 2500
# for better resolution of figures
pct_counts_mt_cutoff = 5
n_genes_cutoff = 2500
total_counts_cutoff = 5000


# In[142]:


adata_list_postqc = []
for adata in adata_list_nodoub:
    adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff]
    adata = adata[adata.obs.total_counts < total_counts_cutoff]
    adata = adata[adata.obs.pct_counts_mt < pct_counts_mt_cutoff]
    adata_list_postqc.append(adata)


# In[144]:


adata_list_postqc


# In[90]:


df = pd.concat([x.obs for x in adata_list_postqc])
df = df.sort_values('batch')
for value in qc_metrics:
    g = generate_qc_figs(df, value)
    g.savefig(os.path.join(qc_output_dir, value + '_firstroundqc.png'), dpi = 300, bbox_inches = 'tight')


# In[96]:


for sample in samples:
    print(sample+'...')
    for value in qc_metrics:
        sns.displot(df[df.batch == sample][value])
        plt.show()


# In[133]:


# set the second round cutoffs
n_genes_cutoff = 1500
total_counts_cutoff = 2000
pct_counts_mt_cutoff = 1


# In[134]:


for i,adata in enumerate(adata_list_postqc):
    adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff]
    adata = adata[adata.obs.total_counts < total_counts_cutoff]
    adata = adata[adata.obs.pct_counts_mt < pct_counts_mt_cutoff]
    print(adata.shape)
    adata_list_postqc[i] = adata


# In[145]:


adata_list_postqc


# In[146]:


df = pd.concat([x.obs for x in adata_list_postqc])
df = df.sort_values('batch')
for value in qc_metrics:
    g = generate_qc_figs(df, value)
    g.savefig(os.path.join(qc_output_dir, value + '_afterqc.png'), dpi = 300, bbox_inches = 'tight')


# In[149]:


output_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/scrna/data/velocyto_cellbender_processed/after_qc/indv_samples"


# In[148]:


for i,adata in enumerate(adata_list_postqc):
    adata.write(os.path.join(output_dir, samples[i]+"-filtered_combined.h5ad"))


# In[ ]:




