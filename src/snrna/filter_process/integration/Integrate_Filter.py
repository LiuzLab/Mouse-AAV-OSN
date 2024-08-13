import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import os
import rapids_singlecell as rsc # only use this if you have a gpu; if no GPU then use sc equivalent
import scvi

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")

import logging
#reduce logging amount
# Configure logging
logging.basicConfig(level=logging.WARNING)

input_dir = "data/velocyto_cellbender_processed/after_qc/indv_samples"
output_dir = "data/velocyto_cellbender_processed/after_qc/integrated"

# Process all screens
samples = [f'screen-{i}' for i in range(1, 5)]

adata_list = []
for sample in samples:
    #check which sample we're using
    print(f'loading {sample}' )
    # read in object
    adata = sc.read_h5ad(os.path.join(input_dir, sample + '-filtered_combined.h5ad'))
    adata_list.append(adata)

# set the batches
for i,adata in enumerate(adata_list):
    adata.obs['batch'] = samples[i]
    print(adata_list[i].obs['batch'])

#merge raw, unnormalized counts
adata = sc.concat(adata_list)

# Make obs and vars unique
adata.obs_names_make_unique()
adata.var_names_make_unique()

aav = ['AAV1', 'AAV7', 'AAVDJ8', 'AAVRH10']
sum(adata.var_names.isin(aav))

adata.write(os.path.join(output_dir, 'integrated_raw.h5ad'))

# perform basic filtering and QC
sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 1)

sum(adata.var_names.isin(aav))

adata

scvi.model.SCVI.setup_anndata(adata, layer = "spliced",
                             batch_key="batch",
                             #continuous_covariate_keys=['pct_counts_mt', 'total_counts']
                             )
model = scvi.model.SCVI(adata)
model.train(accelerator = 'gpu')

# catch doublets
solo = scvi.external.SOLO.from_scvi_model(model)
solo.train(accelerator = "gpu",
            )
df = solo.predict()
df['prediction'] = solo.predict(soft = False)

print(df)
print(df.groupby('prediction').count())

df['dif'] = df.doublet - df.singlet
sns.displot(df[df.prediction == 'doublet'], x = 'dif')
plt.axvline(1)
plt.show()

doublet_score_cutoff = 1
doublets = df[(df.prediction == 'doublet') & (df.dif > doublet_score_cutoff)]
doublets.shape

print(adata.shape)
adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
adata.obs
adata = adata[~adata.obs.doublet]
print(adata.shape)

qc_metrics = ['total_counts', 'log1p_total_counts',
              'n_genes', 'n_genes_by_counts', 
              'pct_counts_mt',
              'pct_counts_in_top_20_genes']
df2 = pd.concat([x.obs for x in adata])

for value in qc_metrics:
    sns.displot(df2[value])
    plt.show()

adata.write(os.path.join(output_dir, 'integrated-nodoublets.h5ad'))

adata.shape

bdata = adata.copy()

n_genes_cutoff = 1000
total_counts_cutoff = 1500
pct_counts_mt_cutoff = .6

bdata = bdata[bdata.obs.n_genes_by_counts < n_genes_cutoff]
bdata = bdata[bdata.obs.total_counts < total_counts_cutoff]
bdata = bdata[bdata.obs.pct_counts_mt < pct_counts_mt_cutoff]
print(bdata.shape)

bdata.write(os.path.join(output_dir, 'integrated-allfilter.h5ad'))
