# Make adjustments to annotations - JJ 7/29/2024

# load in the packages
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import os
import rapids_singlecell as rsc
import scvi # unneeded
import scvelo as scv
import palantir as pltr
from sklearn_ann.kneighbors.annoy import AnnoyTransformer  # currently unneeded
import math
import textwrap
import cellrank
from scipy import stats
sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300)
import warnings
import logging

# Ignore all warnings
warnings.filterwarnings('ignore')

# Set logging level to ERROR to suppress info and warning messages
logging.getLogger().setLevel(logging.ERROR)

# Filter out the specific DeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning, message="is_categorical_dtype is deprecated")
input_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/test/scrna/data/velocyto_cellbender_processed/analysis/"

output_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/test/scrna/results"
cr_output_dir = os.path.join(output_dir, "cellrank", "figs")
os.makedirs(cr_output_dir, exist_ok = True)
cr_output_table_dir = os.path.join(output_dir, 'cellrank', 'table')
os.makedirs(cr_output_table_dir, exist_ok = True)
# Process all screens
samples = [f'screen-{i}' for i in range(1, 5)]
# pull out gene list
aav = ['AAV1', 'AAV7', 'AAVDJ8', 'AAVRH10']
kwargs = dict(frameon=False, size=10, linewidth=1.5,
              legend_loc = 'on data', legend_fontsize = 6
               , cmap = "Spectral_r")

plot_kwargs = {'dpi':300, 'bbox_inches':'tight', 'facecolor':'white'}
## Load in object
file_name = "fullobject_palantir_dpt.h5ad"
ad = sc.read_h5ad(os.path.join(input_dir,file_name))

# add the original object into the raw attribute of the object (prior to HVG filtering)
ad_raw = sc.read_h5ad(os.path.join('/home/johnathan/projects/arrenkiel_scrnaseq/test/scrna/data/velocyto_cellbender_processed/after_qc/annotated', 'annotated_object_allgenes.h5ad'))
ad.shape
ad.raw = ad_raw.copy()

## Fix the celltype annotations
ad.obs['celltype_subs'] = ad.obs['celltype']
new_categories = ['MVC2', 'BowG2']
ad.obs['celltype_subs'] = ad.obs.celltype_subs.cat.add_categories(new_categories)
ad.obs.loc[ad.obs['celltype'] == 'isus', 'celltype'] = 'BowG'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype'] = 'MVC'

ad.obs.loc[ad.obs['leiden_3'] == '59', 'celltype_subs'] = 'BowG2'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype_subs'] = 'MVC2'
sc.pl.umap(ad, color = ['celltype', 'Muc2', 'Mybl1', 'Lgr5'], **kwargs)
with plt.rc_context({"figure.figsize": (12, 12), "figure.dpi": (300)}):
    sc.pl.umap(ad, color = ['celltype', 'Bpifb9b', 'Bpifb9a', 'Bpifb5',
                        'Trpm5', 'Ascl3',
                        'Obp2a', 'Obp2b'], **kwargs,
               ncols = 2)
result = ad.uns['dealeiden_2']
groups = result['names'].dtype.names
markers = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'logfoldchanges', 'pvals_adj']})
# Extract results for clusters 25 and 22
clusters_of_interest = ['40']

# Initialize an empty list to store the DataFrames for each cluster
selected_markers = []

# Loop through each cluster of interest and filter the DataFrame
for cluster in clusters_of_interest:
    cluster_df = markers.loc[:, [col for col in markers.columns if col.startswith(cluster + '_')]]
    selected_markers.append(cluster_df)

# Combine the DataFrames for the selected clusters
selected_markers_df = pd.concat(selected_markers, axis=1)

# Display the results
print(selected_markers_df.head(20))
gene_exp = ad.raw.to_adata()[:, list(markers.values())].to_df()

# 4. Calculate mean expression per cell type
mean_exp = gene_exp.groupby(ad.obs['celltype']).mean()

# 5. Hierarchical clustering to order cell types and genes
cell_type_order = mean_exp.index[sns.clustermap(mean_exp).dendrogram_row.reordered_ind]
gene_order = mean_exp.columns[sns.clustermap(mean_exp.T).dendrogram_row.reordered_ind]

# 6. Reorder the data
mean_exp_ordered = mean_exp.loc[cell_type_order, gene_order]

# 7. Calculate fraction of cells expressing each gene
fraction_expressed = (gene_exp > 0).groupby(ad.obs['celltype']).mean()
fraction_expressed_ordered = fraction_expressed.loc[cell_type_order, gene_order]

# 8. Plot improved dot plot
fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(mean_exp_ordered, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Mean expression'})

# Overlay dots
for i, cell_type in enumerate(mean_exp_ordered.index):
    for j, gene in enumerate(mean_exp_ordered.columns):
        fraction = fraction_expressed_ordered.loc[cell_type, gene]
        ax.scatter(j + 0.5, i + 0.5, s=fraction * 500, color='black', alpha=0.5)

ax.set_xlabel('Genes')
ax.set_ylabel('Cell Types')
ax.set_title('Improved Dot Plot')
plt.tight_layout()
plt.show()
sc.tl.rank_genes_groups(ad, 'celltype', method='wilcoxon')

# Get the top markers for each group
marker_dict = {}
for group in ad.obs['celltype'].unique():
    markers = sc.get.rank_genes_groups_df(ad, group=group, key='rank_genes_groups')
    top_markers = markers.head(30)['names'].tolist()  # Get top 5 markers
    marker_dict[group] = top_markers

# Convert to a DataFrame for easier viewing
marker_df = pd.DataFrame.from_dict(marker_dict, orient='index')
marker_df.columns = [f'Marker_{i+1}' for i in range(30)]

# Display the marker table
print(marker_df)
marker_df.loc['MV2']
marker_df.loc['MV2']
marker_df.loc['dSus']
marker_df.loc['iSus']
marker_df.loc['GBC']
ad.obs['celltype_subs'] = ad.obs['celltype']
new_categories = ['MVC2', 'BowG2']
ad.obs['celltype_subs'] = ad.obs.celltype_subs.cat.add_categories(new_categories)
ad.obs.loc[ad.obs['celltype'] == 'MV2', 'celltype'] = 'BowG'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype'] = 'MVC'

ad.obs.loc[ad.obs['leiden_3'] == '59', 'celltype_subs'] = 'BowG2'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype_subs'] = 'MVC2'
# repeat the process for dSus 
# the evidence is far too weak
ad.obs['celltype_subs'] = ad.obs.celltype_subs.cat.add_categories(['mSus2'])
ad.obs.loc[ad.obs['celltype'] == 'dSus', 'celltype'] = 'mSus'
ad.obs.loc[ad.obs['celltype_subs'] == 'dSus', 'celltype_subs'] = 'mSus2'
# drop the unused categories
ad.obs['celltype'] = ad.obs['celltype'].cat.remove_unused_categories()
ad.obs['celltype_subs'] = ad.obs['celltype_subs'].cat.remove_unused_categories()
# save the corrected cell types and new "subcluster" annotations
ad.write_h5ad(os.path.join(input_dir, file_name))