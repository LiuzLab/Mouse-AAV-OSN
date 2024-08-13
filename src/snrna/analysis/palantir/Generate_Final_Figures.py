#!/usr/bin/env python
# coding: utf-8

# # Notebook for Generating Final Figures

# ## Load in the Packages

# In[1]:


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
import random


# In[2]:


# Set the random seed
np.random.seed(42)
random.seed(42)


# In[3]:


sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300)


# In[4]:


import warnings
import logging

# Ignore all warnings
warnings.filterwarnings('ignore')

# Set logging level to ERROR to suppress info and warning messages
logging.getLogger().setLevel(logging.ERROR)

# Filter out the specific DeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning, message="is_categorical_dtype is deprecated")


# In[5]:


input_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/test/scrna/data/velocyto_cellbender_processed/analysis/"

output_dir = "/home/johnathan/projects/arrenkiel_scrnaseq/test/scrna/results"
cr_output_dir = os.path.join(output_dir, "cellrank", "figs")
os.makedirs(cr_output_dir, exist_ok = True)
cr_output_table_dir = os.path.join(output_dir, 'cellrank', 'table')
os.makedirs(cr_output_table_dir, exist_ok = True)


# In[6]:


# Process all screens
samples = [f'screen-{i}' for i in range(1, 5)]
# pull out gene list
aav = ['AAV1', 'AAV7', 'AAVDJ8', 'AAVRH10']


# In[7]:


kwargs = dict(frameon=False, size=10, linewidth=1.5,
              legend_loc = 'on data', legend_fontsize = 6
               , cmap = "Spectral_r")


# In[8]:


plot_kwargs = {'dpi':300, 'bbox_inches':'tight', 'facecolor':'white'}


# In[9]:


def chi_squared_test(a_probs, b_probs):
    # Ensure both arrays have the same length
    assert len(a_probs) == len(b_probs), "HBC and GBC probability arrays must have the same length"
    
    # Calculate chi-squared statistic
    chi2_stat = np.sum((a_probs - b_probs)**2 / b_probs)
    
    # Degrees of freedom
    df = len(a_probs) - 1
    
    # Calculate p-value
    p_value = 1 - stats.chi2.cdf(chi2_stat, df)
    
    return chi2_stat, df, p_value


# ## Load in object

# In[10]:


file_name = "fullobject_palantir_dpt.h5ad"


# In[11]:


ad = sc.read_h5ad(os.path.join(input_dir,file_name))


# In[10]:


# run some checks for the identity of this cluster
print(ad.obs.celltype.value_counts())
print(ad.obs.celltype_subs.value_counts())
print(ad.shape)
print(ad.raw.to_adata().shape)


# ## Load in genes of interest

# In[12]:


genes_of_interest = ['Ano2',
'Olfm1',
'Chga',
'Kcnk10',
'Cnga2',
'Sox11',
'Dpysl5', 'Dpysl2', 'Hdac',
'Lhx2',
'Ncam2']
olf_genes = ['Or6b2', 'Or6b2b', 'Or6b3', 'Or10j27', 'Or12k8']
osn_maturation_genes = ["Ascl1", "Neurog1", "Neurog2", "Tbx21", "Emx2", "Sox2", "Sox11", "Ngn1", "Lhx2", "Lhx9", "Ebf1", "Ebf2", "Mash1", "Otx2", "Gap43", "Neurod1", "Dlx5", "Dlx6", "Pax6", "Zic1", "Zic2", "Foxg1", "Irx3", "Bdnf", "Gdnf"
]
axon_targeting_genes= ["Reln", "Shh", "Wnt3a", "Dcx", "Nrg1", "Nrp1", "EphA4", "Bdnf", "Slit1", "Robo1", "Robo2", "PSA-NCAM", "Prok2", "Tbr2", "Sema3A", "Nrp2", "Gad1", "EphA5", "EphA6", "Efna5", "Ntn1", "Dcc", "Cdh1", "Cdh2", "Pcdh", "L1CAM",] 
axon_targeting_genes2= ["Tnr", "Cntn2", "Ctnnb1", "Sema7A", "Robo3", "Itga3", "Itgb1", "PlxnA1", "Robo1", "Robo2", "PSA-NCAM", "Cdk5", "Sema3C", "Dscam", "Lhx2", "Slitrk1", "Slitrk2", "Gfra1", "Lhx9", "Gli1", "Gli2", "Lrig1", "Erbb4", "Nrn1", "Slit3", "Bcl11b", "FOXP2", "Nrp2", "EphB2", "Atf5", "Axin2", "Nrp1", "Tnc", "Zic2", "Or2p2", "Or2ab1", "Or9e1", "Or2f1"
]
other_genes_of_interest = ["Nrp2", "PlxnA4", "NrCAM", "Fzd3", "Celsr3", "Dscam", "Slit2", "Sox11", "Gfra1", "Robo4", "Arhgap21", "Cux1", "Cux2", "Spinophilin", "Kalirin", "Neurod2", "Pcdh10", "Map2", "Shank3", "Dlg4", "Bdnf", "Camk2a", "Arpc3", "Mef2c", "Notch1", "Omp", "Gnai2", "Gucy1b2", "Adcy3", "Pax6", "Cyp2g1", "Slc8a1"
]


# In[13]:


ad


# ## Add original raw object back to AnnData object - Skip if after 7/31/2024

# In[24]:


ad_raw = sc.read_h5ad(os.path.join('/home/johnathan/projects/arrenkiel_scrnaseq/test/scrna/data/velocyto_cellbender_processed/after_qc/annotated', 'annotated_object_allgenes.h5ad'))


# In[25]:


ad.shape


# In[26]:


ad.raw = ad_raw.copy()


# In[27]:


ad.raw.to_adata()


# In[30]:


markers.values()


# ## Fix the celltype annotations - Skip if running after 7/31/2024

# In[48]:


ad.obs['celltype_subs'] = ad.obs['celltype']


# In[49]:


new_categories = ['MVC2', 'BowG2']
ad.obs['celltype_subs'] = ad.obs.celltype_subs.cat.add_categories(new_categories)


# In[50]:


ad.obs.loc[ad.obs['celltype'] == 'isus', 'celltype'] = 'BowG'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype'] = 'MVC'

ad.obs.loc[ad.obs['leiden_3'] == '59', 'celltype_subs'] = 'BowG2'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype_subs'] = 'MVC2'


# In[84]:


sc.pl.umap(ad, color = ['celltype', 'Muc2', 'Mybl1', 'Lgr5'], **kwargs)


# In[ ]:


with plt.rc_context({"figure.figsize": (12, 12), "figure.dpi": (300)}):
    sc.pl.umap(ad, color = ['celltype', 'Bpifb9b', 'Bpifb9a', 'Bpifb5',
                        'Trpm5', 'Ascl3',
                        'Obp2a', 'Obp2b'], **kwargs,
               ncols = 2)


# In[46]:


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


# In[ ]:


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


# In[38]:


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


# In[39]:


marker_df.loc['MV2']


# In[ ]:


marker_df.loc['MV2']


# In[81]:


marker_df.loc['dSus']


# In[129]:


marker_df.loc['iSus']


# In[128]:


marker_df.loc['GBC']


# In[ ]:


ad.obs['celltype_subs'] = ad.obs['celltype']


# In[ ]:


new_categories = ['MVC2', 'BowG2']
ad.obs['celltype_subs'] = ad.obs.celltype_subs.cat.add_categories(new_categories)


# In[ ]:


ad.obs.loc[ad.obs['celltype'] == 'MV2', 'celltype'] = 'BowG'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype'] = 'MVC'

ad.obs.loc[ad.obs['leiden_3'] == '59', 'celltype_subs'] = 'BowG2'
ad.obs.loc[ad.obs['leiden_3'] == '60', 'celltype_subs'] = 'MVC2'


# In[89]:


# repeat the process for dSus 
# the evidence is far too weak
ad.obs['celltype_subs'] = ad.obs.celltype_subs.cat.add_categories(['mSus2'])
ad.obs.loc[ad.obs['celltype'] == 'dSus', 'celltype'] = 'mSus'
ad.obs.loc[ad.obs['celltype_subs'] == 'dSus', 'celltype_subs'] = 'mSus2'


# In[90]:


# drop the unused categories
ad.obs['celltype'] = ad.obs['celltype'].cat.remove_unused_categories()
ad.obs['celltype_subs'] = ad.obs['celltype_subs'].cat.remove_unused_categories()


# In[91]:


# save the corrected cell types and new "subcluster" annotations
ad.write_h5ad(os.path.join(input_dir, file_name))


# # Generate the final UMAP, Dotplot, and Barplot figures

# In[13]:


kwargs['legend_fontsize'] = 4
final_umap = sc.pl.umap(ad, color = 'celltype',title = '', **kwargs, return_fig = True)
final_umap.savefig('scrna/results/final_umap_withiSus.png', **plot_kwargs)

plt.show()
plt.close()


# In[200]:


kwargs['legend_fontsize'] = 4
final_umap = sc.pl.umap(ad, color = 'celltype_subs',title = 'Annotated UMAP', **kwargs, return_fig = True)
final_umap.savefig('scrna/results/final_umap_withiSus_subclusters.png', **plot_kwargs)

plt.show()
plt.close()


# In[276]:


import scanpy as sc
import matplotlib.pyplot as plt

# Generate the UMAP plot
umap = sc.pl.umap(ad, color='celltype',
                  frameon=False, 
                  legend_fontsize=6, return_fig=True, show=True,
                  title='Annotated UMAP')
ax = umap.axes[0]

# Remove the original legend
ax.get_legend().remove()

# Extract handles and labels from the original legend
handles, labels = ax.get_legend_handles_labels()

# Create a dictionary mapping abbreviated labels to full names
label_mapping = {
    'ACC': 'Airway Ciliated Cell (ACC)',
    'BowG': 'Bowman\'s Glands (BowG)',
    'Gob': 'Goblet Cells (Gob)',
    'MVC': 'Microvillar Cells (MVC)',
    'Mes': 'Mesenchymal Cells (Mes)',
    'iOSN': 'Immature olfactory sensory neurons (iOSN)',
    'mOEC': 'Mature olfactory ensheathing cells (mOEC)',
    'mSus' : 'Mature sustentacular cells (mSus)',
    'VE' : 'Vascular endothelial Cells (VE)',
    'SCC' : 'Solitary chemosensory Cell (SCC)',
    'GBC' : 'Globose basal cell (GBC)',
    'HBC' : 'Horizontal basal cell (HBC)',
    'mOSN': 'Mature olfactory sensory neurons (mOSN)',
    'Immune' : 'Immune cells (Immune)',
    'RE' : 'Respiratory Epithelium (RE)',
    'UnSec' : 'Undefined secretory cells (UnSec)',
    'LE' : 'Lymphatic epithelium (LE)',
    'iSus' : 'Immature sustentacular cells (iSus)'
}

# Create custom labels using the dictionary
custom_labels = [label_mapping.get(label, label) for label in labels]

# Create a new legend with the custom labels and desired number of columns
ax.legend(handles, custom_labels, ncol=2, title='Cell Types', 
          bbox_to_anchor=(0.5, -2), loc='center',
          fontsize=22, title_fontsize=32)

# Set other legend properties
ax.get_legend().set_frame_on(True)

# Save the plot with the new legend
umap.savefig(os.path.join('scrna/results/', 'annotated_umap_custom_legend.png'), dpi=300,
             bbox_inches='tight')

plt.show()


# In[202]:


# Group by 'sample' and 'celltype' and count occurrences
celltype_counts = ad.obs.groupby(['batch', 'celltype']).size()

# Optionally, unstack to get a more readable format
celltype_counts_unstacked = celltype_counts.unstack(fill_value=0)

print(celltype_counts_unstacked)


# In[203]:


# Ensure the categories in the order of the columns of `celltype_counts_unstacked`
new_order = celltype_counts_unstacked.columns.tolist()
print(new_order)


# In[204]:


ad.obs.celltype.cat.categories


# In[205]:


ad.uns['celltype_colors']


# In[206]:


# Perform Chi-Squared GOF Test
chi2, p, dof, expected = stats.chi2_contingency(celltype_counts_unstacked.T)
print(f'Chi-Squared Test Statistic: {chi2}')
print(f'p-value: {p}')
print(f'Degrees of Freedom: {dof}')
print('Expected Counts:')
print(expected)


# In[207]:


# Melt the DataFrame for seaborn
celltype_counts_melted = celltype_counts_unstacked.reset_index().melt(id_vars='batch', value_name='count', var_name='celltype')

# Reorder the categories in the order of value counts
new_order = celltype_counts_unstacked.columns.tolist()

# Reorder the categories in the AnnData object
ad.obs['celltype'] = ad.obs['celltype'].cat.reorder_categories(new_order, ordered=True)

# Create a color palette that matches the reordered cell type categories
color_pal = {col: ad.uns['celltype_colors'][i] for i, col in enumerate(ad.obs['celltype'].cat.categories)}


# Create a bar plot
plt.figure(figsize=(14, 8))
sns.barplot(data=celltype_counts_melted, x='batch', y='count', hue='celltype', palette = color_pal)

# Add p-value and chi-squared test statistic to the plot
chi2_text = f'Chi-Squared: {chi2:.2f}\np-value: {p:.2e}\ndf: {dof}'
plt.annotate(chi2_text, xy=(1.1, 1), xycoords='axes fraction', fontsize=16,
             xytext=(-70, 10), textcoords='offset points',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

plt.xticks(rotation=45)
plt.title('Cell Type Counts per Sample', fontsize = 32)
plt.xlabel('Sample', fontsize = 22)
plt.ylabel('Count', fontsize = 22)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
plt.legend(title='Cell Type', title_fontsize = 22, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize = 16)
plt.tight_layout()
plt.savefig('scrna/results/celltype_proportion_barplot_bysample.png', **plot_kwargs)
plt.show()


# In[208]:


#6/20/2024 -> Remake this Figure with NO grid lines
# Calculate the value counts
cell_type_counts = ad.obs['celltype'].value_counts()

# Get the colors for the cell types in the order of cell_type_counts
colors = [color_pal[cell_type] for cell_type in cell_type_counts.index]

# Create the bar plot
plt.figure(figsize=(10, 6))
bar_plot = cell_type_counts.plot(kind='bar', color=colors, width=0.7)

# Annotate the bars with the actual numbers
for p in bar_plot.patches:
    bar_plot.annotate(
        str(p.get_height()), 
        (p.get_x() + p.get_width() / 2., p.get_height()), 
        ha='center', va='center', 
        xytext=(2, 17), 
        textcoords='offset points',
        rotation = 45,
        fontsize = 14
    )

plt.ylim(0, max(ad.obs.celltype.value_counts())+500)

# Remove gridlines
plt.grid(False)

# Customize the plot
plt.title('Cell Type Distribution', fontsize = 22)
plt.xlabel('Cell Type', fontsize = 22)
plt.ylabel('Number of Cells', fontsize = 22)
plt.xticks(rotation=45, fontsize = 12)  # Rotate x-axis labels if they are too long
plt.savefig(os.path.join('scrna/results/celltype_proportion_barplot.png'),**plot_kwargs)
plt.show()


# In[209]:


# Create a bar plot
plt.figure(figsize=(14, 8))
ax = sns.barplot(data=celltype_counts_melted, x='batch', y='count', hue='celltype', palette=color_pal)

# Add numbers above the bars at 45-degree angle
for p in ax.patches:
    ax.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 9), 
                textcoords = 'offset points', 
                rotation=45)

plt.xticks(rotation=90)
plt.title('Cell Type Counts per Sample')
plt.xlabel('Sample')
plt.ylabel('Count')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()


# In[210]:


from scipy import stats

# Perform Chi-Squared GOF Test
chi2, p, dof, expected = stats.chi2_contingency(celltype_counts_unstacked.T)
print(f'Chi-Squared Test Statistic: {chi2}')
print(f'p-value: {p}')
print(f'Degrees of Freedom: {dof}')
print('Expected Counts:')
print(expected)


# In[211]:


# Set up the figure for 2x2 subplots
batches = ad.obs['batch'].unique()
n_batches = len(batches)
n_rows = 2
n_cols = 2

fig, axs = plt.subplots(n_rows, n_cols, figsize=(14, 14))
axs = axs.flatten()  # Flatten the 2D array of axes to make indexing easier

# Generate scatter plots for each batch
for i, batch in enumerate(batches):
    if i < n_rows * n_cols:  # Ensure we don't exceed the number of subplots
        batch_data = ad[ad.obs['batch'] == batch]
        scatter = axs[i].scatter(
            batch_data.obs['n_genes'],
            batch_data.obs['total_counts'],
            c=batch_data.obs['pct_counts_mt'],
            cmap='viridis'
        )
        axs[i].set_title(f'Batch: {batch}')
        axs[i].set_ylabel('Total Counts')
        axs[i].set_xlabel('Number of Genes by Counts')
        
        # Remove gridlines
        axs[i].grid(False)

# Remove any unused subplots
for j in range(i+1, n_rows * n_cols):
    fig.delaxes(axs[j])

# Adjust spacing between subplots
plt.tight_layout()
plt.subplots_adjust(hspace=0.3, wspace=0.3, bottom=0.15)

# Add a colorbar below the subplots
cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.02])  # [left, bottom, width, height]
cbar = fig.colorbar(scatter, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Percentage of Mitochondrial Counts (%)')

plt.savefig('scrna/results/final_qc_figs.png', **plot_kwargs)

plt.show()


# In[ ]:


markers = {
    'ACC' : 'Cdhr3',
    'BowG' : 'Muc5b',
    'Gob' : 'Agr2', # or LTF, Muc5ac
    'MVC' : 'Ascl3',
    'Mes' : 'Carmn',
    'iOSN' : 'Sox11',
    'mOEC' : 'Alx3',
    'mOSN' : 'Ano2', #Olfm1 works too
    'mSus' : 'Muc2',
    'VE' : 'Flt1',
    'SCC' : 'Trpm5',
    'GBC' : 'Hes6',
    'HBC' : 'Krt5',
    'Immune' : 'Cd74',
    'LE' : 'Ccl21a',
    'RE' : 'Bpifb4',
    'UnSec' : 'Tmsb4x',
    #'BowG2' : 'Bpifb9b', #leiden 40
    #'MVC2' : 'Obp2a', #leiden 39
    #mSus2 : "???", 
    'iSus' : 'Bpifa1'
}


# In[ ]:


sc.tl.dendrogram(ad, groupby = 'celltype')


# In[ ]:


marker_dot_plot = sc.pl.dotplot(ad, groupby = 'celltype', var_names = markers, log = True, standard_scale = 'group',
              cmap = 'Blues', dendrogram= True, return_fig = True)
marker_dot_plot.savefig('scrna/results/final_marker_dotplot.png', **plot_kwargs)
plt.show()
plt.close()


# In[ ]:


sc.pl.embedding(ad, basis = 'draw_graph_fa', color = ['celltype', 'palantir_pseudotime', 'dpt_pseudotime'],
                **kwargs)


# # Start Palantir Analysis here

# ### Work with full object here

# In[42]:


ad_og = ad.copy()
ad = ad[ad.obs['celltype'] != 'ACC']
print(ad.obs.celltype.value_counts())


# In[46]:


# Step 1: Extract the original color palette from the anndata object
original_celltype_colors = ad_og.uns['celltype_colors']
original_cell_types = ad_og.obs['celltype'].cat.categories

# Step 2: Create a dictionary mapping cell types to colors
color_palette = {cell_type: color for cell_type, color in zip(original_cell_types, original_celltype_colors)}

# Step 3: Subset the anndata object
#already subsetted

# Step 4: Ensure the subsetted anndata object uses the same color palette
# Extract cell types present in the subset
subset_cell_types = ad.obs['celltype'].cat.categories

# Create a custom palette for the subset based on the original color mapping
subset_custom_palette = [color_palette[cell_type] for cell_type in subset_cell_types]

# Add the custom palette into original object
ad.uns['celltype_colors'] = subset_custom_palette


# In[43]:


rsc.pp.neighbors(ad, n_neighbors= 30, use_rep = 'X_scvi')
rsc.tl.umap(ad)
sc.pl.umap(ad, color = 'celltype', **kwargs)
rsc.tl.mde(ad, n_neighbors=30, n_pcs = 30)
sc.pl.embedding(ad, basis = 'mde', color = 'celltype')


# In[44]:


# Run diffusion maps
pca_rd = pltr.utils.run_pca(ad)


# In[45]:


# Variance explained by each PC
variance_explained = ad.uns['pca']['variance_ratio']

# Calculate cumulative variance explained
cumulative_variance_explained = np.cumsum(variance_explained)

# Find the number of PCs that explain 85% of the variance
num_pcs_85 = np.argmax(cumulative_variance_explained >= 0.85) + 1

print(f"Number of PCs explaining 85% of the variance: {num_pcs_85}")

# Plot the cumulative variance explained
plt.figure(figsize=(8, 5))
plt.plot(cumulative_variance_explained, marker='o')
plt.axhline(y=0.85, color='r', linestyle='--')
plt.axvline(x=num_pcs_85 - 1, color='r', linestyle='--')
plt.xlabel('Number of Principal Components')
plt.ylabel('Cumulative Variance Explained')
plt.title('Cumulative Variance Explained by Principal Components')
plt.show()


# In[68]:


dff_map = pltr.utils.run_diffusion_maps(ad, n_components = 30)


# In[69]:


ms_data = pltr.utils.determine_multiscale_space(ad)


# In[70]:


imputed_X = pltr.utils.run_magic_imputation(ad, n_jobs = -1)


# In[71]:


pltr.plot.plot_diffusion_components(ad, cmap = 'Spectral_r')
plt.show()
plt.close()


# In[72]:


terminal_states = pd.Series(
    [#"mOSN1", 
     #'mOSN2',
     #'mOSN3',
     'mOSN',
     #'mOSN5',
     #'mOSN6', 
     "BowG", 
     #"Gob", 
     "mSus", 
     "mOEC" ,
     "MVC",
 #    "ACC",
     ],
    index=[
      #  "AACAGGGGTCCCGTGA",  #mOSN1
      #  'AATTCCTTCCTGGGTG', #mOSN2
        #'CCGCGTTAGTCAGAGC', #mOSN3
        'CAGCACGAGAGATGGT', #mOSN4
        #'CACCAAGAGCTTAGGA', #mOSN5
        #'ACCATTTCTCTCCCTA', #mOSN6
        "GCTGAATCGTAGTATC",  # Bowman's Gland
         # "ATGAGTCTCCACCGGG", #gob
           "GAAGTAAAGGAGCTAG", #mature Sus
           "CAGGCCACTTTTAAAA",  #mOEC
          "AGCCAGCAGAAGGGAT", #MVC
 #          "CGCATGGAGAAGGATG", #ACC
    ]
)


# In[73]:


pltr.plot.highlight_cells_on_umap(ad, terminal_states)
plt.show()
plt.close()


# In[63]:


early_cell = pltr.utils.early_cell(ad, 'HBC')


# In[74]:


print(early_cell)


# In[75]:


#CGACTTGTCCATACTT, GGCACTTGTGACCCGC = GBC
#ACAGCGGGTATCAGGG, CCCTCCCTACTCAGTT = HBC
#GGACACTTCCAGCCTT # this for iSus/Trans?
pr_res = pltr.core.run_palantir(
    ad, early_cell,  #"CCCTCCCTACTCAGTT", 
    terminal_states=terminal_states, 
    #max_iterations = 50,
    knn = 15,
    use_early_cell_as_start= True,
)


# In[76]:


pseudotime_ptr_results = pltr.plot.plot_palantir_results(ad, #embedding_basis= 'X_mde',
                                                         s = 5, cmap = 'Spectral_r')
plt.show()
plt.close()


# In[77]:


hbc_index = np.where(ad.obs['celltype'] == 'HBC')[0]
hbc_barcodes = ad.obs_names[hbc_index]


# In[78]:


gbc_index = np.where(ad.obs['celltype'] == 'GBC')[0]
gbc_barcodes = ad.obs_names[gbc_index]


# In[79]:


iSus_index = np.where(ad.obs['celltype'] == 'iSus')[0]
iSus_barcodes = ad.obs_names[iSus_index]


# In[80]:


fate_probs = ad.obsm['palantir_fate_probabilities']


# In[81]:


# Select HBC fate probabilities
hbc_fate_probs = fate_probs.iloc[hbc_index]

# Now we can melt the DataFrame
hbc_fate_probs_melted = hbc_fate_probs.reset_index().melt(id_vars='barcode', var_name='Cell Type', value_name='Probability')

# Select GBC fate probabilities
gbc_fate_probs = fate_probs.iloc[gbc_index]

# Now we can melt the DataFrame
gbc_fate_probs_melted = gbc_fate_probs.reset_index().melt(id_vars='barcode', var_name='Cell Type', value_name='Probability')

hbc_fate_probs_melted['group'] = 'HBC'
gbc_fate_probs_melted['group'] = 'GBC'


# In[82]:


iSus_fate_probs = fate_probs.iloc[iSus_index]
iSus_fate_probs_melted = iSus_fate_probs.reset_index().melt(id_vars='barcode', var_name='Cell Type', value_name='Probability')
iSus_fate_probs_melted['group'] = 'iSus'


# In[83]:


# Set up the plot
plt.figure(figsize=(16, 8))  # Increased figure size for better visibility

# Create the violin plot with modifications
sns.violinplot(x='Cell Type', y='Probability', data=iSus_fate_probs_melted,
               inner='box',  # This will show the box plot inside the violin
               cut=0,  # This extends the violin plot to the full range of the data
               scale='width',  # This makes all violins the same width
               width=0.8,  # Adjust this to make violins wider or narrower
               palette='Set1',
               linewidth=2)  # Change color palette (you can try others like 'Set2', 'deep', etc.) (husl is good)

# Add individual points
sns.swarmplot(x='Cell Type', y='Probability', data=iSus_fate_probs_melted,
              color='black', edgecolor='black', size=3, alpha=0.5)

# Customize the plot
plt.title('Palantir Fate Probabilities for iSus', fontsize=16)
plt.xlabel('Cell Type', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.ylim(-0.05, 1.05)  # Adjust y-axis limits if needed

# Improve layout
plt.tight_layout()
#plt.savefig(os.path.join(pltir_output_dir, 'trans_average_fate_probabilities.png'), **plot_kwargs)

# Show the plot
plt.show()
plt.close()


# In[84]:


# Combine the data
combined_data = pd.concat([hbc_fate_probs_melted, gbc_fate_probs_melted])

# Set up the plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

color_palette = sns.color_palette(palette='Set1')

# Plot HBC data
sns.violinplot(x='Cell Type', y='Probability', data=hbc_fate_probs_melted, ax=ax1, 
               palette=color_palette, inner='box', cut=0, scale='width', width=0.8)
sns.swarmplot(x='Cell Type', y='Probability', data=hbc_fate_probs_melted, ax=ax1,
              color='black', edgecolor='white', size=3, alpha=0.5)
ax1.set_title('Palantir Fate Probabilities for HBCs', fontsize=16)
ax1.set_ylim(-0.05, 1.05)

# Plot GBC data
sns.violinplot(x='Cell Type', y='Probability', data=gbc_fate_probs_melted, ax=ax2, 
               palette=color_palette, inner='box', cut=0, scale='width', width=0.8)
sns.swarmplot(x='Cell Type', y='Probability', data=gbc_fate_probs_melted, ax=ax2,
              color='black', edgecolor='white', size=3, alpha=0.5)
ax2.set_title('Palantir Fate Probabilities for GBCs', fontsize=16)
ax2.set_ylim(-0.05, 1.05)

# Customize both plots
for ax in [ax1, ax2]:
    ax.set_xlabel('Cell Type', fontsize=14)
    ax.set_ylabel('Probability', fontsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
plt.tight_layout()
#plt.savefig(os.path.join(pltir_output_dir, 'hbc_gbc_sidebyside_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[85]:


# Combine the data
combined_data = pd.concat([hbc_fate_probs_melted, gbc_fate_probs_melted])

# Set up the plot
plt.figure(figsize=(20, 10))
color_palette = sns.color_palette("Set1")

# Create a custom order for cell types and groups
cell_types = hbc_fate_probs_melted['Cell Type'].unique()
groups = ['HBC', 'GBC']

# Create the violin plot
ax = sns.violinplot(x='Cell Type', y='Probability', hue='group', data=combined_data,
                    split=True, inner="box", palette=color_palette[:2],
                    cut=0, scale='width', width=0.8)

# Add swarm plots with manual offset
for i, group in enumerate(groups):
    sns.swarmplot(x='Cell Type', y='Probability', data=combined_data[combined_data['group'] == group],
                  color='black', edgecolor='white', size=3, alpha=0.5,
                  ax=ax, dodge=True, 
                  zorder=1)  # Ensure points are on top

# Adjust swarmplot positions
for i, col in enumerate(ax.collections[len(cell_types)*2:]):  # Skip box plot artists
    if i % 2 == 0:  # HBC points
        col.set_offsets(col.get_offsets() - np.array([0.2, 0]))
    else:  # GBC points
        col.set_offsets(col.get_offsets() + np.array([0.2, 0]))

# Customize the plot
plt.title('Palantir Fate Probabilities for HBCs and GBCs', fontsize=32)
plt.xlabel('Cell Type', fontsize=22)
plt.ylabel('Probability', fontsize=22)
plt.ylim(-0.05, 1.2)
plt.xticks(rotation=45, ha='right', fontsize = 22)
y_ticks = np.arange(0, 1.2, 0.2)  # Create an array of y-ticks from 0 to 1 with a step of 0.1
plt.yticks(y_ticks, fontsize = 22)
plt.legend(fontsize=22, bbox_to_anchor=(1.06, 1), loc='upper right', edgecolor = 'black')
plt.grid(False)

# Function to add significance stars
def get_stars(p):
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    else:
        return "ns"

# Perform t-tests and add significance bars with stars
for i, cell_type in enumerate(cell_types):
    hbc_probs = hbc_fate_probs_melted[hbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    gbc_probs = gbc_fate_probs_melted[gbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    
    t_stat, p_value = stats.ttest_ind(hbc_probs, gbc_probs)
    
    y_max = max(hbc_probs.max(), gbc_probs.max())
    y_bar = y_max + 0.05
    
    plt.plot([i-0.2, i+0.2], [y_bar, y_bar], color='black', linewidth=3)
    stars = get_stars(p_value)
    plt.text(i, y_bar + 0.02, stars, ha='center', va='bottom', fontsize=22) # significant stars
    #plt.text(i, y_bar + 0.06, f'p={p_value:.3f}', ha='center', va='bottom', fontsize=22) # p-values

# Calculate average probabilities for each fate
avg_probs_hbc = hbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()
avg_probs_gbc = gbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()

# Perform chi-squared test between HBC and GBC
chi2, df, p_value = chi_squared_test(avg_probs_hbc.values, avg_probs_gbc.values)

# Add text box with chi-squared test results
textstr = '\n'.join((
    f'χ²={chi2:.2f}',
    f'df={df}',
    f'p={p_value:.4e}'
))

# Position the text box in figure coords
props = dict(boxstyle='round', facecolor='white', alpha=1, edgecolor = 'black')
plt.text(0.95, 0.95, textstr, transform=plt.gcf().transFigure, fontsize=22,
         verticalalignment='top', bbox=props)

# Adjust layout and show the plot
#plt.tight_layout()
plt.subplots_adjust(top=0.85)  # Make room for the text box
#plt.savefig(os.path.join(output_dir, 'final_minimalsubset_hbcvgbc_compared_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[86]:


# Combine the data
combined_data = pd.concat([hbc_fate_probs_melted, iSus_fate_probs_melted])

# Set up the plot
plt.figure(figsize=(20, 10))
color_palette = sns.color_palette("Set1")

# Create a custom order for cell types and groups
cell_types = hbc_fate_probs_melted['Cell Type'].unique()
groups = ['HBC', 'iSus']

# Create the violin plot
ax = sns.violinplot(x='Cell Type', y='Probability', hue='group', data=combined_data,
                    split=True, inner="box", palette=color_palette[:2],
                    cut=0, scale='width', width=0.8)

# Add swarm plots with manual offset
for i, group in enumerate(groups):
    sns.swarmplot(x='Cell Type', y='Probability', data=combined_data[combined_data['group'] == group],
                  color='black', edgecolor='white', size=3, alpha=0.5,
                  ax=ax, dodge=True, 
                  zorder=1)  # Ensure points are on top

# Adjust swarmplot positions
for i, col in enumerate(ax.collections[len(cell_types)*2:]):  # Skip box plot artists
    if i % 2 == 0:  # HBC points
        col.set_offsets(col.get_offsets() - np.array([0.2, 0]))
    else:  # GBC points
        col.set_offsets(col.get_offsets() + np.array([0.2, 0]))

# Customize the plot
plt.title('Palantir Fate Probabilities for HBCs and iSus', fontsize=32)
plt.xlabel('Cell Type', fontsize=22)
plt.ylabel('Probability', fontsize=22)
plt.ylim(-0.05, 1.2)
plt.xticks(rotation=45, ha='right', fontsize = 22)
plt.yticks(y_ticks, fontsize = 22)
plt.legend(fontsize=22, bbox_to_anchor=(1.06, 1), loc='upper right', edgecolor = 'black')
plt.grid(False)

# Function to add significance stars
def get_stars(p):
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    else:
        return "ns"

# Perform t-tests and add significance bars with stars
for i, cell_type in enumerate(cell_types):
    hbc_probs = hbc_fate_probs_melted[hbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    iSus_probs = iSus_fate_probs_melted[iSus_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    
    t_stat, p_value = stats.ttest_ind(hbc_probs, iSus_probs)
    
    y_max = max(hbc_probs.max(), iSus_probs.max())
    y_bar = y_max + 0.05
    
    plt.plot([i-0.2, i+0.2], [y_bar, y_bar], color='black', linewidth=3)
    stars = get_stars(p_value)
    plt.text(i, y_bar + 0.02, stars, ha='center', va='bottom', fontsize=22)
    #plt.text(i, y_bar + 0.06, f'p={p_value:.3f}', ha='center', va='bottom', fontsize=22)

# Calculate average probabilities for each fate
avg_probs_hbc = hbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()
avg_probs_iSus = iSus_fate_probs_melted.groupby('Cell Type')['Probability'].mean()

# Perform chi-squared test between HBC and GBC
chi2, df, p_value = chi_squared_test(avg_probs_hbc.values, avg_probs_iSus.values)

# Add text box with chi-squared test results
textstr = '\n'.join((
    f'χ²={chi2:.2f}',
    f'df={df}',
    f'p={p_value:.4e}'
))

# Position the text box in figure coords
props = dict(boxstyle='round', facecolor='white', alpha=1, edgecolor = 'black')
plt.text(0.90, 0.95, textstr, transform=plt.gcf().transFigure, fontsize=22,
         verticalalignment='top', bbox=props)

# Adjust layout and show the plot
plt.tight_layout()
plt.subplots_adjust(top=0.85)  # Make room for the text box
#plt.savefig(os.path.join(output_dir, 'final_minimalsubset_hbcvisus_compared_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[87]:


# Combine the data
combined_data = pd.concat([gbc_fate_probs_melted, iSus_fate_probs_melted])

# Set up the plot
plt.figure(figsize=(20, 10))
color_palette = sns.color_palette("Set1")

# Create a custom order for cell types and groups
cell_types = gbc_fate_probs_melted['Cell Type'].unique()
groups = ['GBC', 'iSus']

# Create the violin plot
ax = sns.violinplot(x='Cell Type', y='Probability', hue='group', data=combined_data,
                    split=True, inner="box", palette=color_palette[:2],
                    cut=0, scale='width', width=0.8)

# Add swarm plots with manual offset
for i, group in enumerate(groups):
    sns.swarmplot(x='Cell Type', y='Probability', data=combined_data[combined_data['group'] == group],
                  color='black', edgecolor='white', size=3, alpha=0.5,
                  ax=ax, dodge=True, 
                  zorder=1)  # Ensure points are on top

# Adjust swarmplot positions
for i, col in enumerate(ax.collections[len(cell_types)*2:]):  # Skip box plot artists
    if i % 2 == 0:  # HBC points
        col.set_offsets(col.get_offsets() - np.array([0.2, 0]))
    else:  # GBC points
        col.set_offsets(col.get_offsets() + np.array([0.2, 0]))

# Customize the plot
plt.title('Palantir Fate Probabilities for GBCs and iSus', fontsize=32)
plt.xlabel('Cell Type', fontsize=22)
plt.ylabel('Probability', fontsize=22)
plt.ylim(-0.05, 1.2)
plt.xticks(rotation=45, ha='right', fontsize = 22)
plt.yticks(y_ticks, fontsize = 22)
plt.legend(fontsize=22, bbox_to_anchor=(1.06, 1), loc='upper right', edgecolor = 'black')
plt.grid(False)

# Function to add significance stars
def get_stars(p):
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    else:
        return "ns"

# Perform t-tests and add significance bars with stars
for i, cell_type in enumerate(cell_types):
    gbc_probs = gbc_fate_probs_melted[gbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    iSus_probs = iSus_fate_probs_melted[iSus_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    
    t_stat, p_value = stats.ttest_ind(hbc_probs, iSus_probs)
    
    y_max = max(gbc_probs.max(), iSus_probs.max())
    y_bar = y_max + 0.05
    
    plt.plot([i-0.2, i+0.2], [y_bar, y_bar], color='black', linewidth=3)
    stars = get_stars(p_value)
    plt.text(i, y_bar + 0.02, stars, ha='center', va='bottom', fontsize=22)
#    plt.text(i, y_bar + 0.06, f'p={p_value:.3f}', ha='center', va='bottom', fontsize=22)

# Calculate average probabilities for each fate
avg_probs_gbc = gbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()
avg_probs_iSus = iSus_fate_probs_melted.groupby('Cell Type')['Probability'].mean()

# Perform chi-squared test between HBC and GBC
chi2, df, p_value = chi_squared_test(avg_probs_gbc.values, avg_probs_iSus.values)

# Add text box with chi-squared test results
textstr = '\n'.join((
    f'χ²={chi2:.2f}',
    f'df={df}',
    f'p={p_value:.4e}'
))

# Position the text box in figure coords
props = dict(boxstyle='round', facecolor='white', alpha=1, edgecolor = 'black')
plt.text(0.90, 0.95, textstr, transform=plt.gcf().transFigure, fontsize=22,
         verticalalignment='top', bbox=props)

# Adjust layout and show the plot
#plt.tight_layout()
plt.subplots_adjust(top=0.85)  # Make room for the text box
#plt.savefig(os.path.join(output_dir, 'final_minimalsubset_gbcvisus_compared_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[88]:


expression_threshold = 1e-5

# Initialize a dictionary to store the counts
expression_counts = {}

# Iterate through each gene in the list and count the cells expressing it
for gene in aav:
    if gene in ad.raw.to_adata().var_names:
        # Create a boolean mask for expression > 0
        mask = ad.raw.to_adata()[:, gene].X > expression_threshold
        # Sum the mask to get the count of cells expressing the gene
        count = mask.sum()
        # Store the count in the dictionary
        expression_counts[gene] = count
    else:
        print(f"Gene {gene} not found in the dataset.")

# Convert the dictionary to a DataFrame for easy viewing
expression_counts_df = pd.DataFrame.from_dict(expression_counts, orient='index', columns=['Count'])

# Print the DataFrame
print(expression_counts_df)

# Optionally, visualize the counts
expression_counts_df.plot(kind='bar', legend=False)
plt.title('Frequency of Cells Expressing Each AAV Gene')
# Remove gridlines
plt.grid(False)
plt.xlabel('Gene')
plt.ylabel('Number of Cells')
#plt.savefig(os.path.join(folder, 'aav_total_count.png'), dpi = 300)
plt.show()


# In[89]:


sum(expression_counts_df['Count'])


# In[ ]:


# Initialize a DataFrame to store the counts
expression_counts = pd.DataFrame(index=ad.obs['celltype'].unique(), columns=aav, data=0)

# Iterate through each gene in the list
for gene in aav:
    if gene in ad.raw.to_adata().var_names:
        # Create a boolean mask for expression > 0
        mask = ad.raw.to_adata()[:, gene].X > expression_threshold
        # Convert the mask to a DataFrame and add the cell type information
        mask_df = pd.DataFrame(mask.toarray(), index=ad.obs_names, columns=[gene])
        mask_df['celltype'] = ad.obs['celltype'].values
        # Group by cell type and sum the boolean values
        count_df = mask_df.groupby('celltype')[gene].sum()
        # Update the expression_counts DataFrame
        expression_counts[gene] = count_df
    else:
        print(f"Gene {gene} not found in the dataset.")

# Fill NaN values with 0
expression_counts = expression_counts.fillna(0)

# Print the DataFrame
print(expression_counts)

# Optionally, visualize the counts for each gene stratified by cell type
expression_counts.plot(kind='bar', stacked=True, figsize=(12, 8))
# Remove gridlines
plt.grid(False)
plt.title('Frequency of Cells Expressing Each AAV Gene by Cell Type')
plt.xlabel('Cell Type')
plt.ylabel('Number of Cells')
plt.legend(title='Gene')
plt.xticks(rotation=90)
#plt.savefig(os.path.join(folder, 'aav_per_cell_count.png'), dpi = 300)
plt.show()


# ### Work with subsetted object

# #### DPT/PAGA First

# In[284]:


ad.obs.celltype.cat.categories


# In[348]:


celltypes_to_keep = ['HBC', 'GBC', 'iOSN', 'mOSN', 'mSus', 'iSus', 'BowG', 'Gob', 'MVC', 'mOEC']
celltypes_to_keep.remove('Gob')
#celltypes_to_keep.append('mOEC')
celltypes_to_keep.remove('mOEC')


# In[349]:


print(celltypes_to_keep)


# In[286]:


ad


# In[350]:


ad2 = ad[ad.obs.celltype.isin(celltypes_to_keep)].copy()
# drop immune, LE, VE, RE, Mes, UnSec
# drop mOEC, ACC, SCC due to issues with pseudotime


# In[351]:


# Step 1: Extract the original color palette from the anndata object
original_celltype_colors = ad.uns['celltype_colors']
original_cell_types = ad.obs['celltype'].cat.categories

# Step 2: Create a dictionary mapping cell types to colors
color_palette = {cell_type: color for cell_type, color in zip(original_cell_types, original_celltype_colors)}

# Step 3: Subset the anndata object
#already subsetted

# Step 4: Ensure the subsetted anndata object uses the same color palette
# Extract cell types present in the subset
subset_cell_types = ad2.obs['celltype'].cat.categories

# Create a custom palette for the subset based on the original color mapping
subset_custom_palette = [color_palette[cell_type] for cell_type in subset_cell_types]

# Add the custom palette into original object
ad2.uns['celltype_colors'] = subset_custom_palette


# In[352]:


ad2.uns['neighbors']


# In[357]:


rsc.pp.neighbors(ad2, n_neighbors= 30, n_pcs = 30, use_rep = 'X_scvi') #normally use X_scvi for latent space
sc.tl.umap(ad2, method = 'rapids')
sc.pl.umap(ad2, color = 'celltype', **kwargs)
rsc.tl.mde(ad2, n_neighbors=15, n_pcs = 15)


# In[359]:


ad2_mde = sc.pl.embedding(ad2, basis = 'mde', color = 'celltype', **kwargs, return_fig = True, title = '')
ad2_mde.savefig(os.path.join(output_dir, 'final_minimalsubsetMDEprojection.png'), **plot_kwargs)


# In[360]:


rsc.tl.draw_graph(ad2)
sc.pl.draw_graph(ad2, color = 'celltype', **kwargs)


# In[361]:


sc.tl.paga(ad2, groups = 'celltype')


# In[362]:


sc.pl.paga(ad2, color='celltype', threshold = 0.1,
                        #edge_width_scale=0.5,
                        frameon=False,
                        fontsize=6,
                        show=False,
                        save=False)


# In[363]:


rsc.tl.draw_graph(ad2, init_pos = 'paga')


# In[365]:


paga_compare = sc.pl.paga_compare(ad2, basis='X_umap', color='celltype', 
                                  edge_width_scale=0.5,
                                  frameon=False,
                                  fontsize=6, 
                                  show=False,
                                  save=False,
                                  title = '')


# In[366]:


root_cell1 = np.flatnonzero(ad2.obs["celltype"] == "HBC")[0]
root_cell2 = np.where(ad2.obs['celltype'] == 'HBC')[0][0]
print(root_cell1)
print(root_cell2)


# In[367]:


ad2.uns['iroot'] = root_cell1


# In[368]:


rsc.tl.diffmap(ad2, n_comps = 15)
rsc.pp.neighbors(ad2, n_neighbors=15, use_rep='X_diffmap')


# In[369]:


sc.tl.dpt(ad2, n_dcs = 10)


# In[370]:


sc.pl.draw_graph(ad2, color = ['celltype', 'dpt_pseudotime', 'dpt_order'],
                 legend_loc = 'on data', legend_fontsize = 6,
                  cmap = 'Spectral_r', ncols = 2)
plt.show()
plt.close()


# In[371]:


sc.pl.embedding(ad2, basis = 'mde', color = ['celltype', 'dpt_pseudotime', 'dpt_order'],
                 legend_loc = 'on data', legend_fontsize = 6,
                  cmap = 'Spectral_r', ncols = 2)
plt.show()
plt.close()


# #### Next work with palantir

# In[372]:


# Run diffusion maps
pca_rd = pltr.utils.run_pca(ad2)
dff_map = pltr.utils.run_diffusion_maps(ad2, n_components = 30)


# In[373]:


ms_data = pltr.utils.determine_multiscale_space(ad2)


# In[374]:


imputed_X = pltr.utils.run_magic_imputation(ad2, n_jobs = -1)


# In[375]:


pltr.plot.plot_diffusion_components(ad2, cmap = 'Spectral_r', embedding_basis = 'X_mde')
plt.show()
plt.close()


# In[380]:


terminal_states = pd.Series(
    [#"mOSN1", 
     #'mOSN2',
     #'mOSN3',
     'mOSN',
     #'mOSN5',
     #'mOSN6', 
     "BowG", 
     #"Gob", 
     "mSus", 
    # "mOEC" ,
     "MVC",
 #    "ACC",
     ],
    index=[
      #  "AACAGGGGTCCCGTGA",  #mOSN1
      #  'AATTCCTTCCTGGGTG', #mOSN2
        #'CCGCGTTAGTCAGAGC', #mOSN3
        'CAGCACGAGAGATGGT', #mOSN4
        #'CACCAAGAGCTTAGGA', #mOSN5
        #'ACCATTTCTCTCCCTA', #mOSN6
        "GCTGAATCGTAGTATC",  # Bowman's Gland
         # "ATGAGTCTCCACCGGG", #gob
           "GAAGTAAAGGAGCTAG", #mature Sus
     #      "CAGGCCACTTTTAAAA",  #mOEC
          "AGCCAGCAGAAGGGAT", #MVC
 #          "CGCATGGAGAAGGATG", #ACC
    ]
)


# In[45]:


ad.obs.loc[ad.obs['celltype'] == 'mOSN'].tail(20)


# In[381]:


pltr.plot.highlight_cells_on_umap(ad2, terminal_states, embedding_basis = 'X_mde')
plt.show()
plt.close()


# In[378]:


pltr.utils.early_cell(ad2, 'HBC')


# In[382]:


#CGACTTGTCCATACTT, GGCACTTGTGACCCGC = GBC
#ACAGCGGGTATCAGGG, CCCTCCCTACTCAGTT, AGAGCCCTCCCTATTA, ATGGGTCTCCCGTTGT = HBC
pr_res = pltr.core.run_palantir(
    ad2, "ATGGGTCTCCCGTTGT", 
    terminal_states=terminal_states, 
    #max_iterations = 50,
    #num_waypoints = 2000, 
    knn = 30,
    use_early_cell_as_start= True,
)


# In[383]:


output_dir


# In[384]:


pseudotime_ptr_results = pltr.plot.plot_palantir_results(ad2, embedding_basis= 'X_mde',
                                                         s = 5, cmap = 'Spectral_r')
plt.savefig(os.path.join(output_dir, 'final_minimalsubsetMDEprojection_palantir.png'), **plot_kwargs)
plt.show()
plt.close()


# In[385]:


hbc_index = np.where(ad2.obs['celltype'] == 'HBC')[0]
hbc_barcodes = ad2.obs_names[hbc_index]


# In[386]:


gbc_index = np.where(ad2.obs['celltype'] == 'GBC')[0]
gbc_barcodes = ad2.obs_names[gbc_index]


# In[387]:


iSus_index = np.where(ad2.obs['celltype'] == 'iSus')[0]
iSus_barcodes = ad2.obs_names[iSus_index]


# In[388]:


fate_probs = ad2.obsm['palantir_fate_probabilities']


# In[389]:


# Select HBC fate probabilities
hbc_fate_probs = fate_probs.iloc[hbc_index]

# Now we can melt the DataFrame
hbc_fate_probs_melted = hbc_fate_probs.reset_index().melt(id_vars='barcode', var_name='Cell Type', value_name='Probability')

# Select GBC fate probabilities
gbc_fate_probs = fate_probs.iloc[gbc_index]

# Now we can melt the DataFrame
gbc_fate_probs_melted = gbc_fate_probs.reset_index().melt(id_vars='barcode', var_name='Cell Type', value_name='Probability')

hbc_fate_probs_melted['group'] = 'HBC'
gbc_fate_probs_melted['group'] = 'GBC'


# In[390]:


iSus_fate_probs = fate_probs.iloc[iSus_index]
iSus_fate_probs_melted = iSus_fate_probs.reset_index().melt(id_vars='barcode', var_name='Cell Type', value_name='Probability')
iSus_fate_probs_melted['group'] = 'iSus'


# In[391]:


# Extract the color palette from the anndata object
celltype_colors = ad2.uns['celltype_colors']
cell_types_in_order = ad2.obs['celltype'].cat.categories  # Ensure this matches the cell types in the correct order

# Create a dictionary mapping cell types to colors
color_palette = {cell_type: color for cell_type, color in zip(cell_types_in_order, celltype_colors)}

# Create a custom color palette for seaborn
custom_palette = [color_palette[cell_type] for cell_type in hbc_fate_probs_melted['Cell Type'].unique()]

# Set up the plot
plt.figure(figsize=(16, 8))

# Create the violin plot with modifications
sns.violinplot(x='Cell Type', y='Probability', data=hbc_fate_probs_melted,
               inner='box',
               cut=0,
               scale='width',
               width=0.8,
               palette=custom_palette,
               linewidth=2)

# Add individual points
sns.swarmplot(x='Cell Type', y='Probability', data=hbc_fate_probs_melted,
              color='black', edgecolor='black', size=3, alpha=0.5)

# Customize the plot
plt.title('Palantir Fate Probabilities for HBC', fontsize=32)
plt.xlabel('Cell Type', fontsize=32)
plt.ylabel('Probability', fontsize=32)
plt.yticks(fontsize = 22)
plt.xticks(rotation=45, ha='right', fontsize = 22)
plt.ylim(-0.05, 1.05)
plt.grid(False)

# Improve layout
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'hbc_average_fate_probabilities.png'))

# Show the plot
plt.show()
plt.close()


# In[412]:


sc.pl.embedding(ad2, basis = 'X_mde', color = ['celltype','dpt_pseudotime', 'palantir_pseudotime', 'palantir_entropy'], legend_loc = None, ncols = 2,
                title = ['','', '',''], frameon=False)


# In[399]:


# Set up the plot
plt.figure(figsize=(16, 8))  # Increased figure size for better visibility

# Create the violin plot with modifications
sns.violinplot(x='Cell Type', y='Probability', data=iSus_fate_probs_melted,
               inner='box',  # This will show the box plot inside the violin
               cut=0,  # This extends the violin plot to the full range of the data
               scale='width',  # This makes all violins the same width
               width=0.8,  # Adjust this to make violins wider or narrower
               palette=custom_palette,
               linewidth=2)  # Change color palette (you can try others like 'Set2', 'deep', etc.) (husl is good)

# Add individual points
sns.swarmplot(x='Cell Type', y='Probability', data=iSus_fate_probs_melted,
              color='black', edgecolor='black', size=3, alpha=0.5)

# Customize the plot
plt.title('Palantir Fate Probabilities for iSus', fontsize=16)
plt.xlabel('Cell Type', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.ylim(-0.05, 1.05)  # Adjust y-axis limits if needed

# Improve layout
plt.tight_layout()
#plt.savefig(os.path.join(pltir_output_dir, 'trans_average_fate_probabilities.png'), **plot_kwargs)

# Show the plot
plt.show()
plt.close()


# In[401]:


# Combine the data
combined_data = pd.concat([hbc_fate_probs_melted, gbc_fate_probs_melted])

# Set up the plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

color_palette = sns.color_palette(palette='Set1')

# Plot HBC data
sns.violinplot(x='Cell Type', y='Probability', data=hbc_fate_probs_melted, ax=ax1, 
               palette=custom_palette, inner='box', cut=0, scale='width', width=0.8)
sns.swarmplot(x='Cell Type', y='Probability', data=hbc_fate_probs_melted, ax=ax1,
              color='black', edgecolor='white', size=3, alpha=0.5)
ax1.set_title('Palantir Fate Probabilities for HBCs', fontsize=16)
ax1.set_ylim(-0.05, 1.05)

# Plot GBC data
sns.violinplot(x='Cell Type', y='Probability', data=gbc_fate_probs_melted, ax=ax2, 
               palette=custom_palette, inner='box', cut=0, scale='width', width=0.8)
sns.swarmplot(x='Cell Type', y='Probability', data=gbc_fate_probs_melted, ax=ax2,
              color='black', edgecolor='white', size=3, alpha=0.5)
ax2.set_title('Palantir Fate Probabilities for GBCs', fontsize=16)
ax2.set_ylim(-0.05, 1.05)

# Customize both plots
for ax in [ax1, ax2]:
    ax.set_xlabel('Cell Type', fontsize=14)
    ax.set_ylabel('Probability', fontsize=14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
plt.tight_layout()
#plt.savefig(os.path.join(pltir_output_dir, 'hbc_gbc_sidebyside_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[395]:


def chi_squared_test(a_probs, b_probs):
    # Ensure both arrays have the same length
    assert len(a_probs) == len(b_probs), "HBC and GBC probability arrays must have the same length"
    
    # Calculate chi-squared statistic
    chi2_stat = np.sum((a_probs - b_probs)**2 / b_probs)
    
    # Degrees of freedom
    df = len(a_probs) - 1
    
    # Calculate p-value
    p_value = 1 - stats.chi2.cdf(chi2_stat, df)
    
    return chi2_stat, df, p_value


# In[409]:


# Combine the data
combined_data = pd.concat([hbc_fate_probs_melted, gbc_fate_probs_melted])

# Set up the plot
plt.figure(figsize=(20, 10))
color_palette = sns.color_palette("Set1")

# Create a custom order for cell types and groups
cell_types = hbc_fate_probs_melted['Cell Type'].unique()
groups = ['HBC', 'GBC']

# Create the violin plot
ax = sns.violinplot(x='Cell Type', y='Probability', hue='group', data=combined_data,
                    split=True, inner="box", palette=color_palette[:2],
                    cut=0, scale='width', width=0.8)

# Add swarm plots with manual offset
for i, group in enumerate(groups):
    sns.swarmplot(x='Cell Type', y='Probability', data=combined_data[combined_data['group'] == group],
                  color='black', edgecolor='white', size=3, alpha=0.5,
                  ax=ax, dodge=True, 
                  zorder=1)  # Ensure points are on top

# Adjust swarmplot positions
for i, col in enumerate(ax.collections[len(cell_types)*2:]):  # Skip box plot artists
    if i % 2 == 0:  # HBC points
        col.set_offsets(col.get_offsets() - np.array([0.2, 0]))
    else:  # GBC points
        col.set_offsets(col.get_offsets() + np.array([0.2, 0]))

# Customize the plot
plt.title('Palantir Fate Probabilities for HBCs and GBCs', fontsize=32)
plt.xlabel('Cell Type', fontsize=22)
plt.ylabel('Probability', fontsize=22)
plt.ylim(-0.05, 1.2)
plt.xticks(rotation=45, ha='right', fontsize = 22)
y_ticks = np.arange(0, 1.2, 0.2)  # Create an array of y-ticks from 0 to 1 with a step of 0.1
plt.yticks(y_ticks, fontsize = 22)
plt.legend(fontsize=22, bbox_to_anchor=(1.06, 1.1), loc='upper right', edgecolor = 'black')
plt.grid(False)

# Function to add significance stars
def get_stars(p):
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    else:
        return "ns"

# Perform t-tests and add significance bars with stars
for i, cell_type in enumerate(cell_types):
    hbc_probs = hbc_fate_probs_melted[hbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    gbc_probs = gbc_fate_probs_melted[gbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    
    t_stat, p_value = stats.ttest_ind(hbc_probs, gbc_probs)
    
    y_max = max(hbc_probs.max(), gbc_probs.max())
    y_bar = y_max + 0.05
    
    plt.plot([i-0.2, i+0.2], [y_bar, y_bar], color='black', linewidth=3)
    stars = get_stars(p_value)
    plt.text(i, y_bar + 0.02, stars, ha='center', va='bottom', fontsize=22) # significant stars
    #plt.text(i, y_bar + 0.06, f'p={p_value:.3f}', ha='center', va='bottom', fontsize=22) # p-values

# Calculate average probabilities for each fate
avg_probs_hbc = hbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()
avg_probs_gbc = gbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()

# Perform chi-squared test between HBC and GBC
chi2, df, p_value = chi_squared_test(avg_probs_hbc.values, avg_probs_gbc.values)

# Add text box with chi-squared test results
textstr = '\n'.join((
    f'χ²={chi2:.2f}',
    f'df={df}',
    f'p={p_value:.2e}'
))

# Position the text box in figure coords
props = dict(boxstyle='round', facecolor='white', alpha=1, edgecolor = 'black')
plt.text(0.90, 0.78, textstr, transform=plt.gcf().transFigure, fontsize=22,
         verticalalignment='top', bbox=props)

# Adjust layout and show the plot
#plt.tight_layout()
plt.subplots_adjust(top=0.85)  # Make room for the text box
plt.savefig(os.path.join(output_dir, 'final_minimalsubset_hbcvgbc_compared_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[410]:


# Combine the data
combined_data = pd.concat([hbc_fate_probs_melted, iSus_fate_probs_melted])

# Set up the plot
plt.figure(figsize=(20, 10))
color_palette = sns.color_palette("Set1")

# Create a custom order for cell types and groups
cell_types = hbc_fate_probs_melted['Cell Type'].unique()
groups = ['HBC', 'iSus']

# Create the violin plot
ax = sns.violinplot(x='Cell Type', y='Probability', hue='group', data=combined_data,
                    split=True, inner="box", palette=color_palette[:2],
                    cut=0, scale='width', width=0.8)

# Add swarm plots with manual offset
for i, group in enumerate(groups):
    sns.swarmplot(x='Cell Type', y='Probability', data=combined_data[combined_data['group'] == group],
                  color='black', edgecolor='white', size=3, alpha=0.5,
                  ax=ax, dodge=True, 
                  zorder=1)  # Ensure points are on top

# Adjust swarmplot positions
for i, col in enumerate(ax.collections[len(cell_types)*2:]):  # Skip box plot artists
    if i % 2 == 0:  # HBC points
        col.set_offsets(col.get_offsets() - np.array([0.2, 0]))
    else:  # GBC points
        col.set_offsets(col.get_offsets() + np.array([0.2, 0]))

# Customize the plot
plt.title('Palantir Fate Probabilities for HBCs and iSus', fontsize=32)
plt.xlabel('Cell Type', fontsize=22)
plt.ylabel('Probability', fontsize=22)
plt.ylim(-0.05, 1.2)
plt.xticks(rotation=45, ha='right', fontsize = 22)
plt.yticks(y_ticks, fontsize = 22)
plt.legend(fontsize=22, bbox_to_anchor=(1.06, 1.06), loc='upper right', edgecolor = 'black')
plt.grid(False)

# Function to add significance stars
def get_stars(p):
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    else:
        return "ns"

# Perform t-tests and add significance bars with stars
for i, cell_type in enumerate(cell_types):
    hbc_probs = hbc_fate_probs_melted[hbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    iSus_probs = iSus_fate_probs_melted[iSus_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    
    t_stat, p_value = stats.ttest_ind(hbc_probs, iSus_probs)
    
    y_max = max(hbc_probs.max(), iSus_probs.max())
    y_bar = y_max + 0.05
    
    plt.plot([i-0.2, i+0.2], [y_bar, y_bar], color='black', linewidth=3)
    stars = get_stars(p_value)
    plt.text(i, y_bar + 0.02, stars, ha='center', va='bottom', fontsize=22)
    #plt.text(i, y_bar + 0.06, f'p={p_value:.3f}', ha='center', va='bottom', fontsize=22)

# Calculate average probabilities for each fate
avg_probs_hbc = hbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()
avg_probs_iSus = iSus_fate_probs_melted.groupby('Cell Type')['Probability'].mean()

# Perform chi-squared test between HBC and GBC
chi2, df, p_value = chi_squared_test(avg_probs_hbc.values, avg_probs_iSus.values)

# Add text box with chi-squared test results
textstr = '\n'.join((
    f'χ²={chi2:.2f}',
    f'df={df}',
    f'p={p_value:.2e}'
))

# Position the text box in figure coords
props = dict(boxstyle='round', facecolor='white', alpha=1, edgecolor = 'black')
plt.text(0.89, 0.75, textstr, transform=plt.gcf().transFigure, fontsize=22,
         verticalalignment='top', bbox=props)

# Adjust layout and show the plot
#plt.tight_layout()
plt.subplots_adjust(top=0.85)  # Make room for the text box
plt.savefig(os.path.join(output_dir, 'final_minimalsubset_hbcvisus_compared_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[408]:


# Combine the data
combined_data = pd.concat([gbc_fate_probs_melted, iSus_fate_probs_melted])

# Set up the plot
plt.figure(figsize=(20, 10))
color_palette = sns.color_palette("Set1")

# Create a custom order for cell types and groups
cell_types = gbc_fate_probs_melted['Cell Type'].unique()
groups = ['GBC', 'iSus']

# Create the violin plot
ax = sns.violinplot(x='Cell Type', y='Probability', hue='group', data=combined_data,
                    split=True, inner="box", palette=color_palette[:2],
                    cut=0, scale='width', width=0.8)

# Add swarm plots with manual offset
for i, group in enumerate(groups):
    sns.swarmplot(x='Cell Type', y='Probability', data=combined_data[combined_data['group'] == group],
                  color='black', edgecolor='white', size=3, alpha=0.5,
                  ax=ax, dodge=True, 
                  zorder=1)  # Ensure points are on top

# Adjust swarmplot positions
for i, col in enumerate(ax.collections[len(cell_types)*2:]):  # Skip box plot artists
    if i % 2 == 0:  # HBC points
        col.set_offsets(col.get_offsets() - np.array([0.2, 0]))
    else:  # GBC points
        col.set_offsets(col.get_offsets() + np.array([0.2, 0]))

# Customize the plot
plt.title('Palantir Fate Probabilities for GBCs and iSus', fontsize=32)
plt.xlabel('Cell Type', fontsize=32)
plt.ylabel('Probability', fontsize=32)
plt.ylim(-0.05, 1.2)
plt.xticks(rotation=45, ha='right', fontsize = 22)
plt.yticks(y_ticks, fontsize = 22)
plt.legend(fontsize=32, bbox_to_anchor=(1, 1.06), loc='upper right', edgecolor = 'black')
plt.grid(False)

# Function to add significance stars
def get_stars(p):
    if p <= 0.0001:
        return "****"
    elif p <= 0.001:
        return "***"
    elif p <= 0.01:
        return "**"
    elif p <= 0.05:
        return "*"
    else:
        return "ns"

# Perform t-tests and add significance bars with stars
for i, cell_type in enumerate(cell_types):
    gbc_probs = gbc_fate_probs_melted[gbc_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    iSus_probs = iSus_fate_probs_melted[iSus_fate_probs_melted['Cell Type'] == cell_type]['Probability']
    
    t_stat, p_value = stats.ttest_ind(hbc_probs, iSus_probs)
    
    y_max = max(gbc_probs.max(), iSus_probs.max())
    y_bar = y_max + 0.05
    
    plt.plot([i-0.2, i+0.2], [y_bar, y_bar], color='black', linewidth=5)
    stars = get_stars(p_value)
    plt.text(i, y_bar + 0.02, stars, ha='center', va='bottom', fontsize=32)
#    plt.text(i, y_bar + 0.06, f'p={p_value:.3f}', ha='center', va='bottom', fontsize=22)

# Calculate average probabilities for each fate
avg_probs_gbc = gbc_fate_probs_melted.groupby('Cell Type')['Probability'].mean()
avg_probs_iSus = iSus_fate_probs_melted.groupby('Cell Type')['Probability'].mean()

# Perform chi-squared test between HBC and GBC
chi2, df, p_value = chi_squared_test(avg_probs_gbc.values, avg_probs_iSus.values)

# Add text box with chi-squared test results
textstr = '\n'.join((
    f'χ²={chi2:.2f}',
    f'df={df}',
    f'p={p_value:.2e}'
))

# Position the text box in figure coords
props = dict(boxstyle='round', facecolor='white', alpha=1, edgecolor = 'black')
plt.text(0.82, 0.70, textstr, transform=plt.gcf().transFigure, fontsize=32,
         verticalalignment='top', bbox=props)

# Adjust layout and show the plot
#plt.tight_layout()
plt.subplots_adjust(top=0.85)  # Make room for the text box
plt.savefig(os.path.join(output_dir, 'final_minimalsubset_gbcvisus_compared_average_fate_probabilities.png'), **plot_kwargs)

plt.show()
plt.close()


# In[413]:


ad2.raw.to_adata()


# In[414]:


input_dir


# In[415]:


print(ad.uns['celltype_colors'])
print(ad2.uns['celltype_colors'])


# In[416]:


ad2.write(os.path.join(input_dir, 'minimalsubsetobject_palantir_dpt.h5ad'))


# # Stop here

# In[94]:


masks = pltr.presults.select_branch_cells(ad2, q=.005, eps=.005)


# In[96]:


terminal_states.values.tolist()


# In[99]:


ad2.obsm['branch_masks']


# In[100]:


gene_trends = pltr.presults.compute_gene_trends(
    ad2,
    expression_key="MAGIC_imputed_data",
)


# In[101]:


genes_of_interest = ['Ano2',
'Olfm1',
'Chga',
'Kcnk10',
'Cnga2',
'Sox11',
'Dpysl5', 'Dpysl2', 'Hdac',
'Lhx2',
'Ncam2']
olf_genes = ['Or6b2', 'Or6b2b', 'Or6b3', 'Or10j27', 'Or12k8']
osn_maturation_genes = ["Ascl1", "Neurog1", "Neurog2", "Tbx21", "Emx2", "Sox2", "Sox11", "Ngn1", "Lhx2", "Lhx9", "Ebf1", "Ebf2", "Mash1", "Otx2", "Gap43", "Neurod1", "Dlx5", "Dlx6", "Pax6", "Zic1", "Zic2", "Foxg1", "Irx3", "Bdnf", "Gdnf"
]
axon_targeting_genes= ["Reln", "Shh", "Wnt3a", "Dcx", "Nrg1", "Nrp1", "EphA4", "Bdnf", "Slit1", "Robo1", "Robo2", "PSA-NCAM", "Prok2", "Tbr2", "Sema3A", "Nrp2", "Gad1", "EphA5", "EphA6", "Efna5", "Ntn1", "Dcc", "Cdh1", "Cdh2", "Pcdh", "L1CAM",] 
axon_targeting_genes2= ["Tnr", "Cntn2", "Ctnnb1", "Sema7A", "Robo3", "Itga3", "Itgb1", "PlxnA1", "Robo1", "Robo2", "PSA-NCAM", "Cdk5", "Sema3C", "Dscam", "Lhx2", "Slitrk1", "Slitrk2", "Gfra1", "Lhx9", "Gli1", "Gli2", "Lrig1", "Erbb4", "Nrn1", "Slit3", "Bcl11b", "FOXP2", "Nrp2", "EphB2", "Atf5", "Axin2", "Nrp1", "Tnc", "Zic2", "Or2p2", "Or2ab1", "Or9e1", "Or2f1"
]
other_genes_of_interest = ["Nrp2", "PlxnA4", "NrCAM", "Fzd3", "Celsr3", "Dscam", "Slit2", "Sox11", "Gfra1", "Robo4", "Arhgap21", "Cux1", "Cux2", "Spinophilin", "Kalirin", "Neurod2", "Pcdh10", "Map2", "Shank3", "Dlg4", "Bdnf", "Camk2a", "Arpc3", "Mef2c", "Notch1", "Omp", "Gnai2", "Gucy1b2", "Adcy3", "Pax6", "Cyp2g1", "Slc8a1"
]


# In[102]:


genes = olf_genes + genes_of_interest + other_genes_of_interest +  axon_targeting_genes + axon_targeting_genes2 + osn_maturation_genes
print(len(genes))
genes = list(set(genes))
print(len(genes))


# In[103]:


# Filter genes to keep only those present in the Scanpy object
genes = [gene for gene in genes if gene in ad.var_names]

print(f"Number of genes after filtering: {len(genes)}")

