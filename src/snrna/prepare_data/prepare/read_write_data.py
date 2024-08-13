import os
import scanpy as sc
import pandas as pd
import anndata as ad

def read_files(sample, type = 'spliced'):
    input_dir = '/home/johnathan/projects/arrenkiel_scrnaseq/align_count/output/starsolo/'
    
    file_dir = input_dir + 'ASAP-' + sample + '/Solo.out/Velocyto/filtered/'
    
    barcodes = pd.read_csv(file_dir + 'barcodes.tsv', header=None, sep='\t')
    genes = pd.read_csv(file_dir + 'features.tsv', header=None, sep='\t')
    mtx = sc.read_mtx(file_dir + type+'.mtx')
    adata = ad.AnnData.transpose(mtx)
    
    return barcodes, genes, adata        

def setup_files(sample, type = 'spliced'):
    barcodes, genes, adata = read_files(sample, type)
    adata.obs_names = barcodes[0].values
    adata.var_names = genes[1].values
    return adata

def save_files(adata, filename, output_dir):
    output_path = os.path.join(output_dir, filename)
    adata.write(output_path)
