from prepare.read_write_data import setup_files, save_files
import os

file_types = ['spliced', 'unspliced']

for ft in file_types:
    file_list = []
    vars = ['screen-1', 'screen-2', 'screen-3', 'screen-4']

    for var in vars:
        adata = setup_files(var, ft)
        file_list.append(adata)

    #output_dir = '/mnt/atlas_local/jonathan/home/projects/arrenkiel_scrnaseq/scrna/data/before_qc/'
    output_dir = '/home/johnathan/projects/arrenkiel_scrnaseq/scrna/data/velocyto/before_qc/'
    os.makedirs(output_dir, exist_ok=True)

    for i, file in enumerate(file_list):
        save_files(file, vars[i] +'_'+ ft + '.h5ad', output_dir)
