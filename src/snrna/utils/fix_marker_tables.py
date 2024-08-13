import pandas as pd
import os
import argparse

# Function to modify column names
def modify_column_names(file_path):
    df = pd.read_csv(file_path)
    new_column_names = {}
    for col in df.columns:
        if col.endswith('_n'):
            new_column_names[col] = col[:-2] + '_name'
        elif col.endswith('_l'):
            new_column_names[col] = col[:-2] + '_log_fold_change'
        elif col.endswith('_p'):
            new_column_names[col] = col[:-2] + '_p_val_adj'
    df.rename(columns=new_column_names, inplace=True)
    df.to_csv(file_path, index=False)

# Define the clustering resolutions
#leiden_resolutions = ['leiden_0.1', 'leiden_0.2', 'leiden_0.5', 'leiden_1']
#leiden_resolutions = ['leiden005', 'leiden01', 'leiden02', 'leiden05', 'leiden1']
#neighbor_keys = ['scvi_mde_', 'tfor_harm_']
leiden_resolutions = ['0.1', '0.5', '1', '1.5', '2', '3', '8', '10']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify column names in cluster marker CSV files.")
    parser.add_argument('directory', type=str, help='Directory containing the CSV files.')
    args = parser.parse_args()

    # Modify column names for each CSV file
    # Loop over each resolution, compute markers, and save to CSV
 #   for key in neighbor_keys:
for res in leiden_resolutions:
    filename = os.path.join(args.directory, f'cluster_markers_leiden_{res}.csv')
    modify_column_names(filename)
    print(f'Modified column names for {filename}')
