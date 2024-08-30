import os
import shutil
import argparse

def standardize_files(input_dir, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Process each subdirectory in the input directory
    for subdir in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue

        # Create corresponding subdirectory in output
        out_subdir = os.path.join(output_dir, subdir)
        os.makedirs(out_subdir, exist_ok=True)

        # Process files in the subdirectory
        for filename in os.listdir(subdir_path):
            if filename.endswith('.txt.gz') or filename.endswith('_README.txt'):
                # Remove prefix from filename
                new_filename = filename.split('_geo_submission_', 1)[-1]
                if new_filename == filename:  # If no prefix was found
                    new_filename = filename.split('geo_submission_', 1)[-1]
                
                # Copy file to output directory with new name
                shutil.copy2(
                    os.path.join(subdir_path, filename),
                    os.path.join(out_subdir, new_filename)
                )

        print(f"Processed {subdir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standardize GEO submission files")
    parser.add_argument("input_dir", help="Path to the input directory containing all screen folders")
    parser.add_argument("output_dir", help="Path to the output directory for standardized files")
    
    args = parser.parse_args()
    
    standardize_files(args.input_dir, args.output_dir)
    print("Preprocessing completed successfully!")