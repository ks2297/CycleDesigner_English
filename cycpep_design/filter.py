import pandas as pd
import os
import shutil

def filter(csv_file, to_rosetta):
    file_path = csv_file
    output_path = to_rosetta
    os.makedirs(output_path, exist_ok=True)

    # read csv
    df = pd.read_csv(file_path)

    # Filter data that meets the criteria
    # filtered_df = df[(df['plddt'] > 80) & (df['iptm'] > 0.5) & (df['ipae'] < 0.35) & (df['rmsd'] < 3.5)]

    # Save filtered data to a new CSV file
    # filtered_df.to_csv('filtered_data.csv', index=False)

    # Filter data that meets the criteria
    filtered_df = df[(df['plddt'] > 80) & (df['iptm'] > 0.5) & (df['ipae'] < 0.35) & (df['rmsd'] < 3.5)]

    # Iterate through filtered rows and copy files
    for _, row in filtered_df.iterrows():
        # Build the full file path
        src_path = os.path.join(
            row['folder'], f"{row['protein_id']}_{row['backbone']}_{row['seqid']}_unrelaxed_{row['file']}.pdb"
        )
        
        if os.path.exists(src_path):
            # Copy the file to the destination folder
            dest_path = os.path.join(output_path, os.path.basename(src_path))
            shutil.copy(src_path, dest_path)
            # print(f"Copied file: {src_path} to {dest_path}")
        else:
            print(f"File does not exist, cannot copy: {src_path}")

    return "done"
