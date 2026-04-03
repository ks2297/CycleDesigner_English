import os

def check_subfolders_for_pdb_files(parent_folder):
    # Iterate through each subfolder in the parent folder
    for subfolder in os.listdir(parent_folder):
        subfolder_path = os.path.join(parent_folder, subfolder)
        
        # Make sure it is a directory
        if os.path.isdir(subfolder_path):
            contains_unrelaxed_rank_file = False
            
            # Check all files in this subfolder
            for filename in os.listdir(subfolder_path):
                if "unrelaxed_rank" in filename and filename.endswith(".pdb"):
                    contains_unrelaxed_rank_file = True
                    break
            
            # If no matching file is found, print the subfolder name
            if not contains_unrelaxed_rank_file:
                print(f"Subfolder '{subfolder}' does not contain an 'unrelaxed_rank' PDB file")

# Example usage
parent_folder_path = "/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T"
check_subfolders_for_pdb_files(parent_folder_path)
