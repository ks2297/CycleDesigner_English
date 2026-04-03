import json
import os
import numpy as np

def fasta_read(fasta_file):
    with open(fasta_file, 'r') as f:
        list_info = f.readlines()
        peptide_length = len(list_info[1].strip())  # Calculate the peptide chain length
        return peptide_length

def interface_plddt_read(interface_file):
    with open(interface_file, 'r') as f:
        list_info = f.readlines()
        chainA = [int(item) - 1 for item in list_info[0].strip().split()]  # Convert residue numbers to 0-indexed
        chainB = [int(item) - 1 for item in list_info[1].strip().split()]  # Convert residue numbers to 0-indexed
        return chainA, chainB

# Set paths
path_fasta = '/home/yons/lab/zch/rfd2HF/inputdir/1104_H_25T/1sfi_0_1.fasta'
path_json = '/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T/1sfi_0_1/1sfi_0_1_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json'
path_interface = '/path/to/interface_directory'

# List to store results
result_data = []

for fasta_file in os.listdir(path_fasta):
    # Read peptide chain length and interface information
    peptide_length = fasta_read(os.path.join(path_fasta, fasta_file))
    info = fasta_file.split('.')[0]
    chainA, chainB = interface_plddt_read(os.path.join(path_interface, f"{info.split('_')[0]}.txt"))
    
    # Read PAE data from the JSON file
    json_file = os.path.join(path_json, info, f'rank_001.json')  # Assuming the rank_001 model
    with open(json_file) as f:
        confidences = json.load(f)
        pae_matrix = np.array(confidences['pae'])

        # Calculate interface PAE
        peptide_receptor_pae = sum(pae_matrix[y][x] for y in chainA for x in chainB)
        receptor_peptide_pae = sum(pae_matrix[y][x] for y in chainB for x in chainA)
        num = len(chainA) * len(chainB)
        total_interface_pae = round((peptide_receptor_pae / num + receptor_peptide_pae / num) / 2, 2)

        result_data.append([info, total_interface_pae])

# Output results
for data in result_data:
    print(f"Peptide Name: {data[0]}, Interface PAE: {data[1]}")




import json
import numpy as np

# Read the JSON file
json_file = "your_file.json"
with open(json_file) as f:
    confidences = json.load(f)
    pae_matrix = np.array(confidences['pae'])

# Set sequence lengths
pep_len = ...  # Set the pep_len value
tar_len = ...  # Set the tar_len value

# Calculate iPAE
ipae = (np.mean(pae_matrix[:tar_len, tar_len:]) + np.mean(pae_matrix[tar_len:, :tar_len])) / 2

print("ipae:", ipae)
