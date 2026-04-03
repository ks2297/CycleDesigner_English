import os
import pandas as pd

# Load the CSV file
csv_file_path = '/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T/log_ipae_rmsd.csv'
data = pd.read_csv(csv_file_path)

# Generate identifier column: protein_id_backbone_seqid
data['identifier'] = data['protein_id'].astype(str) + '_' + \
                     data['backbone'].astype(str) + '_' + \
                     data['seqid'].astype(str)

# Specify the folder path containing .fasta files
fasta_folder = '/home/yons/lab/zch/rfd2HF/inputdir/1104_H_25T'  # Replace with the actual path

# List to store sequences
sequences = []

# Iterate through identifiers to find corresponding .fasta files
for identifier in data['identifier']:
    fasta_file_path = os.path.join(fasta_folder, f"{identifier}.fasta")
    if os.path.exists(fasta_file_path):
        # Open and read the last line of the .fasta file
        with open(fasta_file_path, 'r') as file:
            lines = file.readlines()
            # Get the last line, stripping the newline character
            sequence = lines[-1].strip() if lines else ''
    else:
        # If the file does not exist, record as empty
        sequence = ''
    sequences.append(sequence)

# Add the sequence column to the DataFrame
data['sequence'] = sequences

# Save the updated CSV file
output_csv_path = '/home/yons/lab/zch/rfd2HF/outputdir/1104_H_25T/log_ipae_rmsd-1120.csv'
data.to_csv(output_csv_path, index=False)

print(f"Processing complete. File saved to: {output_csv_path}")
