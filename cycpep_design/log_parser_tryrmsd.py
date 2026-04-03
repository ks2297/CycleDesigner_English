import os
import re
import json
import numpy as np
import pandas as pd
from calc_position_aligned_CaRMSD import main  # Renamed from: 计算基于位置对齐的CaRMSD

def extract_plddt_from_file(log_file, fasta_dir, native_dir):
    """
    Extract pLDDT and iPAE metrics from a specified log file.
    """
    path = os.path.dirname(log_file)
    metrics = []
    
    with open(log_file, 'r') as fr:
        lines = fr.readlines()
        for line in lines:
            items = line.strip().split(' ')
            # Filter lines that start with rank_00
            if len(items) >= 3 and items[2].startswith('rank_00'):
                # Initialize record
                m = [path, items[2], -1, -1, -1]
                
                # Extract plddt, ptm, and iptm
                if len(items) >= 4:
                    m[2] = items[3].split('=')[1]
                if len(items) >= 5:
                    m[3] = items[4].split('=')[1]
                if len(items) >= 6:
                    m[4] = items[5].split('=')[1]
                
                # Parse path info to generate protein_id and other fields
                protein_id = path.split('/')[-1].split('_')[0]
                backbone = path.split('_')[1]
                seqid = path.split('_')[2]
                rank = items[2].split('_')[1][2]
                model = items[2].split('_')[-3]
                
                # Extract iPAE from the JSON scores file
                json_file_path = f"{path}/{protein_id}_{backbone}_{seqid}_scores_rank_00{rank}_alphafold2_multimer_v3_model_{model}_seed_000.json"
                fasta_file_path = f"{fasta_dir}/{protein_id}_{backbone}_{seqid}.fasta"
                if os.path.exists(json_file_path):
                    with open(json_file_path) as f:
                        confidences = json.load(f)
                        pae_matrix = np.array(confidences['pae'])
                        with open(fasta_file_path) as fasta_file:
                            for fasta_line in fasta_file:
                                fasta_line = fasta_line.strip()
                                if fasta_line.endswith(":"):
                                    tar_len = len(fasta_line) - 1
                                    pae_interaction1 = np.mean(pae_matrix[:tar_len, tar_len:])
                                    pae_interaction2 = np.mean(pae_matrix[tar_len:, :tar_len])
                                    ipae = round((pae_interaction1 + pae_interaction2) / 2, 2)
                                    ipae = ipae/31
                                    # ipae = (np.mean(pae_matrix[:tar_len, tar_len:]) + np.mean(pae_matrix[tar_len:, :tar_len])) / 2
                                    m.append(ipae)

                # Calculate RMSD against the native structure
                pdb_file_path = f"{path}/{protein_id}_{backbone}_{seqid}_unrelaxed_rank_00{rank}_alphafold2_multimer_v3_model_{model}_seed_000.pdb"
                native_file_path = f"{native_dir}/{protein_id}_{backbone}.pdb"
                if os.path.exists(pdb_file_path) and os.path.exists(native_file_path):
                    rmsd = main(pdb_file_path, native_file_path)
                    m.append(rmsd)


                # # return native_file
                # match = re.search(r"_rank_00(\d+)", pdb_file_name)
                # for root, _, files in os.walk(pdb_dir):
                #     for pdb_file in files:
                #         if pdb_file.endswith(".pdb"):
                #             pdb_file_path = os.path.join(root, pdb_file)
                #             rank_number = extract_rank_number(pdb_file)

                #             # Find the matching native file
                #             native_file = find_native_file(native_dir, pdb_file)
                #             if native_file:
                #                 try:
                #                     # Calculate RMSD
                #                     rmsd_value = main(pdb_file_path, native_file)
                #                     m.append([rmsd_value])
                #                 except Exception as e:
                #                     print(f"Error calculating RMSD for {pdb_file} and {native_file.name}: {e}")

                
                # Merge info and add to the results list
                m = [protein_id, backbone, seqid, rank] + m
                metrics.append(m)
    
    return metrics

def extract_plddt_from_folder(folder, fasta_dir, native_dir):
    """
    Walk through all log files in the specified folder and extract pLDDT and iPAE metrics.
    """
    metrics = []
    
    for root, _, files in os.walk(folder):
        for f in files:
            if f == 'log.txt':
                fullname = os.path.join(root, f)
                metrics += extract_plddt_from_file(fullname, fasta_dir, native_dir)
    
    return metrics

def write_to_file(metrics, output_file):
    """
    Write the extracted pLDDT and iPAE metrics to a CSV file.
    """
    # Adjusted columns to match the structure of metrics
    df = pd.DataFrame(metrics, columns=[
        'protein_id', 'backbone', 'seqid', 'rank', 'folder', 'file', 'plddt', 'ptm', 'iptm', 'ipae', 'rmsd'
    ])
    df.to_csv(output_file, index=False)
    print(f"Data written to {output_file}")
    return "done"
