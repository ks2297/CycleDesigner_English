import os
import csv
import re
from pathlib import Path
from calc_position_aligned_CaRMSD import main  # Renamed from: 计算基于位置对齐的CaRMSD

def find_native_file(native_dir, pdb_file_name):
    """
    Find the corresponding native file based on the PDB ID.
    """
    # prefix = pdb_file_name[:4]
    prefix = pdb_file_name[:6]
    native_file = next(Path(native_dir).glob(f"{prefix}*.pdb"), None)
    return native_file

def extract_rank_number(pdb_file_name):
    """
    Extract the rank number.
    """
    match = re.search(r"_rank_00(\d+)", pdb_file_name)
    return match.group(1) if match else None

def process_files_and_write_csv(pdb_dir, native_dir, output_file):
    # pdb_dir = "/home/yons/lab/zch/rfd2HF/outputdir/1023/1023_30T"
    # native_dir = "/home/yons/inputs/pdbfixer/rfd_native/fixed"
    # output_file = "rmsd_results.csv"

    results = []

    # Iterate through all subdirectories and PDB files in the PDB directory
    for root, _, files in os.walk(pdb_dir):
        for pdb_file in files:
            if pdb_file.endswith(".pdb"):
                pdb_file_path = os.path.join(root, pdb_file)
                rank_number = extract_rank_number(pdb_file)

                # Find the matching native file
                native_file = find_native_file(native_dir, pdb_file)
                if native_file:
                    try:
                        # Calculate RMSD
                        rmsd_value = main(pdb_file_path, native_file)
                        results.append([pdb_file, rank_number,native_file.name, rmsd_value])
                    except Exception as e:
                        print(f"Error calculating RMSD for {pdb_file} and {native_file.name}: {e}")

    # Write results to CSV file
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PDB File", "Rank", "Native File", "RMSD"])
        writer.writerows(results)

    print(f"RMSD calculation complete. Results written to {output_file}")

if __name__ == "__main__":
    process_files_and_write_csv()
