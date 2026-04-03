import os
import re

def extract_pdbid_from_filename(fa_file):
    """Extract the PDB ID from an FA filename."""
    filename = os.path.basename(fa_file)
    pdbid = filename.split('_')[0]  # Extract the '1sfi' part
    return pdbid

def read_sequences_from_fa(fa_file):
    """Read all valid sequences from an FA file."""
    sequences = []
    with open(fa_file, 'r') as f:
        content = f.read()
        # Match lines containing only uppercase letter sequences
        all_sequences = re.findall(r'^[A-Z]+$', content, re.MULTILINE)
        
        # Filter out sequences composed entirely of G
        for seq in all_sequences:
            if set(seq) != {'G'}:  # If the sequence is not all G, it is valid
                sequences.append(seq)

    if len(sequences) < 5:
        raise ValueError(f"{fa_file} contains fewer than five valid sequences")
    return sequences[:5]

def read_fasta_sequences(fasta_path):
    """Read chain sequences from a FASTA file."""
    sequences = []
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        # Extract sequence lines that don't start with '>'
        for line in lines:
            if not line.startswith(">"):
                sequences.append(line.strip())
    return sequences

def replace_chain_and_save(pdb_id, original_sequences, new_sequence, output_folder, index, fa_index):
    """Replace the second chain's sequence and save as a new FASTA file."""
    modified_sequences = original_sequences.copy()
    if len(modified_sequences) > 1:
        modified_sequences[1] = new_sequence  # Replace the second chain

    # Generate new FASTA file content
    fasta_content = f">{pdb_id}_{fa_index}_{index}\n" + "\n".join(modified_sequences)
    output_path = os.path.join(output_folder, f"{pdb_id}_{fa_index}_{index}.fasta")
    with open(output_path, "w") as f:
        f.write(fasta_content)
    # print(f"Generated new file: {output_path}")

def process_fa_file(fa_file, fasta_folder, output_folder, fa_index):
    """Process a single FA file and replace the second chain in the matching FASTA file."""
    pdb_id = extract_pdbid_from_filename(fa_file)
    new_sequences = read_sequences_from_fa(fa_file)

    matched_file = f"{pdb_id}.fasta"
    fasta_path = os.path.join(fasta_folder, matched_file)

    if os.path.exists(fasta_path):
        # Read chain sequences from the original FASTA file
        original_sequences = read_fasta_sequences(fasta_path)

        # Replace the second chain with each of the five new sequences and save
        for i, new_seq in enumerate(new_sequences, start=1):
            replace_chain_and_save(pdb_id, original_sequences, new_seq, output_folder, i, fa_index)
    else:
        print(f"No matching FASTA file found: {matched_file}")

def main(fa_folder, fasta_folder):
    """Main function: iterate through FA files and process matching FASTA files."""
    output_folder = os.path.join(fa_folder, "new")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(fa_folder):
        if re.match(r".+_\d\.fa$", filename):  # Match .fa files ending with a digit
            fa_file = os.path.join(fa_folder, filename)
            fa_index = filename.split('_')[-1].split('.')[0]  # Extract the trailing digit
            process_fa_file(fa_file, fasta_folder, output_folder, fa_index)

if __name__ == "__main__":
    fa_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed/fasta/241023_50T"
    fasta_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed/fasta"
    main(fa_folder, fasta_folder)
