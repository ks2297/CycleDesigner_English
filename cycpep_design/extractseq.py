import os
from Bio.PDB import PDBParser, Polypeptide
from Bio.SeqUtils import seq1

def extract_sequences_from_pdb(file_path):
    """Extract sequences from all chains in a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(file_path).split('.')[0], file_path)
    sequences = []

    for model in structure:
        for chain in model:
            sequence = []
            for residue in chain:
                if Polypeptide.is_aa(residue, standard=True):  # Check if it is a standard amino acid
                    aa = seq1(residue.resname)  # Convert three-letter code to one-letter code
                    sequence.append(aa)
            if sequence:
                sequences.append("".join(sequence))

    # Append ":" to the end of the first chain's sequence
    if sequences:
        sequences[0] += ":"

    return sequences

def write_fasta(pdb_id, sequences, output_folder):
    """Write extracted sequences to a FASTA file."""
    fasta_content = f">{pdb_id}\n" + "\n".join(sequences)
    output_path = os.path.join(output_folder, f"{pdb_id}.fasta")
    with open(output_path, "w") as f:
        f.write(fasta_content)

def main(input_folder, output_folder):
    """Iterate through PDB files, extract sequences, and write to FASTA files."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".pdb"):
            pdb_path = os.path.join(input_folder, filename)
            pdb_id = os.path.splitext(filename)[0]
            sequences = extract_sequences_from_pdb(pdb_path)
            write_fasta(pdb_id, sequences, output_folder)

if __name__ == "__main__":
    input_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed"
    output_folder = "/home/yons/inputs/pdbfixer/rfd_native/fixed/fasta"
    main(input_folder, output_folder)
