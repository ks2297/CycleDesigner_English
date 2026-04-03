import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def fix_missing_residues_and_gaps(input_pdb, output_pdb):
    # Load the PDB file
    fixer = PDBFixer(filename=input_pdb)

    # Find and fill in missing residues
    # print("Finding and fixing missing amino acid residues and sequence gaps...")
    fixer.findMissingResidues()
    
    if fixer.missingResidues:
        print(f"missing residues: {fixer.missingResidues}")
    else:
        print("No missing residues found.")

    # Find and add all missing atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Handle nonstandard residues (gaps); replace them if present
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Add missing hydrogen atoms
    # print("Adding missing hydrogen atoms...")
    fixer.addMissingHydrogens()

    # Save the repaired PDB file
    print(f"fixed pdb saved to{output_pdb}")
    with open(output_pdb, 'w') as out_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        # print("Usage: python fix_pdb.py input_pdb output_pdb")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]

    fix_missing_residues_and_gaps(input_pdb, output_pdb)
