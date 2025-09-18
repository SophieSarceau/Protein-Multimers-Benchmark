import os
from tqdm import tqdm
from Bio.PDB import MMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.PDB.Polypeptide import is_aa, three_to_one
import matplotlib.pyplot as plt
import shutil

cif_dir = './cif_files'
output_dir = './processed_files'
os.makedirs(output_dir, exist_ok=True)

parser = MMCIFParser(QUIET=True)
total_lengths = []
pdb_count = 0

def extract_chain_sequences(structure):
    """Extracts protein sequences from a PDB structure."""
    chain_seqs = {}
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if is_aa(residue, standard=True):
                    try:
                        seq += three_to_one(residue.get_resname())
                    except KeyError:
                        # Skip unknown or modified amino acids
                        continue
            if seq:
                chain_seqs[chain.id] = seq
    return chain_seqs

for cif_file in tqdm(os.listdir(cif_dir), desc="Processing CIF files"):
    if not cif_file.endswith('.cif'):
        continue

    pdb_id = cif_file.split('.')[0]
    cif_path = os.path.join(cif_dir, cif_file)

    try:
        structure = parser.get_structure(pdb_id, cif_path)
        chain_seqs = extract_chain_sequences(structure)

        # Calculate total length for the entire PDB
        total_len = sum(len(seq) for seq in chain_seqs.values())

        # Skip if total length is 0 or exceeds the threshold
        if total_len == 0 or total_len > 5120:
            continue

        # Create a subdirectory for the PDB ID
        pdb_out_dir = os.path.join(output_dir, pdb_id)
        os.makedirs(pdb_out_dir, exist_ok=True)

        # Copy the original CIF file
        shutil.copy(cif_path, os.path.join(pdb_out_dir, f"{pdb_id}.cif"))

        # Save sequences to a FASTA file
        fasta_path = os.path.join(pdb_out_dir, f"{pdb_id}.fasta")
        records = [SeqRecord(Seq(seq), id=f"{pdb_id}_{chain_id}", description="") for chain_id, seq in chain_seqs.items()]
        SeqIO.write(records, fasta_path, "fasta")

        # Collect statistics
        total_lengths.append(total_len)
        pdb_count += 1
        
    except Exception as e:
        print(f"Error processing {cif_file}: {e}")

print(f"Total PDBs processed: {pdb_count}")
if total_lengths:
    stats = (
        f"Total sequence length stats: "
        f"min={min(total_lengths)}, "
        f"max={max(total_lengths)}, "
        f"mean={sum(total_lengths)/len(total_lengths):.2f}"
    )
    print(stats)

# Plotting the distribution of total sequence lengths
plt.figure(figsize=(6, 5))
plt.hist(total_lengths, bins=50, color='skyblue', edgecolor='black')
plt.xlabel('Total Sequence Length per PDB')
plt.ylabel('Number of PDBs')
plt.title('Distribution of Total Sequence Lengths')
plt.tight_layout()
plt.savefig('./total_sequence_length_distribution.png', dpi=300)
plt.show()
