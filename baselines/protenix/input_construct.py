import os
import json
from collections import Counter
from Bio import SeqIO
from tqdm import tqdm

def read_fasta(filepath):
    """Reads a FASTA file and returns a list of sequences (including duplicates)."""
    return [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]

def construct_json_per_fasta(casp_root_dir):
    # Go through each CASP subfolder (e.g., CASP15, CASP16)
    for casp_subdir in tqdm(os.listdir(casp_root_dir)):
        # Only select folder instead of files
        if not os.path.isdir(os.path.join(casp_root_dir, casp_subdir)):
            continue
        casp_subdir_path = os.path.join(casp_root_dir, casp_subdir)
        target_id = casp_subdir
        fasta_path = os.path.join(casp_subdir_path, f"{target_id}.fasta")
        seqs = read_fasta(fasta_path)
        if seqs:
            seq_counter = Counter(seqs)
            protein_entry = {
                "name": target_id,
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": seq,
                            "count": count
                        }
                    } for seq, count in seq_counter.items()
                ]
            }
            # Save JSON next to the fasta file
            json_filename = target_id + ".json"
            json_path = os.path.join(casp_subdir_path, json_filename)
            with open(json_path, 'w') as out_f:
                json.dump([protein_entry], out_f, indent=4)

if __name__ == "__main__":
    casp_root_dir = "./CASP16"  # Or "./CASP" for all CASP folders
    construct_json_per_fasta(casp_root_dir)
