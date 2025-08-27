import os
import json
from tqdm import tqdm
from collections import Counter
from Bio import SeqIO

CASP_ID = 16
CASP_DIR = f"./CASP{CASP_ID}"

def parse_fasta_biopython(fasta_path):
    # Use Biopython to parse FASTA and return a list of sequences (as strings)
    return [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]

def main():
    all_json_objs = []
    for target in tqdm(os.listdir(CASP_DIR)):
        target_dir = os.path.join(CASP_DIR, target)
        if not os.path.isdir(target_dir):
            continue
        fasta_file = os.path.join(target_dir, f"{target}.fasta")
        if not os.path.isfile(fasta_file):
            continue
        seqs = parse_fasta_biopython(fasta_file)
        seq_counts = Counter(seqs)
        json_obj = {
            "name": f"CASP{CASP_ID}_{target}",
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": seq,
                        "count": count
                    }
                } for seq, count in seq_counts.items()
            ],
            "dialect": "alphafoldserver",
            "version": 1
        }
        all_json_objs.append(json_obj)

    # Save combined JSON
    combined_json_path = f"./CASP{CASP_ID}/CASP{CASP_ID}_all_targets.json"
    with open(combined_json_path, "w") as out_f:
        json.dump(all_json_objs, out_f, indent=4)
    print(f"Wrote combined JSON to {combined_json_path}")

if __name__ == "__main__":
    main()
