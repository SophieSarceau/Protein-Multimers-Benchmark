import json
import os
import re
from pathlib import Path
from typing import List
from Bio import SeqIO
from tqdm import tqdm

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")  # keep uncommon letters if present

def chain_labels(n: int) -> List[str]:
    """Return n Excel-like labels: A..Z, AA..AZ, BA.. etc., but AF3 usually expects single letters.
    We’ll still support >26 sequences just in case."""
    labels = []
    i = 0
    while len(labels) < n:
        x = i
        s = ""
        while True:
            s = chr(ord('A') + (x % 26)) + s
            x = x // 26 - 1
            if x < 0:
                break
        labels.append(s)
        i += 1
    return labels

def fasta_to_af3_json(name: str, records: List[SeqIO.SeqRecord], model_seeds: List[int]):
    all_chains = chain_labels(len(records))

    # Group records by sequence
    seq_to_chains = {}
    for chain_id, rec in zip(all_chains, records):
        seq = str(rec.seq)
        if not seq:
            raise ValueError(f"{name}: empty sequence for record {rec.id}")
        if seq not in seq_to_chains:
            seq_to_chains[seq] = []
        seq_to_chains[seq].append(chain_id)

    sequences_block = []
    for seq, chain_ids in seq_to_chains.items():
        # AF3 allows str | list[str] for "id". Use a list if multiple chains, else a string.
        protein_id = chain_ids if len(chain_ids) > 1 else chain_ids[0]
        sequences_block.append({
            "protein": {
                "id": protein_id,
                "sequence": seq
            }
        })

    return {
        "name": name,
        "sequences": sequences_block,
        "modelSeeds": model_seeds,
        "dialect": "alphafold3",
        "version": 1
    }

def main(casp_folder: str, seeds: List[int], overwrite: bool):
    folders = [f for f in os.listdir(casp_folder) if os.path.isdir(os.path.join(casp_folder, f))]
    for folder in tqdm(folders, desc="Processing folders"):
        target_id = folder
        target_fasta = os.path.join(casp_folder, target_id, f"{target_id}.fasta")
        if not os.path.isfile(target_fasta):
            print(f"Skipping {target_id}: no FASTA file found at {target_fasta}")
            continue
        out_json = os.path.join(casp_folder, target_id, f"{target_id}_input.json")
        if os.path.isfile(out_json) and not overwrite:
            print(f"Skipping {target_id}: output JSON {out_json} already exists (use overwrite=True to replace)")
            continue
        records = list(SeqIO.parse(target_fasta, "fasta"))
        try:
            af3_json_data = fasta_to_af3_json(name=target_id, records=records, model_seeds=seeds)
            with open(out_json, "w") as f:
                json.dump(af3_json_data, f, indent=4)
        except ValueError as e:
            print(f"Error processing {target_id}: {e}")

if __name__ == "__main__":
    casp_folder = './CASP16'
    seeds = [1]
    overwrite = True
    main(casp_folder=casp_folder, seeds=seeds, overwrite=overwrite)
