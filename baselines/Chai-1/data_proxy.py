import os
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict


def read_folder(folder_path):
    """
    Reads all folder names in the specified directory.
    :param folder_path: Path to the directory containing folders.
    :return: List of folder names.
    """
    if not os.path.isdir(folder_path):
        raise ValueError(f"The path {folder_path} is not a valid directory.")

    folder_names = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]

    return folder_names


def parse_fasta_with_biopython(fasta_path):
    """
    Parses a FASTA file using BioPython and extracts chain sequences.
    :param fasta_path: Path to the FASTA file.
    :return: Dictionary with chain IDs as keys and Bio.Seq sequences as values.
    """
    seqs = {}

    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Extract chain ID (use the full record.id or just the first part)
            seq = record.seq
            if seq not in seqs.keys():
                seqs[seq] = 1
            else:
                seqs[seq] += 1
    except Exception as e:
        print(f"Error reading {fasta_path}: {e}")
        return {}

    return seqs


def analyze_sequence_composition(seqs):
    """
    Analyzes the composition of sequences.
    :param seqs: Dictionary with sequences as keys and their counts as values.
    :return: Dictionary with sequence composition statistics.
    """
    if not seqs:
        return {
            "total_chains": 0,
            "composition": "empty",
            "chain_formula": "empty"
        }

    total_chains = sum(seqs.values())
    unique_sequences = len(seqs)

    # Sort counts in descending order for consistent naming
    counts = sorted(seqs.values(), reverse=True)

    # Determine composition type
    if unique_sequences == 1:
        # All chains have the same sequence (homo-oligomer)
        composition = "homo"
        chain_formula = f"A{total_chains}"
    else:
        # Multiple different sequences (hetero-oligomer)
        composition = "hetero"

        # Generate chain formula like A2B1C1
        chain_labels = []
        for i, count in enumerate(counts):
            chain_letter = chr(ord('A') + i)  # A, B, C, D, ...
            chain_labels.append(f"{chain_letter}{count}")

        chain_formula = "".join(chain_labels)

    # Additional analysis
    analysis = {
        "total_chains": total_chains,
        "composition": composition,
        "chain_formula": chain_formula,
        "sequence_counts": counts,
    }

    return analysis


def rewrite_fasta_with_chain_labels(fasta_path, target_id, seqs, output_path=None):
    """
    Rewrites a FASTA file with standardized chain labels.
    :param fasta_path: Path to the input FASTA file.
    :param target_id: Target ID to use in the new headers.
    :param seqs: Dictionary with sequences as keys and their counts as values.
    :param analysis: Analysis results from analyze_sequence_composition.
    :param output_path: Path for the output file. If None, overwrites the input file.
    :return: Dictionary mapping old headers to new headers.
    """
    if output_path is None:
        output_path = fasta_path

    # Sort sequences by count (descending) to assign chain letters consistently
    sorted_seqs = sorted(seqs.items(), key=lambda x: x[1], reverse=True)
    seq_to_chain = {}

    # Assign chain letters based on the analysis
    for i, (seq, count) in enumerate(sorted_seqs):
        chain_letter = chr(ord('A') + i)
        seq_to_chain[str(seq)] = chain_letter

    header_mapping = {}
    new_records = []
    chain_counters = {}

    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq)
            chain_letter = seq_to_chain[seq]

            # Track count for this chain type
            if chain_letter not in chain_counters:
                chain_counters[chain_letter] = 0
            chain_counters[chain_letter] += 1

            # Create new header
            new_header = f"protein|{target_id}-{chain_letter}{chain_counters[chain_letter]}"
            header_mapping[record.id] = new_header

            # Create new record
            record.id = new_header
            record.description = ""
            new_records.append(record)

        # Write the new FASTA file with sequences in one line
        with open(output_path, 'w') as output_file:
            for record in new_records:
                output_file.write(f">{record.id}\n{str(record.seq)}\n")
        print(f"Rewritten FASTA saved to: {output_path}")

    except Exception as e:
        print(f"Error processing {fasta_path}: {e}")
        return {}

    return header_mapping


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read folder names from a specified directory and analyze FASTA files.")
    parser.add_argument("--folder_path", type=str, help="Path to the directory containing folders.")

    args = parser.parse_args()

    folders = read_folder(args.folder_path)
    for folder in folders:
        print(f"Processing folder: {folder}")
        target_id = folder
        fasta_path = os.path.join(args.folder_path, folder, target_id+".fasta")
        seqs = parse_fasta_with_biopython(fasta_path)
        analysis = analyze_sequence_composition(seqs)
        print(f"Analysis for {target_id}:")
        print(f"Total chains: {analysis['total_chains']}")
        print(f"Composition: {analysis['composition']}")
        print(f"Chain formula: {analysis['chain_formula']}")
        print(f"Sequence counts: {analysis['sequence_counts']}")
        output_path = os.path.join(args.folder_path, folder, target_id + ".fasta")
        header_mapping = rewrite_fasta_with_chain_labels(fasta_path, target_id, seqs, output_path=output_path)
        print(f"Header mapping: {header_mapping}")
        print("-" * 40)
