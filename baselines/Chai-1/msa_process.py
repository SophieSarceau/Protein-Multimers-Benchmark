import os
import sys
import argparse
from Bio import SeqIO


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


def aln_to_a3m(aln_file, a3m_file, query_id="A"):
    """
    Convert alignment file (.aln) to A3M format (.a3m)
    
    Args:
        aln_file (str): Path to input alignment file
        a3m_file (str): Path to output A3M file
        query_id (str): ID for the query sequence (first sequence)
    """
    sequences = []
    
    # Read alignment file
    with open(aln_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line:  # Skip empty lines
                sequences.append(line)
    
    if not sequences:
        raise ValueError("No sequences found in alignment file")
    
    # Write A3M format
    with open(a3m_file, 'w') as f:
        for i, seq in enumerate(sequences):
            if i == 0:
                # First sequence is the query
                f.write(f">{query_id}\n")
            else:
                # Subsequent sequences get generic IDs
                f.write(f">seq_{i}\n")
            
            # Convert alignment format to A3M format
            # In A3M format, lowercase letters represent insertions relative to query
            # Uppercase letters represent match states
            # Gaps in query are removed, gaps in other sequences become lowercase
            
            if i == 0:
                # Query sequence: remove gaps, keep uppercase
                query_seq = seq.replace('-', '')
                f.write(f"{query_seq}\n")
            else:
                # Other sequences: convert based on query gaps
                a3m_seq = convert_to_a3m(sequences[0], seq)
                f.write(f"{a3m_seq}\n")


def convert_to_a3m(query_with_gaps, target_with_gaps):
    """
    Convert a target sequence to A3M format based on query sequence gaps
    
    Args:
        query_with_gaps (str): Query sequence with gaps
        target_with_gaps (str): Target sequence with gaps
    
    Returns:
        str: Target sequence in A3M format
    """
    if len(query_with_gaps) != len(target_with_gaps):
        raise ValueError("Query and target sequences must have same length")
    
    a3m_seq = ""
    
    for q_char, t_char in zip(query_with_gaps, target_with_gaps):
        if q_char == '-':
            # Gap in query - this position is an insertion in target
            if t_char != '-':
                a3m_seq += t_char.lower()
            # If both have gaps, skip this position
        else:
            # Match state position in query
            if t_char == '-':
                a3m_seq += '-'
            else:
                a3m_seq += t_char.upper()
    
    return a3m_seq


def add_metadata_to_a3m(a3m_file, metadata_list=None):
    """
    Add metadata headers to A3M sequences (optional enhancement)
    
    Args:
        a3m_file (str): Path to A3M file to modify
        metadata_list (list): List of metadata strings for each sequence
    """
    if metadata_list is None:
        return
    
    with open(a3m_file, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    seq_count = 0
    
    for line in lines:
        if line.startswith('>'):
            if seq_count < len(metadata_list) and metadata_list[seq_count]:
                # Replace header with metadata
                new_lines.append(f">{metadata_list[seq_count]}\n")
            else:
                new_lines.append(line)
            seq_count += 1
        else:
            new_lines.append(line)
    
    with open(a3m_file, 'w') as f:
        f.writelines(new_lines)


def get_best_msa_combinations(msa_sort_file, top_n=5):
    """
    Read the msa_sort_file and select the best N msa combinations based on scores
    
    Args:
        msa_sort_file (str): Path to the all.sort file containing MSA combinations and scores
        top_n (int): Number of top combinations to return (default: 5)
    
    Returns:
        list: List of tuples containing (combination_name, score) for the best combinations
    """
    if not os.path.exists(msa_sort_file):
        print(f"Warning: MSA sort file not found: {msa_sort_file}")
        return []

    combinations = []

    try:
        with open(msa_sort_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:  # Skip empty lines
                    parts = line.split()
                    if len(parts) >= 2:
                        combination_name = parts[0]
                        try:
                            score = float(parts[1])
                            combinations.append((combination_name, score))
                        except ValueError:
                            print(f"Warning: Invalid score format in line: {line}")
                            continue
    except Exception as e:
        print(f"Error reading MSA sort file {msa_sort_file}: {e}")
        return []

    # Sort by score in descending order (higher scores are better)
    combinations.sort(key=lambda x: x[1], reverse=True)
    combinations = combinations[:top_n]  # Limit to top N combinations
    # Return only the combination names, not scores
    combinations = [name for name, score in combinations]

    # Return top N combinations
    return combinations


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


# Example usage
if __name__ == "__main__":
    # python msa_process.py --folder_path ./benchmark/CASP15/ --msa_path /mnt/rna01/zhengxz/DMFold
    parser = argparse.ArgumentParser(description="Convert alignment file to A3M format")
    parser.add_argument("--folder_path", type=str, required=True, help="Path to folder containing alignment files")
    parser.add_argument("--msa_path", type=str, default="msa", help="Subdirectory for MSA files")
    parser.add_argument("--query_id", type=str, default="A", help="ID for the query sequence")

    args = parser.parse_args()

    folder_path = args.folder_path
    msa_path = args.msa_path

    folders = read_folder(folder_path)
    for folder in folders:
        print(f"Processing folder: {folder}")
        # Create a msa directory if it doesn't exist
        msa_dir = os.path.join(folder_path, folder, "msa")
        os.makedirs(msa_dir, exist_ok=True)
        target_id = folder
        fasta_file_path = os.path.join(folder_path, folder, f"{target_id}.fasta")
        seqs = parse_fasta_with_biopython(fasta_file_path)
        if len(seqs.keys()) >=2: # For Heteromer
            ori_aln_path = os.path.join(msa_path, folder, "AF2Models")
            msa_sort_file = os.path.join(ori_aln_path, "all.sort")
            if os.path.exists(msa_sort_file):
                # Read the msa_sort_file and select the best five msa combinations
                msa_combinations = get_best_msa_combinations(msa_sort_file, top_n=5)
                print(f"Best MSA combinations for {target_id}: {msa_combinations}")
            # Copy the aln file to the msa directory
            for msa_combination in msa_combinations:
                aln_file_path = os.path.join(ori_aln_path, msa_combination, "multi2.2-v4", "seq", "msas", "paired.aln")
                if os.path.exists(aln_file_path):
                    a3m_file_path = os.path.join(msa_dir, f"{msa_combination}.a3m")
                    aln_to_a3m(aln_file_path, a3m_file_path, query_id=args.query_id)
                else:
                    print(f"Warning: MSA file not found for combination {msa_combination} in {ori_aln_path}")
            # Use chai-lab a3m-to-pqt command to convert A3M files to PQT format
            os.system(f"chai-lab a3m-to-pqt {msa_dir}")
        else: # For Homomer
            ori_aln_path = os.path.join(msa_path, folder, "AF2Models")
            msa_file = os.path.join(ori_aln_path, "all.list")
            # Read all combinations from the all.list file
            if os.path.exists(msa_file):
                with open(msa_file, 'r') as f:
                    msa_combinations = [line.strip() for line in f if line.strip()]
                print(f"All MSA combinations for {target_id}: {msa_combinations}")
            for msa_combination in msa_combinations:
                a3m_file_path = os.path.join(ori_aln_path, msa_combination, "multi2.2-v4", "seq", "msas", "A", "alphafold2.a3m")
                if os.path.exists(a3m_file_path):
                    # Copy the A3M file to the msa directory
                    new_a3m_file_path = os.path.join(msa_dir, f"{msa_combination}.a3m")
                    # use cp command to copy the file
                    os.system(f"cp {a3m_file_path} {new_a3m_file_path}")
                else:
                    print(f"Warning: A3M file not found for combination {msa_combination} in {ori_aln_path}")
            # Use chai-lab a3m-to-pqt command to convert A3M files to PQT format
            os.system(f"chai-lab a3m-to-pqt {msa_dir}")

    print("MSA processing completed.")
    sys.exit(0)
