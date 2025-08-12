import os
from Bio import SeqIO
from tqdm import tqdm

def read_msa_file(msa_path):
    """
    Reads all .a3m files in the directory of msa_file,
    uses the first sequence (query_0) as key, msa_file path as value.
    Returns a dict: {sequence: msa_file_path}
    """
    msa_dict = {}
    msa_files = [f for f in os.listdir(msa_path) if f.endswith(".a3m")]
    for file in msa_files:
        with open(os.path.join(msa_path, file), "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if records:
                query_seq = str(records[0].seq)
                msa_dict[query_seq] = os.path.join(msa_path, file)
    return msa_dict

def main(input_dir):
    # Read all folders in the input directory
    folders = [f for f in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, f))]
    for folder in tqdm(folders):
        target_id = folder
        # Read the sequences from the FASTA files in the folder
        fasta_file = os.path.join(input_dir, folder, target_id+".fasta")
        with open(fasta_file, "r") as f:
            sequences = list(SeqIO.parse(f, "fasta"))
        # Read MSA file
        msa_dict = read_msa_file(os.path.join(input_dir, folder, "msa", target_id, "msa_resmsa_seq_0"))
        # Construct the input FASTA file
        input_dir_path = os.path.join(input_dir, folder, "input")
        os.makedirs(input_dir_path, exist_ok=True)
        input_fasta_path = os.path.join(input_dir_path, f"{target_id}.fasta")
        with open(input_fasta_path, "w") as out_f:
            chain_id = ord('A')
            for seq_record in sequences:
                seq_str = str(seq_record.seq)
                msa_path = msa_dict.get(seq_str)
                header = f">{chr(chain_id)}|protein|{msa_path}"
                out_f.write(f"{header}\n{seq_str}\n")
                chain_id += 1

if __name__ == "__main__":
    # Define the path to the input directory
    input_dir = "./CASP16"
    main(input_dir)
