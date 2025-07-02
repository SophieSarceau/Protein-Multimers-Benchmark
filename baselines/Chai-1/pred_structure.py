import os
import sys
import argparse
from tqdm import tqdm


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


if __name__ == "__main__":
    # python pred_structure.py --folder_path ./benchmark/CASP15/
    parser = argparse.ArgumentParser(description="Read folder names from a specified directory.")
    parser.add_argument("--folder_path", type=str, required=True, help="Path to the directory containing folders.")

    args = parser.parse_args()

    folders = read_folder(args.folder_path)
    for folder in tqdm(folders):
        target_id = folder
        fasta_path = os.path.join(args.folder_path, folder, target_id + ".fasta")
        msa_path = os.path.join(args.folder_path, folder, "msa")
        # For non-msa
        no_msa_output_path = os.path.join(args.folder_path, folder, "output")
        msa_output_path = os.path.join(args.folder_path, folder, "output_msa")
        if not os.path.exists(no_msa_output_path):
            os.makedirs(no_msa_output_path)
        if not os.path.exists(msa_output_path):
            os.makedirs(msa_output_path)
        print(f"Processing folder: {folder}")
        # Run the command
        os.system(f"chai-lab fold {fasta_path} {no_msa_output_path}")
        os.system(f"chai-lab fold {fasta_path} {msa_output_path} --msa-directory {msa_path}")
        print(f"Processed {folder} successfully.")
