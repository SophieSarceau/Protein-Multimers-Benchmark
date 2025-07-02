import os
import sys
import pandas as pd
import argparse
from tqdm import tqdm
import re


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


def parse_tmscore_output(output_file):
    """
    Parse USAlign output file to extract sequence length, RMSD, and TM-score.
    
    Args:
        output_file (str): Path to the USAlign output file
    
    Returns:
        dict: Dictionary containing seq_len, rmsd, and tmscore
    """
    if not os.path.exists(output_file):
        return {"seq_len": None, "rmsd": None, "tmscore": None}
    
    seq_len = None
    rmsd = None
    tmscore = None
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Extract Length of Structure_1 (sequence length)
        seq_len_match = re.search(r'Length of Structure_1:\s*(\d+)\s*residues', content)
        if seq_len_match:
            seq_len = int(seq_len_match.group(1))
        
        # Extract RMSD
        rmsd_match = re.search(r'RMSD=\s*([\d.]+)', content)
        if rmsd_match:
            rmsd = float(rmsd_match.group(1))
        
        # Extract TM-score (normalized by length of Structure_2)
        tmscore_match = re.search(r'TM-score=\s*([\d.]+)\s*\(normalized by length of Structure_2', content)
        if tmscore_match:
            tmscore = float(tmscore_match.group(1))
    
    except Exception as e:
        print(f"Error parsing {output_file}: {e}")
        return {"seq_len": None, "rmsd": None, "tmscore": None}
    
    return {"seq_len": seq_len, "rmsd": rmsd, "tmscore": tmscore}


if __name__ == "__main__":
    # python cal_tmscore.py --folder_path ./benchmark/CASP15/
    parser = argparse.ArgumentParser(description="Calculate TM-scores for predicted structures.")
    parser.add_argument("--folder_path", type=str, required=True, help="Path to the directory containing folders.")

    args = parser.parse_args()

    # Create a list to store results instead of DataFrame
    results_list = []

    folders = read_folder(args.folder_path)
    for folder in tqdm(folders):
        target_id = folder
        output_path = os.path.join(args.folder_path, folder, "output")
        msa_output_path = os.path.join(args.folder_path, folder, "output_msa")

        # Initialize row data
        row_data = {"target_id": target_id, "seq_len": None, 
                   "no_msa_rmsd": None, "no_msa_tmscore": None,
                   "msa_rmsd": None, "msa_tmscore": None}

        # Process no-MSA results
        if os.path.exists(output_path):
            gt_path = os.path.join(args.folder_path, folder, target_id + ".pdb")
            pred_path = os.path.join(output_path, "pred.model_idx_0.cif")
            spu_path = os.path.join(output_path, "spu")
            tmscore_output = os.path.join(output_path, "tmscore.txt")

            # Run USAlign
            os.system(f"./US-align/USalign {pred_path} {gt_path} -mm 1 -o {spu_path} > {tmscore_output}")

            # Parse results
            no_msa_results = parse_tmscore_output(tmscore_output)
            row_data["seq_len"] = no_msa_results["seq_len"]
            row_data["no_msa_rmsd"] = no_msa_results["rmsd"]
            row_data["no_msa_tmscore"] = no_msa_results["tmscore"]

        # Process MSA results
        if os.path.exists(msa_output_path):
            gt_path = os.path.join(args.folder_path, folder, target_id + ".pdb")
            pred_path = os.path.join(msa_output_path, "pred.model_idx_0.cif")
            spu_path = os.path.join(msa_output_path, "spu")
            tmscore_output = os.path.join(msa_output_path, "tmscore.txt")

            # Run USAlign
            os.system(f"./US-align/USalign {pred_path} {gt_path} -mm 1 -o {spu_path} > {tmscore_output}")

            # Parse results
            msa_results = parse_tmscore_output(tmscore_output)
            if row_data["seq_len"] is None:  # If we didn't get seq_len from no_msa
                row_data["seq_len"] = msa_results["seq_len"]
            row_data["msa_rmsd"] = msa_results["rmsd"]
            row_data["msa_tmscore"] = msa_results["tmscore"]

        # Add row to results list
        results_list.append(row_data)

    # Create DataFrame from list of dictionaries
    results = pd.DataFrame(results_list)

    # Save results to xlsx file
    output_xlsx = os.path.join(args.folder_path, "tmscore_results.xlsx")
    results.to_excel(output_xlsx, index=False)
    print(f"Results saved to {output_xlsx}")

    # Display summary statistics
    print("\nSummary Statistics:")
    print(f"Total targets processed: {len(results)}")
    print(f"No-MSA results available: {results['no_msa_tmscore'].notna().sum()}")
    print(f"MSA results available: {results['msa_tmscore'].notna().sum()}")

    if results['no_msa_tmscore'].notna().any():
        print(f"Average No-MSA TM-score: {results['no_msa_tmscore'].mean():.4f}")
    if results['msa_tmscore'].notna().any():
        print(f"Average MSA TM-score: {results['msa_tmscore'].mean():.4f}")
