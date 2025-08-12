import os
import sys
import pandas as pd
import argparse
from tqdm import tqdm
import re
import json

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
        seq_len_match = re.search(r'Length of Structure_1:\s*(\d+)\s*residues', content)
        if seq_len_match:
            seq_len = int(seq_len_match.group(1))
        rmsd_match = re.search(r'RMSD=\s*([\d.]+)', content)
        if rmsd_match:
            rmsd = float(rmsd_match.group(1))
        tmscore_match = re.search(r'TM-score=\s*([\d.]+)\s*\(normalized by length of Structure_2', content)
        if tmscore_match:
            tmscore = float(tmscore_match.group(1))
    except Exception as e:
        print(f"Error parsing {output_file}: {e}")
        return {"seq_len": None, "rmsd": None, "tmscore": None}
    return {"seq_len": seq_len, "rmsd": rmsd, "tmscore": tmscore}

if __name__ == "__main__":
    # python cal_tmscore.py --folder_path ./CASP16/
    folder_path = './CASP16/'

    results_list = []
    folders = read_folder(folder_path)
    for folder in tqdm(folders):
        try:
            target_id = folder
            output_path = os.path.join(folder_path, folder, "output", "boltz_results_input", "predictions", target_id)
            # protenix may not have MSA, so only one output
            row_data = {"target_id": target_id, "seq_len": None, "rmsd": None, "tmscore": None}
            if os.path.exists(output_path):
                gt_path = os.path.join(folder_path, folder, target_id + ".pdb")
                pred_path = os.path.join(output_path, target_id+"_model_0.cif")
                spu_path = os.path.join(output_path, "spu")
                tmscore_output = os.path.join(output_path, "tmscore.txt")
                # Run USAlign
                os.system(f"./US-align/USalign {pred_path} {gt_path} -mm 1 -o {spu_path} > {tmscore_output}")
                # Parse results
                results = parse_tmscore_output(tmscore_output)
                row_data["seq_len"] = results["seq_len"]
                row_data["rmsd"] = results["rmsd"]
                row_data["tmscore"] = results["tmscore"]
            results_list.append(row_data)
        except Exception as e:
            print(f"Error processing folder {folder}: {e}")
            continue
    results = pd.DataFrame(results_list)
    output_xlsx = os.path.join(folder_path, "tmscore_results.xlsx")
    results.to_excel(output_xlsx, index=False)
    print(f"Results saved to {output_xlsx}")
    print("\nSummary Statistics:")
    print(f"Total targets processed: {len(results)}")
    print(f"Results available: {results['tmscore'].notna().sum()}")
    if results['tmscore'].notna().any():
        print(f"Average TM-score: {results['tmscore'].mean():.4f}")
