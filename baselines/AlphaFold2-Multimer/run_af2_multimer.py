import os
import sys
from tqdm import tqdm

def get_folders(path):
    return [f.path for f in os.scandir(path) if f.is_dir()]

def main(casp_folder, data_dir):
    folders = get_folders(casp_folder)
    original_cwd = os.getcwd()
    # Define the directory where run_alphafold.sh is located
    script_dir = os.path.join(original_cwd, 'ParallelFold')

    for folder in tqdm(folders):
        target_id = os.path.basename(folder)

        # Create absolute paths for directories and files
        output_dir = os.path.join(original_cwd, casp_folder, target_id, 'af2_pred')
        fasta_path = os.path.join(original_cwd, casp_folder, target_id, f'{target_id}.fasta')
        abs_data_dir = os.path.join(original_cwd, data_dir)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # Check if there is already a prediction
        pred_file = os.path.join(output_dir, target_id, 'ranking_debug.json')
        if os.path.exists(pred_file):
            print(f'Skipping {target_id} as prediction already exists.')
            continue
        print(f'Processing {target_id}...')
        # Change into the script's directory before running it
        os.chdir(script_dir)

        # Use absolute paths for the script's arguments
        cmd = f'./run_alphafold.sh \
                -d {abs_data_dir} \
                -o {output_dir} \
                -p multimer \
                -i {fasta_path} \
                -t 2021-11-01'
        os.system(cmd)

        # Change back to the original directory
        os.chdir(original_cwd)

if __name__ == "__main__":
    # casp_folder = "./CASP16"
    casp_folder = "./smp1"
    data_dir = "./ParallelFold/data"
    main(casp_folder, data_dir)
