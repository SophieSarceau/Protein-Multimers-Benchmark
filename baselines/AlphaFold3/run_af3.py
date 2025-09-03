import os
import sys
from tqdm import tqdm

# Directory for model weights
MODEL_DIR = './alphafold3/model_weights'

# Directory for databases
DB_DIR = '/mnt/dna01/library-seq/alphafold3-database'

# Paths to database files
UNIREF90_DATABASE_PATH = os.path.join(DB_DIR, 'uniref90_2022_05.fa')
MGNIFY_DATABASE_PATH = os.path.join(DB_DIR, 'mgy_clusters_2022_05.fa')
PDB_DATABASE_PATH = os.path.join(DB_DIR, 'mmcif_files')
SMALL_BFD_DATABASE_PATH = os.path.join(DB_DIR, 'bfd-first_non_consensus_sequences.fasta')
UNIPROT_CLUSTER_ANNOT_DATABASE_PATH = os.path.join(DB_DIR, 'uniprot_all_2021_04.fa')
SEQRES_DATABASE_PATH = os.path.join(DB_DIR, 'pdb_seqres_2022_09_28.fasta')

# Paths to HMMER binaries
JACKHMMER_BINARY_PATH = './alphafold3/hmmer/bin/jackhmmer'
NHMMER_BINARY_PATH = './alphafold3/hmmer/bin/nhmmer'
HMMALIGN_BINARY_PATH = './alphafold3/hmmer/bin/hmmalign'
HMMSEARCH_BINARY_PATH = './alphafold3/hmmer/bin/hmmsearch'
HMMBUILD_BINARY_PATH = './alphafold3/hmmer/bin/hmmbuild'

def construct_af3_command(input_json: str, output_dir: str, model_dir: str, db_dir: str):
    cmd = (
        f"python ./alphafold3/run_alphafold.py "
        f"--json_path={input_json} "
        f"--output_dir={output_dir} "
        f"--model_dir={model_dir} "
        f"--db_dir={db_dir} "
        f"--uniref90_database_path={UNIREF90_DATABASE_PATH} "
        f"--mgnify_database_path={MGNIFY_DATABASE_PATH} "
        f"--pdb_database_path={PDB_DATABASE_PATH} "
        f"--small_bfd_database_path={SMALL_BFD_DATABASE_PATH} "
        f"--uniprot_cluster_annot_database_path={UNIPROT_CLUSTER_ANNOT_DATABASE_PATH} "
        f"--seqres_database_path={SEQRES_DATABASE_PATH} "
        f"--jackhmmer_binary_path={JACKHMMER_BINARY_PATH} "
        f"--nhmmer_binary_path={NHMMER_BINARY_PATH} "
        f"--hmmalign_binary_path={HMMALIGN_BINARY_PATH} "
        f"--hmmsearch_binary_path={HMMSEARCH_BINARY_PATH} "
        f"--hmmbuild_binary_path={HMMBUILD_BINARY_PATH} "
    )
    return cmd

def main(casp_folder: str, out_folder: str):
    os.makedirs(out_folder, exist_ok=True)
    folders = [f for f in os.listdir(casp_folder) if os.path.isdir(os.path.join(casp_folder, f))]
    for folder in tqdm(folders):
        target_id = folder
        target_json_input = os.path.join(casp_folder, target_id, f"{target_id}_input.json")
        out_folder_target = os.path.join(out_folder, target_id)
        os.makedirs(out_folder_target, exist_ok=True)
        if not os.path.isfile(target_json_input):
            print(f"Skipping {target_id}: no input JSON file found at {target_json_input}")
            continue
        cmd = construct_af3_command(
            input_json=target_json_input,
            output_dir=out_folder_target,
            model_dir=MODEL_DIR,
            db_dir=DB_DIR
        )
        print(f"Running AF3 for {target_id}...")
        ret = os.system(cmd)
        if ret != 0:
            print(f"AF3 failed for {target_id} with return code {ret}")
        else:
            print(f"AF3 completed successfully for {target_id}")
        print("-" * 40)

if __name__ == "__main__":
    casp_folder = './CASP16'
    out_folder = './local_pred2'
    main(casp_folder=casp_folder, out_folder=out_folder)
