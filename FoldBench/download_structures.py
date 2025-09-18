import os
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial
import warnings
from Bio.PDB import MMCIFParser, MMCIFIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning

def read_pdb_ids(csv_file_path):
    df = pd.read_csv(csv_file_path)
    pdb_ids = df['pdb_id'].tolist()
    # Remove -assembly1 suffix if present
    pdb_ids = [pdb_id.replace('-assembly1', '') for pdb_id in pdb_ids]
    # Remove duplicates
    pdb_ids = list(set(pdb_ids))

    return pdb_ids

def download_and_unzip(pdb_id, download_dir):
    """Downloads, unzips, and cleans a single CIF file."""
    base_url = "https://files.rcsb.org/download/{}-assembly1.cif.gz"
    url = base_url.format(pdb_id)
    output_path_gz = os.path.join(download_dir, f"{pdb_id}.cif.gz")
    output_path_cif = os.path.join(download_dir, f"{pdb_id}.cif")
    
    # Use wget to download. -q for quiet, -O to specify output file.
    os.system(f"wget -q {url} -O {output_path_gz}")
    
    # Unzip the file if download was successful
    if os.path.exists(output_path_gz):
        os.system(f"gunzip -f {output_path_gz}")

    # Clean the CIF file using BioPython if it exists
    if os.path.exists(output_path_cif):
        try:
            # Suppress warnings about discontinuous chains, etc.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', PDBConstructionWarning)
                parser = MMCIFParser()
                structure = parser.get_structure(pdb_id, output_path_cif)

            # Identify and remove ligands, water, and ions (HETATMs)
            for model in structure:
                chains_to_remove = []
                for chain in model:
                    residues_to_remove = []
                    for residue in chain:
                        # The first element of residue.id is the hetero-field.
                        # ' ' for standard residues (amino acids/nucleotides)
                        # 'W' for water
                        # 'H_' + name for other hetero-atoms
                        if residue.id[0] != ' ':
                            residues_to_remove.append(residue.id)
                    
                    for res_id in residues_to_remove:
                        chain.detach_child(res_id)
                    
                    # If a chain becomes empty after removing residues, mark it for removal
                    if not list(chain.get_residues()):
                        chains_to_remove.append(chain.id)

                for chain_id in chains_to_remove:
                    model.detach_child(chain_id)

            # Save the cleaned structure back to the CIF file
            io = MMCIFIO()
            io.set_structure(structure)
            io.save(output_path_cif)

        except Exception as e:
            print(f"Could not process file {pdb_id}.cif. Error: {e}")

def download_cif_files(pdb_ids, download_dir='./cif_files'):
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    # Use multiprocessing Pool to download files in parallel
    num_processes = min(cpu_count(), 8)  # Limit to 8 processes to avoid overloading the server

    # Create a partial function to pass the download_dir to the worker
    worker_func = partial(download_and_unzip, download_dir=download_dir)

    with Pool(processes=num_processes) as pool:
        # Use tqdm to display a progress bar
        list(tqdm(pool.imap_unordered(worker_func, pdb_ids), total=len(pdb_ids), desc="Downloading CIF files"))

if __name__ == "__main__":
    csv_file_path = './interface_protein_protein.csv'
    pdb_ids = read_pdb_ids(csv_file_path)
    print(f"Total PDB IDs to download: {len(pdb_ids)}")
    download_cif_files(pdb_ids, download_dir='./cif_files')
