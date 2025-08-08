import os
from tqdm import tqdm

def process_msa(input_json, out_dir):
    """
    Process the MSA using Protenix.
    
    Args:
        input_json (str): Path to the input JSON file.
        out_dir (str): Directory to save the output.
    """
    os.system(f"protenix msa --input {input_json} --out_dir {out_dir}")

def main(casp_root_dir):
    # Get the list of folders
    casp_folders = [f for f in os.listdir(casp_root_dir) if os.path.isdir(os.path.join(casp_root_dir, f))]
    for casp_folder in tqdm(casp_folders):
        input_json = os.path.join(casp_root_dir, casp_folder, f"{casp_folder}.json")
        out_dir = os.path.join(casp_root_dir, casp_folder, "msa")
        process_msa(input_json, out_dir)

if __name__ == "__main__":
    casp_root_dir = "./CASP16"
    main(casp_root_dir)
