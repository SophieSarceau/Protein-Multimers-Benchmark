import os
from tqdm import tqdm

def main(casp_root_path):
    # Read all folders in the input directory
    folders = [f for f in os.listdir(casp_root_path) if os.path.isdir(os.path.join(casp_root_path, f))]
    for folder in tqdm(folders):
        # Perform structure prediction for each folder
        target_id = folder
        input_path = os.path.join(casp_root_path, target_id, "input")
        output_path = os.path.join(casp_root_path, target_id, "output")
        # Run structure prediction
        try:
            os.system(f"boltz predict {input_path} --out_dir {output_path} --override")
        except Exception as e:
            print(f"Error occurred while predicting structure for {target_id}: {e}")

if __name__ == "__main__":
    casp_root_path = './CASP16'
    main(casp_root_path)
