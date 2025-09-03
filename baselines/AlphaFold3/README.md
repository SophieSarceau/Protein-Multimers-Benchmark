## AlphaFold3 Evaluation

### 1. Construct Input JSON
Run `python generate_input_json.py`.

### 2. Run AlphaFold3 (Online Server)
#### 2.1 Job Submission
Submit the generated input json files to AlphaFold3 server.

#### 2.2 TM-Score Calculation
Run `python cal_server_tmscore.py`.

### 3. Run AlphaFold3 (Local)
#### 3.1 Clone AF3
```bash
git clone https://github.com/google-deepmind/alphafold3.git
cd alphafold3
conda create --name af3 python=3.11
conda activate af3
pip install -r requirements.txt
pip install -e .
```

#### 3.2 Build Data
```bash
cd src/alphafold3
python build_data.py
```

#### 3.3 Download Hmmer
```bash
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar -xvf hmmer.tar.gz
cd hmmer-*/
./configure --prefix=/path/to/your/install/dir
make
make install
```

#### 3.4 Construct Input JSON
Run `python gen_local_input_json.py`.

#### 3.5 Run AF3
Run `python run_af3.py`.

#### 3.6 TM-Score Calculation
Run `python cal_local_tmscore.py`.
