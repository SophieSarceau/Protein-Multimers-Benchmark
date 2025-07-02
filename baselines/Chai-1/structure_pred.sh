#!/bin/bash

python data_proxy.py --folder_path ./benchmark/CASP15/

python msa_process.py --folder_path ./benchmark/CASP15/ --msa_path ../DMFold

python pred_structure.py --folder_path ./benchmark/CASP15/

python cal_tmscore.py --folder_path ./benchmark/CASP15/
