## AlphaFold2-Multimer Evaluation

### 1. Install ParaFold
Please refer to https://github.com/hermannschwaerzlerUIBK/ParallelFold.

### 2. Run AlphaFold2-Multimer
You need `gcc-13.2.0` to support `GLIBCXX_3.4.29`.
Then run `python run_af2_multimer.py`

### 3. Evaluate the results
Run `python cal_tmscore.py` to calculate TM-score.
