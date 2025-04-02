# Task2: Using short regulatory sequence to predict cell type specific transcription factor bindings

# Briefing

Text:
We aim to test how well genomic language models predict transcription factor (TF) binding, defined by binary ChIP-seq data in tissue-specific and TF-specific context. Specifically, how the models' sequence features generalise to unseen tissues and unseen TFs. Input: 
1. DNA sequence (from reference genome) under ChIP peaks.
2. Tissue-specific features obtained from bulk RNA seq.
3. TF-specific features (e.g., from motifs).

The task is divided into 4 smaller tasks. 
- a. For a fixed TF, predict in unseen tissues. 
- b. For a fixed tissue, predict in unseen TFs.
- c. When both the TF & tissue are in the training set, but the specific combination is unseen. 
- d. When the TF & tissue are both unseen. 

## Overview of the data & cartoon for the sub-tasks.
* Number of ChIP samples. 
![Alt text](./data/figs/ct_tf_counts.svg?raw=true "Distribution of ChIP samples compiled from ENCODE")
* Cartoon for the 4 tasks
![Alt text](./data/figs/cartoon.png?raw=true "task2") 



Features:
- Sequence dimension: relative small - mostly around 1kb.
- Sample size: number of the binding region is large. Numbers of tissues & TFs are small - 19 tissues and 312 TFs.
- Output dimension: Depends on the task. The full dimension is (n_sequence, n_tissue, n_TF).
- Prior knowedge on output dimension: TF binding motif, TF RNA expression.


# Data
- Bulk rna-seq. File: ad_rna_bulk.h5ad. Entry: tpm
- Chip-seq. 962 samples in 19 tissues. File: ad_chip_chromsubset.h5ad. Entry: binary
- (Optional) TF position weight matrices from JASPAR. File: 20250327110737_JASPAR2024_combined_matrices_532533_meme.txt. 
- (For writing fasta) Reference genome. File: hg38.fa 

## Preprocessing details & intermediate data
- Raw ChIP data (ad_chip.h5ad) is filtered for shared tissues with RNA data (ad_chip_filtered.h5ad). See preprocess.ipynb.
- ad_chip_filtered.h5ad is further subsetted (ad_chip_chromsubset.h5ad). See pp_subset.ipynb. 
    - ChIP regions are subsetted to chr19-22.
    - TFs without motif annotation in JASPAR are removed. 


# Demo data load
* Specify the data scope & splitting in config.json (Example: ./tasks/taska/trial1/config.json).
    * "all_ct"/"all_tf": Subset data by celltype/TF. Leave blank for all. 

```python
from helpers import * 

DATA_DIR = "/mnt/project/task2_Seq2CellxTF"

data = DataLoader(
    trial_path="./tasks/taska/trial1",
    rna_path = f"{DATA_DIR}/data/ad_rna_bulk.h5ad",
    chip_path = f"{DATA_DIR}/data/ad_chip_chromsubset.h5ad",
    tf_path = f"{DATA_DIR}/genomes/JASPAR_human_TFs_meme/20250327110737_JASPAR2024_combined_matrices_532533_meme.txt",
    rna_func = None,
    tf_func = None,
    fasta_ref = f"{DATA_DIR}/genomes/hg38/hg38.fa",
    )

ds = data.read_ds(key="train")
print(ds['Y_chip'].shape,       # (n_sample, n_region)
      len(ds['X_seq']),         # (n_region,)
      ds['X_rna'].shape,        # (n_sample, n_feature)
      len(ds['X_tf']))          # (n_sample)
                                # n_sample = n_TF * n_tissues 
```
* Region & sample indices are stored in ```data._idx_seqs``` & ```data._keys_cttf```. 

# Evaluation 
Submit as csv files (```columns=data._idx_seqs['test']```, ```index=data._keys_cttf['test']```) [here](http://10.64.155.14:5013/). 



# Env
- See task1_SNP2GEX for creating environments on the CPOS server. 
#### Dependencies for helpers.py 
```
conda install bioconda::bedtools -y
conda install --channel conda-forge --channel bioconda pybedtools

pip3 install pyfaidx scanpy 
```