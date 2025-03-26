# Task2: Using short regulatory sequence to predict cell type specific transcription factor bindings

# Briefing

Text:

Features:
- Input dimension: relative small - mostly around 1kb (51bp-44793bp).
- Sample size: number of the binding region is large
- Ourput dimension: huge - over 2000 (and can be more)
- Prior knowedge on output dimension: TF binding motif, TF RNA expression.


# Data
- Bulk RNA seq (ad_rna_bulk.h5ad) encompassing 600 samples in 111 tissues. Entry: tpm
- ChIP (ad_chip_filtered.h5ad) encompassing 1701 samples in 111 tissues. Entry: binary
- ChIP sequences (chip_seqs.fasta)

- See preprocess.ipynb for preprocessing and other details.

# Demo data load
```
import scanpy as sc 
data_dir = "/ssd/users/cfx/DNA-AI-hackathon/task2_Seq2CellxTF/data"
rna_bulk = sc.read_h5ad(f"{data_dir}/ad_rna_bulk.h5ad")
chip_bulk = sc.read_h5ad(f"{data_dir}/ad_chip_filtered.h5ad")
```

