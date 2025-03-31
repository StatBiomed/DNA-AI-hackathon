# Task1: Using genetic variants to predict personalised gene expression
## Task 1a: for gene level expression
### Object
Develop a sequence-based model to predict personalized gene expression by utilizing the DNA sequences. Your model should capture how DNA sequence variations influence gene expression levels.
### Evaluation Metrics
Model performance will be assessed through three key dimensions:

* **Across-genes**: Pearson correlation between true and predicted genome-wide expression profiles ​across held-out genes per test individual.

* **Across-individuals per seen gene**: Per-gene correlation between predictions and ground truth ​across held-out individuals for genes seen during training.

* **Across-individuals per unseen gene**: Per-gene correlation between predictions and ground truth ​across held-out individuals for unseen genes from training.

![Fig.1](imgs/Fig1.jpg)

### Data
Root directory: `/home/smli/ssd/miniHackathon/` 

The dataset contains one subfoler and one CSV file:
* **fasta/** directory: 
    - Contains ​two FASTA files per individual (named <sampleID>_allele_[1/2].fasta)
    - Format per file:
        * Each record contains ​32kb flanking sequence around a gene's transcription start site (TSS)
        * ​Header syntax:
        ```bash
        >ENSG00000006282 chr:17;strand:+;start:48604419;end:48636419
        GGAGGGGGGCCTACTCTCTGACCCTGGCAAATCTTGGAGAAGGGGTTCATAGGTACAGAT
        TTCTGAGGGGGGTCCCTGGCTCCCACCAAAGGCACCCAGACAGCTCTCCATAGCTGCATC
        CCCTCCTGGTTCCTGGTCCCCTGCCACCCATCCCCACATCACCATGCCCTTCACTAGAG...
        ```
        Annotations include chromosome, strand, genomic coordinates.
* **partitions.csv**:
    - Contains gene expression levels and dataset splits:
    
        | Column name | Description |
        | --- | --- |
        | `sample` | individual identifier | 
        | `gene` | ENSG gene IDs |
        | `symbol` | gene symbols |
        | `chr` | chromosome |
        | `TSS` | transcriptional start site |
        | `strand` | DNA strand of the gene |
        | `rpkm`, `log_TPM`, `raw count`, `log_count` | gene expression levels |
        | `labels` | Group identifiers (`seen gene seen individuals`/`seen gene unseen individuals`/`unseen gene unseen individuals`) | 
        | `split` | Dataset partition (`train`/`test`) |


### Demo data load
We provide a dataloader for loading the raw sequences and the truth labels from the fasta file and the partition CSV file. The example code shown below:

```python
import Dataset
from torch.utils.data import DataLoader
import os
DATA_ROOT = "/home/smli/ssd/miniHackathon/"

# set the `split` as "test" when you are evaluating.
# sample_list = [], or gene_list = [] means all individuals and all genes of the specific split category will be used. if you want to select subset of them, pass the list to the corresponding argument.
# you can also set target_name="raw count"/"log_count"/"rpkm" as you want

trainDataset = Dataset.SampleGeneExpressionDataset(split="train",csv_file=os.path.join(DATA_ROOT,"partitions.csv"),consensus_root=os.path.join(DATA_ROOT,"fasta"),sample_list=[],gene_list=[],target_name="log_TPM")
# set the batch size or the sampler as you want.
trainDataLoader = DataLoader(trainDataset, batch_size=1,shuffle=True)
# iterate the dataset, modify it as needed to utilize gLM embeddings!
for cur_seq_1, cur_seq_2, sample_target,sample_name,gene_name in trainDataLoader:
    print(f"cur_seq_1: {cur_seq_1}")
    print(f"cur_seq_2: {cur_seq_2}")
    print(f"sample_target:{sample_target}")
    print(f"sample_name:{sample_name}")
    print(f"gene_name: {gene_name}")
    break
```
Alternatively, you can also implement your own dataloader.

### Evaluation

You can upload your predictions and perform evaluation at our [ranking borad](http://10.64.155.14:5011) 

* **Usage**: make sure your uploaded CSV file at least has three columns: `sample`, `gene`, `log_TPM`, where `log_TPM` stores the predictions.

## Task 1b: for promoter-level activities
# Data

# Demo data load

