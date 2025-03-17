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
We provide a dataloader for loading the raw sequences or the one-hot encoded matrix from the fasta file, and the target from the partition CSV file as the batch. 

TO BE UPDATED SOON..

Alternatively, you can also implement your own dataloader.

## Task 1b: for promoter-level activities
# Data

# Demo data load

