---
permalink: /
layout: page
title: About
list_title: AI Hackathon on Genomic Language Model
---


Welcome to the [AI Hackathon on Genomic Language Model](https://statbiomed.github.io/DNA-AI-hackathon/) 
hosted by the StatBiomed lab at the University of Hong Kong.

Materials are stored in this GitHub repo: 
[https://github.com/StatBiomed/DNA-AI-hackathon](https://github.com/StatBiomed/DNA-AI-hackathon)

### Schedule

- Full-day intense hackathon: April 2nd (**8.30am-8.30pm**) briefing at _Seminar Rm 1A_, G/F, 5 Sassoon Rd, 2025

- Tutorial: March 26th, 2025, **2-3.30pm** at _Seminar Rm 1A_, G/F, 5 Sassoon Rd, by Dr Shumin Li (& Fangxin Cai / Ruiyan Hou)

- Seminar: March 25th, 2025, **4-5pm** at _Lecture Theatre 1_, 21 Sassoon Rd, by Dr Shumin Li


### Candidate models

1. Evo2 ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.02.18.638918v1)): 
   7B or 40B with trained on 9.3 trillion DNA base pairs (in 1 million token 
   contexts)

2. NucleotideTransformer ([paper](https://www.nature.com/articles/s41592-024-02523-z))
   Model: ranging from 50 million up to 2.5 billion parameters and integrating 
   information from 3,202 human genomes and 850 genomes from diverse species

3. Caduceus ([arXiv](https://arxiv.org/abs/2403.03234))
   MambaDNA-based model

4. GET ([paper](https://www.nature.com/articles/s41586-024-08391-z))
   general expression transformer


### Tasks

#### Task1: SNP2GEX

See details in [./task1_SNP2GEX/](https://github.com/StatBiomed/DNA-AI-hackathon/tree/main/task1_SNP2GEX/)

In this task, we will use individuals's genomic sequence (with both common and 
rare variants) to predict its cis gene expression, at both gene level expression (task 1a)
or promoter activities (task 1b).

The success of this task will have a huge impact on precision medicine, 
including the study of the functional effects of rare variants and somatic 
variants, which lack power in common eQTL studies.

Multiple labs have been working on this, but it remains an extremely challenging
task. We are benchmarking multiple state-of-the-art genomic language models.

#### Task2: Seq2CellxTF

See details [./task2_Seq2CellxTF/](https://github.com/StatBiomed/DNA-AI-hackathon/tree/main//task2_Seq2CellxTF/)

In this task, we will use the reference sequence and its regulatory regions to
predict the cell type-specific transcription factor bindings.

As you can see, the input data is limited, mainly the relatively short sequence,
optionally with its cis-gene expression (maybe chromatin open accessibility in 
the future), but the output is extremely large, not only for many TFs but also in diverse 
cell conditions.

This would be an ideal task to demonstrate the zero-shot capability of 
"foundation model"


### Sever
- Hostname: hpcf3.cpos.hku.hk (access only within HKU network)
- 4 computing nodes, each with 4 x L40 GPU cards, 120 thread CPU and 1TB of ram Total 10TB of SSD storage
- Duration is from 31 Mar to 6 April (~ 1 week)
- Jobs submission via PBS scheduler
- Project user home folders (/home2) and project folder (/mnt/project) creation,  with 10TB shared quota assigned

### Sponsors

* HKU's [CPOS](https://cpos.hku.hk/) for the GPU support!
* HKUMed & School of Biomedical Sciences
