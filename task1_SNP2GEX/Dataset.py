from enformer_pytorch import Enformer,seq_indices_to_one_hot,str_to_one_hot
import os
import torch
from tqdm import tqdm
import pyfaidx
import pandas as pd


class SampleGeneExpressionDataset(torch.utils.data.Dataset):
    def __init__(self,split,csv_file,ref_path=ref_PATH,\
    consensus_root=fasta_ROOT,seq_len=32000,device='cuda',target_name='log_TPM',\
    selected_gene=[],selected_samples=[],return_ref=False):
        cur_split_file = pd.read_csv(csv_file)
        f1 = cur_split_file['split']==split
        f2 = cur_split_file['sample']=='ref'
        self.consensus_root = consensus_root
        self.target_name = target_name
        self.device = device
        self.dictionary = {'A': 0, 'T': 1, 'C': 2,'G':3}
        self.seq_len = 32000
        self.return_ref = return_ref

        if(len(selected_samples)>0):
            f3 = cur_split_file['sample'].isin(selected_samples)
            self.sample_list = selected_samples
        else:
            print("use all samples for the current dataset")
            f3 = ~f2
            self.sample_list = list(cur_split_file['sample'].unique())

        self.sample_df = cur_split_file[f1&f3].reset_index()
        if(len(selected_gene)!=0):
            self.gene_list = sorted(selected_gene)
            self.sample_df = self.sample_df[self.sample_df['gene'].isin(self.gene_list)].reset_index()
        else:
            self.gene_list = sorted(list(self.sample_df['gene'].unique()))
        
        self.sample_fasta = self.loadSampleFasta()


    def loadSampleFasta(self):

        '''

        return sample_chr fasta dictionary

        '''
        non_ref_sample_data =self.sample_df
        sample_fasta = dict()
        not_existed_samples = []

        for index, cur_sample in tqdm(enumerate(self.sample_list), total=len(self.sample_list),desc="Loading sample fasta files"):
            cur_index= cur_sample

            if(cur_index in sample_fasta):
                continue
            
            cur_fn = [cur_index+'_allele_1.fasta',cur_index+'_allele_2.fasta']
            cur_fasta_data_list = []
         
            for temp_fn in cur_fn:
                cur_fasta_path = os.path.join(self.consensus_root,temp_fn)
                if(not os.path.exists(cur_fasta_path)):
                    not_existed_samples.append(cur_index)
                    continue
                curFastaData = pyfaidx.Fasta(cur_fasta_path)
                cur_fasta_data_list.append(curFastaData)
            sample_fasta[cur_index] = cur_fasta_data_list
        
        self.sample_df = self.sample_df[~self.sample_df['sample'].isin(not_existed_samples)]
        print("sample file not found",not_existed_samples)
        return sample_fasta

    def __len__(self):
        return self.sample_df.shape[0]

    def __getitem__(self,idx):

        row = self.sample_df.iloc[idx]
        #print("idx",idx,'row',row)
        cur_sample = row['sample']
        sample_target = row[self.target_name]
        cur_sample_fasta_file = self.sample_fasta[cur_sample]
        cur_gene_name = row["gene"]
        cur_strand = row["strand"]

        # allele 1 information
        cur_sample_fasta_file_1 = cur_sample_fasta_file[0]
        original_cur_sample_seq_1 =  str(cur_sample_fasta_file_1[cur_gene_name])
       
        # allele 2 information
        cur_sample_fasta_file_2 = cur_sample_fasta_file[1]
        original_cur_sample_seq_2 =  str(cur_sample_fasta_file_2[cur_gene_name])


        if(cur_strand=='-'):

            pass
        if(not self.return_ref):
            return original_cur_sample_seq_1, original_cur_sample_seq_2, sample_target,cur_sample,cur_gene_name
        
    
    