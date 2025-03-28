import pybedtools
import pandas as pd 
import scanpy as sc 
import re 
from pyfaidx import Fasta
from tqdm import tqdm 
import numpy as np 
import anndata 
import json 
import gc 
from typing import Callable, Optional

BASE_DIR = "/ssd/users/cfx"

class DataLoader:
    def __init__(self, 
                 rna_path = "./data/ad_rna_bulk.h5ad",
                 chip_path = "./data/ad_chip_chromsubset.h5ad",
                 trial_path = "./data/tasks/task1/trial1",
                 tf_path = f"{BASE_DIR}/genomes/JASPAR_human_TFs_meme/20250327110737_JASPAR2024_combined_matrices_532533_meme.txt",
                 rna_func: Optional[ Callable[[anndata.AnnData], pd.DataFrame]] = None,
                 tf_func : Optional[ Callable[[pd.DataFrame], pd.DataFrame]] = None,
                 fasta_ref = f"{BASE_DIR}/genomes/hg38/hg38.fa",
                 **tf_func_kwargs):
        """_summary_

        Args:
            rna_path (str, optional):
            chip_path (str, optional):
            trial_path (str, optional):
            tf_path (_type_, optional): Path to TF tokens. If None, self.read_ds() won't contain "X_tf"
            rna_func (Optional[ Callable[[anndata.AnnData], pd.DataFrame]], optional): 
                    takes self.rna, \
                    returns df of (sample, features), same indices as self.rna.obs['term_name']. \
                Defaults to take 20 pcs of log1p(tpm).\
            tf_func (Optional[ Callable[[pd.DataFrame], pd.DataFrame]], optional): \
                    takes df with columns ['tf', 'pwm']\
                    Ret: df with tf as index, should have 'feature' column. \
                Defaults to PWM.
            fasta_ref (_type_, optional): reference genome. Defaults to f"{BASE_DIR}/genomes/hg38/hg38.fa".
        """
        self.rna = sc.read_h5ad(rna_path)
        self.chip = sc.read_h5ad(chip_path)
        self.trial_path = trial_path
        self.tf_path = tf_path
        self.fasta_ref = fasta_ref
        
        self._idx_seqs = {"train":[], "test":[], "val":[], "all":[]}
        self._keys_cttf = {"train":[], "test":[], "val":[], "all":[]}

        with open(f"{trial_path}/config.json") as json_data:
            self.split_config = json.load(json_data)

        self.subset_ds_by_config()
        
        if not ("chip_grouped" in self.__dict__):
            self.pp_chip()
        if not ("rna_grouped" in self.__dict__):
            self.pp_rna(rna_func=rna_func)
        if not ("tf_grouped" in self.__dict__ or tf_path is None):
            self.pp_tf(tf_func=tf_func, **tf_func_kwargs)
        
        self.set_ds_keys()
        pass
    
    def subset_ds_by_config(self):
        ## fixed tf, cross celltypes 
        if len(self.split_config["celltype_TF"]['all_ct'])==0:
            all_ct = list(set(self.chip.obs["Biosample term name"]))
        else:
            all_ct = self.split_config["celltype_TF"]['all_ct']
        if len(self.split_config["celltype_TF"]['all_tf'])==0:
            all_tf = list(set(self.chip.obs["Target of assay"]))
        else:
            all_tf = self.split_config["celltype_TF"]['all_tf']
            
        print(f"subset chip to {all_ct}, {all_tf}")
        self.chip = self.chip[self.chip.obs['Target of assay'].isin(all_tf),:]
        self.chip = self.chip[self.chip.obs['Biosample term name'].isin(all_ct),:]
        self.chip = self.chip[:,self.chip.X.sum(axis=0).flatten()>0]
        self.chip = self.chip[self.chip.X.sum(axis=1).flatten()>0,:]
        print(f"chip shape : {self.chip.shape}")
        
        ## write fasta for subsetted seqs 
        print(f"write chip region fastas to {self.trial_path}/chip_seqs.fasta")
        chip_regions = pd.DataFrame(index = self.chip.var.index, columns=["chrom", "start", "stop"])
        chip_regions['chrom'] = chip_regions.index.str.split(":").str[0]
        chip_regions['start'] = chip_regions.index.str.split(":").str[1].str.split("-").str[0]
        chip_regions['stop'] = chip_regions.index.str.split(":").str[1].str.split("-").str[1]
        seq_pos_2_fasta(seqs=chip_regions, 
                fasta_ref=self.fasta_ref, 
                save_path=f"{self.trial_path}/chip_seqs.fasta")
        return

    
    def set_ds_keys(self):
        for k in ["train", "test", "val"]:
            self._set_seq_keys(key=k)
            self._set_cttf_keys(key=k)
        return
    def _set_seq_keys(self, key):
        ## seq idx from config
        fasta = Fasta(f"{self.trial_path}/chip_seqs.fasta")
        k_all = []          # chrXXX
        for k in fasta.keys():
            k_all.append(k.split(":")[0])
        k_all = list(set(k_all))
         
        if key=="test" or key == "val":
            k_seqs = self.split_config["seq"][key]
        elif key=="train":
            k_rest = self.split_config["seq"]["test"]+self.split_config["seq"]["val"]
            k_seqs = [k for k in k_all if not (k in k_rest)]
        elif key=="all":
            k_seqs = k_all
        
        self._idx_seqs[key] = []
        for i, k in enumerate(fasta.keys()):
            chr_id = k.split(":")[0]
            if chr_id in k_seqs:
                self._idx_seqs[key].append(i)
        return
    def _set_cttf_keys(self, key):
        ## _cttf keys
        k_all = self.chip.obs["_cttf"].to_list()
        k_all = list(set(k_all))
        
        if key=="test" or key == "val":
            k_tfct = self.split_config["celltype_TF"][key]
        elif key=="train":
            k_rest = self.split_config["celltype_TF"]["test"]+self.split_config["celltype_TF"]["val"]
            k_tfct = [k for k in k_all if not (k in k_rest)]
        elif key=="all":
            k_tfct = k_all
        
        self._keys_cttf[key] = k_tfct
        return 
    
    def pp_tf(self, tf_func=None, **tf_func_kwargs):
        '''
        tf_func: 
            input: df with columns ['tf', 'pwm']
            Ret: df with tf index, should have 'feature' column
        '''
        @staticmethod 
        def _tf_func(df_tf:pd.DataFrame):
            df_out = df_tf.drop_duplicates(subset=['tf'])
            df_out.index = df_out['tf']
            df_out['feature'] = df_out['pwm']
            return df_out
        
        if tf_func is None:
            tf_pwms = read_JASPAR_pwms(self.tf_path)
            tf_func = _tf_func 
            tf_func_kwargs["df_tf"] = tf_pwms
        
        self.tf_grouped = tf_func(**tf_func_kwargs)
        return
    
    def pp_rna(self, rna_func = None):
        '''
        rna_func: 
            input: self.rna
            Ret: df of (sample, features), same indices as self.rna.obs['term_name']
        '''
        @staticmethod
        def _rna_func(rna:anndata.AnnData) -> pd.DataFrame:
            sc.pp.log1p(rna)
            sc.pp.pca(rna, n_comps=20)
            rna_pc = pd.DataFrame(rna.obsm['X_pca'], index=self.rna.obs.index)
            return rna_pc 
        if rna_func is None:
            rna_func = _rna_func
        
        k_all = self.rna.obs['term_name'].to_list()
        self.rna.obs.index = k_all 
        df_rna_new = rna_func(self.rna)
        self.rna_grouped = df_rna_new.groupby(self.rna.obs.index).mean()
        return
    
    def pp_chip(self):
        k_all = self.chip.obs['Biosample term name'].astype(str) + "_" + self.chip.obs['Target of assay'].astype(str)
        self.chip.obs["_cttf"] = k_all
        self.chip.obs.index = self.chip.obs["_cttf"]
        
        ## group duplicated rows & columns 
        self.chip_grouped = self.chip.to_df().groupby(self.chip.obs.index).any()
        return 
    
    def read_ds(self, key="test") -> dict:
        seqs = self.read_seqs(key=key)
        X_rna = self.read_rna(key=key)
        Y_chip = self.read_chip(key=key)
        
        ds = {"X_seq": seqs, 
              "X_rna": X_rna,
              "Y_chip": Y_chip}
        
        if not (self.tf_path is None):
            X_tf = self.read_tf(key=key)
            ds['X_tf'] = X_tf
        
        return ds 
    
    
    def read_tf(self, key="test"):
        k_tf = [x.split("_")[1] for x in self._keys_cttf[key]]
        X_tf = [self.tf_grouped.loc[tf,"feature"] for tf in k_tf]   
        return X_tf 
    
    def read_rna(self, key="test"):
        k_ct = [x.split("_")[0] for x in self._keys_cttf[key]]
        X_rna = self.rna_grouped.loc[k_ct,:].to_numpy()
        return X_rna
    
    def read_chip(self, key="test"):
        k_tfct = self._keys_cttf[key]
        
        X = self.chip_grouped.loc[k_tfct,:]
        if len(self._idx_seqs[key]) == 0:
            raise Exception("empty indices, run read_seqs() first")
        X = X.iloc[:,self._idx_seqs[key]].to_numpy()
        return X 
    
    def read_seqs(self, key="test"):
        fasta = Fasta(f"{self.trial_path}/chip_seqs.fasta")
        idx_seqs = self._idx_seqs[key]
        seqs = []
        for i in tqdm(idx_seqs, desc="reading fasta"):
            seqs.append(fasta[i][:].seq)
        return seqs 
    

## writing fasta
def seq_pos_2_fasta(seqs:pd.DataFrame, fasta_ref=f"{BASE_DIR}/genomes/hg38/hg38.fa", save_path=""):
    '''
    seqs: pd dataframe of columns ["chrom", 'start', 'stop'] 
    '''
    a = pybedtools.BedTool.from_dataframe(seqs)
    a = a.sequence(fi=fasta_ref)
    a.save_seqs(save_path)    
    return

## reading jaspar pwms 
class MOTIF:
    def __init__(self, name, tf=None, width=0) -> None:
        self.name=name
        self.tf = tf
        self.width = width
        self.pwm_couter=0
        self.pwm = np.zeros((self.width, 4))
        pass
    
    def read_pwm(self, l):
        w = list(filter(None, re.split(r'[\s]+', l)))
        arr = np.array(w, dtype=float)
        self.pwm[self.pwm_couter, :] = arr
        self.pwm_couter +=1
        return

def read_JASPAR_pwms(meme_file):
    start_line_num, width = 0, 0
    motifs = []
    with open(meme_file, 'r') as handle:
        for i, line in enumerate(handle):
            match_motif = re.search("^MOTIF", line)
            start_pwm = re.search("^letter-probability matrix:", line)
            end_pwm = re.search("^URL", line)
            
            if match_motif:
                keys = re.split(r'[\s]', line)
                motif = keys[1]
                tf = keys[2].split('.')[-1]
            
            if start_pwm:
                width = re.search(r'w= [(0-9)]+', line)
                width = int(line[width.start()+3: width.end()])
                mot = MOTIF(motif, tf, width)
                start_line_num = i
            
            pwm_lines = i>start_line_num and i<=start_line_num+width
            if pwm_lines:
                mot.read_pwm(line)
            
            if end_pwm:
                motifs.append(mot)

    motif_name = [mot.name for mot in motifs]
    tf_name = [mot.tf for mot in motifs]
    pwm = [mot.pwm for mot in motifs]
    df = pd.DataFrame.from_dict({'motif':motif_name ,'tf': tf_name, 'pwm': pwm})
    # df['tf']=df['tf'].str.split('::')
    df = df.explode('tf')
    df = df.sort_values('tf').reset_index(drop=True)
    return df 

def pad_pwm_df(df, pad_val=0):
    '''pad pwms to same length'''
    max_len = list(map(lambda x: x.shape[0], df['pwm']))
    max_len = np.max(max_len)
    print(f"pwms max_len={max_len}")
    pwms = list(map(lambda x: np.pad(x, ((0,max_len-x.shape[0]), (0,0)), constant_values=pad_val),  df['pwm']))
    df['pwm'] = pwms
    pwms = np.asarray(pwms)
    return pwms
