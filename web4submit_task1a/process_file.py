import sys
import numpy as np
import pandas as pd
import scipy.stats as st
import logging
import os
logging.basicConfig(level=logging.DEBUG)


# online deloyment
DATA_ROOT = "/ssd/users/smli/miniHackathon/"
#local test
#DATA_ROOT = "/Users/shuminli/Documents/DNA-AI-hackathon/"

def process_file(filepath):
    # Example: Process the file (assuming it's a CSV)
    # data = pd.read_csv(filepath)
    
    import datetime
    now = datetime.datetime.now()
    date_time = now.strftime("%Y-%m-%d#%H:%M:%S")
    
    #res_out = np.genfromtxt(filepath, dtype=float)

    # force the users to upload a CSV file with columns of gene, sample, prediction
    res_out = pd.read_csv(filepath)

    df = pd.read_csv(os.path.join(DATA_ROOT,"partitions.csv"))

    test_labels = ['seen genes seen individuals',"seen genes unseen individuals","unseen genes"]

    res = dict()
    for cur_test_label in test_labels:
        
        idx = df['labels'] == cur_test_label
        X_obs_data = df[idx]

        # align two dataframe
        f1 = res_out["gene"].isin(X_obs_data["gene"].values)
        f2 = res_out["sample"].isin(X_obs_data["sample"].values)
        res_out_subset = res_out[f1&f2]
        res_out_ordered = X_obs_data[['gene', 'sample']].merge(res_out_subset, on=['gene', 'sample'], how='left').dropna(subset="log_TPM")
 
        matrix_obs  = X_obs_data.pivot(index='sample', columns='gene', values='log_TPM')
        matrix_res  = res_out_ordered.pivot(index='sample', columns='gene', values='log_TPM')
       
        if(matrix_obs.shape[0]!=matrix_res.shape[0] or matrix_obs.shape[1]!=matrix_res.shape[1]):
            
            R_cross_gene = np.nan
            R_cross_sample = np.nan
            MAE = np.nan
        else:
            R_cross_gene = np.round(st.pearsonr(matrix_obs, matrix_res, axis=1)[0].mean(),5)
            R_cross_sample = np.round(st.pearsonr(matrix_obs, matrix_res, axis=0)[0].mean(),5)
            MAE = np.round(np.mean(np.abs(matrix_obs - matrix_res)),5)

        cur_test_label = cur_test_label.replace("individuals",'indiv')

        res[f"{cur_test_label} PCC (cross gene)"] = R_cross_gene
        res[f"{cur_test_label} PCC (cross indiv)"] = R_cross_sample
        res[f"{cur_test_label} MAE"] = MAE

    if(os.path.exists('results.csv')):
        df_all = pd.read_csv('results.csv')
    else:
        df_all = pd.DataFrame({"1.Rank":[],"UploadTime":[],
                               "seen genes seen indiv PCC (cross gene)":[],\
                                "seen genes seen indiv PCC (cross indiv)":[],"Train MAE":[],\
                                   "seen genes unseen indiv PCC (cross gene)":[],\
                                    "seen genes unseen indiv PCC (cross indiv)":[],\
                                        "unseen genes all indiv PCC (cross gene)":[],\
                                            "unseen genes all indiv PCC (cross indiv)":[]})
    
    rank_select_col = "seen genes unseen indiv PCC (cross indiv)"
    _rank = np.sum(df_all[rank_select_col].values >= res[rank_select_col]) + 1
    res['1.rank'] = _rank
    res['UploadTime'] = date_time
    
    df_res = pd.DataFrame([res])
    # test about the dataframe
    df_res.to_csv("df_res.csv",index=False)
    
    return df_res.to_csv(index=False)
    # Perform some operations (for simplicity, we just return the first 5 rows)
    # return data.head().to_csv(index=False)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python process_file.py <file_path>")
        sys.exit(1)
    filepath = sys.argv[1]
    result = process_file(filepath)
    print(result)
