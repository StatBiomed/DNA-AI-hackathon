import sys 
import pandas as pd 
import sklearn.metrics 
import datetime
import numpy as np 
import os 

BASE_DIR = "/ssd/users/cfx/"
Groundtruth_base = f"{BASE_DIR}/DNA-AI-hackathon/task2_Seq2CellxTF/tasks"

def metrics(gt, y):
    auroc = sklearn.metrics.roc_auc_score(y_true=gt, y_score=y)
    auprc = sklearn.metrics.average_precision_score(y_true=gt, y_score=y)
    f1 = sklearn.metrics.f1_score(y_true=gt, y_pred=(y>0.5))
    for v in [auroc, auprc, f1]:
        v = 0 if np.isnan(v) else v 
    return auroc, auprc, f1 

def process_file(filepath, task="a"):
    now = datetime.datetime.now()
    date_time = now.strftime("%Y-%m-%d#%H:%M:%S")
    gt = pd.read_csv(f"{Groundtruth_base}/task{task}/ground_truth.csv", index_col=0)
    res_out = pd.read_csv(filepath, index_col=0)                 ## (n_sample, n_region)
    shared_idx_row = gt.index.intersection(res_out.index)
    shared_idx_col = gt.columns.intersection(res_out.columns)
    
    gt = gt.loc[shared_idx_row, shared_idx_col].to_numpy()
    res_out = res_out.loc[shared_idx_row, shared_idx_col].to_numpy()
    
    aurocs_cross_region, auprcs_cross_region, f1_cross_region = np.zeros((gt.shape[0],)), np.zeros((gt.shape[0],)), np.zeros((gt.shape[0],))
    for i in range(gt.shape[0]):
        _auroc, _auprc, _f1 = metrics(gt[i,:], res_out[i,:])
        aurocs_cross_region[i] = _auroc
        auprcs_cross_region[i] = _auprc
        f1_cross_region[i] = _f1 
    
    m_cross_r = {"aurocs_cross_region": aurocs_cross_region.mean(),
                 "auprcs_cross_region": auprcs_cross_region.mean(),
                 "f1_cross_region": f1_cross_region.mean(),
                 }
    
    if gt.shape[0] > 1:
        aurocs_cross_sample, auprcs_cross_sample, f1_cross_sample = np.zeros((gt.shape[1],)), np.zeros((gt.shape[1],)), np.zeros((gt.shape[1],))
        for i in range(gt.shape[1]):
            _auroc, _auprc, _f1 = metrics(gt[:,i], res_out[:,i])
            aurocs_cross_sample[i] = _auroc
            auprcs_cross_sample[i] = _auprc
            f1_cross_sample[i] = _f1 
        m_cross_s = {"aurocs_cross_sample" : auprcs_cross_sample.mean(), 
                     "auprcs_cross_sample" : auprcs_cross_sample.mean(), 
                     "f1_cross_sample" : f1_cross_sample.mean()
                     }
    else:               ## only one sample 
        m_cross_s = {"aurocs_cross_sample" : "NA", 
                     "auprcs_cross_sample" : "NA", 
                     "f1_cross_sample" : "NA"
                     }
        
    res = {**m_cross_r, **m_cross_s}
    res = {k:np.round(res[k],5) for k in res.keys()}
    
    results_file = f"results/results_{taskId}.csv"
    
    if(os.path.exists(results_file)):
        df_all = pd.read_csv(results_file)
    else:
        df_all = pd.DataFrame({k:[] for k in list(res.keys())+["1.Rank", "UploadTime"]})
    
    rank_select_col = "aurocs_cross_region"
    
    _rank = np.sum(df_all[rank_select_col].values >= res[rank_select_col]) + 1      ## should change rank of all other rows
    res['1.rank'] = _rank
    res['UploadTime'] = date_time
    
    df_res = pd.DataFrame([res])
    # test about the dataframe
    # df_res.to_csv("df_res.csv",index=False)
    
    return df_res.to_csv(index=False)
    



if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python process_file.py <file_path> <taskId>")
        sys.exit(1)
    filepath = sys.argv[1]
    taskId = sys.argv[2]
    result = process_file(filepath, taskId)
    print(result)
