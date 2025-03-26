import sys
import numpy as np
import pandas as pd
import scipy.stats as st

def process_file(filepath):
    # Example: Process the file (assuming it's a CSV)
    # data = pd.read_csv(filepath)
    
    import datetime
    now = datetime.datetime.now()
    date_time = now.strftime("%Y-%m-%d#%H:%M:%S")
    
    res_out = np.genfromtxt(filepath, dtype=float)
    df = pd.read_csv("/ssd/users/smli/miniHackathon/partitions.csv")
    
    idx = df['labels'] == 'seen genes seen individuals'
    X_obs = df[idx]['log_TPM'].values.reshape(100, -1)
    X_pre = res_out[idx].reshape(100, -1)
    
    R_gene = st.pearsonr(X_obs, X_pre, axis=1)[0]
    R_samp = st.pearsonr(X_obs, X_pre, axis=0)[0]
    MAE = np.mean(np.abs(X_obs - X_pre))
    
    df_all = pd.read_csv('results.csv')
    _rank = np.sum(df_all['Train_PCC_gene'].values >= np.round(np.mean(R_gene), 5)) + 1
    
    df_res = pd.DataFrame({
        '1.Rank': [_rank],
        'UploadTime': [date_time],
        'Train_PCC_gene': [np.round(np.mean(R_gene), 5)],
        'Train_PCC_sample': [np.round(np.mean(R_samp), 5)],
        'Train_MAE': [np.round(MAE, 5)]
    })
    
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
