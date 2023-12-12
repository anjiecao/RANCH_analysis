import pandas as pd
import ipdb
from pathlib import Path
import os

def pooled_sd(group):
    n = 20  # Size of each subset
    ssd = np.sum((n - 1) * group['n_samples_sd']**2)
    count = group['n_samples_sd'].count()
    return np.sqrt(ssd / (n * count - count))

# set directory to here
os.chdir(Path(__file__).absolute().parent)

# Load data
EIG_data = pd.read_csv("sim_data/summarized_output_EIG.csv")

grouped_EIG_data = EIG_data.groupby(['trial_id', 'param_id']).agg({
    'n_samples_mean': 'mean',  
    'n_samples_sd': pooled_sd  
}).reset_index()

# Load trial and param csvs
trial_info = pd.read_csv("sim_data/trial_info.csv")
param_info = pd.read_csv("sim_data/eig.csv")

# Merge EIG data with trial and param info based on id columns
grouped_EIG_data = grouped_EIG_data.merge(trial_info, on="trial_id", how="left") \
                   .merge(param_info, on="param_id", how="left")


