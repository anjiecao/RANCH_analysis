import pandas as pd
import numpy as np
import ipdb
from pathlib import Path
import os
import argparse

# set directory to here
#os.chdir(Path(__file__).absolute().parent)

def pooled_sd(group):
    group = group.to_frame()
    n = 20  # Size of each subset
    #print("before ssd")
    ssd = np.sum((n - 1) * group['sample_n_std']**2)
    # print(ssd)
    count = group['sample_n_std'].count()
    return np.sqrt(ssd / (n * count - count))



def parse_data(args):

    sim_data = pd.read_csv(args.sim_data_path)
    # Load trial and param csvs
    trial_info = pd.read_csv(args.trial_info_path)
    param_info = pd.read_csv(args.param_info_path)
    

    grouped_sim_data = sim_data.groupby(['trial_id', 'param_id']).agg({
        'sample_n_mean': 'mean',  
        'sample_n_std': pooled_sd
    }).reset_index()

    grouped_sim_data = grouped_sim_data.merge(trial_info, on="trial_id", how="left") \
                   .merge(param_info, on="param_id", how="left")
    
    grouped_sim_data.to_csv(args.output_path, index = False)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("sim_data_path", type=str, help="Path to csv with the simulation results")
    parser.add_argument("trial_info_path", type=str, help="Path to trial_info csv")
    parser.add_argument("param_info_path", type=str, help="Path to param_info csv")
    parser.add_argument("output_path", type=str, help="Path to the output file")
    args = parser.parse_args()
    parse_data(args)

