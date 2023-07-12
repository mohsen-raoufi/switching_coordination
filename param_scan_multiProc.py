import pandas as pd
import matplotlib.pyplot as plt
import SwitchingCoordination as sc
import numpy as np
from tqdm import tqdm

import multiprocessing 
num_processors = multiprocessing.cpu_count()
print("Number of processors: ", num_processors)

from joblib import Parallel, delayed

# Define the method to run multi_processing on MonteCarlo Rep.s
multi_processing_method = "joblib" # "single_process" # "multiprocessing" # "joblib"

from importlib import reload
reload(sc)
# %load_ext autoreload
# %autoreload 2

import time

# simulation parameters
N = 3
ref_time = 5.0
noise_std = 0.5
avg_frequency = 0.0
write_file = False # True

# monte-carlo parameters
n_mc_reps = 32

# make an empty out_data to fill in later
out_data = {}
out_data_list = []

param_scan_dict = {"switchingRate": {"range": np.logspace(0,1,2), "log": True},
                   "N": {"range": np.linspace(3,10,15), "log": False}}

# initialize a parameter dictionary
params = sc.InitParams(N=-1,switchingRate=-1,
                            refTime=ref_time,noiseStd=noise_std,
                            avgFrequency=avg_frequency, writeFile=write_file,showAnimation=False)

start = time.perf_counter()

# first parameters loop
for i_switching_rate, switching_rate in enumerate(tqdm(param_scan_dict['switchingRate']['range'])):
    params['switchingRate'] = switching_rate

    # second parameter loop
    for i_N, N in enumerate(tqdm(param_scan_dict['N']['range'])):
        params['N'] = int(N)

        ## # Normal for loop with single processor
        if multi_processing_method=="single_process":
            for mc_iter in np.arange(n_mc_reps):
            
        #     #perform a single simulation
                out_data_tmp = sc.SingleSimulation(params)
                # make an empty (temporary) dict to put all the data (params + output) into
                tmp_dict = {}
                # put the params into the dict
                for key, val in dict.items(params):
                    tmp_dict[key] = val
                # put the time and order arrays into the dict
                tmp_dict['t'] = np.array(out_data_tmp[0]['t'])
                tmp_dict['order'] = np.array(out_data_tmp[0]['order'])
                # tmp_dict["mc_iter"] = mc_iter

                # append it to the list 
                out_data_list.append(tmp_dict)

        
        # # parallel processing using "multiprocessing" 
        elif multi_processing_method=="multiprocessing":
            with multiprocessing.Pool(processes=num_processors) as pool:
                results = pool.map(sc.SingleSimulation, [params] * n_mc_reps)

                for out_data_tmp in results:
                # out_data_tmp = results
                    tmp_dict = {}
            #         # put the params into the dict
                    for key, val in dict.items(params):
                        tmp_dict[key] = val
                    # put the time and order arrays into the dict
                    tmp_dict['t'] = np.array(out_data_tmp[0]['t'])
                    tmp_dict['order'] = np.array(out_data_tmp[0]['order'])
                    # tmp_dict["mc_iter"] = mc_iter

                    # append it to the list 
                    out_data_list.append(tmp_dict)
            

        ## # parallel processing using "joblib"
        elif multi_processing_method=="joblib":
            results = Parallel(n_jobs=num_processors)(delayed(sc.SingleSimulation)(params) for _ in range(n_mc_reps))
            for out_data_tmp in results:
                tmp_dict = {}
        #         # put the params into the dict
                for key, val in dict.items(params):
                    tmp_dict[key] = val
                # put the time and order arrays into the dict
                tmp_dict['t'] = np.array(out_data_tmp[0]['t'])
                tmp_dict['order'] = np.array(out_data_tmp[0]['order'])
                # tmp_dict["mc_iter"] = mc_iter

                # append it to the list 
                out_data_list.append(tmp_dict)
    

finish = time.perf_counter()

print(f'Finished in {round(finish-start, 2)} second(s)')

# convert it to a pd.df
out_data_df = pd.DataFrame(out_data_list)