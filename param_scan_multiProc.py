import pandas as pd
import matplotlib.pyplot as plt
import SwitchingCoordination as sc
import numpy as np
from tqdm import tqdm
import pickle

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

<<<<<<< HEAD
param_scan_dict = {"switchingRate": {"range": np.logspace(0,1,2), "log": True},
                   "N": {"range": np.linspace(3,10,15), "log": False}}
=======
param_scan_dict = {"switchingRate": {"range": np.logspace(-3.0, 2.0, 2), "log": True},
                   "N": {"range": np.linspace(3,15,2), "log": False}}
>>>>>>> 40245842f60940cc3a9ecb3e5bac115232d5e440

# initialize a parameter dictionary
params = sc.InitParams(N=-1,switchingRate=-1,
                            refTime=ref_time,noiseStd=noise_std,
                            avgFrequency=avg_frequency, writeFile=write_file,showAnimation=False)

# save the params merged with param_scan_dict to a pickle file for later use
params_to_save = params | param_scan_dict
fileName = 'scan_'+str(n_mc_reps)+'nMC_'+'switchRate'

with open(f'{fileName}_params.pkl', 'wb') as handle:
    pickle.dump(params_to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

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
                 #tmp_dict["mc_iter"] = mc_iter

                # append it to the list 
                out_data_list.append(tmp_dict)
    

finish = time.perf_counter()

print(f'Finished in {round(finish-start, 2)} second(s)')

# convert it to a pd.df
out_data_df = pd.DataFrame(out_data_list)

# add the mean order for each single simulation
out_data_df['meanOrder'] = [np.mean(x) for x in out_data_df['order']]

#save DataFrame to pickle file
out_data_df.to_pickle(f'{fileName}_data.pkl')

# calculate the mean of the average order over monte-carlo reps.
avg_order_over_reps =  np.empty(shape=(len(param_scan_dict['N']['range']), len(param_scan_dict['switchingRate']['range'])))
for i_N, N in enumerate(param_scan_dict['N']['range']):
    for i_switching_rate, switching_rate in enumerate(param_scan_dict['switchingRate']['range']):
        out_data_tmp = out_data_df.loc[out_data_df['switchingRate'] == param_scan_dict['switchingRate']['range'][i_switching_rate]]
        out_data_tmp = out_data_tmp.loc[out_data_tmp['N'] == param_scan_dict['N']['range'][i_N]]
        avg_order_tmp = np.mean(out_data_tmp['meanOrder'][:])
        avg_order_over_reps[i_N, i_switching_rate] = avg_order_tmp

'''
N_0 = param_scan_dict['N']['range'][0]
N_end = param_scan_dict['N']['range'][-1]
rate_0 = param_scan_dict['switchingRate']['range'][0]
rate_end = param_scan_dict['switchingRate']['range'][-1]
mean_order = out_data_df['meanOrder'][0]

# plot order parameter versus time
plt.figure(figsize=(10,4))
plt.title(f'N = {N_0}-{N_end}, switching rate = {rate_0}-{rate_end} (log)')
plt.plot(out_data_df.loc[0]['t'],(out_data_df.loc[0]['order']),'.-')
plt.axhline(out_data_df['meanOrder'][0], linestyle='--', color='r', label='mean order')
plt.xlabel('time')
plt.ylabel('order')
plt.legend()
plt.text(0.2,0.9,f'Finished in {round(finish-start, 2)} second(s), mean order: {mean_order}')
'''