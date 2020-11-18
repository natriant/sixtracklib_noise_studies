# You need to give as argument the name of the file with the tbt data.
# example: python job_005b_cmpt_tune_1particle.pyimport NAFFlib ./output/tbt.pkl

import sys
import pickle
import numpy as np
import NAFFlib as pnf
import matplotlib.pyplot as plt
import pandas as pd

deltas = np.linspace(-8e-3, 8e-3, 20)

Qx_list = []
Qy_list = []

for i in deltas:
    my_dict = pickle.load( open('./output/tbt{}.pkl'.format(i), 'rb'))
    x_data = []
    y_data = []
    for turn in range(len(my_dict['x'])):
        x_data.append(my_dict['x'][turn][0])
        y_data.append(my_dict['y'][turn][0])
    
    # Compute the tune
    signal_x = x_data
    signal_y = y_data

    Qx_list.append(pnf.get_tune(np.array(signal_x)))
    Qy_list.append(pnf.get_tune(np.array(signal_y)))

# dictionary of lists
my_dict = {'dp':list(deltas) , 'Qx': Qx_list, 'Qy': Qy_list}   
df = pd.DataFrame(my_dict)  
    
# saving the dataframe  
df.to_csv('tracking_dp_QxQy.csv')


