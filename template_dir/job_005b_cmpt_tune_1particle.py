# You need to give as argument the name of the file with the tbt data.
# example: python job_005b_cmpt_tune_1particle.pyimport NAFFlib ./output/tbt.pkl

import sys
import pickle
import numpy as np
import NAFFlib as pnf
import matplotlib.pyplot as plt

deltas = np.linspace(-4e-4, 4e-4, 10)
Qx_list = []
Qy_list = []

for i in deltas:
    my_dict = pickle.load( open('./output/tbt_{}.pkl'.format(i), 'rb'))
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

# Plot the expected chroma
qpx = qpy = 2
qx0 = 0.13
qy0 = 0.18
expected_x = qpx*deltas + qx0
expected_y = qpy*deltas + qy0

# H plane
fig, ax = plt.subplots()
ax.plot(deltas, Qx_list, '.-', c='C0', label='simulation')
ax.plot(deltas, expected_x, '--', c='r', label='analytical')
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=10)
plt.ylabel('Qχ', fontsize=15)
plt.xlabel('dp/p', fontsize=15)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('./output/Qx_vs_dpp.png')
plt.show()


# V plane
fig, ax = plt.subplots()
ax.plot(deltas, Qy_list, '.-', c='C1', label='simulation')
ax.plot(deltas, expected_y, '--', c='r', label='analytical')
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=10)
plt.ylabel('Qy', fontsize=15)
plt.xlabel('dp/p', fontsize=15)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('./output/Qy_vs_dpp.png')
plt.show()

