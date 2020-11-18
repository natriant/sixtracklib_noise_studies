'''
Compare the chromatic behavior from PTC and Sixtracklib.
The data from both libraries should be in csv format.
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

params = {'legend.fontsize': 37,
          'figure.figsize': (11.5, 9.5),
          'axes.labelsize': 37,
          'axes.titlesize': 37,
          'xtick.labelsize': 37,
          'ytick.labelsize': 37,
          'image.cmap': 'jet',
          'lines.linewidth': 2,
          'lines.markersize': 8,
          'font.family': 'sans-serif'}

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)

ptc_data = pd.read_csv('sps_ptc_QpxQpy5e-1_noklofklod_b3b5b7.csv')
sixtracklib_data = pd.read_csv('sixtracklib_sps_QpxQpy5e-1_noklofklod_b3b5b7.csv')
print(ptc_data.keys())
print(sixtracklib_data.keys())

fig, ax = plt.subplots(1, 1)
# data from ptc
ax.plot(ptc_data['DP'], ptc_data['QX0'], '-o', c='C0')
ax.plot(ptc_data['DP'], ptc_data['QY0'], '-o', c='C1')
# data from sixtracklib
ax.plot(sixtracklib_data['dp'], sixtracklib_data['Qx'], '-^', c='C0')
ax.plot(sixtracklib_data['dp'], sixtracklib_data['Qy'], '-^', c='C1')

# plot for legend
ax.plot(np.arange(-1, -0.5), np.arange(-1, -0.5), '-', linewidth=3, c='C0', label='Qx')
ax.plot(np.arange(-1, -0.5), np.arange(-1, -0.5), '-', linewidth=3, c='C1', label='Qy')
# set axis limits
ax.set_xlim(-9e-3, 9e-3)
ax.set_ylim(0.05, 0.21)

# set labels
ax.set_xlabel('dp/p [1]')
ax.set_ylabel(r'$\mathrm{Q_x,Q_y}$')
ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3))
ax.grid(linestyle='--')
plt.legend(loc=7)

plt.tight_layout()
#plt.savefig('sps_270GeV_QpyQpx5e-1_noklofklod_b3b5b7_PTCvsSixtracklib.png')
plt.show()
