import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt

import simulation_parameters as pp

# plotting parameters
params = {'legend.fontsize': 25,
          'figure.figsize': (12.5, 10.5),
          'axes.labelsize': 25,
          'axes.titlesize': 25,
          'xtick.labelsize': 25,
          'ytick.labelsize': 25,
          'image.cmap': 'jet',
          'lines.linewidth': 1,
          'lines.markersize': 8,
          'font.family': 'sans-serif'}

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)

# load the lattice parameters, needed for the calculation of the closed orbit, from the twiss file
infile = open(pp.input_dir+'/twiss_sanity_check.pkl', 'rb')
twiss_for_CO = pickle.load(infile) # type : dict
beta_y = twiss_for_CO['betay_init']
beta_y_CC1 = twiss_for_CO['betay_cc1']
muy = twiss_for_CO['muy_cc1']

# parameters of the CC kick
V_cc1 = pp.cravity1_voltage
f_cc1 = 400e6 #Hz
ps_cc1 = pp.cravity1_phase 

# general parameters
Qy = pp.Qy
E_0 = 26e9 # if you want the exact number use madx
clight = 299792458 # speed of light [m/s]

k = 2 * np.pi * f_cc1 / clight # wavenumber of the cavity

# estimate the expected closed orbit distortion from the standard CC
n_particles = pp.n_macroparticles
start, stop = -0.6, 0.6
step = (stop-start)/n_particles
initial_sigmas = np.arange(start, stop, step)

# vertical kick from the standard CC
deg2rad=np.pi/180.
delta_py_CC1 = V_cc1 * np.sin((ps_cc1 -90.)*deg2rad + k * np.array(initial_sigmas))/E_0
# closed orbit at the start of the lattice
y_co_CC = (np.sqrt(beta_y*beta_y_CC1)) * np.array(delta_py_CC1)*np.cos(2*np.pi*muy - np.pi*Qy)/ (2*np.sin(np.pi*Qy)) 


# load the data from the simulation 
df_1 =  pd.read_pickle(pp.output_dir+'/tbt.pkl') 

# plot
plt.scatter(initial_sigmas, y_co_CC, label = 'theoretical prediction', linewidth = 6)
plt.scatter(df_1['sigma'][pp.n_turns_max - 1], df_1['y'][pp.n_turns_max - 1], label ='sixtracklib', linewidth=0.5)
plt.xlabel('sigma [m]')
plt.ylabel(r'$y_{co}(s)$ [m] ')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('sanity_check_rf_multipole.png')
