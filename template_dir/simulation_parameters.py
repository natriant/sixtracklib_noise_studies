import numpy as np

from pysixtrack.particles import Particles


input_dir = './input/'
output_dir = './output/'
plot_dir = './png/'

Qx = 26.130
Qy = 26.180
Qpx = 0.0
Qpy = 0.0
axx_val = 500.0
ayy_val = 500.0
# MAD-X parameters dictionary
madx_settings = {'QH':Qx, 'QV':Qy, 'QPH':Qpx, 'QPV':Qpy, 'axx_val':axx_val, 'ayy_val':ayy_val}
use_aperture = True #False #
seq_name = 'sps'

mass = Particles.pmass
p0c = 269.99e9 #25.92e9
harmonic_number = 4620
V_RF = 0. #2.37*1e6 #2.0*1e6
lag_RF_deg = 180.

bunchlength_rms = 1.55e-3 #0.155 #0.22 # meters
JohoParameter = 4 # longitudinal form factor of binomial distribution (inf = Gaussian)
z_max = np.sqrt((JohoParameter+1.)/2.) * 2 * bunchlength_rms
neps_x = 2e-6 #2.5e-6
neps_y = 2e-6 #2.5e-6 / 1e6
n_macroparticles = 10000 #100000 
number_of_particles = 1 #1e11 
macro_size = number_of_particles/n_macroparticles


# tracking parameters
#track_with = 'pysixtrack'
track_with = 'sixtracklib'
n_turns_max = 100000 #10000
turns_between_print = 100
turns_to_print = range(turns_between_print,n_turns_max+1,turns_between_print)

# parameters for cc1
cravity1_voltage = 3e6 #0 #3e6 [V]
cravity1_phase = 90.
cravity1_ks0L_from_turn = lambda turn: np.interp(turn, [0,1,1e12],
        cravity1_voltage / p0c * np.array([0,1,1]))

# parameters for cc2
cravity2_voltage = 3e6 #1e6 # [V]
cravity2_phase = 270. #90. 
cravity2_ks0L_from_turn = lambda turn: np.interp(turn, [0,1,1e12],
        cravity2_voltage / p0c * np.array([0,1,1])) # remember to rump up when you have 1 cavity

