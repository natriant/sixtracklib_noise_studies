import numpy as np
import pickle

from pysixtrack.particles import Particles

from lib.GenerateNoiseKicks import *


input_dir = './input/'
output_dir = './output/'
plot_dir = './png/'

Qx = 26.130
Qy = 26.180
Qpx = 0.5
Qpy = 0.5
ayy_val = 0.0 #1e5 #float(%ayy_val) # 1/m
axy_val = 0.0
# MAD-X parameters dictionary
madx_settings = {'QH':Qx, 'QV':Qy, 'QPH':Qpx, 'QPV':Qpy, 'ayy_val':ayy_val, 'axy_val':axy_val}
use_aperture = True #False #
seq_name = 'sps'

mass = Particles.pmass
p0c = 269.99e9 # [eV]
harmonic_number = 4620
V_RF = 5.088*1e6 # [V]
lag_RF_deg = 180.

bunchlength_rms = 0.0155 #float(%sigma_z_val) # meters
JohoParameter = 4 # longitudinal form factor of binomial distribution (inf = Gaussian)
z_max = np.sqrt((JohoParameter+1.)/2.) * 2 * bunchlength_rms
neps_x = 2e-6 # [m]
neps_y = 2e-6 # [m]
n_macroparticles = 20000 #100000 
number_of_particles = 1 #1e11 
macro_size = number_of_particles/n_macroparticles


# tracking parameters
#track_with = 'pysixtrack'
track_with = 'sixtracklib'
n_turns_max = 500000
turns_between_print = 10000
turns_to_print = range(turns_between_print,n_turns_max+1,turns_between_print)


# parameters for cc1
cravity1_voltage = 0. #1e6 #3e6 [V]
cravity1_phase = 90.
cravity1_ks0L_from_turn = lambda turn: np.interp(turn, [0,1,1e12],
        cravity1_voltage / p0c * np.array([0,1,1]))
cravity1_phase_from_turn = lambda turn: cravity1_phase

# parameters for cc2
cravity2_voltage = 1e6 #1e6 # [V]
cravity2_phase = 0. * np.ones(n_turns_max)  # 270. when two CCs are ON operating in opposite phase 

###### flags for the type of noise ################
rad2deg=180./np.pi # In Sixtracklib phase is set in deg
noise_type = 'PN' # 'AN', 'BOTH', PN: Phase noise, AN: Amplitude, BOTH: AN+PN

# White noise kicks or colored noise
#stdNoise = 1e-8 # fix it rad^2/Hz or V/Hz
#noiseKicks = white_noise(stdNoise, n_turns_max) # rad for PN or 1 for AN
#assert len(noiseKicks) == n_turns_max
#if noise_type == 'PN':
#    cravity2_phase += rad2deg*noiseKicks
#if noise_type = 'AN':

    
# Load the phase noise kicks sequence generated from Coast3-Setting3
path_to_data = '.' #path to the kicks' file
with open(path_to_data+'/{}_realNoise_v1.pkl'.format(noise_type), 'rb') as f:
    noiseKicks = rad2deg*np.array(pickle.load(f))[:n_turns_max]  # in deg
    assert len(noiseKicks) == n_turns_max
cravity2_phase += noiseKicks

# # Modulate crab cavity phase
# A_sec = 0e-12 #75e-12 # amplitude of the CC phase modulation in seconds
# mod_period_sec = 4.5e-3 #  period of the CC phase modulation in seconds
# f_rf = 400e6 # frequency of the CC RF in Hz (from Mad-x) hardcoded
# f_rev = 43.45e3 # revolution frequency of SPS in Hz
# mod_signal = rad2deg*modulated_rf_phase(A_sec, mod_period_sec, f_rf, len(noiseKicks), f_rev) # amplitude in rad
# cravity2_phase += mod_signal

cravity2_ks0L_from_turn = lambda turn: np.interp(turn, [0,200,1e12],
        cravity2_voltage / p0c * np.array([0,1,1])) # remember to rump up (200 turns) when you have 1 cavity, np.interp(turn, [0, 200, 1e12]
cravity2_phase_from_turn = lambda turn: cravity2_phase[turn]
