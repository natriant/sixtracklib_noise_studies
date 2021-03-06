import numpy as np
import pickle
import os
import datetime
from argparse import ArgumentParser

import pysixtrack
import sixtracklib

from lib.EmittanceCalculation import calculate_emittances
from lib.GenerateNoiseKicks import modulated_rf_phase
import simulation_parameters as pp
from pysixtrack.particles import Particles

###### flags for the type of noise ################
noise_type = 'PN' # 'AN', 'BOTH', PN: Phase noise, AN: Amplitude, BOTH: AN+PN

# Load the sequency from the phase kicks from Coast3-Setting3
path_to_data = '/home/natriant/sixtracklib_cc_test/' #path to the kicks' file
if noise_type == 'AN' or noise_type == 'PN':
    with open(path_to_data+'%path'+'/{}_realNoise_v1.pkl'.format(noise_type), 'rb') as f:
        noiseKicks = pickle.load(f)
        noiseKicks = np.array(noiseKicks)
        print(len(noiseKicks)) # amplitude in rad


# Modulate RF phase
A_sec = %A_sec #75e-12 # amplitude of the RF phase modulation in seconds
mod_period_sec = 4.5e-3 #  period of the RF phase modulation in seconds
f_rf = 400e6 # frequency of the CC RF in Hz (from Mad-x) hardcoded
f_rev = 43.45e3 # revolution frequency of SPS in Hz (doubts about this value)

mod_signal = modulated_rf_phase(A_sec, mod_period_sec, f_rf, pp.n_turns_max+1, f_rev) # amplitude in rad

#####################
parser = ArgumentParser()
parser.add_argument("--device", help="set opencl device")
args = parser.parse_args()

os.makedirs(pp.output_dir, exist_ok=True)

with open(pp.input_dir + 'line.pkl', 'rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid), keepextra=True)

elements = sixtracklib.Elements()
elements.append_line(line)

n_part = pp.n_macroparticles
circumference = line.get_length()

# directory to save the final distribution
parts_distribution_dict = {'x': [], 'px':[], 'y': [], 'py': [], 'sigma': [], 'delta': []}
# directory to save the tbt emittances
tbt_dict = {'turn':[], 'time':[], 'intensity':[], 'neps_x':[], 'neps_y':[], 'std_sigma':[]}

time_cum = 0

if pp.track_with == 'sixtracklib':

    ps = sixtracklib.ParticlesSet().fromfile('input/sixtracklib.particles')

    t_start = datetime.datetime.now()
    print('%s: start tracking %d turns'%(str(t_start)[:-7], pp.n_turns_max))

    if args.device is None:
        job = sixtracklib.TrackJob(elements, ps)
    else:
        job = sixtracklib.TrackJob(elements, ps, device=args.device)


    # Collect elements
    job.collect()
    cravity1_id = line.element_names.index('cravity.1')
    cravity1 = job.beam_elements_buffer.get_object(cravity1_id)
    assert cravity1 is not None

    cravity2_id = line.element_names.index('cravity.2')
    cravity2 = job.beam_elements_buffer.get_object(cravity2_id)
    assert cravity2 is not None


    for turn in range(1, pp.n_turns_max+1):
        
        # Update elements 
        cravity1.set_ksl(pp.cravity1_ks0L_from_turn(turn), 0)
        cravity1.set_ps(pp.cravity1_phase, 0)


        if noise_type == 'AN':
            if measured_noise:
                cravity2.set_ksl(pp.cravity2_ks0L_from_turn(turn)+noiseKicks[turn-1]*pp.cravity2_voltage/pp.p0c, 0)
                cravity2.set_ps(pp.cravity2_phase, 0)
            else:
                cravity2.set_ksl(pp.cravity2_ks0L_from_turn(turn)+noiseKicks[turn-1], 0)
                cravity2.set_ps(pp.cravity2_phase, 0)
        if noise_type == 'PN':
            cravity2.set_ksl(pp.cravity2_ks0L_from_turn(turn), 0)
            rad2deg=180./np.pi # convert rad to degrees
            #if measured_noise: # no scaling factor in the phase
            cravity2.set_ps(pp.cravity2_phase+noiseKicks[turn-1]+mod_signal[turn-1]*rad2deg, 0)
            #else:
            #   cravity2.set_ps(pp.cravity2_phase+noiseKicks[turn-1]*pp.p0c*rad2deg/pp.cravity2_voltage, 0)
        '''        
        if noise_type == 'BOTH':
            rad2deg=180./np.pi # convert rad to degrees
            if measured_noise:
                cravity2.set_ksl(pp.cravity2_ks0L_from_turn(turn)+noiseKicks_an[turn-1]*pp.cravity2_voltage/pp.p0c, 0)
                cravity2.set_ps(pp.cravity2_phase+noiseKicks_pn[turn-1]*rad2deg, 0)
            else:
                cravity2.set_ksl(pp.cravity2_ks0L_from_turn(turn)+noiseKicks_AN[turn-1], 0)
                cravity2.set_ps(pp.cravity2_phase+noiseKicks_PN[turn-1]*pp.p0c*rad2deg/pp.cravity2_voltage, 0)
        '''

        job.push_beam_elements()

        # Tracking
        job.track_until(turn)
        time_cum += circumference / (ps.particles[0].beta0[0]*pysixtrack.Particles.clight)

        job.collect_particles()
        res = ps.particles[0]

        job.push_particles() # maybe that's not needed in the presence of real CC

        indx_alive = np.where(res.state)
        x = res.x[indx_alive]
        px = res.px[indx_alive]
        y = res.y[indx_alive]
        py = res.py[indx_alive]
        sigma = res.sigma[indx_alive]
        delta = res.delta[indx_alive]
        betagamma = res.beta0*res.gamma0
        neps_x, neps_y = calculate_emittances(x, px, y, py, sigma, delta, betagamma[0])
        n_mp_alive = len(indx_alive[0])
        intensity = pp.macro_size*n_mp_alive
        tbt_dict['time'].append(time_cum)
        tbt_dict['turn'].append(turn)
        tbt_dict['neps_x'].append(neps_x)
        tbt_dict['neps_y'].append(neps_y)
        tbt_dict['intensity'].append(intensity)
        tbt_dict['std_sigma'].append(np.std(sigma))
        '''
        For long tracking the tbt coordinates are not damped due to large volume of the file, which results to long running time.
        '''


        if turn in pp.turns_to_print or turn == pp.n_turns_max:       
            print('%s: completed turn %d of %d (%d%%)'%(str(datetime.datetime.now())[:-7], 
                turn, pp.n_turns_max, int(100*(turn+1)/pp.n_turns_max)) )
            with open(pp.output_dir + 'tbt.pkl', 'wb') as fid:
                pickle.dump(tbt_dict, fid)
            fid.close()


        if turn == pp.n_turns_max: # save the final distribution
            parts_distribution_dict['x'].append(x)
            parts_distribution_dict['y'].append(y)
            parts_distribution_dict['px'].append(px)
            parts_distribution_dict['py'].append(py)
            parts_distribution_dict['sigma'].append(sigma)
            parts_distribution_dict['delta'].append(delta)

            with open(pp.output_dir + 'final_distribution.pkl', 'wb') as ff:
                pickle.dump(parts_distribution_dict, ff)
            ff.close()


    t_stop = datetime.datetime.now()
    (t_stop-t_start).total_seconds()
    simulation_info = {'n_turns': turn, 
                       'n_macroparticles': pp.n_macroparticles, 
                       'NoiseType': noise_type,
                       'WhiteNoise': white_noise,
                       'CreateNoiseKicks': create_noise_kicks,
                       'PeakedNoise': peaked_noise,
                       'simulation_time_in_seconds': (t_stop-t_start).total_seconds()}
    with open(pp.output_dir + 'simulation_info.pkl', 'wb') as fid:
            pickle.dump(simulation_info, fid)
