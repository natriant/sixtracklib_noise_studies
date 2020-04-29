import numpy as np
import pickle
import os
import datetime
from argparse import ArgumentParser

import pysixtrack
import sixtracklib

from lib.EmittanceCalculation import calculate_emittances
import simulation_parameters as pp
from pysixtrack.particles import Particles


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

tbt_dict = {'turn':[], 'time':[], 'intensity':[], 'neps_x':[], 'neps_y':[],
            'x': [], 'px':[], 'y': [], 'py': [], 'sigma': [], 'delta': []}
#tbt_dict = {'turn':[], 'time':[], 'intensity':[], 'neps_x':[], 'neps_y':[], 'std_sigma':[]}

time_cum = 0


if pp.track_with == 'sixtracklib':

    ps = sixtracklib.ParticlesSet().fromfile('input/sixtracklib.particles')

    t_start = datetime.datetime.now()
    print('%s: start tracking %d turns'%(str(t_start)[:-7], pp.n_turns_max))

    if args.device is None:
        job = sixtracklib.TrackJob(elements, ps)
    else:
        job = sixtracklib.TrackJob(elements, ps, device=args.device)

    # to update elements
    job.collect()
    cravity1_id = line.element_names.index('cravity.1')
    cravity1 = job.beam_elements_buffer.get_object(cravity1_id)
    assert cravity1 is not None

    cravity2_id = line.element_names.index('cravity.2')
    cravity2 = job.beam_elements_buffer.get_object(cravity2_id)
    assert cravity2 is not None


    #1. CREATE WHITE PHASE NOISE 
    # stdPhaseNoise = 1e-8
    # phaseKicks = np.random.normal(0, stdPhaseNoise, pp.n_turns_max) # white noise
    #2. CREATE PEAKED NOISE
    # A. Noise parameters
    phi_0 = 1e-8 # amplitude of noise
    Delta_psi = 0.32 # the peak of the spectrum
    # B. Parameters for ksi 
    mean = 0.0
    std = 0.02 # the rms width of the noise spectrum 
    psi_t = 0
    psi_t_list = [] # list to append the phase of the noise signal
    # C. create the phase of the noise signal
    for i in range(0, pp.n_turns_max):
        psi_t_list.append(psi_t)
        ksi = np.random.normal(mean, std) # different seed on each turn
        psi_t = psi_t + 2*np.pi*Delta_psi + 2*np.pi*ksi
    # D. Construct the noise signal
    phi_noise = phi_0*np.cos(psi_t_list)

    # amplitude noise
    #stdAmpNoise = 1e-8
    #ampKicks = np.random.normal(0, stdAmpNoise, pp.n_turns_max)

    for turn in range(1, pp.n_turns_max+1):

        # update elements
        cravity1.set_ksl(pp.cravity1_ks0L_from_turn(turn), 0)        
        cravity1.set_ps(pp.cravity1_phase, 0)

        cravity2.set_ksl(pp.cravity2_ks0L_from_turn(turn), 0)
        cravity2.set_ps(pp.cravity2_phase, 0)

    
        job.track_until(turn)

        time_cum += circumference / (ps.particles[0].beta0[0]*pysixtrack.Particles.clight)

        job.collect_particles()
        res = ps.particles[0]

        # Uncomment for amplitude noise
        #res.py += ampKicks[turn-1]*np.sin(2*np.pi*400.789e6/(res.beta0*pysixtrack.Particles.clight)*res.sigma)
        # Uncommnet for phase noise
        res.py += phi_noise[turn-1]*np.cos(2*np.pi*400.789e6/(res.beta0*pysixtrack.Particles.clight)*res.sigma) # phase noise
        job.push_particles()

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
        #tbt_dict['std_sigma'].append(np.std(sigma))
        '''
        For long tracking the tbt coordinates are not damped due to large volume of the file, which results to long running time.
        '''
        tbt_dict['x'].append(x)
        tbt_dict['y'].append(y)
        tbt_dict['px'].append(px)
        tbt_dict['py'].append(py)
        tbt_dict['sigma'].append(sigma)
        tbt_dict['delta'].append(delta)

        if turn in pp.turns_to_print or turn == pp.n_turns_max:       
            print('%s: completed turn %d of %d (%d%%)'%(str(datetime.datetime.now())[:-7], 
                turn, pp.n_turns_max, int(100*(turn+1)/pp.n_turns_max)) )
            with open(pp.output_dir + 'tbt.pkl', 'wb') as fid:
                pickle.dump(tbt_dict, fid)

    t_stop = datetime.datetime.now()
    (t_stop-t_start).total_seconds()
    simulation_info = {'n_turns': turn, 
                       'n_macroparticles': pp.n_macroparticles, 
                       'simulation_time_in_seconds': (t_stop-t_start).total_seconds()}
    with open(pp.output_dir + 'simulation_info.pkl', 'wb') as fid:
            pickle.dump(simulation_info, fid)
