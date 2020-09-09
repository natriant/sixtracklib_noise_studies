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

###### flags for choos the type of noise ################
phase_noise = True
amplitude_noise = False

white_noise = True
peaked_noise = False  # for now available only for phase noise


##### create the noise #########

if white_noise:
    print('white noise')
    stdNoise = 1e-8
    noiseKicks = np.random.normal(0, stdNoise, pp.n_turns_max)

if peaked_noise: # only phase noise for now  
    # A. Noise parameters
    phi_0 = 1e-8 # amplitude of noise
    Delta_psi = 0.18 # the peak of the spectrum
    print('peaked noise at {}'.format(Delta_psi))
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
    noiseKicks = phi_0*np.cos(psi_t_list)

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

#tbt_dict = {'turn':[], 'time':[], 'intensity':[], 'neps_x':[], 'neps_y':[],
#            'x': [], 'px':[], 'y': [], 'py': [], 'sigma': [], 'delta': []}
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


    for turn in range(1, pp.n_turns_max+1):
        
        job.track_until(turn)
        time_cum += circumference / (ps.particles[0].beta0[0]*pysixtrack.Particles.clight)

        job.collect_particles()
        res = ps.particles[0]

        if amplitude_noise:
            res.py += noiseKicks[turn-1]*np.sin(2*np.pi*400.789e6/(res.beta0*pysixtrack.Particles.clight)*res.sigma)
        if phase_noise:
            res.py += noiseKicks[turn-1]*np.cos(2*np.pi*400.789e6/(res.beta0*pysixtrack.Particles.clight)*res.sigma)

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
        tbt_dict['std_sigma'].append(np.std(sigma))
        '''
        For long tracking the tbt coordinates are not damped due to large volume of the file, which results to long running time.
        '''
        #tbt_dict['x'].append(x)
        #tbt_dict['y'].append(y)
        #tbt_dict['px'].append(px)
        #tbt_dict['py'].append(py)
        #tbt_dict['sigma'].append(sigma)
        #tbt_dict['delta'].append(delta)

        if turn in pp.turns_to_print or turn == pp.n_turns_max:       
            print('%s: completed turn %d of %d (%d%%)'%(str(datetime.datetime.now())[:-7], 
                turn, pp.n_turns_max, int(100*(turn+1)/pp.n_turns_max)) )
            with open(pp.output_dir + 'tbt.pkl', 'wb') as fid:
                pickle.dump(tbt_dict, fid)
            fid.close()
            
    t_stop = datetime.datetime.now()
    (t_stop-t_start).total_seconds()
    simulation_info = {'n_turns': turn, 
                       'n_macroparticles': pp.n_macroparticles, 
                       'AmplitudeNoise': amplitude_noise,
                       'PhaseNoise': phase_noise,
                       'WhiteNoise': white_noise,
                       'PeakedNoise': peaked_noise,
                       'simulation_time_in_seconds': (t_stop-t_start).total_seconds()}
    with open(pp.output_dir + 'simulation_info.pkl', 'wb') as fid:
            pickle.dump(simulation_info, fid)
