# Lattice is prepared using job000_prepare_lattice.py.
# However, the intial conditions of the particle are created here.
import sys
import numpy as np
import pickle

import pysixtrack
import sixtracklib

from lib.GenerateMatchedDistribution import generate_longitudinal_distribution
from lib.GenerateMatchedDistribution import generate_transverse_distribution
import simulation_parameters as pp

with open(pp.input_dir + 'line.pkl', 'rb') as fid:
    line = pysixtrack.Line.from_dict(pickle.load(fid), keepextra=True)

with open(pp.input_dir + 'particle_on_CO.pkl', 'rb') as fid:
    partCO = pysixtrack.Particles.from_dict(pickle.load(fid))

# get beta functions from twiss table
with open(pp.input_dir + 'twiss_at_start.pkl', 'rb') as fid:
    twiss_at_start = pickle.load(fid)

with open(pp.input_dir + 'twiss_summary.pkl', 'rb') as fid:
    twiss_summary = pickle.load(fid)
   
n_macroparticles = 1  
n_turns = 1000

# Set initial transverse and longitudinal offset from the close orbit
Dx_wrt_CO, Dpx_wrt_CO, Dy_wrt_CO, Dpy_wrt_CO = 1e-4, 0, 1e-4, 0
Dsigma_wrt_CO, Ddelta_wrt_CO = 200e-3, 0.0 #, float(sys.argv[1])

print('delta = {}'.format(Ddelta_wrt_CO))

# Save a sixtracklib particle set
ps = sixtracklib.ParticlesSet()
p = ps.Particles(num_particles=n_macroparticles)

part = partCO.copy()
part.x += Dx_wrt_CO
part.px += Dpx_wrt_CO
part.y += Dy_wrt_CO
part.py += Dpy_wrt_CO
part.sigma += Dsigma_wrt_CO
part.delta += Ddelta_wrt_CO
part.partid = 1
part.state = 1
part.elemid = 0
part.turn = 0

p.from_pysixtrack(part, 0) # 0 is the particle index 
print('Initial condition:')
print(ps.particles[0]) # print the initial condition

# Prepare for tracking
print('Prepare for tracking')
elements = sixtracklib.Elements()
elements.append_line(line)

n_part = pp.n_macroparticles
circumference = line.get_length()

tbt_dict = {'turn':[], 'time':[], 'x': [], 'px':[], 'y': [], 'py': [], 'sigma': [], 'delta': []}

# flags for noise type 
dpy_kick = False
real_CC_noise = True

# flags for type of noise
white_noise = False
peaked_noise = True

if white_noise:
    print('white noise')
    stdNoise = 1e-8
    noiseKicks = np.random.normal(0, stdNoise, n_turns)

if peaked_noise: # only phase noise for now  
    # A. Noise parameters
    phi_0 = 1e-6 # amplitude of noise
    Delta_psi = 0.27 # the peak of the spectrum
    print('peaked noise at {}'.format(Delta_psi))
    # B. Parameters for ksi 
    mean = 0.0
    std = 0.2 # the rms width of the noise spectrum 
    psi_t = 0
    psi_t_list = [] # list to append the phase of the noise signal
    # C. create the phase of the noise signal
    for i in range(0, n_turns):
        psi_t_list.append(psi_t)
        ksi = np.random.normal(mean, std) # different seed on each turn
        psi_t = psi_t + 2*np.pi*Delta_psi + 2*np.pi*ksi
    # D. Construct the noise signal
    noiseKicks = phi_0*np.cos(psi_t_list)



job = sixtracklib.TrackJob(elements, ps)

# if you want to update elements, check job002..
time_cum = 0
for turn in range(1, n_turns+1):
    job.push_beam_elements()
    job.track_until(turn)
    time_cum += circumference / (ps.particles[0].beta0[0]*pysixtrack.Particles.clight)


    job.collect_particles()
    res = ps.particles[0]

    if dpy_kick:
        print('dpy_kick')
        # Uncomment for amplitude noise
        #res.py += ampKicks[turn-1]*np.sin(2*np.pi*400.789e6/(res.beta0*pysixtrack.Particles.clight)*res.sigma)
        # Uncommnet for phase noise
        res.py += noiseKicks[turn-1]*np.cos(2*np.pi*400.789e6/(res.beta0*pysixtrack.Particles.clight)*res.sigma) # phase noise
    if real_CC_noise: # only phase noise for now in CC2
        print('real CC kick')
        # First you need to update the elements
        cravity1_id = line.element_names.index('cravity.1')
        cravity1 = job.beam_elements_buffer.get_object(cravity1_id)
        assert cravity1 is not None
        cravity2_id = line.element_names.index('cravity.2')
        cravity2 = job.beam_elements_buffer.get_object(cravity2_id)
        assert cravity2 is not None


        # Now make the changes you want
        # CC1
        cravity1.set_ksl(pp.cravity1_ks0L_from_turn(turn), 0)
        cravity1.set_ps(pp.cravity1_phase, 0)
        # CC2
        cravity2.set_ksl(pp.cravity2_ks0L_from_turn(turn), 0)
        rad2deg=180./np.pi # convert rad to degrees
        cravity2.set_ps(pp.cravity2_phase+noiseKicks[turn-1]*pp.p0c*rad2deg/pp.cravity2_voltage, 0)
        #cravity2.set_ps(pp.cravity2_phase, 0)
        job.push_beam_elements()


    job.push_particles()
    indx_alive = np.where(res.state)
    x = res.x[indx_alive]
    px = res.px[indx_alive]
    y = res.y[indx_alive]
    py = res.py[indx_alive]
    sigma = res.sigma[indx_alive]
    delta = res.delta[indx_alive]
    betagamma = res.beta0*res.gamma0
    n_mp_alive = len(indx_alive[0])
    tbt_dict['time'].append(time_cum)
    tbt_dict['turn'].append(turn)
    tbt_dict['x'].append(x)
    tbt_dict['y'].append(y)
    tbt_dict['px'].append(px)
    tbt_dict['py'].append(py)
    tbt_dict['sigma'].append(sigma)
    tbt_dict['delta'].append(delta)


with open(pp.output_dir + 'tbt.pkl', 'wb') as fid:
                pickle.dump(tbt_dict, fid)
