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
Dsigma_wrt_CO, Ddelta_wrt_CO = 0, float(sys.argv[1])

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

job = sixtracklib.TrackJob(elements, ps)

# if you want to update elements, check job002..
time_cum = 0
for turn in range(1, n_turns+1):
    job.push_beam_elements()
    job.track_until(turn)
    time_cum += circumference / (ps.particles[0].beta0[0]*pysixtrack.Particles.clight)


    job.collect_particles()
    res = ps.particles[0]
    
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


with open(pp.output_dir + 'tbt_{}.pkl'.format(sys.argv[1]), 'wb') as fid:
                pickle.dump(tbt_dict, fid)
