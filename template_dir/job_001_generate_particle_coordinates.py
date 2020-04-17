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
   

# generate longitudinal distribution
gamma_transition = twiss_summary['gammatr']
circumference = twiss_summary['length']
harmonic_number = pp.harmonic_number
rf_voltage = pp.V_RF
rf_phase = pp.lag_RF_deg / 180. * np.pi
z_max = pp.z_max
JohoParameter = pp.JohoParameter 
n_macroparticles = pp.n_macroparticles
Dsigma_wrt_CO, Ddelta_wrt_CO = generate_longitudinal_distribution(partCO, 
            gamma_transition, circumference, harmonic_number, rf_voltage, 
            rf_phase, z_max, JohoParameter, n_macroparticles)

neps_x, neps_y = pp.neps_x, pp.neps_y
Dx_wrt_CO, Dy_wrt_CO, \
    Dpx_wrt_CO, Dpy_wrt_CO = generate_transverse_distribution(partCO,twiss_at_start,
                            neps_x,neps_y,Ddelta_wrt_CO,n_macroparticles)

with open(pp.input_dir + 'initial_distribution_wrt_CO.pkl', 'wb') as fid:
    pickle.dump({
                'Dx_wrt_CO': Dx_wrt_CO,
                'Dy_wrt_CO': Dy_wrt_CO,
                'Dpx_wrt_CO': Dpx_wrt_CO,
                'Dpy_wrt_CO': Dpy_wrt_CO,
                'Dsigma_wrt_CO': Dsigma_wrt_CO,
                'Ddelta_wrt_CO': Ddelta_wrt_CO,
                }, fid)


# save pysixtrack particles
particles = pysixtrack.Particles(p0c   = pp.p0c, 
                                 x     = partCO.x + Dx_wrt_CO, 
                                 px    = partCO.px + Dpx_wrt_CO, 
                                 y     = partCO.y + Dy_wrt_CO, 
                                 py    = partCO.py + Dpy_wrt_CO, 
                                 sigma = partCO.sigma + Dsigma_wrt_CO,
                                 delta = partCO.delta + Ddelta_wrt_CO)

with open(pp.input_dir + 'pysixtrack.particles', 'wb') as fid:
    pickle.dump(particles, fid)


# save a sixtracklib particle set
ps = sixtracklib.ParticlesSet()
p = ps.Particles(num_particles=n_macroparticles)

for i_part in range(n_macroparticles):
    part = partCO.copy()
    part.x += Dx_wrt_CO[i_part]
    part.px += Dpx_wrt_CO[i_part]
    part.y += Dy_wrt_CO[i_part]
    part.py += Dpy_wrt_CO[i_part]
    part.sigma += Dsigma_wrt_CO[i_part]
    part.delta += Ddelta_wrt_CO[i_part]
    part.partid = i_part
    part.state = 1
    part.elemid = 0
    part.turn = 0

    p.from_pysixtrack(part, i_part)

ps.to_file(pp.input_dir + 'sixtracklib.particles')

