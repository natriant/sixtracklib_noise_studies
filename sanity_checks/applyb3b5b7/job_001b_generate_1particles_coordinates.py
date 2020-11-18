import numpy as np
import pickle
import sys
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

n_macroparticles = 1

# Set initial transverse and longitudinal offset from the close orbit
# With dispersive contribution
Dx_wrt_CO, Dpx_wrt_CO, Dy_wrt_CO, Dpy_wrt_CO = 9e-4+float(sys.argv[1])*(-0.48342603), 0+float(sys.argv[1])*(-0.01996213), 2*7e-4, 0

# Without dispersive contribution
#Dx_wrt_CO, Dpx_wrt_CO, Dy_wrt_CO, Dpy_wrt_CO = 9e-4, 0, 2*7e-4, 0
#Dsigma_wrt_CO, Ddelta_wrt_CO = 0.05, float(sys.argv[1])


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

ps.to_file(pp.input_dir + 'sixtracklib.particles')

