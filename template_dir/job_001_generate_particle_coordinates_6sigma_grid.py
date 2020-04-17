import numpy as np
import pickle
import matplotlib.pyplot as plt

import pysixtrack
import sixtracklib

from lib.GenerateMatchedDistribution import generate_longitudinal_distribution
from lib.GenerateMatchedDistribution import generate_transverse_distribution
import simulation_parameters as pp


# Plotting parameters
params = {'legend.fontsize': 20,
          'figure.figsize': (9.5, 8.5),
          'axes.labelsize': 27,
          'axes.titlesize': 23,
          'xtick.labelsize': 27,
          'ytick.labelsize': 27,
          'image.cmap': 'jet',
          'lines.linewidth': 1,
          'lines.markersize': 5,
          'font.family': 'sans-serif'}

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)


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

#Dx_wrt_CO, Dy_wrt_CO, \
#    Dpx_wrt_CO, Dpy_wrt_CO = generate_transverse_distribution(partCO,twiss_at_start,
#                            neps_x,neps_y,Ddelta_wrt_CO,n_macroparticles)


#### Create the initial condition
steps = np.sqrt(n_macroparticles)
# Α. Compute the geometric emittance
m0 = pp.mass # [eV]
E_rest = m0
E0 = np.sqrt(pp.p0c**2 + E_rest**2)
gamma_0 = E0/E_rest
beta_0 = np.sqrt(1-1/(gamma_0**2))
ex_geom = pp.neps_x/(beta_0*gamma_0)  
ey_geom = pp.neps_y/(beta_0*gamma_0)  

# Β. Compute the rms beam size
sigma_x = np.sqrt(ex_geom*twiss_at_start['betx'])
sigma_y = np.sqrt(ey_geom*twiss_at_start['bety'])

# C. Define the upper limit of the distibution
xmax, ymax = 6*sigma_x, 6*sigma_y

# D. Change to action phase variables just to create the linspace
Jxmin, Jymin = 10**(-13), 10**(-13) # for zero you cannot calculate the tune
Jxmax, Jymax = xmax**2/(2*twiss_at_start['betx']), ymax**2/(2*twiss_at_start['bety'])

Jx = np.linspace(10**(-13), Jxmax, steps)
Jy = np.linspace(10**(-13), Jymax, steps)

# E. Return to x-y
x = np.sqrt(Jx*2*twiss_at_start['betx'])
y = np.sqrt(Jy*2*twiss_at_start['bety'])

# F. meshgrid
xx, yy = np.meshgrid(x, y)

Dx_wrt_CO = xx.flatten()
Dy_wrt_CO = yy.flatten()
Dpx_wrt_CO = np.zeros(n_macroparticles)
Dpy_wrt_CO = np.zeros(n_macroparticles)

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
    #part.px += Dpx_wrt_CO[i_part]
    part.y += Dy_wrt_CO[i_part]
    #part.py += Dpy_wrt_CO[i_part]
    #part.sigma += Dsigma_wrt_CO[i_part]
    #part.delta += Ddelta_wrt_CO[i_part]
    part.partid = i_part
    part.state = 1
    part.elemid = 0
    part.turn = 0

    p.from_pysixtrack(part, i_part)

ps.to_file(pp.input_dir + 'sixtracklib.particles')


plt.scatter(ps.particles[0].x, ps.particles[0].y)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.xlim(-0.001, 0.003)
plt.ylim(-0.001, 0.005) 
plt.tight_layout()
plt.savefig(pp.input_dir+'initial_distribution.png')

