import pickle
from math import *
import numpy as np
import NAFFlib as pnf
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import simulation_parameters as pp
import scipy.stats as ss

# Plotting parameters
params = {'legend.fontsize': 20,
          'figure.figsize': (9.5, 8.5),
          'axes.labelsize': 22,
          'axes.titlesize': 23,
          'xtick.labelsize': 22,
          'ytick.labelsize': 22,
          'image.cmap': 'jet',
          'lines.linewidth': 3,
          'lines.markersize': 5,
          'font.family': 'sans-serif'}

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)

# Plot the initial distribution 
initial_distribution = pickle.load(open('./input/initial_coordinates.pkl', 'rb'))
fig, ax = plt.subplots(1,1,figsize=(9,8))
ax.plot(np.array(initial_distribution['x'])*1e3, np.array(initial_distribution['y'])*1e3, '.')
ax.set_xlabel('x (mm)', fontsize=20)
ax.set_ylabel('y (mm)', fontsize=20)
plt.grid()
plt.tight_layout()
plt.savefig('./input/initial_condition.png')


tbt_data = pickle.load(open('./output/tbt.pkl', 'rb'))
twiss = pickle.load(open('input/twiss_at_start.pkl', 'rb'))
plane_of_interest = 'y' # 'y' # type:string

n_turns = tbt_data['turn'][-1]
n_particles = len(tbt_data['x'][0])
Qx0 = pp.Qx
Qy0 = pp.Qy

# We need to group the data in a list for each particle.
# (u, pu) <--> (x, px) or (y, py)
u_data = {}
pu_data = {}
for particle in range(n_particles):
    u_data[particle] = []
    pu_data[particle] = []

# maybe even 100 turns are enough
for particle in range(n_particles):
    for turn in range(n_turns):
        u_data[particle].append(tbt_data['{}'.format(plane_of_interest)][turn][particle])
        pu_data[particle].append(tbt_data['p{}'.format(plane_of_interest)][turn][particle])

# Remove any lost particles as NAFF will crash
lost_particles = []
Q_list = []

for particle in range(n_particles):
    if np.isnan(u_data[particle]).any() or np.isnan(pu_data[particle]).any():
        lost_particles.append(particle)
        print('particle {} lost'.format(particle))
    else:        
        signal = u_data[particle]
        Q_list.append(pnf.get_tune(np.array(signal)))


# Use the correct twiss parameters
if plane_of_interest == 'x':
    beta = twiss['betx']
    alpha = twiss['alfx']
else:
    beta = twiss['bety']
    alpha = twiss['alfy']


mean, var  = ss.distributions.expon.fit(Q_list)
print('mean={}, variance={}'.format(mean, var))

x = np.linspace(mean-5*var, mean+5*var, len(Q_list))
fitted_data = ss.distributions.expon.pdf(x, mean, var)

# string with distribution info
textstr = '\n'.join((
    r'$\mu=%.2f$' % (mean, ),
    r'$\sigma=%.5f$' % (var, )))

# Histogram
f = plt.figure(figsize=(9.5,8.5))
ax = f.add_subplot(111)
ax.hist(Q_list, color = 'C0')
ax.plot(x, fitted_data, 'r-')
ax.set_ylabel(r'$\rho (\nu_b)$')
ax.set_xlabel('Betatron tune')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=24,
        verticalalignment='top', bbox=props)

plt.grid()
plt.tight_layout()

savefig = True
if savefig:
    plt.savefig('tune_distribution_ayy1e4.png')
print('ok')
plt.show()
