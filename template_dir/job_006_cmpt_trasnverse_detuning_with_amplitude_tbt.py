import pickle
from math import *
import numpy as np
import NAFFlib as pnf
import matplotlib.pyplot as plt
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

# Convert to normalised coordiantes 
u_norm = initial_distribution[plane_of_interest]/sqrt(beta)
print('xmax', max(initial_distribution[plane_of_interest]))
pu_norm = initial_distribution['p{}'.format(plane_of_interest)]*sqrt(beta) + initial_distribution[plane_of_interest]*alpha/sqrt(beta)
# Compute the initial action

J_initial = (u_norm**2 + pu_norm**2)/2
print('J{}_min= {}, J{}_max={}'.format(plane_of_interest, min(J_initial), plane_of_interest, max(J_initial)))

# As the detuning is linear with action (Jx,y) we can compute the total tune spread by substracting the tune of the particle at the minimum action from the one at the maximum action.
# Find minimum and maximum action in the array
J_max = np.amax(J_initial)
J_min = np.amin(J_initial)
# Find the indeces of the maximum and minimum action
index_Jmax = np.where(J_initial == np.amax(J_max))[0][0]
index_Jmin = np.where(J_initial == np.amin(J_min))[0][0]

# Total tune spread
DQ = Q_list[index_Jmax] - Q_list[index_Jmin]
print('total tune spread ={}'.format(DQ))

# Compute the detuning
if plane_of_interest == 'x':
    mytune = Qx0-26. # we are interested only in the fractional part of the tune
else:
    mytune = Qy0-26.
detuning = np.array([i-mytune for i in Q_list])
print('min DQ', min(detuning))

# Linear fit, the slope of the fit should be equal with the corresponding detuning coefficient.
[m_pn, b_pn], cov_pn = np.polyfit(2*J_initial, detuning, 1, cov=True)

fig, ax = plt.subplots(1,1,figsize=(8,7))
ax.plot(np.array(2*J_initial)*1e7, detuning*1e3, 'o', c='b', label = r'$\Delta Q_{}={:.2f}e-3$'.format(plane_of_interest, DQ*1e3))
plt.plot(np.array(2*J_initial)*1e7, (m_pn*2*J_initial+b_pn)*1e3, c = 'r', linewidth=2, label=r'$\alpha_{}={:.2f}$'.format(plane_of_interest, m_pn))
ax.set_xlabel('2J{}'.format(plane_of_interest) + r'$\cdot 10^{-7}$')
ax.set_ylabel('(Q{} -Q{}0)'.format(plane_of_interest, plane_of_interest)+ r'$\cdot 10^{-3}$')
#ax.set_ylim(0.0,1e-4)
plt.tight_layout()
plt.grid()
plt.legend()

savefig = True
if savefig:
    plt.savefig('tune_shift.png')
print('ok')

