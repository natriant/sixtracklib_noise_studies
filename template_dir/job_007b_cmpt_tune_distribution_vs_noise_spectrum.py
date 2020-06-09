import pickle
from math import *
import numpy as np
import NAFFlib as pnf
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import simulation_parameters as pp
import pandas as pd
from lib.FitDistribution import *

# Plotting parameters
params = {'legend.fontsize': 18,
          'figure.figsize': (9.5, 8.5),
          'axes.labelsize': 18,
          'axes.titlesize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'image.cmap': 'jet',
          'lines.linewidth': 3,
          'lines.markersize': 5,
          'font.family': 'sans-serif'}

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)


tbt_data = pickle.load(open('./output/tbt_ayy1e5.pkl', 'rb'))
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

# Load data, type: pandas
data = pd.Series(Q_list)

# Plot for comparison
plt.figure(figsize=(12,8))
ax = data.plot(kind='hist',  bins=50, color = 'C1', normed=True, alpha=0.5)
# save plot limits
#dataYLim = ax.get_ylim()

# Find best fit distribution
best_distribution, best_fit_name, best_fit_params = best_fit_distribution(data, 200, ax)
best_dist = getattr(st, best_fit_name)
print('The best fit is a {} distribution'.format(best_fit_name))

param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
param_str = ', '.join(['{}={}'.format(k,v) for k,v in zip(param_names, best_fit_params)])

print('param_names',param_names)
print('param_Strings', param_str)


# Compute the mean and the variance of the distribution
if best_fit_name == 'expon':    
    mean, var = best_distribution.stats(loc=best_fit_params[0], scale=best_fit_params[1], moments='mv')
elif best_fit_name == 'gamma':
    mean, var = best_distribution.stats(best_fit_params[0], loc=best_fit_params[1], scale=best_fit_params[2], moments='mv')
elif best_fit_name =='beta':
     mean, var = best_distribution.stats(best_fit_params[0], best_fit_params[1], loc=best_fit_params[2], scale=best_fit_params[3], moments='mv')
else:
    print('computing the mean and variance of this distibution is not yet implemented')

print('mu={}, var={}, simga={}'.format(mean, var, np.sqrt(var)))

# Update plots
#ax.set_ylim(dataYLim)
plt.close()

# Make PDF with best params
pdf = make_pdf(best_dist, best_fit_params)

# Display
plt.figure(figsize=(12,8))
ax = pdf.plot(lw=2, label='PDF', legend=True)
ax = data.plot(kind='hist', bins=50, normed=True, alpha=0.5, label='Data', legend=True, ax=ax)

param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
param_str = ', '.join(['{}={:.3f}'.format(k,v) for k,v in zip(param_names, best_fit_params)])
dist_str = '{}({})'.format(best_fit_name, param_str)

ax.set_title('Best fit distribution \n {} \n mean={:.3f}, sigma={:.5f}'.format(dist_str, mean, np.sqrt(var) ))

#ax.set_title('Best fit distribution:{} \n mu={:.3f}, sigma={:.5f}'.format(best_fit_name, mean, np.sqrt(var)))
ax.set_xlabel('Betatron tune')
ax.set_ylabel(r'$\rho (\nu_b)$')
#plt.savefig('tune_distribution_ayy5e3.png')
plt.show()
plt.close()

######################################################################################
# Define the sampling parameters
T = 1 # sampling rate
N = 1000 # number of turns
f = np.linspace(0, 1/T, N)

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
for i in range(0, N):
    psi_t_list.append(psi_t)
    ksi = np.random.normal(mean, std) # different seed on each turn
    psi_t = psi_t + 2*np.pi*Delta_psi + 2*np.pi*ksi


# D. Construct the noise signal
noiseKicks = phi_0*np.cos(psi_t_list)
# E. Obtain the spectrum of the noise signall
fft = np.fft.fft(noiseKicks)

fig, ax1 = plt.subplots()

#ax1 = pdf.plot(lw=2, label='PDF', legend=True)
ax1 = data.plot(kind='hist', bins=50, color='C1', normed=True, alpha=0.5, label='Data', legend=True)#, ax=ax)

ax2 = ax1.twinx() # instantiate a second axes that shares the same x-axis
ax2.plot(f[:N // 2], np.abs(fft)[:N // 2] * 1 / N, color='k')#, label=i) # 1 / N is a normalization factor
ax2.set_ylim(0.0, 1e-8)
ax2.set_ylabel('Noise amplitude (arbitrary units)')


param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
param_str = ', '.join(['{}={:.3f}'.format(k,v) for k,v in zip(param_names, best_fit_params)])
dist_str = '{}({})'.format(best_fit_name, param_str)

ax1.set_xlim(0.0, 0.5)
ax1.set_title('Best fit distribution \n {} \n mean={:.3f}, sigma={:.5f}'.format(dist_str, mean, np.sqrt(var) ))

#ax.set_title('Best fit distribution:{} \n mu={:.3f}, sigma={:.5f}'.format(best_fit_name, mean, np.sqrt(var)))
ax1.set_xlabel('Betatron tune')
ax1.set_ylabel(r'$\rho (\nu_b)$')
#plt.savefig('tune_distribution_ayy5e3.png')
plt.show()







