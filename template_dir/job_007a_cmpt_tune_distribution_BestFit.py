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


tbt_data = pickle.load(open('./output/tbt_ayy5e3.pkl', 'rb'))
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
ax = data.plot(kind='hist', bins=50, normed=True, alpha=0.5)
# Save plot limits
dataYLim = ax.get_ylim()

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
ax.set_ylim(dataYLim)

# Make PDF with best params
pdf = make_pdf(best_dist, best_fit_params)

# Display
plt.figure(figsize=(12,8))
ax = pdf.plot(lw=2, label='PDF', legend=True)
data.plot(kind='hist', bins=50, normed=True, alpha=0.5, label='Data', legend=True, ax=ax)

param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
param_str = ', '.join(['{}={:.3f}'.format(k,v) for k,v in zip(param_names, best_fit_params)])
dist_str = '{}({})'.format(best_fit_name, param_str)

ax.set_title('Best fit distribution \n {} \n mean={:.3f}, sigma={:.5f}'.format(dist_str, mean, np.sqrt(var) ))

#ax.set_title('Best fit distribution:{} \n mu={:.3f}, sigma={:.5f}'.format(best_fit_name, mean, np.sqrt(var)))
ax.set_xlabel('Betatron tune')
ax.set_ylabel(r'$\rho (\nu_b)$')
#plt.savefig('tune_distribution_ayy5e3.png')
plt.show()
