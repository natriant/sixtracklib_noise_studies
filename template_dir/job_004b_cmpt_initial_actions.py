import pickle
from math import * 
import numpy as np
import simulation_parameters as pp

initial_distribution = pickle.load(open('./input/initial_coordinates.pkl', 'rb'))
twiss = pickle.load(open('input/twiss_at_start.pkl', 'rb'))
plane_of_interest = 'y' # 'y' # type:string


# Use the correct twiss parameters
if plane_of_interest == 'x':
    beta = twiss['betx']
    alpha = twiss['alfx']
    e_norm = pp.neps_x
else:
    beta = twiss['bety']
    alpha = twiss['alfy']
    e_norm = pp.neps_y

# Convert to normalised coordiantes 
u_norm = initial_distribution[plane_of_interest]/sqrt(beta)
print('xmax', max(initial_distribution[plane_of_interest]))
pu_norm = initial_distribution['p{}'.format(plane_of_interest)]*sqrt(beta) + initial_distribution[plane_of_interest]*alpha/sqrt(beta)
# Compute the initial action

J_initial = (u_norm**2 + pu_norm**2)/2
print('J{}_min= {} m, J{}_max={} m'.format(plane_of_interest, min(J_initial), plane_of_interest, max(J_initial)))

print('rms J{} = {} m'.format(plane_of_interest, np.std(J_initial)))

# Sanity check: The rms J should be equal with the geometric emittance.
p0c = pp.p0c
E_rest = pp.mass
E_0 = np.sqrt(p0c**2+E_rest**2)
gamma_0 =  E_0/E_rest # gamma realtivistic of the reference particle  
beta_0 = np.sqrt(1-1/gamma_0**2) # beta realtivistic of the reference particle
e_geom = e_norm/(gamma_0*beta_0)
print('e{} geometric = {} m'.format(plane_of_interest, e_geom))

