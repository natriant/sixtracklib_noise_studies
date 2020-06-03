import pickle
from math import * 
import numpy as np

initial_distribution = pickle.load(open('./input/initial_coordinates.pkl', 'rb'))
twiss = pickle.load(open('input/twiss_at_start.pkl', 'rb'))
plane_of_interest = 'y' # 'y' # type:string


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

print('rms J{} = {}'.format(plane_of_interest, np.std(J_initial)))
