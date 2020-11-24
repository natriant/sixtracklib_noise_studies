import pickle
from math import * 
import numpy as np
import simulation_parameters as pp

initial_distribution = pickle.load(open('./input/initial_coordinates.pkl', 'rb'))
twiss = pickle.load(open('input/twiss_sanity_check.pkl', 'rb'))


def cmpt_normalised_coordinates(u, up, beta, alpha):
    # (u, up)--> (x,xp) or (y, yp), beta, alpha optic functions 
    u_n = u / np.sqrt(beta)
    up_n = alpha*u/np.sqrt(beta) + np.sqrt(beta) * up
    return u_n, up_n

def cmpt_actions(u_n, up_n):
    J = (1/2)*(u_n**2+up_n**2)
    return J


# Optics at CC2
beta_y = twiss['betay_cc2']
beta_x = twiss['betax_cc2']
alpha_y = twiss['alphay_cc2']
alpha_x = twiss['alphax_cc2']

# normalised emittance
e_norm_x, e_norm_y = pp.neps_x, pp.neps_y

# Coordinates
x, px = initial_distribution['x'], initial_distribution['px']
y, py = initial_distribution['y'], initial_distribution['py']

# Normalised coordinates 
x_n, px_n = cmpt_normalised_coordinates(x, px, beta_x, alpha_x)
y_n, py_n = cmpt_normalised_coordinates(y, py, beta_y, alpha_y)


# Compute actions
Jx_init = cmpt_actions(x_n, px_n)
Jy_init = cmpt_actions(y_n, py_n)

print(f'Initial <Jx> = {np.mean(Jx_init)}')
print(f'Initial <Jy> = {np.mean(Jy_init)}')

# Compute the geometric emittances for Sanity check. It should be: <Jx> = ex_geom and <Jy> = ey_geom
p0c = pp.p0c
E_rest = pp.mass
E_0 = np.sqrt(p0c**2+E_rest**2)
gamma_0 =  E_0/E_rest # gamma realtivistic of the reference particle  
beta_0 = np.sqrt(1-1/gamma_0**2) # beta realtivistic of the reference particle
e_geom_y = e_norm_y/(gamma_0*beta_0)
e_geom_x = e_norm_x/(gamma_0*beta_0)
print(f'e_geom_x = {e_geom_x}, e_geom_y = {e_geom_y}')
print('The Jx is different as the dispersion is not taken into account here (?)')

# Test <J^2> vs <J>^2
print(f'<Jy^2> = {np.mean(Jy_init**2)}, <Jy>^2 = {np.mean(Jy_init)**2}')
print(f'<Jx^2> = {np.mean(Jx_init**2)}, <Jx>^2 = {np.mean(Jx_init)**2}')

# rms J, np.std(J) = np.sqrt(var(j) = np.sqrt(<J^2> - <J>^2>)
print(f'np.std(Jy) = {np.std(Jy_init)}')
print(f'np.sqrt(Var((Jy)) = {np.sqrt(np.mean(Jy_init**2)- np.mean(Jy_init)**2)}')



