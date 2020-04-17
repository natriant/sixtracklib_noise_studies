import sixtracklib
import pickle

ps = sixtracklib.ParticlesSet().fromfile('./input/sixtracklib.particles')

x = ps.particles[0].x
px = ps.particles[0].px
y = ps.particles[0].y
py = ps.particles[0].py
sigma = ps.particles[0].sigma
delta = ps.particles[0].delta

coordinates = {'x': x, 'px':px, 'y':y, 'py':py, 'sigma':sigma, 'delta':delta}

pickle.dump(coordinates, open('./input/initial_coordinates.pkl', 'wb'))
