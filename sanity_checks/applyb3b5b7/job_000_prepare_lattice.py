import pickle
import os
import numpy as np
import matplotlib.pylab as plt

from cpymad.madx import Madx
import pysixtrack
from pysixtrack.particles import Particles
import pysixtrack.be_beamfields.tools as bt

import simulation_parameters as pp


os.makedirs(pp.input_dir, exist_ok=True)

mad = Madx()
mad.options.echo = False
mad.options.info = False
mad.warn = False
mad.chdir('./madx')
mad.call('sps_thin_crabcavity.madx')

for parameter in pp.madx_settings:
    setting = pp.madx_settings[parameter]
    mad.input(f'{parameter} = {setting};')

mad.use(pp.seq_name)

# Include b3b5b7 in MBA and MBB
mad.call('./sps/cmd/sps_setMultipoles_upto7.cmd')
mad.input('exec, set_Multipoles_270GeV;')
mad.call('./sps/cmd/sps_assignMultipoles_upto7.cmd')
mad.input('exec, AssignMultipoles;')
mad.command.readtable(file='err.out', table='errors')
errors = mad.table.errors

# Tune and Chromaticity matching
mad.call('./sps/cmd/sps_matching.cmd')
mad.input('exec, SPS_matchtunes(QH, QV);')
mad.input('exec, SPS_setchroma_Q26(QPH, QPV);')
mad.input('acta.31637, harmon=%d;'%pp.harmonic_number)
mad.input('exec, match_chroma(QPH ,QPV);')


# Power the octupoles
#mad.input('exec, match_octupoles(ayy_val, axy_val);') # use this line, if the input is the ayy, axy coefficients and then matching to klof, klod
mad.input('klof=1.0;')
mad.input('klod=1.0;')


line = pysixtrack.Line.from_madx_sequence(mad.sequence.sps, install_apertures=pp.use_aperture, exact_drift=True)
line.apply_madx_errors(errors)

print(f'MBA after multiple errors: {line.elements[47].knl}')

mad.call('./sps/cmd/sps_ptc_QxQyvsdeltap.cmd')
mad.input('exec, plot_QxQyVSdeltap;')
#quit()

# twiss
twtable = mad.twiss()
with open(pp.input_dir + 'twiss_summary.pkl', 'wb') as fid:
    pickle.dump(twtable.summary, fid)

# twiss for sanity check
names = twtable.name
index_init = np.where(names =='mystart:1')
index_cc1 = np.where(names =='cravity.1:1')

print(twtable.keys())
ddx_init = twtable.dx[index_init]
ddpx_init = twtable.dpx[index_init]
print(ddx_init, ddpx_init)
#quit()

betay_init = twtable.bety[index_init]
betay_cc1 = twtable.bety[index_cc1]
muy_cc1 = twtable.muy[index_cc1]
index_cc2 = np.where(names =='cravity.2:1')
betay_cc2 = twtable.bety[index_cc2]
muy_cc2 = twtable.muy[index_cc2]
my_twiss_dict = {'betay_init': betay_init, 'betay_cc1': betay_cc1, 'muy_cc1':muy_cc1, 'betay_cc2': betay_cc2, 'muy_cc2':muy_cc2}
with open(pp.input_dir + 'twiss_sanity_check.pkl', 'wb') as fid:
        pickle.dump(my_twiss_dict, fid)

# Generate line
#line = pysixtrack.Line.from_madx_sequence(mad.sequence.sps, 
#    install_apertures=pp.use_aperture)

# enable RF
i_cavity = line.element_names.index('acta.31637')
line.elements[i_cavity].voltage = pp.V_RF
line.elements[i_cavity].lag = pp.lag_RF_deg

with open(pp.input_dir + 'line.pkl', 'wb') as fid:
    pickle.dump(line.to_dict(keepextra=True), fid)

part_on_CO = line.find_closed_orbit(guess=[twtable['x'][0], twtable['px'][0],
                                           twtable['y'][0], twtable['py'][0], 
                                           0., 0.], p0c=pp.p0c, method='get_guess')

# Save particle on CO
with open(pp.input_dir + 'particle_on_CO.pkl', 'wb') as fid:
    pickle.dump(part_on_CO.to_dict(), fid)

# Save twiss at start ring
with open(pp.input_dir + 'twiss_at_start.pkl', 'wb') as fid:
    pickle.dump({
        'alfx': twtable.alfx[0],
        'alfy': twtable.alfy[0],
        'betx': twtable.betx[0],
        'bety': twtable.bety[0], 
        'dx':   twtable.dx[0],
        'dy':   twtable.dy[0],
        'dpx':  twtable.dpx[0],
        'dpy':  twtable.dpy[0]}, fid)


