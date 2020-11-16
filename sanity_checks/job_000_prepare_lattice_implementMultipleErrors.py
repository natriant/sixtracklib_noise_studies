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
mad.chdir('madx')
mad.call('sps_thin_crabcavity.madx')
for parameter in pp.madx_settings:
    setting = pp.madx_settings[parameter]
    mad.input(f'{parameter} = {setting};')

mad.use(pp.seq_name)

# Include b3b5b7 in MBA and MBB
# Sanity check, select example MBA, index 47

mba = mad.sequence.sps.expanded_elements[47]

# Actually apply the errors
mad.call('./sps/cmd/sps_setMultipoles_upto7.cmd')
mad.input('exec, set_Multipoles;')
mad.call('./sps/cmd/sps_assignMultipoles_upto7.cmd')
mad.input('exec, AssignMultipoles;')

mad.command.readtable(file='err.out', table='errors')
errors = mad.table.errors
pysixtrack_elements = pysixtrack.Line.from_madx_sequence(mad.sequence.sps, exact_drift=True)
pysixtrack_elements.apply_madx_errors(errors)

print(f'MBA before multiple errors: {mba}')
print(f'knl = {mba.knl}')

print(f'MBA after multiple errors: {pysixtrack_elements.elements[47].knl}')


