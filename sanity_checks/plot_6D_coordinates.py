import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt

# Load the 6D coordinates, from Sixtracklib tracking
studies_list = ['oldSixtracklib', 'newSixtracklib']

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
fig.suptitle('6D coordinates')

for index, study in enumerate(studies_list):
    mylinestyle = '-'
    if index == 1:
        mylinestyle='--'
    tbt = pkl.load(open(f'tbt_{study}.pkl', 'rb'))
    print(tbt.keys())

    ax1.plot(tbt['turn'], 1e3*np.array(tbt['x']), linestyle=mylinestyle)
    ax2.plot(tbt['turn'], 1e4*np.array(tbt['px']), linestyle=mylinestyle)
    ax3.plot(tbt['turn'], 1e3*np.array(tbt['y']), linestyle=mylinestyle)
    ax4.plot(tbt['turn'], 1e4*np.array(tbt['py']), linestyle=mylinestyle)
    ax5.plot(tbt['turn'], 1e2*np.array(tbt['sigma']), linestyle=mylinestyle)
    ax6.plot(tbt['turn'], 1e4*np.array(tbt['delta']), linestyle=mylinestyle, label=study)

ax1.set_ylabel('x [mm]')
ax2.set_ylabel('px [1e-4]')
ax3.set_ylabel('y [mm]')
ax4.set_ylabel('py [1e-4]')
ax5.set_ylabel('sigma [cm]')
ax6.set_ylabel('delta [1e-4]')
ax5.set_xlabel('Turns')
ax6.set_xlabel('Turns')

# Put a legend below current axis
plt.legend(loc='upper center', bbox_to_anchor=(-0.3, 4.2),
          fancybox=True, ncol=5)

plt.subplots_adjust(wspace=0.4, hspace=0.4)

savefig = False
if savefig:
    ax1.set_xlim(0, 150)
    ax2.set_xlim(0, 150)
    ax3.set_xlim(0, 150)
    ax4.set_xlim(0, 150)
    ax5.set_xlim(0, 150)
    ax6.set_xlim(0, 150)
    plt.savefig('oldvsnewSixtracklib')
else:
    plt.show()
