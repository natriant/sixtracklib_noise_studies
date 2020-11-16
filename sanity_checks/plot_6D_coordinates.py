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

    ax1.plot(tbt['turn'], tbt['x'], linestyle=mylinestyle)
    ax2.plot(tbt['turn'], tbt['px'], linestyle=mylinestyle)
    ax3.plot(tbt['turn'], tbt['y'], linestyle=mylinestyle)
    ax4.plot(tbt['turn'], tbt['py'], linestyle=mylinestyle)
    ax5.plot(tbt['turn'], tbt['sigma'], linestyle=mylinestyle)
    ax6.plot(tbt['turn'], tbt['delta'], linestyle=mylinestyle)

plt.show()