import pickle
import matplotlib.pyplot as plt
import numpy as np

params = {'legend.fontsize': 18,
          'figure.figsize': (8.5, 6.5),
          'axes.labelsize': 18,
          'axes.titlesize': 18,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'image.cmap': 'jet',
          'lines.linewidth': 1,
          'lines.markersize': 5,
          'font.family': 'sans-serif'}


plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)


tbt_list = ['rmsksi8e-2']
path_to_tbt = './output/'


T = 1 # sampling rate
N = 1000 # number of turns
f = np.linspace(0,1/T, N) # frequencies


for i in tbt_list:
    tbt= pickle.load(open(path_to_tbt+'tbt_{}.pkl'.format(i), 'rb'))
    my_signal = []
    for turn in range(0, N):
        my_signal.append(tbt['y'][turn][0]) 
    fft = np.fft.fft(my_signal)
    plt.plot(f[:N // 2], np.abs(fft)[:N // 2] * 1 / N, label=i) # 1 / N is a normalization factor


    filtered_freq_fft = fft.copy()
    for index, ff in enumerate(f):
        if np.abs(ff) > 0.18-0.05 and np.abs(ff) < 0.18+0.05 :
            print('filtered{}, index{}'.format(ff, index))
            filtered_freq_fft[index] = 0

    plt.plot(f[:N // 2], np.abs(filtered_freq_fft)[:N // 2] * 1 / N, label=i) # 1 / N is a normalization factor

plt.legend()
plt.grid()
plt.ylabel("Amplitude")
plt.xlabel("Frequency")
plt.yscale('log')
plt.tight_layout()
plt.show()

