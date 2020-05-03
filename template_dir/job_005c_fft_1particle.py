import pickle
import matplotlib.pyplot as plt
import numpy as np

tbt_list = ['noNoise', 'WhiteNoise', 'PeakedNoise0_32_ksi0_002']
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


plt.legend()
plt.grid()
plt.ylabel("Amplitude")
plt.xlabel("Frequency")
plt.yscale('log')
plt.tight_layout()
plt.show()

