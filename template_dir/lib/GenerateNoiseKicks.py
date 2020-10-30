import numpy as np


def seconds2rad(my_t, my_f_rf):
    phi = my_t*2*np.pi*my_f_rf
    return phi


def rad2sec(my_rad, my_f_rf):
    sec = my_rad/(2*np.pi*my_f_rf)
    return sec


def modulated_rf_phase(my_A_sec, my_mod_period_s, my_f_rf, my_turns_max, my_f_rev):
    # returns the signal of the modulated RF, here the CC. The amplitude in rad
    A_rad = seconds2rad(my_A_sec, my_f_rf)  # amplitude of the RF modulation in rad
    mod_freq_Hz = 1 / my_mod_period_s  # frequency of the RF phase modulation in Hz

    turns = np.arange(my_turns_max)
    time = turns/my_f_rev  # seconds

    # Create the sin signal
    my_signal = A_rad*np.sin(2*np.pi*mod_freq_Hz*time)
    return my_signal


