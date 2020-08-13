import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fdtd_1d.constants import BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY

class Case:
    ''' creates a bunch of informations about the file of interest in "saved_data" '''

    def __init__(self, filename, zero_paddling):
        self.path = "./fdtd_1d/saved_data/"+filename
        self.raw_data = None
        self.abs_fft = None
        self.abs_fft_sqrd = None
        self.omega = None
        self.normalized_abs_fft = None
        self.normalized_abs_fft_sqrd = None
        self._set_raw_data()
        self._recreate_fft(zero_pad=zero_paddling)

    def _set_raw_data(self):
        df = pd.read_csv(self.path, skiprows=[0], sep=',', header=None)
        df = df.T
        self.raw_data = np.array(df[0], float)[0:30000]

    def _recreate_fft(self, zero_pad=0):
        df = pd.read_csv(self.path, sep=',', header=None, nrows=1)
        dt = float(df[1])
        timestep_duration = int(df[3]) + zero_pad
        self.omega = 2*np.pi*np.linspace(0, 1/(2*dt), timestep_duration//2)
        data = np.pad(self.raw_data, (0, zero_pad), 'constant', constant_values=0)
        data_fft = np.fft.fft(data)
        self.abs_fft = 2/timestep_duration * np.abs(data_fft[0:timestep_duration//2])
        self.abs_fft_sqrd = self.abs_fft**2
        abs_fft_max = np.max(self.abs_fft)
        self.normalized_abs_fft = self.abs_fft/abs_fft_max
        abs_fft_sqrd_max = np.max(self.abs_fft_sqrd)
        self.normalized_abs_fft_sqrd = self.abs_fft_sqrd/abs_fft_sqrd_max


#def recreate_fft(data):

#def convert_to_float(data):
'''
df1 = pd.read_csv("./fdtd_1d/saved_data/E_two_lambda_laterpeak.csv", skiprows=[0], sep=',', header=None)
df1 = df1.T
P = np.array(df1[1], float)


fig, axes = plt.subplots()
axes.plot(P, color=ORANGE)
plt.show()
'''
Case0 = Case('just_testing_2.csv', zero_paddling=80000)
Case1 = Case('P_two_lambda_laterpeak_16000_737.csv', zero_paddling=130000)
Case2 = Case('E_two_lambda_laterpeak_16000_737.csv', zero_paddling=130000)

fig, axes = plt.subplots(2, 2)
#axes[0][0].plot(Case1.omega, Case1.abs_fft, color=TEAL)
axes[0][0].plot(Case1.raw_data, color=TEAL)
axes[0][1].set_yscale(value='log')
axes[0][1].plot(Case1.omega, Case1.normalized_abs_fft_sqrd, color=TEAL)
#axes[1][0].plot(Case2.omega, Case2.abs_fft, color=ORANGE)
axes[1][0].plot(Case2.raw_data, color=ORANGE)
axes[1][1].set_yscale(value='log')
axes[1][1].plot(Case2.omega, Case2.normalized_abs_fft_sqrd, color=ORANGE)
plt.show()



