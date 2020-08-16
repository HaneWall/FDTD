from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fdtd_1d.constants import BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY

class Case:
    ''' creates a bunch of informations about the file of interest in "saved_data" '''

    def __init__(self):
        self.path = "./fdtd_1d/saved_data"


class Case_qpm_harmonic(Case):
    ''' reconstructs a bunch of information about a file of interest created by a qpm_harmonic benchmark'''

    def __init__(self, zero_padding, filename):
        super().__init__()
        self.path = self.path+'/qpm_harmonic/'+filename
        self.raw_data = None
        self.abs_fft = None
        self.abs_fft_sqrd = None
        self.omega = None
        self.normalized_abs_fft = None
        self.normalized_abs_fft_sqrd = None
        self._set_raw_data()
        self.zero_padding = zero_padding
        self._recreate_fft(zero_pad=zero_padding)

    def _set_raw_data(self):
        df = pd.read_csv(self.path, skiprows=[0], sep=',', header=None)
        df = df.T
        self.raw_data = np.array(df[0], float)[0:22843]

    def _recreate_fft(self, zero_pad=0):
        df = pd.read_csv(self.path, sep=',', header=None, nrows=1)
        dt = float(df[1])
        timestep_duration = len(self.raw_data) + zero_pad
        self.omega = 2*np.pi*np.linspace(0, 1/(2*dt), timestep_duration//2)
        data = np.pad(self.raw_data, (0, zero_pad), 'constant', constant_values=0)
        data_fft = np.fft.fft(data)
        self.abs_fft = 2/timestep_duration * np.abs(data_fft[0:timestep_duration//2])
        self.abs_fft_sqrd = self.abs_fft**2
        abs_fft_max = np.max(self.abs_fft)
        self.normalized_abs_fft = self.abs_fft/abs_fft_max
        abs_fft_sqrd_max = np.max(self.abs_fft_sqrd)
        self.normalized_abs_fft_sqrd = self.abs_fft_sqrd/abs_fft_sqrd_max


class Case_qpm_harmonic_length(Case):
    ''' reconstructs a bunch of information about a file of interest created by a qpm_harmonic benchmark'''

    def __init__(self, zero_padding, filename):
        super().__init__()
        self.path = self.path + '/qpm_harmonic_length/' + filename
        self.raw_data = None
        self.abs_fft = None
        self.abs_fft_sqrd = None
        self.omega = None
        self.normalized_abs_fft = None
        self.normalized_abs_fft_sqrd = None
        self._set_raw_data()
        self.zero_padding = zero_padding
        self._recreate_fft(zero_pad=zero_padding)

    def _set_raw_data(self):
        df = pd.read_csv(self.path, skiprows=[0], sep=',', header=None)
        df = df.T
        self.raw_data = np.array(df[0], float)[0:22900]

    def _recreate_fft(self, zero_pad=0):
        df = pd.read_csv(self.path, sep=',', header=None, nrows=1)
        dt = float(df[1])
        timestep_duration = len(self.raw_data) + zero_pad
        self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * dt), timestep_duration // 2)
        data = np.pad(self.raw_data, (0, zero_pad), 'constant', constant_values=0)
        data_fft = np.fft.fft(data)
        self.abs_fft = 2 / timestep_duration * np.abs(data_fft[0:timestep_duration // 2])
        self.abs_fft_sqrd = self.abs_fft ** 2
        abs_fft_max = np.max(self.abs_fft)
        self.normalized_abs_fft = self.abs_fft / abs_fft_max
        abs_fft_sqrd_max = np.max(self.abs_fft_sqrd)
        self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max


#Case0 = Case_qpm_harmonic_length(filename='P_testing_benchmark_obj.csv', zero_padding=80000)
Case0 = Case_qpm_harmonic_length(filename='E_1775_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv', zero_padding=0)

Case1 = Case_qpm_harmonic_length(filename='E_1480_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv', zero_padding=0)
Case2 = Case_qpm_harmonic_length(filename='E_1185_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv', zero_padding=0)
Case3 = Case_qpm_harmonic_length(filename='E_890_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv', zero_padding=0)
Case4 = Case_qpm_harmonic_length(filename='E_595_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv', zero_padding=0)
Case5 = Case_qpm_harmonic_length(filename='E_5_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv', zero_padding=0)
Case6 = Case_qpm_harmonic_length(filename='E_300_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv', zero_padding=0)


def merge_data(paths, positions):
    collected_df = pd.DataFrame()
    for path, pos in zip(paths, positions):
        indv_df = pd.read_csv(path, sep=',', header=None, skiprows=[0])
        indv_df = indv_df.T
        collected_df[pos] = indv_df[0]
    data_matrix = collected_df.to_numpy()
    im = plt.imshow(np.abs(data_matrix), cmap='magma', aspect='auto')
    plt.colorbar(im, orientation='horizontal')
    plt.show()


fig, axes = plt.subplots(2, 2)
#axes[0][0].plot(Case1.omega, Case1.abs_fft, color=TEAL)
axes[0][1].set_yscale(value='log')
#axes[1][0].plot(Case2.omega, Case2.abs_fft, color=ORANGE)
axes[0][0].plot(Case0.raw_data, color=GREY)
axes[0][1].plot(Case0.omega, Case0.normalized_abs_fft_sqrd, color=GREY)
axes[0][0].plot(Case1.raw_data, color=ORANGE)
axes[0][1].plot(Case1.omega, Case1.normalized_abs_fft_sqrd, color=ORANGE)
axes[0][0].plot(Case2.raw_data, color=TEAL)
axes[1][1].set_yscale(value='log')
axes[0][1].plot(Case2.omega, Case2.normalized_abs_fft_sqrd, color=TEAL)
axes[0][1].plot(Case3.omega, Case3.normalized_abs_fft_sqrd, color=CYAN)
axes[0][0].plot(Case3.raw_data, color=CYAN)
axes[0][1].plot(Case4.omega, Case4.normalized_abs_fft_sqrd, color=BLUE)
axes[0][0].plot(Case4.raw_data, color=BLUE)
axes[0][1].plot(Case6.omega, Case6.normalized_abs_fft_sqrd, color=MAGENTA)
axes[0][0].plot(Case6.raw_data, color=MAGENTA)
axes[0][1].plot(Case5.omega, Case5.normalized_abs_fft_sqrd, color=RED)
axes[0][0].plot(Case5.raw_data, color=RED)
plt.show()


#positions = [5 + 59*i for i in range(50)]
#path_names = ['fdtd_1d/saved_data/qpm_harmonic_length/E_'+str(pos)+'_two_lambda_27000_courant_1_timestep_peak_16000_analyze_shg_smaller_peak_trying_to_absorb.csv' for pos in positions]
#merge_data(path_names, positions)




