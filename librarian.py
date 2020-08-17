from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from fdtd_1d.constants import c0, BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY

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


class QPM_Length_benchmark(Case):
    def __init__(self, dir_name, zero_padding=0):
        super().__init__()
        self.path += '/qpm_harmonic_length/' + dir_name
        self.dt = None
        self.dx = None
        self.timesteps = None
        self.merged_data = None
        self.abs_fft = []
        self.abs_fft_sqrd = []
        self.normalized_abs_fft = []
        self.normalized_abs_fft_sqrd = []

        self.first_har = 1.77306e15
        self.second_har = 2 * self.first_har
        self.zero_padding = zero_padding
        self.positions = []

        self._set_grid_information()
        self._set_observer_positions()
        self._set_merge_data()
        self._recreate_fft()
        self.relative_width = (np.array(self.positions) - 5) * self.dx


    def _set_observer_positions(self):
        for file in os.listdir(self.path):
            string = os.path.basename(file)
            string = string[2:]
            self.positions.append(int(string[:-4]))

    def _set_grid_information(self):
        first_file = self.path + '/E_5.csv'
        df = pd.read_csv(first_file, sep=',', header=None, nrows=1)
        self.dt = float(df[1])
        self.dx = self.dt * c0
        self.timesteps = int(df[3])

    # TODO ERROR IN MERGING DATA
    def _set_merge_data(self):
        collected_df = pd.DataFrame()
        for pathfile, pos in zip(os.listdir(self.path), self.positions):
            indv_df = pd.read_csv(self.path+'/'+str(pathfile), sep=',', header=None, skiprows=[0])
            indv_df = indv_df.T
            collected_df[pos] = indv_df[0]
        self.merged_data = collected_df.to_numpy()

    def show_trace(self):
        im = plt.imshow(self.merged_data**2, cmap='magma', aspect='auto')
        plt.colorbar(im, orientation='horizontal')
        plt.show()

    def _recreate_fft(self):
        x, y = np.shape(self.merged_data)
        padded_data = []
        timestep_duration = self.timesteps + self.zero_padding
        self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timestep_duration // 2) / self.first_har
        for obs_ind in range(x):
            padded_data.append(np.pad(self.merged_data[obs_ind], (0, self.zero_padding), 'constant', constant_values=0))
            data_fft = np.fft.fft(padded_data[obs_ind])
            self.abs_fft.append(2/timestep_duration * np.abs(data_fft))
            self.abs_fft_sqrd.append(np.array(self.abs_fft[obs_ind]) ** 2)
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[obs_ind])
            self.normalized_abs_fft_sqrd.append(np.array(self.abs_fft_sqrd[obs_ind])/abs_fft_sqrd_max)

class Case_qpm_harmonic_length_old(Case):
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
        self.first_har = 1.77306e15
        self.second_har = 2 * self.first_har
        self._set_raw_data()
        self.zero_padding = zero_padding
        self._recreate_fft(zero_pad=zero_padding)


    def _set_raw_data(self):
        df = pd.read_csv(self.path, skiprows=[0], sep=',', header=None)
        df = df.T
        self.raw_data = np.array(df[0], float)

    def _recreate_fft(self, zero_pad=0):
        df = pd.read_csv(self.path, sep=',', header=None, nrows=1)
        dt = float(df[1])
        timestep_duration = len(self.raw_data) + zero_pad
        self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * dt), timestep_duration // 2) / self.first_har
        data = np.pad(self.raw_data, (0, zero_pad), 'constant', constant_values=0)
        data_fft = np.fft.fft(data)
        self.abs_fft = 2 / timestep_duration * np.abs(data_fft[0:timestep_duration // 2])
        self.abs_fft_sqrd = self.abs_fft ** 2
        abs_fft_max = np.max(self.abs_fft)
        self.normalized_abs_fft = self.abs_fft / abs_fft_max
        abs_fft_sqrd_max = np.max(self.abs_fft_sqrd)
        self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max


test = QPM_Length_benchmark(dir_name='testing_new_database2')
test.show_trace()


'''Case0 = Case_qpm_harmonic_length(filename='P_testing_benchmark_obj.csv', zero_padding=80000)
Case0 = Case_qpm_harmonic_length_old(filename='E_5_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv', zero_padding=0)
Case1 = Case_qpm_harmonic_length_old(filename='E_269_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv', zero_padding=0)
Case2 = Case_qpm_harmonic_length_old(filename='E_533_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv', zero_padding=0)
Case3 = Case_qpm_harmonic_length_old(filename='E_797_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv', zero_padding=0)
Case4 = Case_qpm_harmonic_length_old(filename='E_1061_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv', zero_padding=0)
Case5 = Case_qpm_harmonic_length_old(filename='E_1325_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv', zero_padding=0)
Case6 = Case_qpm_harmonic_length_old(filename='E_1589_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv', zero_padding=0)
'''
'''
def show_trace(paths, positions):
    collected_df = pd.DataFrame()
    for path, pos in zip(paths, positions):
        indv_df = pd.read_csv(path, sep=',', header=None, skiprows=[0])
        indv_df = indv_df.T
        collected_df[pos] = indv_df[0]
    data_matrix = collected_df.to_numpy()
    im = plt.imshow(np.abs(data_matrix), cmap='magma', aspect='auto')
    plt.colorbar(im, orientation='horizontal')
    plt.show()
'''





'''
positions = [5 + 88*i for i in range(50)]
print(positions)
path_names = ['fdtd_1d/saved_data/qpm_harmonic_length/E_'+str(pos)+'_two_lambda_but3long_30000_courant_1_timestep_peak_8000_analyze_shg_even_smaller_peak.csv' for pos in positions]
show_trace(path_names, positions)

'''

'''fig, axes = plt.subplots(ncols=2, nrows=1)
#axes[0][0].plot(Case1.omega, Case1.abs_fft, color=TEAL)
axes[1].set_yscale(value='log')
#axes[1][0].plot(Case2.omega, Case2.abs_fft, color=ORANGE)
axes[0].plot(Case0.raw_data, color=GREY)
axes[1].plot(Case0.omega, Case0.normalized_abs_fft_sqrd, color=GREY)
axes[0].plot(Case1.raw_data, color=ORANGE)
axes[1].plot(Case1.omega, Case1.normalized_abs_fft_sqrd, color=ORANGE)
axes[0].plot(Case2.raw_data, color=TEAL)
#axes[1][1].set_yscale(value='log')
axes[1].plot(Case2.omega, Case2.normalized_abs_fft_sqrd, color=TEAL)
axes[1].plot(Case3.omega, Case3.normalized_abs_fft_sqrd, color=CYAN)
axes[0].plot(Case3.raw_data, color=CYAN)
axes[1].plot(Case4.omega, Case4.normalized_abs_fft_sqrd, color=BLUE)
axes[0].plot(Case4.raw_data, color=BLUE)
axes[1].plot(Case5.omega, Case5.normalized_abs_fft_sqrd, color=RED)
axes[0].plot(Case5.raw_data, color=RED)
axes[1].plot(Case6.omega, Case6.normalized_abs_fft_sqrd, color=MAGENTA)
axes[0].plot(Case6.raw_data, color=MAGENTA)
'''
#im = axes.imshow(np.abs(merge_data(path_names, positions), cmap='magma', aspect='auto'))
#plt.colorbar(im, orientation='horizontal')
#plt.show()





