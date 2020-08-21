from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from fdtd_1d.constants import c0, BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm
from matplotlib import ticker
from werkzeug.utils import cached_property


color_spec = [BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY]

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

class Lorentz_Slab_benchmark(Case):


    def __init__(self, dir_name):
        super().__init__()
        self.path += '/harmonic_lorentz_slab/' + dir_name
        self.grids = os.listdir(self.path)
        self.dir_name = dir_name
        self._load_data()


    def _load_data(self):
        self.grid_informations = np.transpose(np.load(self.path+'/info.npy'))
        self.dx = self.grid_informations[0]
        self.timesteps = self.grid_informations[1]
        self.N_lambdas = self.grid_informations[2]
        self.length_media = self.grid_informations[3]
        self.theo_amplitude_merged = np.load(self.path+'/theory_ampl.npy')
        self.theo_phasenunterschied_merged = np.load(self.path+'/theory_phase.npy')
        self.exp_amplitude_merged = np.load(self.path+'/exp_ampl.npy')
        self.exp_phase_merged = np.load(self.path+'/exp_phase.npy')
        self.width_in_dx = np.load(self.path+'/width.npy')


    def visualize(self):
        max_resolution_grid = int(np.argmax(self.N_lambdas))

        fig, axes = plt.subplots(2, 2)
        axes[0][0].grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)
        axes[0][1].grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)
        axes[1][0].grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)
        axes[1][1].grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)
        axes[0][0].set_xlabel('width in m', fontsize=14)
        axes[0][1].set_xlabel('width in m', fontsize=14)
        axes[1][0].set_xlabel('width in m', fontsize=14)
        axes[1][1].set_xlabel('width in m', fontsize=14)
        axes[0][0].set_ylabel('transmitted amplitude ' + r'$E_{z,tr}$', fontsize=14)
        axes[1][0].set_ylabel('Phasenunterschied', fontsize=14)
        axes[1][1].set_ylabel(r'$d(\phi_{exp},\phi_{theo})$', fontsize=14)
        axes[0][1].set_ylabel(r'$E_{tr,theo}$ / $E_{tr,FDTD}$', fontsize=14)

        axes[0][0].plot(self.width_in_dx[max_resolution_grid][0:int(self.length_media[max_resolution_grid])-1] * self.dx[max_resolution_grid],
                        self.theo_amplitude_merged[max_resolution_grid][0:int(self.length_media[max_resolution_grid])-1], label='Theorie', color=ORANGE)
        axes[1][0].plot(self.width_in_dx[max_resolution_grid][0:int(self.length_media[max_resolution_grid])-1] * self.dx[max_resolution_grid],
                        self.theo_phasenunterschied_merged[max_resolution_grid][0:int(self.length_media[max_resolution_grid])-1], color=ORANGE,
                        label='Theorie')

        for grid in range(len(self.dx)):
            axes[0][0].plot(self.width_in_dx[grid][0:int(self.length_media[grid])-1] * self.dx[grid], self.exp_amplitude_merged[grid][0:int(self.length_media[grid])-1], color=color_spec[grid], linestyle='dashed', label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambdas[grid]))
            axes[0][1].plot(self.width_in_dx[grid][0:int(self.length_media[grid])-1] * self.dx[grid], self.theo_amplitude_merged[grid][0:int(self.length_media[grid])-1]/self.exp_amplitude_merged[grid][0:int(self.length_media[grid])-1], color=color_spec[grid], linestyle='dashed', label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambdas[grid]))
            axes[1][0].plot(self.width_in_dx[grid][0:int(self.length_media[grid])-1] * self.dx[grid], self.exp_phase_merged[grid][0:int(self.length_media[grid])-1], color= color_spec[grid], linestyle='dashed', label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambdas[grid]))
            axes[1][1].plot(self.width_in_dx[grid][0:int(self.length_media[grid])-1] * self.dx[grid], np.abs(self.exp_phase_merged[grid][0:int(self.length_media[grid])-1] - self.theo_phasenunterschied_merged[grid][0:int(self.length_media[grid])-1]), color=color_spec[grid], label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambdas[grid]))

        axes[0][0].legend(loc='best')
        axes[0][1].legend(loc='best')
        axes[1][0].legend(loc='best')
        axes[1][1].legend(loc='best')

        plt.show()

class QPM_Length_benchmark(Case):


    def __init__(self, dir_name, zero_padding=0):
        super().__init__()
        self.path += '/qpm_harmonic_length/' + dir_name
        self.dt = None
        self.dx = None
        self.timesteps = None
        self.merged_data = None
        self.windowed_merged_data = None

        self.padded_data = []
        self.abs_fft = []
        self.abs_fft_sqrd = []
        self.normalized_abs_fft = []
        self.normalized_abs_fft_sqrd = []

        self.omega = None
        self.first_har = 1.77306e15
        self.second_har = 2 * self.first_har
        self.zero_padding = zero_padding
        self.positions = []
        self.fft_window = []

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
        self.positions = np.sort(self.positions)

    def _set_grid_information(self):
        first_file = self.path + '/E_5.csv'
        df = pd.read_csv(first_file, sep=',', header=None, nrows=1)
        self.dt = float(df[1])
        self.dx = self.dt * c0
        self.timesteps = int(df[3])

    def _set_merge_data(self):
        collected_df = pd.DataFrame()
        for pos, ind in zip(self.positions, range(len(self.positions))):
            indv_df = pd.read_csv(self.path+'/E_'+str(pos)+'.csv', sep=',', header=None, skiprows=[0])
            indv_df = indv_df.T
            collected_df[ind] = indv_df[0]
        self.merged_data = np.transpose(collected_df.to_numpy())

    def set_fft_limits(self, past_from_max, future_from_max):
        self.windowed_merged_data = np.zeros(shape=(len(self.positions), (past_from_max + future_from_max + 1)))
        for ind in range(len(self.positions)):
            abs_obs = np.abs(self.merged_data[ind])
            observer_max_ts = np.argmax(abs_obs)
            self.fft_window.append([observer_max_ts-past_from_max, observer_max_ts+future_from_max])
            self.windowed_merged_data[ind] = self.padded_data[ind][self.fft_window[ind][0]:self.fft_window[ind][1]+1]
        self._recreate_fft()

    def show_trace(self):
        fig, axes = plt.subplots()
        im = axes.imshow(np.transpose(self.merged_data)**2, cmap='magma', aspect='auto')
        if self.fft_window:
            for obs_ind in range(len(self.positions)):
                fft_rect = Rectangle(xy=(obs_ind-0.5, self.fft_window[obs_ind][1]), height=
                                     -(self.fft_window[obs_ind][1]-self.fft_window[obs_ind][0]), width=1,
                                     facecolor='none', edgecolor='white')
                axes.add_patch(fft_rect)
        plt.colorbar(im, orientation='horizontal')
        plt.show()

    def _recreate_fft(self):
        self.padded_data = []
        if not self.fft_window:                 # fft_window empty
            x, y = np.shape(self.merged_data)
            timestep_duration = self.timesteps + self.zero_padding
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timestep_duration // 2) / self.first_har
            for obs_ind in range(x):
                self.padded_data.append(np.pad(self.merged_data[obs_ind], (0, self.zero_padding), 'constant', constant_values=0))
                data_fft = np.fft.fft(self.padded_data[obs_ind])
                self.abs_fft.append(2/timestep_duration * np.abs(data_fft[0:timestep_duration // 2]))
                self.abs_fft_sqrd.append(np.array(self.abs_fft[obs_ind]) ** 2)
                abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[obs_ind])
                self.normalized_abs_fft_sqrd.append(np.array(self.abs_fft_sqrd[obs_ind])/abs_fft_sqrd_max)

        else:
            self.abs_fft = []
            self.abs_fft_sqrd = []
            self.normalized_abs_fft = []
            self.normalized_abs_fft_sqrd = []

            x, y = np.shape(self.merged_data)
            timestep_duration = self.fft_window[0][1] - self.fft_window[0][0] + self.zero_padding
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timestep_duration // 2) / self.first_har
            for obs_ind in range(x):
                self.padded_data.append(
                    np.pad(self.windowed_merged_data[obs_ind], (0, self.zero_padding), 'constant', constant_values=0))
                data_fft = np.fft.fft(self.padded_data[obs_ind])
                self.abs_fft.append(2 / timestep_duration * np.abs(data_fft[0:timestep_duration // 2]))
                self.abs_fft_sqrd.append(np.array(self.abs_fft[obs_ind]) ** 2)
                abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[obs_ind])
                self.normalized_abs_fft_sqrd.append(np.array(self.abs_fft_sqrd[obs_ind]) / abs_fft_sqrd_max)

    def visualize_over_frequencies(self):
        x = self.relative_width
        y = self.omega
        X, Y = np.meshgrid(x, y)
        fig, axes = plt.subplots()
        c = axes.contourf(X, Y, np.transpose(self.normalized_abs_fft_sqrd), levels=20, locator=ticker.LogLocator(), cmap='magma', norm=LogNorm())
        #fig, axes = plt.subplots()
        #im = axes.imshow(self.normalized_abs_fft_sqrd, cmap='magma', aspect='auto')
        axes.set_xlabel('Distanz in m')
        axes.set_ylabel('Harmonische')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.colorbar(c, orientation='horizontal')
        plt.show()


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


lorentz_test = Lorentz_Slab_benchmark(dir_name='new_data_files')
lorentz_test.visualize()








