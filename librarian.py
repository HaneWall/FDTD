from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import time
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

class Load_QPM_length(Case):

    def __init__(self, dir_name, zero_padding=0):
        super().__init__()
        self.path += '/qpm_harmonic_length/' + dir_name
        self.zero_padding = zero_padding
        self._load_data()
        self.windowed_data = None
        self.padded_observed_data = None
        self.first_har = 1.77e15
        self.second_har = 2 * self.first_har

    def _load_data(self):
        start_time = time.time()
        self.observed_data = np.load(self.path + '/E_data.npy')
        self.no_observer, self.timesteps = np.shape(self.observed_data)
        self.grid_information = np.load(self.path + '/info.npy')
        self.relative_observer_pos = np.load(self.path + '/relative_ind.npy')
        self.dt, self.dx, self.no_lambdas, self.peak_ts, self.pulse_duration = self.grid_information
        print("loaded in --- %s seconds ---" % (time.time() - start_time))

    def zero_pad(self, number_of_zeros):
        self.padded_observed_data = np.pad(self.observed_data[0:self.no_observer][:], [(0, 0), (0, number_of_zeros)], constant_values=0)

    def set_fft_limits(self, past_from_max, future_from_max):
        self.windowed_data = np.zeros(shape=(self.no_observer, (past_from_max + future_from_max + 1)))

        #TODO: this is pretty bad, because the phase is almost never symmetric
        abs_obs = np.abs(self.padded_observed_data[0:self.no_observer][:])
        observer_max_ts = np.argmax(abs_obs[0:self.no_observer][:], axis=1)
        begin = np.array(observer_max_ts-past_from_max)
        end = np.array(observer_max_ts+future_from_max + 1)
        for obs in range(self.no_observer):
            self.windowed_data[obs][:] = self.padded_observed_data[obs][begin[obs]:end[obs]]


    def fft(self):
        if self.windowed_data is not None:
            x, y = np.shape(self.windowed_data)
            padded_length = y + self.zero_padding
            timesteps = padded_length
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps // 2)
            self.fft_matrix = np.fft.fft(self.windowed_data[0:self.no_observer][:], n=padded_length)
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps//2))
            for obs in range(self.no_observer):
                self.abs_fft[obs][:] = 2/timesteps * np.abs(self.fft_matrix[obs][0:timesteps//2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

        elif self.padded_observed_data is not None:
            timesteps = self.timesteps + self.zero_padding
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps//2)
            self.fft_matrix = np.fft.fft(self.padded_observed_data[0:self.no_observer][:])
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps // 2))
            for obs in range(self.no_observer):
                self.abs_fft = 2 / timesteps * np.abs(self.fft_matrix[obs][0:timesteps//2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

        else:
            timesteps = self.timesteps
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps // 2)
            self.fft_matrix = np.fft.fft(self.observed_data[0:self.no_observer][:])
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps // 2))
            for obs in range(self.no_observer):
                self.abs_fft = 2 / timesteps * np.abs(self.fft_matrix[obs][0:timesteps//2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

    def visualize_over_frequencies(self):
        x = self.relative_observer_pos * self.dx
        y = self.omega /self.first_har
        X, Y = np.meshgrid(x, y)
        fig, axes = plt.subplots()
        c = axes.contourf(X, Y, np.transpose(self.normalized_abs_fft_sqrd), levels=20, locator=ticker.LogLocator(),
                          cmap='magma', norm=LogNorm())
        # fig, axes = plt.subplots()
        # im = axes.imshow(self.normalized_abs_fft_sqrd, cmap='magma', aspect='auto')
        axes.set_xlabel('Distanz in m')
        axes.set_ylabel(r'Harmonische')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.colorbar(c, orientation='horizontal')
        plt.show()

    def visualize_2nd_over_length(self):
        index_of_second = np.argmin(np.abs(self.omega-self.second_har))
        shg_amplitude = np.zeros(self.no_observer)
        for obs in range(self.no_observer):
            shg_amplitude[obs] = self.abs_fft[obs][index_of_second]

        fig, axes = plt.subplots()
        axes.plot(self.relative_observer_pos * self.dx, shg_amplitude)
        axes.set_xlim([0, 2.55e-05])
        axes.grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.show()

    def show_trace(self):
        fig, axes = plt.subplots()
        im = axes.imshow(np.transpose(self.observed_data**2), cmap='magma', aspect='auto')
        '''if self.fft_window:
            for obs_ind in range(len(self.positions)):
                fft_rect = Rectangle(xy=(obs_ind - 0.5, self.fft_window[obs_ind][1]), height=
                -(self.fft_window[obs_ind][1] - self.fft_window[obs_ind][0]), width=1,
                                     facecolor='none', edgecolor='white')
                axes.add_patch(fft_rect)'''
        plt.colorbar(im, orientation='horizontal')
        plt.show()



qpm_test = Load_QPM_length('npy_format_test')
qpm_test.zero_pad(7000)
qpm_test.set_fft_limits(past_from_max=6000, future_from_max=6000)
qpm_test.fft()
qpm_test.visualize_over_frequencies()
qpm_test.visualize_2nd_over_length()








