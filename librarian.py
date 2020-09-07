
import matplotlib.pyplot as plt
import os
import time
import numpy as np
from fdtd_1d.constants import c0, eps0, mu0, BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY
from fdtd_1d.analytical_solutions import Qpm_Length_Analytical
from matplotlib.colors import LogNorm
from matplotlib import ticker
from scipy.signal import get_window, hilbert


color_spec = [BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY]

class Case:
    ''' creates a bunch of informations about the file of interest in "saved_data" '''

    def __init__(self):
        self.path = "./fdtd_1d/saved_data"


class Load_QPM_end_P(Case):

    def __init__(self, dir_name, zero_padding=0):
        super().__init__()
        self.path += '/qpm_harmonic/' + dir_name
        self.zero_padding = zero_padding
        self._load_data()
        self.windowed_data = None
        self.padded_observed_data = None
        self.first_har = 1.77e15
        self.second_har = 2 * self.first_har

    def _load_data(self):
        start_time = time.time()
        self.observed_data = np.load(self.path + '/E_P_data.npy')
        self.no_observer, self.timesteps = np.shape(self.observed_data)
        self.grid_information = np.load(self.path + '/info.npy')
        self.relative_observer_pos = np.load(self.path + '/relative_ind.npy')
        self.dt, self.dx, self.no_lambdas, self.peak_ts, self.pulse_duration = self.grid_information
        print("loaded in --- %s seconds ---" % (time.time() - start_time))

    def zero_pad(self, number_of_zeros):
        self.padded_observed_data = np.pad(self.observed_data[0:self.no_observer][:], [(0, 0), (0, number_of_zeros)],
                                           constant_values=0)

    def set_fft_limits(self, past_from_max, future_from_max):
        self.windowed_data = np.zeros(shape=(self.no_observer, (past_from_max + future_from_max + 1)))

        if self.padded_observed_data is not None:
            # print(np.shape(self.padded_observed_data))
            x, y = np.shape(self.padded_observed_data)
            self.analytical_signal = np.empty(shape=(self.no_observer, y), dtype='complex_')
            self.envelope_amplitude = np.zeros(shape=(self.no_observer, y))
            for obs in range(self.no_observer):
                self.analytical_signal[obs] = hilbert(self.padded_observed_data[obs][:])
                self.envelope_amplitude[obs] = np.abs(self.analytical_signal[obs][:])

        else:
            x, y = np.shape(self.observed_data)
            self.analytical_signal = np.empty(shape=(self.no_observer, y), dtype='complex_')
            self.envelope_amplitude = np.zeros(shape=(self.no_observer, y))
            for obs in range(self.no_observer):
                self.analytical_signal[obs] = hilbert(self.observed_data[obs][:])
                self.envelope_amplitude[obs] = np.abs(self.analytical_signal[obs][:])

        self.max_timestep = np.zeros(self.no_observer)
        self.max_timestep = np.argmax(self.envelope_amplitude[0:self.no_observer][:], axis=1)
        begin = np.array(self.max_timestep - past_from_max)
        # print(begin)
        end = np.array(self.max_timestep + future_from_max + 1)

        if self.padded_observed_data is not None:
            for obs in range(self.no_observer):
                self.windowed_data[obs][:] = self.padded_observed_data[obs][begin[obs]:end[obs]]

        else:
            for obs in range(self.no_observer):
                self.windowed_data[obs][:] = self.observed_data[obs][begin[obs]:end[obs]]

    def fft(self):
        if self.windowed_data is not None:
            # print(np.shape(self.padded_observed_data))
            x, y = np.shape(self.windowed_data)
            padded_length = y + self.zero_padding
            timesteps = padded_length
            # window = np.hamming(y)
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps // 2)
            self.fft_matrix = np.fft.fft(self.windowed_data[0:self.no_observer][:], n=padded_length)
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps // 2))
            for obs in range(self.no_observer):
                self.abs_fft[obs][:] = 2 / timesteps * np.abs(self.fft_matrix[obs][0:timesteps // 2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

        elif self.padded_observed_data is not None:
            timesteps = self.timesteps + self.zero_padding
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps // 2)
            self.fft_matrix = np.fft.fft(self.padded_observed_data[0:self.no_observer][:])
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps // 2))
            for obs in range(self.no_observer):
                self.abs_fft[obs][:] = 2 / timesteps * np.abs(self.fft_matrix[obs][0:timesteps // 2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

        else:
            timesteps = self.timesteps
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps // 2)
            self.fft_matrix = np.fft.fft(self.observed_data[0:self.no_observer][:])
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps // 2))
            for obs in range(self.no_observer):
                self.abs_fft[obs][:] = 2 / timesteps * np.abs(self.fft_matrix[obs][0:timesteps // 2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

    def show_spectrum(self):
        normed_omega = self.omega / self.first_har

        fig, axes = plt.subplots(nrows=1, ncols=2)

        axes[0].set_xlim([0, 4])
        axes[1].set_xlim([0, 4])
        axes[0].set_yscale('log')
        axes[1].set_yscale('log')
        axes[0].grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)
        axes[1].grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)

        axes[0].plot(normed_omega, self.normalized_abs_fft_sqrd[1][:], color=TEAL, label='P')
        axes[1].plot(normed_omega, self.normalized_abs_fft_sqrd[0][:], color=TEAL, label='E')
        axes[0].legend(loc='best')
        axes[1].legend(loc='best')

        plt.show()

    def show_P(self):
        fig, axes = plt.subplots()
        axes.plot(self.observed_data[1])

        plt.show()

    def show_E(self):
        fig, axes = plt.subplots()
        axes.plot(self.observed_data[0])

        plt.show()


class Load_Soliton(Case):

    def __init__(self, dir_name):
        super().__init__()
        self.path += '/soliton/' + dir_name
        self._load_data()
        self._do_hilbert_transforms()

    def _load_data(self):
        start_time = time.time()
        self.observed_data = np.load(self.path + '/Int_Pos_E.npy')
        self.dx = 25e-09
        self.peak_timestep, self.pulse_duration = np.load(self.path + '/info.npy')
        self.intensities = np.load(self.path + '/intensities.npy')
        self.propagation_distance = np.load(self.path + '/propagations.npy')
        self.mu_m = self.dx*np.arange(2001) - 1000*self.dx
        print("loaded in --- %s seconds ---" % (time.time() - start_time))

    def _do_hilbert_transforms(self):
        self.abs_observed_data = np.empty_like(self.observed_data)
        self.hilbert_transforms = np.empty_like(self.observed_data)
        for int_ind in range(self.intensities.size):
            for x_ind in range(self.propagation_distance.size):
                self.abs_observed_data[int_ind][x_ind][:] = np.abs(self.observed_data[int_ind][x_ind][:])
                self.hilbert_transforms[int_ind][x_ind][:] = np.abs(hilbert(self.observed_data[int_ind][x_ind][:]))

    def plot_hilbert_transforms(self):
        x, y, z = np.shape(self.observed_data)
        fig, axes = plt.subplots(nrows=1, ncols=x)
        if x == 1:
            for ind in range(y):
                axes.plot(self.hilbert_transforms[0][y][:], color=color_spec[y])
            plt.show()

        else:
            for ind_x in range(x):
                axes[ind_x].grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
                axes[ind_x].set_yticks([])
                axes[ind_x].ticklabel_format(scilimits=(0, 0))
                for ind_y in range(y):
                    if ind_x == 0:
                        if ind_y == 0:
                            axes[ind_x].plot(self.mu_m, self.abs_observed_data[ind_x][ind_y][:]/np.max(self.abs_observed_data[ind_x][ind_y][:]), color='white', alpha = 0.4)
                        axes[ind_x].fill_between(self.mu_m, y1=(self.hilbert_transforms[ind_x][ind_y][:] / np.max(self.hilbert_transforms[ind_x][0][:])) - ind_y, y2=-ind_y, color=color_spec[ind_y], alpha=0.5, label='{:.3e}'.format(self.propagation_distance[ind_y] + 1000*self.dx))
                    else:
                        if ind_y == 0:
                            axes[ind_x].plot(self.mu_m, self.abs_observed_data[ind_x][ind_y][:]/np.max(self.abs_observed_data[ind_x][ind_y][:]), color='white', alpha = 0.4)
                        axes[ind_x].fill_between(self.mu_m, y1=(self.hilbert_transforms[ind_x][ind_y][:] / np.max(self.hilbert_transforms[ind_x][0][:])) - ind_y, y2=-ind_y, color=color_spec[ind_y], alpha=0.5)

                    plt.figlegend(loc='upper center', ncol=y, fancybox=False, shadow=False, edgecolor='black', title='Distanz in m')
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
        self.observed_data_mono = np.load(self.path + '/E_data_mono.npy')
        self.no_observer, self.timesteps = np.shape(self.observed_data)
        self.grid_information = np.load(self.path + '/info.npy')
        self.relative_observer_pos = np.load(self.path + '/relative_ind.npy')
        self.dt, self.dx, self.no_lambdas, self.peak_ts, self.pulse_duration = self.grid_information
        print("loaded in --- %s seconds ---" % (time.time() - start_time))

    def zero_pad(self, number_of_zeros):
        self.padded_observed_data = np.pad(self.observed_data[0:self.no_observer][:], [(0, 0), (number_of_zeros, number_of_zeros)], constant_values=0)
        self.padded_observed_data_mono = np.pad(self.observed_data_mono[0:self.no_observer][:], [(0, 0), (number_of_zeros, number_of_zeros)], constant_values=0)

    def set_fft_limits(self, past_from_max, future_from_max):
        self.windowed_data = np.zeros(shape=(self.no_observer, (past_from_max + future_from_max + 1)))
        self.windowed_data_mono = np.zeros(shape=(self.no_observer, (past_from_max + future_from_max + 1)))

        if self.padded_observed_data is not None:
            x, y = np.shape(self.padded_observed_data)
            self.analytical_signal = np.empty(shape=(self.no_observer, y), dtype='complex_')
            self.envelope_amplitude = np.zeros(shape=(self.no_observer, y))
            self.analytical_signal_mono = np.empty(shape=(self.no_observer, y), dtype='complex_')
            self.envelope_amplitude_mono = np.zeros(shape=(self.no_observer, y))
            for obs in range(self.no_observer):
                self.analytical_signal[obs] = hilbert(self.padded_observed_data[obs][:])
                self.analytical_signal_mono[obs] = hilbert(self.padded_observed_data_mono[obs][:])
                self.envelope_amplitude[obs] = np.abs(self.analytical_signal[obs][:])
                self.envelope_amplitude_mono[obs] = np.abs(self.analytical_signal_mono[obs][:])

        else:
            x, y = np.shape(self.observed_data)
            self.analytical_signal = np.empty(shape=(self.no_observer, y), dtype='complex_')
            self.analytical_signal_mono = np.empty(shape=(self.no_observer, y), dtype='complex_')
            self.envelope_amplitude = np.zeros(shape=(self.no_observer, y))
            self.envelope_amplitude_mono = np.zeros(shape=(self.no_observer, y))
            for obs in range(self.no_observer):
                self.analytical_signal[obs] = hilbert(self.observed_data[obs][:])
                self.analytical_signal_mono[obs] = hilbert(self.observed_data_mono[obs][:])
                self.envelope_amplitude[obs] = np.abs(self.analytical_signal[obs][:])
                self.envelope_amplitude_mono[obs] = np.abs(self.analytical_signal_mono[obs][:])


        self.max_timestep = np.zeros(self.no_observer)
        self.max_timestep = np.argmax(self.envelope_amplitude[0:self.no_observer][:], axis=1)
        self.max_timestep_mono = np.zeros(self.no_observer)
        self.max_timestep_mono = np.argmax(self.envelope_amplitude_mono[0:self.no_observer][:], axis=1)

        begin = np.array(self.max_timestep - past_from_max)
        end = np.array(self.max_timestep + future_from_max + 1)

        begin_mono = np.array(self.max_timestep_mono - past_from_max)
        end_mono = np.array(self.max_timestep_mono + future_from_max + 1)

        if self.padded_observed_data is not None:
            for obs in range(self.no_observer):
                self.windowed_data[obs][:] = self.padded_observed_data[obs][begin[obs]:end[obs]]
                self.windowed_data_mono[obs][:] = self.padded_observed_data_mono[obs][begin_mono[obs]:end_mono[obs]]

        else:
            for obs in range(self.no_observer):
                self.windowed_data[obs][:] = self.observed_data[obs][begin[obs]:end[obs]]
                self.windowed_data_mono[obs][:] = self.observed_data_mono[obs][begin_mono[obs]:end_mono[obs]]


    def fft(self):

        if self.windowed_data is not None:
            x, y = np.shape(self.windowed_data)
            self.timesteps_wind = y
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), self.timesteps_wind // 2)
            self.fft_matrix = 2 / self.timesteps_wind * np.fft.fft(self.windowed_data[0:self.no_observer][:])
            self.fft_matrix_mono = 2 / self.timesteps_wind * np.fft.fft(self.windowed_data_mono[0:self.no_observer][:])
            self.abs_fft = np.zeros(shape=(self.no_observer, self.timesteps_wind//2))
            self.abs_fft_mono = np.zeros(shape=(self.no_observer, self.timesteps_wind // 2))
            for obs in range(self.no_observer):
                self.abs_fft[obs][:] = np.abs(self.fft_matrix[obs][0:self.timesteps_wind//2])
                self.abs_fft_mono[obs][:] = np.abs(self.fft_matrix_mono[obs][0:self.timesteps_wind // 2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            self.abs_fft_sqrd_mono = self.abs_fft_mono[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

        elif self.padded_observed_data is not None:
            timesteps = self.timesteps + self.zero_padding
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps//2)
            self.fft_matrix = np.fft.fft(self.padded_observed_data[0:self.no_observer][:])
            self.fft_matrix_mono = np.fft.fft(self.padded_observed_data_mono[0:self.no_observer][:])
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps // 2))
            self.abs_fft_mono = np.zeros(shape=(self.no_observer, timesteps // 2))
            for obs in range(self.no_observer):
                self.abs_fft[obs][:] = 2 / timesteps * np.abs(self.fft_matrix[obs][0:timesteps//2])
                self.abs_fft_mono[obs][:] = 2 / timesteps * np.abs(self.fft_matrix_mono[obs][0:timesteps // 2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            self.abs_fft_sqrd_mono = self.abs_fft_mono[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

        else:
            timesteps = self.timesteps
            self.omega = 2 * np.pi * np.linspace(0, 1 / (2 * self.dt), timesteps // 2)
            self.fft_matrix = 2 / timesteps * np.fft.fft(self.observed_data[0:self.no_observer][:])
            self.fft_matrix_mono = 2 / timesteps * np.fft.fft(self.observed_data_mono[0:self.no_observer][:])
            self.abs_fft = np.zeros(shape=(self.no_observer, timesteps // 2))
            self.abs_fft_mono = np.zeros(shape=(self.no_observer, timesteps // 2))
            for obs in range(self.no_observer):
                self.abs_fft[obs][:] =  np.abs(self.fft_matrix[obs][0:timesteps//2])
                self.abs_fft_mono[obs][:] = np.abs(self.fft_matrix[obs][0:timesteps // 2])
            self.abs_fft_sqrd = self.abs_fft[:][:] ** 2
            self.abs_fft_sqrd_mono = self.abs_fft_mono[:][:] ** 2
            abs_fft_sqrd_max = np.max(self.abs_fft_sqrd[:][:])
            self.normalized_abs_fft_sqrd = self.abs_fft_sqrd / abs_fft_sqrd_max

    def visualize_over_frequencies(self):
        x = self.relative_observer_pos * self.dx
        y = self.omega /self.first_har
        X, Y = np.meshgrid(x, y)
        fig, axes = plt.subplots()
        c = axes.contourf(X, Y, np.transpose(self.normalized_abs_fft_sqrd), levels=20, locator=ticker.LogLocator(),
                          cmap='magma', norm=LogNorm())
        axes.set_xlabel('Distanz in m')
        axes.set_ylabel(r'Harmonische')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.colorbar(c, orientation='horizontal')
        plt.show()

    def visualize_n_over_length(self, number_of_harmonic):
        intensity = 5*10e12
        E_0 = np.sqrt(2*np.sqrt(mu0/eps0)*intensity)
        index_of_second = np.argmin(np.abs(self.omega-number_of_harmonic*self.first_har))
        index_of_main = np.argmin(np.abs(self.omega-self.first_har))
        first_amplitude = self.abs_fft[0][index_of_main]

        shg_amplitude = np.zeros(self.no_observer)
        shg_amplitude_mono = np.zeros(self.no_observer)
        for obs in range(self.no_observer):
            shg_amplitude[obs] = self.abs_fft[obs][index_of_second]
            shg_amplitude_mono[obs] = self.abs_fft_mono[obs][index_of_second]

        #analytical = Qpm_Length_Analytical(initial_values=[E_0, 0])
        analytical = Qpm_Length_Analytical(initial_values=[E_0, 0])
        A_1, A_2, z = analytical.integrate(end_z=2.55e-5, dz=0.0001e-6)

        analytical_mono = Qpm_Length_Analytical(initial_values=[E_0, 0], mode='mono')
        A_1_mono, A_2_mono, z_mono = analytical_mono.integrate(end_z=2.55e-5, dz=0.0001e-6)

        fig, axes = plt.subplots()
        for i in range(int(self.no_lambdas)*2 + 1):
            axes.axvline(x=i*738*self.dx, color='black', linestyle='dashed', alpha=1)
        axes.grid(True, linestyle=(0, (1, 5)), color=GREY, linewidth=1)
        axes.plot(z, np.abs(A_2)*231/(1.38*E_0), color=ORANGE)
        axes.plot(z, np.abs(A_2_mono)*231/(1.38*E_0), color=ORANGE, label='Theorie')
        #axes.plot(z, np.abs(A_2)/(1.38*first_amplitude), color=ORANGE)
        #axes.plot(z, np.abs(A_2_mono) / (1.38 * first_amplitude), color=ORANGE)
        axes.plot(self.relative_observer_pos * self.dx, shg_amplitude / (1.38 * first_amplitude), linewidth=1.5, color=TEAL,
                  linestyle='dashed', label='FDTD')
        axes.plot(self.relative_observer_pos * self.dx, shg_amplitude_mono / (1.38 * first_amplitude), linewidth=1.5,
                  color=TEAL,
                  linestyle='dashed')
        axes.set_xlim([0, 2.55e-05])
        axes.set_xlabel('Propagationsdistanz in m')
        axes.set_ylabel('Normalisierte zweite Harmonische')
        plt.figlegend(loc='upper center', ncol=2, fancybox=False, shadow=False, edgecolor='black')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.show()

    def visualize_windowed_data(self):
        fig, axes = plt.subplots()
        x, y = np.shape(self.windowed_data)
        for obs in range(600):
            axes.plot(np.arange(y), self.windowed_data[obs][:])
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


class Load_Lorentz_Slab(Case):

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
        axes[0][0].ticklabel_format(scilimits=(0, 0))
        axes[0][1].ticklabel_format(scilimits=(0, 0))
        axes[1][0].ticklabel_format(scilimits=(0, 0))
        axes[1][1].ticklabel_format(scilimits=(0, 0))
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



#soliton_test = Load_Soliton('late_night')
#soliton_test.plot_hilbert_transforms()


#lorentz_test = Load_Lorentz_Slab('new_data_files_per_terminal')
#lorentz_test.visualize()
# 51750 sec

'''
qpm_test_end_P = Load_QPM_end_P('P_end_new')
qpm_test_end_P.zero_pad(30000)
qpm_test_end_P.set_fft_limits(past_from_max=18000, future_from_max=18000)
qpm_test_end_P.show_E()
qpm_test_end_P.fft()
qpm_test_end_P.show_spectrum()
'''

qpm_test = Load_QPM_length('mono_and_not_mono')
qpm_test.zero_pad(37000)
qpm_test.set_fft_limits(past_from_max=18000, future_from_max=18000)
qpm_test.fft()
qpm_test.visualize_windowed_data()
#qpm_test.visualize_over_frequencies()
qpm_test.visualize_n_over_length(number_of_harmonic=2)









