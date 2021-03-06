from .backend import backend as bd
import numpy as np
import os
import time
import fdtd_1d as f
import matplotlib.pyplot as plt
from .constants import c0, BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY
from werkzeug.utils import cached_property
from multiprocessing import Pool






color_spec = [BLUE, CYAN, TEAL, ORANGE, RED, MAGENTA, GREY]

# Implementation of a dielectric slab
def theory_dielectric_slab_complex(grid):
    d = ((grid.materials[0].position[-1] - grid.materials[0].position[0] + 1)) * grid.dx
    omega = grid.sources[0].omega
    n_real = np.sqrt((np.abs(grid.materials[0].epsilon_complex(omega)) + grid.materials[0].epsilon_real(omega)) / 2)
    kappa = np.sqrt((np.abs(grid.materials[0].epsilon_complex(omega)) - grid.materials[0].epsilon_real(omega)) / 2)
    n = n_real + 1j*kappa
    k0 = grid.sources[0].omega / c0
    k = n*k0
    q = ((n - 1) ** 2) / ((n + 1) ** 2) * np.exp(2j * k * d)
    e_inc = grid.sources[0].ampl
    e_tr = e_inc * (2 / (n + 1)) * (2 * n / (n + 1)) * (1 / (1 - q)) * np.exp(1j * (k - k0) * d)
    theo_amplitude = np.abs(e_tr)
    theo_phasenunterschied = np.angle(e_tr)
    return theo_amplitude, theo_phasenunterschied

def theory_dielectric_slab(grid):
    d = ((grid.materials[0].position[-1] - grid.materials[0].position[0] + 1)) * grid.dx
    n = np.sqrt(grid.eps[grid.materials[0].position[1]])
    k0 = grid.sources[0].omega / c0
    k = n * k0
    q = ((n - 1) ** 2) / ((n + 1) ** 2) * np.exp(2j * k * d)
    e_inc = grid.sources[0].ampl
    e_tr = e_inc * (2 / (n + 1)) * (2 * n / (n + 1)) * (1 / (1 - q)) * np.exp(1j * (k - k0) * d)
    theo_amplitude = np.abs(e_tr)
    theo_phasenunterschied = np.angle(e_tr)
    return theo_amplitude, theo_phasenunterschied



class benchmark:
    '''Parent class for all benchmarks. Declares path/directory to store data'''


    def __init__(self, name, benchmark_type):
        self.name = name
        self.benchmark_type = benchmark_type
        self.dir_path = None
        self.grids = []

    def allocate_directory(self):

        self.dir_path = os.path.join(os.path.dirname(__file__), 'saved_data/'+self.benchmark_type+'/'+self.name)
        os.mkdir(path=self.dir_path)

    def store_obs_data(self):
        self.allocate_directory()
        for grid in self.grids:
            grid.store_obs_data(benchmark=True, benchmark_name=self.name)



class Harmonic_Slab_Lorentz_Setup(benchmark):


    def __init__(self,  name, dx, length_grid_in_dx, length_media_in_dx, start_index_media, wavelength, ampl, conductivity, eps_inf, gamma, w0, chi_1, chi_2, chi_3, timesteps, courant=1):
        super().__init__(name=name, benchmark_type='harmonic_lorentz_slab')
        self.start_media = start_index_media
        self.dx = dx
        self.indices = []
        self.grids = []

        for grid in range(len(self.dx)):
            self.indices.append([start_index_media + 2 + i for i in np.arange(0, length_media_in_dx[grid] - 1)])
            self.grids.append('grid_'+str(grid))

        self.eps_inf = eps_inf
        self.Nx = length_grid_in_dx
        self.length_media = np.array(length_media_in_dx)
        self.lamb = wavelength
        self.timesteps = timesteps
        self.courant = courant
        self.N_lambda = np.zeros(len(self.dx))                                                             #[[] for _ in range(len(self.dx))]
        self.ampl = ampl
        self.conductivity = conductivity
        self.gamma = gamma
        self.w0 = w0

        self.chi_1 = chi_1
        self.chi_2 = chi_2
        self.chi_3 = chi_3

        self.eps_real = None
        self.eps_imag = None
        self.eps_complex = None
        self.n_real = None

        self.wo_phase_merged = np.zeros(len(self.dx))
        self.theo_phasenunterschied_merged = np.zeros(shape=(len(self.dx), np.max(self.length_media)-1))
        self.theo_amplitude_merged = np.zeros(shape=(len(self.dx), np.max(self.length_media)-1))
        self.exp_phase_merged = np.zeros(shape=(len(self.dx), np.max(self.length_media)-1))
        self.exp_amplitude_merged = np.zeros(shape=(len(self.dx), np.max(self.length_media)-1))


    def _grid_wo_slab(self):
        position_src = self.start_media - 1
        for grid in range(len(self.dx)):
            position_obs = self.Nx[grid] - 3
            end_mur = self.Nx[grid] - 1
            wo_grid = f.Grid(self.Nx[grid], dx=self.dx[grid], courant=self.courant)
            if wo_grid.courant == 1:
                wo_grid[position_src] = f.ActivatedSinus(name='SinsquaredActivated', wavelength=self.lamb,
                                                         carrier_wavelength=(self.lamb * 30), phase_shift=0,
                                                         amplitude=self.ampl, tfsf=True)
            else:
                wo_grid[position_src] = f.ActivatedSinus(name='SinsquaredActivated', wavelength=self.lamb,
                                                         carrier_wavelength=(self.lamb * 30), phase_shift=0,
                                                         amplitude=self.ampl, tfsf=False)
            wo_grid[position_obs] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=self.timesteps[grid]-200)

            if wo_grid.courant == 0.5:
                wo_grid[0] = f.LeftSideGridBoundary()
                wo_grid[end_mur] = f.RightSideGridBoundary()
            else:
                wo_grid[0] = f.LeftSideMur()
                wo_grid[end_mur] = f.RightSideMur()

            wo_grid.run_timesteps(self.timesteps[grid], vis=False)
            self.wo_phase_merged[grid] = wo_grid.local_observers[0].phase

    def _grids_w_slab(self):
        position_src = self.start_media - 1
        for grid in range(len(self.dx)):
            position_obs = self.Nx[grid] - 3
            end_mur = self.Nx[grid] - 1
            for ind_media, ind_array in zip(self.indices[grid], range(self.length_media[grid])):
                # Step 1: init grid
                w_grid = 'slab' + str(ind_media)
                w_grid = f.Grid(nx=self.Nx[grid], dx=self.dx[grid], courant=self.courant)

                # Step 2: init media
                w_grid[self.start_media:ind_media] = f.LorentzMedium(name='media', permeability=1, eps_inf=self.eps_inf, conductivity=self.conductivity, gamma=self.gamma, chi_1=self.chi_1, chi_2=self.chi_2, chi_3=self.chi_3, w0=self.w0)

                # Step 3: init source
                if w_grid.courant == 1:
                    w_grid[position_src] = f.ActivatedSinus(name='SinsquaredActivated', wavelength=self.lamb,
                                                             carrier_wavelength=(self.lamb * 30), phase_shift=0,
                                                             amplitude=self.ampl, tfsf=True)
                else:
                    w_grid[position_src] = f.ActivatedSinus(name='SinsquaredActivated', wavelength=self.lamb,
                                                             carrier_wavelength=(self.lamb * 30), phase_shift=0,
                                                             amplitude=self.ampl, tfsf=False)

                # Step 4: add observer
                w_grid[position_obs] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=self.timesteps[grid]-200)

                # Step 5: add boundaries
                if w_grid.courant == 0.5:
                    w_grid[0] = f.LeftSideGridBoundary()
                    w_grid[end_mur] = f.RightSideGridBoundary()
                else:
                    w_grid[0] = f.LeftSideMur()
                    w_grid[end_mur] = f.RightSideMur()

                # Step 6: run simulation
                w_grid.run_timesteps(timesteps=self.timesteps[grid], vis=False)

                # Step 7: misc
                self.exp_amplitude_merged[grid][ind_array] = w_grid.local_observers[0].amplitude
                self.exp_phase_merged[grid][ind_array] = w_grid.local_observers[0].phase
                self.theo_amplitude_merged[grid][ind_array] = theory_dielectric_slab_complex(w_grid)[0]
                self.theo_phasenunterschied_merged[grid][ind_array] = theory_dielectric_slab_complex(w_grid)[1]

                # if list self.eps_real is empty:
                if self.eps_real is None:
                    self.eps_real = w_grid.materials[0].epsilon_real(w_grid.sources[0].omega)
                    self.eps_imag = w_grid.materials[0].epsilon_imag(w_grid.sources[0].omega)
                    self.eps_complex = w_grid.materials[0].epsilon_complex(w_grid.sources[0].omega)
                    self.n_real = np.sqrt((np.abs(self.eps_complex) + self.eps_real)/2)



    def _visualize(self):
        max_resolution_grid = int(np.argmax(self.N_lambda))

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

        axes[0][0].plot(np.array((np.array(self.indices[max_resolution_grid])) - self.start_media) * self.dx[max_resolution_grid],
                        self.theo_amplitude_merged[max_resolution_grid][0:self.length_media[max_resolution_grid] - 1], label='Theorie', color=ORANGE)
        axes[1][0].plot(np.array((np.array(self.indices[max_resolution_grid])) - self.start_media) * self.dx[max_resolution_grid],
                        self.theo_phasenunterschied_merged[max_resolution_grid][0:self.length_media[max_resolution_grid] - 1], color=ORANGE,
                        label='Theorie')

        for grid in range(len(self.dx)):
            axes[0][0].plot(np.array((np.array(self.indices[grid])) - self.start_media) * self.dx[grid], self.exp_amplitude_merged[grid][0:self.length_media[grid] - 1], color=color_spec[grid], linestyle='dashed', label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambda[grid]))
            axes[0][1].plot(np.array((np.array(self.indices[grid])) - self.start_media) * self.dx[grid], self.theo_amplitude_merged[grid][0:self.length_media[grid] - 1]/self.exp_amplitude_merged[grid][0:self.length_media[grid] - 1], color=color_spec[grid], linestyle='dashed', label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambda[grid]))
            axes[1][0].plot(np.array((np.array(self.indices[grid])) - self.start_media) * self.dx[grid], self.get_exp_phasedifference[grid][0:self.length_media[grid] - 1], color= color_spec[grid], linestyle='dashed', label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambda[grid]))
            axes[1][1].plot(np.array((np.array(self.indices[grid])) - self.start_media) * self.dx[grid], np.abs(self.get_exp_phasedifference[grid][0:self.length_media[grid] - 1] - self.theo_phasenunterschied_merged[grid][0:self.length_media[grid] - 1]), color=color_spec[grid], label=r'$N_{\lambda}=$' + '{0:.3}'.format(self.N_lambda[grid]))

        axes[0][0].legend(loc='best')
        axes[0][1].legend(loc='best')
        axes[1][0].legend(loc='best')
        axes[1][1].legend(loc='best')

        plt.show()

    @cached_property
    def get_exp_phasedifference(self):
        phase_diff = np.zeros(shape=(len(self.dx), np.max(self.length_media)-1))
        for grid in range(len(self.dx)):
            phase_diff[grid][0:self.length_media[grid] - 1] = -np.array(self.exp_phase_merged[grid][0:self.length_media[grid] - 1]) + self.wo_phase_merged[grid]
            mask = (phase_diff[grid][:]) > np.pi
            phase_diff[grid][mask] -= 2*np.pi
        return phase_diff

    def _set_N_lambda(self):
        for grid in range(len(self.dx)):
            self.N_lambda[grid] = self.lamb/(self.dx[grid]*self.n_real)

    def store_obs_data(self):
        self.allocate_directory()
        file_grid_informations = os.path.join(self.dir_path, 'info.npy')
        file_width_in_dx = os.path.join(self.dir_path, 'width.npy')
        file_theory_ampl = os.path.join(self.dir_path, 'theory_ampl.npy')
        file_theory_phase = os.path.join(self.dir_path, 'theory_phase.npy')
        file_exp_ampl = os.path.join(self.dir_path, 'exp_ampl.npy')
        file_exp_phase = os.path.join(self.dir_path, 'exp_phase.npy')

        grid_informations = np.zeros(shape=(len(self.grids), 4))
        width_in_dx = np.zeros(shape=(len(self.grids), np.max(self.length_media)))
        for grid in range(len(self.grids)):
            width_in_dx[grid][0:self.length_media[grid] - 1] = np.array(self.indices[grid]) - self.start_media
            grid_informations[grid][0] = self.dx[grid]
            grid_informations[grid][1] = self.timesteps[grid]
            grid_informations[grid][2] = self.N_lambda[grid]
            grid_informations[grid][3] = self.length_media[grid]

        np.save(file_grid_informations, arr=grid_informations)
        np.save(file_width_in_dx, arr=width_in_dx)
        np.save(file_theory_ampl, arr=self.theo_amplitude_merged)
        np.save(file_theory_phase, arr=self.theo_phasenunterschied_merged)
        np.save(file_exp_ampl, arr=self.exp_amplitude_merged)
        np.save(file_exp_phase, arr=self.get_exp_phasedifference)

    def run_benchmark(self):
        start_time = time.time()
        self._grid_wo_slab()
        self._grids_w_slab()
        self._set_N_lambda()
        print("computed in --- %s seconds ---" % (time.time() - start_time))
        self._visualize()


class QPM_Length_multiple_chi_2(benchmark):
    def __init__(self, number_of_lambdas, timesteps, name, peak_timestep, pulse_duration, number_of_distributed_observer, chi2, multi=True):
        super().__init__(name=name, benchmark_type='qpm_harmonic_length_multiple')
        self.no_of_lambdas = number_of_lambdas
        self.half_qpm = 737
        self.chi_2 = chi2
        self.multi = multi
        self.dx = 4e-09
        self.dt = 0.5*self.dx/c0
        self.timesteps = timesteps
        self.nx = self.no_of_lambdas * (2*self.half_qpm + 1) + 10           # 10 is just a buffer to place boundarys and src
        self.start_media = 5
        self.ending_indices = np.array([self.start_media + i*self.half_qpm for i in range(self.no_of_lambdas*2 + 1)])
        self.peak_timestep = peak_timestep
        self.pulse_duration = pulse_duration
        self.no_observer = number_of_distributed_observer
        self.obs_distance = (self.nx - 10)//self.no_observer
        self.obs_positions = np.array([5 + i*self.obs_distance for i in range(self.no_observer)])

    def _allocate_memory(self):
        self.grid_information = np.array([self.dt, self.dx, self.no_of_lambdas, self.peak_timestep, self.pulse_duration])
        self.relative_observer_pos = self.obs_positions - self.start_media

    def store_obs_data(self):
        start_time = time.time()
        self.allocate_directory()
        self._allocate_memory()

        file_grid_info = os.path.join(self.dir_path, 'info.npy')
        file_relative_pos = os.path.join(self.dir_path, 'relative_ind.npy')
        file_data = os.path.join(self.dir_path, 'E_data.npy')
        file_chi_info = os.path.join(self.dir_path, 'chi2.npy')

        np.save(file_grid_info, arr=self.grid_information)
        np.save(file_relative_pos, arr=self.relative_observer_pos)
        np.save(file_data, arr=self.observed_data)
        np.save(file_chi_info, arr=self.chi_2)
        print("saved in --- %s seconds ---" % (time.time() - start_time))

    def _create_grids(self, chi2):
        print('process initiated')
        qpm_grid = f.Grid(nx=self.nx, dx=self.dx, benchmark='qpm_harmonic_length_multiple', courant=0.5)

        for indices in range(len(self.ending_indices) - 1):
            if indices % 2 == 0:
                    qpm_grid[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                        name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0],
                        chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

            else:
                    qpm_grid[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                        name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[chi2, 0, 0],
                        chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

        qpm_grid[3] = f.GaussianImpulseWithFrequency(name='Varin', Intensity=5*10e12, wavelength=1.064e-06, pulse_duration=self.pulse_duration, peak_timestep=self.peak_timestep, tfsf=False)

        for pos in self.obs_positions:
            qpm_grid[int(pos)] = f.E_FFTObserver(name='Varin', first_timestep=0, second_timestep=self.timesteps - 1)


        qpm_grid[0] = f.LeftSideGridBoundary()
        qpm_grid[self.nx - 1] = f.RightSideGridBoundary()

        # step 6: run simulation
        qpm_grid.run_timesteps(timesteps=self.timesteps, vis=False)

        observed_grid_data = bd.zeros(shape=(self.no_observer, self.timesteps))
        for obs in range(self.no_observer):
            observed_grid_data[obs] = qpm_grid.local_observers[obs].observed_E

        return observed_grid_data

    def run_benchmark(self):
        start_time = time.time()
        processpool = Pool()
        self.observed_data = bd.zeros(shape=(len(self.chi_2), self.no_observer, self.timesteps))
        if self.multi:
            self.observed_data = bd.stack(processpool.map(self._create_grids, self.chi_2))
        else:
            for ind_chi, chi in zip(range(len(self.chi_2)), self.chi_2):
                self.observed_data[ind_chi] = self._create_grids(self.chi_2[ind_chi])
        print("computed in --- %s seconds ---" % (time.time() - start_time))


class QPM_Length(benchmark):

    def __init__(self, number_of_lambdas, timesteps, name, peak_timestep, pulse_duration, number_of_distributed_observer):
        super().__init__(name=name, benchmark_type='qpm_harmonic_length')
        self.no_of_lambdas = number_of_lambdas
        self.half_qpm = 737
        self.dx = 4e-09
        self.timesteps = timesteps
        self.nx = self.no_of_lambdas * (2*self.half_qpm + 1) + 10           # 10 is just a buffer to place boundarys and src
        self.start_media = 5
        self.ending_indices = np.array([self.start_media + i*self.half_qpm for i in range(self.no_of_lambdas*2 + 1)])
        self.peak_timestep = peak_timestep
        self.pulse_duration = pulse_duration
        self.no_observer = number_of_distributed_observer
        self.obs_distance = (self.nx - 10)//self.no_observer
        self.obs_positions = np.array([5 + i*self.obs_distance for i in range(self.no_observer)])

    def _create_grid(self):
        qpm_grid = f.Grid(nx=self.nx, dx=self.dx, benchmark='qpm_harmonic_length', courant=0.5)
        self.grids.append(qpm_grid)

        for indices in range(len(self.ending_indices) - 1):
            if indices % 2 == 0:
                    qpm_grid[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                        name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0],
                        chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

            else:
                    qpm_grid[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                        name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[-30.e-12, 0, 0],
                        chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

        qpm_grid[3] = f.GaussianImpulseWithFrequency(name='Varin', Intensity=5*10e12, wavelength=1.064e-06, pulse_duration=self.pulse_duration, peak_timestep=self.peak_timestep, tfsf=False)

        for pos in self.obs_positions:
            qpm_grid[int(pos)] = f.E_FFTObserver(name='Varin', first_timestep=0, second_timestep=self.timesteps - 1)


        qpm_grid[0] = f.LeftSideGridBoundary()
        qpm_grid[self.nx - 1] = f.RightSideGridBoundary()

        #qpm_grid[0] = f.LeftSideMur()
        #qpm_grid[self.nx - 1] = f.RightSideMur()

        # step 6: run simulation
        qpm_grid.run_timesteps(timesteps=self.timesteps, vis=False)

    def _create_grid_mono(self):
        qpm_grid_mono = f.Grid(nx=self.nx, dx=self.dx, benchmark='qpm_harmonic_length', courant=0.5)
        self.grids.append(qpm_grid_mono)

        for indices in range(len(self.ending_indices) - 1):
            if indices % 2 == 0:
                    qpm_grid_mono[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                        name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0],
                        chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

            else:
                    qpm_grid_mono[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                        name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0],
                        chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

        qpm_grid_mono[3] = f.GaussianImpulseWithFrequency(name='Varin', Intensity=5*10e12, wavelength=1.064e-06, pulse_duration=self.pulse_duration, peak_timestep=self.peak_timestep, tfsf=False)

        for pos in self.obs_positions:
            qpm_grid_mono[int(pos)] = f.E_FFTObserver(name='Varin', first_timestep=0, second_timestep=self.timesteps - 1)


        qpm_grid_mono[0] = f.LeftSideGridBoundary()
        qpm_grid_mono[self.nx - 1] = f.RightSideGridBoundary()

        #qpm_grid[0] = f.LeftSideMur()
        #qpm_grid[self.nx - 1] = f.RightSideMur()

        # step 6: run simulation
        qpm_grid_mono.run_timesteps(timesteps=self.timesteps, vis=False)


    def _allocate_memory(self):
        self.observed_data = np.zeros(shape=(self.no_observer, self.timesteps))
        self.observed_data_mono = np.zeros(shape=(self.no_observer, self.timesteps))
        observer_object_list = np.array(self.grids[0].local_observers)
        observer_object_list_mono = np.array(self.grids[1].local_observers)
        for obs_ind in range(self.no_observer):
            self.observed_data[obs_ind][:] = observer_object_list[obs_ind].observed_E[:]
            self.observed_data_mono[obs_ind][:] = observer_object_list_mono[obs_ind].observed_E[:]
        self.grid_information = np.array([self.grids[0].dt, self.grids[0].dx, self.no_of_lambdas, self.peak_timestep, self.pulse_duration])
        self.relative_observer_pos = self.obs_positions - self.start_media

    def store_obs_data(self):
        start_time = time.time()
        self.allocate_directory()
        self._allocate_memory()

        file_grid_info = os.path.join(self.dir_path, 'info.npy')
        file_relative_pos = os.path.join(self.dir_path, 'relative_ind.npy')
        file_data = os.path.join(self.dir_path, 'E_data.npy')
        file_data_mono = os.path.join(self.dir_path, 'E_data_mono.npy')

        np.save(file_grid_info, arr=self.grid_information)
        np.save(file_relative_pos, arr=self.relative_observer_pos)
        np.save(file_data, arr=self.observed_data)
        np.save(file_data_mono, arr=self.observed_data_mono)
        print("saved in --- %s seconds ---" % (time.time() - start_time))

    def run_benchmark(self):
        start_time = time.time()
        self._create_grid()
        self._create_grid_mono()
        print("computed in --- %s seconds ---" % (time.time() - start_time))


class QPM_end_P(benchmark):
    ''' reproduces paper QuasiPhaseMatching from Varin's Paper '''

    def __init__(self, number_of_lambdas, timesteps, name, peak_timestep, pulse_duration):
        super().__init__(name=name, benchmark_type='qpm_harmonic')
        self.no_of_lambdas = number_of_lambdas
        # Varin parameters
        self.half_qpm = 737
        self.dx = 4e-09
        self.timesteps = timesteps
        self.nx = number_of_lambdas * (2*self.half_qpm + 1) + 10
        self.start_media = 5
        self.ending_indices = np.array([self.start_media + i*self.half_qpm for i in range(self.no_of_lambdas*2 + 1)])
        self.peak_timestep = peak_timestep
        self.pulse_duration = pulse_duration


    def _create_grid(self):
        self.position_P_obs = self.nx - 9
        self.position_E_obs = self.nx - 8

        # step 1: init grid
        qpm_grid = f.Grid(nx=self.nx, dx=self.dx, benchmark='qpm_harmonic', courant=0.5)
        self.grids.append(qpm_grid)

        for indices in range(len(self.ending_indices) - 1):
            if indices % 2 == 0:
                qpm_grid[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                    name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0],
                    chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

            else:
                qpm_grid[self.ending_indices[indices]:self.ending_indices[indices + 1]] = f.LorentzMedium(
                    name='Varin', permeability=1, eps_inf=1.0, chi_1=[2.42, 9.65, 1.46], chi_2=[-30.e-12, 0, 0],
                    chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

        qpm_grid[3] = f.GaussianImpulseWithFrequency(name='Varin', Intensity=10e12, wavelength=1.064e-06,
                                                     pulse_duration=self.pulse_duration,
                                                     peak_timestep=self.peak_timestep, tfsf=True)


        qpm_grid[self.position_E_obs] = f.E_FFTObserver(name='Varin', first_timestep=0, second_timestep=self.timesteps - 1)
        qpm_grid[self.position_P_obs] = f.P_FFTObserver(name='Varin', first_timestep=0, second_timestep=self.timesteps - 1)

        qpm_grid[0] = f.LeftSideGridBoundary()
        qpm_grid[self.nx - 1] = f.RightSideGridBoundary()

        #qpm_grid[0] = f.LeftSideMur()
        #qpm_grid[self.nx - 1] = f.RightSideMur()

        # step 6: run simulation
        qpm_grid.run_timesteps(timesteps=self.timesteps, vis=True)

    def _allocate_memory(self):
        self.observed_data = np.zeros(shape=(2, self.timesteps))
        self.observed_data[0][:] = self.grids[0].local_observers[0].observed_E[:]
        self.observed_data[1][:] = self.grids[0].local_observers[1].observed_P[:]
        self.grid_information = np.array([self.grids[0].dt, self.grids[0].dx, self.no_of_lambdas, self.peak_timestep,
                                          self.pulse_duration])


    def store_obs_data(self):
        start_time = time.time()
        self.allocate_directory()
        self._allocate_memory()

        self.relative_observer_pos = np.array([self.position_E_obs - self.start_media, self.position_P_obs - self.start_media])

        file_grid_info = os.path.join(self.dir_path, 'info.npy')
        file_relative_pos = os.path.join(self.dir_path, 'relative_ind.npy')
        file_data = os.path.join(self.dir_path, 'E_P_data.npy')

        np.save(file_grid_info, arr=self.grid_information)
        np.save(file_relative_pos, arr=self.relative_observer_pos)
        np.save(file_data, arr=self.observed_data)
        print("saved in --- %s seconds ---" % (time.time() - start_time))

    def run_benchmark(self):
        start_time = time.time()
        self._create_grid()
        print("computed in --- %s seconds ---" % (time.time() - start_time))


class Soliton(benchmark):

    ''' tries to show the majestic combined effect of GVD and SPM '''

    def __init__(self, name, central_wavelength, pulse_duration, intensities, x_to_snapshot, peak_timestep, frame_width_in_dx, dx, multi=True):
        super().__init__(name=name, benchmark_type='soliton')
        self.central_wavelength = central_wavelength
        self.pulse_duration = pulse_duration
        self.intensities = np.array(intensities)
        self.x_to_snapshot = x_to_snapshot
        self.peak_timestep = peak_timestep
        self.name = name
        self.dx = dx
        self.nx = int(x_to_snapshot[-1]/self.dx) + 5020
        self.frame_width = frame_width_in_dx
        self.multi = multi


    def _allocate_memory(self):
        self.grid_information = np.array([self.peak_timestep, self.pulse_duration, self.dx, self.frame_width])
        self.used_propagations = np.array(self.x_to_snapshot)
        self.used_intensities = np.array(self.intensities)

    def _create_grids(self, intensity):
        soliton_grid = f.Grid(self.nx, dx=self.dx, courant=0.5) #12mm 550000

        # Step 1: init media
        soliton_grid[5:self.nx-10] = f.CentroRamanMedium(name='Varin', chi_1=[0.69617, 0.40794, 0.89748], w0=[2.7537e16, 1.6205e16, 1.9034e14], chi_3=[1.94e-22, 0, 0], alpha=[0.7, 0, 0], wr=[8.7722e13, 0, 0], gamma_K=[0, 0, 0], gamma_R=[3.1250e13, 0, 0], permeability=1, conductivity=0, eps_inf=1)

        # Step 2: init src
        soliton_grid[3] = f.SechEnveloped(name='Varin', wavelength=1.5e-06, pulse_duration=self.pulse_duration, Intensity=intensity, peak_timestep=self.peak_timestep, tfsf=False)

        # Step 3: init frame
        soliton_grid[5:(6+self.frame_width)] = f.MovingFrame(x_to_snapshot=self.x_to_snapshot, central_wavelength=self.central_wavelength)

        # Step 4: init boundaries
        soliton_grid[0] = f.LeftSideGridBoundary()
        soliton_grid[self.nx - 1] = f.RightSideGridBoundary()

        soliton_grid.local_observers[0]._allocate_memory()
        timesteps = soliton_grid.local_observers[0].timesteps_to_store[-1] + 5000

        # Step 5: start benchmark
        soliton_grid.run_timesteps(timesteps, vis=False)

        # due to multiprocessing it is more convenient to store data by process (nobody loves waiting)
        observed_process_data = soliton_grid.local_observers[0].stored_data
        return observed_process_data

    def store_obs_data(self):
        start_time = time.time()
        self._allocate_memory()
        self.allocate_directory()

        file_grid_info = os.path.join(self.dir_path, 'info.npy')
        file_observed_data = os.path.join(self.dir_path, 'Int_Pos_E.npy')
        file_used_propagations = os.path.join(self.dir_path, 'propagations.npy')
        file_used_intensities = os.path.join(self.dir_path, 'intensities.npy')

        np.save(file_grid_info, arr=self.grid_information)
        np.save(file_used_intensities, arr=self.used_intensities)
        np.save(file_observed_data, arr=self.observed_data)
        np.save(file_used_propagations, arr=self.used_propagations)
        print("stored in --- %s seconds ---" % (time.time() - start_time))

    def run_benchmark(self):
        start_time = time.time()
        if self.multi:
            processpool = Pool()
            self.observed_data = np.array(processpool.map(self._create_grids, self.intensities))
        else:
            self.observed_data = np.zeros(shape=(np.size(self.intensities), len(self.x_to_snapshot), self.frame_width + 1))
            for intensity_index, intensity in zip(range(np.size(self.intensities)), self.intensities):
                self.observed_data[intensity_index] = self._create_grids(intensity)
        print("computed in --- %s seconds ---" % (time.time() - start_time))


class TiO2_Si02_Dielectric_Mirror_Setup:

    def __init__(self, N_lambda_media, wavelength_guided_for, wavelength, ampl, timesteps, number_of_layer_pairs, vary_layers=False, vary_inc_wavelength=False):
        self.N_lambda = N_lambda_media
        self.layer_number = [0 + i for i in range(number_of_layer_pairs)]
        self.number_of_layer = number_of_layer_pairs
        self.ampl = ampl
        self.timesteps = timesteps
        self.lamb = wavelength
        self.lamb_guided = wavelength_guided_for
        self.vary_layers = vary_layers
        self.vary_lambda = vary_inc_wavelength
        self.ti_n = 2.519                                                               # refractiveindex.info
        self.si_n = 1.453                                                               # refractiveindex.info
        self.dx = self.lamb/(self.ti_n * self.N_lambda)
        self.d_ti = int(self.lamb_guided / (self.ti_n * 4 * self.dx))                   # huge problem
        self.d_si = int(self.lamb_guided / (self.si_n * 4 * self.dx))                   # huge problem
        self.Nx = self.number_of_layer*(self.d_si + self.d_ti) + 14
        self.starting_locations_ti = [8 + i*(self.d_ti + self.d_si) for i in self.layer_number]
        self.starting_locations_si = np.array(self.starting_locations_ti) + self.d_ti
        self.position_src = 6
        self.position_obs = 3
        self.incident_wavelengths = [self.lamb + (i / 101) * self.lamb_guided for i in np.arange(0, 101, 1)]
        self.refl_ampl = []
        self.theory_R = []

    def _construct_non_vary_grid(self):
        end_mur = self.Nx - 1
        grid_non_vary = f.Grid(dx=self.dx, nx=self.Nx)
        grid_non_vary[0] = f.LeftSideMur()
        grid_non_vary[end_mur] = f.RightSideMur()
        grid_non_vary[self.position_src] = f.ActivatedSinus(name='Laser', wavelength=self.lamb, carrier_wavelength=6000e-09, tfsf=True, amplitude=self.ampl, phase_shift=0)
        grid_non_vary[self.position_obs] = f.QuasiHarmonicObserver(name='Observer', first_timestep=self.timesteps - 200)

        for layer_pair in self.layer_number:
            ti_start = self.starting_locations_ti[layer_pair]
            si_start = self.starting_locations_si[layer_pair]
            grid_non_vary[ti_start:si_start] = f.NonDispersiveMedia('TiO2', permittivity=self.ti_n**2, permeability=1, conductivity=0)
            grid_non_vary[si_start:(si_start+self.d_si)] = f.NonDispersiveMedia('SiO2', permittivity=self.si_n**2, permeability=1, conductivity=0)

        grid_non_vary.run_timesteps(self.timesteps)
        self.refl_ampl.append(grid_non_vary.local_observers[0].amplitude)

    def _construct_vary_grid(self):
        for layer in self.layer_number:
            Nx = (layer + 1) * (self.d_si + self.d_ti) + 12
            end_mur = Nx - 1
            grid_vary = f.Grid(dx=self.dx, nx=Nx)
            grid_vary[end_mur] = f.RightSideMur()
            grid_vary[0] = f.LeftSideMur()
            grid_vary[self.position_obs] = f.QuasiHarmonicObserver(name='Observer', first_timestep=self.timesteps - 300)
            grid_vary[self.position_src] = f.ActivatedSinus(name='Laser', wavelength=self.lamb, carrier_wavelength=6000e-09, tfsf=True, amplitude=self.ampl, phase_shift=0)

            for layers in range(0, layer + 1):
                ti_start = self.starting_locations_ti[layers]
                si_start = self.starting_locations_si[layers]
                grid_vary[ti_start:si_start] = f.NonDispersiveMedia('TiO2', permittivity=self.ti_n**2, permeability=1, conductivity=0)
                grid_vary[si_start:(si_start + self.d_si)] = f.NonDispersiveMedia('SiO2', permittivity=self.si_n**2, permeability=1, conductivity=0)

            grid_vary.run_timesteps(self.timesteps, vis=False)
            self.refl_ampl.append(grid_vary.local_observers[0].amplitude)
            self.theory_R.append(((self.si_n**(2*(layer+1))-self.ti_n**(2*(layer+1)))/(self.si_n**(2*(layer+1))+self.ti_n**(2*(layer+1))))**2)

    def _construct_vary_wavelength_grid(self):
        for wavelength in self.incident_wavelengths:
            end_mur = self.Nx - 1
            grid_vary_lamb = f.Grid(dx=self.dx, nx=self.Nx)
            grid_vary_lamb[0] = f.LeftSideMur()
            grid_vary_lamb[end_mur] = f.RightSideMur()
            grid_vary_lamb[self.position_src] = f.ActivatedSinus('Laser', wavelength=wavelength, carrier_wavelength=6000e-09, tfsf=True, amplitude=self.ampl, phase_shift=0)
            grid_vary_lamb[self.position_obs] = f.QuasiHarmonicObserver(name='Observer', first_timestep=self.timesteps - 200)

            for layer_pair in self.layer_number:
                ti_start = self.starting_locations_ti[layer_pair]
                si_start = self.starting_locations_si[layer_pair]
                grid_vary_lamb[ti_start:si_start] = f.NonDispersiveMedia('TiO2', permittivity=self.ti_n**2, permeability=1, conductivity=0)
                grid_vary_lamb[si_start:(si_start + self.d_si)] = f.NonDispersiveMedia('SiO2', permittivity=self.si_n**2, permeability=1, conductivity=0)

            grid_vary_lamb.run_timesteps(self.timesteps, vis=False)
            self.refl_ampl.append(grid_vary_lamb.local_observers[0].amplitude)

    def _visualize_vary_grid(self):
        fig, axes = plt.subplots(1, 1)
        axes.plot(np.array(self.layer_number)+1, np.array(self.refl_ampl)**2, linestyle='dashed', color='blue', marker='o', alpha=0.5, label='FDTD')
        axes.plot(np.array(self.layer_number)+1, np.array(self.theory_R), linestyle='dashed', color='red', marker='o', alpha=0.5, label='model')
        axes.legend(loc='best')
        axes.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes.set_xlabel('Anzahl der Paare', fontsize=14)
        axes.set_ylabel(r'Reflektionsgrad $\mathcal{R}$', fontsize=14)
        plt.show()

    def _visualize_vary_lamb(self):
        fig, axes = plt.subplots(1, 1)
        axes.plot(np.array(self.incident_wavelengths)*10**9, np.array(self.refl_ampl)**2, linestyle='dashed', color='blue', marker='o',
                  alpha=0.5, label='FDTD')
        axes.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes.set_xlabel(r'$\lambda$ in nm', fontsize=14)
        axes.set_ylabel(r'Reflektionsgrad $\mathcal{R}$', fontsize=14)
        plt.show()

    def run_benchmark(self):
        if self.vary_layers:
            self._construct_vary_grid()
            self._visualize_vary_grid()

        elif self.vary_lambda:
            self._construct_vary_wavelength_grid()
            self._visualize_vary_lamb()

        else:
            self._construct_non_vary_grid()


class Harmonic_Slab_Setup:

    def __init__(self, dx,  length_grid_in_dx, length_media_in_dx, start_index_media, wavelength, epsilon, ampl, timesteps):
        self.indices = [start_index_media + 2 + i for i in np.arange(0, length_media_in_dx - 1)]
        self.start_media = start_index_media
        self.dx = dx
        self.eps = epsilon
        self.Nx = length_grid_in_dx
        self.length_media = length_media_in_dx
        self.lamb = wavelength
        self.ampl = ampl
        self.wo_phase = []
        self.theo_phasenunterschied = []
        self.theo_amplitude = []
        self.exp_phase = []
        self.exp_amplitude = []
        self.timesteps = timesteps

    def _grid_wo_slab(self):
        position_src = self.start_media - 1
        position_obs = self.Nx - 3
        end_mur = self.Nx - 1
        wo_grid = f.Grid(self.Nx, dx=self.dx)
        wo_grid[position_src] = f.ActivatedSinus(name='SinsquaredActivated', wavelength=self.lamb, carrier_wavelength=10000.e-09, phase_shift=0, amplitude=self.ampl, tfsf=True)
        wo_grid[position_obs] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=self.timesteps-200)
        wo_grid[0] = f.LeftSideMur()
        wo_grid[end_mur] = f.RightSideMur()
        wo_grid.run_timesteps(self.timesteps, vis=False)
        self.wo_phase.append(wo_grid.local_observers[0].phase)

    def _grids_w_slab(self):
        position_src = self.start_media - 1
        position_obs = self.Nx - 3
        end_mur = self.Nx - 1
        for ind in self.indices:
            # Step 1: init grid
            w_grid = 'slab' + str(ind)
            w_grid = f.Grid(nx=self.Nx, dx=self.dx)

            # Step 2: init media
            w_grid[self.start_media:ind] = f.NonDispersiveMedia(name='media', permeability=1, permittivity=self.eps, conductivity=0)

            # Step 3: init source
            w_grid[position_src] = f.ActivatedSinus(name='SinsquaredActivated', wavelength=self.lamb, carrier_wavelength=10000.e-09, phase_shift=0, amplitude=self.ampl, tfsf=True)

            # Step 4: add observer
            w_grid[position_obs] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=self.timesteps-200)

            # Step 5: add boundaries
            w_grid[0] = f.LeftSideMur()
            w_grid[end_mur] = f.RightSideMur()

            # Step 6: run simulation
            w_grid.run_timesteps(timesteps=self.timesteps, vis=False)

            # Step 7: misc
            self.exp_amplitude.append(w_grid.local_observers[0].amplitude)
            self.exp_phase.append(w_grid.local_observers[0].phase)
            self.theo_amplitude.append(theory_dielectric_slab_complex(w_grid)[0])
            self.theo_phasenunterschied.append(theory_dielectric_slab_complex(w_grid)[1])

    def _visualize(self):
        fig, axes = plt.subplots(2, 2)
        fig.suptitle(r'$\epsilon_r=$ {epsilon}'.format(epsilon=self.eps) + r'   $N_{\lambda_{media}}=$' + '{0:.3}'.format(self.lamb/(self.dx*np.sqrt(self.eps))), fontsize=20)
        axes[0][0].plot(np.array(self.indices) - self.start_media, np.array(self.theo_amplitude), label='theorie', color='blue', marker='o', alpha=0.5)
        axes[0][0].grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes[0][0].plot(np.array(self.indices) - self.start_media, np.array(self.exp_amplitude), label='FDTD', linestyle='dashed', color='red', marker='s', alpha=0.5)
        axes[0][0].legend(loc='best')
        axes[0][0].set_xlabel('Breite des Mediums in ' + r'$\Delta_x$', fontsize=14)
        axes[0][0].set_ylabel('Transmittierte Amplitude ' + r'$Ez_{tr}$', fontsize=14)
        axes[0][0].set_xlim([0, self.length_media + 1])
        axes[0][1].plot(np.array(self.indices) - self.start_media, np.array(self.theo_amplitude) / np.array(self.exp_amplitude), color='black')
        axes[0][1].set_ylabel(r'$E_{tr,theo}$ / $E_{tr,FDTD}$', fontsize=14)
        axes[0][1].set_xlabel('Breite des Mediums in ' + r'$\Delta_x$', fontsize=14)
        axes[0][1].grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes[0][1].set_xlim([0, self.length_media + 1])
        axes[1][0].set_ylabel('Phasenunterschied', fontsize=14)
        axes[1][0].plot(np.array(self.indices) - self.start_media, self.theo_phasenunterschied, label='theorie', color='blue', alpha=0.5)
        axes[1][0].set_xlabel('Breite des Mediums in ' + r'$\Delta_x$', fontsize=14)
        axes[1][0].grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes[1][0].plot(np.array(self.indices) - self.start_media, -np.array(self.exp_phase) + self.wo_phase, color='red', linestyle='dashed',
                        label='FDTD', alpha=0.5)
        axes[1][0].set_xlim([0, self.length_media + 1])
        axes[1][0].legend()
        axes[1][1].set_xlabel('Breite des Mediums in ' + r'$\Delta_x$', fontsize=14)
        axes[1][1].set_ylabel(r'$d(\phi_{exp},\phi_{theo})$', fontsize=14)
        axes[1][1].plot(np.array(self.indices) - self.start_media,
                        np.abs(-np.array(self.exp_phase) + self.wo_phase - np.array(self.theo_phasenunterschied)), color='black')
        axes[1][1].grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes[1][1].set_xlim([0, self.length_media + 1])
        plt.show()

    def run_benchmark(self):
        self._grid_wo_slab()
        self._grids_w_slab()
        self._visualize()






