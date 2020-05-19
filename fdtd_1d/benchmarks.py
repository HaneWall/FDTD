import numpy as np
import fdtd_1d as f
import matplotlib.pyplot as plt
from .constants import c0

# Implementation of a dielectric slab

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
        self.d_ti = int(self.lamb_guided / (self.ti_n * 4 * self.dx))
        self.d_si = int(self.lamb_guided / (self.si_n * 4 * self.dx))
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
            self.theo_amplitude.append(theory_dielectric_slab(w_grid)[0])
            self.theo_phasenunterschied.append(theory_dielectric_slab(w_grid)[1])

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