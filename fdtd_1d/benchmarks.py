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