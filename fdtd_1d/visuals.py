import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as ani
from matplotlib.patches import Rectangle
import numpy as np

from .observer import QuasiHarmonicObserver
from .observer import FFTObserver
from .constants import mu0, eps0, c0



def visualize(Grid):
    fig, axes = plt.subplots(2, 1, dpi=100, figsize=(12,8))
    axes[0].plot(Grid.Ez)
    axes[0].xaxis.grid(linestyle='dotted')
    axes[0].set_xlabel('Cell')
    axes[0].set_ylabel(r'$E_z$', fontsize=12, rotation=0)
    axes[0].set_title('time passed: {0:.3e}s,  timesteps passed: {timesteps}'.format(Grid.time_passed,
                                                                                timesteps=Grid.timesteps_passed))
    axes[0].set_ylim([-1.5, 1.5])
    axes[1].plot(Grid.Hy)
    axes[1].xaxis.grid(linestyle='dotted')
    axes[1].set_xlabel('Cell')
    axes[1].set_ylabel(r'$H_y$', fontsize=12, rotation=0)
    for mat in Grid.materials:
        media_repr_0 = Rectangle(xy=(mat.position[0]-0.5, -1.4), height=2.8, width=(mat.position[-1] - mat.position[0] + 1),
                               color='grey', fill=True, alpha=mat.eps * 0.12)
        axes[0].add_patch(media_repr_0)
        axes[0].annotate(
            s=r'$\epsilon_r$={0:.2f}'.format(mat.eps) + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity),
            xy=(media_repr_0.get_x() + 0.1, media_repr_0.get_y() + 0.2), color='black')
        media_repr_1 = Rectangle(xy=(mat.position[0] - 0.5, -1.4), height=2.8, width=(mat.position[-1] - mat.position[0] + 1),
                                 color='grey', fill=True, alpha=mat.eps * 0.12)
        axes[1].add_patch(media_repr_1)
        if mat.model=='Lorentz':
            s = r'$\epsilon(\omega)$' + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity)
        else:
            s = r'$\epsilon_r$={0:.2f}'.format(mat.eps) + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity)
        axes[1].annotate(s, xy=(media_repr_1.get_x() + 0.1, media_repr_1.get_y() + 0.2), color='black')
    for src in Grid.sources:
        src_repr = Rectangle(xy=(src.position - 0.5, -src.ampl), height=2 * src.ampl, width=1, color='red', alpha=0.3)
        axes[0].add_patch(src_repr)

    for obs in Grid.local_observers:
        if isinstance(obs, QuasiHarmonicObserver):
            obs_repr = Rectangle(xy=(obs.position - 0.5, -1.4), height=2.8, width=1, color='green', alpha=0.3)
        elif isinstance(obs, FFTObserver):
            obs_repr = Rectangle(xy=(obs.position - 0.5, -1.4), height=2.8, width=1, color='indigo', alpha=0.3)
        axes[0].add_patch(obs_repr)

    fig.tight_layout()
    plt.show()

def visualize_permittivity(Grid):
    number_of_plots = len(Grid.materials)
    list = np.arange(0, number_of_plots, 1)
    omega = np.arange(1e12, 1e18, 1e12)
    fig, axes = plt.subplots(1, number_of_plots, dpi=100)

    if number_of_plots == 1:
        axes.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes.semilogx(omega, Grid.materials[0].epsilon_real(omega), label=r'$\epsilon_{real}$')
        axes.semilogx(omega, Grid.materials[0].epsilon_imag(omega), label=r'$\epsilon_{imag}$')
        axes.legend(loc='best')
        axes.set_xlabel(r'$\omega$')

    else:
        for ax, mat in zip(axes, list):
            ax.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
            ax.semilogx(omega, Grid.materials[mat].epsilon_real(omega), label=r'$\epsilon_{real}$')
            ax.semilogx(omega, Grid.materials[mat].epsilon_imag(omega), label=r'$\epsilon_{imag}$')
            ax.legend(loc='best')
            ax.set_xlabel(r'$\omega$')
    plt.show()

def visualize_fft(grid):
    fft_objects = []
    for elem in grid.local_observers:
        if isinstance(elem, FFTObserver):
            fft_objects.append(elem)

    number_of_plots = len(fft_objects)
    fig, axes = plt.subplots(1, number_of_plots, dpi=100)
    list = np.arange(0, number_of_plots, 1)

    if number_of_plots == 1:
        frequency_x = 2*np.pi*np.linspace(0, 1/(2*grid.dt), (fft_objects[0].second_timestep - fft_objects[0].first_timestep)//2)
        axes.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
        axes.plot(frequency_x, 2/fft_objects[0].timestep_duration * np.abs(fft_objects[0].fft[0:fft_objects[0].timestep_duration//2]), marker='o', alpha=0.8, linestyle='-.')
        axes.set_title('Position {}'.format(fft_objects[0].position))
        axes.set_xlabel('Frequenz ' + r'$\omega$')
        axes.set_ylabel('Amplitude')

    else:
        for ax, fft_obs in zip(axes, list):
            frequency_x = 2*np.pi*np.linspace(0, 1/(2*grid.dt), fft_objects[fft_obs].timestep_duration//2)
            ax.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
            ax.plot(frequency_x, 2/fft_objects[fft_obs].timestep_duration * np.abs(fft_objects[fft_obs].fft[0:fft_objects[fft_obs].timestep_duration//2]), marker='o', alpha=0.8, linestyle='-.')
            ax.set_title('Position {}'.format(fft_objects[fft_obs].position))
            ax.set_xlabel('Frequenz '+r'$\omega$')
            ax.set_ylabel('Amplitude')
    plt.show()

class AnimateTillTimestep(ani.TimedAnimation):

    def __init__(self, grid_obj, final_timestep):
        self.grid = grid_obj
        fig_ani, self.axes_ani = plt.subplots(2, 1, dpi=100, figsize=(12, 10))
        self.final_timestep = final_timestep
        self.grid.timesteps = self.final_timestep

        self.axes_ani[0].xaxis.grid(linestyle='dotted')
        self.axes_ani[0].set_xlabel('Cell')
        self.axes_ani[0].set_ylabel(r'$E_z$', fontsize=12, rotation=0)
        self.axes_ani[1].xaxis.grid(linestyle='dotted')
        self.axes_ani[1].set_xlabel('Cell')
        self.axes_ani[1].set_ylabel(r'$H_y$', fontsize=12, rotation=0)

        self.line_0 = Line2D([], [])
        self.line_1 = Line2D([], [])

        self.axes_ani[0].add_line(self.line_0)
        self.axes_ani[1].add_line(self.line_1)

        self.axes_ani[0].set_xlim([0, self.grid.nx - 1])
        self.axes_ani[1].set_xlim([0, self.grid.nx - 1])
        self.axes_ani[0].set_ylim([-1.5, 1.5])
        self.axes_ani[1].set_ylim([-1.5/np.sqrt(mu0/eps0), 1.5/np.sqrt(mu0/eps0)])

        for mat in self.grid.materials:
            media_repr_0 = Rectangle(xy=(mat.position[0] - 0.5, -1.4), height=2.8, width=(mat.position[-1] - mat.position[0] + 1),
                                   color='grey', fill=True, alpha=mat.eps * 0.12)
            self.axes_ani[0].add_patch(media_repr_0)
            if mat.model == 'Lorentz':
                s = r'$\epsilon(\omega)$' + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity)
            else:
                s = r'$\epsilon_r$={0:.2f}'.format(mat.eps) + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity)

            self.axes_ani[0].annotate(s, xy=(media_repr_0.get_x() + 0.1, media_repr_0.get_y() + 0.2), color='black')
            media_repr_1 = Rectangle(xy=(mat.position[0] - 0.5 , -1.4), height=2.8, width=(mat.position[-1] - mat.position[0] + 1),
                                     color='grey', fill=True, alpha=mat.eps * 0.12)
            self.axes_ani[1].add_patch(media_repr_1)

        for src in self.grid.sources:
            src_repr = Rectangle(xy=(src.position - 0.5, -src.ampl), height=2 * src.ampl, width=1, color='red',
                                 alpha=0.3)
            self.axes_ani[0].add_patch(src_repr)

        for obs in self.grid.local_observers:
            if isinstance(obs, QuasiHarmonicObserver):
                obs_repr = Rectangle(xy=(obs.position - 0.5, -1.4), height=2.8, width=1, color='green', alpha=0.3)
            elif isinstance(obs, FFTObserver):
                obs_repr = Rectangle(xy=(obs.position - 0.5, -1.4), height=2.8, width=1, color='indigo', alpha=0.3)
            self.axes_ani[0].add_patch(obs_repr)

        ani.TimedAnimation.__init__(self, fig_ani, blit=True, interval=5, repeat=False)

    def _draw_frame(self, framedata):
        self.grid.timesteps_passed = framedata
        self.grid.update()
        self.axes_ani[0].set_title('time passed: {0:.3e}s,  timesteps passed: {timesteps}'.format(self.grid.time_passed,
                                                                                                  timesteps=self.grid.timesteps_passed))
        self.line_0.set_data(np.arange(self.grid.nx), self.grid.Ez)
        self.line_1.set_data(np.arange(self.grid.nx), self.grid.Hy)

        self._drawn_artists = [self.line_0, self.line_1]

    def new_frame_seq(self):
        return iter(range(self.final_timestep + 1))

    def _init_draw(self):
        lines = [self.line_0, self.line_1]
        for l in lines:
            l.set_data([], [])

    def create_animation(self):
        plt.show()