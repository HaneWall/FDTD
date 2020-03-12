import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as ani
from matplotlib.patches import Rectangle
import numpy as np
from .constants import c0

def visualize(Grid):
    fig, axes = plt.subplots(2, 1, dpi=100, figsize=(12,8))
    axes[0].plot(Grid.E)
    axes[0].xaxis.grid(linestyle='dotted')
    axes[0].set_xlabel('Cell')
    axes[0].set_ylabel('E', fontsize=12, rotation=0)
    axes[0].set_title('time passed: {0:.3e}s,  timesteps passed: {timesteps}'.format(Grid.time_passed,
                                                                                timesteps=Grid.timesteps_passed))
    axes[0].set_ylim([-1.5, 1.5])
    axes[1].plot(Grid.B)
    axes[1].xaxis.grid(linestyle='dotted')
    axes[1].set_xlabel('Cell')
    axes[1].set_ylabel('B', fontsize=12, rotation=0)
    for mat in Grid.materials:
        media_repr_0 = Rectangle(xy=(mat.position[0], -1.4), height=2.8, width=(mat.position[-1] - mat.position[0]),
                               color='grey', fill=True, alpha=0.3)
        axes[0].add_patch(media_repr_0)
        axes[0].annotate(
            s=r'$\epsilon_r$={0:.2f}'.format(mat.eps) + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity),
            xy=(media_repr_0.get_x() + 0.1, media_repr_0.get_y() + 0.2), color='black')
        media_repr_1 = Rectangle(xy=(mat.position[0], -1.4), height=2.8, width=(mat.position[-1] - mat.position[0]),
                                 color='grey', fill=True, alpha=0.3)
        axes[1].add_patch(media_repr_1)
        axes[1].annotate(
            s=r'$\epsilon_r$={0:.2f}'.format(mat.eps) + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity),
            xy=(media_repr_1.get_x() + 0.1, media_repr_1.get_y() + 0.2), color='black')
    for src in Grid.sources:
        src_repr = Rectangle(xy=(src.position - 0.5, -src.ampl), height=2 * src.ampl, width=1, color='red', alpha=0.3)
        axes[0].add_patch(src_repr)
    fig.tight_layout()
    plt.show()


class AnimateTillTimestep:

    def __init__(self, grid_obj, final_timestep):
        self.grid = grid_obj
        self.fig_ani, self.axes_ani = plt.subplots(2, 1, dpi=100, figsize=(12,8))
        self.final_timestep = final_timestep
        self.animation = None
        self.grid.timesteps = self.final_timestep
        self.axes_ani[0].xaxis.grid(linestyle='dotted')
        self.axes_ani[0].set_xlabel('Cell')
        self.axes_ani[0].set_ylabel('E', fontsize=12, rotation=0)
        self.axes_ani[1].xaxis.grid(linestyle='dotted')
        self.axes_ani[1].set_xlabel('Cell')
        self.axes_ani[1].set_ylabel('B', fontsize=12, rotation=0)
        self.line_0 = Line2D(np.arange(self.grid.nz), self.grid.E)
        self.line_1 = Line2D(np.arange(self.grid.nz), self.grid.B)
        self.axes_ani[0].add_line(self.line_0)
        self.axes_ani[1].add_line(self.line_1)
        self.axes_ani[0].set_xlim([0, self.grid.nz -1])
        self.axes_ani[1].set_xlim([0, self.grid.nz - 1])
        self.axes_ani[0].set_ylim([-1.5, 1.5])
        self.axes_ani[1].set_ylim([-1.5 / c0, 1.5 / c0])

        for mat in self.grid.materials:
            media_repr_0 = Rectangle(xy=(mat.position[0], -1.4), height=2.8, width=(mat.position[-1] - mat.position[0]),
                                   color='grey', fill=True, alpha=0.3)
            self.axes_ani[0].add_patch(media_repr_0)
            self.axes_ani[0].annotate(
                s=r'$\epsilon_r$={0:.2f}'.format(mat.eps) + '\n' + r'$\sigma$={0:.2f}'.format(mat.conductivity),
                xy=(media_repr_0.get_x() + 0.1, media_repr_0.get_y() + 0.2), color='black')

            media_repr_1 = Rectangle(xy=(mat.position[0], -1.4), height=2.8, width=(mat.position[-1] - mat.position[0]),
                                     color='grey', fill=True, alpha=0.3)
            self.axes_ani[1].add_patch(media_repr_1)

        for src in self.grid.sources:
            src_repr = Rectangle(xy=(src.position - 0.5, -src.ampl), height=2 * src.ampl, width=1, color='red',
                                 alpha=0.3)
            self.axes_ani[0].add_patch(src_repr)

    def step(self, i):
        self.grid.timesteps_passed = i
        self.grid.update()
        self.axes_ani[0].set_title('time passed: {0:.3e}s,  timesteps passed: {timesteps}'.format(self.grid.time_passed,
                                                                                                  timesteps=self.grid.timesteps_passed))
        self.line_0.set_data(np.arange(self.grid.nz), self.grid.E)
        self.line_1.set_data(np.arange(self.grid.nz), self.grid.B)

        return self.line_0, self.line_1

    def create_animation(self):
        # fps = 1000 / interval
        self.animation = ani.FuncAnimation(self.fig_ani, func=self.step, frames=self.final_timestep + 1, blit=False, interval=100, repeat=False)
        self.fig_ani.tight_layout()
        plt.show()