from .constants import c0, eps0
import numpy as np
from .visuals import visualize, AnimateTillTimestep

def curl_E(field, cell):
    return field[cell + 1] - field[cell]


def curl_B(field, cell):
    return field[cell] - field[cell - 1]


class Grid:

    def __init__(self, nz, dz):
        self.mu = np.ones(nz)               # permeability - free space / vacuum
        self.eps = np.ones(nz)              # permittivity - free space / vacuum
        self.conductivity = np.zeros(nz)
        self.nz = nz                        # nz: number of cells in z-direction
        self.dz = dz                        # dz: width of one cell in m
        self.E = np.zeros(nz)
        self.B = np.zeros(nz)
        self.courant = 0.5                  # 1 = magic time step ( Taflove - numerical error is minimal ), ≤1! 0.5 for absorbing boundary conditions
        self.dt = dz * self.courant / c0
        self.timesteps = None
        self.timesteps_passed = 0
        self.sources = []                   # saving source.py-objects
        self.materials = []                 # saving material.py-objects
        self.boundaries = []                # saving boundary.py-objects


    @property
    def time_passed(self):
        return self.timesteps_passed * self.dt

    # overwriting input 'obj[] = ..' operator
    def __setitem__(self, key, placing_obj):
        # placing_obj is representing materials, sources or boundaries
        # note that each obj owns a different declaration of _place_into_grid
        placing_obj._place_into_grid(grid=self, index=key)


    def run_time(self, simulate_t):
        # simulate_t in s
        self.timesteps = int(simulate_t / self.dt)

        for time_step in range(1, self.timesteps + 1):          # range output: [..)
            self.update()
            self.timesteps_passed += 1

        visualize(self)

    def run_timesteps(self, timesteps):
        # discrete timesteps
        self.timesteps = timesteps

        for time_step in range(1, self.timesteps + 1):
            self.update()
            self.timesteps_passed += 1

        visualize(self)

    def animate_timesteps(self, timesteps):
        vid = AnimateTillTimestep(grid_obj=self, final_timestep=timesteps)
        vid.create_animation()

    def cc_minus(self, cell):                      # conductivity coefficient for non dispersive media 1 - ..
        coefficient_minus = (1 - (self.conductivity[cell] * self.dt) / (2 * eps0 * self.eps[cell]))
        return coefficient_minus

    def cc_plus(self, cell):
        coefficient_plus = (1 + (self.conductivity[cell] * self.dt) / (2 * eps0 * self.eps[cell]))
        return coefficient_plus

    def step_E(self, cell):
        self.E[cell] = self.cc_minus(cell) / self.cc_plus(cell) * self.E[cell] - c0 * self.courant / (self.eps[cell] * self.cc_plus(cell)) * curl_B(self.B, cell)

    def step_B(self, cell):
        self.B[cell] = self.B[cell] - self.courant / c0 * curl_E(self.E, cell)

    # main algorithm
    def update(self):
        # note that step_E() is dependent on object

        # saving E - boundaries
        for bound in self.boundaries:
            bound.save_E()

        # uüpdating E - field, index [1,last cell]
        for index in range(1, self.nz):
            self.step_E(index)

        # updating E - Sources
        for source in self.sources:
            source.step_E()

        # updating E - boundaries
        for bound in self.boundaries:
            bound.step_E()

        # saving B - boundaries
        for bound in self.boundaries:
            bound.save_B()

        # updating B - field, index [first cell, last cell - 1]
        for index in range(0, self.nz - 1):
            self.step_B(index)

        # updating B - Sources (under construction)

        # updating B - boundaries
        for bound in self.boundaries:
            bound.step_B()


