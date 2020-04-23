from .constants import c0, eps0, mu0
import numpy as np
import pandas as pd
from .visuals import visualize, AnimateTillTimestep

def curl_Ez(field, cell):
    return field[cell + 1] - field[cell]


def curl_Hy(field, cell):
    return field[cell] - field[cell - 1]

class Grid:

    def __init__(self, nx, dx):
        self.mu = np.ones(nx)               # permeability - free space / vacuum
        self.eps = np.ones(nx)              # permittivity - free space / vacuum
        self.conductivity = np.zeros(nx)
        self.nx = nx                        # nx: number of cells in x-direction
        self.dx = dx                        # dx: width of one cell in m
        self.Ez = np.zeros(nx)
        self.Hy = np.zeros(nx)
        self.courant = 1                   # 1 = magic time step ( Taflove - numerical error is minimal )
        self.dt = dx * self.courant / c0
        self.timesteps = None
        self.timesteps_passed = 0
        self.sources = []                   # saving source.py-objects
        self.materials = []                 # saving material.py-objects
        self.boundaries = []                # saving boundary.py-objects
        self.local_observers = []           # saving time-data of the field

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

        #visualize(self)

    def animate_timesteps(self, timesteps):
        vid = AnimateTillTimestep(grid_obj=self, final_timestep=timesteps)
        vid.create_animation()

    def get_observed_signals(self):
        dict = {'name': [], 'position': [], 'first timestep': [], 'second timestep': [], 'amplitude': [], 'phase': []}
        for observer in self.local_observers:
            dict['name'].append(observer.observer_name)
            dict['position'].append(observer.position)
            dict['first timestep'].append(observer.first_timestep)
            dict['second timestep'].append(observer.second_timestep)
            dict['amplitude'].append(observer.amplitude)
            dict['phase'].append(observer.phase)
        tab_sig = pd.DataFrame.from_dict(dict)
        print(tab_sig.to_string())


    def ca(self, cell):                            # Taflove convention
        return (1 - (self.conductivity[cell] * self.dt) / (2 * eps0 * self.eps[cell])) / (1 + (self.conductivity[cell] * self.dt) / (2 * eps0 * self.eps[cell]))

    def cb(self, cell):                            # Taflove convention
        return 1 / (1 + (self.conductivity[cell] * self.dt) / (2 * eps0 * self.eps[cell]))

    def step_Ez(self, cell):
        self.Ez[cell] = self.ca(cell) * self.Ez[cell] + self.dt / (eps0 * self.eps[cell] * self.dx) * self.cb(cell) * curl_Hy(self.Hy, cell)

    def step_Hy(self, cell):
        self.Hy[cell] = self.Hy[cell] + self.dt/(mu0 * self.mu[cell] * self.dx) * curl_Ez(self.Ez, cell)

    # main algorithm
    def update(self):
        # note that steps are dependent on object
        # saving Hy - boundaries
        for bound in self.boundaries:
            bound.save_Hy()

        # updating Hy - field, index [first cell, last cell - 1]
        for index in range(0, self.nx - 1):
            self.step_Hy(index)

        # updating Hy - Sources
        for source in self.sources:
            source.step_Hy()

        # updating Hy - boundaries
        for bound in self.boundaries:
            bound.step_Hy()

        # saving Ez - boundaries
        for bound in self.boundaries:
            bound.save_Ez()

        # updating E - field, index [1,last cell]
        for index in range(1, self.nx):
            self.step_Ez(index)

        # updating E - Sources
        for source in self.sources:
            source.step_Ez()

        # updating E - boundaries
        for bound in self.boundaries:
            bound.step_Ez()

        # saving local points in order to extract phase and amplitude data
        for observer in self.local_observers:
            observer.save_Ez()




