import numpy as np
import pandas as pd

from .type_management import ndarray, Number, float64, int64
from .constants import c0, eps0, mu0
from .visuals import visualize, AnimateTillTimestep, visualize_permittivity, visualize_fft
from .observer import QuasiHarmonicObserver, E_FFTObserver, P_FFTObserver
from werkzeug.utils import cached_property
import time
from tqdm import tqdm
#from progressbar import ProgressBar
from numba import njit

#pbar = ProgressBar()

class Grid:

    def __init__(self, nx, dx, courant = 1, benchmark=None, name=None):
        self.mu = np.ones(nx)               # permeability - free space / vacuum
        self.eps = np.ones(nx)              # permittivity - free space / vacuum
        self.conductivity = np.zeros(nx)
        self.nx = nx                        # nx: number of cells in x-direction
        self.dx = dx                        # dx: width of one cell in m
        self.Ez = np.zeros(nx)
        self.Hy = np.zeros(nx)
        self.J_p = np.zeros(nx)
        self.P = np.zeros(nx)
        self.courant = courant                  # 1 = magic time step ( Taflove - numerical error is minimal )
        self.dt = dx * courant / c0
        self.timesteps = None
        self.timesteps_passed = 0
        self.grid_name = name               # needed if benchmarks create multiple grids (distinguish them)
        self.benchmark_type = benchmark     # some objects have to know which kind of benchmark is processed to work properly


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


    def curl_Ez(self) -> ndarray:
        return self.Ez[1:] - self.Ez[:-1]


    def curl_Hy(self) -> ndarray:
        return self.Hy[1:] - self.Hy[:-1]


    def run_time(self, simulate_t, vis=True):
        # simulate_t in s
        self.timesteps = int(simulate_t / self.dt)

        for time_step in tqdm(range(1, self.timesteps + 1)):          # range output: [..)
            self.update()
            self.timesteps_passed += 1
        if vis:
            visualize(self)

    def run_timesteps(self, timesteps, vis=True):
        # discrete timesteps
        self.timesteps = timesteps

        for time_step in tqdm(range(1, self.timesteps + 1)):
            self.update()
            self.timesteps_passed += 1
        if vis:
            visualize(self)

    def animate_timesteps(self, timesteps):
        vid = AnimateTillTimestep(grid_obj=self, final_timestep=timesteps)
        vid.create_animation()

    def visualize_permittivity(self):
        visualize_permittivity(self)

    def visualize_fft_observed(self):
        visualize_fft(self)

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

    @cached_property
    def ca(self):                            # Taflove convention
        return (1 - (self.conductivity * self.dt) / (2 * eps0 * self.eps)) / (1 + (self.conductivity * self.dt) / (2 * eps0 * self.eps))

    @cached_property
    def cb(self):                            # Taflove convention
        return 1 / (1 + (self.conductivity * self.dt) / (2 * eps0 * self.eps))


    def step_Ez(self):
        self.Ez[1:self.nx] = self.ca[1:self.nx] * self.Ez[1:self.nx] + self.dt / (eps0 * self.eps[1:self.nx]) * self.cb[1:self.nx] * (self.curl_Hy()/self.dx - self.J_p[1:self.nx])


    def step_Hy(self):
        self.Hy[0:self.nx-1] = self.Hy[0:self.nx-1] + self.dt/(mu0 * self.mu[0:self.nx-1] * self.dx) * self.curl_Ez()

    # main algorithm

    def update(self):
        # note that steps are dependent on object
        # updating polarisation P

        # saving local points in order to extract phase and amplitude data
        for observer in self.local_observers:
            observer.save_Ez()
            observer.save_P()

        for mat in self.materials:
            #start_time = time.time()
            mat.step_P_tilde()
            #print("computed P tilde in --- %s seconds ---" % (time.time() - start_time))
            #start_time = time.time()
            mat.step_J_p()
            #print("computed J_p in --- %s seconds ---" % (time.time() - start_time))
            #start_time = time.time()
            mat.step_P()
            #print("computed P in --- %s seconds ---" % (time.time() - start_time))
            #start_time = time.time()
            mat.step_G()
            #print("computed G in --- %s seconds ---" % (time.time() - start_time))
            #start_time = time.time()
            mat.step_Q()
            #print("computed Q in --- %s seconds ---" % (time.time() - start_time))


        # updating Hy - boundaries
        for bound in self.boundaries:
            bound.save_Hy()

        # updating Hy - field, index [first cell, last cell - 1]
        #start_time = time.time()
        self.step_Hy()
        #print("computed H in --- %s seconds ---" % (time.time() - start_time))

        # updating Hy - Sources
        for source in self.sources:
            source.step_Hy()

        # updating Hy - boundaries
        for bound in self.boundaries:
            bound.step_Hy()

        # saving E - boundaries
        for bound in self.boundaries:
            bound.save_Ez()

        # updating E - field, index [1,last cell]
        #start_time = time.time()
        self.step_Ez()
        #print("computed E in --- %s seconds ---" % (time.time() - start_time))

        # updating E - Sources
        for source in self.sources:
            source.step_Ez()

        # updating E - boundaries
        for bound in self.boundaries:
            bound.step_Ez()







