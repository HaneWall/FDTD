from fdtd_1d.utilities import get_amplitude_and_phase
import numpy as np
import csv
import os.path
from werkzeug.utils import cached_property

class ParentObserver:

    def __init__(self):
        self.position = None
        self.grid = None
        self.observer_name = None

    def _place_into_grid(self, grid, index):
        if isinstance(index, int):
            self.grid = grid
            self.grid.local_observers.append(self)
            self.position = index
        elif isinstance(index, slice):
            raise KeyError('Not supporting slicing for observer!')

    def save_Ez(self):
        pass

    def save_P(self):
        pass

# differentiate between quasi harmonic situations (2 points required - minimal memory usage) e.g. ramped up signal / ONE frequency
# and fft (array for one wavelength) e.g. wave packages

class QuasiHarmonicObserver(ParentObserver):
    '''
    ramping up SinusoidalImpulse via ActivatedSinus -> waiting for steadystate
        -> save two points at the observers positions in time domain with T/4 away from each other -> reconstruct sinusoidal in time domain
        -> get Amplitude + Phase with amplitude > 0 and based on A cos (omega * t + phi)
    '''

# Note that phase and amplitude information is based on A*cos(wt + phi)

    def __init__(self, name, first_timestep):
        super().__init__()
        self.observer_name = name
        self.first_timestep = first_timestep
        self.observedE = []

    @property
    def second_timestep(self):
        return self.first_timestep + int(self.grid.sources[0].period / (4 * self.grid.dt))

    @property
    def phase(self):
        return get_amplitude_and_phase(grid=self.grid, first_timestep=self.first_timestep,
                                       second_timestep=self.second_timestep, data=self.observedE)[1]

    @property
    def amplitude(self):
        return get_amplitude_and_phase(grid=self.grid, first_timestep=self.first_timestep,
                                       second_timestep=self.second_timestep, data=self.observedE)[0]

    def save_Ez(self):
        if self.grid.timesteps_passed == self.first_timestep:
            self.observedE.append(self.grid.Ez[self.position])

        elif self.grid.timesteps_passed == self.second_timestep:
            self.observedE.append(self.grid.Ez[self.position])

class E_FFTObserver(ParentObserver):
    '''
    stores an array of E at specific position from timestep1 to timestep2
    '''

    def __init__(self, name, first_timestep, second_timestep):
        super().__init__()
        self.observer_name = name
        self.first_timestep = first_timestep
        self.second_timestep = second_timestep
        self.observed_E = []
        self.Ez_fft = []

    @cached_property
    def timestep_duration(self):
        return self.second_timestep - self.first_timestep

    @cached_property
    def fft(self):
        self.Ez_fft = np.fft.fft(self.observed_E)
        return self.Ez_fft

    # saving Ez in order to fft
    def save_Ez(self):
        if self.grid.timesteps_passed in range(self.first_timestep, self.second_timestep + 1):
            self.observed_E.append(self.grid.Ez[self.position])

    # store Ez (and Ez_fft) in order to analyze data wo computing simulation again
    def store_Ez_data(self, filename):
        if self.grid.benchmark_type is None:
            filepath_0 = os.path.join(os.path.dirname(__file__), 'saved_data')

        else:
            filepath_0 = os.path.join(os.path.dirname(__file__), 'saved_data/'+self.grid.benchmark_type)

        filepath_1 = os.path.join(filepath_0, filename)
        with open(filepath_1, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['dt', self.grid.dt, 'timestep_duration', self.timestep_duration])
            writer.writerow(self.observed_E)
            #writer.writerow(self.Ez_fft)



class P_FFTObserver(ParentObserver):
    '''
    stores an array of P at specific position from timestep1 to timestep2
    '''

    def __init__(self, name, first_timestep, second_timestep):
        super().__init__()
        self.observer_name = name
        self.first_timestep = first_timestep
        self.second_timestep = second_timestep
        self.observed_P = []
        self.P_fft = []

    @cached_property
    def timestep_duration(self):
        return self.second_timestep - self.first_timestep

    @cached_property
    def fft(self):
        self.P_fft = np.fft.fft(self.observed_P)
        return self.P_fft

    def save_P(self):
        if self.grid.timesteps_passed in range(self.first_timestep, self.second_timestep + 1):
            self.observed_P.append(self.grid.P[self.position])

    def store_P_data(self, filename):
        if self.grid.benchmark_type is None:
            filepath_0 = os.path.join(os.path.dirname(__file__), 'saved_data')

        else:
            filepath_0 = os.path.join(os.path.dirname(__file__), 'saved_data/' + self.grid.benchmark_type)

        filepath_1 = os.path.join(filepath_0, filename)
        with open(filepath_1, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['dt', self.grid.dt, 'timestep_duration', self.timestep_duration])
            writer.writerow(self.observed_P)
            #writer.writerow(self.P_fft)