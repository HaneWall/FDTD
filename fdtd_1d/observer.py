from fdtd_1d.utilities import get_amplitude_and_phase, numerical_group_velocity
from fdtd_1d.constants import c0
import numpy as np
import csv
import os
from werkzeug.utils import cached_property

def _create_array_from_slice(slice):
    arr = []
    for cell in range(slice.start, slice.stop):
        arr.append(cell)
    return arr


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
            self.grid = grid
            self.grid.local_observers.append(self)
            arr = _create_array_from_slice(index)
            self.position = arr

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
        self.type = 'E'
        self.first_timestep = first_timestep
        self.second_timestep = second_timestep
        self.observed_E = np.zeros(self.second_timestep-self.first_timestep + 1)
        self.Ez_fft = None

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
            self.observed_E[self.grid.timesteps_passed-self.first_timestep] = (self.grid.Ez[self.position])

class P_FFTObserver(ParentObserver):
    '''
    stores an array of P at specific position from timestep1 to timestep2
    '''

    def __init__(self, name, first_timestep, second_timestep):
        super().__init__()
        self.observer_name = name
        self.type = 'P'
        self.first_timestep = first_timestep
        self.second_timestep = second_timestep
        self.observed_P = np.zeros(self.second_timestep - self.first_timestep + 1)
        self.P_fft = None

    @cached_property
    def timestep_duration(self):
        return self.second_timestep - self.first_timestep

    @cached_property
    def fft(self):
        self.P_fft = np.fft.fft(self.observed_P)
        return self.P_fft

    def save_P(self):
        if self.grid.timesteps_passed in range(self.first_timestep, self.second_timestep + 1):
            self.observed_P[self.grid.timesteps_passed-self.first_timestep] = self.grid.P[self.position]

class MovingFrame(ParentObserver):

    def __init__(self, x_to_snapshot, central_wavelength):
        super().__init__()
        self.x_store = x_to_snapshot
        self.central_wavelength = central_wavelength
        self.central_omega = 2*np.pi*c0/central_wavelength
        self.i = 0
        self.stored_data = None

    @cached_property
    def t_start(self):
        delay_distance_in_dx = self.position[0] + int((self.position[-1] - self.position[0])/2) - self.grid.sources[0].position
        peak_timestep = self.grid.sources[0].peak_timestep
        return (peak_timestep + int((delay_distance_in_dx * c0) / (self.corrected_group_velocity * self.grid.courant)))

    @cached_property
    def corrected_group_velocity(self):
        return numerical_group_velocity(central_wavelength=self.central_wavelength,
                                        group_velocity=self.grid.materials[0].group_velocity(omega=self.central_omega),
                                        dx=self.grid.dx, courant=self.grid.courant,
                                        n_real=self.grid.materials[0].n_real(self.central_omega))


    def _allocate_memory(self):
        self.stored_data = np.empty(shape=(len(self.x_store), len(self.position)))
        self.position_frame = np.array(self.position)
        self.timesteps_to_store = np.array([self.t_start + int((x * c0)/(self.corrected_group_velocity * self.grid.dx * self.grid.courant)) for x in self.x_store])

    def _new_position(self):
        movement = int((self.grid.timesteps_passed - self.t_start) * self.corrected_group_velocity * self.grid.courant / c0)
        self.position_frame = np.array(self.position) + movement

    def save_Ez(self):
        if self.stored_data is None:
            self._allocate_memory()

        if self.grid.timesteps_passed in self.timesteps_to_store:
            self._new_position()
            self.stored_data[self.i][:] = self.grid.Ez[self.position_frame[0]:(self.position_frame[-1] + 1)]
            self.i += 1
        else:
            pass

