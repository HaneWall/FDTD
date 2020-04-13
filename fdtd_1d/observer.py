import numpy as np
from .grid import Grid


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

# differentiate between quasi harmonic situations (2 points required - minimal memory usage) e.g. ramped up signal / ONE frequency
# and fft (array for one wavelength) e.g. wave packages

class QuasiHarmonicObserver(ParentObserver):
    '''ramping up SinusoidalImpulse via ActivatedSinus -> waiting for steadystate
        -> save two points at the observers positions in time domain with T/4 away from each other -> reconstruct sinusoidal in time domain
        -> get Amplitude + Phase'''

# Note that phase and amplitude information is based on A*cos(wt + phi)

    def __init__(self, name, first_timestep):
        super().__init__()
        self.observer_name = name
        self.first_timestep = first_timestep
        self.observedE = []
        self.signed_phase = None
        self.signed_amplitude = None
        self.phase = None
        self.amplitude = None

    @property
    def second_timestep(self):
        return self.first_timestep + int(self.grid.sources[0].period / (4 * self.grid.dt))

    def save_E(self):
        if self.grid.timesteps_passed == self.first_timestep:
            self.observedE.append(self.grid.E[self.position])

        elif self.grid.timesteps_passed == self.second_timestep:
            self.observedE.append(self.grid.E[self.position])


    def _set_signed_phase(self):
        self.signed_phase = np.arctan2(self.observedE[1] * np.cos(self.grid.sources[0].omega * self.grid.dt * self.first_timestep) - self.observedE[0]*np.cos(self.grid.sources[0].omega * self.grid.dt * self.second_timestep),
                      self.observedE[1] * np.sin(self.grid.sources[0].omega * self.grid.dt * self.first_timestep) - self.observedE[0]*np.sin(self.grid.sources[0].omega * self.grid.dt * self.second_timestep))

    def _set_signed_amplitude(self):
        self.signed_amplitude = self.observedE[0] / (np.cos(self.grid.sources[0].omega * self.grid.dt * self.first_timestep) * np.cos(self.signed_phase) - np.sin(self.grid.sources[0].omega * self.grid.dt * self.first_timestep)*np.sin(self.signed_phase))

    # in order to get more beautiful phase and amplitude data (unsigned):

    def set_amplitude_phase(self):
        self._set_signed_phase()
        self._set_signed_amplitude()

        if self.signed_amplitude > 0:
            self.amplitude = self.signed_amplitude
            self.phase = self.signed_phase

        elif self.signed_amplitude < 0 and self.signed_phase >= 0:
            self.amplitude = -self.signed_amplitude
            self.phase = self.signed_phase - np.pi

        elif self.signed_amplitude < 0 and self.signed_phase < 0:
            self.amplitude = -self.signed_amplitude
            self.phase = self.signed_phase + np.pi