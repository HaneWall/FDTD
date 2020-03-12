import numpy as np
from .constants import c0
from .grid import Grid


class ParentSource:

    def __init__(self):
        self.position = None
        self.grid = None
        self.source_name = None

    def _place_into_grid(self, grid, index):
        if isinstance(index, int):
            self.grid = grid
            self.grid.sources.append(self)
            self.position = index
        elif isinstance(index, slice):
            raise KeyError('Not supporting slicing for sources!')

# defining subclasses - some sources

class GaussianImpulse(ParentSource):
    '''creates a bell-shaped curve'''
    # currently only supports E-Sources

    @property
    def sigma(self):
        return self.fwhm / (self.grid.dz * 2.355)       # approximating 2.255 = 2*sqrt(2*ln(2))

    def __init__(self, name, amplitude, peak_timestep, fwhm):
        super().__init__()
        self.peak_timestep = peak_timestep
        self.source_name = name
        self.fwhm = fwhm
        self.ampl = amplitude

    # as soft source
    def step_E(self):
        self.grid.E[self.position] += self.ampl * np.exp(-0.5 * ((self.peak_timestep -
                                                      self.grid.timesteps_passed) / self.sigma) ** 2)


class SinusoidalImpulse(ParentSource):
    '''creates oscillating source'''
    # currently only supports E-Sources

    def __init__(self, name, wavelength, phase_shift, amplitude):
        super().__init__()
        self.source_name = name
        self.lamb = wavelength
        self.phase = phase_shift
        self.ampl = amplitude


    @property
    def omega(self):
        return 2*np.pi*c0 / self.lamb

    @property
    def period(self):
        return self.lamb / c0


    # as hard source - note that reflection information is forfeited in order to create perfect shape
    # for soft source change '=' to '+=' and vice versa
    def step_E(self):
        self.grid.E[self.position] += self.ampl * np.sin(self.omega * self.grid.time_passed + self.phase)

class EnvelopeSinus(ParentSource):
    '''creates an enveloped oscillation (ampl * gaussian_impulse * sin(wt + phase))'''
    # currently only supports E-Sources

    def __init__(self, name, wavelength, phase_shift, amplitude, fwhm, peak_timestep):
        super().__init__()
        self.source_name = name
        self.lamb = wavelength
        self.phase = phase_shift
        self.ampl = amplitude
        self.fwhm = fwhm
        self.peak_timestep = peak_timestep

    @property
    def sigma(self):
        return self.fwhm / (self.grid.dz * 2.355)  # approximating 2.355 = 2*sqrt(2*ln(2))

    @property
    def omega(self):
        return 2 * np.pi * c0 / self.lamb

    @property
    def period(self):
        return self.lamb / c0


    # as hard source - note that reflection information is forfeited in order to create perfect shape
    # for soft source change '=' to '+=' and vice versa
    def step_E(self):
        self.grid.E[self.position] += self.ampl * np.exp(-0.5 * ((self.peak_timestep -
                                                                  self.grid.timesteps_passed)/self.sigma) ** 2) * np.sin(self.omega * self.grid.time_passed + self.phase)