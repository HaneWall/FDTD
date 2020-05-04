import numpy as np
from .constants import c0, mu0, eps0
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
        return self.fwhm / (self.grid.dx * 2.355)       # approximating 2.255 = 2*sqrt(2*ln(2))

    def __init__(self, name, amplitude, peak_timestep, fwhm, tfsf=False):
        super().__init__()
        self.tfsf = tfsf
        self.peak_timestep = peak_timestep
        self.source_name = name
        self.fwhm = fwhm
        self.ampl = amplitude

    def step_Ez(self):
        self.grid.Ez[self.position] += self.ampl * np.exp(-0.5 * (((self.peak_timestep - 0.5) -
                                                                   (self.grid.timesteps_passed + 0.5)) / self.sigma) ** 2)
    def step_Hy(self):
        if not self.tfsf:
            pass
        else:
            self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * np.exp(-0.5 * ((self.peak_timestep -
                                                      self.grid.timesteps_passed) / self.sigma) ** 2)

class SinusoidalImpulse(ParentSource):
    '''creates oscillating source'''
    # currently only supports E-Sources

    def __init__(self, name, wavelength, phase_shift, amplitude,tfsf=False):
        super().__init__()
        self.tfsf = tfsf
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
    def step_Ez(self):
        self.grid.Ez[self.position] += self.ampl * np.sin(self.omega * (self.grid.time_passed + self.grid.dt) + self.phase)

    def step_Hy(self):
        if not self.tfsf:
            pass
        else:
            self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * np.sin(self.omega * self.grid.time_passed + self.phase)

class EnvelopeSinus(ParentSource):
    '''creates an enveloped oscillation (ampl * gaussian_impulse * sin(wt + phase))'''
    # currently only supports E-Sources

    def __init__(self, name, wavelength, phase_shift, amplitude, fwhm, peak_timestep, tfsf=False):
        super().__init__()
        self.tfsf = tfsf
        self.source_name = name
        self.lamb = wavelength
        self.phase = phase_shift
        self.ampl = amplitude
        self.fwhm = fwhm
        self.peak_timestep = peak_timestep

    @property
    def sigma(self):
        return self.fwhm / (self.grid.dx * 2.355)  # approximating 2.355 = 2*sqrt(2*ln(2))

    @property
    def omega(self):
        return 2 * np.pi * c0 / self.lamb

    @property
    def period(self):
        return self.lamb / c0


    def step_Ez(self):
        self.grid.Ez[self.position] += self.ampl * np.exp(-0.5 * (((self.peak_timestep - 0.5) -
                                                                   (self.grid.timesteps_passed + 0.5))/self.sigma) ** 2) * np.sin(self.omega * (self.grid.time_passed + self.grid.dt) + self.phase)

    def step_Hy(self):
        if not self.tfsf:
            pass
        else:
            self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * np.exp(-0.5 * ((self.peak_timestep -
                                                  self.grid.timesteps_passed)/self.sigma) ** 2) * np.sin(self.omega * self.grid.time_passed + self.phase)

class ActivatedSinus(ParentSource):
    '''sin squared as activation function'''

    def __init__(self, name, wavelength, carrier_wavelength, amplitude, phase_shift, tfsf=False):
        super().__init__()
        self.tfsf = tfsf
        self.name = name
        self.lamb = wavelength
        self.carrier_lamb = carrier_wavelength
        self.phase = phase_shift
        self.ampl = amplitude

    @property
    def carrier_omega(self):
        return 2 * np.pi * c0 / self.carrier_lamb

    @property
    def omega(self):
        return 2 * np.pi * c0 / self.lamb

    @property
    def period(self):
        return self.lamb / c0

    def step_Ez(self):
        if self.carrier_omega * self.grid.time_passed < np.pi / 2:
            self.grid.Ez[self.position] += self.ampl * (np.sin(self.carrier_omega * (self.grid.time_passed + self.grid.dt)))**2 * np.sin(self.omega * (self.grid.time_passed + self.grid.dt) + self.phase)

        else:
            self.grid.Ez[self.position] += self.ampl * np.sin(self.omega * (self.grid.time_passed + self.grid.dt) + self.phase)

    def step_Hy(self):
        if not self.tfsf:
            pass

        else:
            if self.carrier_omega * self.grid.time_passed < np.pi / 2:
                self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * (np.sin(self.carrier_omega * self.grid.time_passed))**2 * np.sin(self.omega * self.grid.time_passed + self.phase)

            else:
                self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * np.sin(self.omega * self.grid.time_passed + self.phase)