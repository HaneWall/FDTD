import numpy as np
from werkzeug.utils import cached_property
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

    @cached_property
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
        self.grid.Ez[self.position] += self.ampl * np.exp(-0.5 * (((self.peak_timestep) -
                                                                   (self.grid.timesteps_passed)) / self.sigma) ** 2)
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


    @cached_property
    def omega(self):
        return 2*np.pi*c0 / self.lamb

    @cached_property
    def period(self):
        return self.lamb / c0


    # as hard source - note that reflection information is forfeited in order to create perfect shape
    # for soft source change '=' to '+=' and vice versa
    def step_Ez(self):
        self.grid.Ez[self.position] += self.ampl * np.sin(self.omega * self.grid.time_passed + self.phase)

    def step_Hy(self):
        if not self.tfsf:
            pass
        else:
            self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * np.sin(self.omega * self.grid.time_passed + self.phase)


class GaussianImpulseWithFrequency(ParentSource):
    def __init__(self, name, Intensity, peak_timestep, pulse_duration, wavelength, tfsf=False, phase=0):
        super().__init__()
        self.tfsf = tfsf
        self.peak_timestep = peak_timestep
        self.source_name = name
        self.lamb = wavelength
        self.phase = phase
        self.pulse_duration = pulse_duration
        self.intensity = Intensity

    @cached_property
    def omega(self):
        return 2 * np.pi * c0 / self.lamb

    @cached_property
    def ampl(self):
        return np.sqrt(2*np.sqrt(mu0/eps0)*self.intensity)

    def step_Ez(self):
        self.grid.Ez[self.position] += self.ampl * np.exp(-2*np.log(2)*(((self.peak_timestep - self.grid.timesteps_passed)*self.grid.dt)/self.pulse_duration)**2) \
                                       * np.sin(self.omega * self.grid.time_passed + self.phase)

    def step_Hy(self):
        if not self.tfsf:
            pass
        else:
            self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * np.exp(-2*np.log(2)*(((self.peak_timestep - self.grid.timesteps_passed)*self.grid.dt)/self.pulse_duration)**2) \
                                       * np.sin(self.omega * self.grid.time_passed + self.phase)

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

    @cached_property
    def sigma(self):
        return self.fwhm / (self.grid.dx * 2*np.sqrt(2*np.log(2)))  # approximating 2.355 = 2*sqrt(2*ln(2))

    @cached_property
    def omega(self):
        return 2 * np.pi * c0 / self.lamb

    @cached_property
    def period(self):
        return self.lamb / c0


    def step_Ez(self):
        self.grid.Ez[self.position] += self.ampl * np.exp(-0.5 * (((self.peak_timestep - self.grid.timesteps_passed)/self.sigma) ** 2)) * np.sin(self.omega * self.grid.time_passed + self.phase)

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

    @cached_property
    def carrier_omega(self):
        return 2 * np.pi * c0 / self.carrier_lamb

    @cached_property
    def omega(self):
        return 2 * np.pi * c0 / self.lamb

    @cached_property
    def period(self):
        return self.lamb / c0

    def step_Ez(self):
        if self.carrier_omega * self.grid.time_passed < np.pi / 2:
            self.grid.Ez[self.position] += self.ampl * (np.sin(self.carrier_omega * self.grid.time_passed))**2 * np.sin(self.omega * self.grid.time_passed  + self.phase)

        else:
            self.grid.Ez[self.position] += self.ampl * np.sin(self.omega * self.grid.time_passed  + self.phase)

    def step_Hy(self):
        if not self.tfsf:
            pass

        else:
            if self.carrier_omega * self.grid.time_passed < np.pi / 2:
                self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * (np.sin(self.carrier_omega * self.grid.time_passed))**2 * np.sin(self.omega * self.grid.time_passed + self.phase)

            else:
                self.grid.Hy[self.position - 1] += -np.sqrt(eps0/mu0) * self.ampl * np.sin(self.omega * self.grid.time_passed + self.phase)