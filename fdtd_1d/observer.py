from fdtd_1d.utilities import get_amplitude_and_phase

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
        -> get Amplitude + Phase with amplitude > 0 and based on A cos (omega * t + phi)'''

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
