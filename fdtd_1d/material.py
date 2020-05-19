from .grid import Grid
from .constants import eps0



def _create_array_from_slice(slice):
    arr = []
    for cell in range(slice.start, slice.stop):
        arr.append(cell)
    return arr


class Vacuum:

    def __init__(self):
        self.eps = 1                    # permittivity
        self.mu = 1                     # permeability
        self.grid = None
        self.position = None
        self.material_name = 'vacuum'
        self.conductivity = 0           # sigma
        self.model = None

    def _place_into_grid(self, grid, index):
        if isinstance(index, int):
            self.grid = grid
            self.grid.materials.append(self)
            self.position = index
            self.grid.eps[index] = self.eps
            self.grid.mu[index] = self.mu
            self.grid.conductivity[index] = self.conductivity
            self.grid.all_E_mats.add(index)

        elif isinstance(index, slice):          # note that [1:5] means cell INDICES 1, 2, 3, 4 are getting manipulated
            self.grid = grid
            self.grid.materials.append(self)
            arr = _create_array_from_slice(index)
            self.position = arr
            for cell in arr:
                self.grid.eps[cell] = self.eps
                self.grid.mu[cell] = self.mu
                self.grid.conductivity[cell] = self.conductivity
                self.grid.all_E_mats.add(cell)

        else:
            raise KeyError('Currently not supporting these kind of keys! Use slice or simple index.')


# defining subclasses - convenient materials for investigations


class NonDispersiveMedia(Vacuum):

    def __init__(self, name, permittivity, permeability, conductivity):
        super().__init__()
        self.material_name = name
        self.eps = permittivity
        self.mu = permeability
        self.conductivity = conductivity

    def epsilon_real(self):
        return self.eps

    def epsilon_imag(self):
        return 0

    def step_J_p(self, index):
        pass

class LorentzMedium(Vacuum):
    ''' electron cloud behaves like harmonic oscillator '''

    def __init__(self, name, permeability, conductivity, eps_inf, gamma, w0, chi_1):
        super().__init__()
        self.material_name = name
        self.model = 'Lorentz'
        self.eps = eps_inf
        self.mu = permeability
        self.conductivity = conductivity
        self.gamma = gamma
        self.w_0 = w0
        self.chi_1 = chi_1

    @property
    def a(self):
        return (self.grid.dt * self.w_0**2) / (1 + self.gamma/2 * self.grid.dt)

    @property
    def b(self):
        return (1 - self.gamma/2 * self.grid.dt) / (1 + self.gamma/2 * self.grid.dt)

    def epsilon_real(self, omega):
        return self.eps + (self.chi_1 * self.w_0**2) / ((self.w_0**2 - omega**2)**2 + self.gamma**2 * omega**2)

    def epsilon_imag(self, omega):
        return (-self.chi_1 * self.w_0**2 * self.gamma * omega) / ((self.w_0**2 - omega**2)**2 + self.gamma**2 * omega**2)

    def step_J_p(self, index):
        self.grid.J_p[index] = self.b * self.grid.J_p[index] + self.a * (eps0 * self.chi_1 * self.grid.Ez[index] - self.grid.P[index])

class CustomMedia(Vacuum):
    # UNDER CONSTRUCTION!!!

    def __init__(self, name, permittivity, permeability, conductivity, dispersion_model):
        super().__init__()
        self.material_name = name
        self.eps = permittivity
        self.mu = permeability
        self.conductivity = conductivity
        self.model = dispersion_model

# ... for example class Silicon(Vacuum): ...?!
