from .grid import Grid
from werkzeug.utils import cached_property
import numpy as np
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

        elif isinstance(index, slice):          # note that [1:5] means cell INDICES 1, 2, 3, 4 are getting manipulated
            self.grid = grid
            self.grid.materials.append(self)
            arr = _create_array_from_slice(index)
            self.position = arr
            for cell in arr:
                self.grid.eps[cell] = self.eps
                self.grid.mu[cell] = self.mu
                self.grid.conductivity[cell] = self.conductivity

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

    def epsilon_real(self, omega):
        return self.eps

    def epsilon_imag(self, omega):
        return 0

    def epsilon_complex(self, omega):
        return self.eps

    def step_P(self):
        pass

    def step_J_p(self):
        pass

class LorentzMedium(Vacuum):
    ''' electron cloud behaves like one or multiple harmonic oscillator '''
    # Todo: implement chi in such a way, that you can decribe chi_n -> variable E_vec length (inf order)

    def __init__(self, name, permeability, conductivity, eps_inf, gamma, w0, chi_1, chi_2, chi_3):
        super().__init__()
        self.material_name = name
        self.model = 'Lorentz'
        self.eps = eps_inf
        self.mu = permeability
        self.conductivity = conductivity
        self.gamma = np.array(gamma)
        self.w_0 = np.array(w0)
        self.chi_1 = chi_1
        self.chi_2 = chi_2
        self.chi_3 = chi_3
        self.J_p_k = None
        self.P_k = None
        self.start_of_media = None

    def _allocate_J_p_k(self):
        self.J_p_k = np.zeros(shape=(len(self.position), len(self.chi_1)))

    def _allocate_P_k(self):
        self.P_k = np.zeros(shape=(len(self.position), len(self.chi_1)))



    @cached_property
    def a(self):
        a = np.zeros(len(self.w_0))
        for w_k, gamma_k, ind in zip(self.w_0, self.gamma, range(len(self.w_0))):
            a[ind] = (self.grid.dt * w_k**2) / (1 + gamma_k/2 * self.grid.dt)
        return a

    @cached_property
    def b(self):
        b = np.zeros(len(self.w_0))
        for gamma_k, ind in zip(self.gamma, range(len(self.w_0))):
            b[ind] = (1 - gamma_k/2 * self.grid.dt) / (1 + gamma_k/2 * self.grid.dt)
        return b

    @cached_property
    def chi_matrix(self):
        chi_m = np.array([self.chi_1, self.chi_2, self.chi_3])
        return chi_m

    def epsilon_real(self, omega):
        eps_real = np.real(self.epsilon_complex(omega))
        return eps_real

    def epsilon_imag(self, omega):
        eps_imag = np.imag(self.epsilon_complex(omega))
        return eps_imag

    def epsilon_complex(self, omega):
        eps_complex = self.eps
        for w_k, gamma_k, chi_1_k in zip(self.w_0, self.gamma, self.chi_1):
            eps_complex += chi_1_k * (w_k**2) / (w_k**2 - omega**2 - 1j * gamma_k * omega)
        return eps_complex

    def step_P(self):
        if self.J_p_k is None:
            self._allocate_J_p_k()
            self._allocate_P_k()

        self.P_k[0:len(self.position)] = self.P_k[0:len(self.position)] + self.grid.dt * self.J_p_k[0:len(self.position)]
        self.grid.P[self.position[0]:(self.position[-1] + 1)] = np.sum(self.P_k, axis=1)

    def step_J_p(self):
        if self.J_p_k is None:
            self._allocate_J_p_k()
            self._allocate_P_k()

        # define relative_index
        E_matrix = np.zeros(shape=(len(self.position), 3))
        for ind, pos in zip(range(len(self.position)), self.position):
            E_matrix[ind] = [self.grid.Ez[pos]**i for i in range(1, 4)]

        self.J_p_k[0:len(self.position)] = self.b * self.J_p_k[0:len(self.position)] + \
                                           self.a * (eps0 * np.matmul(E_matrix, self.chi_matrix)[0:len(self.position)] - self.P_k[0:len(self.position)])

        self.grid.J_p[self.position[0]:(self.position[-1] + 1)] = np.sum(self.J_p_k, axis=1)

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
