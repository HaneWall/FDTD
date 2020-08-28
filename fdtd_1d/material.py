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

    def step_P(self):
        pass

    def step_J_p(self):
        pass

    def step_Q(self):
        pass

    def step_G(self):
        pass

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
        a[:] = (self.grid.dt * self.w_0[:] ** 2) / (1 + self.gamma[:]/2 * self.grid.dt)
        return a


    @cached_property
    def b(self):
        b = np.zeros(len(self.w_0))
        b[:] = (1 - self.gamma[:]/2 * self.grid.dt) / (1 + self.gamma[:]/2 * self.grid.dt)
        return b

    @cached_property
    def chi_matrix(self):
        chi_m = np.array([self.chi_1, self.chi_2, self.chi_3])
        return chi_m


    def P_Tilde(self):
        E_matrix = np.zeros(shape=(len(self.position), 3))
        for ind, pos in zip(range(len(self.position)), self.position):
            E_matrix[ind] = [self.grid.Ez[pos] ** i for i in range(1, 4)]

        return eps0 * np.matmul(E_matrix, self.chi_matrix)

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


        self.J_p_k[0:len(self.position)] = self.b * self.J_p_k[0:len(self.position)] + \
                                           self.a * (self.P_Tilde()[0:len(self.position)] - self.P_k[0:len(self.position)])

        self.grid.J_p[self.position[0]:(self.position[-1] + 1)] = np.sum(self.J_p_k, axis=1)


class CentroRamanMedium(Vacuum):


    def __init__(self, name, permeability, conductivity, eps_inf, gamma_K, gamma_R, w0, chi_1, chi_3, alpha, wr):
        super().__init__()
        self.material_name = name
        self.model = 'Raman'
        self.eps = eps_inf
        self.mu = permeability
        self.conductivity = conductivity
        self.gamma_K = np.array(gamma_K)
        self.gamma_R = np.array(gamma_R)
        self.w_0 = np.array(w0)
        self.w_R = np.array(wr)
        self.chi_1 = chi_1
        self.chi_3 = chi_3
        self.J_p_k = None
        self.P_k = None
        self.G_k = None
        self.Q_k = None


        self.first_order_term = None
        self.third_order_term = None

        self.alpha = np.array(alpha)
        self.start_of_media = None

    @cached_property
    def a(self):
        a = np.zeros(len(self.w_0))
        a[:] = (self.grid.dt * self.w_0[:] ** 2) / (1 + self.gamma_K[:]/2 * self.grid.dt)
        return a

    @cached_property
    def b(self):
        b = np.zeros(len(self.w_0))
        b[:] = (1 - self.gamma_K[:]/2 * self.grid.dt) / (1 + self.gamma_K[:]/2 * self.grid.dt)
        return b

    @cached_property
    def c(self):
        c = np.zeros(len(self.w_R))
        c[:] = (self.w_R[:] ** 2 * self.grid.dt) / (1 + self.grid.dt * self.gamma_R[:])
        return c

    @cached_property
    def d(self):
        d = np.zeros(len(self.w_R))
        d[:] = (1 - self.gamma_R[:] * self.grid.dt) / (1 + self.gamma_R[:] * self.grid.dt)
        return d

    @cached_property
    def chi_matrix(self):
        chi_m = np.array([self.chi_1, self.chi_3])
        return chi_m

    @cached_property
    def raman_factor(self):
        pre_raman = np.tile(np.transpose((1 - self.alpha) * self.chi_matrix[1]), (len(self.position), 1))
        return pre_raman

    @property
    def E_matrix(self):
        E_matrix = np.zeros(shape=(len(self.position), 3))
        for ind, pos in zip(range(len(self.position)), self.position):
            E_matrix[ind] = [self.grid.Ez[pos] ** i for i in [1, 2, 3]]
        return E_matrix

    @property
    def P_Tilde(self):
        if self.first_order_term is None:
            self._allocate_first_order_term()
            self._allocate_third_order_term()

        E_1_vector = np.transpose(self.E_matrix)[0]
        E_1_matrix = np.transpose(np.tile(E_1_vector, (len(self.chi_1), 1)))
        E_3_vector = np.transpose(self.E_matrix)[2]

        self.first_order_term = np.outer(np.transpose(self.chi_matrix[0]), E_1_vector)
        self.third_order_term = np.outer(np.transpose(self.alpha * self.chi_matrix[1]), E_3_vector)


        raman_term = self.raman_factor * self.Q_k * E_1_matrix


        return eps0 * (np.transpose(self.first_order_term) + np.transpose(self.third_order_term) + raman_term)

    def _allocate_J_p_k(self):
        self.J_p_k = np.zeros(shape=(len(self.position), len(self.chi_1)))

    def _allocate_P_k(self):
        self.P_k = np.zeros(shape=(len(self.position), len(self.chi_1)))

    def _allocate_G_k(self):
        self.G_k = np.zeros(shape=(len(self.position), len(self.chi_1)))

    def _allocate_Q_k(self):
        self.Q_k = np.zeros(shape=(len(self.position), len(self.chi_1)))

    def _allocate_first_order_term(self):
        self.first_order_term = np.zeros(shape=(len(self.chi_1), len(self.position)))

    def _allocate_third_order_term(self):
        self.third_order_term = np.zeros(shape=(len(self.chi_1), len(self.position)))


    def step_P(self):
        if self.J_p_k is None:
            self._allocate_P_k()
            self._allocate_J_p_k()
            self._allocate_G_k()
            self._allocate_Q_k()

        self.P_k[0:len(self.position)] = self.P_k[0:len(self.position)] + self.grid.dt * self.J_p_k[0:len(self.position)]
        self.grid.P[self.position[0]:(self.position[-1] + 1)] = np.sum(self.P_k, axis=1)

    def step_J_p(self):
        if self.J_p_k is None:
            self._allocate_P_k()
            self._allocate_J_p_k()
            self._allocate_G_k()
            self._allocate_Q_k()

        self.J_p_k[0:len(self.position)] = self.b * self.J_p_k[0:len(self.position)] + \
                                           self.a * (self.P_Tilde[0:len(self.position)] - self.P_k[0:len(self.position)])

        self.grid.J_p[self.position[0]:(self.position[-1] + 1)] = np.sum(self.J_p_k, axis=1)

    def step_G(self):
        if self.J_p_k is None:
            self._allocate_P_k()
            self._allocate_J_p_k()
            self._allocate_G_k()
            self._allocate_Q_k()

        self.G_k[0:len(self.position)] = self.d *self.G_k[0:len(self.position)] + self.c * \
                                         (np.transpose(np.tile(np.transpose(self.E_matrix)[1], (len(self.chi_1), 1))) - self.Q_k[0:len(self.position)])

    def step_Q(self):
        if self.J_p_k is None:
            self._allocate_P_k()
            self._allocate_J_p_k()
            self._allocate_G_k()
            self._allocate_Q_k()

        self.Q_k[0:len(self.position)] = self.Q_k[0:len(self.position)] + self.grid.dt * self.G_k[0:len(self.position)]

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
