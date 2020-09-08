from .grid import Grid
from werkzeug.utils import cached_property
import numpy as np
import time
from .backend import backend as bd
from .backend import NumpyBackend
from .constants import eps0, c0
from .type_management import ndarray


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

    def step_P_tilde(self):
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
        self.gamma = bd.array(gamma)
        self.w_0 = bd.array(w0)
        self.chi_1 = chi_1
        self.chi_2 = chi_2
        self.chi_3 = chi_3
        self.J_p_k = None
        self.P_k = None
        self.P_tilde = None
        self.start_of_media = None

    def _allocate_P_tilde(self):
        self.P_tilde = bd.zeros((len(self.position), len(self.chi_1)))

    def _allocate_E_arrays(self):
        self.E_1 = bd.zeros(len(self.position))
        self.E_2 = bd.zeros(len(self.position))
        self.E_3 = bd.zeros(len(self.position))
        self.E_matrix = bd.empty((len(self.position), 3))

    def _allocate_J_p_k(self):
        self.J_p_k = bd.zeros((len(self.position), len(self.chi_1)))

    def _allocate_P_k(self):
        self.P_k = bd.zeros((len(self.position), len(self.chi_1)))


    @cached_property
    def a(self):
        a = bd.zeros(len(self.w_0))
        a = (self.grid.dt * self.w_0 ** 2) / (1 + self.gamma/2 * self.grid.dt)
        return a


    @cached_property
    def b(self):
        b = bd.zeros(len(self.w_0))
        b = (1 - self.gamma/2 * self.grid.dt) / (1 + self.gamma/2 * self.grid.dt)
        return b

    @cached_property
    def chi_matrix(self):
        chi_m = bd.array([self.chi_1, self.chi_2, self.chi_3])
        return chi_m

    '''@property
    def P_Tilde(self):
        E_matrix = np.zeros(shape=(len(self.position), 3))
        for ind, pos in zip(range(len(self.position)), self.position):
            E_matrix[ind] = [self.grid.Ez[pos] ** i for i in range(1, 4)]

        return eps0 * np.matmul(E_matrix, self.chi_matrix)'''

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

    def n_real(self, omega):
        return np.sqrt((np.abs(self.epsilon_complex(omega)) + self.epsilon_real(omega)) / 2)

    def n_kappa(self, omega):
        return np.sqrt((np.abs(self.epsilon_complex(omega)) - self.epsilon_real(omega)) / 2)

    def group_velocity(self, omega):
        n_1 = self.n_real(omega + 0.000001*omega)
        n_2 = self.n_real(omega - 0.000001*omega)
        diff_n = (n_1 - n_2)/(0.000002*omega)
        return (c0)/(self.n_real(omega) + omega*diff_n)

    def step_P_tilde(self):
        #start_time = time.time()
        if self.J_p_k is None:
            self._allocate_E_arrays()
            self._allocate_P_tilde()
            self._allocate_J_p_k()
            self._allocate_P_k()

        '''self.E_1[:] = self.grid.Ez[self.position[0]:(self.position[-1] + 1)]
        self.E_2[:] = self.E_1[:] ** 2
        self.E_3[:] = self.E_1[:] ** 3
        self.E_matrix = np.transpose(np.array([self.E_1, self.E_2, self.E_3]))'''
        #start_time = time.time()
        self.E_1 = self.grid.Ez[self.position[0]:(self.position[-1] + 1)]
        self.E_matrix = bd.stack((self.E_1, self.E_1 ** 2, self.E_1 ** 3), axis=1)
        #print("computed matrix in --- %s seconds ---" % (time.time() - start_time))
        #start_time = time.time()
        self.P_tilde = eps0 * bd.matmul(self.E_matrix, self.chi_matrix)
        #print("computed mamul in --- %s seconds ---" % (time.time() - start_time))

    def step_P(self):
        self.P_k = self.P_k + self.grid.dt * self.J_p_k
        self.grid.P[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.P_k, axis=1)
        '''if isinstance(bd, NumpyBackend):
            self.grid.P[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.P_k, axis=1)
        else:
            self.grid.P[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.P_k, dim=1)'''


    def step_J_p(self):
        '''self.J_p_k[0:len(self.position)][:] = self.b * self.J_p_k[0:len(self.position)] + \
                                           self.a * (self.P_tilde[0:len(self.position)] - self.P_k[0:len(self.position)])'''

        self.J_p_k = self.b * self.J_p_k + self.a * (self.P_tilde - self.P_k)

        self.grid.J_p[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.J_p_k, axis=1)
        '''if isinstance(bd, NumpyBackend):
            self.grid.J_p[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.J_p_k, axis=1)
        else:
            self.grid.J_p[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.J_p_k, dim=1)'''

class CentroRamanMedium(Vacuum):


    def __init__(self, name, permeability, conductivity, eps_inf, gamma_K, gamma_R, w0, chi_1, chi_3, alpha, wr):
        super().__init__()
        self.material_name = name
        self.model = 'Raman'
        self.eps = eps_inf
        self.mu = permeability
        self.conductivity = conductivity
        self.gamma_K = bd.array(gamma_K)
        self.gamma_R = bd.array(gamma_R)
        self.w_0 = bd.array(w0)
        self.w_R = bd.array(wr)
        self.chi_1 = chi_1
        self.chi_3 = chi_3
        self.J_p_k = None
        self.P_k = None
        self.G_k = None
        self.Q_k = None
        self.E_1 = None
        self.E_2 = None
        self.E_3 = None

        self.alpha = bd.array(alpha)
        self.start_of_media = None

    @cached_property
    def a(self):
        a = bd.zeros(len(self.w_0))
        a = (self.grid.dt * self.w_0 ** 2) / (1 + self.gamma_K/2 * self.grid.dt)
        return a

    @cached_property
    def b(self):
        b = bd.zeros(len(self.w_0))
        b = (1 - self.gamma_K/2 * self.grid.dt) / (1 + self.gamma_K/2 * self.grid.dt)
        return b

    @cached_property
    def c(self):
        c = bd.zeros(len(self.w_R))
        c = (self.w_R ** 2 * self.grid.dt) / (1 + self.grid.dt * self.gamma_R)
        return c

    @cached_property
    def d(self):
        d = bd.zeros(len(self.w_R))
        d = (1 - self.gamma_R * self.grid.dt) / (1 + self.gamma_R * self.grid.dt)
        return d

    @cached_property
    def chi_matrix(self):
        chi_m = bd.array([self.chi_1, self.chi_3])
        return chi_m

    def _allocate_J_p_k(self):
        self.J_p_k = bd.empty((len(self.position), len(self.chi_1)))

    def _allocate_P_k(self):
        self.P_k = bd.empty((len(self.position), len(self.chi_1)))

    def _allocate_P_tilde(self):
        self.P_Tilde = bd.empty((len(self.position), len(self.chi_1)))

    def _allocate_G_k(self):
        self.G_k = bd.empty((len(self.position), len(self.chi_1)))

    def _allocate_Q_k(self):
        self.Q_k = bd.empty((len(self.position), len(self.chi_1)))

    def _allocate_E_arrays(self):
        self.E_1 = bd.zeros(len(self.position))
        self.E_2 = bd.zeros(len(self.position))
        self.E_3 = bd.zeros(len(self.position))

    def epsilon_real(self, omega):
        eps_real = np.real(self.epsilon_complex(omega))
        return eps_real

    def epsilon_imag(self, omega):
        eps_imag = np.imag(self.epsilon_complex(omega))
        return eps_imag

    def epsilon_complex(self, omega):
        eps_complex = self.eps
        for w_k, gamma_k, chi_1_k in zip(self.w_0, self.gamma_K, self.chi_1):
            eps_complex += chi_1_k * (w_k ** 2) / (w_k ** 2 - omega ** 2 - 1j * gamma_k * omega)
        return eps_complex

    def n_real(self, omega):
        return np.sqrt((np.abs(self.epsilon_complex(omega)) + self.epsilon_real(omega)) / 2)

    def n_kappa(self, omega):
        return np.sqrt((np.abs(self.epsilon_complex(omega)) - self.epsilon_real(omega)) / 2)

    def group_velocity(self, omega):
        n_1 = self.n_real(omega + 0.000001 * omega)
        n_2 = self.n_real(omega - 0.000001 * omega)
        diff_n = (n_1 - n_2) / (0.000002 * omega)
        return (c0) / (self.n_real(omega) + omega * diff_n)


    def step_P_tilde(self):
        if self.J_p_k is None:
            self._allocate_P_k()
            self._allocate_J_p_k()
            self._allocate_Q_k()
            self._allocate_G_k()
            self._allocate_E_arrays()

        #start_time = time.time()
        self.E_1 = self.grid.Ez[self.position[0]:(self.position[-1] + 1)]
        self.E_2 = self.E_1 ** 2
        self.E_3 = self.E_1 ** 3
        #print("got E-fields in --- %s seconds ---" % (time.time() - start_time))

        first_term = bd.outer(self.E_1, bd.array(self.chi_1))
        second_term = bd.outer(self.E_3, (self.chi_matrix[1] * self.alpha))
        third_term_1_alpha = (1 - self.alpha) * self.chi_matrix[1]
        third_term = third_term_1_alpha * (self.Q_k * self.E_1[:, np.newaxis])
        self.P_Tilde = eps0 * (first_term + second_term + third_term)

    def step_P(self):
        self.P_k = self.P_k + self.grid.dt * self.J_p_k
        '''if isinstance(bd, NumpyBackend):
            self.grid.P[self.position[0]:(self.position[-1] + 1)] = np.sum(self.P_k, axis=1)
        else:
            self.grid.P[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.P_k, dim=1)'''


    def step_J_p(self):
        self.J_p_k = self.b * self.J_p_k + self.a * (self.P_Tilde - self.P_k)
        self.grid.J_p[self.position[0]:(self.position[-1] + 1)] = bd.sum(self.J_p_k, axis=1)


    def step_G(self):
        self.G_k = self.d * self.G_k + self.c * (self.E_2[:, np.newaxis] - self.Q_k)

    def step_Q(self):
        self.Q_k = self.Q_k + self.grid.dt * self.G_k

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
