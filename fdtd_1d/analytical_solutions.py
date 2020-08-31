
from scipy.integrate import odeint, complex_ode
from fdtd_1d.constants import c0
import numpy as np


class Qpm_Length_Analytical:

    def __init__(self, initial_values, mode='qpm'):
        self.omega_1 = 1.77e15
        self.omega_2 = 2 * 1.77e15
        #self.n_1 = np.abs(1.5745775436250355 + 1j * 1.5745775436250355)
        #self.n_2 = np.abs(1.6381145184887687 + 1j * 1.6381145184887687)
        self.n_1 = 2.226788917202639
        self.n_2 = 2.3166437687670887

        self.k_1 = self.n_1 * self.omega_1 / c0
        self.k_2 = self.n_2 * self.omega_2 / c0
        self.delta_k = 2 * self.k_1 - self.k_2
        self.lambda_qpm = 5.9e-6
        self.k_m = 2 * np.pi / self.lambda_qpm
        self.mode = mode
        self.d_eff = 15e-15
        self.A = np.array(initial_values)

    def _allocate_storage(self, steps):
        self.solution = np.empty(shape=(steps, 2), dtype='complex')

    def qpm_coupled_system(self, z, A):
        A_1 = A[0]
        A_2 = A[1]

        if self.mode == 'qpm':
            d_m = np.sign(np.sin(self.k_m * z)) * self.d_eff
        else:
            d_m = self.d_eff

        da1 = (2j * self.omega_1 * d_m) / (self.n_1 * c0) * A_2 * np.conj(A_1) * np.exp(-1j * self.delta_k * z)
        da2 = (1j * self.omega_2 * d_m) / (self.n_2 * c0) * A_1 ** 2 * np.exp(1j * self.delta_k * z)
        return [da1, da2]

    def integrate(self, end_z, dz):
        steps = int(end_z/dz)
        self.z = np.arange(0, stop=end_z, step=dz)
        self._allocate_storage(steps)
        f = complex_ode(self.qpm_coupled_system)
        f.set_initial_value(self.A)
        i = 0

        while f.successful() and f.t <= end_z:
            self.solution[i] = f.integrate(f.t + dz)
            i += 1

        self.A_1_solution = np.transpose(self.solution)[0]
        self.A_2_solution = np.transpose(self.solution)[1]
        return [self.A_1_solution, self.A_2_solution, self.z]


