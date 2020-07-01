from .constants import c0
from werkzeug.utils import cached_property
import numpy as np


class Boundary:

    def __init__(self):
        self.position = None
        self.grid = None

    def _place_into_grid(self, grid, index):
        if isinstance(index, int):
            self.grid = grid
            self.grid.boundaries.append(self)
            self.position = index
        elif isinstance(index, slice):
            raise KeyError('Not supporting slicing for boundaries!')


class LeftSideGridBoundary(Boundary):
    '''
    incoming waves from the right hand side are getting absorbed, only works for courant = 0.5
    '''

    def __init__(self):
        super().__init__()
        self.arr_Ez = [0]
        self.arr_Hy = [0]

    def save_Ez(self):
        self.arr_Ez.append(self.grid.Ez[self.position + 1])     # 1st : [0] -> [0, Ez[pos+1]]

    def save_Hy(self):
        self.arr_Hy.append(self.grid.Hy[self.position + 1])

    def step_Ez(self):
        self.grid.Ez[self.position] = self.arr_Ez.pop(0)        # 2nd : [0, Ez[pos+1]] -> [Ez[pos+1]]

    def step_Hy(self):
        self.grid.Hy[self.position] = self.arr_Hy.pop(0)


class RightSideGridBoundary(Boundary):
    '''
    incoming waves from the left hand side are getting absorbed, only works for courant = 0.5
    '''

    def __init__(self):
        super().__init__()
        self.arr_Ez = [0]
        self.arr_Hy = [0]

    def save_Ez(self):
        self.arr_Ez.append(self.grid.E[self.position - 1])

    def save_Hy(self):
        self.arr_Hy.append(self.grid.B[self.position - 1])

    def step_Ez(self):
        self.grid.Ez[self.position] = self.arr_Ez.pop(0)

    def step_Hy(self):
        self.grid.Hy[self.position] = self.arr_Hy.pop(0)

class LeftSideMur(Boundary):
    '''
    1st Order Mur
    incoming waves from the right hand side are getting absorbed
    '''

    def __init__(self):
        super().__init__()
        self.prev_Ez = [0, 0]
        self.prev_Hy = [0, 0]

    @cached_property
    def c_mat(self):
        return c0 / np.sqrt(self.grid.eps[self.position] * self.grid.mu[self.position])

    @cached_property
    def radiation_coeff(self):
        u = self.c_mat * self.grid.dt / self.grid.dx
        return (u - 1) / (u + 1)

    def save_Ez(self):
        self.prev_Ez[0] = self.grid.Ez[self.position]
        self.prev_Ez[1] = self.grid.Ez[self.position + 1]

    def step_Ez(self):
        self.grid.Ez[self.position] = self.prev_Ez[1] + self.radiation_coeff * (self.grid.Ez[self.position + 1] - self.prev_Ez[0])

    def save_Hy(self):
        pass
        #self.prev_Hy[0] = self.grid.Hy[self.position]
        #self.prev_Hy[1] = self.grid.Hy[self.position + 1]


    def step_Hy(self):
        pass
        #self.grid.Hy[self.position] = self.prev_Hy[1] + self.radiation_coeff * (self.grid.Hy[self.position + 1] - self.prev_Hy[0])

class RightSideMur(Boundary):
    '''
    1st Order Mur
    incoming waves from the left hand side are getting absorbed
    '''

    def __init__(self):
        super().__init__()
        self.prev_Ez = [0, 0]
        self.prev_Hy = [0, 0]

    @cached_property
    def c_mat(self):
        return c0 / np.sqrt(self.grid.eps[self.position] * self.grid.mu[self.position])

    @cached_property
    def radiation_coeff(self):
        u = self.c_mat * self.grid.dt / self.grid.dx
        return (u - 1) / (u + 1)

    def save_Ez(self):
        self.prev_Ez[0] = self.grid.Ez[self.position - 1]
        self.prev_Ez[1] = self.grid.Ez[self.position]

    def step_Ez(self):
        self.grid.Ez[self.position] = self.prev_Ez[0] + self.radiation_coeff * (self.grid.Ez[self.position - 1] - self.prev_Ez[1])

    def save_Hy(self):
        self.prev_Hy[0] = self.grid.Hy[self.position - 1]
        self.prev_Hy[1] = self.grid.Hy[self.position]

    def step_Hy(self):
        self.grid.Hy[self.position] = self.prev_Hy[0] + self.radiation_coeff * (self.grid.Hy[self.position - 1] - self.prev_Hy[1])