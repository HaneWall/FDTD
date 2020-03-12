'''
UNDER CONSTRUCTION!
Only works properly for courant = 0.5
'''

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
        self.arr_E = [0]
        self.arr_B = [0]

    def save_E(self):
        self.arr_E.append(self.grid.E[self.position + 1])     # 1st : [0] -> [0, E[pos+1]]

    def save_B(self):
        self.arr_B.append(self.grid.B[self.position + 1])

    def step_E(self):
        self.grid.E[self.position] = self.arr_E.pop(0)        # 2nd : [0, E[pos+1]] -> [E[pos+1]]

    def step_B(self):
        self.grid.B[self.position] = self.arr_B.pop(0)


class RightSideGridBoundary(Boundary):
    '''
    incoming waves from the left hand side are getting absorbed, only works for courant = 0.5
    '''

    def __init__(self):
        super().__init__()
        self.arr_E = [0]
        self.arr_B = [0]

    def save_E(self):
        self.arr_E.append(self.grid.E[self.position - 1])

    def save_B(self):
        self.arr_B.append(self.grid.B[self.position - 1])

    def step_E(self):
        self.grid.E[self.position] = self.arr_E.pop(0)

    def step_B(self):
        self.grid.B[self.position] = self.arr_B.pop(0)