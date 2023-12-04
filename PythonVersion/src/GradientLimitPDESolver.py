'''
 # @ author: Xu Ran | Bao Chenyu
 # @ date: 2023-10-28 00:21:58
 # @ license: MIT
 # @ language: python 3
 # @ description: 2D gradient limit pde solver
 '''

import numpy as np
import scipy as scp
from BackgroundGrid import UniformSquareGrid
from Boundary import *

class GradientLimitPDESolver:
    
    # gradient limit pde solver
    # gradient limit pde is a pde like:
    # $\nabla^2 u = G(x,y,u)$
    # on cartesian coordinate system it can be written as:
    # $(\nabla^2-\Lambda)u = f \to Au = f$
    
    def __init__(self) -> None:
        self.reset()
        pass
    
    def reset(self) -> None:
        self.background_grid: UniformSquareGrid = UniformSquareGrid()
        self.boundary_list: list[Boundary] = []
        self.boundary_values: np.ndarray = np.zeros(len(self.background_grid.known_index_array))
        self.result: np.ndarray = np.zeros((self.background_grid.n, self.background_grid.n))
        pass
    
    def addBoundary(self, boundary: Boundary) -> None:
        self.boundary_list.append(boundary)
        pass
    
    def setBoundaryValues(self, boundary_values) -> None:
        self.boundary_values = np.zeros(len(self.background_grid.known_index_array)) + boundary_values
        pass
    
    def __calSourceCoeff(self) -> np.ndarray:
        # calculate source coefficient $\Lambda$
        source_coeff = np.zeros((self.background_grid.n, self.background_grid.n))
        for i_row in range(self.background_grid.n):
            for j_col in range(self.background_grid.n):
                x_i, y_j = self.background_grid.position(i_row, j_col)
                for boundary in self.boundary_list:
                    source_coeff[i_row, j_col] += boundary.calSumPsiJ(x_i, y_j)
                    pass
                pass
            pass
        return source_coeff.T.flatten() # psi_n J_n
        # ! here .T is a trick
        pass
    
    def __calSource(self) -> np.ndarray:
        source = np.zeros((self.background_grid.n, self.background_grid.n))
        for i_row in range(self.background_grid.n):
            for j_col in range(self.background_grid.n):
                x_i, y_j = self.background_grid.position(i_row, j_col)
                for boundary in self.boundary_list:
                    source[i_row, j_col] += boundary.calSumPsiI(x_i, y_j)
                    pass
                pass
            pass
        return -source.T.flatten() # -psi_n I_n
        # ! here .T is a trick
        pass
    
    def solve(self) -> None:
        A = self.background_grid.laplacian_operator - scp.sparse.diags(self.__calSourceCoeff())
        f = self.__calSource()
        result = np.zeros(self.background_grid.n**2)
        result[self.background_grid.known_index_array] = self.boundary_values
        f_unknown = (f - A.dot(result))[self.background_grid.unknown_index_array]
        result[self.background_grid.unknown_index_array] = scp.sparse.linalg.spsolve(
            A[self.background_grid.unknown_index_array][:, self.background_grid.unknown_index_array],
            f_unknown
        )
        self.result = np.array(result).reshape(self.background_grid.n, self.background_grid.n).T
        # ! here .T is a trick
        pass
    
    pass