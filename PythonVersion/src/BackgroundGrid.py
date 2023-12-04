'''
 # @ author: Xu Ran | Bao Chenyu
 # @ date: 2023-10-27 23:05:54
 # @ license: MIT
 # @ language: python 3
 # @ description: 2D background grid
 '''
 
import numpy as np
from scipy import sparse

diag_coeff = -4.0
offdiag_coeff = 1.0

class UniformSquareGrid:
    
    # x0: left bottom of the square
    # y0: left bottom of the square
    # edge_len: edge length of the square
    # n: number of nodes on one edge of the square
    
    def __init__(self) -> None:
        self.reset()
        pass
    
    def reset(self) -> None:
        self.x0 = 0.0
        self.y0 = 0.0
        self.edge_len = 1.0
        self.n = 2
        self.setParameters((self.x0, self.y0), self.edge_len, self.n)
        pass
    
    def setParameters(self, xy0, edge_len, n) -> None:
        self.x0, self.y0 = xy0
        self.edge_len = edge_len
        self.n = n
        self.__calParameters()
        self.__calLaplacianOperator()
        pass
    
    def position(self, i, j) -> tuple:
        x_i = self.x0 + i * self.delta
        y_j = self.y0 + j * self.delta
        return x_i, y_j
        pass
    
    def __calParameters(self) -> None:
        self.delta = self.edge_len / (self.n - 1)
        self.__getKnownIndexArray()
        self.__getUknownIndexArray()
        pass
    
    def __getKnownIndexArray(self) -> None:
        # known_index_array is a 1D array
        # the index of known_index_array is the boundary node index
        self.known_index_array = np.zeros(4 * (self.n - 1), dtype = int)
        self.known_index_array[0 : self.n] = np.arange(0, self.n)
        for i_row in range(1, self.n - 1):
            self.known_index_array[(i_row-1)*2 + self.n] = i_row * self.n
            self.known_index_array[(i_row-1)*2 + self.n + 1] = (i_row + 1) * self.n - 1
            pass
        self.known_index_array[-self.n : ] = np.arange((self.n-1)*self.n, self.n*self.n)
        pass
    
    def __getUknownIndexArray(self) -> None:
        self.unknown_index_array = np.setdiff1d(np.arange(0, self.n*self.n), self.known_index_array)
        pass
    
    def __calLaplacianOperator(self) -> None:
        # Laplacian operator
        # use sparse matrix to save memory
        # the matrix is a square matrix with shape (n*n, n*n)
        # 1. diag part
        diag = np.zeros(self.n*self.n) + diag_coeff / self.delta**2
        diag_row = np.arange(0, self.n*self.n)
        diag_col = np.arange(0, self.n*self.n)
        # 2. upper neighbour to diag
        upper_neighbour = np.zeros(self.n*self.n - 1) + offdiag_coeff / self.delta**2
        upper_neighbour_row = np.arange(0, self.n*self.n - 1)
        upper_neighbour_col = np.arange(1, self.n*self.n)
        # 3. lower neighbour to diag
        lower_neighbour = np.zeros(self.n*self.n - 1) + offdiag_coeff / self.delta**2
        lower_neighbour_row = np.arange(1, self.n*self.n)
        lower_neighbour_col = np.arange(0, self.n*self.n - 1)
        # 4. upper far to diag
        upper_far = np.zeros(self.n*self.n - self.n) + offdiag_coeff / self.delta**2
        upper_far_row = np.arange(0, self.n*self.n - self.n)
        upper_far_col = np.arange(self.n, self.n*self.n)
        # 5. lower far to diag
        lower_far = np.zeros(self.n*self.n - self.n) + offdiag_coeff / self.delta**2
        lower_far_row = np.arange(self.n, self.n*self.n)
        lower_far_col = np.arange(0, self.n*self.n - self.n)
        # 6. assemble
        data = np.concatenate((diag, upper_neighbour, lower_neighbour, upper_far, lower_far))
        row = np.concatenate((diag_row, upper_neighbour_row, lower_neighbour_row, upper_far_row, lower_far_row))
        col = np.concatenate((diag_col, upper_neighbour_col, lower_neighbour_col, upper_far_col, lower_far_col))
        # 7. generate sparse matrix
        self.laplacian_operator = sparse.coo_matrix((data, (row, col)), shape=(self.n*self.n, self.n*self.n)).tocsr()
        pass
    
    pass