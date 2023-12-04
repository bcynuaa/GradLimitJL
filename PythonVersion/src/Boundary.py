'''
 # @ author: Xu Ran | Bao Chenyu
 # @ date: 2023-10-27 20:56:09
 # @ license: MIT
 # @ language: python 3
 # @ description: 2D boundary
 '''

import numpy as np
from scipy import integrate

dim = 2
epsilon = 1e-10
epsilon2 = 1e-10

class Boundary:
    
    # ---------------------------------------------------------------------------------------------
    
    def __init__(self) -> None:
        self.reset()
        pass
    
    def allocate(self, n_nodes) -> None:
        self.n_nodes = n_nodes
        self.n_edges = n_nodes
        self.xy = np.zeros((n_nodes, dim))
        self.nodal_space = np.zeros(n_nodes)
        self.nodal_psi = np.zeros(n_nodes)
        self.linear_psi = np.zeros(n_nodes)
        pass
    
    def reset(self):
        self.allocate(0)
        pass
    
    # ---------------------------------------------------------------------------------------------
    
    def setNodesNumber(self, n_nodes) -> None:
        self.allocate(n_nodes)
        pass
    
    def setNodesPosition(self, xy) -> None:
        self.xy = np.zeros((self.n_nodes, dim)) + xy
        pass
    
    def setNodalSpace(self, nodal_space) -> None:
        self.nodal_space = np.zeros(self.n_nodes) + nodal_space
        pass
    
    def setNodalPsi(self, nodal_psi) -> None:
        self.nodal_psi = np.zeros(self.n_nodes) + nodal_psi
        pass
    
    def setLinearPsi(self, linear_psi) -> None:
        self.linear_psi = np.zeros(self.n_nodes) + linear_psi
        pass
        
    
    # ---------------------------------------------------------------------------------------------
    
    def calSumPsiJNodalPart(self, x, y) -> float:
        nodal_part = 0.0
        for i_node in range(self.n_nodes):
            x_i, y_i = self.xy[i_node]
            nodal_psi_i = self.nodal_psi[i_node]
            Ji = 1 / ( (x-x_i)**2 + (y-y_i)**2 + epsilon2 )
            nodal_part += nodal_psi_i * Ji
            pass
        return nodal_part
        pass
    
    def calSumPsiJLinearPart(self, x, y) -> float:
        linear_part = 0.0
        for i_edge in range(self.n_edges):
            x_i, y_i = self.xy[i_edge]
            x_j, y_j = self.xy[(i_edge+1)%self.n_edges]
            linear_psi_i = self.linear_psi[i_edge]
            edge_len = np.sqrt( (x_j-x_i)**2 + (y_j-y_i)**2 )
            # def funcForIntegral(t):
            #     xi = x_i + t * (x_j - x_i)
            #     eta = y_i + t * (y_j - y_i)
            #     return 1 / ((x-xi)**2 + (y-eta)**2 + epsilon2)
            #     pass
            # integral_part = integrate.quad(funcForIntegral, 0, 1)[0]
            # linear_part += linear_psi_i * integral_part
            # ! why following code is wrong?
            # ? wtf, it's right now
            if np.abs(x_j-x_i) < epsilon:
                def funcForIntegral(eta) -> float:
                    k_inv = (x_j - x_i) / (y_j - y_i)
                    xi = k_inv * (eta - y_i) + x_i
                    return 1 / ( (x-xi)**2 + (y-eta)**2 + epsilon2 )
                    pass
                integral_part = integrate.quad(funcForIntegral, y_i, y_j)[0] / (y_j - y_i)
                pass
            else:
                def funcForIntegral(xi) -> float:
                    k = (y_j - y_i) / (x_j - x_i)
                    eta = k * (xi - x_i) + y_i
                    return 1 / ( (x-xi)**2 + (y-eta)**2 + epsilon2 )
                    pass
                integral_part = integrate.quad(funcForIntegral, x_i, x_j)[0] / (x_j - x_i)
                pass
            linear_part += linear_psi_i * integral_part
            pass
        return linear_part
        pass
    
    # calculate $\sum_{n=1}^N \psi_n J_n$
    def calSumPsiJ(self, x, y) -> float:
        nodal_part = self.calSumPsiJNodalPart(x, y)
        linear_part = self.calSumPsiJLinearPart(x, y)
        return nodal_part + linear_part
        pass
    
    # ---------------------------------------------------------------------------------------------
    
    def calSumPsiINodalPart(self, x, y) -> float:
        nodal_part = 0.0
        for i_node in range(self.n_nodes):
            x_i, y_i = self.xy[i_node]
            nodal_space_i = self.nodal_space[i_node]
            nodal_psi_i = self.nodal_psi[i_node]
            Ji = 1 / ( (x-x_i)**2 + (y-y_i)**2 + epsilon2 ) * nodal_space_i
            nodal_part += nodal_psi_i * Ji
            pass
        return nodal_part
        pass
    
    def calSumPsiILinearPart(self, x, y) -> float:
        linear_part = 0.0
        for i_edge in range(self.n_edges):
            x_i, y_i = self.xy[i_edge]
            x_j, y_j = self.xy[(i_edge+1)%self.n_edges]
            s_i = self.nodal_space[i_edge]
            s_j = self.nodal_space[(i_edge+1)%self.n_edges]
            linear_psi_i = self.linear_psi[i_edge]
            edge_len = np.sqrt( (x_j-x_i)**2 + (y_j-y_i)**2 )
            integral_part = 0.0
            # def funcForIntegral(t):
            #     xi = x_i + t * (x_j - x_i)
            #     eta = y_i + t * (y_j - y_i)
            #     s = s_i + t * (s_j - s_i)
            #     return s / ((x-xi)**2 + (y-eta)**2 + epsilon2)
            #     pass
            # integral_part = integrate.quad(funcForIntegral, 0, 1)[0]
            # linear_part += linear_psi_i * integral_part
            # ! why following code is wrong?
            # ? wtf, it's right now
            if np.abs(x_j-x_i) < epsilon:
                def funcForIntegral(eta) -> float:
                    k_inv = (x_j - x_i) / (y_j - y_i)
                    xi = k_inv * (eta - y_i) + x_i
                    s = (s_j - s_i) / (y_j - y_i) * (eta - y_i) + s_i
                    return 1 / ( (x-xi)**2 + (y-eta)**2 + epsilon2 ) * s
                    pass
                integral_part = integrate.quad(funcForIntegral, y_i, y_j)[0] / (y_j - y_i)
                pass
            else:
                def funcForIntegral(xi) -> float:
                    k = (y_j - y_i) / (x_j - x_i)
                    eta = k * (xi - x_i) + y_i
                    s = (s_j - s_i) / (x_j - x_i) * (xi - x_i) + s_i
                    return 1 / ( (x-xi)**2 + (y-eta)**2 + epsilon2 ) * s
                    pass
                integral_part = integrate.quad(funcForIntegral, x_i, x_j)[0] / (x_j - x_i)
                pass
            linear_part += linear_psi_i  * integral_part
            pass
        return linear_part
        pass
    
    def calSumPsiI(self, x, y) -> float:
        nodal_part = self.calSumPsiINodalPart(x, y)
        linear_part = self.calSumPsiILinearPart(x, y)
        return nodal_part + linear_part
        pass
    
    # ---------------------------------------------------------------------------------------------
    
    pass

class CircleBoundary(Boundary):
    
    def __init__(self) -> None:
        super().__init__()
        pass
    
    def setCircle(self, center, radius, n_nodes) -> None:
        self.setNodesNumber(n_nodes)
        theta = np.linspace(0, 2*np.pi, n_nodes, endpoint=False)
        self.setNodesPosition(center + radius * np.array([np.cos(theta), np.sin(theta)]).T)
        self.center = center
        self.radius = radius
        pass
    
    pass

class SquareBoundary(Boundary):
    
    def __init__(self) -> None:
        super().__init__()
        pass
    
    def setSquare(self, xy0, edge_len, n_nodes_each_edge) -> None:
        self.setNodesNumber(4 * (n_nodes_each_edge-1) )
        dx = edge_len / (n_nodes_each_edge-1)
        self.x0, self.y0 = xy0
        # edge 1
        self.xy[0:n_nodes_each_edge-1, 0] = self.x0 + np.arange(n_nodes_each_edge-1) * dx
        self.xy[0:n_nodes_each_edge-1, 1] = self.y0
        # edge 2
        self.xy[n_nodes_each_edge-1:2*n_nodes_each_edge-2, 0] = self.x0 + edge_len
        self.xy[n_nodes_each_edge-1:2*n_nodes_each_edge-2, 1] = self.y0 + np.arange(n_nodes_each_edge-1) * dx
        # edge 3
        self.xy[2*n_nodes_each_edge-2:3*n_nodes_each_edge-3, 0] = self.x0 + edge_len - np.arange(n_nodes_each_edge-1) * dx
        self.xy[2*n_nodes_each_edge-2:3*n_nodes_each_edge-3, 1] = self.y0 + edge_len
        # edge 4
        self.xy[3*n_nodes_each_edge-3:4*n_nodes_each_edge-4, 0] = self.x0
        self.xy[3*n_nodes_each_edge-3:4*n_nodes_each_edge-4, 1] = self.y0 + edge_len - np.arange(n_nodes_each_edge-1) * dx
        # other properties
        self.edge_len = edge_len
        self.n_nodes_each_edge = n_nodes_each_edge
        pass
    
    pass