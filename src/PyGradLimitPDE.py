'''
 # @ author: Xu Ran | Bao Chenyu
 # @ date: 2023-10-30 18:44:25
 # @ license: MIT
 # @ language: python 3
 # @ description: python interface for GradLimitPDE
 '''

import os
import numpy as np
import juliacall

# read in current folder's "GradLimitPDE.jl" using juliacall
GradientLimitPDEJulia = juliacall.Main.include(os.path.join(os.path.dirname(__file__), "GradientLimitPDE.jl"))
JuliaVector = juliacall.Main.Vector

Boundary = GradientLimitPDEJulia.Boundary
createCircleBoundary = GradientLimitPDEJulia.createCircleBoundary
createSquareBoundary = GradientLimitPDEJulia.createSquareBoundary
calBoundaryContributions = GradientLimitPDEJulia.calBoundaryContributions

UniformSquareGrid = GradientLimitPDEJulia.UniformSquareGrid
positionAt = GradientLimitPDEJulia.positionAt
positionsAt = GradientLimitPDEJulia.positionsAt

def GradientLimitPDE(background_grid, boundarys):
    return GradientLimitPDEJulia.GradientLimitPDE(background_grid, boundarys[0], boundarys[1])
    pass

def solve(gradient_limit_pde):
    return np.array(GradientLimitPDEJulia.solve(gradient_limit_pde))
    pass