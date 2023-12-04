# GradLimitJL

A repository for gradient limit pde in julia.

Xu Ran is a co-mate in research group of my research leader. This repository helps him solve a simple PDE problem in mesh-generating.

The original version was written in python of which the efficiency is not good enough. So I rewrite it in julia and makes it callable in python within `juliacall`.

In conclusion, this repository solves a possion equation with the form:

$$
\nabla^2 s = f
$$

where $f$ is the boundary mesh source term and $s$ is the mesh size scalar distribution all over the domain. This equation was solved on background grid with finite difference method.