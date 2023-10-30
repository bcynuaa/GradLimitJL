# -----
 # @ # author: Xu Ran | Bao Chenyu
 # @ # date: 2023-10-30 13:15:50
 # @ # license: MIT
 # @ # language: julia
 # @ # description: 2D gradient limit pde module
 # -----

module GradientLimitPDEModule

include("Boundary.jl");
include("BackgroundGrid.jl");

using SparseArrays;

using .BoundaryModule;
using .BackgroundGridModule;

struct GradientLimitPDE
    background_grid::UniformSquareGrid
    boundarys::Vector{Boundary}
end

"""
# Solve the gradient limit pde

∇²u=Mu+f → Au=f
"""
function solve(gradient_limit_pde::GradientLimitPDE)::Matrix{Float64}
    sum_psi_J_at_nodes::Matrix{Float64} = zeros(
        gradient_limit_pde.background_grid.n,
        gradient_limit_pde.background_grid.n
    );
    sum_psi_I_s_at_nodes::Matrix{Float64} = zeros(
        gradient_limit_pde.background_grid.n,
        gradient_limit_pde.background_grid.n
    );
    for i_row = 2: gradient_limit_pde.background_grid.n-1
        for j_col = 2: gradient_limit_pde.background_grid.n-1
            x_pos, y_pos = positionAt(i_row, j_col, gradient_limit_pde.background_grid);
            for boundary in gradient_limit_pde.boundarys
                sum_psi_J, sum_psi_I = calBoundaryContributions(x_pos, y_pos, boundary);
                sum_psi_J_at_nodes[i_row, j_col] += sum_psi_J;
                sum_psi_I_s_at_nodes[i_row, j_col] += sum_psi_I;
            end
        end
    end
    M_diag::Vector{Float64} = vec(sum_psi_J_at_nodes);
    A_mat::SparseMatrixCSC{Float64, Int64} = gradient_limit_pde.background_grid.laplacian_operator - spdiagm(0 => M_diag);
    f_vec::Vector{Float64} = -vec(sum_psi_I_s_at_nodes);
    result::Vector{Float64} = zeros(gradient_limit_pde.background_grid.n^2);
    result[gradient_limit_pde.background_grid.known_index_s] .= gradient_limit_pde.background_grid.boundary_values;
    f_vec .-= A_mat * result;
    result[gradient_limit_pde.background_grid.unknown_index_s] .= begin
    A_mat[gradient_limit_pde.background_grid.unknown_index_s, gradient_limit_pde.background_grid.unknown_index_s] \ f_vec[gradient_limit_pde.background_grid.unknown_index_s];
    end
    return reshape(result, (gradient_limit_pde.background_grid.n, gradient_limit_pde.background_grid.n));
end

export GradientLimitPDE, solve;
export Boundary, UniformSquareGrid;
export createCircleBoundary, createSquareBoundary;

end # module GradientLimitPDEModule