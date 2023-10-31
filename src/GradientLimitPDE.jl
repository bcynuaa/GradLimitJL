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
    inner_boundary::Boundary
    outer_boundary::Boundary
end

function GradientLimitPDE(
    background_grid::UniformSquareGrid,
    boundarys::Vector{Boundary}
)::GradientLimitPDE
    return GradientLimitPDE(
        background_grid,
        boundarys[1],
        boundarys[2]
    );
end

function isAntiClockwise(x1, y1, x2, y2, x3, y3)::Bool
    return (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) > 0;
end


function nodeInsideBoundary(x, y, boundary::Boundary)::Bool
    for i_node in 1: boundary.n_nodes
        i_node_next = mod1(i_node+1, boundary.n_nodes);
        if isAntiClockwise(
            x, y,
            boundary.xy[i_node, 1], boundary.xy[i_node, 2],
            boundary.xy[i_node_next, 1], boundary.xy[i_node_next, 2]
        ) != true
            return false;
        end
    end
    return true;
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
    known_index_s::Vector{Int64} = gradient_limit_pde.background_grid.known_index_s;
    additional_known_index_s::Vector{Int64} = Int64[];
    additional_known_value_s::Vector{Float64} = Float64[];
    for i_row = 2: gradient_limit_pde.background_grid.n-1
        for j_col = 2: gradient_limit_pde.background_grid.n-1
            x_pos, y_pos = positionAt(i_row, j_col, gradient_limit_pde.background_grid);
            # inner_boundary part
            sum_psi_J, sum_psi_I = calBoundaryContributions(x_pos, y_pos, gradient_limit_pde.inner_boundary);
            sum_psi_J_at_nodes[i_row, j_col] += sum_psi_J;
            sum_psi_I_s_at_nodes[i_row, j_col] += sum_psi_I;
            # outer_boundary part
            sum_psi_J, sum_psi_I = calBoundaryContributions(x_pos, y_pos, gradient_limit_pde.outer_boundary);
            sum_psi_J_at_nodes[i_row, j_col] += sum_psi_J;
            sum_psi_I_s_at_nodes[i_row, j_col] += sum_psi_I;
            if nodeInsideBoundary(x_pos, y_pos, gradient_limit_pde.inner_boundary)
                push!(additional_known_index_s, i_row + (j_col-1)*gradient_limit_pde.background_grid.n);
                push!(additional_known_value_s, sum_psi_I_s_at_nodes[i_row, j_col] / sum_psi_J_at_nodes[i_row, j_col]);
            end
        end
    end
    known_index_s = vcat(known_index_s, additional_known_index_s);
    known_index_values = vcat(gradient_limit_pde.background_grid.boundary_values, additional_known_value_s);
    unknown_index_s = setdiff(1: gradient_limit_pde.background_grid.n^2, known_index_s);
    M_diag::Vector{Float64} = vec(sum_psi_J_at_nodes);
    A_mat::SparseMatrixCSC{Float64, Int64} = gradient_limit_pde.background_grid.laplacian_operator - spdiagm(0 => M_diag);
    f_vec::Vector{Float64} = -vec(sum_psi_I_s_at_nodes);
    result::Vector{Float64} = zeros(gradient_limit_pde.background_grid.n^2);
    result[known_index_s] .= known_index_values;
    f_vec .-= A_mat * result;
    result[unknown_index_s] .= begin
    A_mat[unknown_index_s, unknown_index_s] \ f_vec[unknown_index_s];
    end
    return reshape(result, (gradient_limit_pde.background_grid.n, gradient_limit_pde.background_grid.n));
end

export GradientLimitPDE, solve;
export Boundary, UniformSquareGrid;
export createCircleBoundary, createSquareBoundary, calBoundaryContributions;

end # module GradientLimitPDEModule