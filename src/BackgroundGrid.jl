# -----
 # @ # author: Xu Ran | Bao Chenyu
 # @ # date: 2023-10-29 19:17:49
 # @ # license: MIT
 # @ # language: julia
 # @ # description: 2d background grid module
 # -----

module BackgroundGridModule

using SparseArrays;

const diag_value::Float64 = -4.0;
const offdiag_value::Float64 = 1.0;

struct UniformSquareGrid
    x_0::Float64 # left bottom x
    y_0::Float64 # left bottom y
    side_len::Float64
    n::Int64 # n nodes on each side
    boundary_values::Vector{Float64}
    delta::Float64
    x_s::Vector{Float64}
    y_s::Vector{Float64}
    unknown_index_s::Vector{Int64}
    known_index_s::Vector{Int64}
    laplacian_operator::SparseMatrixCSC{Float64, Int64}
end

function index2DToIndex1D(i::Int64, j::Int64, n::Int64)::Int64
    return (i-1) * n + j;
end

function UniformSquareGrid(xy0::Vector{Float64}, side_len::Float64, n::Int64, boundary_values)::UniformSquareGrid
    delta::Float64 = side_len / (n-1);
    boundary_values = zeros(4 * (n-1)) .+ boundary_values;
    x_s::Vector{Float64} =  LinRange(0, side_len, n) .+ xy0[1];
    y_s::Vector{Float64} =  LinRange(0, side_len, n) .+ xy0[2];
    unknown_index_s = index2DToIndex1D.(Vector(2: n-1), Vector(2: n-1)', n);
    unknown_index_s = reshape(unknown_index_s', (n-2)^2);
    # known index is the complement of unknown index
    all_index_s = index2DToIndex1D.(Vector(1: n), Vector(1: n)', n);
    all_index_s = reshape(all_index_s', n^2);
    known_index_s = setdiff(all_index_s, unknown_index_s);
    laplacian_operator::SparseMatrixCSC{Float64, Int64} = spdiagm(
        -1 => ones(n^2-1) .* offdiag_value,
        0 => ones(n^2) .* diag_value,
        1 => ones(n^2-1) .* offdiag_value,
        n => ones(n^2-n) .* offdiag_value,
        -n => ones(n^2-n) .* offdiag_value
    );
    laplacian_operator ./= delta^2;
    return UniformSquareGrid(
        xy0[1],
        xy0[2],
        side_len,
        n,
        boundary_values,
        delta,
        x_s,
        y_s,
        unknown_index_s,
        known_index_s,
        laplacian_operator
    );
end

function positionAt(i::Int64, j::Int64, grid::UniformSquareGrid)::Vector{Float64}
    return [grid.x_s[i], grid.y_s[j]];
end

function positionsAt(grid::UniformSquareGrid)::Matrix{Vector{Float64}}
    n::Int64 = grid.n;
    positions_mat::Matrix{Vector{Float64}} = positionAt.(Vector(1: n), Vector(1: n)', Ref(grid));
    return positions_mat;
end

export UniformSquareGrid, positionAt, positionsAt;

end # module BackgroundGridModule