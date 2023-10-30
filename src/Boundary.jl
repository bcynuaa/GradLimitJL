# -----
 # @ # author: Xu Ran | Bao Chenyu
 # @ # date: 2023-10-29 19:15:47
 # @ # license: MIT
 # @ # language: julia
 # @ # description: 2d boundary struct
 # -----

module BoundaryModule

const dim::Int64 = 2;
const epsilon::Float64 = 1e-10;

function distance(x_1, y_1, x_2, y_2)::Float64
    return sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2);
end

function distance2(x_1, y_1, x_2, y_2)::Float64
    return (x_1 - x_2)^2 + (y_1 - y_2)^2;
end

struct Boundary
    n_nodes::Int64
    n_edges::Int64
    xy::Matrix{Float64}
    edge_lens::Vector{Float64}
    nodal_sources::Vector{Float64}
    nodal_psi_s::Vector{Float64}
    linear_psi_s::Vector{Float64}
end

function Boundary(
    n_nodes::Int64,
    xy::Matrix{Float64},
    nodal_sources,
    nodal_psi_s,
    linear_psi_s
)::Boundary
    n_edges::Int64 = n_nodes;
    edge_lens::Vector{Float64} = zeros(n_edges);
    for i_node::Int64 = 1: n_nodes;
        i_node_next::Int64 = mod1(i_node + 1, n_nodes);
        edge_lens[i_node] = distance(
            xy[i_node, 1],
            xy[i_node, 2],
            xy[i_node_next, 1],
            xy[i_node_next, 2]
        );
    end
    nodal_sources::Vector{Float64} = zeros(n_nodes) .+ nodal_sources;
    nodal_psi_s::Vector{Float64} = zeros(n_nodes) .+ nodal_psi_s;
    linear_psi_s::Vector{Float64} = zeros(n_edges) .+ linear_psi_s;
    return Boundary(
        n_nodes,
        n_edges,
        xy,
        edge_lens,
        nodal_sources,
        nodal_psi_s,
        linear_psi_s
    );
end

function createCircleBoundary(
    n_nodes::Int64,
    radius::Float64,
    center_xy::Vector{Float64},
    nodal_sources,
    nodal_psi_s,
    linear_psi_s
)::Boundary
    theta::Vector{Float64} = Vector(range(0, stop=2*pi, length=n_nodes+1)[1:n_nodes]);
    xy::Matrix{Float64} = zeros(n_nodes, dim);
    xy[:, 1] = center_xy[1] .+ radius .* cos.(theta);
    xy[:, 2] = center_xy[2] .+ radius .* sin.(theta);
    n_edges::Int64 = n_nodes;
    edge_lens::Vector{Float64} = zeros(n_nodes) .+ 2*radius*sin(pi/n_nodes);
    nodal_sources::Vector{Float64} = zeros(n_nodes) .+ nodal_sources;
    nodal_psi_s::Vector{Float64} = zeros(n_nodes) .+ nodal_psi_s;
    linear_psi_s::Vector{Float64} = zeros(n_nodes) .+ linear_psi_s;
    return Boundary(
        n_nodes,
        n_edges,
        xy,
        edge_lens,
        nodal_sources,
        nodal_psi_s,
        linear_psi_s
    );
end

function createSquareBoundary(
    n_nodes_per_side::Int64,
    side_len::Float64,
    left_bottom_xy::Vector{Float64},
    nodal_sources,
    nodal_psi_s,
    linear_psi_s
)::Boundary
    n_nodes::Int64 = 4 * (n_nodes_per_side-1);
    n_edges::Int64 = n_nodes;
    xy::Matrix{Float64} = zeros(n_nodes, dim);
    edge_lens::Vector{Float64} = zeros(n_nodes);
    nodal_sources::Vector{Float64} = zeros(n_nodes) .+ nodal_sources;
    nodal_psi_s::Vector{Float64} = zeros(n_nodes) .+ nodal_psi_s;
    linear_psi_s::Vector{Float64} = zeros(n_nodes) .+ linear_psi_s;
    per_edge_len::Float64 = side_len / (n_nodes_per_side - 1);
    edge_lens .= per_edge_len;
    # allocate nodes anticlockwise around square from left bottom
    # side 1
    xy[1:n_nodes_per_side-1, 1] .= LinRange(0., side_len-per_edge_len, n_nodes_per_side-1);
    xy[1:n_nodes_per_side-1, 2] .= 0.;
    # side 2
    xy[n_nodes_per_side:2*(n_nodes_per_side-1), 1] .= side_len;
    xy[n_nodes_per_side:2*(n_nodes_per_side-1), 2] .= LinRange(0., side_len-per_edge_len, n_nodes_per_side-1);
    # side 3
    xy[2*(n_nodes_per_side-1)+1:3*(n_nodes_per_side-1), 1] .= LinRange(side_len, per_edge_len, n_nodes_per_side-1);
    xy[2*(n_nodes_per_side-1)+1:3*(n_nodes_per_side-1), 2] .= side_len;
    # side 4
    xy[3*(n_nodes_per_side-1)+1:4*(n_nodes_per_side-1), 1] .= 0.;
    xy[3*(n_nodes_per_side-1)+1:4*(n_nodes_per_side-1), 2] .= LinRange(side_len-per_edge_len, 0., n_nodes_per_side-1);
    # shift to left bottom
    xy[:, 1] .+= left_bottom_xy[1];
    xy[:, 2] .+= left_bottom_xy[2];
    return Boundary(
        n_nodes,
        n_edges,
        xy,
        edge_lens,
        nodal_sources,
        nodal_psi_s,
        linear_psi_s
    );
end

# actually, 2 points' precision is enough
# const gauss_points::Vector{Float64} = [-1/sqrt(3), 1/sqrt(3)] ./ 2 .+ 1/2;
# const gauss_weights::Vector{Float64} = [1., 1.] ./ 2;
const gauss_points::Vector{Float64} = [-0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526] ./ 2 .+ 1/2;
const gauss_weights::Vector{Float64} = [0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385] ./ 2;

function calBoundaryContributions(
    x, y,
    boundary::Boundary
)::Tuple{Float64, Float64}
    # nodal part
    distances2::Vector{Float64} = distance2.(boundary.xy[:, 1], boundary.xy[:, 2], x, y);
    nodal_psi_J_s::Vector{Float64} = 1 ./ (distances2 .+ epsilon) .* boundary.nodal_psi_s;
    nodal_psi_I_s::Vector{Float64} = nodal_psi_J_s .* boundary.nodal_sources;
    nodal_sum_psi_J::Float64 = sum(nodal_psi_J_s);
    nodal_sum_psi_I::Float64 = sum(nodal_psi_I_s);
    # linear part
    linear_psi_J_s::Vector{Float64} = zeros(boundary.n_edges);
    linear_psi_I_s::Vector{Float64} = zeros(boundary.n_edges);
    for i_edge::Int64 = 1: boundary.n_edges
        # get temp variables
        i_edge_next::Int64 = mod1(i_edge + 1, boundary.n_edges);
        x_1::Float64, y_1::Float64 = boundary.xy[i_edge, :];
        x_2::Float64, y_2::Float64 = boundary.xy[i_edge_next, :];
        s_1::Float64 = boundary.nodal_sources[i_edge];
        s_2::Float64 = boundary.nodal_sources[i_edge_next];
        # calculate linear integral using gauss quadrature
        dr::Vector{Float64} = [x_2 - x_1, y_2 - y_1];
        r_at_gauss_points::Matrix{Float64} = [x_1, y_1]' .+ gauss_points .* dr';
        r_to_xy_at_gauss_points::Matrix{Float64} = r_at_gauss_points .- [x, y]';
        ds::Float64 = s_2 - s_1;
        s_at_gauss_points::Vector{Float64} = gauss_points .* ds .+ s_1;
        r2_inv_at_gauss_points::Vector{Float64} = 1 ./ (sum(r_to_xy_at_gauss_points.^2, dims=2)[:,] .+ epsilon);
        linear_psi_J_s[i_edge] = gauss_weights' * r2_inv_at_gauss_points .* boundary.linear_psi_s[i_edge];
        linear_psi_I_s[i_edge] = gauss_weights' * (r2_inv_at_gauss_points .* s_at_gauss_points) .* boundary.linear_psi_s[i_edge];
    end
    linear_sum_psi_J::Float64 = sum(linear_psi_J_s);
    linear_sum_psi_I::Float64 = sum(linear_psi_I_s);
    sum_psi_J::Float64 = nodal_sum_psi_J + linear_sum_psi_J;
    sum_psi_I::Float64 = nodal_sum_psi_I + linear_sum_psi_I;
    return sum_psi_J, sum_psi_I;
end

export  Boundary,
        createCircleBoundary,
        createSquareBoundary;
export calBoundaryContributions;

end # module Boundary