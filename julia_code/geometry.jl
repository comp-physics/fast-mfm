using NearestNeighbors
using DataStructures
using Distances
using StaticArrays

# takes in the N point locations in d-dim space as an d \times N matrix
function create_coloring(x, ρ, distance)
    d = size(x, 1)
    N = size(x, 2)
    out = Vector{Int}[]
    x_list = Vector{SVector{d, Float64}}(undef, N)

    for k = 1 : N
        x_list[k] = SVector(Tuple(x[:, k]))
    end

    if typeof(distance) == Euclidean
        tree = KDTree(x_list, distance)
    else
        tree = BallTree(x_list, distance)
    end

    blocking = inrange(tree, x_list, ρ)

    assigned = falses(N)
    
    while count(assigned) < N
        nonblocked = trues(N) .&& .!assigned
        push!(out, Int[])
        while count(nonblocked) > 0
            # find the first index that is not blocked
            new_ind = findfirst(nonblocked)
            # add the new index to the current color
            push!(out[end], new_ind)
            # block all indices that are blocked by the new index
            nonblocked[blocking[new_ind]] .= false
            # note index as assigned
            assigned[new_ind] = true
        end
    end
    return out
end

function create_rhs(coloring)
    rhs = zeros(sum(length.(coloring)), length(coloring))
    for k = 1 : length(coloring)
        rhs[coloring[k], k] .= 1.0
    end
    return rhs
end

function deconstruct_measurement(measurement, coloring, x, distance)
    d = size(x, 1)
    N = size(x, 2)
    n_colors = length(coloring)
    @assert size(measurement, 2) == n_colors
    @assert size(measurement, 1) == N
    out = Vector{Int}[]
    x_list = Vector{SVector{d, Float64}}(undef, N)

    for k = 1 : N
        x_list[k] = SVector(Tuple(x[:, k]))
    end

    tree_function(x) = BallTree(x, distance)

#     function tree_function(x) 
#         @show typeof(x)
#         BallTree(x, distance)
#     end

    out = zeros(N, N)
    for (c_index, c) = enumerate(coloring)
        @show c
        new_color = zeros(N, length(c))
        tree = tree_function(x_list[c])
        nearest = nn(tree, x_list)[1]
        for k = 1 : N
            new_color[k, nearest[k]] = measurement[k, c_index]
        end
        out[:, c] .= new_color
    end
    return out
end
