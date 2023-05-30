using LinearAlgebra
using NearestNeighbors
using StaticArrays
using DataStructures

# This function subdivides an input set of indices into two 
# and forms the resulting wavelet basis vector.
# It then adds these vectors to basis_vectors,
# adds their diameter and center to 
function subdivide!(index_set, level, x, diameters, centers, levels, basis_vectors)
    n = length(index_set)
    N = length(x)
    # if n = 1, we are done (on the finest scale)
    if n == 1
        return 0
    # if n > 1, we go through another round of splitting
    else 
        level = level + 1 
        diameter = maximum(x[index_set]) - minimum(x[index_set])
        center = (maximum(x[index_set]) + minimum(x[index_set])) / 2

        # choosing a pivot to divide the array
        set_plus = index_set[findall(k -> x[k] ≤ center, index_set)]
        set_minus = index_set[findall(k -> x[k] > center, index_set)]

        # constructs the new basis vector, such that it is orthogonal to
        # the next higher scale
        basis_vector = zeros(N)
        basis_vector[set_plus] .= 1.0 / length(set_plus)
        basis_vector[set_minus] .= -1.0 / length(set_minus)
        basis_vector = basis_vector / norm(basis_vector)
        push!(diameters, diameter)
        push!(centers, center)
        push!(levels, level)
        push!(basis_vectors, basis_vector)
        subdivide!(set_plus, level, x, diameters, centers, levels, basis_vectors)
        subdivide!(set_minus, level, x, diameters, centers, levels, basis_vectors)
    end
end


# Input is vector of physical locations. 
# uses regular subdivision and thus assumes roughly uniform mesh. 
function create_multiscale_basis(x::Vector{Float64})
    N = length(x)
    diameters = Float64[]
    centers = Float64[]
    levels = Int[]
    basis_vectors = Vector{Float64}[]
    #creating the first, coarsest basis function on level 0 
    push!(diameters, maximum(x) - minimum(x))
    push!(centers, (maximum(x) + minimum(x)) / 2)
    push!(levels, 0)
    push!(basis_vectors, ones(N) / norm(ones(N)))
    subdivide!(collect(1 : N), 0, x, diameters, centers, levels, basis_vectors) 
    P = sortperm(levels)
    return reduce(hcat, basis_vectors[P]), levels[P], diameters[P], centers[P]
end

# computes the multicolor ordering within a level. 
# Centers and diameters are the vectors of centers and diameters of the 
# indices within the scale under consideration. 
# right now hardcoded dimension 1.
function multicolor_ordering(ρ, x, diameters, indices, distance)
    ρh = ρ * maximum(diameters)
    d = 1
    N = length(x)
    out = Vector{Int}[]
    x_list = Vector{SVector{d, Float64}}(undef, N)

    for k = 1 : N
        x_list[k] = SVector((x[k],))
    end

    if typeof(distance) == Euclidean
        tree = KDTree(x_list, distance)
    else
        tree = BallTree(x_list, distance)
    end

    blocking = inrange(tree, x_list, ρh)

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
    return [indices[color] for color in out]
end

function farthest_multicolor_ordering(ρ, x, diameters, indices, distance)
    # using the old input variable names for ρh
    ρh = ρ * maximum(diameters)
    d = 1

    tree_function(x) = BallTree(x, distance)

    N = length(x)
    out = Vector{Int}[]
    input_array = Vector{SVector{d, Float64}}(undef, N)

    for k = 1 : N
        input_array[k] = SVector((x[k],))
    end


    # Old function header
    # function construct_multicolor_ordering(input_array::AbstractVector, ρh::Real, tree_function)
    # Swapping out general RT for Float64
    # RT = real_type(eltype(input_array))
    RT = Float64

    # Vector (colors) of Vectors of supernodes  
    out = Vector{Int}[]

    assigned = falses(length(input_array))
    # While not all nodes are assigned to a color
    tree = tree_function(input_array)
    heap = MutableBinaryMaxHeap(fill(typemax(RT), length(input_array)))
    while !isempty(heap)
        # add new color
        push!(out, Int[])
        # reset heap values
        # Presently a little inefficient since we revisit every entry of the multicolor 
        # ordering. Still doesn't change the asymptotic complexity
        for k = 1 : length(input_array)
            # only assigns a value to points that have not already been assigned and thus
            # removed from the heap
            !assigned[k] && (heap[k] = typemax(RT)) 
        end
        # do furthest point sampling, finish if either the heap is empty, or there are no
        # sufficiently distant points left 
        while !isempty(heap) && (first(heap) == typemax(RT) || first(heap) > ρh)
            # get the id of the new pivot
            ~, top_id = top_with_handle(heap)
            # remove the new pivot from the heap
            pop!(heap); assigned[top_id] = true; push!(out[end], top_id)
            # add the new pivot to the 
            # obtain the affected nodes and their distance to the pivot 
            number_affected_nodes = length(inrange(tree, input_array[top_id], ρh))
            affected_nodes, distances = knn(tree, input_array[top_id], number_affected_nodes)
            # Update the distances of nodes that could be affedcted. 
            for (node, dist) in zip(affected_nodes, distances)
                # update distance if node is still in heap
                !assigned[node] && update!(heap, node, dist)
            end
            
        end
    end
    return [indices[color] for color in out]
end


function construct_measurement_matrix(x, ρ, distance, max_level=typemax(Int))
    N = length(x) 
    basis, levels, diameters, centers = create_multiscale_basis(x)
    coloring = Vector{Int}[]
    # checking that levels is ordered and consecutive
    @assert unique(levels) == collect(levels[1] : levels[end])
    for level in unique(levels)
        # find the set of indices corresponding to the level
        level_indices = findall(l -> l == level, levels)
        # Swapped out the ordering algorithm. The farthest point way of selecting colors doesn't seem helpful
        append!(coloring, multicolor_ordering(ρ, centers[level_indices], diameters[level_indices], level_indices, distance))
        # append!(coloring, farthest_multicolor_ordering(ρ, centers[level_indices], diameters[level_indices], level_indices, distance))
    end
    
    # We make sure that the dofs of the same color appear consecutively
    ordering = reduce(vcat, coloring) 
    reverse_ordering = similar(ordering)
    reverse_ordering[ordering] .= 1 : length(ordering)
    basis = basis[:, ordering] 
    # checking that the multicolor ordering does not mix levels
    # display(levels)
    # display(levels[ordering])
    # display(ordering)
    # display(coloring)
    @assert levels == levels[ordering]
    diameters = diameters[ordering]
    centers = centers[ordering]
    coloring = [reverse_ordering[c] for c in coloring]
    #check that the colors in the reordered coloring are consecutive 
    for c in coloring
        @assert c == collect(c[1] : c[end])
    end

    measurement_matrix = zeros(N, 0)
    for color in coloring
        measurement_matrix = hcat(measurement_matrix, sum(basis[:, color], dims=2))
    end

    last_index_of_level = findlast(k -> levels[k] ≤ max_level, 1 : length(levels))
    last_color = findlast(c -> minimum(c) ≤ last_index_of_level, coloring)

    return measurement_matrix[:, 1 : min(last_color, end)], basis, levels, diameters, centers, coloring[1 : min(last_color, end)]
end

function deconstruct_measurement(input_column, color_indices, x, centers, distance)
    d = 1
    N = length(x)
    out = Vector{Int}[]
    x_list = Vector{SVector{d, Float64}}(undef, N)
    center_list = Vector{SVector{d, Float64}}(undef, N)

    for k = 1 : N
        x_list[k] = SVector((x[k],))
        center_list[k] = SVector((centers[k],))
    end
    # can possibly use other trees, doesn't matter for now
    tree_function(x) = BallTree(x, distance)
    out = zeros(N, length(color_indices))
    tree = tree_function(center_list[color_indices])
    # obtaining the nearest colored index for each index
    nearest = nn(tree, x_list)[1]
    for k = 1 : N
        out[k, nearest[k]] = input_column[k]
    end
    return out
end

function LU_reconstruction(M, Ocol, Orow, basis, coloring, x, centers, distance)
    N = length(x)
    L = zeros(N, 0)
    U = zeros(0, N)
    D = zeros(0)

    # iterating over the colors
    for (color_index, color) in enumerate(coloring)
        new_cols = deconstruct_measurement(Ocol[:, color_index], color, x, centers, distance)
        new_rows = deconstruct_measurement(Orow[:, color_index], color, x, centers, distance)

        new_D = (basis[:, color]' * new_cols + new_rows' * basis[:, color]) / 2
        
        @assert new_D ≈ diagm(diag(new_D))
        new_D = diagm(1 ./ diag(new_D))
        append!(D, diag(new_D))

        L = hcat(L, new_cols)
        U = vcat(U, new_rows')

        Ocol -=  new_cols  * new_D * new_rows' * M
        Orow -=  new_rows  * new_D * new_cols' * M
    end
    return L, U, diagm(D)
end

# function reconstruct(M, O, coloring, x, centers, distance)
#     N = length(x)
#     Q = zeros(N, 0)
#     R = zeros(N, 0)
# 
# 
#     # iterating over the colors
#     for (color_index, color) in enumerate(coloring)
#         new_cols = deconstruct_measurement(O[:, color_index], color, x, centers, distance)
#         Q = hcat(Q, new_cols) 
# 
# 
# 
# 
# end

