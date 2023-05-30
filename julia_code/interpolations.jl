using Distances
function exponential_interpolation(A, x, x_int, l)
    N_old = length(x)
    N_new = length(x_int)
    x2d = reduce(hcat, [x[i]; x[j]] for i = 1 : N_old, j = 1 : N_old)

    x_int2d = reduce(hcat, [x_int[i]; x_int[j]] for i = 1 : N_new, j = 1 : N_new)

    Θ = exp.(-pairwise(Euclidean(), x2d)
 / l)
    # the "data" of the interpolation problem. 
    y = A[:]
    Θinv_y = Θ \ y

    # the offsets of the batches
    lims = Vector(1 : div(N_new ^ 2, 10) : (N_new ^ 2))
    push!(lims, N_new)

    out = Vector{Float64}(undef, N_new ^ 2)

    for k = 1 : (length(lims) - 1)
        println("Interpolating batch $(k)out of $(length(lims) - 1)")
        Θ_test_train = exp.(- pairwise(Euclidean(), x_int2d[:, lims[k]: lims[k + 1] - 1], x2d) / l)
        out[lims[k]: lims[k + 1] - 1] = Θ_test_train * Θinv_y
    end

    return reshape(out, N_new, N_new)
end 