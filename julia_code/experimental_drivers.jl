using Distances 
using Random
using LinearAlgebra
include("factorized_recovery.jl")
# function that estimates the best cut-off level based on a test set 
# obtained from 

function sparse_recovery(A, x, distance, ρ, max_level)
    opnorm_A = opnorm(A)
    return sparse_recovery(A, x, distance, ρ, max_level, opnorm_A)
end


function sparse_recovery(A, x, distance, ρ, max_level, opnorm_A)
    M,  basis, levels, diameters, centers, coloring = construct_measurement_matrix(x, ρ, distance, max_level)

    Ocol = A * M
    Orow = A' * M

    # Seems to work when either the input matrix is psd (presumably, symmetry is the main concern) or when the  basis used is the standard one
    L, U, D = LU_reconstruction(M, Ocol, Orow, basis, coloring, x, centers, distance)
 
    nm = norm(L * D * U - A) / norm(A)
    op_nm = opnorm(L * D * U - A) / opnorm_A
    n_meas =  size(M, 2)
    return L * D * U, nm, op_nm, n_meas
end


function app_opt_l_for_rho(A, x, distance, ρ, n_test)
    # Obtaining the geometric quantities that don't depend on rho
    ~,  basis, levels, diameters, centers, ~ = construct_measurement_matrix(x, 2.0, distance)

    level_list = minimum(levels) : maximum(levels)

    test_set = randn(size(A, 1), n_test) 
    test_set = test_set / qr(test_set).R
    opnorm_A = opnorm(A)

    # finding the best cutoff level
    err = Inf
    best_level = 0
    for level in level_list[2 : end]   
        @show level
        M,  basis, levels, diameters, centers, coloring = construct_measurement_matrix(x, ρ, distance, level)
        Ocol = A * M
        Orow = A' * M

        # Seems to work when either the input matrix is psd (presumably, symmetry is the main concern) or when the  basis used is the standard one
        L, U, D = LU_reconstruction(M, Ocol, Orow, basis, coloring, x, centers, distance)
        # new_err = opnorm(L * D * U - A)
        new_err = norm(L * D * U * test_set - A * test_set)
        if new_err < err
            @show best_level = level
            @show err = new_err
        else
            break
        end
        # @show size(M, 2)
        # @show norm(L * D * U - A) / norm(A)
        # @show opnorm(L * D * U - A) / opnorm(A) 
    end

    M,  basis, levels, diameters, centers, coloring = construct_measurement_matrix(x, ρ, distance, best_level)
    Ocol = A * M
    Orow = A' * M

    # Seems to work when either the input matrix is psd (presumably, symmetry is the main concern) or when the  basis used is the standard one
    L, U, D = LU_reconstruction(M, Ocol, Orow, basis, coloring, x, centers, distance)
 
    nm = norm(L * D * U - A) / norm(A)
    op_nm = opnorm(L * D * U - A) / opnorm_A
    n_meas =  size(M, 2)
    return nm, op_nm, n_meas, best_level
end

function app_opt_l_for_rho(A, x, distance, ρ, n_test, m_test)
    # Obtaining the geometric quantities that don't depend on rho
    ~,  basis, levels, diameters, centers, ~ = construct_measurement_matrix(x, 2.0, distance)

    level_list = minimum(levels) : maximum(levels)

    col_test_set = randn(size(A, 1), n_test) 
    col_test_set = col_test_set / qr(col_test_set).R

    row_test_set = randn(size(A, 2), m_test) 
    row_test_set = row_test_set / qr(row_test_set).R
    opnorm_A = opnorm(A)

    # finding the best cutoff level
    err = Inf
    best_level = 0
    for level in level_list[2 : end]   
        @show level
        M,  basis, levels, diameters, centers, coloring = construct_measurement_matrix(x, ρ, distance, level)
        Ocol = A * M
        Orow = A' * M

        # Seems to work when either the input matrix is psd (presumably, symmetry is the main concern) or when the  basis used is the standard one
        L, U, D = LU_reconstruction(M, Ocol, Orow, basis, coloring, x, centers, distance)
        # new_err = opnorm(L * D * U - A)
        new_err = sqrt(norm(L * D * U * col_test_set - A * col_test_set)^2 + norm(row_test_set' * L * D * U - row_test_set' * A )^2)
        if new_err < err
            @show best_level = level
            @show err = new_err
        else
            println("Breaking:")
            @show new_err
            break
        end
        # @show size(M, 2)
        # @show norm(L * D * U - A) / norm(A)
        # @show opnorm(L * D * U - A) / opnorm(A) 
    end

    M,  basis, levels, diameters, centers, coloring = construct_measurement_matrix(x, ρ, distance, best_level)
    Ocol = A * M
    Orow = A' * M

    # Seems to work when either the input matrix is psd (presumably, symmetry is the main concern) or when the  basis used is the standard one
    L, U, D = LU_reconstruction(M, Ocol, Orow, basis, coloring, x, centers, distance)
 
    nm = norm(L * D * U - A) / norm(A)
    op_nm = opnorm(L * D * U - A) / opnorm_A
    n_meas =  size(M, 2)
    return nm, op_nm, n_meas, best_level
end

# Uses the randomized methods described in HMT
# to construct random approximation to matrix
# A is input matrix
# rk is the target rank
# opn is used as an imput argument to avoid repeatedly 
# computing the opnorm
function rand_lr(A, rk)
    opn = opnorm(A)
    return rand_lr(A, rk, opn)
end 

# same as above, but with the third 
function rand_lr(A, rk, opn)
    Ω = randn(size(A, 1), rk)
    # rk matvecs
    Y = A * Ω
    Q = Matrix(qr(Y).Q)
    # rk matTvecs
    Approx = Q * (Q' * A)
    nm = norm(Approx - A) / norm(A)
    op_nm = opnorm(Approx - A) / opn
    n_meas = rk
    return Approx, nm, op_nm, n_meas
end 
