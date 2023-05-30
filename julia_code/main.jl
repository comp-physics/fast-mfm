using CompressingSolvers
using DelimitedFiles
using Distances 
using Random
using GLMakie
using LinearAlgebra

include("factorized_recovery.jl")
include("experimental_drivers.jl")
matrix = "D2121"
# matrix = "D1121"
# matrix = ""
# matrix = "100N-1Nvx-2Nvy"
# matrix = "2000N-3Nvx-2Nvy"
# matrix = "2000N-1Nvx-1Nvy"
if matrix[1] == 'D'
    A = readdlm("data/matrices/turbulent-channel-flow-park/$(matrix).csv", ',') 
    x = Matrix(readdlm("data/matrices/turbulent-channel-flow-park/y_faces.csv", ',')')[:]
else 
    A = readdlm("data/matrices/multivortex/case-$(matrix)/$(matrix)-D.csv", ',')
    x = readdlm("data/matrices/multivortex/case-$(matrix)/$(matrix)-x1_faces.csv", ',')[:]
end

# Writing the matrix path 
if matrix[1] == 'D'
    path = "julia_D/matrices/$(matrix)"
else 
    path = "julia_D/matrices/$(matrix)-D"
end



# Parameters of the algorithm.
# \rho is the oversampling radius,
# while max_measure is the maximal 
# number of measurements allowed
# probably want to replace the latter with a 
# parameter l that signifies the finest level
# of the multiresolution basis to be used. 
ρ = 2.5
max_level = 6


# The turbulent channel flow uses nonperiodic boundary conditions
if matrix[1] == 'D'
    distance = Euclidean()
else 
    # h is only needed for the periodic case
    h = x[2] - x[1]
    distance = PeriodicEuclidean(maximum(x) - minimum(x) + h)
end

A_peeling, nm, op_nm, n_meas = sparse_recovery(A, x, distance, ρ, max_level)

@show n_meas
@show nm
@show op_nm


open("$(path)_peeling_n_meas_$(n_meas).csv", "w") do io
    writedlm(io, A_peeling, ',')
end

A_rand_lr = rand_lr(A, n_meas)[1]

open("$(path)_rand_lr_n_meas_$(n_meas).csv", "w") do io
    writedlm(io, A_rand_lr, ',')
end



A_svd = svd(A)
A_svd = A_svd.U[:, 1:n_meas] * diagm(A_svd.S[1 : n_meas]) * A_svd.Vt[1 : n_meas, :]

open("$(path)_svd_n_meas_$(n_meas).csv", "w") do io
    writedlm(io, A_svd, ',')
end