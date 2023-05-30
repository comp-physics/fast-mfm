using DelimitedFiles
using Distances 
using GLMakie
using LinearAlgebra

######################################################################
# This file loops over rho from 1 to 20, and computes the optimal
# level of cutoff (with respect to the operator norm) 
# It then returns the optimal cutoff levels
######################################################################

include("factorized_recovery.jl")
include("experimental_drivers.jl")

# matrix = "D2121"
# matrix = "D1121"
# matrix = ""
# matrix = "100N-1Nvx-2Nvy"
matrix = "2000N-3Nvx-2Nvy"
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
    path = "julia_D/$(matrix)A_boussinesque.csv"
else 
    path = "julia_D/$(matrix)-D_boussinesque.csv"
end

nm_list = Float64[]
op_nm_list = Float64[]
n_meas_list = Int[]

A_boussinesque = diagm(A * ones(size(A, 1)))

push!(nm_list, norm(A_boussinesque - A) / norm(A))
push!(op_nm_list, opnorm(A_boussinesque - A) / opnorm(A_boussinesque))
push!(n_meas_list, 1)

open(path, "w") do io
    write(io, 
    "# boussineque approximation
    # norm_list, opnorm_list, n_meas_list\n")
    writedlm(io, hcat(nm_list, op_nm_list, n_meas_list), ',')
end
