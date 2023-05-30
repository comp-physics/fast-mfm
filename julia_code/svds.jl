using DelimitedFiles
using Distances 
using GLMakie
using LinearAlgebra
using Interpolations

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

interpol = false

# Writing the matrix path 
if matrix[1] == 'D'
    if interpol == false
        base_path = "julia_D/$(matrix)"
    else
        x = interpolate(x, BSpline(Linear()))[1 : (1 / interpol) : end]
        A = interpolate(A, BSpline(Linear()))[1 : (1 / interpol) : end, 1 : (1 / interpol) : end]
        base_path = "julia_D/$(matrix)_interpolated_$(interpol)"
    end
else 
    base_path = "julia_D/$(matrix)-D"
end

path = "$(base_path)_svd.csv"


op_nm_list = svdvals(A) / opnorm(A)
n_meas_list = 1 : size(A, 1)

open(path, "w") do io
    write(io, 
    "# svdvals of A, normalized by opnorm(A) 
    # opnorm_list, n_meas_list\n")
    writedlm(io, hcat(op_nm_list, n_meas_list), ',')
end


########################################
# Random low rank approximation
########################################
op_nm_list = Float64[]
nm_list = Float64[]
n_meas_list = Int[]

opn = opnorm(A)
for rk = vcat(Vector(1 : 20), 30:10:min(size(A, 1), min(200, size(A, 1))), Vector(250 : 50 : min(500, size(A, 1))))
    @show rk
    ~, nm, op_nm, n_meas = rand_lr(A, rk, opn)
    push!(op_nm_list, op_nm)
    push!(nm_list, nm)
    push!(n_meas_list, n_meas)
end

# Writing the matrix path 
path = "$(base_path)_rand_lr.csv"

open(path, "w") do io
    write(io, 
    "# random low rank approx of A 
    # nm_list, opnorm_list, n_meas_list\n")
    writedlm(io, hcat(nm_list, op_nm_list, n_meas_list), ',')
end