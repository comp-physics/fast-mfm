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

n_test = 2
m_test = 2

matrix = "D2121"
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
        path = "julia_D/$(matrix)_sweeping.csv"
    else
        x = interpolate(x, BSpline(Linear()))[1 : (1 / interpol) : end]
        A = interpolate(A, BSpline(Linear()))[1 : (1 / interpol) : end, 1 : (1 / interpol) : end]
        path = "julia_D/$(matrix)_interpolated_$(interpol)_sweeping.csv"
    end
else 
    path = "julia_D/$(matrix)-D_sweeping.csv"
end


# The turbulent channel flow uses nonperiodic boundary conditions
if matrix[1] == 'D'
    distance = Euclidean()
else 
    # h is only needed for the periodic case
    h = x[2] - x[1]
    distance = PeriodicEuclidean(maximum(x) - minimum(x) + h)
end

nm_list = Float64[]
op_nm_list = Float64[]
n_meas_list = Int[]
max_level_list = Int[]
ρ_list = Float64[]


~,  basis, levels, diameters, centers, ~ = construct_measurement_matrix(x, 2.0, distance)
opnorm_A = opnorm(A)

for ρ in 1 : 1.0: 20
    @show ρ
    for max_level = 1 : maximum(levels)
        @show max_level
        ~, nm, op_nm, n_meas = sparse_recovery(A, x, distance, ρ, max_level, opnorm_A)
        push!(nm_list, nm)
        push!(op_nm_list, op_nm)
        push!(n_meas_list, n_meas)
        push!(ρ_list, ρ)
        push!(max_level_list, max_level)
    end
end

open(path, "w") do io
    write(io, 
    "# accuracy of sparse recovery for different rho and max_level
    # rho_list, norm_list, opnorm_list, n_meas_list, max_level_list\n")
    writedlm(io, hcat(ρ_list, nm_list, op_nm_list, n_meas_list, max_level_list), ',')
end
