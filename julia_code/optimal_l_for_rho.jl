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
        path = "julia_D/$(matrix)_opt_l_for_rho_n_test_$(n_test)_m_test_$(m_test).csv"
    else
        
        # ongoing attempt at using exponential interpolation
        # x_int = interpolate(x, BSpline(Linear()))[1 : (1 / interpol) : end]
        # A_interpol = exponential_interpolation(A, x, x_int, 10.0 * minimum(diff(x)))

        # x = x_int
        # A = A_int


        x = interpolate(x, BSpline(Linear()))[1 : (1 / interpol) : end]
        A = interpolate(A, BSpline(Linear()))[1 : (1 / interpol) : end, 1 : (1 / interpol) : end]
        path = "julia_D/$(matrix)_interpolated_$(interpol)_opt_l_for_rho_n_test_$(n_test)_m_test_$(m_test).csv"
    end
else 
    path = "julia_D/$(matrix)-D_opt_l_for_rho_n_test_$(n_test)_m_test_$(m_test).csv"
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
best_level_list = Int[]
ρ_list = 1 : 0.5 : 30

for ρ in ρ_list
    nm, op_nm, n_meas, best_level = app_opt_l_for_rho(A, x, distance, ρ, n_test, m_test)
    push!(nm_list, nm)
    push!(op_nm_list, op_nm)
    push!(n_meas_list, n_meas)
    push!(best_level_list, best_level)
end

open(path, "w") do io
    write(io, 
    "# Optimal cut-off level determined based on n_test = $(n_test) and m_test = $(m_test) additional matrix-vector products 
    # rho_list, norm_list, opnorm_list, n_meas_list, best_level_list\n")
    writedlm(io, hcat(ρ_list, nm_list, op_nm_list, n_meas_list, best_level_list), ',')
end
