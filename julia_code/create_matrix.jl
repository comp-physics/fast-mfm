using DelimitedFiles
using ArgParse
include("geometry.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "x_file"
        arg_type = String 
        help = "The file name for x"
        required=true
    "rho" 
        help = "The separation radius"
        arg_type = Float64
        required = true
    "meas_file" 
        help = "The file name ofr the measurements"
        arg_type = String
        required = true
    "output_file_name"
        arg_type = String 
        help = "the input file name"
        required = true
end
parsed_args = parse_args(s)
ρ = parsed_args["rho"]
x_file = parsed_args["x_file"]
meas_file = parsed_args["meas_file"]
out_file = parsed_args["output_file_name"]

x = Matrix(readdlm(x_file, ',')')
measurement = Matrix(readdlm(meas_file, ','))
# coloring = create_coloring(x, ρ, Euclidean())
coloring = create_coloring(x, ρ, PeriodicEuclidean((2 * x[end] - x[end - 1],)))

display(measurement)
@show size(measurement)

A = deconstruct_measurement(measurement, coloring, x, Euclidean())

display(A)

open(out_file, "w") do io
    writedlm(io, A, ',')
end



