using DelimitedFiles
using ArgParse
include("geometry.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "input_file_name"
        arg_type = String 
        help = "the input file name"
        required=true
    "rho" 
        help = "The physical separation within a color"
        arg_type = Float64
        required = true
    "output_file_name"
        arg_type = String 
        help = "the input file name"
        required = true
end
parsed_args = parse_args(s)
ρ = parsed_args["rho"]
in_file = parsed_args["input_file_name"]
out_file = parsed_args["output_file_name"]

x = Matrix(readdlm(in_file, ',')')
# coloring = create_coloring(x, ρ, Euclidean())
coloring = create_coloring(x, ρ, PeriodicEuclidean((2 * x[end] - x[end - 1],)))
rhs = create_rhs(coloring)
open(out_file, "w") do io
    writedlm(io, rhs, ',')
end



