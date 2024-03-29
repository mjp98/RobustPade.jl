module RobustPade

using LinearAlgebra, PaddedViews, ToeplitzMatrices
using Polynomials
using TaylorSeries

export robustpade, robustpade_table, robustpade_size

include("util.jl")
include("pade.jl")

end
