module RobustPade

using LinearAlgebra, PaddedViews, ToeplitzMatrices
using Polynomials
using TaylorSeries

export robustpade

include("util.jl")
include("pade.jl")

end
