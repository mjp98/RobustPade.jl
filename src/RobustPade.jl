module RobustPade

using LinearAlgebra, PaddedViews, ToeplitzMatrices
using Polynomials
using TaylorSeries

export robustpade, robustpade_coefficients

include("util.jl")
include("pade.jl")

end
