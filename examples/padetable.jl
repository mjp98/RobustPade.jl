using RobustPade, Plots, ColorSchemes

include("plotutil.jl")

# Examples from P. Gonnet, S. Guettel, and L. N. Trefethen, "Robust Pade approximation via SVD", SIAM Rev., 55:101-117, 2013.

M = 20
N = 20

# Colorscheme
basecolors = ColorSchemes.colorschemes[:glasbey_bw_minc_20_minl_30_n256][:]
padecolors = repeat(basecolors, (2 * ((N + 1) * (M + 1)) ÷ length(basecolors)) + 1)

# Figure 2 from SIAM Rev.
Fig2 = []
Fig2_f = [exp, cos, z -> (z^5 - 1) / (z^5 + 1), z -> log(5 + z^5)]
for f in Fig2_f
    push!(Fig2, plotpadetable(padetable(f, M, N), padecolors))
end
plot(Fig2...; layout=(2, 2))

# Figure 3 from SIAM Rev.
Fig3 = []
Fig3_f = [z -> 1 + z + z^8 + z^20 + z^30]
for f in Fig3_f
    push!(Fig3, plotpadetable(padetable(f, M, N), padecolors))
end
plot(Fig3...)

# Figure 4 from SIAM Rev.
Fig4 = []
n = 30;
ε = 1e-6;
Fig4_coeffs = ones(n) .+ ε * randn(n)
for tol in [1e-8, 1e-6, 10e-5]
    push!(Fig4, plotpadetable(padetable(Fig4_coeffs, M, N; tol), padecolors))
end
plot(Fig4...; layout=(1, 3))
