using GenericLinearAlgebra, TaylorSeries, RobustPade

m, n = 20, 20
coeffs = taylor_expand(sin, big(0.0); order=m + n + 1).coeffs
r = robustpade(coeffs, m, n)
x = big(π) / 12
r(x) - sin(x) |> abs
