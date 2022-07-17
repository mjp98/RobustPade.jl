
"""
    robustpade(
        coeffs::AbstractVector{T},
        m::Integer,
        n::Integer;
        tol::Real=eps(float(real(T)))
    )

computes the (m,n) Pade approximant to a function with Taylor coefficients 'coeffs' using SVD following Parchon et al. SIAM Review.

Adapted from the Chebfun implementation at https://github.com/chebfun/chebfun/blob/master/padeapprox.m (modified BSD license)


    robustpade(
        f::Function,
        m::Integer,
        n::Integer;
        x=0.
        tol::Real=100eps()
)

computes the (m,n) Pade approximant to a function f using TaylorSeries.taylor_expand to computed the Taylor coefficients.

"""
function robustpade(f::Function,m::Integer,n::Integer,x=0.;tol::Real = 100eps())
    taylorexpansion = taylor_expand(f, x; order=m+n+1)
    robustpade(taylorexpansion.coeffs, m,n; tol)
end

function robustpade(coeffs::AbstractVector,args...)
    robustpade(float.(coeffs), args...)
end

function robustpade(coeffs::AbstractVector{T}, m::Integer, n::Integer; tol::Real=epsreal(T)) where {T<:RealOrComplexFloat}

    @assert m >= 0 "m must be a non-negative integer."
    @assert n >= 0 "n must be a non-negative integer."
    @assert tol > 0 "tol must be a positive real."

    if length(coeffs) < m+n+1
        mn1,lc = m+n+1, length(coeffs)
        @warn "$mn1 coefficients required to determine ($m,$n) approximant. Padding coefficient vector of length $lc with zeros."
    end

    col = PaddedView(0, coeffs, (1:m+n+1,))
    row = PaddedView(0, first(coeffs, 1), (1:n+1,))

    if issubarrayzero(col, 1:m+1, tol)
        a = zeros(T, 1)
        b = ones(T, 1)
        return Polynomial(a)//Polynomial(b)
    else
        # Compute absolute tolerance.
        ts = tol * norm(col)
        # Form Toeplitz matrix.
        Z = Toeplitz(col, row)
        # Do diagonal hopping across block.
        m, n = robustpade_hop!(m, n, ts, Z)
        # Hopping finished. Now compute b and a.
        if iszero(n)
            a = first(col, m + 1)
            b = ones(T, 1)
        else
            @views C = Z[m+2:m+n+1, 1:n+1]
            _, _, V = svd0!(Matrix(C))
            # Null vector gives b.
            @views b = V[:, end]

            # Do final computation via reweighted QR for better zero preservation.
            @. b = abs(b) + sqrt(epsreal(T))
            Q, _ = qr!((C * Diagonal(b))')
            # Compensate for reweighting.
            @views b .*= Q[:, end]
            b ./= norm(b)
            # Multiplying gives a.
            @views a = Z[1:m+1, 1:n+1] * b
        end
        # Discard zeros
        i = firstnonzeroindex(b, tol)
        a = a[i:lastnonzeroindex(a, tol)]
        b = b[i:lastnonzeroindex(b, tol)]
        # Normalise
        a ./= first(b)
        b ./= first(b)
        return Polynomial(a)//Polynomial(b)
    end
end

function robustpade_hop!(m, n, ts, Z)
    # Special case n == 0.
    n == 0 && return m, n
    # Form Toeplitz matrix
    @views C = Z[(m+2):(m+n+1), 1:(n+1)]
    # Compute numerical rank. e.g. ρ = sum(z->z>ts,svdvals(C));
    ρ = rank(C, ts)
    # Compute rank-deficit.
    Δ = n - ρ
    # Break if full-rank.
    Δ == 0 && return m, n   #  m-Δ < 0 && return m,n
    # Decrease m, n if rank-deficient.
    return robustpade_hop!(m - Δ, n - Δ, ts, Z)
end
