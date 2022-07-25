#  Contains code that is based in part on Chebfun v5's chebfun/padeappox.m,
# which is distributed with the following license:

# Copyright (c) 2015, The Chancellor, Masters and Scholars of the University
# of Oxford, and the Chebfun Developers. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Oxford nor the names of its
#       contributors may be used to endorse or promote products derived from
#       this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""
    robustpade(
        f::Function,
        m::Integer,
        n::Integer,
        x=0.;
        tol::Real=100eps()
    )

computes the (m,n) Pade approximant to a function f using TaylorSeries.taylor_expand to compute the Taylor expansion.

"""
function robustpade(f::Function,m::Integer,n::Integer,x=0.;kwargs...)
    taylorexpansion = taylor_expand(f, x; order=m+n+1)
    robustpade(taylorexpansion, m,n;kwargs...)
end

robustpade(p::Taylor1,args...;kwargs...) = robustpade(p.coeffs, args...;kwargs...)
robustpade(p::Polynomial,args...;kwargs...) = robustpade(p.coeffs, args...;kwargs...)

# Ensure float coefficients
function robustpade(coeffs::AbstractVector{<:Union{Integer,Complex{<:Integer}}},m::Integer, n::Integer;kwargs...)
    robustpade(float.(coeffs), m, n;kwargs...)
end

"""
    robustpade(
        coeffs::AbstractVector{T},
        m::Integer,
        n::Integer;
        tol::Real=eps(float(real(T)))
    )

computes the (m,n) Pade approximant to a function with Taylor coefficients 'coeffs' using SVD following Parchon et al. SIAM Review.

Adapted from the Chebfun implementation at https://github.com/chebfun/chebfun/blob/master/padeapprox.m (modified BSD license)

"""
function robustpade(coeffs::AbstractVector{T}, m::Integer, n::Integer; tol::Real=epsreal(T)) where {T<:RealOrComplexFloat}
    a,b = robustpade_coefficients(coeffs,m,n;tol)
    return Polynomial(a)//Polynomial(b)
end

function robustpade_coefficients(coeffs::AbstractVector{T}, m::Integer, n::Integer; tol::Real=epsreal(T)) where {T<:RealOrComplexFloat}

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
        return a,b
    else
        # Compute absolute tolerance.
        ts = tol * norm(col)
        # Form Toeplitz matrix.
        Z = Toeplitz(col, row)
        # Do diagonal hopping across block.
        m, n = robustpade_hop!(m, n, Z, ts)
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
        return a,b
    end
end

function robustpade_hop!(m, n, Z, ts)
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
    return robustpade_hop!(m - Δ, n - Δ, Z, ts)
end

function robustpade_size(f,m,n,args...;kwargs...)
    r = robustpade(f,m,n,args...;kwargs...)
    return degree(r.num),degree(r.den)
end
