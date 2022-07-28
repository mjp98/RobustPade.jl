
#  Contains code that is based in part on Chebfun v5's chebfun/tests/misc/test_padeappox.m,
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

using RobustPade, Polynomials, TaylorSeries, LinearAlgebra
using Test

function polesandresidues(f)
    _,d = residues(f)
    p = [x for x in keys(d)]
    r = [first(d[x]) for x in p]
    return p, r
end

@testset "RobustPade.jl" begin
    # adapted from https://github.com/chebfun/chebfun/blob/master/tests/misc/test_padeapprox.m

    tol = 1e-10
    @testset "from function" begin
        m,n=10,10
        f = x -> (x^4 - 3) / ((x + 3.2) * (x - 2.2))
        r = robustpade(f, 10, 10)
        @test degree(r.num) == 4
        @test degree(r.den) == 2
        pol, res = polesandresidues(r)
        @test norm(sort(pol) - [-3.2; 2.2]) < tol
    end
    @testset "from coefficient vector" begin
        c = [1; 1im]
        a,b = robustpade(c, 0, 1)
        @test a[1] ≈ 1 atol = tol
        @test b[2] ≈ -1im atol = tol
    end
    @testset begin
        f = x -> x / (1 - x)
        r = robustpade(f, 5, 6)
        @test degree(r.num) == 1
        @test degree(r.den) == 1
        @test norm(r.num.coeffs - [0; 1]) < tol
        @test norm(r.den.coeffs - [1; -1]) < tol
        pol, res = polesandresidues(r)
        @test norm(pol - [1]) < tol
        @test norm(res - [-1]) < tol
    end
    @testset "padding" begin
        r = robustpade(Polynomial([1,2,3]), 5, 0)
        @test degree(r.num) == 2
        @test degree(r.den) == 0

    end
    @testset "polynomial" begin
        f = x -> 1 + x + x^4
        m,n = robustpade_size(f, 5, 0)
        @test m==4
        @test n==0
    end

    @testset "zero" begin
        a,b = robustpade(zeros(10), 5, 0)
        @test a == [0]
        @test b == [1]
        r = robustpade(zero, 5, 0)
        @test r(0) == 0
    end
end
