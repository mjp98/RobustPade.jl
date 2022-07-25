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
        r = robustpade(c, 0, 1)
        a = r.num.coeffs
        b = r.den.coeffs
        @test a[1] ≈ 1 atol = tol
        @test b[2] ≈ -1im atol = tol
    end
    @testset "from Polynomial" begin
        c = Polynomial([1; 1im])
        r = robustpade(c, 0, 1)
        a = r.num.coeffs
        b = r.den.coeffs
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
        r = robustpade([1,2,3], 5, 0)
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
        a,b = robustpade_coefficients(zeros(10), 5, 0)
        @test a == [0]
        @test b == [1]
        r = robustpade(zero, 5, 0)
        @test r(0) == 0
    end
end
