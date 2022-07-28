var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = RobustPade","category":"page"},{"location":"#RobustPade","page":"Home","title":"RobustPade","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for RobustPade.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [RobustPade]","category":"page"},{"location":"#RobustPade.robustpade-Tuple{Polynomials.Polynomial, Vararg{Any}}","page":"Home","title":"RobustPade.robustpade","text":"robustpade(\n    f::Polynomial,\n    m::Integer,\n    n::Integer;\n    kwargs...\n)\n\ncomputes the (m,n) Pade approximant corresponding to the Taylor expansion given by a f::Polynomials.Polynomial\n\n\n\n\n\n","category":"method"},{"location":"#RobustPade.robustpade-Tuple{TaylorSeries.Taylor1, Vararg{Any}}","page":"Home","title":"RobustPade.robustpade","text":"robustpade(\n    f::Taylor1,\n    m::Integer,\n    n::Integer;\n    kwargs...\n)\n\ncomputes the (m,n) Pade approximant corresponding to the Taylor expansion given by a f::TaylorSeries.Taylor1\n\n\n\n\n\n","category":"method"},{"location":"#RobustPade.robustpade-Union{Tuple{T}, Tuple{AbstractVector{T}, Integer, Integer}} where T<:Union{AbstractFloat, Complex{<:AbstractFloat}}","page":"Home","title":"RobustPade.robustpade","text":"robustpade(\n    coeffs::AbstractVector{T},\n    m::Integer,\n    n::Integer;\n    tol::Real=eps(float(real(T)))\n) where T\n\ncomputes the (m,n) Pade approximant to a function with Taylor coefficients coeffs::AbstractVector using SVD following Parchon et al. SIAM Review.\n\nAdapted from the Chebfun implementation at https://github.com/chebfun/chebfun/blob/master/padeapprox.m (modified BSD license)\n\n\n\n\n\n","category":"method"},{"location":"#RobustPade.robustpade-Union{Tuple{T}, Tuple{Function, Integer, Integer}, Tuple{Function, Integer, Integer, T}} where T<:Number","page":"Home","title":"RobustPade.robustpade","text":"robustpade(\n    f::Function,\n    m::Integer,\n    n::Integer,\n    x::T=0.;\n    tol::Real=100eps(float(real(T)))\n) where T<:Number\n\ncomputes the (m,n) Pade approximant to a function f using TaylorSeries.taylor_expand to compute the Taylor expansion at x.\n\n\n\n\n\n","category":"method"},{"location":"#RobustPade.robustpade_size-Tuple{AbstractVector, Integer, Integer}","page":"Home","title":"RobustPade.robustpade_size","text":"robustpade_size(f::AbstractVector,m::Integer,n::Integer;kwargs...)\n\ncomputes the degrees of the (m,n) Pade approximant corresponding to the Taylor expansion coefficients f\n\n\n\n\n\n","category":"method"},{"location":"#RobustPade.robustpade_size-Tuple{Any, Integer, Integer, Vararg{Any}}","page":"Home","title":"RobustPade.robustpade_size","text":"robustpade_table(f,m::Integer,n::Integer,args...;kwargs...)\n\ncomputes the degrees of the (m,n) Pade approximant corresponding to f\n\n\n\n\n\n","category":"method"},{"location":"#RobustPade.robustpade_table-Tuple{Any, Integer, Integer, Vararg{Any}}","page":"Home","title":"RobustPade.robustpade_table","text":"robustpade_table(f,M::Integer,N::Integer,args...;kwargs...)\n\ncomputes the (0:M,0:N) Pade table corresponding to f\n\n\n\n\n\n","category":"method"}]
}
