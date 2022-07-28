# Helper functions
const RealOrComplexFloat = Union{AbstractFloat,Complex{<:AbstractFloat}}
epsreal(z) = eps(float(real(z)))

# Try to match MATLAB svd(A,0) described at https://uk.mathworks.com/help/matlab/ref/double.svd.html

function svd0!(A)
    m, n = size(A)
    m > n && return svd!(A; full=false)
    svd!(A; full=true)
end

issubarrayzero(x, idx, tol) = @views norm(x[idx], Inf) <= tol * norm(x, Inf)

# chop leading and trailing zeros

function findfirst_defaultlast(f,x)
    i = findfirst(f,x)
    isnothing(i) ? lastindex(x) : i
end
function findlast_defaultlast(f,x)
    i = findlast(f,x)
    isnothing(i) ? lastindex(x) : i
end
function firstnonzeroindex(x, tol)
    findfirst_defaultlast(z -> abs(z) > tol,x)
end
function lastnonzeroindex(x, tol)
    findlast_defaultlast(z -> abs(z) > tol,x)
end
