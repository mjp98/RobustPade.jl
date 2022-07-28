# RobustPade.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjp98.github.io/RobustPade.jl/dev)
[![Build Status](https://github.com/mjp98/RobustPade.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mjp98/RobustPade.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mjp98/RobustPade.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mjp98/RobustPade.jl)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Julia implementation of the Pade approximation algorithm of P. Gonnet, S. Guettel, and L. N. Trefethen, "Robust Pade approximation via SVD", SIAM Rev., 55:101-117, 2013.

## Installation

This package is unregistered, so may be added by 

```julia
(v1.7) pkg> add https://github.com/mjp98/RobustPade.jl
```

This package supports Julia v1.7 and later.

## Usage

```julia
julia> using RobustPade
```
### Example

Construct a Pade approximant of ``\exp`` of type ``(2,2)``, returning a `Polynomials.RationalFunction`

```julia
julia> r = robustpade(exp,2, 2)
(1.0 + 0.4999999999999998*x + 0.08333333333333322*x^2) // (1.0 - 0.5000000000000002*x + 0.08333333333333347*x^2)
```
When the first argument is a function, `TaylorSeries.taylor_expand` is called to generate the taylor expansion from which the Pade approximation is computed.

## Similar packages
Pade approximation algorithms are available in:

- [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl) method is unstable - see this [issue](https://github.com/JuliaMath/Polynomials.jl/issues/161)
- SymPy's methods may be called, as described in this [comment](https://github.com/JuliaMath/Polynomials.jl/issues/161#issuecomment-456744016). The packages [Wynn.jl](https://github.com/J-Revell/Wynn.jl) and [Pade.jl](https://github.com/J-Revell/Pade.jl) provide some methods for calculating epsilon tables using SymPy.

For more general rational approximation algorithms, see

 - [Remez.jl](https://github.com/simonbyrne/Remez.jl)
 - [BaryRational.jl](https://github.com/macd/BaryRational.jl)
