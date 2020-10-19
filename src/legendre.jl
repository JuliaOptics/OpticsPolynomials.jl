# this file depends on jacobi.jl from JuliaOptics/OpticsPolynomials that
# should be found in the same directory

export legendre, legendre_series, legendre_sum

"""
    legendre(n, x)

Compute the Legendre polynomial of order n at point x.

The Legendre polynomials are a special case of the Jacobi polynomials with
α, β = 0.

See also: [`legendre_series`](@ref), [`legendre_sum`](@ref)
"""
function legendre(n, x)
    return jacobi(n, 0, 0, x)
end

"""
    legendre_series(ns, α, β, x)

Compute a series of Legendre polynomials of orders n.
Returns an array with shape (size(x)..., length(ns)).
That is, the _final_ dimension contains the modes and the first dimension(s)
are spatial.

See also: [`legendre`](@ref), [`legendre_sum`](@ref)
"""
function legendre_series(ns, x)
    return jacobi_series(ns, 0., 0., x)
end
"""
    legendre_sum(ns, weights, x)

Compute a sum of Legendre polynomial of order n weighted by weights.

See also: [`legendre`](@ref), [`legendre_series`](@ref)
"""
function legendre_sum(ns, weights, x)
    return jacobi_sum(ns, weights, 0., 0., x)
end
