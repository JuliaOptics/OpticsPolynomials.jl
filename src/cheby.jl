# this file depends on jacobi.jl from JuliaOptics/OpticsPolynomials that
# should be found in the same directory

export cheby1, cheby1_series, cheby1_sum, cheby2, cheby2_series, cheby2_sum

"""
    cheby1(n, x)

Compute the Chebyshev polynomial of the first kind of order n at point x.

This family of Chebyshev polynomials are a special case of the Jacobi polynomials
with α, β = -1/2.

See also: [`cheby1_series`](@ref), [`cheby1_sum`](@ref)
"""
function cheby1(n, x)
    return jacobi(n, -0.5, -0.5, x)
end

"""
    cheby1_series(ns, α, β, x)

Compute a series of Chebyshev polynomials of the first kind of orders n.
Returns an array with shape (size(x)..., length(ns)).
That is, the _final_ dimension contains the modes and the first dimension(s)
are spatial.

See also: [`cheby1`](@ref), [`cheby1_sum`](@ref)
"""
function cheby1_series(ns, x)
    return jacobi_series(ns, -0.5, -0.5, x)
end
"""
    cheby1_sum(ns, weights, x)

Compute a sum of Chebyshev polynomial of the first kind of order n weighted by weights.

See also: [`cheby1`](@ref), [`cheby1_series`](@ref)
"""
function cheby1_sum(ns, weights, x)
    return jacobi_sum(ns, weights, -0.5, -0.5, x)
end

"""
    cheby2(n, x)

Compute the Chebyshev polynomial of the second kind of order n at point x.

This family of Chebyshev polynomials are a special case of the Jacobi polynomials
with α, β = 1/2.

See also: [`cheby2_series`](@ref), [`cheby2_sum`](@ref)
"""
function cheby2(n, x)
    return jacobi(n, 0.5, 0.5, x)
end

"""
    cheby2_series(ns, α, β, x)

Compute a series of Chebyshev polynomials of the second kind of orders n.
Returns an array with shape (size(x)..., length(ns)).
That is, the _final_ dimension contains the modes and the first dimension(s)
are spatial.

See also: [`cheby1`](@ref), [`cheby1_sum`](@ref)
"""
function cheby2_series(ns, x)
    return jacobi_series(ns, 0.5, 0.5, x)
end
"""
    cheby1_sum(ns, weights, x)

Compute a sum of Chebyshev polynomial of the second kind of order n weighted by weights.

See also: [`cheby1`](@ref), [`cheby1_series`](@ref)
"""
function cheby2_sum(ns, weights, x)
    return jacobi_sum(ns, weights, 0.5, 0.5, x)
end
