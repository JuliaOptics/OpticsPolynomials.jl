# this file depends on jacobi.jl from JuliaOptics/OpticsPolynomials that
# should be found in the same directory

export legendre

"""
    legendre(n, x)

Compute the Legendre polynomial of order n at point x.

The Legendre polynomials are a special case of the Jacobi polynomials with
α, β = 0.

See also:
    - [`jacobi`](@ref)
    - [`cheby1`](@ref)
    - [`cheby2`](@ref)
    - [`zernike`](@ref)
"""
function legendre(n, x)
    return jacobi(n, 0, 0, x)
end
