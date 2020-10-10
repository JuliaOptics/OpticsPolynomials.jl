import .Jacobi

"""
    legendre(n, x)

Compute the Legendre polynomial of order n at point x.

The Legendre polynomials are a special case of the Jacobi polynomials with
α, β = 0.  The implementation here reflects that.
"""
function legendre(n, x)
    return Jacobi.jacobi(n, 0, 0, x)
end
