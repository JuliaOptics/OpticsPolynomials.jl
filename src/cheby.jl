import .Jacobi

"""
    cheby1(n, x)

Compute the Chebyshev polynomial of the first kind of order n at point x.

This family of Chebyshev polynomials are a special case of the Jacobi polynomials
with α, β = -1/2.
"""
function cheby1(n, x)
    return Jacobi.jacobi(n, -0.5, -0.5, x)
end

"""
    cheby2(n, x)

Compute the Chebyshev polynomial of the second kind of order n at point x.

This family of Chebyshev polynomials are a special case of the Jacobi polynomials
with α, β = 1/2.
"""
function cheby2(n, x)
    return Jacobi.jacobi(n, 0.5, 0.5, x)
end
