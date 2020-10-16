# Chebyshev

The Chebyshev polynomials are two families of orthogonal polynomials over the univariate domain ``x \in [-1,1]``.  They are special cases of the Jacobi polynomials for ``\alpha = \beta = -1/2`` for the Chebyshev polynomials of the first kind, and ``\alpha = \beta = 1/2`` for the Chebyshev polynomials of the second kind.  The implementation in this library uses a recurrence relation for the Jacobi polynomials to compute them.  It is fast and stable to high order $n$.

## Usage

The Chebyshev usage follows from the jacobi one, in that the input is an array or a scalar with substantially higher performance when invoked on arrays.


## Core Functions

```@docs
OpticsPolynomials.cheby1
OpticsPolynomials.cheby2
```
