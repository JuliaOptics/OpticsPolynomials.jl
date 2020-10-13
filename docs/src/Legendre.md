# Legendre

The Legendre polynomials are a family of orthogonal polynomials over the univariate domain $x \in [-1,1]$.  They are a special case of the Jacobi polynomials for $\alpha = \beta = 0$  The implementation in this library uses a recurrence relation for the Jacobi polynomials to compute them.  It is fast and stable to high order $n$.

## Usage

The Legendre usage follows from the jacobi one, in that the input is an array or a scalar with substantially higher performance when invoked on arrays.


## Core Functions

```@docs
OpticsPolynomials.legendre
```
