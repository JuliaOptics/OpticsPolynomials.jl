# Jacobi

The Jacobi polynomials are a family of orthogonal polynomials over the univariate domain ``x \in [-1,1]`` with respect to weighting parameters ``\alpha, \beta``.  The implementation in this library uses a recurrence relation to compute them.  It is fast and stable to high order $n$.  The recurrence relation is:

```math
a \cdot P_n^{(\alpha,\beta)} = b \cdot x \cdot P_{n-1}^{(\alpha,\beta)} - c \cdot P_{n-2}^{(\alpha,\beta)}
```

for some scalar ``a,b,c``.  Consult "Handbook of Mathematical Functions," Abramowitz and Stegun or Wikipedia for more information.

## Usage

The input is an array or a scalar.  Performance is substantially higher when invoked on arrays.


## Core Functions

```@docs
OpticsPolynomials.jacobi
```

## Utilities

```@docs
OpticsPolynomials.jacobi_weight
```
