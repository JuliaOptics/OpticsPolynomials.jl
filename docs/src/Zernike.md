# Zernike

The Zernike polynomials are a family of orthogonal polynomials over the unit disk, $\rho \in [0,1]$, $\theta \in [0, 2\pi]$.  They were originally defined in "Beugungstheorie des schneidenver-fahrens und seiner verbesserten form,der phasenkontrastmethode", F. Zernike (1934).  A more typical citation is "Zernike polynomials and atmospheric turbulence," Noll (1976).

The implementation in this library utilizes the property that, given a Jacobi polynomial of order ``n`` with weights ``\alpha,\beta \equiv P_{n}^{(\alpha,\beta)}(x)``, the radial component of the Zernike polynomials for argument ``\rho`` is ``P_{(n-m)/2}^{(0,|m|)}(\2rho^2 -1)``.  This is combined with a simple ``\sin(m\theta)`` or cosine equivalent, depending on the sign of ``m``.  The implementation here reflects that.

## Usage

The Zernike usage follows from the jacobi one, in that the input is an array or a scalar with substantially higher performance when invoked on arrays.

```@example simple
using GridCreation
x, y = mkCartVecs(1/8, 16); # -1,1
r, t = cartVecsToPolarGrid(x,y);



## Core Functions
```@docs
OpticsPolynomials.zernike
OpticsPolynomials.zernike_norm
```

## Index Conversions

This package includes a suite of functions for converting between various indexing schemes.  All functions work for arbitrary order.

### Noll

```@docs
OpticsPolynomials.zernike_noll_to_nm
OpticsPolynomials.zernike_nm_to_noll
```

### Fringe

```@docs
OpticsPolynomials.zernike_fringe_to_nm
OpticsPolynomials.zernike_nm_to_fringe
```

### ANSI Z80.28:2004

```@docs
OpticsPolynomials.zernike_ansi_j_to_nm
OpticsPolynomials.zernike_nm_to_ansi_j
```

## Utilities

```@docs
OpticsPolynomials.zernike_zero_separation
```
