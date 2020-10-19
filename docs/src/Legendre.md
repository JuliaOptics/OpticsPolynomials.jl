# Legendre

The Legendre polynomials are orthogonal over the univariate domain ``x \in [-1,1]``.  They are special cases of the Jacobi polynomials for ``\alpha = \beta = -0``.  The implementation in this library uses a recurrence relation for the Jacobi polynomials to compute them.  It is fast and stable to high order $n$.

## Usage

The Legendre usage follows from the jacobi one, in that the input is an array or a scalar with substantially higher performance when invoked on arrays.  See the [Jacobi](./Jacobi.md) documentation for more detailed information on this matter.

The key API is mirrored, substituting `jacobi` for `legendre`to produce `legendre`, `legendre_series`, and `legendre_sum`.

```@example simplest
using OpticsPolynomials
using Plots

x = collect(-1:0.01:1);
y = legendre(5,x);
plot(x,y)
png("legendre-5");
```
![](legendre-5.png)

```@example compare12
using OpticsPolynomials
using Plots

x = collect(-1:0.01:1);
y = legendre_series([1,2,3,4,5],x);
plot(x,y)
png("legendre-series");
```

![](legendre-series.png)

## Core Functions

```@docs
legendre
legendre_series
legendre_sum
```
