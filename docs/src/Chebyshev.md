# Chebyshev

The Chebyshev polynomials are two families of orthogonal polynomials over the univariate domain ``x \in [-1,1]``.  They are special cases of the Jacobi polynomials for ``\alpha = \beta = -1/2`` for the Chebyshev polynomials of the first kind, and ``\alpha = \beta = 1/2`` for the Chebyshev polynomials of the second kind.  The implementation in this library uses a recurrence relation for the Jacobi polynomials to compute them.  It is fast and stable to high order $n$.

## Usage

The Chebyshev usage follows from the jacobi one, in that the input is an array or a scalar with substantially higher performance when invoked on arrays.  See the [Jacobi](./Jacobi.md) documentation for more detailed information on this matter.

The key API is mirrored, substituting `jacobi` for `cheby1` or `cheby2` to produce `cheby1`, `cheby1_series`, and `cheby1_sum` and associated `cheby2` functions.

```@example simplest
using OpticsPolynomials
using Plots

x = collect(-1:0.01:1);
y = cheby1(5,x);
plot(x,y)
png("cheby1-5");
```
![](cheby1-5.png)

The visual difference between the polynomials of the first and second kind can be seen:

```@example compare12
using OpticsPolynomials
using Plots

x = collect(-1:0.01:1);
y = cheby1_series([1,2,3,4,5],x);
plot(x,y)
png("cheby1-series");

y = cheby2_series([1,2,3,4,5],x);
plot(x,y)
png("cheby2-series");
```

### First Kind

![](cheby1-series.png)

### Second Kind

![](cheby2-series.png)

## Core Functions

```@docs
OpticsPolynomials.cheby1
OpticsPolynomials.cheby1_series
OpticsPolynomials.cheby1_sum
OpticsPolynomials.cheby2
OpticsPolynomials.cheby2_series
OpticsPolynomials.cheby2_sum
```
