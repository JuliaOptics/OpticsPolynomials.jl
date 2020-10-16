using LoopVectorization
using MuladdMacro

export jacobi, jacobi_weight
"""
    jacobi_weight(α, β, x)

Compute the weighting function (1-x)^α * (1-x)^β .
"""
function jacobi_weight(α, β, x)
	return (1-x)^α * (1-x)^β
end

"""
	calcac/b(n, α, β, x)

Two functions calcac and calcb compute the three coefficients in the recurrence
relation for the Jacobi polynomial.

It is broken into two functions to separate a and c, which do not depend on x
and can be calculated once for a vector of x, from b which does depend on x and
must be calculated for each element.

n is the order of the polynomial.

α,β are the two weighting parameters.

x is the argument or point at which the polynomial will be evaluated.
"""
function calcac(n, α, β)
    a = (2n) * (n + α + β) * (2n + α + β - 2)
    c = 2 * (n + α - 1) * (n + β - 1) * (2n + α + β)
    return a, c
end

@inline function calcb(n, α, β, x)
	# the parens emphasize the terms, not PEMDAS
    b1 = @muladd (2n + α + β - 1)
    b2 = @muladd (2n + α + β)
    b2 = b2 * (b2 - 2)
    b = @muladd b1 * (b2 * x + α^2 - β^2)
    return b
end

"""
    jacobi(n, α, β, x)

Compute the Jacobi polynomial of order n and weights α,β at point x.

This function uses a recurrence relation and is numerical stable to very high
order.  The computation time is linear w.r.t. the order n.  jacobi(n, a, b, x)
is also linear w.r.t. the size of argument x.  x should be passed as an array,
and jacobi should not be called with dot syntax, for best performance.
"""
function jacobi(n, α, β, x)
    if n == 0
        return ones(size(x))
    elseif n == 1
        return (α + 1) .+ (α + β + 2) .* ((x.-1)./2)
    end
	# Pnm2 = 1, skip because it is waste calculation
	# do scalar division on a to avoid vector division
	# the cod ebelow, therefore, visually does not
	# match the math, though it is correct
    Pnm1 = (α + 1) .+ (α + β + 2) .* ((x.-1)./2)
	a, c = calcac(2, α, β)
	inva = 1 / a
    Pn = ((calcb.(2, α, β, x) .* Pnm1) .- c) .* inva
    if n == 2
        return Pn
    end
    # Pnm2 = similar(Pnm1)
    for i = 3:n
		# Pnm2, Pnm1 = Pnm1, Pn
		Pnm2, Pnm1, Pn = Pnm1, Pn, Pnm2  # Pnm2 avoids a realloc
        a, c = calcac(n, α, β)
        inva = 1 / a
        @avx unroll=4 for i in eachindex(Pn)
            Pn[i] = (calcb(n, α, β, x[i]) * Pnm1[i] - c * Pnm2[i]) * inva
		end
        # Pn = ((calcb.(n, α, β, x) .* Pnm1) .- (c .* Pnm2)) .* inva
	end
    return Pn
end
