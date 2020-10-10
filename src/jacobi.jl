
module Jacobi

export jacobi, jacobi_weight
"""
    jacobi_weight(α, β, x)

Compute the weighting function (1-x)^α * (1-x)^β .
"""
function jacobi_weight(α, β, x)
	return (1-x)^α * (1-x )^β
end

"""
    abc(n, α, β, x)

Compute the three coefficients in the recurrence relation for the Jacobi polynomial.

n is the order of the polynomial.

α,β are the two weighting parameters.

x is the argument or point at which the polynomial will be evaluated.
"""
function abc(n, α, β, x)
	# the parens emphasize the terms, not PEMDAS
	a = (2n) * (n + α + β) * (2n + α + β - 2)
	b1 = (2n + α + β -1)
	b2 = (2n + α + β)
	b2 = b2 * (b2 - 2)
	b = b1 * (b2 * x + α^2 - β^2)
	c = 2 * (n + α - 1) * (n + β - 1) * (2n + α + β)
	return a, b, c
end

"""
    jacobi(n, α, β, x)

Compute the Jacobi polynomial of order n and weights α,β at point x.

This function uses a recurrence relation and is numerical stable to very high
order.  The computation time is linear w.r.t. the order n.  jacobi.(n, a, b, x)
is also linear w.r.t. the size of argument x.
"""
function jacobi(n, α, β, x)
	if n == 0
		return 1
	elseif n == 1
		return (α + 1) + (α + β + 2) * ((x-1)/2)
	end
	# this uses a loop to avoid recursion
	# we init P of n-2, n-1, and n=2.
	# the latter is the first recursive term.
	# then loop from 3 up to n, updating
	# Pnm2, Pnm1, a, b, c, Pn
	Pnm2 = jacobi(0, α, β, x)
	Pnm1 = jacobi(1, α, β, x)
	a, b, c = abc(2, α, β, x)
	Pn = ((b * Pnm1) - (c * Pnm2)) / a
	if n == 2
		return Pn
	end
	for i = 3:n
		Pnm2 = Pnm1
		Pnm1 = Pn
		a, b, c = abc(n, α, β, x)
		Pn = ((b * Pnm1) - (c * Pnm2)) / a
	end
	return Pn
end

end
