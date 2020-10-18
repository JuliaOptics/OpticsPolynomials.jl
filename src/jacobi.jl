using LoopVectorization
using MuladdMacro

export jacobi, jacobi_weight
"""
    jacobi_weight(α, β, x)

Compute the weighting function (1-x)^α * (1-x)^β .
"""
function jacobi_weight(α, β, x)
	return (1 - x)^α * (1 - x)^β
end

# calcac_startb computes the a and c terms of the recurrence relation for the
# jacobi polynomials.  It also computes three pieces of b which do not depend on
# x.  This is not a formal docstring because documeter complains if I document
# the function and do not export it.

function calcac_startb(n, α, β)
    a = (2n) * (n + α + β) * (2n + α + β - 2)
    c = 2 * (n + α - 1) * (n + β - 1) * (2n + α + β)
    b1 = (2n + α + β - 1)
    b2 = (2n + α + β)
    b2 = b2 * (b2 - 2)
    b3 = α^2 - β^2
    return a, c, b1, b2, b3
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
        return (α + 1) .+ (α + β + 2) .* ((x .- 1) ./ 2)
    end
	Pnm1 = (α + 1) .+ (α + β + 2) .* ((x .- 1) ./ 2)
	a, c, b1, b2, b3 = calcac_startb(2, α, β)
	inva = 1 / a
	# .- c .* Pnm2, but Pnm2 == 1, so drop it
	Pn = similar(Pnm1)
    Pn .= @muladd (b1 .* (b2 .* x .+ b3) .* Pnm1 .- c) .* inva
    if n == 2
        return Pn
    end
	Pnm2 = similar(Pn)
	for i = 3:n
		# assign Pnm2 to Pn to reuse a buffer
		Pnm2, Pnm1, Pn = Pnm1, Pn, Pnm2
		a, c, b1, b2, b3 = calcac_startb(i, α, β)
        inva = 1 / a
		@avx unroll = 4 for i in eachindex(Pn)
			Pn[i] = (b1 * (b2 * x[i] + b3) * Pnm1[i] - c * Pnm2[i]) * inva
		end
	end
	return Pn
end

"""
    jacobi_series(ns, α, β, x)

Compute a series of jacobi polynomials of orders n.  Returns an array with shape
(length(ns), size(x)...).

See also: [`jacobi`](@ref), [`jacobi_sum`](@ref)
"""
function jacobi_series(ns, α, β, x)
    # allocate the output buffer
    out = Array{eltype(x), ndims(x)+1}(undef, length(ns), size(x)...)
    lowest_plane = 1
    if ns[1] == 0
        # set the first element of out to ones, if the first n the user wants
        # is n=0
        fill!(selectdim(out, 1, lowest_plane), 1.)
        lowest_plane += 1
    end
    # seed the recurrence relation
	Pnm1 = (α + 1) .+ (α + β + 2) .* ((x .- 1) ./ 2)
	a, c, b1, b2, b3 = calcac_startb(2, α, β)
	inva = 1 / a

    # unlike the regular jacobi function, we do not need to allocate
    # a scratch buffer, since we can use the lowest plane of the output
    # bufffer with impunity.  This only works because there is no weighting
    # so the lowest plane is guaranteed to contain a "pure" Jacobi polynomial
    # the exception is one buffer for Pnm2 to start up since we don't know if
    # the user opted to keep it

    # .- c .* Pnm2, but Pnm2 == 1, so drop it
    Pn = @muladd (b1 .* (b2 .* x .+ b3) .* Pnm1 .- c) .* inva
    Pnm2 = similar(Pn)

    if ns[lowest_plane] == 1
        view = selectdim(out, 1, lowest_plane)
        view .= Pnm1
        lowest_plane += 1
    end
	if ns[lowest_plane] == 2
        view = selectdim(out, 1, lowest_plane)
        view .= Pn
        lowest_plane += 1
    end

    # now do the loop that avoids recurrence
	for i = 3:n
        # Pnm2, Pnm1, Pn are purposfully not part of out
		Pnm2, Pnm1, Pn = Pnm1, Pn, Pnm2
		a, c, b1, b2, b3 = calcac_startb(i, α, β)
        inva = 1 / a
		@avx unroll = 4 for i in eachindex(Pn)
			Pn[i] = (b1 * (b2 * x[i] + b3) * Pnm1[i] - c * Pnm2[i]) * inva
        end
        for j in eachindex(ns)
			if ns[j] == i
                # this order is relevant to the sum
                view = selectdim(out, 1, lowest_plane)
                view .= Pn
                lowest_plane += 1
                break
			end
		end
	end
	return out
end

"""
    jacobi_sum(ns, weights, α, β, x)

Compute a sum of Jacobi polynomial of order n weighted by weights.

See [`jacobi`](@ref) for more details.  α, β are the spatial weights and x the
argument, passed as an array.
"""
function jacobi_sum(ns, weights, α, β, x)
	# start by allocating the output array, be a little clever
	# around the initialization value
	out = similar(x)
	if ns[1] == 0
		fill!(out, weights[1])
	else
		fill!(out, 0)
	end

	# Pnm1 = P_1, Pnm2= P_0, == 1
	# this block of setup is the same as before
	Pnm1 = (α + 1) .+ (α + β + 2) .* ((x.-1)./2)
	a, c, b1, b2, b3 = calcac_startb(2, α, β)
	inva = 1 / a
	# .- c .* Pnm2, but Pnm2 == 1, so drop it
    Pn = @muladd (b1 .* (b2 .* x .+ b3) .* Pnm1 .- c) .* inva

	# add P_2 if we can.  Two iffs is faster than
	# a searchsorted
	if ns[1] == 1
		out .+= weights[1] .* Pnm1
	elseif ns[2] == 1
		out .+= weights[2] .* Pnm1
	end
	if ns[1] == 2
		out .+= weights[1] .* Pn
	elseif ns[2] == 2
		out .+= weights[2] .* Pn
	end

	Pnm2 = similar(Pn);
	for i = 3:ns[end]
		# assign Pnm2 to Pn to reuse a buffer
		Pnm2, Pnm1, Pn = Pnm1, Pn, Pnm2
		a, c, b1, b2, b3 = calcac_startb(i, α, β)
        inva = 1 / a
		@avx unroll=4 for i in eachindex(Pn)
			Pn[i] = (b1 * (b2 * x[i] + b3) * Pnm1[i] - c * Pnm2[i]) * inva
		end
		# i is the order that was just computed
		for j in eachindex(ns)
			if ns[j] == i
				# this order is relevant to the sum
				out .+= weights[j] .* Pn
				break
			end
		end
	end
	return out
end
