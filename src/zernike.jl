# this file depends on jacobi.jl from JuliaOptics/OpticsPolynomials that
# should be found in the same directory

export zernike
export zernike_series
export zernike_norm
export zernike_nm_to_fringe
export zernike_nm_to_ansi_j
export zernike_ansi_j_to_nm
export zernike_noll_to_nm
export zernike_fringe_to_nm
export zernike_zero_separation


function kronecker(i,j)
    return i==j ? 1 : 0
end

"""
    zernike_norm(n, m)

Norm of Zernike polynomial of radial order n, azimuthal order m.

The norm is the average squared distance to zero.  By multiplying a zernike
value by the norm, the term is given unit stdev or RMS.
"""
function zernike_norm(n::Integer, m::Integer)
    return √(2 * (n+1)) / (1 + kronecker(m, 0))
end

"""
    zernike_nm_to_fringe(n, m)

Map (n,m) ANSI indices to a single fringe index.
"""
function zernike_nm_to_fringe(n::Integer, m::Integer)
    term1 = (1 + (n + abs(m))/2)^2
    term2 = 2 * abs(m)
    term3 = (1 + sign(m)) / 2
    return Int(trunc(term1 - term2 - term3)) + 1
end

"""
    zernike_nm_to_ansi_j(n, m)

Map (n,m) ANSI indices to a single ANSI j index.

See also:
    - [`zernike_ansi_j_to_nm`](@ref) (reciprocal of this function)
"""
function zernike_nm_to_ansi_j(n::Integer, m::Integer)
    return (n * (n + 2) + m) ÷ 2
end

"""
    zernike_ansi_to_ansi_j(n, m)

Map (n,m) ANSI indices to a single ANSI j index.

See also:
    - [`zernike_nm_to_ansi_j`](@ref) (reciprocal of this function)
"""
function zernike_ansi_j_to_nm(j::Integer)
    n = (-3 + √(9 + 8j))÷2
    m = 2j - n * (n + 2)
    return Int(n), Int(m)
end

"""
    zernike_noll_to_nm(j)

Map j Noll index to ANSI (n,m) indices.
"""
function zernike_noll_to_nm(j::Integer)
    n = ceil(Int, (-1 + √(1 + 8j))/2 - 1)
    if n == 0
        m = 0
    else
        nseries = (n+1) * (n+2) ÷ 2
        residual = j - nseries

        if isodd(j)
            sign = -1
        else
            sign = 1
        end

        if isodd(n)
            ms = [1,1]
        else
            ms = [0]
        end

        for i=0:(n÷2)-1
            push!(ms, ms[end]+2)
            push!(ms, ms[end])
        end

        m = ms[end+residual] * sign
    end
    return n, Int(m)
end

"""
    zernike_fringe_to_nm(j)

Map j Fringe index to ANSI (n,m) indices.
"""
function zernike_fringe_to_nm(j::Integer)
    m_n = 2 * (ceil(Int, √j) - 1)
    g_s = (m_n ÷ 2)^2 + 1
    n = m_n ÷ 2 + (j-g_s)÷2
    m = (m_n - n) * (1 - mod(j-g_s, 2) * 2)
    return n, m
end

"""
    zernike_zero_separation(n)

Minimum zero separation of Zernike polynomial of radial order n.  Useful for
computing sample count requirements.
"""
function zernike_zero_separation(n::Integer)
    return 1 / n^2
end

"""
    zernike(n, m, ρ, θ[; norm])

Zernike polynomial of radial order n and azimuthal order m, evaluated at the
point (ρ, θ).  No normalization is required of (ρ, θ), though the polynomials
are orthogonal only over the unit disk.

norm is a boolean flag indicating whether the result should be orthonormalized
(scaled to unit RMS) or not.

The zernike polynomials' radial basis is a special case of the Jacobi
polynomials under the transformation n_jacobi = (n-m)/2, α=0, β=|m|, x=2ρ^2-1.
"""
function zernike(n::Integer, m::Integer, ρ, θ; norm::Bool=true)
    x = 2 .* ρ.^2 .- 1
    am = abs(m)
    n_j = (n - am) ÷ 2
    out = jacobi(n_j, 0, am, x)
    if m != 0
        if m < 0
            out .*= (ρ.^am .* sin.(m.*θ))
        else
            out .*= (ρ.^am .* cos.(m.*θ))
        end

    end
	if norm
		out .*= zernike_norm(n,m)
	end
    return out
end


function zernike_series(nms, ρ, θ, norm::Bool=true)
    # 1.  Calculate a lookup table of Jacobi polynomials for the unique elements combinations of $|m|$ and $n_j$.
    # 2.  Calculate a lookup table of $\rho^{|m|}$ for each unique $|m|$
    # 3.  Calculate a lookup table of $\sin{m\theta}$ and $\cos{m\theta}$ for each unique $m$.
    # 4.  For each input $n, m$ look up the appropriate elements and compute the Zernike polynomial.
	x = 2 .* ρ.^2 .- 1
    n_j_list = [((n - abs(m)) ÷ 2, abs(m)) for (n,m) in nms]
    abs_m_type = typeof(n_j_list[1][2])
    dtype = eltype(ρ)
    jacobi_sequences_mnj = Dict{abs_m_type,Array{typeof(n_j_list[1][1]), 1}}()
	for (nj, am) in n_j_list
		a = get!(jacobi_sequences_mnj, am) do
			Array{abs_m_type, 1}()
		end
		push!(a, nj)
	end
	for key in keys(jacobi_sequences_mnj)
		# sort may not be needed.  Inherently sorted
		# but not sure if unique guarantees not shuffling
		# input
		sort!(unique!(jacobi_sequences_mnj[key]))
	end
    jacobi_sequences_mnj = sort(jacobi_sequences_mnj)
    # jacobi sequences_mnj now maps |m| => max n_j
    # this is step 1 at the top of the function
    jacobi_sequences = Dict{abs_m_type,Array{dtype,ndims(ρ)+1}}()
    for k in keys(jacobi_sequences_mnj)
        njs = jacobi_sequences_mnj[k]
        jacobi_sequences[k] = jacobi_series(njs, 0., dtype(k), x)
    end

    # jacobi_sequences now maps |m| => arrays with n_j along the last dimension
    # powers_of_rho is now ρ^|m| for each unique |m|
    plane_idx = ndims(ρ)+1
    sines = Dict{abs_m_type, typeof(ρ)}()
    cosines = Dict{abs_m_type, typeof(ρ)}()
    powers_of_rho = Dict{abs_m_type, typeof(ρ)}()
    out = Array{eltype(ρ), ndims(ρ)+1}(undef, size(ρ)..., length(nms))
    for i in eachindex(nms)
        n, m = nms[i]
        absm = abs(m)

        out_view = selectdim(out, plane_idx, i)
		n_j = (n-abs(m))÷2
		jac_plane = first(searchsorted(jacobi_sequences_mnj[absm], n_j))
        jacobi_piece = selectdim(jacobi_sequences[absm], plane_idx, jac_plane)
        if m == 0
            if norm
                norm_val = zernike_norm(n, m)
                out_view .= norm_val .* jacobi_piece
            else
                out_view .= jacobi_piece
            end
        else
            if m < 0
				# these evaluate the sin/cosine only as needed,
				# the do block is only called if get! would error
                azpiece = get!(sines, m) do
					sin.(m .* θ)
				end
            else
                azpiece = get!(cosines, m) do
					cos.(m .* θ)
				end
            end
            radialpiece = get!(powers_of_rho, absm) do
				ρ .^ absm
			end

            if norm
                norm_val = zernike_norm(n, m)
                out_view .= norm_val .* jacobi_piece .* radialpiece .* azpiece
            else
                out_view .= jacobi_piece .* radialpiece .* azpiece
            end
        end
    end
    return out
end

