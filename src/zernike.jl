# this file depends on jacobi.jl from JuliaOptics/OpticsPolynomials that
# should be found in the same directory

export zernike
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
    num = √(2 * (n+1)) / (1 + kronecker(m, 0))
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
    x = 2ρ.^2 .- 1
    n_j = (n - m) ÷ 2
    am = abs(m)
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
