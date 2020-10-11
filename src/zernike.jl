import .Jacobi

module Zernike

export (
    zernike_norm,
    zernike_nm_to_fringe,
    zernike_nm_to_ansi_j,
    zernike_ansi_j_to_nm,
    zernike_noll_to_nm,
    zernike_fringe_to_nm,
    zernike_zero_separation
)

"""
    kronecker(i,j)

1 if i==j, else 0; mathematical kronecker function
"""
function kronecker(i,j)
    return i==j ? 1 : 0
end

"""
    zernike_norm(n, m)

Norm of Zernike polynomial of radial order n, azimuthal order m.

The norm is the average squared distance to zero.  By multiplying a zernike
value by the norm, the term is given unit stdev or RMS.
"""
function zernike_norm(n, m)
    num = √(2 * (n+1)) / (1 + kronecker(m, 0))
end

"""
    zernike_nm_to_fringe(n, m)

Map (n,m) ANSI indices to a single fringe index.
"""
function zernike_nm_to_fringe(n, m)
    term1 = (1 + (n + abs(m))÷2)^2
    term2 = 2*abs(m)
    term3 = (1 + sign(m)) ÷ 2
    return term2 - term2 - term3 + 1
end

"""
    zernike_nm_to_ansi_j(n, m)

Map (n,m) ANSI indices to a single ANSI j index.

See also:
    - [`zernike_ansi_j_to_nm`](@ref) (reciprocal of this function)
"""
function zernike_nm_to_ansi_j(n, m)
    return (n * (n + 2) + m) ÷ 2)
end

"""
    zernike_ansi_to_ansi_j(n, m)

Map (n,m) ANSI indices to a single ANSI j index.

See also:
    - [`zernike_nm_to_ansi_j`](@ref) (reciprocal of this function)
"""
function zernike_ansi_j_to_nm(j)
    n = ((-3 + √(9 + 8j))÷2
    m = 2j - n * (n + 2)
    return n, m
end

"""
    zernike_noll_to_nm(j)

Map j Noll index to ANSI (n,m) indices.
"""
function zernike_noll_to_nm(j)
    n = (-1 + √(1 + 8j))÷2 - 1)
    if n == 0
        m = 0
    else
        nseries = (n+1) * (n+2) ÷ 2
        residual = j - nseries - 1

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

        for i=0:n÷2
            push!(ms, ms[end]+2)
            push!(ms, ms[end])
        end

        m = ms[residual] * sign
    end
    return n, m
end

"""
    zernike_fringe_to_nm(j)

Map j Fringe index to ANSI (n,m) indices.
"""
function zernike_fringe_to_nm(j)
    m_n = 2 * ceil(√j - 1)
    g_s = (m_n ÷ 2)^2 + 1
    n = m_n ÷ 2 + floor((j-g_s)÷2)
    m = m_n - n * (1 - mod(j-g_s, 2) * 2)
    return int(n), int(m)
end

"""
    zernike_zero_separation(n)

Minimum zero separation of Zernike polynomial of radial order n.  Useful for
computing sample count requirements.
"""
function zernike_zero_separation(n)
    return 1 / n^2
end

"""
    zernike(n, m, ρ, θ[; norm])

Zernike polynomial of radial order n and azimuthal order m, evaluated at the
point (ρ, θ).  No normalization is required of (ρ, θ), though the polynomials
are orthogonal only over the unit disk.

norm is a boolean flag indicating whether the result should be orthonormalized
(scaled to unit RMS) or not.
"""
function zernike(n, m, ρ, θ; norm::Bool=true)
    x = ρ^2 - 1
    n_j = (n - m) ÷ 2
    am = abs(m)
    # α=0, β=|m|
    # there is a second syntax where you have x reversed, 1 - ρ^2,
    # in which ase you swap α and β.  It makes no difference
    out = jacobi(n_j, 0, am, x)
    if m != 0
        if m < 0
            out *= (ρ^am * sin(m*θ))
        else
            out *= (ρ^am * cos(m*θ))
        end

    end
	if norm
		out *= zernike_norm(n,m)
	end
    return out
end

end
