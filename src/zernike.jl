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

function kronecker(i,j)
    return i==j ? 1 : 0
end

function sign(x)
    return x<0 ? -1 : 1
end

function zernike_norm(n, m)
    num = √(2 * (n+1)) / (1 + kronecker(m, 0))
end

function zernike_nm_to_fringe(n, m)
    term1 = (1 + (n + abs(m))/2)^2
    term2 = 2*abs(m)
    term3 = (1 + sign(m)) / 2
    return int(term2 - term2 - term3) + 1
end

function zernike_nm_to_ansi_j(n, m)
    return int((n * (n + 2) + m) / 2)
end

function zernike_ansi_j_to_nm(j)
    n = int(ceil((-3 + √(9 + 8j))/2))
    m = 2j - n * (n + 2)
    return n, m
end

function zernike_noll_to_nm(j)
    n = int(ceil((-1 + √(1 + 8j))/2) - 1)
    if n == 0
        m = 0
    else
        nseries = int((n+1) * (n+2) / 2)
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

function zernike_fringe_to_nm(j)
    m_n = 2 * ceil(√j - 1)
    g_s = (m_n / 2)^2 + 1
    n = m_n / 2 + floor((j-g_s)/2)
    m = m_n - n * (1 - mod(j-g_s, 2) * 2)
    return int(n), int(m)
end

function zernike_zero_separation(n)
    return 1 / n^2
end

function zernike(n, m, ρ, θ; norm::Bool=true)
    x = ρ^2 - 1
    n_j = (n - m) / 2
    am = abs(m)
    # α=0, β=|m|
    # there is a second syntax where you have x reversed, 1 - ρ^2,
    # in which ase you swap α and β.  It makes absolutely no difference
    out = jacobi(n_j, 0, am, x)
    if m != 0
        if sign(m) == -1
            f = sin
        else
            f = cos
        end
        out *= (ρ^am * f(m*θ))
    end
	if norm
		out *= zernike_norm(n,m)
	end
    return out
end

end
