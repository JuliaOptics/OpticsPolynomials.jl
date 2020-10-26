using GridCreation


@testset "zernike zero separation" begin
    for i in 1:100
        sep = zernike_zero_separation(i)
        @test 1/sep ≈ Float64(i)^2
    end
end

@testset "zernike fringe and nm round trip" begin
    for i in 1:100
        n, m = zernike_fringe_to_nm(i)
        j = zernike_nm_to_fringe(n,m)
        @test j == i
    end
end

@testset "zernike ANSI round trips" begin
    for i in 1:100
        n, m = zernike_ansi_j_to_nm(i)
        j = zernike_nm_to_ansi_j(n,m)
        @test j == i
    end
end


@testset "zernike gives correct values" begin
    for i in 1:100
        # a visual test proves the correctness of the zernike function in
        # general.  The values could only really be tested in a regression
        # fashion, generally, which is undesireable and assumes another
        # implementation is de-facto correct.

        # the test here, therefore, verifies only the analytic value of one
        # at the edge for the maxima of the azumithal term.
        n, m = zernike_noll_to_nm(i)
        out = zernike(n, m, [1.], [0.], norm=false)
        @test (out ≈ [1.]) || (out ≈ [0.])
        out = zernike(n, m, [1.], [0], norm=true)
        norm = [zernike_norm(n,m)]
        @test (out ≈ norm) || (out ≈ [0.])
    end
end

@testset "zernike series contains proper modes" begin
    x, y = mkCartVecs(1/8, 16)
    r, t = cartVecsToPolarGrid(x,y)
    modes = 1:11
    nms = zernike_noll_to_nm.(modes)
    out = zernike_series(nms, r, t)
    for mode in modes
        n, m = zernike_noll_to_nm(mode)
        truth = zernike(n, m, r, t)
        @test out[:, :, mode] ≈ truth
    end
end


@testset "zernike sum contains proper modes" begin
    x, y = mkCartVecs(1/8, 16)
    r, t = cartVecsToPolarGrid(x,y)
    nms = zernike_noll_to_nm.(1:11)
    weights = zeros(11)
    weights[5]=1.
    out = zernike_sum(nms, weights, r, t)
    truth = zernike(2, -2, r, t)
    @test out ≈ truth
end
