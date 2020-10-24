
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
        @test out ≈ [1.]
        out = zernike(n, m, [1.], [0.], norm=true)
        @test out ≈ [zernike_norm(n, m)]
    end
end
