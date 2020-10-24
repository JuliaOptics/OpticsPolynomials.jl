
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
