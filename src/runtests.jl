using Test

using OpticsPolynomials

@testset "jacobi correct values" begin
    # this only has to go as high as 3 to show that a three term recurrence
    # relation behaves correctly, but as a run for record we go a bit higher
    for n in 1:10
        # design of jacobi polynomials is such that edge = 1
        @test jacobi(n, 0, 0, 1) ≈ 1
    end
end
