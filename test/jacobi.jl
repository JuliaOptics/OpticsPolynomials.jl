
@testset "jacobi correct values edge,origin" begin
    # this only has to go as high as 3 to show that a three term recurrence
    # relation behaves correctly, but as a run for record we go a bit higher
    for n in 1:10
        # design of jacobi polynomials is such that edge = 1
        @test jacobi(n, 0., 0., [1.]) ≈ [1.]
        correct = [-1.]
        if iseven(n)
            correct = [1.]
        end
        @test jacobi(n, 0., 0., [-1.]) ≈ correct
        if !iseven(n)
            @test jacobi(n, 0., 0., [0.]) ≈ [0.]
        end
    end
end

@testset "jacobi sum correct values" begin
    terms = [1,3,5,7,9]
    weights = ones(size(terms)...)
    @test jacobi_sum(terms, weights, 0., 0., [0.]) ≈ [0.]
    terms = [1,2,3,4,5]
    weights = ones(size(terms)...)
    # -1 accumulates 3 times and +1 two times on left edge => expect -1.
    # 1 accumulates 5 times on right edge => expect 5.
    @test jacobi_sum(terms, weights, 0., 0., [-1., 1.]) ≈ [-1., 5.]
end

@testset "jacobi series correct values contiguous" begin
    terms = [1,2,3,4,5]
    ledge = [-1.,1.,-1.,1.,-1.]
    center = [0.,-0.5,0.,0.375,0.]
    redge = [1.,1.,1.,1.,1.]
    series = jacobi_series(terms, 0., 0., [-1., 0., 1.])
    @test series[1,:] ≈ ledge
    @test series[2,:] ≈ center
    @test series[3,:] ≈ redge
end

@testset "jacobi series correct values sparse" begin
    terms = [1,3,5]
    ledge = [-1.,-1.,-1.]
    center = [0.,0.,0.]
    redge = [1.,1.,1.]
    series = jacobi_series(terms, 0., 0., [-1., 0., 1.])
    @test series[1,:] ≈ ledge
    @test series[2,:] ≈ center
    @test series[3,:] ≈ redge
end

