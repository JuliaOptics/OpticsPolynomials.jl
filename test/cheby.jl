
@testset "cheby1 correct values origin" begin
    # legendre == jacobi(a,b=0,0); same function body
    for n in 1:10
        if !iseven(n)
            @test cheby1(n, [0.]) ≈ [0.]
        end
    end
end

@testset "cheby2 correct values origin" begin
    # legendre == jacobi(a,b=0,0); same function body
    for n in 1:10
        if !iseven(n)
            @test cheby2(n, [0.]) ≈ [0.]
        end
    end
end

@testset "cheby1 sum correct values origin" begin
    terms = [1,3,5,7,9]
    weights = ones(size(terms)...)
    @test cheby1_sum(terms, weights, [0.]) ≈ [0.]
end

@testset "cheby2 sum correct values origin" begin
    terms = [1,3,5,7,9]
    weights = ones(size(terms)...)
    @test cheby2_sum(terms, weights, [0.]) ≈ [0.]
end

@testset "cheby1 series correct value origin" begin
    terms = [1,3,5]
    out = [0., 0., 0.,]
    series = cheby1_series(terms, [0.,])
    @test series ≈ out'
end

@testset "cheby2 series correct value origin" begin
    terms = [1,3,5]
    out = [0., 0., 0.,]
    series = cheby2_series(terms, [0.,])
    @test series ≈ out'
end

