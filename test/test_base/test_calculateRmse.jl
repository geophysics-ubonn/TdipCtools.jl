@testset "calculateRmse" begin
    f = ones(Float64, 10)
    d = ones(Float64, 10)
    Cinv = Matrix{Float64}(I, length(d), length(d))
    @test calculateRmse(f, d, Cinv) == 0.0
 
    f = ones(Float64, 10)
    d = zeros(Float64, 10)
    Cinv = Matrix{Float64}(I, length(d), length(d))
    @test calculateRmse(f, d, Cinv) == 1.0
 
    f = ones(Float64, 10) .+ 1
    d = zeros(Float64, 10)
    Cinv = Matrix{Float64}(I, length(d), length(d))
    @test calculateRmse(f, d, Cinv) > 1.0
 
    f = ones(Float64, 10) .- 0.5
    d = zeros(Float64, 10)
    Cinv = Matrix{Float64}(I, length(d), length(d))
    @test calculateRmse(f, d, Cinv) < 1.0
 end