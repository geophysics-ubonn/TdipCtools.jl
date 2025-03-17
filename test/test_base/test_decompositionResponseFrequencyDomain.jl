@testset "decompositionResponseFrequencyDomain" begin
    omega = Float64[0.0, 1.0]
    gamma = Float64[0.5, 0.5]
    tau_grid = Float64[1.0, 1.0]
    R0 = 1.0
    response = decompositionResponseFrequencyDomain(omega, gamma, tau_grid, R0)
    @test response[1] == 1.0 - 0.0im
    @test response[2] == 0.5 - 0.5im
end