@testset "decompositionResponseTimeDomain" begin
    tau_grid = Float64[0.1, 1.0]
    timesteps = Float64[0.1, 1.0]
    m = Float64[0.0, 0.0]
    G = createMatrixForwardOperator(timesteps, tau_grid, 1.0)
    @test G * exp.(m) == decompositionResponseTimeDomain(m, tau_grid, timesteps, 1.0)
end