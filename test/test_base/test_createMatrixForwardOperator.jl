@testset "createMatrixForwardOperator" begin
    timesteps = 10.0 .^ (range(-1, 0, 2))
    tau_grid = 10.0 .^ (range(-1, 0, 2))
    G = createMatrixForwardOperator(timesteps, tau_grid)
    gamma = Float64[1.0, 1.0]
    transient = G * gamma
    reference_transient = (debyeResponseTimeDomain(timesteps, 1.0, 0.1)
                           .+
                           debyeResponseTimeDomain(timesteps, 1.0, 1.0))
    difference = reference_transient - transient
    @test difference == Float64[0.0, 0.0]
 
    timesteps = 10.0 .^ (range(-1, 0, 2))
    tau_grid = 10.0 .^ (range(-1, 0, 2))
    G = createMatrixForwardOperator(timesteps, tau_grid)
    gamma = Float64[1.0, 0.1]
    transient = G * gamma
    reference_transient = (debyeResponseTimeDomain(timesteps, 1.0, 0.1)
                           .+
                           debyeResponseTimeDomain(timesteps, 0.1, 1.0))
    difference = reference_transient - transient
    @test difference == Float64[0.0, 0.0]
 end