using BenchmarkTools

@testset "debyeDecomposition" begin
   test_relaxation_time = 0.5
   timesteps = 10.0 .^ (range(-1, 0, 20))
   d = debyeResponseTimeDomain(timesteps, 1.0, test_relaxation_time)
   Cinv = Matrix{Float64}(I, length(d), length(d)) * 1e5
   tau_grid = createTauGrid(timesteps, 1.0, 25)
   G = createMatrixForwardOperator(timesteps, tau_grid, 1.5)

   m0 = -6 * ones(length(tau_grid))
   @time m, rmse, lambda = debyeDecomposition(m0, d, G, Cinv, max_iter=1000, lambda=10.0)
   @test abs(tau_grid[argmax(m)] - test_relaxation_time) < 1e-1
   @test m != m0
end