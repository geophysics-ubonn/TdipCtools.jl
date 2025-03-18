using BenchmarkTools

function bench()
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 20))
    std = 0.001
    noise = Normal(0, std)
    noise_realization = rand(noise, length(timesteps))
    d = debyeResponseTimeDomain(timesteps, 1.0, test_relaxation_time) .+ noise_realization

    Cinv = Matrix{Float64}(I, length(d), length(d)) * 1e5
    tau_grid = createTauGrid(timesteps, 1.0, 25)
    G = createMatrixForwardOperator(timesteps, tau_grid)

    m0 = -6 * log.(ones(length(tau_grid)))

    (m, rmse, lambda) = occamDebyeDecomposition(m0, d, G, Cinv)
    return m, rmse, lambda
end

@info "Benchmarking occamDebyeDecomposition"
Random.seed!(100)
for i = 1:10
    @time bench()
end