Nsamples = 100

@testset "occamDebyeDecomposition medium error" begin
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 50))
    # create single term Debye transient
    R0 = 1.0
    chargeability = 0.1
    gamma = R0 * chargeability

    # noise
    rel = 0.01
    Random.seed!(100)
    for i = 1:Nsamples
        d = debyeResponseTimeDomain(timesteps, gamma, test_relaxation_time)
        std = d * rel .+ 0.001
        # synthetic data
        noise_realization = randn(length(timesteps)) .* std
        d = d .+ noise_realization

        # inversion
        Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std .^ -2)
        tau_grid = createTauGrid(timesteps, 1.5, 25)
        G = createMatrixForwardOperator(timesteps, tau_grid)
        m0 = -6 * ones(length(tau_grid))

        (m, rmse, lambda) = occamDebyeDecomposition(m0, d, G, Cinv)

        @test m != m0
        @test abs(rmse - 1.1) < 1e-1
        @test abs(log10(tau_grid[argmax(m)]) - log10(test_relaxation_time)) < 5e-1

        # plots
        # plot(tau_grid, exp.(m), xscale=:log10)
        # savefig("plots/rtd$i.png")
        # scatter(timesteps, d)
        # plot!(timesteps, decompositionResponseTimeDomain(m,
        #     tau_grid, timesteps, gamma))
        # savefig("plots/$i.png")
    end
end

@testset "occamDebyeDecomposition large error" begin
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 50))
    # create single term Debye transient
    R0 = 1.0
    chargeability = 0.1
    gamma = R0 * chargeability

    # noise
    rel = 0.05
    Random.seed!(100)
    for i = 1:Nsamples
        d = debyeResponseTimeDomain(timesteps, gamma, test_relaxation_time)
        std = d * rel .+ 0.001
        # synthetic data
        noise_realization = randn(length(timesteps)) .* std
        d = d .+ noise_realization

        # inversion
        Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std .^ -2)
        tau_grid = createTauGrid(timesteps, 1.5, 25)
        G = createMatrixForwardOperator(timesteps, tau_grid)
        m0 = -6 * ones(length(tau_grid))

        (m, rmse, lambda) = occamDebyeDecomposition(m0, d, G, Cinv)

        @test m != m0
        @test abs(rmse - 1.1) < 1e-1
        @test abs(log10(tau_grid[argmax(m)]) - log10(test_relaxation_time)) < 5e-1

        # plots
        # plot(tau_grid, exp.(m), xscale=:log10)
        # savefig("plots/rtd$i.png")
        # scatter(timesteps, d)
        # plot!(timesteps, decompositionResponseTimeDomain(m,
        #     tau_grid, timesteps, gamma))
        # savefig("plots/$i.png")
    end
end

@testset "test error" begin
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 50))
    # create single term Debye transient
    R0 = 1.0
    chargeability = 0.1
    gamma = R0 * chargeability

    # noise
    rel = 0.01
    Random.seed!(100)

    d = debyeResponseTimeDomain(timesteps, gamma, test_relaxation_time)
    std = d * rel .+ 0.001
    # synthetic data
    noise_realization = randn(length(timesteps)) .* std
    d = d .+ noise_realization

    # inversion
    Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std .^ -2)
    tau_grid = createTauGrid(timesteps, 1.5, 25)
    G = createMatrixForwardOperator(timesteps, tau_grid)
    m0 = -6 * ones(length(tau_grid))

    (m, rmse, lambda) = occamDebyeDecomposition(m0, d, G, Cinv)

    omegas = Float64[0.0]
    std_R0 = 10.0
    errorZ = estimateStandardDeviation(m, G, Cinv, lambda, omegas, tau_grid, std_R0)

    @test real(errorZ[1]) == std_R0
    @test imag(errorZ[1]) == 0.0
    println(errorZ)
end