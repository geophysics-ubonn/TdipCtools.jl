@testset "debyeManager" begin
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 20))
    # create single term Debye transient
    R0 = 1.0
    chargeability = 0.1
    gamma = R0 * chargeability

    std = 0.001
    Random.seed!(1000)

    noise = Normal(0, std)
    noise_realization = rand(noise, length(timesteps))
    d = debyeResponseTimeDomain(timesteps, gamma, test_relaxation_time)
    d = d .+ noise_realization
    Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std^-2)

    manager = initializeManager(d, timesteps, Cinv, R0=R0)
    @test typeof(manager) == Manager
    debyeManager!(manager, lambda=10.0)

    @test abs(argmin(abs.(test_relaxation_time .- manager.tau_grid)) - argmax(exp.(manager.m))) < 3

    setFrequenciesManager!(manager, Float64[1.0, 10.0])
    @test manager.frequencies == Float64[1.0, 10.0]

    spectrumManager(manager)
    # println(spectrum)
end

@testset "occamManager" begin
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 20))
    # create single term Debye transient
    R0 = 1.0
    chargeability = 0.1
    gamma = R0 * chargeability

    std = 0.001
    Random.seed!(1000)

    noise = Normal(0, std)
    noise_realization = rand(noise, length(timesteps))
    d = debyeResponseTimeDomain(timesteps, gamma, test_relaxation_time)
    d = d .+ noise_realization
    Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std^-2)

    manager = initializeManager(d, timesteps, Cinv, R0=R0)
    occamManager!(manager)
    @test abs(argmin(abs.(test_relaxation_time .- manager.tau_grid)) - argmax(exp.(manager.m))) < 3
end

@testset "debyeTomoManager" begin
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 20))
    # create single term Debye transient
    R0 = 1.0
    chargeability = 0.1
    gamma = R0 * chargeability
    std = 0.001

    Random.seed!(1000)

    tomo_d = []
    tomo_Cinv = []
    tomo_t = []
    tomo_r0 = Float64[]

    for i = 1:100
        noise = Normal(0, std)
        noise_realization = rand(noise, length(timesteps))
        d = debyeResponseTimeDomain(timesteps, gamma, test_relaxation_time)
        d = d .+ noise_realization
        Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std^-2)
        push!(tomo_d, d)
        push!(tomo_Cinv, Cinv)
        push!(tomo_t, timesteps)
        push!(tomo_r0, 1.0)
    end

    tomo_d = cat(tomo_d..., dims=2)
    tomo_Cinv = cat(tomo_Cinv..., dims=3)
    tomo_t = cat(tomo_t..., dims=2)

    tomo_manager = initializeTomoManager(tomo_d, tomo_t, tomo_Cinv, tomo_r0)
    @time debyeTomoManager!(tomo_manager, lambda=10.0)

    freqs = 10.0 .^ (range(-2, 2, 10))
    setFrequenciesTomoManager!(tomo_manager, freqs)
    spectra = spectrumTomoManager(tomo_manager)
    @test typeof(tomo_manager.managers[1].spectrum) == Vector{ComplexF64}

    # for manager = tomo_manager.managers
    #     @test abs(argmin(abs.(test_relaxation_time .- manager.tau_grid)) - argmax(exp.(manager.m))) < 3
    # end
end

@testset "occamTomoManager" begin
    test_relaxation_time = 0.5
    timesteps = 10.0 .^ (range(-1, 0, 50))
    # create single term Debye transient
    R0 = 1.0
    chargeability = 0.1
    gamma = R0 * chargeability
    rel = 0.05

    Random.seed!(1000)

    tomo_d = []
    tomo_Cinv = []
    tomo_t = []
    tomo_r0 = Float64[]

    for i = 1:100
        d = debyeResponseTimeDomain(timesteps, gamma, test_relaxation_time)
        std = d * rel .+ 0.001
        noise_realization = randn(length(timesteps)) .* std
        d = d .+ noise_realization
        Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std .^ -2)
        push!(tomo_d, d)
        push!(tomo_Cinv, Cinv)
        push!(tomo_t, timesteps)
        push!(tomo_r0, 1.0)
    end

    tomo_d = cat(tomo_d..., dims=2)
    tomo_Cinv = cat(tomo_Cinv..., dims=3)
    tomo_t = cat(tomo_t..., dims=2)

    tomo_manager = initializeTomoManager(tomo_d, tomo_t, tomo_Cinv, tomo_r0)
    @time occamTomoManager!(tomo_manager)

    # for manager = tomo_manager.managers
    #     # @test abs(argmin(abs.(test_relaxation_time .- manager.tau_grid)) - argmax(exp.(manager.m))) < 3
    #     @test abs(log10(manager.tau_grid[argmax(manager.m)]) - log10(test_relaxation_time)) < 1.0
    # end
end