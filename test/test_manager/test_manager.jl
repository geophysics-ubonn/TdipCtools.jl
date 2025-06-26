@testset "test_manager" begin
    # create data
    r0 = 1.0
    charg = 0.1
    gamma = r0 * charg
    tau = 0.5
    t = 10.0 .^ (range(-1, 1, 300))
    d = debyeResponseTimeDomain(t, gamma, tau)

    ### TEST INITIALIZATION ###
    # init with single standard deviation
    Cinv_Float = 0.001^-2
    manager = initializeManager(d, t, Cinv_Float, R0=r0)
    # init with Vector
    Cinv_Vector = Float64[fill(0.001^-2, length(d))...]
    mB = initializeManager(d, t, Cinv_Vector, R0=r0)
    # init with Matrix
    Cinv_Matrix = Matrix{Float64}(I, (length(d), length(d))) .* Cinv_Vector
    mC = initializeManager(d, t, Cinv_Matrix, R0=r0)
    # test
    @test (manager.d == mC.d) && (manager.d == mB.d)
    @test (manager.Cinv == mC.Cinv) && (manager.Cinv == mB.Cinv)

    ### TEST DEBYE DECOMPOSITION ###
    @test isnothing(manager.m)
    debyeManager!(manager, lambda=100.0)
    @test !isnothing(manager.m)
    @test !isnothing(manager.tau_grid)
    @test manager.G == createMatrixForwardOperator(t, manager.tau_grid)
    @test manager.lambda == 100.0
    @test abs(manager.tau_grid[argmax(manager.m)] - tau) < 0.5

    ### TEST OCCAM DECOMPOSITION ###
    occam_manager = initializeManager(d, t, 0.0001^(-2), R0=r0)
    occam_managerB = initializeManager(d, t, 0.001^(-2), R0=r0)
    occamManager!(occam_manager)
    occamManager!(occam_managerB)
    @test occam_manager.lambda < occam_managerB.lambda

    ### TEST SETTING FREQUENCIES ###
    frequencies = cat(0, (10.0 .^ (range(-1, 1, 10)))..., dims=1)
    setFrequenciesManager!(manager, frequencies)
    @test manager.frequencies == frequencies

    ### TEST SPECTRUM ###
    spectrum = spectrumManager(manager)
    @test spectrum == manager.spectrum
    @test manager.spectrum[1] == r0

    ### TEST ERROR ###
    errorZ = errorManager!(manager, 1.0)
    @info spectrum
    @info errorZ
end