@testset "initializeManager" begin
    # create data
    r0 = 1.0
    charg = 0.1
    gamma = r0 * charg
    tau = 0.5
    t = 10.0 .^ (range(-1, 1, 30))
    d = debyeResponseTimeDomain(t, gamma, tau)

    # init with single standard deviation
    Cinv_Float = 0.1^-2
    mA = initializeManager(d, t, Cinv_Float, R0=r0)

    # init with Vector
    Cinv_Vector = Float64[fill(0.1^-2, length(d))...]
    mB = initializeManager(d, t, Cinv_Vector, R0=r0)

    # init with Matrix
    Cinv_Matrix = Matrix{Float64}(I, (length(d), length(d))) .* Cinv_Vector
    mC = initializeManager(d, t, Cinv_Matrix, R0=r0)

    # test
    @test (mA.d == mC.d) && (mA.d == mB.d)
    @test (mA.Cinv == mC.Cinv) && (mA.Cinv == mB.Cinv)
end

@testset "initializeTomoManager" begin
    # create data
    r0 = 1.0
    charg = 0.1
    gamma = r0 * charg
    tau = 0.5
    t = 10.0 .^ (range(-1, 1, 30))
    d = debyeResponseTimeDomain(t, gamma, tau)

    # controll single manager
    mC = initializeManager(d, t, 0.1^-2, R0=r0)

    # init tomoManager
    Cinv_Vector = Float64[fill(0.1^-2, length(d))...]

    dtomo = cat(d, d, d, dims=2)
    ttomo = cat(t, t, t, dims=2)
    r0tomo = cat(r0, r0, r0, dims=1)
    Ctomo = cat(Cinv_Vector, Cinv_Vector, Cinv_Vector, dims=2)
    Ctomo_Float = cat(0.1^-2, 0.1^-2, 0.1^-2, dims=1)

    tmanagerA = initializeTomoManager(dtomo, ttomo, Ctomo, r0tomo)
    tmanagerB = initializeTomoManager(dtomo, ttomo, Ctomo_Float, r0tomo)

    @test tmanagerA.managers[1].d == mC.d
    @test tmanagerA.managers[1].timesteps == mC.timesteps
    @test tmanagerA.managers[1].Cinv == mC.Cinv

    @test tmanagerB.managers[1].Cinv == mC.Cinv
end