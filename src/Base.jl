function calculateRmse(
    f::Vector{Float64},
    d::Vector{Float64},
    Cinv::Union{SparseMatrixCSC,Matrix}
)::Float64

    return sqrt(((f - d)' * Cinv * (f - d)) / length(f))
end

function createMatrixForwardOperator(
    timesteps::Vector{Float64},
    tau_grid::Vector{Float64},
    gamma::Float64
)::Matrix{Float64}

    G = zeros(Float64, (length(timesteps), length(tau_grid)))
    for (i, t) in enumerate(timesteps)
        G[i, :] = exp.(-t ./ tau_grid)
    end
    return G # * gamma
end

function createTauGrid(
    timesteps::Vector{Float64},
    extend::Float64,
    points_per_decade::Int64
)::Vector{Float64}

    tau_min = log10(minimum(timesteps)) - extend
    tau_max = log10(maximum(timesteps)) + extend
    ntau = Int64((tau_max - tau_min) * points_per_decade)
    return 10.0 .^ (range(tau_min, tau_max, ntau))
end

function decompositionResponseTimeDomain(
    m::Vector{Float64},
    tau_grid::Vector{Float64},
    timesteps::Vector{Float64},
    gamma::Float64
)::Vector{Float64}

    G = createMatrixForwardOperator(timesteps, tau_grid, gamma)
    return G * exp.(m)
end

function decompositionResponseFrequencyDomain(
    omega::Vector{Float64},
    gamma::Vector{Float64},
    tau_grid::Vector{Float64},
    R0::Float64
)::Vector{ComplexF64}

    Z = zeros(ComplexF64, length(omega))
    for (iw, w) in enumerate(omega)
        Polarization::ComplexF64 = 0 + 0im
        for i = 1:length(gamma)
            Polarization += gamma[i] * (1 - (1 + 1im * w * tau_grid[i])^(-1))
        end
        Z[iw] = R0 - Polarization
    end
    return Z
end

function debyeResponseTimeDomain(
    timesteps::Vector{Float64},
    gamma::Float64,
    tau::Float64
)::Vector{Float64}

    return gamma * exp.(-timesteps / tau)
end

function debyeDecomposition(m0::Vector{Float64},
    d::Vector{Float64},
    G::Matrix{Float64},
    Cinv::Union{SparseMatrixCSC,Matrix};
    lambda::Float64=1e6,
    eta::Float64=1.0,
    tune_eta::Bool=true,
    max_iter::Int64=1000
)::Tuple{Vector{Float64},Float64,Float64}

    m = copy(m0)
    Cinv = sparse(Cinv)
    R = initializeSmoothingOperator(length(m))

    # pre-acclocate memory for arrays
    Hessian = Matrix{Float64}(undef, (length(m), length(m)))
    J = Matrix{Float64}(undef, (length(d), length(m)))
    gradient = Vector{Float64}(undef, length(m))
    f = Vector{Float64}(undef, length(d))
    m_update = Vector{Float64}(undef, length(m))
    trial_m = Vector{Float64}(undef, length(m))

    rmse = inversion!(m, m_update, trial_m, d, f, G, Cinv,
        J, Hessian, gradient, lambda, eta,
        tune_eta, max_iter, R)
    return m, rmse, lambda
end

function evaluateJacobian!(
    J::Matrix{Float64},
    G::Matrix{Float64},
    m::Vector{Float64}
)::Nothing

    @assert size(J) == size(G)
    @assert size(J)[2] == length(m)

    # evaluate jacobian matrix for current model
    for i = 1:size(J)[1]
        J[i, :] .= view(G, i, :) .* exp.(m)
    end
    return nothing
end

function initializeSmoothingOperator(
    M::Int64
)::SparseMatrixCSC

    # initialize smoothing operator
    R = spzeros(Float64, (M - 1, M))
    for i in 1:M-1
        R[i, i] = 1.0
        R[i, i+1] = -1.0
    end
    return R' * R
end

function inversion!(
    m::Vector{Float64},
    m_update::Vector{Float64},
    trial_m::Vector{Float64},
    d::Vector{Float64},
    f::Vector{Float64},
    G::Matrix{Float64},
    Cinv::Union{SparseMatrixCSC,Matrix},
    J::Matrix{Float64},
    Hessian::Matrix{Float64},
    gradient::Vector{Float64},
    lambda::Float64,
    eta::Float64,
    tune_eta::Bool,
    max_iter::Int64,
    R::SparseMatrixCSC
)::Float64

    # initialize regularization tuning
    f .= G * exp.(m)
    rmse = calculateRmse(d, f, Cinv)

    for iteration = 1:max_iter
        # evaluate forward response
        f .= G * exp.(m)
        rmse = calculateRmse(f, d, Cinv)

        # perform single gauss newton update
        gaussNewtonUpdate!(m, m_update, d, f, G, Cinv, J, Hessian, gradient, lambda, R)

        # finish and return inversion result of MAP has been found
        if (norm(m_update) < 1e-4)
            @debug "Inversion finished after $iteration iterations."
            return rmse
        end

        # line search for optimal step length
        if tune_eta
            eta = lineSearchStepLength!(m, m_update, trial_m, d, f, Cinv, eta, G, lambda * R)
        end

        # m = m .+ eta * m_update
        axpy!(eta, m_update, m) # faster for large models. ovewrites m!
    end
    @debug "Debye decomposition not converged after max_iter=$max_iter iterations."
    return rmse
end

function lineSearchStepLength!(
    m::Vector{Float64},
    m_update::Vector{Float64},
    trial_m::Vector{Float64},
    d::Vector{Float64},
    f::Vector{Float64},
    Cinv::Union{SparseMatrixCSC,Matrix},
    eta::Float64,
    G::Matrix{Float64},
    R::SparseMatrixCSC
)::Float64

    # calculate rmse for trial step lengths
    trial_etas = Float64[0.0, 0.5, 1.0] * eta
    trial_rmse = Float64[0.0, 0.0, 0.0]

    for i = 1:3
        trial_eta = trial_etas[i]
        trial_m .= m + trial_eta * m_update
        f .= G * exp.(trial_m)
        rmse = (d - f)' * Cinv * (d - f) + trial_m' * R * trial_m
        trial_rmse[i] = rmse
    end

    # fit a parabola to the trial rmse vs trial eta and find the minimum
    A = ones(Float64, (3, 3))
    A[:, 1] = trial_etas .^ 2
    A[:, 2] = trial_etas
    a, b, _ = A \ trial_rmse
    eta_opt = -b / (2 * a)

    # keep eta between 0 and 1 for stability
    if eta_opt < 0.0 || eta_opt > 1.0
        return eta
    end
    return eta_opt
end

function gaussNewtonUpdate!(
    m::Vector{Float64},
    m_update::Vector{Float64},
    d::Vector{Float64},
    f::Vector{Float64},
    G::Matrix{Float64},
    Cinv::Union{SparseMatrixCSC,Matrix},
    J::Matrix{Float64},
    Hessian::Matrix{Float64},
    gradient::Vector{Float64},
    lambda::Float64,
    R::SparseMatrixCSC
)::Nothing

    # evaluate Gauss Newton model update
    evaluateJacobian!(J, G, m)
    Hessian .= J' * Cinv * J + lambda * R
    gradient .= J' * Cinv * (d - f) - lambda * R * m
    m_update .= Hessian \ gradient
    return nothing
end

function occamDebyeDecomposition(
    m0::Vector{Float64},
    d::Vector{Float64},
    G::Matrix{Float64},
    Cinv::Union{SparseMatrixCSC,Matrix};
    max_iter::Int64=1000
)::Tuple{Vector{Float64},Float64,Float64}

    m = copy(m0)
    Cinv = sparse(Cinv)
    f = G * exp.(m)
    rmse = calculateRmse(d, f, Cinv)
    rmse0 = rmse
    target_rmse = 1.1
    R = initializeSmoothingOperator(length(m))

    # pre-acclocate memory for arrays
    Hessian = Matrix{Float64}(undef, (length(m), length(m)))
    J = Matrix{Float64}(undef, (length(d), length(m)))
    gradient = Vector{Float64}(undef, length(m))
    f = Vector{Float64}(undef, length(d))
    m_update = Vector{Float64}(undef, length(m))
    trial_m = Vector{Float64}(undef, length(m))

    # set initial lambda according to Newman and Alaumbaugh 1997
    # kemna 2000 p 67
    evaluateJacobian!(J, G, m)
    lambda = mean(diag(J' * Cinv * J))

    # In the first loop, inversions are performed and lambda is adjusted
    # to reach the target data fit
    for iteration = 1:max_iter
        try
            rmse = inversion!(m, m_update, trial_m, d, f, G, Cinv,
                J, Hessian, gradient, lambda, 1.0, true, 100, R)
        catch e
            @warn "Inversion stopped due to ", e
            break
        end

        rmse_reached = rmse < (target_rmse + 1e-3)
        rmse_not_improved = abs(rmse - rmse0) < 1e-4
        if rmse_reached || rmse_not_improved
            break
        end

        # update lambda. Decrease if underfitted and increase of overfitted
        if rmse > target_rmse
            lambda = lambda / 1.1
        else
            lambda = lambda * 1.1
        end
        rmse0 = rmse
    end

    # In the second loop, inversions are performed and lambda is increased
    # until the rmse increases significantly.
    lambda0 = lambda
    rmse0 = maximum([rmse, target_rmse])
    for iteration = 1:max_iter
        try
            rmse = inversion!(m, m_update, trial_m, d, f, G, Cinv,
                J, Hessian, gradient, lambda, 1.0, true, 100, R)
        catch e
            @warn "Inversion stopped due to ", e
            break
        end

        rmse_reached = (rmse - rmse0) > 1e-2
        lambda_upper_bound_reached = lambda0 * 1000 < lambda
        if rmse_reached || lambda_upper_bound_reached
            break
        end
        lambda = lambda * 1.5
    end
    return m, rmse, lambda
end
