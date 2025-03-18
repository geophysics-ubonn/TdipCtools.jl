mutable struct Manager
    d::Vector{Float64}
    timesteps::Vector{Float64}
    Cinv::Matrix{Float64}
    m::Union{Vector{Float64},Nothing}
    G::Union{Matrix{Float64},Nothing}
    tau_grid::Union{Vector{Float64},Nothing}
    frequencies::Union{Vector{Float64},Nothing}
    lambda::Union{Float64,Nothing}
    rmse::Union{Float64,Nothing}
    R0::Union{Float64,Nothing}
end

mutable struct TomoManager
    managers::Array{Manager}
    size::Int
end

function initializeManager(
    d::Vector{Float64},
    timesteps::Vector{Float64},
    Cinv::Matrix{Float64};
    R0::Float64=1.0
)::Manager

    manager = Manager(d, timesteps, Cinv, fill(nothing, 7)...)
    manager.R0 = R0
    return manager
end

function debyeManager!(
    manager::Manager;
    lambda::Float64=1e4,
)::Nothing

    manager.tau_grid = createTauGrid(manager.timesteps, 1.5, 30)
    manager.G = createMatrixForwardOperator(manager.timesteps, manager.tau_grid)
    m = -6 * ones(Float64, length(manager.tau_grid))
    (manager.m, manager.rmse, manager.lambda) = debyeDecomposition(m, manager.d, manager.G, manager.Cinv,
        lambda=lambda, max_iter=1000)
    return nothing
end

function occamManager!(
    manager::Manager
)::Nothing

    manager.tau_grid = createTauGrid(manager.timesteps, 1.5, 30)
    manager.G = createMatrixForwardOperator(manager.timesteps, manager.tau_grid)
    m = -6 * ones(Float64, length(manager.tau_grid))
    (manager.m, manager.rmse, manager.lambda) = occamDebyeDecomposition(m, manager.d, manager.G,
        manager.Cinv, max_iter=100)
    return nothing
end

function initializeTomoManager(
    d::Matrix{Float64},
    timesteps::Matrix{Float64},
    Cinv::Array{Float64,3}
)::TomoManager

    managers = Manager[]
    for i in 1:size(d, 2)
        push!(managers, Manager(d[:, i], timesteps[:, i],
            Cinv[:, :, i], fill(nothing, 7)...))
    end
    return TomoManager(managers, length(managers))
end

function debyeTomoManager!(
    tomo_manager::TomoManager;
    lambda::Float64=1e4
)::Nothing

    @info "Inverting TomoManager of size $(size(tomo_manager.managers, 1)) using $(Threads.nthreads()) threads."

    Threads.@threads for i in 1:tomo_manager.size
        debyeManager!(tomo_manager.managers[i], lambda=lambda)
    end
    return nothing
end

function spectrumManager(
    manager::Manager
)::Vector{ComplexF64}

    return decompositionResponseFrequencyDomain(2 * pi * manager.frequencies,
                exp.(manager.m), manager.tau_grid, manager.R0)
end

function setFrequenciesManager!(
    manager::Manager,
    frequencies::Vector{Float64}
)::Nothing

    manager.frequencies = frequencies
    return nothing
end

function setFrequenciesTomoManager!(
    tomo_manager::TomoManager,
    frequencies::Vector{Float64},
)::Nothing

    for manager in tomo_manager.managers
        setFrequenciesManager!(manager, frequencies)
    end
end

function occamTomoManager!(
    tomo_manager::TomoManager
)::Nothing

    @info "Inverting TomoManager of size $(size(tomo_manager.managers, 1)) using $(Threads.nthreads()) threads."

    Threads.@threads for i in 1:tomo_manager.size
        occamManager!(tomo_manager.managers[i])
    end
    return nothing
end