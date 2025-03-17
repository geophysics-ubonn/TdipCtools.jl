using Random
using Distributions
using Plots
using LinearAlgebra
using TdipCtools

### Creating synthetic data

timesteps = 10.0 .^ (range(-2, 1, 200))

tau = 0.5
r0 = 1.0
chargeability = 0.1
gamma = r0 * chargeability 

std = 0.001
Random.seed!(1000)

noise = Normal(0, std)
noise_realization = rand(noise, length(timesteps))
d = debyeResponseTimeDomain(timesteps, gamma, tau)
d = d .+ noise_realization

p = plot(timesteps, d, xscale=:log10)
savefig(p, "./examples/plots/synthetic_transient.png")

Cinv = Matrix{Float64}(I, (length(d), length(d))) .* (std^-2)

manager = initializeManager(d, timesteps, Cinv)
debyeManager!(manager, lambda=1.0)

p = plot(manager.tau_grid, exp.(manager.m), xscale=:log10)
savefig(p, "./examples/plots/rtd.png")

freqs = 10.0 .^ (range(-3, 3, 50))
setFrequenciesManager!(manager, freqs)
spectrum = spectrumManager(manager)

theoretical_spectrum = decompositionResponseFrequencyDomain(freqs * 2 * pi, Float64[gamma],
                        Float64[tau], r0)

# println(theoretical_spectrum)

p = plot(freqs, real(spectrum), xscale=:log10)
scatter!(p, freqs, real(theoretical_spectrum))
savefig(p, "./examples/plots/real.png")

p = plot(freqs, imag(spectrum), xscale=:log10)
scatter!(p, freqs, imag(theoretical_spectrum))
savefig(p, "./examples/plots/imag.png")