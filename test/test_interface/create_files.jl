tau = 0.5
timesteps = 10.0 .^(range(-1, 0, 20))
R0 = 2.0
chargeability = 0.1
gamma = R0 * chargeability

std = 0.001
Random.seed!(1000)

noise = Normal(0, std)
noise_realization = rand(noise, length(timesteps))
d = debyeResponseTimeDomain(timesteps, gamma, tau)
d = d .+ noise_realization

Cinv = ones(length(d)) .* std^-2

# write timesteps
open("test_interface/t.txt", "w") do io
    writedlm(io, timesteps)
end

# write transient
open("test_interface/d.txt", "w") do io
    writedlm(io, d)
end

# write cinv
open("test_interface/cinv.txt", "w") do io
    writedlm(io, Cinv)
end

freqs = 10.0 .^ (range(-2, 1, 100))
# write freqs
open("test_interface/freqs.txt", "w") do io
    writedlm(io, freqs)
end