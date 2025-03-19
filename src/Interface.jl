function convert_files(
    timesteps_file::String,
    data_file::String,
    inverse_variances_file::String,
    frequencies_file::String,
    output_file::String
)::Nothing

    @assert isfile(timesteps_file) "Timesteps file not found!"
    @assert isfile(data_file) "Data file not found!"
    @assert isfile(inverse_variances_file) "Inverse variances file not found!"
    @assert isfile(frequencies_file) "Frequencies file not found!"

    timesteps = Float64[readdlm(timesteps_file)...]
    d = Float64[readdlm(data_file)...]
    Cinv = Float64[readdlm(inverse_variances_file)...]
    frequencies = Float64[readdlm(frequencies_file)...]

    manager = initializeManager(d, timesteps, Cinv)
    occamManager!(manager)
    setFrequenciesManager!(manager, frequencies)
    spectrum = spectrumManager(manager)

    open(output_file, "w") do io
        writedlm(io, spectrum)
    end

    return nothing
end