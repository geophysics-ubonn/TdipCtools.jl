function convert_files(
    timesteps_file::String,
    data_file::String,
    inverse_variances_file::String,
    frequencies_file::String,
    output_file::String,
    R0_file::String
)::Nothing

    @assert isfile(timesteps_file) "Timesteps file not found!"
    @assert isfile(data_file) "Data file not found!"
    @assert isfile(inverse_variances_file) "Inverse variances file not found!"
    @assert isfile(frequencies_file) "Frequencies file not found!"
    @assert isfile(R0_file) "R0 file is not found!"

    timesteps = readdlm(timesteps_file)
    d = readdlm(data_file)
    Cinv = readdlm(inverse_variances_file)
    frequencies = readdlm(frequencies_file)
    R0 = readdlm(R0_file)

    @assert size(timesteps) == size(d) "Files not of equal size!"
    @assert size(timesteps) == size(Cinv) "Files not of equal size!"

    if size(timesteps, 2) == 1
        timesteps = Float64[readdlm(timesteps_file)...]
        d = Float64[readdlm(data_file)...]
        Cinv = Float64[readdlm(inverse_variances_file)...]
        frequencies = Float64[readdlm(frequencies_file)...]
        R0 = R0[1]
    
        manager = initializeManager(d, timesteps, Cinv, R0=R0)
        occamManager!(manager)
        setFrequenciesManager!(manager, frequencies)
        spectrum = spectrumManager(manager)
    
        open(output_file, "w") do io
            writedlm(io, spectrum)
        end
    else
        @assert typeof(timesteps) == Matrix{Float64}
        @assert typeof(d) == Matrix{Float64}
        @assert typeof(Cinv) == Matrix{Float64}

        frequencies = Float64[readdlm(frequencies_file)...]
        R0 = Float64[readdlm(R0_file)...]

        _Cinv = []
        for i in 1:size(timesteps,2)
            push!(_Cinv, Matrix{Float64}(I, (size(d, 1), size(d, 1))) .* Cinv[:, i])
        end
        Cinv = cat(_Cinv..., dims=3)

        tomo_manager = initializeTomoManager(d, timesteps, Cinv, R0)
        occamTomoManager!(tomo_manager)
        setFrequenciesTomoManager!(tomo_manager, frequencies)
        spectra = transpose(spectrumTomoManager(tomo_manager))

        open(output_file, "w") do io
            writedlm(io, spectra)
        end
    end

    return nothing
end