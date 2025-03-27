module TdipCtools

using Base.Threads
using LinearAlgebra
using Statistics
using SparseArrays
using DelimitedFiles

greet() = print("Hello World!")

include("Base.jl")
export calculateRmse
export createMatrixForwardOperator
export createTauGrid
export debyeResponseTimeDomain
export decompositionResponseFrequencyDomain
export decompositionResponseTimeDomain
export debyeDecomposition
export decompositionResponseFrequencyDomain
export estimateStandardDeviation
export occamDebyeDecomposition

include("Manager.jl")
export Manager
export debyeManager!
export occamManager!
export initializeManager
export TomoManager
export initializeTomoManager
export debyeTomoManager!
export occamTomoManager!
export setFrequenciesManager!
export setFrequenciesTomoManager!
export spectrumManager
export spectrumTomoManager

include("Interface.jl")
export convert_files

end # module TdipCtools
