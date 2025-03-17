using TdipCtools
using LinearAlgebra
using Plots
using Test
using Random
using Distributions
import Logging

# run Base.jl tests
include("benchmark.jl")
include("test_base/test_debyeDecomposition.jl")
include("test_base/test_calculateRmse.jl")
include("test_base/test_createMatrixForwardOperator.jl")
include("test_base/test_decompositionResponseTimeDomain.jl")
include("test_base/test_occamDebyeDecomposition.jl")
include("test_base/test_decompositionResponseFrequencyDomain.jl")

# run Manager.jl tests
include("test_manager.jl")
