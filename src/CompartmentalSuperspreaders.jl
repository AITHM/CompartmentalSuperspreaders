module CompartmentalSuperspreaders

using CSV
using DataFrames
using Distributions
using DynamicHMC
import ForwardDiff
using LogDensityProblems
import LogDensityProblemsAD: ADgradient
using MCMCDiagnosticTools
using NLopt
using Printf
using Random
using Roots
using SpecialFunctions
using StatsBase
using TransformVariables
using TransformedLogDensities
using UnPack
import CSV.write

include(".\\utils\\structs.jl")
include(".\\utils\\general.jl")
include(".\\utils\\extinction.jl")
include(".\\utils\\printmat.jl")
include(".\\offspring.jl")
export ℓ, ℓ_offspring, Solution, NegBinParameters, ErlangParameters, ZeroInfParameters, MixtureParameters, ClinicalParameters, Variable1Parameters, Variable2Parameters,
NegBinOffspring, ErlangOffspring, ZeroInfOffspring, MixtureOffspring, ClinicalOffspring, Variable1Offspring, Variable2Offspring,
fit_mcmc, fit_mle, fit, fit_offspring_ensemble

end # module CompartmentalSuperspreaders
