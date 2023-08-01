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
export ℓ, ℓ_offspring, Solution, NegBinParameters, SEIR1Parameters, ZeroInfParameters, SEIR2Parameters, ClinicalParameters, SingleTypeParameters, TwoTypeParameters, ThreeTypeParameters, SEIR3Parameters,
NegBinOffspring, SEIR1Offspring, ZeroInfOffspring, SEIR2Offspring, ClinicalOffspring, SingleTypeOffspring, TwoTypeOffspring, ThreeTypeOffspring, SEIR3Offspring,
fit_mcmc, fit_mle, fit, fit_ensemble, predict, compute_scores

end # module CompartmentalSuperspreaders
