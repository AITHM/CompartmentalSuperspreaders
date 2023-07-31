struct Scores
    n_obs::Int64
    n_pars::Int64
    ℓₘₐₓ::Float64
    BIC::Float64
    AIC::Float64
    AICc::Float64
end


Scores(n_obs, n_pars, ℓₘₐₓ) = Scores(n_obs, n_pars, ℓₘₐₓ, compute_scores(ℓₘₐₓ, n_pars, n_obs)...)


function Base.show(io::IO, scores::Scores)
    printblue("# obs : ", scores.n_obs)
    printblue("# pars: ", scores.n_pars)
    printblue("ℓₘₐₓ: ", scores.ℓₘₐₓ)
    printblue("BIC : ", scores.BIC)
    printblue("AIC : ", scores.AIC)
    printblue("AICc: ", scores.AICc)
end



function get_weights(scores::Vector{Scores})
    AICcs = getfield.(scores, :AICc)
    ΔAICc = AICcs .- minimum(AICcs)
    return exp.(-0.5 * ΔAICc) ./ sum(exp.(-0.5 * ΔAICc))
end



#### Parameters ####
abstract type Parameters end

struct SEIR1Parameters <: Parameters
    R::Float64
end


struct NegBinParameters <: Parameters
    R::Float64
    α::Float64
end


struct ZeroInfParameters <: Parameters
    R::Float64
    c::Float64
end


struct SEIR2Parameters <: Parameters
    R::Float64
    c::Float64
    ρ::Float64
end


struct ClinicalParameters <: Parameters
    R::Float64
    σ::Float64
    ρ::Float64
end


struct SingleTypeParameters <: Parameters
    R::Float64
    σ::Float64
end


struct TwoTypeParameters <: Parameters
    R::Float64
    σ::Float64
    c::Float64
    ρ::Float64
end


struct ThreeTypeParameters <: Parameters
    R::Float64
    σ::Float64
    c₁::Float64
    c₂::Float64
    ρ₁::Float64
    ρ₂::Float64
end


struct SEIR3Parameters <: Parameters
    R::Float64
    c₁::Float64
    c₂::Float64
    ρ₁::Float64
    ρ₂::Float64
end


function extract_values(parms::T) where T <: Parameters
    return [getfield(parms, fname) for fname in fieldnames(T)]
end 


function summarize(parms::Vector{T}; quantiles=[0.025, 0.25, 0.5, 0.75, 0.975]) where T <: Parameters
    values = hcat(extract_values.(parms)...)
    return mapslices(x -> quantile(x, quantiles), values, dims=2)
end



struct Solution{T}
    mle::T
    mcmc::Vector{T}
    ess::Vector{Float64}
    R̂::Vector{Float64}
    # summary::Matrix{Float64}
end


# function Solution(mle::T, mcmc::Vector{T}, ess::Vector{Float64}, R̂::Vector{Float64}) where T <: Parameters
#     mle = extract_values(mle)
#     mcmc_summary = summarize(mcmc) 
#     combined = hcat(mle, mcmc_summary, ess, R̂)
#     return Solution(mle, mcmc, ess, R̂, combined)
# end


function combine_sol(sol::Solution)
    T = typeof(sol.mle)
    mle = extract_values(sol.mle)
    mcmc_summary = summarize(sol.mcmc) 
    combined = DataFrame(hcat(mle, mcmc_summary, sol.ess, sol.R̂), ["mle", "2.5%", "25%", "50%", "75%", "97.5%", "ess", "Rhat"])
    return hcat(DataFrame(Model = fill(chop(string(T), tail=10), length(mle)),
                          Parameter = [fname for fname in fieldnames(T)]), 
                          combined)
end


function Base.show(io::IO, sol::Solution)
    T = typeof(sol.mle)
    mle = extract_values(sol.mle)
    mcmc_summary = summarize(sol.mcmc) 
    combined = hcat(mle, mcmc_summary, sol.ess, sol.R̂)
    printmat(combined, rowNames=fieldnames(T), colNames=["mle", "2.5%", "25%", "50%", "75%", "97.5%", "ess", "Rhat"])
end


import DataFrames.DataFrame


function DataFrame(chain::Vector{SEIR1Parameters}; complete=true, α=1, c=1.)
    n = length(chain)
    df = DataFrame(hcat(fill("SEIR1", n), hcat(extract_values.(chain)...)'), vcat(:Model, [fname for fname in fieldnames(SEIR1Parameters)]))
    df[!, :q] = extinction_prob.(chain, α)
    if complete
        df[!, :α] = fill(-1., n)
        df[!, :c] = fill(-1., n)
        df[!, :ρ] = fill(-1., n)
        df[!, :σ] .= -1.
    end
    return df
end


function DataFrame(chain::Vector{NegBinParameters}; complete=true, α=1, c=1.)
    n = length(chain)
    df = DataFrame(hcat(fill("NegBin", n), hcat(extract_values.(chain)...)'), vcat(:Model, [fname for fname in fieldnames(NegBinParameters)]))
    df[!, :q] = extinction_prob.(chain, df.α)
    if complete
        df[!, :c] = fill(-1., n)
        df[!, :ρ] = fill(-1., n)
        df[!, :σ] .= -1.
    end
    return df
end


function DataFrame(chain::Vector{ZeroInfParameters}; complete=true, α=1, c=1.)
    n = length(chain)
    df = DataFrame(hcat(fill("ZeroInf", n), hcat(extract_values.(chain)...)'), vcat(:Model, [fname for fname in fieldnames(ZeroInfParameters)]))
    df[!, :q] = extinction_prob.(chain, α)
    if complete
        df[!, :α] = fill(-1., n)
        df[!, :ρ] = fill(-1., n)
        df[!, :σ] .= -1.
    end
    return df
end


function DataFrame(chain::Vector{SEIR2Parameters}; complete=true, α=1, c=1.)
    n = length(chain)
    df = DataFrame(hcat(fill("SEIR2", n), hcat(extract_values.(chain)...)'), vcat(:Model, [fname for fname in fieldnames(SEIR2Parameters)]))
    df[!, :q] = extinction_prob.(chain, α)
    if complete
        df[!, :α] = fill(-1., n)
        df[!, :σ] .= -1.
    end
    return df
end


function DataFrame(chain::Vector{ClinicalParameters}; complete=true, α=1, c=1.)
    n = length(chain)
    df = DataFrame(hcat(fill("Clinical", n), hcat(extract_values.(chain)...)'), vcat(:Model, [fname for fname in fieldnames(ClinicalParameters)]))
    df[!, :q] = extinction_prob.(chain, α, c)
    if complete
        df[!, :α] = fill(-1., n)
        df[!, :c] = fill(c, n)
        # df[!, :σ] .= -1.
    end
    return df
end


function DataFrame(chain::Vector{SingleTypeParameters}; complete=true, α=1, c=1.)
    n = length(chain)
    df = DataFrame(hcat(fill("SingleType", n), hcat(extract_values.(chain)...)'), vcat(:Model, [fname for fname in fieldnames(SingleTypeParameters)]))
    df[!, :q] = extinction_prob.(chain)
    if complete
        df[!, :α] = fill(-1., n)
        df[!, :c] = fill(c, n)
        df[!, :ρ] = fill(-1., n)
    end
    return df
end


function DataFrame(chain::Vector{TwoTypeParameters}; complete=true, α=1, c=1.)
    n = length(chain)
    df = DataFrame(hcat(fill("TwoType", n), hcat(extract_values.(chain)...)'), vcat(:Model, [fname for fname in fieldnames(TwoTypeParameters)]))
    df[!, :q] = extinction_prob.(chain)
    if complete
        df[!, :α] = fill(-1., n)
    end
    return df
end



#### Problems ####
abstract type Problem end


### Offspring ###
abstract type OffspringProblem <: Problem end

## SEIR1 model
struct SEIR1Offspring <: OffspringProblem
    Z::Vector{Int64}
    n::Vector{Int64}
    α::Int64             # SEIR1 shape parameter
end

SEIR1Offspring(df::DataFrame, α) = SEIR1Offspring(df.Z, df.n, α)
SEIR1Offspring(df::DataFrame, α, c) = SEIR1Offspring(df.Z, df.n, α)

## Negative binomial model
struct NegBinOffspring <: OffspringProblem
    Z::Vector{Int64}
    n::Vector{Int64}
end

NegBinOffspring(df::DataFrame) = NegBinOffspring(df.Z, df.n)
NegBinOffspring(df::DataFrame, α, c) = NegBinOffspring(df.Z, df.n)

## Zero-inflated model
struct ZeroInfOffspring <: OffspringProblem
    Z::Vector{Int64}     # offspring
    n::Vector{Int64}     # frequency of Z
    α::Int64             # SEIR1 shape parameter
end

ZeroInfOffspring(df::DataFrame, α) = ZeroInfOffspring(df.Z, df.n, α)
ZeroInfOffspring(df::DataFrame, α, c) = ZeroInfOffspring(df.Z, df.n, α)

## SEIR2 model
struct SEIR2Offspring <: OffspringProblem
    Z::Vector{Int64}     # offspring
    n::Vector{Int64}     # frequency of Z
    α::Int64             # SEIR1 shape parameter
end

SEIR2Offspring(df::DataFrame, α) = SEIR2Offspring(df.Z, df.n, α)
SEIR2Offspring(df::DataFrame, α, c) = SEIR2Offspring(df.Z, df.n, α)

## Clinical model
struct ClinicalOffspring <: OffspringProblem
    Z::Vector{Int64}     # offspring
    n::Vector{Int64}     # frequency of Z
    α::Int64             # SEIR1 shape parameter
    c::Float64           # clinical fraction
end

ClinicalOffspring(df::DataFrame, α, c) = ClinicalOffspring(df.Z, df.n, α, c)


## SingleType model
struct SingleTypeOffspring <: OffspringProblem
    Z::Vector{Int64}     # offspring
    n::Vector{Int64}     # frequency of Z
end

SingleTypeOffspring(df::DataFrame) = SingleTypeOffspring(df.Z, df.n)
SingleTypeOffspring(df::DataFrame, α, c) = SingleTypeOffspring(df.Z, df.n)


## TwoType model
struct TwoTypeOffspring <: OffspringProblem
    Z::Vector{Int64}     # offspring
    n::Vector{Int64}     # frequency of Z
end

TwoTypeOffspring(df::DataFrame) = TwoTypeOffspring(df.Z, df.n)
TwoTypeOffspring(df::DataFrame, α, c) = TwoTypeOffspring(df.Z, df.n)


## ThreeType model
struct ThreeTypeOffspring <: OffspringProblem
    Z::Vector{Int64}    # offspring
    n::Vector{Int64}    # frequency of Z
end

ThreeTypeOffspring(df::DataFrame) = ThreeTypeOffspring(df.Z, df.n)
ThreeTypeOffspring(df::DataFrame, α, c) = ThreeTypeOffspring(df.Z, df.n)


## SEIR3 model
struct SEIR3Offspring <: OffspringProblem
    Z::Vector{Int64}     # offspring
    n::Vector{Int64}     # frequency of Z
    α::Int64             # SEIR1 shape parameter
end

SEIR3Offspring(df::DataFrame, α) = SEIR3Offspring(df.Z, df.n, α)
SEIR3Offspring(df::DataFrame, α, c) = SEIR3Offspring(df.Z, df.n, α)