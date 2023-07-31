using Distributions

#### Priors ####
# Specify model priors
R_prior = Gamma(2, 1)
α_prior = Exponential(1)
ρ_prior = Uniform()
c_prior = Uniform()
σ_prior = Uniform()


#### Model likelihoods ####
# Single-type + Negative binomial offspring likelihood
function ℓ_offspring(Z, n, R, α)
    loglik = sum(n) * (α * log(α / (α + R)) - loggamma(α)) + Z' * n * log(R / (α + R))
    for i in eachindex(Z)
        loglik += n[i] * (loggamma(Z[i] + α) - logfactorial(Z[i]))
    end
    return loglik
end

# Multi-type offspring likelihood
function ℓ_offspring(Z, n, R, α, c, ρ)
    p = get_p(R, c, ρ, α)
    loglik = -sum(n) * loggamma(α)
    for l in eachindex(Z)
        s = (1 - c) * p[1]^α * (1 - p[1])^Z[l] + c * p[2]^α * (1 - p[2])^Z[l]
        loglik += n[l] * (log(s) + loggamma(Z[l] + α) - logfactorial(Z[l]))
    end
    return loglik
end


function ℓ_offspring(Z, n, R, α, c₁, c₂, ρ₁, ρ₂)
    p = get_p(R, c₁, c₂, ρ₁, ρ₂, α)
    loglik = -sum(n) * loggamma(α)
    for l in eachindex(Z)
        s = (1. - c₁ - c₂) * p[1]^α * (1 - p[1])^Z[l] + c₁ * p[2]^α * (1 - p[2])^Z[l] + c₂ * p[3]^α * (1 - p[3])^Z[l]
        loglik += n[l] * (log(s) + loggamma(Z[l] + α) - logfactorial(Z[l]))
    end
    return loglik
end


# Single-type variable
function ℓ_variable(Z, n, R, σ)
    loglik = sum(n) * ( log(1. + σ) - log(1. - σ) - log(R) )
    for i in eachindex(Z)
        loglik += n[i] * log( (R / (1 + σ + R))^(Z[i] + 1) - (σ * R / (1. + σ + σ * R))^(Z[i] + 1) )
    end
    return loglik
end


# Two-type variable
function ℓ_variable(Z, n, R, σ, c, ρ)
    σ == 1. && return ℓ_offspring(Z, n, R, 2, c, ρ)
    R₁, R₂ = split_R(R, c, ρ)
    loglik = sum(n) * ( log(1. + σ) - log(1. - σ) )
    for i in eachindex(Z)
        loglik += n[i] * log( (1. - c) * 1. / R₁ * ( (R₁ / (1 + σ + R₁))^(Z[i] + 1) - (σ * R₁ / (1. + σ + σ * R₁))^(Z[i] + 1) ) + c * 1. / R₂ * ( (R₂ / (1 + σ + R₂))^(Z[i] + 1) - (σ * R₂ / (1. + σ + σ * R₂))^(Z[i] + 1) ) )
    end
    return loglik
end


# Three-type variable
function ℓ_variable(Z, n, R, σ, c₁, c₂, ρ₁, ρ₂)
    # σ == 1. && return ℓ_offspring(Z, n, R, 2, c, ρ)
    R₁, R₂, R₃ = split_R(R, c₁, c₂, ρ₁, ρ₂)
    loglik = sum(n) * ( log(1. + σ) - log(1. - σ) )
    for i in eachindex(Z)
        loglik += n[i] * log( (1. - c₁ - c₂) * 1. / R₁ * ( (R₁ / (1 + σ + R₁))^(Z[i] + 1) - (σ * R₁ / (1. + σ + σ * R₁))^(Z[i] + 1) ) + 
                        c₁ * 1. / R₂ * ( (R₂ / (1 + σ + R₂))^(Z[i] + 1) - (σ * R₂ / (1. + σ + σ * R₂))^(Z[i] + 1) ) +
                            c₂ * 1. / R₃ * ( (R₃ / (1 + σ + R₃))^(Z[i] + 1) - (σ * R₃ / (1. + σ + σ * R₃))^(Z[i] + 1) ))
    end
    return loglik
end


# SEIR1 model
function ℓ(prob::SEIR1Offspring, R)
    @unpack Z, n, α = prob
    return ℓ_offspring(Z, n, R, α)
end


# Negative Binomial model
function ℓ(prob::NegBinOffspring, R, α)
    @unpack Z, n = prob
    return ℓ_offspring(Z, n, R, α)
end


# Zero-inflated model (one non-infector type + one-infector type)
function ℓ(prob::ZeroInfOffspring, R, c)
    @unpack Z, n, α = prob
    return ℓ_offspring(Z, n, R, α, c, 0.)
end


# SEIR2 model (one subspreader type + one superspreader type)
function ℓ(prob::SEIR2Offspring, R, c, ρ)
    @unpack Z, n, α = prob
    return ℓ_offspring(Z, n, R, α, c, ρ)
end


# Clinical model (one subspreader type + one superspreader type with fixed proportion in each)
function ℓ(prob::ClinicalOffspring, R, σ, ρ)
    @unpack Z, n, α, c = prob
    return ℓ_variable(Z, n, R, σ, c, ρ)
end


# Single-type variable model
function ℓ(prob::SingleTypeOffspring, R, σ)
    @unpack Z, n = prob
    return ℓ_variable(Z, n, R, σ)
end


# Two-type variable model
function ℓ(prob::TwoTypeOffspring, R, σ, c, ρ)
    @unpack Z, n = prob
    return ℓ_variable(Z, n, R, σ, c, ρ)
end


# Three-type variable model
function ℓ(prob::ThreeTypeOffspring, R, σ, c₁, c₂, ρ₁, ρ₂)
    @unpack Z, n = prob
    return ℓ_variable(Z, n, R, σ, c₁, c₂, ρ₁, ρ₂)
end


# SEIR(3) model
function ℓ(prob::SEIR3Offspring, R, c₁, c₂, ρ₁, ρ₂)
    @unpack Z, n, α = prob
    return ℓ_offspring(Z, n, R, α, c₁, c₂, ρ₁, ρ₂)
end


#### Fitting functions ####
# Make type callable
function (prob::SEIR1Offspring)(θ)
    @unpack R = θ               # extract the parameters
    ℓ(prob, R) + logpdf(R_prior, R)
end

# Define transformation for parameters
function transformation(prob::SEIR1Offspring)
    as((R = asℝ₊, ))
end


# Make type callable
function (prob::NegBinOffspring)(θ)
    @unpack R, α = θ           # extract the parameters
    ℓ(prob, R, α) + logpdf(R_prior, R) + logpdf(α_prior, α)
end

# Define transformation for parameters
function transformation(prob::NegBinOffspring)
    as((R = asℝ₊, α = asℝ₊))
end


# Make the type callable
function (prob::ZeroInfOffspring)(θ)
    R, c = θ    # extract parameters
    ℓ(prob, R, c) + logpdf(R_prior, R) + logpdf(c_prior, c)
end

# Define transformation for parameters
function transformation(prob::ZeroInfOffspring)
    as((R = asℝ₊, c = as𝕀))
end


# Make the type callable
function (prob::SEIR2Offspring)(θ)
    R, c, ρ = θ    # extract parameters
    ℓ(prob, R, c, ρ) + logpdf(R_prior, R) + logpdf(c_prior, c) + logpdf(ρ_prior, ρ)
end

# Define transformation for parameters
function transformation(prob::SEIR2Offspring)
    as((R = asℝ₊, c = as𝕀, ρ = as𝕀))
end


# Make the type callable
function (prob::ClinicalOffspring)(θ)
    R, σ, ρ = θ    # extract parameters
    ℓ(prob, R, σ, ρ) + logpdf(R_prior, R) + logpdf(σ_prior, σ) + logpdf(ρ_prior, ρ)
end

# Define transformation for parameters
function transformation(prob::ClinicalOffspring)
    as((R = asℝ₊, σ = as𝕀, ρ = as𝕀))
end


# Make the type callable
function (prob::SingleTypeOffspring)(θ)
    R, σ = θ    # extract parameters
    ℓ(prob, R, σ) + logpdf(R_prior, R) + logpdf(σ_prior, σ)
end

# Define transformation for parameters
function transformation(prob::SingleTypeOffspring)
    as((R = asℝ₊, σ = as𝕀))
end


# Make the type callable
function (prob::TwoTypeOffspring)(θ)
    R, σ, c, ρ = θ    # extract parameters
    ℓ(prob, R, σ, c, ρ) + logpdf(R_prior, R) + logpdf(σ_prior, σ) + logpdf(c_prior, c) + logpdf(ρ_prior, ρ)
end

# Define transformation for parameters
function transformation(prob::TwoTypeOffspring)
    as((R = asℝ₊, σ = as𝕀, c = as𝕀, ρ = as𝕀))
end




## Fitting routines
function fit_mcmc(prob::OffspringProblem; t=transformation(prob), iter=1_000, n_chains=5, reporter=NoProgressReport())
    P = TransformedLogDensity(t, prob)
    ∇P = ADgradient(:ForwardDiff, P)
    results = [mcmc_with_warmup(Random.default_rng(), ∇P, iter, reporter=reporter) for _ in 1:n_chains]
    posterior = TransformVariables.transform.(t, eachcol(pool_posterior_matrices(results)))
    ess, R̂ = ess_rhat(stack_posterior_matrices(results))
    return (posterior=posterior, ess=ess, R̂ = R̂)
end


function fit_mle(prob::SEIR1Offspring)
    @unpack Z, n, α = prob
    opt = Opt(:LN_SBPLX, 1)
    opt.lower_bounds = [0.]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n)])
end


function fit_mle(prob::NegBinOffspring)
    @unpack Z, n = prob
    opt = Opt(:LN_SBPLX, 2)
    opt.lower_bounds = [0., 0.]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 1.])
end


function fit_mle(prob::ZeroInfOffspring)
    @unpack Z, n, α = prob
    opt = Opt(:LN_SBPLX, 2)
    opt.lower_bounds = [0., 0.]
    opt.upper_bounds = [Inf, 1.]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.5])
end


function fit_mle(prob::SEIR2Offspring)
    @unpack Z, n, α = prob
    opt = Opt(:LN_SBPLX, 3)
    opt.lower_bounds = [0., 0., 0.]
    opt.upper_bounds = [Inf, 1., 1.]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.5, 0.3])
end


function fit_mle(prob::ClinicalOffspring)
    @unpack Z, n, α, c = prob
    opt = Opt(:LN_SBPLX, 3)
    opt.lower_bounds = [0., 0., 0.]
    opt.upper_bounds = [Inf, 1., 1.]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.1, 0.3])
end


function fit_mle(prob::SingleTypeOffspring)
    @unpack Z, n = prob
    opt = Opt(:LN_SBPLX, 2)
    opt.lower_bounds = [0., 0.]
    opt.upper_bounds = [Inf, 1.]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.1])
end


function fit_mle(prob::TwoTypeOffspring)
    @unpack Z, n = prob
    opt = Opt(:LN_SBPLX, 4)
    opt.lower_bounds = [0., 0., eps(), eps()]
    opt.upper_bounds = [Inf, 1., 1. - eps(), 1. - eps()]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.1, 0.1, 0.01])
end


function fit_mle(prob::ThreeTypeOffspring)
    @unpack Z, n = prob
    opt = Opt(:LN_SBPLX, 6)
    opt.lower_bounds = [0., 0., eps(), eps(), eps(), eps()]
    opt.upper_bounds = [Inf, 0.01, 1. - eps(), 1. - eps(), 1. - eps(), 1. - eps()]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.0, 0.1, 0.1, 0.2, 0.2])
end


function fit_mle(prob::SEIR3Offspring)
    @unpack Z, n, α = prob
    opt = Opt(:LN_SBPLX, 5)
    opt.lower_bounds = [0., eps(), eps(), eps(), eps()]
    opt.upper_bounds = [Inf, 1. - eps(), 1. - eps(), 1. - eps(), 1. - eps()]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.1, 0.1, 0.2, 0.2])
end



function fit(prob::SEIR1Offspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = SEIR1Parameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [SEIR1Parameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end


function fit(prob::NegBinOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = NegBinParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [NegBinParameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end


function fit(prob::ZeroInfOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = ZeroInfParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [ZeroInfParameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end


function fit(prob::SEIR2Offspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = SEIR2Parameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [SEIR2Parameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end


function fit(prob::ClinicalOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = ClinicalParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [ClinicalParameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end


function fit(prob::SingleTypeOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = SingleTypeParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [SingleTypeParameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end


function fit(prob::TwoTypeOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = TwoTypeParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [TwoTypeParameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end



function fit_ensemble(dataset, α_table, c_table; dir=".\\data\\offspring\\", models=[:NegBin, :TwoType, :Clinical], iter=1_000, n_chains=5, reporter=NoProgressReport())
    # Read in data
    data = CSV.read(dir*dataset*".csv", DataFrame)
    pathogen, location, author, _ = split(dataset, "_")

    # Datatypes
    Offspring = Dict(:NegBin => NegBinOffspring, :TwoType => TwoTypeOffspring, 
                     :Clinical => ClinicalOffspring, :SEIR2 => SEIR2Offspring, :ZeroInf => ZeroInfOffspring, 
                     :SingleType => SingleTypeOffspring, :SEIR1 => SEIR1Offspring)

    # Fit models
    all_sol = Vector{Solution}()
    all_scores = Vector{Scores}()
    for model in models
        sol, scores = fit(Offspring[model](data, α_table[pathogen], c_table[pathogen]), iter=iter, n_chains=n_chains, reporter=reporter)
        push!(all_sol, sol)
        push!(all_scores, scores)
    end

    # Concatenate model results
    parms_summ = vcat(combine_sol.(all_sol)...)
    parms_summ = hcat(DataFrame(Dataset = fill(dataset, nrow(parms_summ))), parms_summ)

    chain = vcat([DataFrame(sol.mcmc, complete=true,  α=α_table[pathogen], c=c_table[pathogen]) for sol in all_sol]...)
    chain = hcat(DataFrame(Dataset = fill(dataset, nrow(chain))), chain)

    # Calculate weights for each model
    weights = get_weights(all_scores)
    score_summ = DataFrame(hcat(fill(dataset, length(models)),
                                models,
                                map(f -> getfield.(all_scores, f), fieldnames(Scores))..., 
                                weights), 
                           vcat([:Dataset, :Model], fieldnames(Scores)..., [:w]))

    # Calculate model fits
    n_tot = sum(data.n)
    Z_max = 150
    Z_array = collect(0:Z_max)
    
    model_fit = vcat([DataFrame(Model = fill(models[idx], Z_max+1), Z=Z_array, n=predict(Z_array, n_tot, all_sol[idx].mle, α_table[pathogen], c_table[pathogen])) for idx in eachindex(models)]...)
    model_fit = hcat(DataFrame(Dataset = fill(dataset, nrow(model_fit))), model_fit)                 
    return parms_summ, score_summ, chain, model_fit
end

function predict(Z::Vector{Int64}, n::Int64, parms::NegBinParameters, α, c)
    return exp.(ℓ_offspring.(Z, 1, parms.R, parms.α)) .* n
end

function predict(Z::Vector{Int64}, n::Int64, parms::SEIR1Parameters, α, c)
    return exp.(ℓ_offspring.(Z, 1, parms.R, α)) .* n
end

function predict(Z::Vector{Int64}, n::Int64, parms::ZeroInfParameters, α, c)
    return exp.(ℓ_offspring.(Z, 1, parms.R, α, parms.c, 0.)) .* n
end

function predict(Z::Vector{Int64}, n::Int64, parms::SEIR2Parameters, α, c)
    return exp.(ℓ_offspring.(Z, 1, parms.R, α, parms.c, parms.ρ)) .* n
end

function predict(Z::Vector{Int64}, n::Int64, parms::ClinicalParameters, α, c)
    return exp.(ℓ_offspring.(Z, 1, parms.R, α, c, parms.ρ)) .* n
end

function predict(Z::Vector{Int64}, n::Int64, parms::SingleTypeParameters, α, c)
    return exp.(ℓ_variable.(Z, 1, parms.R, parms.σ)) .* n
end

function predict(Z::Vector{Int64}, n::Int64, parms::TwoTypeParameters, α, c)
    return exp.(ℓ_variable.(Z, 1, parms.R, parms.σ, parms.c, parms.ρ)) .* n
end