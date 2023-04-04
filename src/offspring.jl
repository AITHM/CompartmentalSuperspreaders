using Distributions

#### Priors ####
# Specify model priors
R_prior = Gamma(2, 1)
α_prior = Exponential(1)
ρ_prior = Uniform()
c_prior = Uniform()


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


# Erlang model
function ℓ(prob::ErlangOffspring, R)
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


# Mixture model (one subspreader type + one superspreader type)
function ℓ(prob::MixtureOffspring, R, c, ρ)
    @unpack Z, n, α = prob
    return ℓ_offspring(Z, n, R, α, c, ρ)
end


# Clinical model (one subspreader type + one superspreader type with fixed proportion in each)
function ℓ(prob::ClinicalOffspring, R, ρ)
    @unpack Z, n, α, c = prob
    return ℓ_offspring(Z, n, R, α, c, ρ)
end



#### Fitting functions ####
# Make type callable
function (prob::ErlangOffspring)(θ)
    @unpack R = θ               # extract the parameters
    ℓ(prob, R) + logpdf(R_prior, R)
end

# Define transformation for parameters
function transformation(prob::ErlangOffspring)
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
function (prob::MixtureOffspring)(θ)
    R, c, ρ = θ    # extract parameters
    ℓ(prob, R, c, ρ) + logpdf(R_prior, R) + logpdf(c_prior, c) + logpdf(ρ_prior, ρ)
end

# Define transformation for parameters
function transformation(prob::MixtureOffspring)
    as((R = asℝ₊, c = as𝕀, ρ = as𝕀))
end


# Make the type callable
function (prob::ClinicalOffspring)(θ)
    R, ρ = θ    # extract parameters
    ℓ(prob, R, ρ) + logpdf(R_prior, R) + logpdf(ρ_prior, ρ)
end

# Define transformation for parameters
function transformation(prob::ClinicalOffspring)
    as((R = asℝ₊, ρ = as𝕀))
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


function fit_mle(prob::ErlangOffspring)
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


function fit_mle(prob::MixtureOffspring)
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
    opt = Opt(:LN_SBPLX, 2)
    opt.lower_bounds = [0., 0.]
    opt.upper_bounds = [Inf, 1.]
    opt.xtol_rel = 1e-4
    opt.max_objective = (x,grad) -> ℓ(prob, x...)
    return NLopt.optimize(opt, [Z' * n / sum(n), 0.3])
end



function fit(prob::ErlangOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = ErlangParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [ErlangParameters(posterior[i]...) for i in eachindex(posterior)]
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


function fit(prob::MixtureOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = MixtureParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [MixtureParameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end


function fit(prob::ClinicalOffspring; iter=1_000, n_chains=5, reporter=NoProgressReport())
    ℓₘₐₓ, xₘₗₑ, ret = fit_mle(prob)
    mle = ClinicalParameters(xₘₗₑ...)
    posterior, ess, R̂ = fit_mcmc(prob, iter=iter, n_chains=n_chains, reporter=reporter)
    mcmc = [ClinicalParameters(posterior[i]...) for i in eachindex(posterior)]
    return Solution(mle, mcmc, ess, R̂), Scores(sum(prob.n), length(xₘₗₑ), ℓₘₐₓ)
end



function fit_offspring_ensemble(dataset, α_table, c_table; dir=".\\data\\offspring\\", iter=1_000, n_chains=5, reporter=NoProgressReport())
    # Read in data
    data = CSV.read(dir*dataset*".csv", DataFrame)
    pathogen, location, author, _ = split(dataset, "_")

    # Fit all models
    n_sol, n_scores = fit(NegBinOffspring(data), iter=iter, n_chains=n_chains, reporter=reporter)
    e_sol, e_scores = fit(ErlangOffspring(data, α_table[pathogen]), iter=iter, n_chains=n_chains, reporter=reporter)
    z_sol, z_scores = fit(ZeroInfOffspring(data, α_table[pathogen]), iter=iter, n_chains=n_chains, reporter=reporter)
    m_sol, m_scores = fit(MixtureOffspring(data, α_table[pathogen]), iter=iter, n_chains=n_chains, reporter=reporter)
    c_sol, c_scores = fit(ClinicalOffspring(data, α_table[pathogen], c_table[pathogen]), iter=iter, n_chains=n_chains, reporter=reporter)

    # Concatenate model results
    all_sol = [n_sol, e_sol, z_sol, m_sol, c_sol]
    parms_summ = vcat(combine_sol.(all_sol)...)
    parms_summ = hcat(DataFrame(Dataset = fill(dataset, nrow(parms_summ))), parms_summ)

    chain = vcat(DataFrame(n_sol.mcmc, complete=true),
                 DataFrame(e_sol.mcmc, complete=true, α=α_table[pathogen]),
                 DataFrame(z_sol.mcmc, complete=true, α=α_table[pathogen]),
                 DataFrame(m_sol.mcmc, complete=true, α=α_table[pathogen]),
                 DataFrame(c_sol.mcmc, complete=true, α=α_table[pathogen], c=c_table[pathogen]))

    chain = hcat(DataFrame(Dataset = fill(dataset, nrow(chain))), chain)

    # Calculate weights for each model
    all_scores = [n_scores, e_scores, z_scores, m_scores, c_scores]
    weights = get_weights(all_scores)
    score_summ = DataFrame(hcat(fill(dataset, 5),
                                [:NegBin, :Erlang, :ZeroInf, :Mixture, :Clinical],
                                map(f -> getfield.(all_scores, f), fieldnames(Scores))..., 
                                weights), 
                           vcat([:Dataset, :Model], fieldnames(Scores)..., [:w]))

    # Calculate model fits
    n_tot = sum(data.n)
    Z_max = 150
    Z_array = collect(0:Z_max)
    model_fit = vcat(DataFrame(Model = fill(:NegBin, Z_max+1), Z = Z_array, n = exp.(ℓ_offspring.(Z_array, 1, n_sol.mle.R, n_sol.mle.α)) .* n_tot),
                     DataFrame(Model = fill(:Erlang, Z_max+1), Z = Z_array, n = exp.(ℓ_offspring.(Z_array, 1, e_sol.mle.R, α_table[pathogen])) .* n_tot),
                     DataFrame(Model = fill(:ZeroInf, Z_max+1), Z = Z_array, n = exp.(ℓ_offspring.(Z_array, 1, z_sol.mle.R, α_table[pathogen], z_sol.mle.c, 0.)) .* n_tot),
                     DataFrame(Model = fill(:Mixture, Z_max+1), Z = Z_array, n = exp.(ℓ_offspring.(Z_array, 1, m_sol.mle.R, α_table[pathogen], m_sol.mle.c, m_sol.mle.ρ)) .* n_tot),
                     DataFrame(Model = fill(:Clinical, Z_max+1), Z = Z_array, n = exp.(ℓ_offspring.(Z_array, 1, c_sol.mle.R, α_table[pathogen], c_table[pathogen], c_sol.mle.ρ)) .* n_tot))
    model_fit = hcat(DataFrame(Dataset = fill(dataset, nrow(model_fit))), model_fit)                 

    return parms_summ, score_summ, chain, model_fit
end

"""


# Loop over all offspring datasets
α = 2
α_table = Dict("SARS-CoV-2" => α,
               "Ebola" => α,
               "Hantavirus" => α,
               "Measles" => α,
               "MERS-CoV" => α,
               "Mpox" => α,
               "SARS-CoV-1" => α,
               "Smallpox" => α,
               "Tuberculosis" => α)

c_table = Dict("SARS-CoV-2" => 0.649,
               "Ebola" => 0.729,
               "Hantavirus" => 1.,
               "Measles" => 0.5,
               "MERS-CoV" => 0.898,
               "Mpox" => 0.898,
               "SARS-CoV-1" => 0.867,
               "Smallpox" => 1.,
               "Tuberculosis" => 0.496)

offspring_datasets =  [chop(dat, tail=4) for dat in readdir(".\\data\\offspring\\") if last(dat, 3) == "csv"]


out = map(x -> fit_offspring_ensemble(x, α_table, c_table), offspring_datasets)
CSV.write(".\\outputs\\offspring\\parm_summary_k5.csv", vcat(getfield.(out, 1)...))
CSV.write(".\\outputs\\offspring\\score_summary_k5.csv", vcat(getfield.(out, 2)...))
CSV.write(".\\outputs\\offspring\\full_chain_k5.csv", vcat(getfield.(out, 3)...))





# Fit to segregated symptomatic v. asymptomatic data
clinical_datasets =  [chop(dat, tail=4) for dat in readdir(".\\data\\offspring\\clinical\\") if last(dat, 3) == "csv"]


out = map(x -> fit_ensemble(x, dir=".\\data\\offspring\\clinical\\"), clinical_datasets)


CSV.write(".\\outputs\\offspring\\clinical\\parm_summary.csv", vcat(getfield.(out, 1)...))
CSV.write(".\\outputs\\offspring\\clinical\\score_summary.csv", vcat(getfield.(out, 2)...))
CSV.write(".\\outputs\\offspring\\clinical\\full_chain.csv", vcat(getfield.(out, 3)...))
"""