"""
    tabulate(x)

    Generate a frequency table from the vector x.
"""
function tabulate(x)
    dct = countmap(x)
    x, n = collect(keys(dct)), collect(values(dct))
    perm = sortperm(x)
    return x[perm], n[perm]
end


"""
    split_R(R, c, ρ)

    Calculate the type-specific reproductive numbers [R₁, R₂] given the population mean, R,
    the superspreader fraction, c, and the relative transmissibility, ρ.
"""
function split_R(R, c, ρ)
    R₂ = R / ((1. - c) * ρ + c)
    R₁ = R₂ * ρ
    return [R₁, R₂]
end


function split_R(R, c₁, c₂, ρ₁, ρ₂)
    R₃ = R / ((1. - c₁ - c₂) * ρ₁ * ρ₂ + c₁ * ρ₂ + c₂)
    R₂ = ρ₂ * R₃
    R₁ = ρ₁ * R₂
    return [R₁, R₂, R₃]
end


"""
    get_p(c, R, ρ, α)

    Calculate the success parameter p for each type given R, c, ρ, and α.
"""
function get_p(R, c, ρ, α)
    g = R / (α * (ρ * (1. - c) + c))
    return vcat((1. + ρ * g)^(-1), 1. / (1. + g))
end


"""
    compute_scores(ℓₘₐₓ, n_pars, n_obs)

    Compute Akaike and Bayesian information criteria scores given maximum likelihood, ℓₘₐₓ, number of parameters, n_pars, and number of observvations, n_obs.
"""
function compute_scores(ℓₘₐₓ, n_pars, n_obs)
    AIC = 2 * n_pars - 2 * ℓₘₐₓ
    AICc = AIC + 2. * n_pars * (n_pars + 1.) / (n_obs - n_pars - 1.) 
    BIC = n_pars * log(n_obs) - 2 * ℓₘₐₓ
    return (BIC=BIC, AIC=AIC, AICc=AICc)
end