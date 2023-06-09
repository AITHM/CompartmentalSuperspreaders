"""
    extinction_prob(parms)

    Calculate the extinction probability of an outbreak given the model parameters.
"""

function extinction_prob(parms::ErlangParameters, α) 
    return find_zero((x) -> (1. + 3. / α * (1. - x))^(-α) - x, 1. / 3.)
end


function extinction_prob(parms::NegBinParameters, α)
    return find_zero((x) -> (1. + 3. / parms.α * (1. - x))^(-parms.α) - x, 1. / 3.)
end


function extinction_prob(parms::ZeroInfParameters, α)
    _, R₀₂ = split_R(3., parms.c, 0.)
    return find_zero((x) -> 1. - parms.c + parms.c * (1. + R₀₂ / α * (1. - x))^(-α) - x, 1. / 3.)
end


function extinction_prob(parms::MixtureParameters, α)
    R₀₁, R₀₂ = split_R(3., parms.c, parms.ρ)
    return find_zero((x) -> (1. - parms.c) * (1. + R₀₁ / α * (1. - x))^(-α) + parms.c * (1. + R₀₂ / α * (1. - x))^(-α) - x, 1. / 3.)
end


function extinction_prob(parms::ClinicalParameters, α, c)
    R₀₁, R₀₂ = split_R(3., c, parms.ρ)
    return find_zero((x) -> (1. - c) * (1. + R₀₁ / α * (1. - x))^(-α) + c * (1. + R₀₂ / α * (1. - x))^(-α) - x, 1. / 3.)
end