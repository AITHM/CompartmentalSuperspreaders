using Distributions

#### Priors ####
# Specify model priors
R_prior = Gamma(2, 1)
α_prior = Exponential(1)
ρ_prior = Uniform()
c_prior = Uniform()