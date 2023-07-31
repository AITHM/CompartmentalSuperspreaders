location_levels = c("Guinea", "Korea", "Zaire", "Beijing", "Batam", "China", "Hong Kong", "India", "Jakarta", "South Korea (a)", "South Korea (b)", "Europe", "West Midlands, UK", "Victoria, Australia")

relabel_locations = function(x, levels=location_levels){
  dplyr::case_when(x == "Guinea" ~ levels[1],
                   x == "Korea" ~ levels[2],
                   x == "Zaire" ~ levels[3],
                   x == "Beijing" ~ levels[4],
                   x == "Batam" ~ levels[5],
                   x == "China" ~ levels[6],
                   x == "HongKong" ~ levels[7],
                   x == "India" ~ levels[8],
                   x == "Jakarta" ~ levels[9],
                   x == "SouthKoreaa" ~ levels[10],
                   x == "SouthKoreab" ~ levels[11],
                   x == "Europe" ~ levels[12],
                   x == "WestMidlands" ~ levels[13],
                   x == "Victoria" ~ levels[14])
}


abbreviated_models = c("NB", "TT", "SEIR(2)", "Clinical", "ST", "SEIR(1)", "Three-type", "SEIR(3)")


abbreviate_models = function(x, levels=abbreviated_models){
  dplyr::case_when(x == "Negative Binomial" ~ levels[1],
                   x == "Two-type" ~ levels[2],
                   x == "SEIR(2)" ~ levels[3],
                   x == "Clinical" ~ levels[4],
                   x == "Single-type" ~ levels[5],
                   x == "SEIR(1)" ~ levels[6],
                   x == "Three-type" ~ levels[7],
                   x == "SEIR(3)" ~ levels[8])
}


make_score_table = function(scores){
  
  scores_txt = scores %>%
    dplyr::select(Pathogen, Location, n_obs, Model, ℓₘₐₓ, BIC, AIC, AICc, w)
  scores_txt$Pathogen = ifelse(duplicated(scores_txt$Pathogen),"", as.character(scores_txt$Pathogen))
  scores_txt$Location = relabel_locations(scores_txt$Location)
  scores_txt$Location = ifelse(duplicated(scores_txt$Location),"",scores_txt$Location)
  scores_txt$n_obs = ifelse(duplicated(scores_txt$n_obs),"",scores_txt$n_obs)
  scores_txt$Model = abbreviate_models(scores_txt$Model)
  colnames(scores_txt) = c("Pathogen", "Location", "$n_{obs}$", "Model", "$\\ell_{\\mathrm{max}}$", "BIC", "AIC", "AIC$_c$", "$w$")
  scores_txt %>% xtable(align = c("l", "l", "l", "r", "c", "r", "r", "r", "r", "r"),
           digits= c(0, 0, 0, 0, 0, -3, -3, -3, -3, 2)) %>%
    print(display = "g", math.style.exponents = TRUE, sanitize.text.function = function(x) {x}, include.rownames=FALSE)
    
}




parameter_levels = c("$R$", "$k_{\\mathrm{NB}}$", "$\\sigma$", "$c$", "$\\rho$")


relabel_parameters = function(x, levels = parameter_levels){
  dplyr::case_when(x == "R" ~ levels[1],
                   x == "α" ~ levels[2],
                   x == "σ" ~ levels[3],
                   x == "c" ~ levels[4],
                   x == "ρ" ~ levels[5])
}


make_params_table = function(params){
  
  params_txt = params %>%
    dplyr::select(Pathogen, Location, Model, Parameter, mle, X2.5., X25., X50., X75., X97.5., ess, Rhat)
  params_txt$Pathogen = ifelse(duplicated(params_txt$Pathogen),"", as.character(params_txt$Pathogen))
  params_txt$Location = relabel_locations(params_txt$Location)
  params_txt$Location = ifelse(duplicated(params_txt$Location),"",params_txt$Location)
  params_txt$Model = abbreviate_models(params_txt$Model)
  params_txt$Parameter = relabel_parameters(params_txt$Parameter)
  colnames(params_txt) = c("Pathogen", "Location", "Model", "Parameter", "MLE", "2.5\\%", "25\\%", "50\\%", "75\\%", "97.5\\%", "ESS", "$\\hat{R}$")
  params_txt %>% xtable(align = c("l", "l", "l", "c", "c", "r", "r", "r", "r", "r", "r", "r", "c"),
                        digits= c(0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, -3, 2)) %>%
    print(display = "g", math.style.exponents = TRUE, sanitize.text.function = function(x) {x}, include.rownames=FALSE)
  
}
