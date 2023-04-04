# Import utils and load required libraries
remove_duplicates = function(df, column){
  current = ""
  for(i in 1:nrow(df)) {
    row = df[i,]
    if(row[column] != current){
      current = row[column]
    }
    else{
      df[i, column]= ""
    }
  }
  return(df)
}

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

max_bold = function(df, bold_column){
  current = df[1, "n_obs"]
  position = 1
  largest = as.numeric(df[1, bold_column])
  for(i in 1:nrow(df)) {
    row = df[i,]
    if(row["n_obs"] != current){
      df[position, bold_column] = bold(df[position, bold_column])
      current = row["n_obs"]
      largest = as.numeric(row[bold_column])
      position = i
    }
    else if(as.numeric(row[bold_column]) > largest){
      largest = as.numeric(row[bold_column])
      position = i
    }
  }
  return(df)
}

make_output_table = function(df, output_fn, table_alignment){
  df = remove_duplicates(df, "Location")
  df = remove_duplicates(df, "Model")
  df = remove_duplicates(df, "Pathogen")
  latex_of_table = xtable(df, type = "latex", auto = TRUE)
  align(latex_of_table) <- table_alignment
  print(latex_of_table, tabular.environment = "tabular", file = output_fn, include.rownames=FALSE)
}




source("outputs/plot_utils.R")


### Data imports ###
data_dir = "data/offspring/"
output_dir = "outputs/offspring/baseline/"   # Edit the final directory to obtain sensitivity results (e.g., baseline -> k5)

# Offspring data
offspring_files = list.files(path=data_dir, pattern ="*.csv")
offspring_datasets = lapply(paste0(data_dir, offspring_files), read.csv)
names(offspring_datasets) = gsub('.{4}$', '', offspring_files)
offspring = bind_rows(offspring_datasets, .id="Dataset")
offspring = cbind(str_split(offspring$Dataset, "_", n=3, simplify=TRUE), offspring)
names(offspring) = c("Pathogen", "Location", "Author", "Dataset", "Z", "n")
offspring$Pathogen = factor(offspring$Pathogen, levels=pathogen_levels)
offspring$Label = relabel_offspring(offspring$Dataset)

# Model predictions
model_fits = read.csv(paste0(output_dir, "model_fit.csv"))
model_fits = cbind(str_split(model_fits$Dataset, "_", n=3, simplify=TRUE), model_fits)
names(model_fits) = c("Pathogen", "Location", "Author", "Dataset", "Model", "Z", "n")
model_fits$Pathogen = factor(model_fits$Pathogen, levels=pathogen_levels)
model_fits$Model = relabel_model(model_fits$Model)
model_fits$Label = relabel_offspring(model_fits$Dataset)

# Model performance
score_summary = read.csv(paste0(output_dir, "score_summary.csv"), stringsAsFactors = TRUE)
score_summary = cbind(str_split(score_summary$Dataset, "_", n=3, simplify=TRUE), score_summary)
names(score_summary)[1:3] = c("Pathogen", "Location", "Author")
score_summary$Label = relabel_offspring(score_summary$Dataset)
score_summary$Model = relabel_model(score_summary$Model)
score_summary$Pathogen = factor(score_summary$Pathogen, levels=pathogen_levels)

# Model parameters
parms_summary = read.csv(paste0(output_dir, "parm_summary.csv"), stringsAsFactors = TRUE)
parms_summary = cbind(str_split(parms_summary$Dataset, "_", n=3, simplify=TRUE), parms_summary)
names(parms_summary)[1:3] = c("Pathogen", "Location", "Author")
parms_summary$Label = relabel_offspring(parms_summary$Dataset)
parms_summary$Model = relabel_model(parms_summary$Model)
parms_summary$Pathogen = factor(parms_summary$Pathogen, levels=pathogen_levels)

# Parameter chains
full_chain = read.csv(paste0(output_dir, "full_chain.csv"))
full_chain = cbind(str_split(full_chain$Dataset, "_", n=3, simplify=TRUE), full_chain)
names(full_chain)[1:3] = c("Pathogen", "Location", "Author")
full_chain$Label = relabel_offspring(full_chain$Dataset)
full_chain$Model = relabel_model(full_chain$Model)
full_chain$Pathogen = factor(full_chain$Pathogen, levels=pathogen_levels)


### Mutate data ###
# Aggregate ≥15 offspring
Z_range = c(0,15)
Zalt_labels = c(seq(0,14), "≥15")
Zalt_ticks = rep("", length(Zalt_labels))
Zalt_ticks[seq(1, length(Zalt_labels), by=5)] = Zalt_labels[seq(1, length(Zalt_labels), by=5)]
offspring$Zalt = factor(cut(offspring$Z, breaks = c(seq(Z_range[1] - 0.5, Z_range[2] - 0.5, by=1), Inf), labels=Zalt_labels), levels=Zalt_labels)
model_fits$Zalt = factor(cut(model_fits$Z, breaks = c(seq(Z_range[1] - 0.5, Z_range[2] - 0.5, by=1), Inf), labels=Zalt_labels), levels=Zalt_labels)
summarize_df = expand.grid(Label = offspring_labels, Zalt = Zalt_labels)

offspring_accum = summarize_df %>%
  left_join(offspring, by = c("Label", "Zalt")) %>%
  group_by(Label, Zalt, Pathogen) %>%
  summarize(count = sum(n), .groups = "drop") %>%
  replace_na(list(count = 0))

model_accum = summarize_df %>%
  left_join(model_fits, by = c("Label", "Zalt")) %>%
  group_by(Model, Label, Zalt, Pathogen) %>%
  summarize(count = sum(n), .groups = "drop") %>%
  replace_na(list(count = 0)) %>%
  filter(Model %in% c("Negative Binomial", "Mixture", "Clinical"))

max.offspring = offspring %>%
  group_by(Dataset) %>%
  summarize(Max.n = max(n))

score_text = score_summary %>%
  filter(Model %in% c("Negative Binomial", "Mixture", "Clinical")) %>%
  mutate(n = as.vector(sapply(max.offspring$Max.n, function(x) x * c(0.9, 0.7, 0.5))),
         wr = paste0("w = ", round(w, 2)),
         aiccs = paste0("AICc = ", signif(AICc, 4)))

score_text$Zalt = 10.5



#### Model fits ####

ggplot(model_accum, aes(x=Zalt)) +
  geom_bar(data=offspring_accum, aes(weight=count, fill=Pathogen), alpha=1.) + 
  geom_point(aes(y = count, col=Model, shape=Model), size=3) +
  geom_line(aes(y = count, group=Model, col=Model, lty=Model), lwd=1) + 
  geom_text(data=score_text, aes(y=n, label=aiccs, col=Model), size=6, family="serif") +
  facet_wrap(Label ~ ., scales="free") +
  xlab("Secondary cases (Z)") +
  ylab("Frequency") +
  scale_x_discrete(labels = Zalt_ticks) +
  scale_y_continuous(label=scales::comma) +
  custom_theme + 
  scale_fill_pathogen() +
  scale_col_model() + 
  scale_shape_model() +
  guides(fill = guide_legend(nrow=4),
         col = guide_legend(nrow=3))


ggsave(paste0(output_dir, "offspring_fits.png"), dpi=600, width=14, height=14)




#### Parameter plots ####

## Plot R for each model * dataset ##
parms_summary %>%
  filter(Parameter=="R") %>%
  ggplot(aes(x = Label, col=Pathogen, shape=Model)) +
  geom_linerange(aes(ymin = X2.5., ymax=X97.5.),position=position_dodge2(width=0.75), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = X25., ymax=X75.), position=position_dodge2(width=0.75), lwd=2, alpha=0.6) +
  geom_point(aes(y=X50.), position=position_dodge2(width=0.75), size=3) +
  geom_hline(yintercept = 1, lty=2, lwd=1) +
  xlab("") +
  ylab("Reproductive number (R)") +
  custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = c(0.2,0.65)) +
  scale_col_pathogen() +
  scale_shape_model() +
  guides(col = guide_legend(nrow=4))


ggsave(paste0(output_dir, "R_estimates.png"), dpi=600, width=12, height=7.5)


### Dispersion parameter ###
parms_summary %>%
  filter(Parameter=="α") %>%
  ggplot(aes(x = Label, col=Pathogen, shape=Model)) +
  geom_linerange(aes(ymin = X2.5., ymax=X97.5.),position=position_dodge2(width=0.75), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = X25., ymax=X75.), position=position_dodge2(width=0.75), lwd=2, alpha=0.6) +
  geom_point(aes(y=X50.), position=position_dodge2(width=0.75), size=3) +
  xlab("") +
  ylab(expression(Dispersion~(k[NB]))) +
  # ylab(Tex("Dispersion parameter, $k_{NB}$")) +
  custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(trans="log10", breaks = c(0.01, 0.03, 0.1, 0.3, 1., 3.)) +
  # scale_y_log10(
  #   # breaks = c(0.01, 0.03, 0.1, 0.3, 1., 3.),
  #   breaks = scales::trans_breaks("log10", function(x) 10^x),
  #   labels = scales::trans_format("log10", scales::math_format(10^.x))
  # ) +
  # coord_trans(y = "log10") +
  annotation_logticks(sides="l") +
  scale_col_pathogen() +
  scale_shape_model()


ggsave(paste0(output_dir, "k_negbin.png"), dpi=600, width=12, height=7.5)



### Superspreader fraction + relative transmissibility ###
parms_summary %>%
  filter(Parameter==c("c", "ρ"), Model %in% c("Mixture")) %>%
  ggplot(aes(x = Label, col=Pathogen, shape=Parameter)) +
  geom_linerange(aes(ymin = X2.5., ymax=X97.5.),position=position_dodge2(width=0.5), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = X25., ymax=X75.), position=position_dodge2(width=0.5), lwd=2, alpha=0.6) +
  geom_point(aes(y=X50.), position=position_dodge2(width=0.5), size=3) +
  xlab("") +
  ylab("Superspreader fraction (c)\nRelative transmissibility (ρ)") +
  custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = c(0.25,0.7)) +
  # scale_y_continuous(trans="log10") +
  scale_y_continuous(labels=scales::percent) +
  scale_col_pathogen() +
  guides(col = guide_legend(nrow=4))


ggsave(paste0(output_dir, "c_rho.png"), dpi=600, width=12, height=7.5)



#### Derived parameters ####
### Extinction probability
full_chain %>%
  filter(Model != "Erlang") %>%
  group_by(Label, Pathogen, Model) %>%
  summarize(lower025 = quantile(q, 0.025),
            lower25 = quantile(q, 0.25),
            median = median(q),
            upper75 = quantile(q, 0.75),
            upper975 = quantile(q, 0.975)) %>%
  ggplot(aes(x = Label, col=Pathogen, shape=Model)) +
  geom_point(aes(y=median), position=position_dodge2(width=0.75), size=3) +
  geom_linerange(aes(ymin = lower025, ymax=upper975), position=position_dodge2(width=0.75), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = lower25, ymax=upper975), position=position_dodge2(width=0.75), lwd=2, alpha=0.6) +
  geom_hline(yintercept=0.20923955891032855, lty=2) +
  annotate("text", x=3, y=0.17, size=6, label="Erlang prediction", family="serif") +
  xlab("") +
  ylab("Extinction probability (q)") +
  scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = seq(0,1,0.1)) +
  custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "right") +
  scale_col_pathogen() +
  scale_shape_model()


ggsave(paste0(output_dir, "q_estimates.png"), dpi=600, width=12, height=7.5)




#### Scores ####
## Plot score by model for each dataset ##
ggplot(score_summary, aes(x=Model, col=Pathogen, shape=Model)) +
  geom_point(aes(y=w), size=4) +
  facet_wrap(Label ~ ., nrow=5) +
  xlab("Model") +
  ylab("Akaike weight (w)") +
  custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position="right") +
  scale_shape_model() +
  scale_col_pathogen()


ggsave(paste0(output_dir, "all_scores.png"), dpi=600, height=12, width=16)







##### Supplementary figures ####



#### Tables ####
### Score table
score_table = score_summary %>% dplyr::select(Pathogen, Location, n_obs, Model, BIC, AIC, AICc, w)
score_table = max_bold(score_table, "w")
score_table = remove_duplicates(score_table, "n_obs")
covid_table = score_table %>% dplyr::filter(Pathogen == "SARS-CoV-2")
make_output_table(covid_table, "Scores_table_covid.tex", "rlllrrrrr")
not_covid_table = score_table %>% dplyr::filter(Pathogen != "SARS-CoV-2")
make_output_table(not_covid_table, "Scores_table_not_covid.tex", "rllrlrrrr")

### Parameters table
parms_summary_table = parms_summary %>% dplyr::select(Pathogen, Location, Model, Parameter, mle, X2.5., X25., X50., X75., X97.5., ess, Rhat) %>%
  dplyr::filter(Model %in% c("Negative Binomial", "Mixture")) 
covid_table = parms_summary_table %>% dplyr::filter(Pathogen == "SARS-CoV-2")
make_output_table(covid_table, "Parameters_table_covid.tex", "llccrrrrrrrrr")
not_covid_table = parms_summary_table %>% dplyr::filter(Pathogen != "SARS-CoV-2")
make_output_table(not_covid_table, "Parameters_table_not_covid.tex", "llccrrrrrrrrr")

#parms_table = parms_summary %>%
  #dplyr::select(Pathogen, Location, Model, Parameter, mle, X2.5., X25., X50., X75., X97.5., ess, Rhat) %>%
  #dplyr::filter(Model %in% c("Negative Binomial", "Mixture")) %>%
  #xtable(include.rownames=FALSE)


#### Posterior plots ###
## Plot  R v. α
full_chain %>%
  filter(Model == "Negative Binomial") %>%
  ggplot(aes(x= R, y= α, col=Pathogen)) +
  geom_point(alpha=0.5) +
  # geom_density_2d(bins=4) +
  geom_vline(xintercept=1, lty=2) +
  facet_wrap(Label ~ ., nrow=5) +
  scale_x_continuous(limits=c(0, 8), expand=c(0,0)) +
  scale_y_continuous(trans="log10", limits = c(0.01, 12.)) +
  ylab("Dispersion parameter (k)") +
  xlab("Reproduction number") +
  custom_theme + scale_col_pathogen()

ggsave(paste0(output_dir, "R_v_k_posterior.png"), dpi=600, height=10, width=16)


## Plot R v. c
full_chain %>%
  filter(Model == "Mixture") %>%
  ggplot(aes(x= R, y= c, col=Pathogen)) +
  geom_point(alpha=0.5) +
  # geom_density_2d(bins=4) +
  geom_vline(xintercept=1, lty=2) +
  facet_wrap(Label ~ ., nrow=5) +
  scale_x_continuous(limits=c(0, 8), expand=c(0,0)) +
  ylab("Superspreader fraction (c)") +
  xlab("Reproduction number") +
  custom_theme + scale_col_pathogen()

ggsave(paste0(output_dir, "R_v_c_posterior.png"), dpi=600, height=10, width=16)


## Plot R v. ρ
full_chain %>%
  filter(Model == "Mixture") %>%
  ggplot(aes(x= R, y=ρ, col=Pathogen)) +
  geom_point(alpha=0.5) +
  # geom_density_2d(bins=4) +
  geom_vline(xintercept=1, lty=2) +
  facet_wrap(Label ~ ., nrow=5) +
  scale_x_continuous(limits=c(0, 8), expand=c(0,0)) +
  ylab("Relative transmissibility (ρ)") +
  xlab("Reproduction number") +
  custom_theme + scale_col_pathogen()

ggsave(paste0(output_dir, "R_v_rho_posterior.png"), dpi=600, height=10, width=16)



## Plot c v. ρ
full_chain %>%
  filter(Model == "Mixture") %>%
  ggplot(aes(x= c, y= ρ, col=Pathogen)) +
  geom_point(alpha=0.5) +
  # geom_density_2d(bins=4) +
  geom_vline(xintercept=1, lty=2) +
  facet_wrap(Label ~ ., nrow=5) +
  ylab("Relative transmissibility (ρ)") +
  xlab("Superspreader fraction (c)") +
  custom_theme + scale_col_pathogen()

ggsave(paste0(output_dir, "c_v_rho_posterior.png"), dpi=600, height=10, width=16)
