# Import utils and load required libraries
source("outputs/plot_utils.R")
source("outputs/table_utils.R")


### Data imports ###
data_dir = "data/offspring/"
output_dir = "outputs/offspring/baseline/"   # Edit the final directory to obtain sensitivity results (e.g., baseline -> k3)


# Offspring data
offspring = read_input(data_dir)

# Model predictions
model_fits = read_output("model_fit")


# Model performance
score_summary = read_output("score_summary")
score_summary = append_ΔAICc(score_summary %>% filter(Model != "Zero-inflated"))



# Model parameters
parms_summary = read_output("parm_summary")

# Parameter chains
full_chain = read_output("full_chain")



#### Model fits ####

plot_fit(offspring, model_fits, score_summary, models=c("Negative Binomial", "Two-type"))  +
  guides(fill = guide_legend(nrow=4),
         col = guide_legend(nrow=2))
ggsave(paste0(output_dir, "offspring_fits.png"), dpi=600, width=14, height=14)


plot_fit(offspring, model_fits, score_summary, models=c("Negative Binomial", "Clinical"))  +
  guides(fill = guide_legend(nrow=4),
         col = guide_legend(nrow=2))
ggsave(paste0(output_dir, "clinical_fits.png"), dpi=600, width=14, height=14)




#### Parameter plots ####

### Dispersion parameter ###

plot_parameter(parms_summary, "α", ylabel=expression(Dispersion~(k[NB])), models=c("Negative Binomial")) +
  scale_y_continuous(trans="log10", breaks = c(0.01, 0.05, 0.1, 0.5, 1., 5.)) +
  annotation_logticks(sides="l") +
  theme(legend.position="right")


ggsave(paste0(output_dir, "k_negbin.png"), dpi=600, width=12, height=7.5)



## Plot R for each model * dataset ##

plot_parameter(parms_summary, "R", ylabel="Reproductive number (R)", models=all_models) +
    geom_hline(yintercept = 1, lty=2, lwd=1) +
  theme(legend.position = "right") 

ggsave(paste0(output_dir, "R_estimates.png"), dpi=600, width=16, height=8)



### Superspreader fraction + relative transmissibility ###

ggplot(parms_summary %>% filter(Parameter %in% c("c", "ρ"), Model %in% c("Two-type", "SEIR(2)")), aes(x=Label, col=Pathogen, shape=Model)) +
  geom_linerange(aes(ymin = X2.5., ymax=X97.5.),position=position_dodge2(width=0.75), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = X25., ymax=X75.), position=position_dodge2(width=0.75), lwd=2, alpha=0.6) +
  geom_point(aes(y=X50.), position=position_dodge2(width=0.75), size=3) +
  facet_grid(Parameter ~ ., labeller=labeller(Parameter = c(c="Superspreader\nfraction (c)", ρ="Relative\ntransmissibility (ρ)"))) +
  xlab("") +
  ylab("") +
  custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = c(0.2,0.65)) +
  scale_col_pathogen() +
  scale_shape_model() +
  theme(legend.position = "right")
  


ggsave(paste0(output_dir, "c_rho.png"), dpi=600, width=12, height=7.5)



plot_parameter(parms_summary, "σ", ylabel="Relative transmissibility (σ)", models=c("Two-type", "Single-type")) +
  theme(legend.position = "right")


ggsave(paste0(output_dir, "sigma.png"), dpi=600, width=12, height=7.5)



#### Derived parameters ####
### Extinction probability

plot_extinction(full_chain, models=all_models)


ggsave(paste0(output_dir, "q_estimates.png"), dpi=600, width=16, height=8)




#### Scores ####
## Plot score by model for each dataset ##

plot_score(score_summary, models=all_models)


ggsave(paste0(output_dir, "all_scores.png"), dpi=600, height=12, width=16)




####### Clinical results #######
clinical_offspring = read_input(dir = "./data/offspring/clinical/")
clinical_offspring$Label = relabel_clinical(clinical_offspring$Dataset)


# Model predictions
clinical_fits = read_output("model_fit", output_dir = "./outputs/offspring/clinical/")
clinical_fits$Label = relabel_clinical(clinical_fits$Dataset)

# Model performance
clinical_scores = read_output("score_summary", output_dir = "./outputs/offspring/clinical/")
clinical_scores$Label = relabel_clinical(clinical_scores$Dataset)



plot_fit(clinical_offspring, clinical_fits, clinical_scores, models=c("Negative Binomial", "Single-type"), Z_max=8, Z_inc=2, labels=clinical_labels, x_text=6) + 
  theme(legend.position="right")

ggsave("./outputs/offspring/clinical/clinical_model_fit_comparison.png", dpi=600, height=6, width=9)


##### Supplementary figures ####

### k=2 ####
# Model predictions
sensk2_fits = read_output("model_fit", output_dir = "./outputs/offspring/k2/")

# Model performance
sensk2_scores = read_output("score_summary", output_dir = "./outputs/offspring/k2/")
sensk2_scores = append_ΔAICc(sensk2_scores)


plot_fit(offspring, sensk2_fits, sensk2_scores, models=c("Negative Binomial", "SEIR(2)"), Z_max=15, Z_inc=5)  +
  guides(fill = guide_legend(nrow=4),
         col = guide_legend(nrow=3))


ggsave("./outputs/offspring/k2/k2_fit_comparison.png", dpi=600, height=6, width=9)

plot_score(sensk2_scores, models=c("Negative Binomial", "SEIR(2)", "SEIR(1)"))


ggsave("./outputs/offspring/k2/k2_all_scores.png", dpi=600, height=12, width=16)


### k=3 ####
# Model predictions
sensk3_fits = read_output("model_fit", output_dir = "./outputs/offspring/k3/")

# Model performance
sensk3_scores = read_output("score_summary", output_dir = "./outputs/offspring/k3/")
sensk3_scores = append_ΔAICc(sensk3_scores %>% filter(Model != "Zero-inflated"))


plot_fit(offspring, sensk3_fits, sensk3_scores, models=c("Negative Binomial", "SEIR(2)"), Z_max=15, Z_inc=5)  +
  guides(fill = guide_legend(nrow=4),
         col = guide_legend(nrow=3))


ggsave("./outputs/offspring/k3/k3_fit_comparison.png", dpi=600, height=6, width=9)

plot_score(sensk3_scores, models=c("Negative Binomial", "SEIR(2)", "SEIR(1)"))

ggsave("./outputs/offspring/k3/k3_all_scores.png", dpi=600, height=12, width=16)



# Check sensitivity of parameter estimates
parms_summary$k = 1
sensk2_params = read_output("parm_summary", output_dir =  "./outputs/offspring/k2/")
sensk2_params$k = 2
sensk3_params = read_output("parm_summary", output_dir =  "./outputs/offspring/k3/")
sensk3_params$k = 3

parms_all = rbind(parms_summary, sensk2_params, sensk3_params)
parms_all$k = factor(parms_all$k)


ggplot(parms_all %>% filter(Parameter %in% c("R", "c", "ρ"), Model == "SEIR(2)"), aes(x=Label, shape=k, col=Pathogen)) +
  geom_linerange(aes(ymin = X2.5., ymax=X97.5.),position=position_dodge2(width=0.75), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = X25., ymax=X75.), position=position_dodge2(width=0.75), lwd=2, alpha=0.6) +
  geom_point(aes(y=X50.), position=position_dodge2(width=0.75), size=3) +
  facet_grid(Parameter ~ Model, scales = "free") +
  ylab("") +
  xlab("") +
  labs(shape="Serial compartments") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "right") +
  scale_col_pathogen()

ggsave("./outputs/offspring/baseline/k_sensitivity.png", dpi=600, height=10, width=10)


ggplot(parms_all %>% dplyr::filter(Model %in% "SEIR(2)", Parameter %in% c("c", "ρ")), aes(x=Label, shape=k, col=Pathogen)) +
  geom_linerange(aes(ymin = X2.5., ymax=X97.5.),position=position_dodge2(width=0.75), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = X25., ymax=X75.), position=position_dodge2(width=0.75), lwd=2, alpha=0.6) +
  geom_point(aes(y=X50.), position=position_dodge2(width=0.75), size=3) +
  facet_grid(Parameter ~ .) +
  ylab("") +
  xlab("") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "right") +
  scale_col_pathogen()


#### Tables ####
### Score table
score_summary %>% filter(Pathogen == "SARS-CoV-2") %>% make_score_table()
score_summary %>% filter(Pathogen != "SARS-CoV-2") %>% make_score_table()

### Parameters table
parms_summary %>% filter(Pathogen == "SARS-CoV-2", Model %in% c("Negative Binomial", "Two-type")) %>% make_params_table()
parms_summary %>% filter(Pathogen != "SARS-CoV-2", Model %in% c("Negative Binomial", "Two-type")) %>% make_params_table()



three_type_scores = read_output("threetype_scores")
three_type_scores = rbind(three_type_scores, (score_summary %>% filter(Model == "Negative Binomial"))[, c(1:11,13)])
three_type_scores = append_ΔAICc(three_type_scores) %>% arrange(Dataset, Model)
three_type_scores = append_w(three_type_scores)

three_type_scores %>% make_score_table()




sensk1_scores = read_output("score_summary", output_dir =  "./outputs/offspring/baseline/")
sensk1_scores$k = 1
sensk2_scores = read_output("score_summary", output_dir =  "./outputs/offspring/k2/")
sensk2_scores$k = 2
sensk3_scores = read_output("score_summary", output_dir =  "./outputs/offspring/k3/")
sensk3_scores$k = 3

all_scores = rbind(sensk1_scores, sensk2_scores, sensk3_scores)
all_scores$k = factor(all_scores$k)


all_scores %>% 
  filter(Model == "SEIR(2)") %>%
  group_by(Dataset) %>%
  mutate(Rel = AICc - min(AICc)) %>%
  ungroup() %>%
  ggplot(aes(x = k,  y=Rel)) +
  geom_col(aes(fill = Pathogen), position=position_dodge2(width=0.75)) +
  facet_wrap(Label ~ ., scales="free") +
  # geom_point(size=2, position=position_dodge2(width=0.75)) +
  xlab("Serial compartments") +
  ylab("Absolute change in AICc") +
  custom_theme +
  scale_col_pathogen()


ggsave(paste0(output_dir, "k_sensitivity_AICc.png"), dpi=600, width=14, height=14)

