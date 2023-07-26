### Library imports ###

# packages <- c("ggplot2", "dplyr", "tidyverse", "ggsci", "xtable", "labeling", "zipfR", "farver")
packages <- c("tidyverse", "ggsci", "xtable")

for (package in packages) {
  if(!require(package, character.only = TRUE)){
    install.packages(package)
  }
  library(package, character.only = TRUE)
}


### Helper functions ###
`%notin%` = Negate(`%in%`)

read_input = function(dir){
  files = list.files(path=dir, pattern ="*.csv")
  input =  lapply(paste0(dir, files), read.csv)
  names(input) = gsub('.{4}$', '', files)
  input = bind_rows(input, .id="Dataset")
  input = cbind(str_split(input$Dataset, "_", n=3, simplify=TRUE), input)
  names(input) = c("Pathogen", "Location", "Author", "Dataset", "Z", "n")
  input$Pathogen = factor(input$Pathogen, levels=pathogen_levels)
  input$Label = relabel_offspring(input$Dataset)
  return(input)
}

read_output = function(output, output_dir="./outputs/offspring/baseline/"){
  out = read.csv(paste0(output_dir, output, ".csv"))
  out = cbind(str_split(out$Dataset, "_", n=3, simplify=TRUE), out)
  names(out)[1:5] = c("Pathogen", "Location", "Author", "Dataset", "Model")
  out$Pathogen = factor(out$Pathogen, levels=pathogen_levels)
  out$Model = relabel_model(out$Model)
  out$Label = relabel_offspring(out$Dataset)
  return(out)
}


append_ΔAICc = function(scores){
  ΔAICc = scores %>% group_by(Label) %>% summarize(Model = Model, ΔAICc = AICc - min(AICc))
  return(left_join(scores, ΔAICc))
}


append_Δl_max = function(scores){
  Δl_max = scores %>% group_by(Label) %>% summarize(Model = Model, Δl_max = ℓₘₐₓ - max(ℓₘₐₓ))
  return(left_join(scores, Δl_max))
}

compartmental_models = c("Two-type", "Clinical", "SEIR(2)", "Single-type", "SEIR(1)")
all_models = c("Negative Binomial", compartmental_models)


### Aesthetic conventions ###
pathogen_levels = c("SARS-CoV-2", "SARS-CoV-1", "Smallpox", "EBV", "MERS-CoV", "Mpox", "Measles", "Tuberculosis")
pathogen_levels = factor(pathogen_levels, levels = pathogen_levels)


scale_fill_pathogen <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c("#ADB6B6", "#FDAF91", "#42B540","#91D1C2B2", "#E18727FF", "#EE4C97FF", "#ADB6B6", "#FED439"), pathogen_levels), 
    drop=T,
    ...
  )
}


scale_col_pathogen <- function(...){
  ggplot2:::manual_scale(
    'col', 
    values = setNames(c("#ADB6B6", "#FDAF91", "#42B540", "#91D1C2B2", "#E18727FF", "#EE4C97FF", "#ADB6B6", "#FED439"), pathogen_levels), 
    ...
  )
}


offspring_labels = c("Batam, 2020 (n=89)", "China, 2020(a) (n=1,178)","China, 2020(b) (n=2,048)","Hong Kong, 2020 (n=290)",  
                     "India, 2020 (n=88,527)", "Jakarta, 2020 (n=1,199)", "Sth. Korea, 2021(a) (n=344)", "Sth. Korea, 2021(b) (n=1,401)", 
                     "Beijing, 2003 (n=33)", "Singapore, 2003 (n=57)", "Europe, 1958-73 (n=32)", "England, 1966 (n=25)", 
                     "Guinea, 2014 (n=152)", "Rep. of Korea, 2015 (n=185)", "Zaire, 1980-84 (n=147)", "Victoria, 2005-15 (n=2,311)")


relabel_offspring = function(x, levels=offspring_labels) {
  case_when(x == "SARS-CoV-2_Batam_Hasan_" ~ levels[1] ,
            x == "SARS-CoV-2_China_Sun_" ~ levels[2],
            x == "SARS-CoV-2_China_Xu_" ~ levels[3],
            x == "SARS-CoV-2_HongKong_Adam_" ~ levels[4],
            x == "SARS-CoV-2_India_Laxminarayan_" ~ levels[5],
            x == "SARS-CoV-2_Jakarta_Hasan_" ~ levels[6],
            x == "SARS-CoV-2_SouthKoreaa_Lim_" ~ levels[7],
            x == "SARS-CoV-2_SouthKoreab_Lim_" ~ levels[8],
            x == "SARS-CoV-1_Beijing_Shen_" ~ levels[9],
            x == "SARS-CoV-1_Singapore_CDC_" ~ levels[10],
            x == "Smallpox_Europe_Fenner_" ~ levels[11],
            x == "Smallpox_WestMidlands_UKGov_" ~ levels[12],
            x == "EBV_Guinea_Althous_" ~ levels[13],
            x == "MERS-CoV_Korea_Chun_" ~ levels[14],
            x == "Mpox_Zaire_Jezek_" ~ levels[15],
            x == "Tuberculosis_Victoria_Melsew_" ~ levels[16]) %>% factor(levels)
}

clinical_labels = c("Asymptomatic, G1 (n=9)", "Symptomatic, G1 (n=30)", "Asymptomatic, G2-5 (n=44)", "Symptomatic, G2-5 (n=73)")


relabel_clinical = function(x, levels=clinical_labels) {
  case_when(x == "SARS-CoV-2_WanzhouG1a_Shi_" ~ levels[1],
            x == "SARS-CoV-2_WanzhouG1s_Shi_" ~ levels[2],
            x == "SARS-CoV-2_WanzhouG25a_Shi_" ~ levels[3],
            x == "SARS-CoV-2_WanzhouG25s_Shi_" ~ levels[4]) %>% factor(levels)
}



clinical_levels = c("Asymptomatic", "Symptomatic")

assign_clinical_status = function(x, levels=clinical_levels){
  case_when(x == "SARS-CoV-2_WanzhouG1a_Shi_" ~ levels[1] ,
            x == "SARS-CoV-2_WanzhouG1s_Shi_" ~ levels[2],
            x == "SARS-CoV-2_WanzhouG25a_Shi_" ~ levels[1],
            x == "SARS-CoV-2_WanzhouG25s_Shi_" ~ levels[2]) %>% factor(levels)
}


scale_fill_clinical <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c("#00468B", "#ED0000"), clinical_levels), 
    drop=T,
    ...
  )
}


cluster_labels = c("Tianjin, 2020\n(c=43, n=135)", "Alberta, 2020-21\n(c=5,032, n=12,338)",
                   "Florida, 2020-21\n(c=147, n=780)", "Georgia, 2020-21\n(c=68, n=324)",
                   "Illinois, 2020-21\n(c=81, n=281)", "Indiana, 2020-21\n(c=105, n=230)",
                   "Manitoba, 2020-21\n(c=1,754, n=3,782)", "Ohio, 2020-21\n(c=122, n=257)",
                   "Ontario, 2020-21\n(c=8,482, n=15,518)", "Pennsylvania, 2020-21\n(c=322, n=564)",
                   "Saskatchewan, 2020-21\n(c=1,211, n=2,063)", "Texas, 2020-21\n(c=369, n=640)",
                   "Wisconsin, 2020-21\n(c=130, n=522)", "Eurasia, 2012-13\n(c=41, n=111)",
                   "Canada, 1998-2001\n(c=49, n=274)", "US, 1997-99\n(c=165, n=336)",
                  "UK, 2010-15\n(c=12,503, n=23,643)", "Netherlands, 2004-15\n(c=4,926, n=8,398)")


relabel_clusters = function(x, levels=cluster_labels){
  case_when(x == "SARS-CoV-2_Tianjin_Zhang_" ~ levels[1],
            x == "SARS-CoV-2_Alberta_Tupper_" ~ levels[2],
            x == "SARS-CoV-2_Florida_Tupper_" ~ levels[3],
            x == "SARS-CoV-2_Georgia_Tupper_" ~ levels[4],
            x == "SARS-CoV-2_Illinois_Tupper_" ~ levels[5],
            x == "SARS-CoV-2_Indiana_Tupper_" ~ levels[6],
            x == "SARS-CoV-2_Manitoba_Tupper_" ~ levels[7],
            x == "SARS-CoV-2_Ohio_Tupper_" ~ levels[8],
            x == "SARS-CoV-2_Ontario_Tupper_" ~ levels[9],
            x == "SARS-CoV-2_Pennsylvania_Tupper_" ~ levels[10],
            x == "SARS-CoV-2_Saskatchewan_Tupper_" ~ levels[11],
            x == "SARS-CoV-2_Texas_Tupper_" ~ levels[12],
            x == "SARS-CoV-2_Wisconsin_Tupper_" ~ levels[13],
            x == "MERS-CoV_Eurasia_Cauchemez_" ~ levels[14],
           x == "Measles_Canada_Blumberg_" ~ levels[15],
           x == "Measles_US_Blumberg_" ~ levels[16],
           x == "Tuberculosis_UK_Brooks-Pollock_" ~ levels[17],
           x == "Tuberculosis_Netherlands_Brooks-Pollock_" ~ levels[18]) %>% factor(levels)
}


model_levels = c("Negative Binomial", "Two-type",  "SEIR(2)", "Zero-inflated", "Clinical", "Single-type", "SEIR(1)")


relabel_model = function(x, levels=model_levels){
  case_when(x == "SEIR1" ~ levels[7],
            x == "NegBin" ~ levels[1],
            x == "ZeroInf" ~ levels[4],
            x == "SEIR2" ~ levels[3],
            x == "Clinical" ~ levels[5],
            x == "SingleType" ~ levels[6],
            x == "TwoType" ~ levels[2],
            x == "Variable1" ~ levels[6],
            x == "Variable2" ~ levels[2]) %>% factor(levels)
}


scale_shape_model <- function(...){
  ggplot2:::manual_scale(
    'shape', 
    values = setNames(c(19,17,3,25,15,18,4), model_levels), 
    ...
  )
}


scale_col_model <- function(...){
  ggplot2:::manual_scale(
    'col', 
    values = setNames(c("#1B1919", "#00468B", "#A73030FF", "#925E9F", "#AD002A", "#008B45FF", "#0099B4"), model_levels), 
    ...
  )
}


custom_theme = theme(
  panel.background = element_rect(fill = "white", colour = "white"),
  axis.line = element_line(linewidth = 0.4, colour = "grey10"),
  text = element_text(size = 20,  family = "serif"),
  legend.key = element_rect(fill = "white", colour = "white"),
  strip.background = element_rect(fill = "white"),
  strip.text = element_text(size = 15, colour = 'black'),
  panel.spacing = unit(1, "lines"), plot.title = element_text(hjust = 0.5),
  legend.position = "bottom"
)


#### Plot types ####

plot_parameter = function(parm_values, parameters, ylabel="", models=c("Negative Binomial")){
parm_values %>%
  filter(Parameter %in% parameters, Model %in% models) %>%
  ggplot(aes(x = Label, col=Pathogen, shape=Model)) +
  geom_linerange(aes(ymin = X2.5., ymax=X97.5.),position=position_dodge2(width=0.75), lwd=2, alpha=0.3) +
  geom_linerange(aes(ymin = X25., ymax=X75.), position=position_dodge2(width=0.75), lwd=2, alpha=0.6) +
  geom_point(aes(y=X50.), position=position_dodge2(width=0.75), size=3) +
  xlab("") +
  ylab(ylabel) +
  custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = c(0.2,0.65)) +
  scale_col_pathogen() +
  scale_shape_model()
}


plot_extinction = function(chain, models=c("Negative Binomial")){
  chain %>%
    filter(Model %in% models) %>%
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
    annotate("text", x=5, y=0.17, size=6, label="Homogeneous limit", family="serif") +
    xlab("") +
    ylab("Extinction probability (q)") +
    scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = seq(0,1,0.1)) +
    custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "right") +
    scale_col_pathogen() +
    scale_shape_model()
}


plot_score = function(scores, models=c("Negative Binomial", "Two-type", "Clinical")){
  ggplot(scores %>% filter(Model %in% models), aes(x=Model, col=Pathogen, shape=Model)) +
    geom_point(aes(y=w), size=4) +
    facet_wrap(Label ~ ., nrow=5) +
    xlab("Model") +
    ylab("Akaike weight (w)") +
    custom_theme + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position="right") +
    scale_shape_model() +
    scale_col_pathogen()
}


plot_pairs = function(chain, par1, par2, model="", xlabel="", ylabel=""){
  chain %>%
    filter(Model == model) %>%
    ggplot(aes(x= .data[[par1]], y= .data[[par2]], col=Pathogen, shape=Model)) +
    geom_point(alpha=0.5) +
    # geom_density_2d(bins=4) +
    # geom_vline(xintercept=1, lty=2) +
    facet_wrap(Label ~ .) +
    ylab(ylabel) +
    xlab(xlabel) +
    custom_theme + scale_col_pathogen()
}


#### Model fit data ####

plot_fit = function(observed, predicted, scores, models=c("Negative Binomial", "Two-type", "Clinical"), labels=offspring_labels, Z_max=15, Z_inc=5, x_text=10.5){
  Z_range = c(0,Z_max)
  Zalt_labels = c(seq(0,Z_max-1), paste0("≥", Z_max))
  Zalt_ticks = rep("", length(Zalt_labels))
  Zalt_ticks[seq(1, length(Zalt_labels), by=Z_inc)] = Zalt_labels[seq(1, length(Zalt_labels), by=Z_inc)]
  summarize_df = expand.grid(Label = labels, Zalt = Zalt_labels)
  
  observed$Zalt = factor(cut(observed$Z, breaks = c(seq(Z_range[1] - 0.5, Z_range[2] - 0.5, by=1), Inf), labels=Zalt_labels), levels=Zalt_labels)
  predicted$Zalt = factor(cut(predicted$Z, breaks = c(seq(Z_range[1] - 0.5, Z_range[2] - 0.5, by=1), Inf), labels=Zalt_labels), levels=Zalt_labels)

  observed_accum = summarize_df %>%
    left_join(observed, by = c("Label", "Zalt")) %>%
    group_by(Label, Zalt, Pathogen) %>%
    summarize(count = sum(n), .groups = "drop") %>%
    replace_na(list(count = 0))
  
  predicted_accum = summarize_df %>%
    left_join(predicted, by = c("Label", "Zalt")) %>%
    group_by(Model, Label, Zalt, Pathogen) %>%
    summarize(count = sum(n), .groups = "drop") %>%
    replace_na(list(count = 0))
  
  max_off = observed %>%
    group_by(Dataset) %>%
    summarize(Max.n = max(n))
  
  text = scores %>%
    arrange(Dataset, Model) %>%
    filter(Model %in% models) %>%
    mutate(n = as.vector(sapply(max_off$Max.n, function(x) x * seq(0.9, 0.9 - 0.2*(length(models)-1), by=-0.2))),
           wr = paste0("w = ", round(w, 2)),
           aiccs = paste0("AICc = ", signif(AICc, 4)))
  text$Zalt = x_text
  
  plt = ggplot(predicted_accum %>% dplyr::filter(Model %in% models), aes(x=Zalt)) +
    geom_bar(data=observed_accum, aes(weight=count, fill=Pathogen), alpha=1.) + 
    geom_point(aes(y = count, col=Model, shape=Model), size=3, alpha=0.7) +
    geom_line(aes(y = count, group=Model, col=Model, lty=Model), lwd=1, alpha=0.7) +
    geom_text(data=text %>% filter(Model %in% models), aes(y=n, label=aiccs, col=Model), size=6, family="serif") +
    facet_wrap(Label ~ ., scales="free") +
    xlab("Secondary cases (Z)") +
    ylab("Frequency") +
    scale_x_discrete(labels = Zalt_ticks) +
    scale_y_continuous(label=scales::comma) +
    custom_theme + 
    scale_fill_pathogen() +
    scale_col_model() + 
    scale_shape_model()
  return(plt)
}
