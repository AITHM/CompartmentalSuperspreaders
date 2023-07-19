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


offspring_labels = c("Batam, 2020 (n=89)", "China, 2020 (n=1,178)","China, 2020 (n=2,048)","Hong Kong, 2020 (n=290)",  
                     "India, 2020 (n=88,527)", "Jakarta, 2020 (n=1,199)", "Sth. Korea, 2021 (n=344)", "Sth. Korea, 2021 (n=1,401)", 
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
  case_when(x == "Erlang" ~ levels[7],
            x == "NegBin" ~ levels[1],
            x == "ZeroInf" ~ levels[4],
            x == "Mixture" ~ levels[3],
            x == "Clinical" ~ levels[5],
            x == "SingleType" ~ levels[6],
            x == "TwoType" ~ levels[2],
            x == "Variable1" ~ levels[6],
            x == "Variable2" ~ levels[2]) %>% factor(levels)
}


scale_shape_model <- function(...){
  ggplot2:::manual_scale(
    'shape', 
    values = setNames(c(19,17,3,18,15,6,4), model_levels), 
    ...
  )
}


scale_col_model <- function(...){
  ggplot2:::manual_scale(
    'col', 
    values = setNames(c("#1B1919", "#00468B", "#42B540", "#ADB6B6", "#AD002A", "#0099B4", "#925E9F"), model_levels), 
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