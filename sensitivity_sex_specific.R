### SENSITIVITY ANALYSIS ###

# EFFECT MODIFICATION: SEX-SPECIFIC ANLYSES #
   # separate models for women and men #

# loading packages
library(tidyverse) # data pre-processing
library(splines) # natural (restricted cubic) splines
library(marginaleffects) # marginal (average) treatment effects
library(sandwich) # robust (sandwich) standard errors
library(ggtext) # fine-tuning the plots
library(cowplot) # arranging multiple plots
library(mice) # pooling estimates across imputations

# loading multiply imputed dataset and the incomplete (original) dataset
UiN_data_weights <- read.csv("N:/durable/Data_analyses/Libor/Social_Gradient_Mental_Health_Lifespan/data/UiN_data_IPAW")[, c("id", "weights")]
imp_UiN_data <- read.csv("N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Data/imp_UiN_data") %>% 
  merge(., UiN_data_weights, by = "id")
UiN_data_analysis <- read.csv("N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Data/UiN_data_analysis") %>% 
  merge(., UiN_data_weights, by = "id", no.dups = F)


# specifying a list of baseline and time-varying covariates
baseline_covariates_t0 <- "age + ethnicity + parental_education + behavioral_monitoring + psych_overcontrol + cold_parenting + parental_alcoholuse + parental_smoking + chronic_condition +"
time_varying_covariates <- c("employment_", "living_situation_", "relationship_", "friends_",
                             "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                             "depression_", "social_support_")
time_varying_covariates_t0 <- paste(paste0(time_varying_covariates, "1"), collapse = " + ")
time_varying_covariates_t1 <- paste(paste0(time_varying_covariates, "2"), collapse = " + ")
time_varying_covariates_t2 <- paste(paste0(time_varying_covariates, "2"), collapse = " + ")


# dropping one time-varying covariate due to convergence issues
time_varying_covariates_adapt <- c("employment_", "living_situation_", "relationship_", "friends_",
                                 "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                                 "depression_")
time_varying_covariates_t2_adapt <- paste(paste0(time_varying_covariates_adapt, "3"), collapse = " + ") 






####################################
### LONELINESS -> HEALTH FACTORS ###
####################################

### EXPOSURE MODELS FOR LONELINESS ###
# modelling treatment assignment/selection mechanisms #

exposure_models_imp_sex <- imp_UiN_data %>% 
  group_by(.imp, sex) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5,
                                           ~ (.x - mean(.x))/sd(.x)),   # standardizing loneliness scores (mean = 0, SD = 1)
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))),
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "gaussian", data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "gaussian", data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "gaussian", data = .x)),
         data = map2(.x = data,
                     .y = exposure_model_t1,
                     ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t2,
                     ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t3,
                     ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))))




### OUTCOME MODELS FOR HEALTH FACTORS ###
# estimating causal effect of loneliness #

### Body Mass Index (BMI) ###

BMI_outcome_models_sex <- exposure_models_imp_sex %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("bmi_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ loneliness_2")), 
                                         family = "gaussian", data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), # Average treatment effect (ATE) using g-computation for pre-specified marginal contrasts
         outcome_model_t2_t3 = map(.x = data,
                                   ~ glm(as.formula(paste0("bmi_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ loneliness_3")), 
                                         family = "gaussian", data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("bmi_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2_adapt, 
                                                           "+ loneliness_4")), 
                                         family = "gaussian", data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average) treatment effects of 1SD loneliness increase compared to mean and -1SD
BMI_outcome_effects_sex <- select(BMI_outcome_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`) 





### Alcohol Use ###

alcohol_outcome_models_sex <- exposure_models_imp_sex %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                                    "+ loneliness_2")), 
                                                  data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                                    "+ loneliness_3")), 
                                                  data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2_adapt, 
                                                                 "+ loneliness_4")), 
                                                  data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

alcohol_outcome_effects_sex <- select(alcohol_outcome_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling estimates according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`) 




### Smoking ###

smoking_outcome_models_sex <- exposure_models_imp_sex %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("smoking_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ loneliness_2")), 
                                         family = "binomial", data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,
                                   ~ glm(as.formula(paste0("smoking_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ loneliness_3")), 
                                         family = "binomial", data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("smoking_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2_adapt, 
                                                           "+ loneliness_4")), 
                                         family = "binomial", data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

smoking_outcome_effects_sex <- select(smoking_outcome_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`)





### Physical activity ###

phys_act_outcome_models_sex <- exposure_models_imp_sex %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                                    "+ loneliness_2")), 
                                                  data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")),
         outcome_model_t2_t3 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                                    "+ loneliness_3")), 
                                                  data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2_adapt, 
                                                                    "+ loneliness_4")), 
                                                  data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

phys_act_outcome_effects_sex <- select(phys_act_outcome_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`) 





### MAIN PLOTS ###

# CAUSAL EFFECT OF LONELINESS ON HEALTH OUTCOMES #
ATE_loneli_health_sex <- bind_rows(BMI_outcome_effects_sex, alcohol_outcome_effects_sex, 
                                   smoking_outcome_effects_sex, phys_act_outcome_effects_sex) %>% 
  ungroup() %>% 
  filter(contrast != "2 - 0") %>% 
  mutate(contrast = relevel(factor(ifelse(contrast == "1 - -1", "-1 SD", "Mean")), ref = "Mean"),
         path = str_sub(effect, nchar(effect)-8, nchar(effect)-4),
         path = case_when(path == "t1_t2" ~ "adolescence \n (17y)",
                          path == "t2_t3" ~ "emerg. adulthood \n (22y)",
                          path == "t3_t4" ~ "young adulthood \n (28y)"),
         health_factor = c(rep("Body mass index", 12), rep("Alcohol use", 12), 
                           rep("Smoking", 12), rep("Physical Exercise", 12)),
         gender = ifelse(sex == "female", "women", "men"))

results_ATE_loneli_health_sex <- ATE_loneli_health_sex %>% 
  mutate(est = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                      " (", formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                      formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")")) %>% 
  select(health_factor, gender, contrast, path, est) %>% 
  pivot_wider(names_from = "path",
              values_from = "est") %>% 
  arrange(health_factor, gender)
write.csv(results_ATE_loneli_health_sex, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_loneli_health_sex")

plot_ATE_loneli_health_sex <- ATE_loneli_health_sex %>% 
  mutate(health_factor = c(rep("**Body Mass Index**<br> (kg/m<sup>2</sup>)", 12), rep("**Alcohol Use**<br> (drinks per 4 weeks)", 12), 
                           rep("**Current Smoking**<br> (probability)", 12), rep("**Physical Exercise**<br> (minutes per week)", 12))) %>% 
  filter(contrast == "Mean") %>% 
  ggplot(aes(x = path, y = estimate, color = gender)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ health_factor, 
             ncol = 2, scales = "free") +
  labs(title = "Causal Effect of Loneliness on Health Outcomes by Gender",
       x = "",
       y = "",
       color = "Gender") +
  theme_minimal() +
  scale_color_manual(values = c("men" = "black", "women" = "darkgrey")) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_ATE_loneli_health_sex.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 9, 
       height = 9,  
       bg="white",
       dpi=900)









####################################
### HEALTH FACTORS -> LONELINESS ###
####################################

# in models for BMI, social support was unlike in other models included also in the outcome model, resulting in singularity problem (due to perfect collinearity with propensity score),
# here, we are explicitly dropping social support 
time_varying_covariates_bmi <- c("employment_", "living_situation_", "relationship_", "friends_",
                                 "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                                 "depression_")
time_varying_covariates_t0_bmi <- paste(paste0(time_varying_covariates_bmi, "1"), collapse = " + ")
time_varying_covariates_t1_bmi <- paste(paste0(time_varying_covariates_bmi, "2"), collapse = " + ")


### Body Mass Index (BMI) ###

BMI_loneli_models_sex <- imp_UiN_data %>% 
  group_by(.imp, sex) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "gaussian",  data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "gaussian",  data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "gaussian",  data = .x)),
         data = map2(.x = data,
                     .y = exposure_model_t1,
                     ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t2,
                     ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t3,
                     ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))),
         outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0_bmi, 
                                                           "+ ns(bmi_2, knots = quantile(bmi_2, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_2, probs = c(0.05, 0.95)))")),
                                         family = "gaussian",  data = .x)),
         effect_t1_t2_plus = map(.x = outcome_model_t1_t2,
                                 ~ avg_comparisons(.x, variables = list(bmi_2 = c(21.75,30)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         effect_t1_t2_minus = map(.x = outcome_model_t1_t2,
                                  ~ avg_comparisons(.x, variables = list(bmi_2 = c(21.75,17)), vcov = "HC3")), # marginal (counter-factual) contrasts
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1_bmi, 
                                                           "+ ns(bmi_3, knots = quantile(bmi_3, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_3, probs = c(0.05, 0.95)))")), 
                                         family = "gaussian",  data = .x)),
         effect_t2_t3_plus = map(.x = outcome_model_t2_t3,
                                 ~ avg_comparisons(.x, variables = list(bmi_3 = c(21.75,30)), vcov = "HC3")),
         effect_t2_t3_minus = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(bmi_3 = c(21.75,17)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(bmi_4, knots = quantile(bmi_4, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_4, probs = c(0.05, 0.95)))")), 
                                         family = "gaussian",  data = .x)),
         effect_t3_t4_plus = map(.x = outcome_model_t3_t4,
                                 ~ avg_comparisons(.x, variables = list(bmi_4 = c(21.75,30)), vcov = "HC3")),
         effect_t3_t4_minus = map(.x = outcome_model_t3_t4,
                                  ~ avg_comparisons(.x, variables = list(bmi_4 = c(21.75,17)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
BMI_loneli_effects_sex <- select(BMI_loneli_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_plus:effect_t3_t4_minus, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`) 





### Current smoking ###

smoke_loneli_models_sex <- imp_UiN_data %>% 
  group_by(.imp, sex) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_1:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("smoking_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "binomial", data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("smoking_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "binomial", data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("smoking_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "binomial", data = .x)),
         data = map2(.x = data,
                     .y = exposure_model_t1,
                     ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t2,
                     ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t3,
                     ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))),
         outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ smoking_2")),
                                         family = "gaussian", data = .x)),
         effect_t1_t2 = map(.x = outcome_model_t1_t2,
                            ~ avg_comparisons(.x, variables = list(smoking_2 = c(0,1)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ smoking_3")), 
                                         family = "gaussian", data = .x)),
         effect_t2_t3 = map(.x = outcome_model_t2_t3,
                            ~ avg_comparisons(.x, variables = list(smoking_3 = c(0,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ smoking_4")), 
                                         family = "gaussian", data = .x)),
         effect_t3_t4 = map(.x = outcome_model_t3_t4,
                            ~ avg_comparisons(.x, variables = list(smoking_4 = c(0,1)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
smoke_loneli_effects_sex <- select(smoke_loneli_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2:effect_t3_t4, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`) 





### Alcohol Use ###

alcohol_loneli_models_sex <- imp_UiN_data %>% 
  group_by(.imp, sex) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("alcohol_use_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "poisson", data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("alcohol_use_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "poisson", data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("alcohol_use_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "poisson", data = .x)),
         data = map2(.x = data,
                     .y = exposure_model_t1,
                     ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t2,
                     ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t3,
                     ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))),
         outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ ns(alcohol_use_2, knots = quantile(alcohol_use_2, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(alcohol_use_2, probs = c(0.00, 0.95)))")),
                                         family = "gaussian",  data = .x)),
         effect_t1_t2_male_heavy = map(.x = outcome_model_t1_t2,
                                       ~ avg_comparisons(.x, variables = list(alcohol_use_2 = c(10,60)), vcov = "HC3")), # if male, heavy use >= 60,
         effect_t1_t2_female_heavy = map(.x = outcome_model_t1_t2,
                                         ~ avg_comparisons(.x, variables = list(alcohol_use_2 = c(10,32)), vcov = "HC3")), # if female, heavy use >= 32
         effect_t1_t2_absti = map(.x = outcome_model_t1_t2,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_2 = c(10,0)), vcov = "HC3")), # marginal (counter-factual) contrasts
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(alcohol_use_3, knots = quantile(alcohol_use_3, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(alcohol_use_3, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian",  data = .x)),
         effect_t2_t3_male_heavy = map(.x = outcome_model_t2_t3,
                                      ~ avg_comparisons(.x, variables = list(alcohol_use_3 = c(10,60)), vcov = "HC3")),
         effect_t2_t3_female_heavy = map(.x = outcome_model_t2_t3,
                                       ~ avg_comparisons(.x, variables = list(alcohol_use_3 = c(10,32)), vcov = "HC3")), 
         effect_t2_t3_absti = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_3 = c(10,0)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(alcohol_use_4, knots = quantile(alcohol_use_4, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(alcohol_use_4, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian",  data = .x)),
         effect_t3_t4_male_heavy = map(.x = outcome_model_t3_t4,
                                       ~ avg_comparisons(.x, variables = list(alcohol_use_4 = c(10,60)), vcov = "HC3")),
         effect_t3_t4_female_heavy = map(.x = outcome_model_t3_t4,
                                         ~ avg_comparisons(.x, variables = list(alcohol_use_4 = c(10,32)), vcov = "HC3")),  
         effect_t3_t4_absti = map(.x = outcome_model_t3_t4,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_4 = c(10,0)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
alcohol_loneli_effects_sex <- select(alcohol_loneli_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_male_heavy:effect_t3_t4_absti, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`) %>% 
  filter(!(sex == "female" & contrast == "60 - 10"),
         !(sex == "male" & contrast == "32 - 10"))




### Physical Activity ###

phys_exer_loneli_models_sex <- imp_UiN_data %>% 
  group_by(.imp, sex) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ MASS::glm.nb(as.formula(paste0("phys_exer_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ MASS::glm.nb(as.formula(paste0("phys_exer_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ MASS::glm.nb(as.formula(paste0("phys_exer_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), weights = weights, data = .x)),
         data = map2(.x = data,
                     .y = exposure_model_t1,
                     ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t2,
                     ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t3,
                     ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))),
         outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ ns(phys_exer_2, knots = quantile(phys_exer_2, probs = 0.50), Boundary.knots = quantile(phys_exer_2, probs = c(0.00, 0.95)))")),
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t1_t2 = map(.x = outcome_model_t1_t2,
                            ~ avg_comparisons(.x, variables = list(phys_exer_2 = c(0,150)), vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(phys_exer_3, knots = quantile(phys_exer_3, probs = 0.50), Boundary.knots = quantile(phys_exer_3, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3 = map(.x = outcome_model_t2_t3,
                            ~ avg_comparisons(.x, variables = list(phys_exer_3 = c(0, 150)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(phys_exer_4, knots = quantile(phys_exer_4, probs = 0.50), Boundary.knots = quantile(phys_exer_4, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4 = map(.x = outcome_model_t3_t4,
                            ~ avg_comparisons(.x, variables = list(phys_exer_4 = c(0, 150)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
phys_exer_loneli_effects_sex <- select(phys_exer_loneli_models_sex, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2:effect_t3_t4, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect, sex) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, sex, contrast, estimate, `2.5 %`, `97.5 %`) 




### MAIN PLOTS ###

# CAUSAL EFFECT OF HEALTH FACTORS ON LONELINESS #
ATE_health_loneli_sex <- bind_rows(BMI_loneli_effects_sex, alcohol_loneli_effects_sex, 
                                   smoke_loneli_effects_sex, phys_exer_loneli_effects_sex) %>% 
  ungroup() %>% 
  mutate(path = str_sub(effect,8,12),
         path = case_when(path == "t1_t2" ~ "adolescence \n (17y)",
                          path == "t2_t3" ~ "emerg. adulthood \n (22y)",
                          path == "t3_t4" ~ "young adulthood \n (28y)"),
         contrast = case_when(contrast == "30 - 21.75" ~ "adiposity",
                              contrast == "21.75 - 17" ~ "underweight",
                              contrast %in% c("60 - 10", "32 - 10") ~ "heavy",
                              contrast == "10 - 0" ~ "abstinence",
                              contrast == "1 - 0" ~ "smoking",
                              contrast %in% c("3 - 0", "150 - 0") ~ "moderate"),
         gender = ifelse(sex == "female", "women", "men"),
         across(c(estimate, `2.5 %`, `97.5 %`), 
                ~ case_when(contrast %in% c("abstinence", "underweight") ~ .x * -1,
                            TRUE ~ .x)))

results_ATE_health_loneli_sex <- ATE_health_loneli_sex %>% 
  mutate(est = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                      " (", formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                      formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")"),
         factor = case_when(contrast %in% c("abstinence", "heavy") ~ "Alcohol use",
                            contrast %in% c("underweight", "adiposity") ~ "Body mass index",
                            contrast == "smoking" ~ "Smoking",
                            contrast == "moderate" ~ "Physical activity")) %>% 
  select(factor, gender, contrast, path, est) %>% 
  pivot_wider(names_from = "path",
              values_from = "est") %>% 
  arrange(factor, contrast)
write.csv(results_ATE_health_loneli_sex, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_health_loneli_sex")



ATE_bmi_loneli_plot_sex <- ATE_health_loneli_sex %>% 
  filter(contrast %in% c("adiposity", "underweight")) %>% 
  mutate(contrast = factor(contrast, levels = c("underweight", "adiposity"))) %>% 
  ggplot(aes(x = path, y = estimate, color = gender, shape = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Body Mass Index",
       x = "",
       y = "",
       shape = "ref.: normative",
       color = "Gender") +
  theme_minimal() +
  scale_color_manual(values = c("men" = "black", "women" = "darkgrey")) +
  scale_shape_manual(values = c("adiposity" = 16, "underweight" = 17)) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(t = 40),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_alcohol_loneli_plot_sex <- ATE_health_loneli_sex %>% 
  filter(contrast %in% c("heavy", "abstinence")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = gender)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Alcohol Use",
       x = "",
       y = "",
       shape = "ref.: occasional/low",
       color = "Gender") +
  theme_minimal() +
  scale_color_manual(values = c("men" = "black", "women" = "darkgrey")) + 
  scale_shape_manual(values = c("heavy" = 16, "abstinence" = 17)) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(t = 40, l = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black")) +
  guides(color = guide_legend(order = 2),  
         shape = guide_legend(order = 1))

ATE_smoke_loneli_plot_sex <- ATE_health_loneli_sex %>% 
  filter(contrast %in% c("smoking")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = gender)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Current smoking",
       x = "",
       y = "",
       color = "Gender",
       shape = "ref.: no smoking") +
  theme_minimal() +
  scale_color_manual(values = c("men" = "black", "women" = "darkgrey")) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(l = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_exer_loneli_plot_sex <- ATE_health_loneli_sex %>% 
  filter(contrast %in% c("moderate")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = gender)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Physical Exercise",
       x = "",
       y = "",
       color = "Gender",
       shape = "ref.: inactive") +
  theme_minimal() +
  scale_color_manual(values = c("men" = "black", "women" = "darkgrey")) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

plot_ATE_health_loneli_sex <- cowplot::plot_grid(ATE_alcohol_loneli_plot_sex, ATE_bmi_loneli_plot_sex, 
                                                 ATE_smoke_loneli_plot_sex, ATE_exer_loneli_plot_sex, 
                                                 ncol = 2) +
  draw_label("Causal Effect of Health Factors on Loneliness by Gender", fontface = 'bold', fontfamily = "Times", size = 16, hjust = 0.5, y = 0.98) +
  draw_label("Loneliness (standardized)", angle = 90, x = 0.03, y = 0.5, vjust = 0.5, fontfamily = "Times", size = 12)
#draw_label("Prospective path", y = 0.03, x = 0.5, hjust = 0.5, fontfamily = "Times", size = 14)

ggsave(filename = "plot_ATE_health_loneli_sex.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 12, 
       height = 9,  
       bg="white",
       dpi=900)

