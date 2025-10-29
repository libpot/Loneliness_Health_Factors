### SENSITIVITY ANALYSIS ###


### UNADJUSTED FOR DEPRESSIVE SYMPTOMS and SOCIAL FACTORS ###
# excluding depressive symptoms and social factors (no. of friends, living situation, relationship status and social support) from the set of time-varying confounders



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

# specifying a list of baseline and time-varying covariates
baseline_covariates_t0 <- "age + sex + ethnicity + parental_education + behavioral_monitoring + psych_overcontrol + cold_parenting + parental_alcoholuse + parental_smoking + chronic_condition +"
time_varying_covariates <- c("employment_", "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_")
time_varying_covariates_t0 <- paste(paste0(time_varying_covariates, "1"), collapse = " + ")
time_varying_covariates_t1 <- paste(paste0(time_varying_covariates, "2"), collapse = " + ")
time_varying_covariates_t2 <- paste(c( "education_3", paste0(time_varying_covariates, "3")), collapse = " + ")



# modelling treatment assignment #
exposure_models_imp_unadj <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5,
                                           ~ (.x - mean(.x))/sd(.x)), 
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), # standardizing loneliness scores (mean = 0, SD = 1)
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "gaussian", weights = weights, data = .x)),
         data_ext = map2(.x = data,
                         .y = exposure_model_t1,
                         ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data_ext = map2(.x = data_ext,
                         .y = exposure_model_t2,
                         ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data_ext = map2(.x = data_ext,
                         .y = exposure_model_t3,
                         ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))))


### Body Mass Index (BMI) ###

BMI_outcome_models_unadj <- exposure_models_imp_unadj %>% 
  mutate(outcome_model_l1_b2 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("bmi_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ loneliness_2")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l1_b2_1sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         effect_l1_b2_2sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), # marginal (counter-factual) contrasts
         outcome_model_l2_b3 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("bmi_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ loneliness_3")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l2_b3_1sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_l2_b3_2sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_l3_b4 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("bmi_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ loneliness_4")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l3_b4_1sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_l3_b4_2sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
BMI_outcome_effects_unadj <- select(BMI_outcome_models_unadj, starts_with("effect")) %>% 
  pivot_longer(cols = effect_l1_b2_1sd:effect_l3_b4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



### Alcohol Use ###

alcohol_outcome_models_unadj <- exposure_models_imp_unadj %>% 
  mutate(outcome_model_l1_b2 = map(.x = data_ext,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                                    "+ loneliness_2")), 
                                                  weights = weights, data = .x)),
         effect_l1_b2_1sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_l1_b2_2sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         outcome_model_l2_b3 = map(.x = data_ext,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                                    "+ loneliness_3")), 
                                                  weights = weights, data = .x)),
         effect_l2_b3_1sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_l2_b3_2sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_l3_b4 = map(.x = data_ext,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                                    "+ loneliness_4")), 
                                                  weights = weights, data = .x)),
         effect_l3_b4_1sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_l3_b4_2sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

alcohol_outcome_effects_unadj <- select(alcohol_outcome_models_unadj, starts_with("effect")) %>% 
  pivot_longer(cols = effect_l1_b2_1sd:effect_l3_b4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



### Smoking ###

smoking_outcome_models_unadj <- exposure_models_imp_unadj %>% 
  mutate(outcome_model_l1_b2 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("smoking_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ loneliness_2")), 
                                         family = "binomial", weights = weights, data = .x)),
         effect_l1_b2_1sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_l1_b2_2sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         outcome_model_l2_b3 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("smoking_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ loneliness_3")), 
                                         family = "binomial", weights = weights, data = .x)),
         effect_l2_b3_1sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_l2_b3_2sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_l3_b4 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("smoking_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ loneliness_4")), 
                                         family = "binomial", weights = weights, data = .x)),
         effect_l3_b4_1sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_l3_b4_2sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

smoking_outcome_effects_unadj <- select(smoking_outcome_models_unadj, starts_with("effect")) %>% 
  pivot_longer(cols = effect_l1_b2_1sd:effect_l3_b4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



### Physical activity ###

phys_act_outcome_models_unadj <- exposure_models_imp_unadj %>% 
  mutate(outcome_model_l1_b2 = map(.x = data_ext,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                                    "+ loneliness_2")), 
                                                  weights = weights, data = .x)),
         effect_l1_b2_1sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_l1_b2_2sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")),
         outcome_model_l2_b3 = map(.x = data_ext,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                                    "+ loneliness_3")), 
                                                  weights = weights, data = .x)),
         effect_l2_b3_1sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_l2_b3_2sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_l3_b4 = map(.x = data_ext,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                                    "+ loneliness_4")), 
                                                  weights = weights, data = .x)),
         effect_l3_b4_1sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_l3_b4_2sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

phys_act_outcome_effects_unadj <- select(phys_act_outcome_models_unadj, starts_with("effect")) %>% 
  pivot_longer(cols = effect_l1_b2_1sd:effect_l3_b4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 

 


### RESULTS ###
ATE_loneli_health_unadj <- bind_rows(BMI_outcome_effects_unadj, alcohol_outcome_effects_unadj, 
                                          smoking_outcome_effects_unadj, phys_act_outcome_effects_unadj) %>% 
  ungroup() %>% 
  mutate(contrast = relevel(factor(ifelse(contrast == "1 - -1", "-1 SD", "Mean")), ref = "Mean"),
         path = str_sub(effect, nchar(effect)-8, nchar(effect)-4),
         path = case_when(path == "l1_b2" ~ "adolescence \n (17y)",
                          path == "l2_b3" ~ "emerg. adulthood \n (22y)",
                          path == "l3_b4" ~ "young adulthood \n (28y)"),
         health_factor = c(rep("Body mass index", 6), rep("Alcohol use", 6), 
                           rep("Smoking", 6), rep("Physical Exercise", 6))) 

results_ATE_loneli_health_unadj <- ATE_loneli_health_unadj %>% 
  mutate(est = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                      " (", formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                      formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")")) %>% 
  select(health_factor, contrast, path, est) %>% 
  pivot_wider(names_from = "path",
              values_from = "est") %>% 
  arrange(health_factor, contrast)
write.csv(results_ATE_loneli_health_unadj, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_loneli_health_unadj")

plot_ATE_loneli_health_unadj <- ATE_loneli_health_unadj %>% 
  filter(contrast == "Mean") %>% 
  ungroup() %>% 
  mutate(health_factor = c(rep("**Body Mass Index**<br> (kg/m<sup>2</sup>)", 3), rep("**Alcohol Use**<br> (drinks per 4 weeks)", 3), 
                           rep("**Current Smoking**<br> (probability)", 3), rep("**Physical Exercise**<br> (minutes per week)", 3))) %>% 
  ggplot(aes(x = path, y = estimate)) +
  geom_point(size = 4, position = position_dodge(width = 0.3), color = "black") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ health_factor, ncol = 2, scales = "free") +
  labs(title = "Causal Effect of Loneliness on Health Outcomes",
       subtitle = "unadjusted for depressive symptoms and social factors",
       x = "",
       y = "") +
  theme_minimal() +
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_ATE_loneli_health_unadj.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 9, 
       height = 9,  
       bg="white",
       dpi=900)

