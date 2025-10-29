######################################################################
### THE DYNAMIC INTERPLAY BETWEEN LONELINESS AND HEALTH BEHAVIORS ###
#####################################################################

              ###############################
              ### STRUCTURAL MEAN MODELS ###
              ##############################



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
time_varying_covariates <- c("employment_", "living_situation_", "relationship_", "friends_",
                             "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                             "depression_", "social_support_")
time_varying_covariates_t0 <- paste(paste0(time_varying_covariates, "1"), collapse = " + ")
time_varying_covariates_t1 <- paste(paste0(time_varying_covariates, "2"), collapse = " + ")
time_varying_covariates_t2 <- paste(c( "education_3", paste0(time_varying_covariates, "3")), collapse = " + ")




####################################
### LONELINESS -> HEALTH FACTORS ###
####################################

### EXPOSURE MODELS FOR LONELINESS ###
# modelling treatment assignment/selection mechanisms #

exposure_models_imp <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                          ~ .x %>% mutate(across(loneliness_1:loneliness_5,
                                           ~ (.x - mean(.x))/sd(.x)),   # standardizing loneliness scores (mean = 0, SD = 1)
                                    across(smoking_3:smoking_5,
                                          ~ ifelse(.x == "smoking", 1, 0)))),
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "gaussian", weights = weights, data = .x)),
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

BMI_outcome_models <- exposure_models_imp %>% 
 mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("bmi_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                         "+ loneliness_2")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                 ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), # Average treatment effect (ATE) using g-computation for pre-specified marginal contrasts
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(loneliness_2 = c(-1,0,1)), by = "loneliness_2", vcov = "HC3")), # marginal (counter-factual) predictions
         outcome_model_t2_t3 = map(.x = data,
                                   ~ glm(as.formula(paste0("bmi_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                         "+ loneliness_3")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(loneliness_3 = c(-1,0,1)), by = "loneliness_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("bmi_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                         "+ loneliness_4")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                              ~ avg_predictions(.x, variables = list(loneliness_4 = c(-1,0,1)), by = "loneliness_4", vcov = "HC3")))
         



# extracting and pooling estimates across imputed datasets: marginal (average) treatment effects of 1SD loneliness increase compared to mean and -1SD
BMI_outcome_effects <- select(BMI_outcome_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



# extracting and pooling estimates across imputed datasets: average predicted values for -1SD, mean and +1SD of loneliness
BMI_outcome_avg_pred <- select(BMI_outcome_models, starts_with("avg_pred")) %>% 
  pivot_longer(cols = avg_pred_t1_t2:avg_pred_t3_t4, names_to = "name", values_to = "means") %>% 
  group_by(name) %>%  nest() %>% 
  mutate(data = map(.x = data, ~ .x %>% mutate(minus = list(.x$means[[1]][1,]))),
         data = map(.x = data, ~ .x %>% mutate(mean = list(.x$means[[1]][2,]))),
         data = map(.x = data, ~ .x %>% mutate(plus = list(.x$means[[1]][3,]))),
         minus = map(.x = data, ~ summary(pool(.x$minus), conf.int = TRUE)),
         mean = map(.x = data, ~ summary(pool(.x$mean), conf.int = TRUE)),
         plus = map(.x = data, ~ summary(pool(.x$plus), conf.int = TRUE))) %>%  
  select(minus, mean, plus) %>% pivot_longer(cols = minus:plus, names_to = "level", values_to = "value") %>% unnest() %>% 
  select(name, level, estimate, `2.5 %`, `97.5 %`) 

# Model-based predictions of health factor (BMI) for different levels of loneliness conditional on average covariate values
# (i.e., predictions for hypothetical representative individual whose personal characteristics are exactly average (numeric) or modal (categorical))
cond_predictions_bmi2 <- as.data.frame(predictions(BMI_outcome_models$outcome_model_t1_t2[[1]], variables = list(loneliness_2 = seq(from = -1.1, to = 1.1, by = 0.1)), newdata = "mean")) %>% 
                              mutate(loneliness = loneliness_2, effect = "T1 on T2")
cond_predictions_bmi3 <- as.data.frame(predictions(BMI_outcome_models$outcome_model_t2_t3[[1]], variables = list(loneliness_3 = seq(from = -1.1, to = 2.1, by = 0.1)), newdata = "mean")) %>% 
                              mutate(loneliness = loneliness_3, effect = "T2 on T3")
cond_predictions_bmi4 <- as.data.frame(predictions(BMI_outcome_models$outcome_model_t3_t4[[1]], variables = list(loneliness_4 = seq(from = -1.1, to = 2.1, by = 0.1)), newdata = "mean")) %>% 
                              mutate(loneliness = loneliness_4, effect = "T3 on T4")

windowsFonts(Times = windowsFont("Times New Roman"))

non_linear_loneli_bmi <- bind_rows(cond_predictions_bmi2, cond_predictions_bmi3, cond_predictions_bmi4) %>% 
  select(estimate, conf.low, conf.high, loneliness, effect) %>% 
  ggplot(aes(x = loneliness, y = estimate, 
             group = effect, fill = effect, 
             ymin = conf.low, ymax = conf.high)) + 
  geom_line(aes(color = effect), size = 1) + 
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(limits = c(-1.1,1.1)) + 
  labs(title = "Predicted Values of Body Mass Index",
       x = "Loneliness (standardized)", 
       y = "Body Mass Index (BMI)",
       color = "Prospective Path",
       fill = "Prospective Path") +
  theme_minimal() +
  theme(text = element_text(family = "Times"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.background = element_blank())

ggsave(filename = "non_linear_loneli_bmi.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 6, 
       height = 4,  
       bg="white",
       dpi=700)




### Alcohol Use ###

alcohol_outcome_models <- exposure_models_imp %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                  "+ loneliness_2")), 
                                                  data = .x, weights = weights)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(loneliness_2 = c(-1,0,1)), by = "loneliness_2", vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                  "+ loneliness_3")), 
                                                  data = .x, weights = weights)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(loneliness_3 = c(-1,0,1)), by = "loneliness_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("alcohol_use_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                 "+ loneliness_4")), 
                                                 data = .x, weights = weights)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                               ~ avg_predictions(.x, variables = list(loneliness_4 = c(-1,0,1)), by = "loneliness_4", vcov = "HC3")))

alcohol_outcome_effects <- select(alcohol_outcome_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling estimates according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



### Smoking ###

smoking_outcome_models <- exposure_models_imp %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("smoking_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                         "+ loneliness_2")), 
                                         family = "binomial", weights = weights, data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(loneliness_2 = c(-1,0,1)), by = "loneliness_2", vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,
                                   ~ glm(as.formula(paste0("smoking_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                         "+ loneliness_3")), 
                                         family = "binomial", weights = weights, data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(loneliness_3 = c(-1,0,1)), by = "loneliness_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("smoking_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                         "+ loneliness_4")), 
                                         family = "binomial", weights = weights, data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                              ~ avg_predictions(.x, variables = list(loneliness_4 = c(-1,0,1)), by = "loneliness_4", vcov = "HC3")))

smoking_outcome_effects <- select(smoking_outcome_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



### Physical activity ###

phys_act_outcome_models <- exposure_models_imp %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                  "+ loneliness_2")), 
                                                  data = .x, weights = weights)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_t1_t2_2sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(loneliness_2 = c(-1,0,1)), by = "loneliness_2", vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                  "+ loneliness_3")), 
                                                  data = .x, weights = weights)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_t2_t3_2sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(loneliness_3 = c(-1,0,1)), by = "loneliness_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ MASS::glm.nb(as.formula(paste0("phys_exer_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                  "+ loneliness_4")), 
                                                   data = .x, weights = weights)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_t3_t4_2sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                              ~ avg_predictions(.x, variables = list(loneliness_4 = c(-1,0,1)), by = "loneliness_4", vcov = "HC3")))

phys_act_outcome_effects <- select(phys_act_outcome_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 




### MAIN PLOTS ###

windowsFonts(Times = windowsFont("Times New Roman"))


# CAUSAL EFFECT OF LONELINESS ON HEALTH OUTCOMES #
ATE_loneli_health <- bind_rows(BMI_outcome_effects, alcohol_outcome_effects, 
                               smoking_outcome_effects, phys_act_outcome_effects) %>% 
  ungroup() %>% 
  filter(contrast != "2 - 0") %>% 
  mutate(contrast = relevel(factor(ifelse(contrast == "1 - -1", "-1 SD", "Mean")), ref = "Mean"),
         path = str_sub(effect, nchar(effect)-8, nchar(effect)-4),
         path = case_when(path == "t1_t2" ~ "adolescence \n (17y)",
                          path == "t2_t3" ~ "emerg. adulthood \n (22y)",
                          path == "t3_t4" ~ "young adulthood \n (28y)"),
         health_factor = c(rep("Body mass index", 6), rep("Alcohol use", 6), 
                           rep("Smoking", 6), rep("Physical Exercise", 6))) 
write.csv(ATE_loneli_health, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_loneli_health_raw")
ATE_loneli_health <- read_csv("N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_loneli_health_raw")


results_ATE_loneli_health <- ATE_loneli_health %>% 
mutate(est = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
             " (", formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
             formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")")) %>% 
  select(health_factor, contrast, path, est) %>% 
  pivot_wider(names_from = "path",
              values_from = "est") %>% 
  arrange(health_factor, contrast)
write.csv(results_ATE_loneli_health, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_loneli_health")


plot_ATE_loneli_bmi <- ATE_loneli_health %>% 
  filter(contrast == "Mean", health_factor == "Body mass index") %>% 
  ggplot(aes(x = path, y = estimate)) +
  geom_point(size = 5, position = position_dodge(width = 0.3), color = "#5773CCFF") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.2, size = 1.5, color = "#5773CCFF") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Body mass index",
       x = "",
       y = expression(kg/m^2)) +
  theme_minimal() +
  theme(strip.text = element_text(size = 15, family = "Times"),
        axis.text = element_text(size = 15, family = "Times", color = "black"),
        text = element_text(size = 15, family = "Times"),
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15, family = "Times"),
        plot.margin = margin(t = 10, b = 10, l = 20, r = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

plot_ATE_loneli_alc <- ATE_loneli_health %>% 
  filter(contrast == "Mean", health_factor == "Alcohol use") %>% 
  ggplot(aes(x = path, y = estimate)) +
  geom_point(size = 5, position = position_dodge(width = 0.3), color = "#5773CCFF") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.2, size = 1.5, color = "#5773CCFF") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Alcohol use",
       x = "",
       y = "drinks per 4 weeks") +
  theme_minimal() +
  theme(strip.text = element_text(size = 15, family = "Times"),
        axis.text = element_text(size = 15, family = "Times", color = "black"),
        text = element_text(size = 15, family = "Times"),
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15, family = "Times"),
        plot.margin = margin(t = 10, b = 10, l = 20, r = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

plot_ATE_loneli_smoke <- ATE_loneli_health %>% 
  filter(contrast == "Mean", health_factor == "Smoking") %>% 
  ggplot(aes(x = path, y = estimate)) +
  geom_point(size = 5, position = position_dodge(width = 0.3), color = "#5773CCFF") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.2, size = 1.5, color = "#5773CCFF") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Smoking",
       x = "",
       y = "probability") +
  theme_minimal() +
  theme(strip.text = element_text(size = 15, family = "Times"),
        axis.text = element_text(size = 15, family = "Times", color = "black"),
        text = element_text(size = 15, family = "Times"),
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15, family = "Times"),
        plot.margin = margin(l = 20, t = 20, r = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))


plot_ATE_loneli_exer <- ATE_loneli_health %>% 
  filter(contrast == "Mean", health_factor == "Physical Exercise") %>% 
  ggplot(aes(x = path, y = estimate)) +
  geom_point(size = 5, position = position_dodge(width = 0.3), color = "#5773CCFF") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.2, size = 1.5, color = "#5773CCFF") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Physical Exercise",
       x = "",
       y = "minutes per week") +
  theme_minimal() +
  theme(strip.text = element_text(size = 15, family = "Times"),
        axis.text = element_text(size = 15, family = "Times", color = "black"),
        text = element_text(size = 15, family = "Times"),
        plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15, family = "Times"),
        plot.margin = margin(l = 20, t = 20, r = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

plot_ATE_loneli_health <- cowplot::plot_grid(plot_ATE_loneli_alc, plot_ATE_loneli_bmi, 
                                             plot_ATE_loneli_smoke, plot_ATE_loneli_exer, 
                                             ncol = 2)
ggsave(filename = "plot_ATE_loneli_health.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 14, 
       height = 12,  
       bg="white",
       dpi=1050)










####################################
### HEALTH FACTORS -> LONELINESS ###
####################################

# in models for BMI, social support was unlike in other models included also in the outcome model, resulting in singularity problem (due to perfect collinearity with propensity score),
# here, we are explicitly dropping social support 
time_varying_covariates_bmi <- c("employment_", "living_situation_", "relationship_", "friends_",
                                 "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                                 "depression_")
time_varying_covariates_t0_bmi <- paste(paste0(time_varying_covariates_bmi, "1"), collapse = " + ")

### Body Mass Index (BMI) ###

BMI_loneli_models <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "gaussian", weights = weights, data = .x)),
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
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t1_t2_plus = map(.x = outcome_model_t1_t2,
                                 ~ avg_comparisons(.x, variables = list(bmi_2 = c(21.75,30)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         effect_t1_t2_minus = map(.x = outcome_model_t1_t2,
                                  ~ avg_comparisons(.x, variables = list(bmi_2 = c(21.75,17)), vcov = "HC3")), # marginal (counter-factual) contrasts
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(bmi_2 = c(17,21.75,30)), by = "bmi_2", vcov = "HC3")), # marginal (counter-factual) predictions
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(bmi_3, knots = quantile(bmi_3, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_3, probs = c(0.05, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_plus = map(.x = outcome_model_t2_t3,
                                 ~ avg_comparisons(.x, variables = list(bmi_3 = c(21.75,30)), vcov = "HC3")),
         effect_t2_t3_minus = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(bmi_3 = c(21.75,17)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(bmi_3 = c(17,21.75,30)), by = "bmi_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(bmi_4, knots = quantile(bmi_4, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_4, probs = c(0.05, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_plus = map(.x = outcome_model_t3_t4,
                                 ~ avg_comparisons(.x, variables = list(bmi_4 = c(21.75,30)), vcov = "HC3")),
         effect_t3_t4_minus = map(.x = outcome_model_t3_t4,
                                  ~ avg_comparisons(.x, variables = list(bmi_4 = c(21.75,17)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                              ~ avg_predictions(.x, variables = list(bmi_4 = c(17,21.75,30)), by = "bmi_4", vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
BMI_loneli_effects <- select(BMI_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_plus:effect_t3_t4_minus, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



### Current smoking ###

smoke_loneli_models <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_1:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("smoking_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "binomial", weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("smoking_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "binomial", weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("smoking_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "binomial", weights = weights, data = .x)),
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
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t1_t2 = map(.x = outcome_model_t1_t2,
                                 ~ avg_comparisons(.x, variables = list(smoking_2 = c(0,1)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(smoking_2 = c(0,1)), by = "smoking_2", vcov = "HC3")), # marginal (counter-factual) predictions
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ smoking_3")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3 = map(.x = outcome_model_t2_t3,
                                 ~ avg_comparisons(.x, variables = list(smoking_3 = c(0,1)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(smoking_3 = c(0,1)), by = "smoking_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ smoking_4")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4 = map(.x = outcome_model_t3_t4,
                                 ~ avg_comparisons(.x, variables = list(smoking_4 = c(0,1)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                              ~ avg_predictions(.x, variables = list(smoking_4 = c(0,1)), by = "smoking_4", vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
smoke_loneli_effects <- select(smoke_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2:effect_t3_t4, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 



### Alcohol Use ###

alcohol_loneli_models <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ MASS::glm.nb(as.formula(paste0("alcohol_use_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ MASS::glm.nb(as.formula(paste0("alcohol_use_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ MASS::glm.nb(as.formula(paste0("alcohol_use_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), weights = weights, data = .x)),
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
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t1_t2_heavy = map(.x = outcome_model_t1_t2,
                                 ~ avg_comparisons(.x, variables = list(alcohol_use_2 = c(10,60)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         effect_t1_t2_absti = map(.x = outcome_model_t1_t2,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_2 = c(10,0)), vcov = "HC3")), # marginal (counter-factual) contrasts
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(alcohol_use_2 = c(0,10,60)), by = "alcohol_use_2", vcov = "HC3")), # marginal (counter-factual) predictions
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(alcohol_use_3, knots = quantile(alcohol_use_3, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(alcohol_use_3, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_heavy = map(.x = outcome_model_t2_t3,
                                 ~ avg_comparisons(.x, variables = list(alcohol_use_3 = c(10,60)), vcov = "HC3")),
         effect_t2_t3_absti = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_3 = c(10,0)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(alcohol_use_3 = c(0,10,60)), by = "alcohol_use_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(alcohol_use_4, knots = quantile(alcohol_use_4, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(alcohol_use_4, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_heavy = map(.x = outcome_model_t3_t4,
                                 ~ avg_comparisons(.x, variables = list(alcohol_use_4 = c(10,60)), vcov = "HC3")),
         effect_t3_t4_absti = map(.x = outcome_model_t3_t4,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_4 = c(10,0)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                              ~ avg_predictions(.x, variables = list(alcohol_use_4 = c(0,10,60)), by = "alcohol_use_4", vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
alcohol_loneli_effects <- select(alcohol_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_heavy:effect_t3_t4_absti, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 

# extracting and pooling estimates across imputed datasets: average predicted values for -1SD, mean and +1SD of loneliness
alcohol_loneli_avg_pred <- select(alcohol_loneli_models, starts_with("avg_pred")) %>% 
  pivot_longer(cols = avg_pred_t1_t2:avg_pred_t3_t4, names_to = "name", values_to = "means") %>% 
  group_by(name) %>%  nest() %>% 
  mutate(data = map(.x = data, ~ .x %>% mutate(abstinence = list(.x$means[[1]][1,]))),
         data = map(.x = data, ~ .x %>% mutate(low = list(.x$means[[1]][2,]))),
         data = map(.x = data, ~ .x %>% mutate(heavy = list(.x$means[[1]][3,]))),
         abstinence = map(.x = data, ~ summary(pool(.x$abstinence), conf.int = TRUE)),
         low = map(.x = data, ~ summary(pool(.x$low), conf.int = TRUE)),
         heavy = map(.x = data, ~ summary(pool(.x$heavy), conf.int = TRUE))) %>%  
  select(abstinence, low, heavy) %>% pivot_longer(cols = c(abstinence, low, heavy), names_to = "level", values_to = "value") %>% unnest() %>% 
  select(name, level, estimate, `2.5 %`, `97.5 %`) 

# Alcohol use
cond_predictions_alc2 <- as.data.frame(predictions(alcohol_loneli_models$outcome_model_t1_t2[[1]], variables = list(alcohol_use_2 = seq(from = 0, to = 81, by = 0.1)), newdata = "mean")) %>% 
  mutate(alcohol_use = alcohol_use_2, effect = "T1 on T2")
cond_predictions_alc3 <- as.data.frame(predictions(alcohol_loneli_models$outcome_model_t2_t3[[1]], variables = list(alcohol_use_3 = seq(from = 0, to = 71, by = 0.1)), newdata = "mean")) %>% 
  mutate(alcohol_use = alcohol_use_3, effect = "T2 on T3")
cond_predictions_alc4 <- as.data.frame(predictions(alcohol_loneli_models$outcome_model_t3_t4[[1]], variables = list(alcohol_use_4 = seq(from = 0, to = 71, by = 0.1)), newdata = "mean")) %>% 
  mutate(alcohol_use = alcohol_use_4, effect = "T3 on T4")

funct_form_alc_loneli <- cond_predictions_alc2 %>% 
  select(estimate, conf.low, conf.high, alcohol_use, effect) %>% 
  ggplot(aes(x = alcohol_use, y = estimate, 
             ymin = conf.low, ymax = conf.high)) + 
  geom_line(aes(color = effect), size = 1) + 
  geom_ribbon(alpha = 0.2) +
  labs(title = "Predicted Values of Loneliness",
       x = "Alcohol Use (standard units per 4 weeks)",
       y = "Loneliness (standardized)",
       color = "Prospective Path",
       fill = "Prospective Path") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "funct_form_alc_loneli.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 7, 
       height = 5,  
       bg="white",
       dpi=900)





### Physical Activity ###

#distr_phys_exer <- ecdf(UiN_data_analysis$phys_exer_3)
#percentile_phys_exer <- distr_phys_exer(150)
#quantile(UiN_data_analysis$phys_exer_2, probs = percentile_phys_exer, na.rm = T) # 3 physical activities at T2 ~= 150 min at T3

phys_exer_loneli_models <- imp_UiN_data %>% 
  group_by(.imp) %>% 
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
         avg_pred_t1_t2 = map(.x = outcome_model_t1_t2,
                              ~ avg_predictions(.x, variables = list(phys_exer_2 = c(0,150)), by = "phys_exer_2", vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(phys_exer_3, knots = quantile(phys_exer_3, probs = 0.50), Boundary.knots = quantile(phys_exer_3, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3 = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(phys_exer_3 = c(0, 150)), vcov = "HC3")),
         avg_pred_t2_t3 = map(.x = outcome_model_t2_t3,
                              ~ avg_predictions(.x, variables = list(phys_exer_3 = c(0,150)), by = "phys_exer_3", vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(phys_exer_4, knots = quantile(phys_exer_4, probs = 0.50), Boundary.knots = quantile(phys_exer_4, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4 = map(.x = outcome_model_t3_t4,
                                  ~ avg_comparisons(.x, variables = list(phys_exer_4 = c(0, 150)), vcov = "HC3")),
         avg_pred_t3_t4 = map(.x = outcome_model_t3_t4,
                              ~ avg_predictions(.x, variables = list(phys_exer_4 = c(0, 150)), by = "phys_exer_4", vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
phys_exer_loneli_effects <- select(phys_exer_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2:effect_t3_t4, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 




### MAIN RESULTS ###
# plots and tables

windowsFonts(Times = windowsFont("Times New Roman"))

# CAUSAL EFFECT OF HEALTH FACTORS ON LONELINESS #
ATE_health_loneli <- bind_rows(BMI_loneli_effects, alcohol_loneli_effects, 
                               smoke_loneli_effects, phys_exer_loneli_effects) %>% 
  ungroup() %>% 
  mutate(path = str_sub(effect,8,12),
         path = case_when(path == "t1_t2" ~ "adolescence \n (17y)",
                          path == "t2_t3" ~ "emerg. adulthood \n (22y)",
                          path == "t3_t4" ~ "young adulthood \n (28y)"),
         contrast = case_when(contrast == "30 - 21.75" ~ "adiposity",
                              contrast == "21.75 - 17" ~ "underweight",
                              contrast == "60 - 10" ~ "heavy",
                              contrast == "10 - 0" ~ "abstinence",
                              contrast == "1 - 0" ~ "smoking",
                              contrast %in% c("3 - 0", "150 - 0") ~ "moderate"),
         health_factor = case_when(contrast %in% c("abstinence", "heavy") ~ "Alcohol use",
                                   contrast %in% c("underweight", "adiposity") ~ "Body mass index",
                                   contrast == "smoking" ~ "Smoking",
                                   contrast == "moderate" ~ "Physical activity"),
         across(c(estimate, `2.5 %`, `97.5 %`), 
                ~ case_when(contrast %in% c("abstinence", "underweight") ~ .x * -1,
                            TRUE ~ .x)))
write.csv(ATE_health_loneli, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_health_loneli_raw")
ATE_health_loneli <- read_csv("N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_health_loneli_raw")

results_ATE_health_loneli <- ATE_health_loneli %>% 
  mutate(est = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                      " (", formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                      formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")")) %>% 
  select(health_factor, contrast, path, est) %>% 
  pivot_wider(names_from = "path",
              values_from = "est") %>% 
  arrange(health_factor, contrast)
write.csv(results_ATE_health_loneli, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_health_loneli")


ATE_bmi_loneli_plot <- ATE_health_loneli %>% 
  filter(contrast %in% c("adiposity", "underweight")) %>% 
  mutate(contrast = factor(contrast, levels = c("underweight", "adiposity"))) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = contrast)) +
  geom_point(size = 5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Body Mass Index",
       x = "",
       y = "",
       color = "ref.: normative",
       shape = "ref.: normative") +
  theme_minimal() +
  scale_color_manual(values = c("adiposity" = "#5773CCFF", "underweight" = "#FFB900FF")) +  
  scale_shape_manual(values = c("adiposity" = 16, "underweight" = 17)) +  
  scale_y_continuous(limits = c(-0.25, 0.42), breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4)) +
  theme(strip.text = element_text(size = 15.5, family = "Times"),
        axis.text = element_text(size = 15.5, family = "Times", color = "black"),
        text = element_text(size = 15.5, family = "Times"),
        plot.title = element_text(face = "bold", size = 17.5, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15.5, family = "Times"), 
        plot.margin = margin(t = 10, l = 15),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_alcohol_loneli_plot <- ATE_health_loneli %>% 
  filter(contrast %in% c("heavy", "abstinence")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = contrast)) +
  geom_point(size = 5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Alcohol Use",
       x = "",
       y = "",
       color = "ref.: moderate",
       shape = "ref.: moderate") +
  theme_minimal() +
  scale_color_manual(values = c("heavy" = "#5773CCFF", "abstinence" = "#FFB900FF")) +  
  scale_shape_manual(values = c("heavy" = 16, "abstinence" = 17)) +  
  scale_y_continuous(limits = c(-0.25, 0.42), breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4)) +
  theme(strip.text = element_text(size = 15.5, family = "Times"),
        axis.text = element_text(size = 15.5, family = "Times", color = "black"),
        text = element_text(size = 15.5, family = "Times"),
        plot.title = element_text(face = "bold", size = 17.5, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15.5, family = "Times"), 
        plot.margin = margin(t = 10, l = 10),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_smoke_loneli_plot <- ATE_health_loneli %>% 
  filter(contrast %in% c("smoking")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast)) +
  geom_point(size = 5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.2, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Smoking",
       x = "",
       y = "",  
       color = "ref.: none") +
  theme_minimal() +
  scale_color_manual(values = c("smoking" = "#5773CCFF")) +  
  scale_y_continuous(limits = c(-0.27, 0.21), breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2)) +
  theme(strip.text = element_text(size = 15.5, family = "Times"),
        axis.text = element_text(size = 15.5, family = "Times", color = "black"),
        text = element_text(size = 15.5, family = "Times"),
        plot.title = element_text(face = "bold", size = 17.5, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15.5, family = "Times"), 
        plot.margin = margin(l = 10, t = 40),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_exer_loneli_plot <- ATE_health_loneli %>% 
  filter(contrast %in% c("moderate")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast)) +
  geom_point(size = 5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.2, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Physical Exercise",
       x = "",
       y = "",
       color = "ref.: inactive") +
  theme_minimal() +
  scale_color_manual(values = c("moderate" = "#5773CCFF")) +  
  scale_y_continuous(limits = c(-0.27, 0.21), breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2)) +
  theme(strip.text = element_text(size = 15.5, family = "Times"),
        axis.text = element_text(size = 15.5, family = "Times", color = "black"),
        text = element_text(size = 15.5, family = "Times"),
        plot.title = element_text(face = "bold", size = 17.5, hjust = 0.5, margin = margin(b = 10)),
        legend.text=element_text(size = 15.5, family = "Times"), 
        plot.margin = margin(l = 20, t = 40),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

plot_ATE_health_loneli <- cowplot::plot_grid(ATE_alcohol_loneli_plot, ATE_bmi_loneli_plot, 
                                             ATE_smoke_loneli_plot, ATE_exer_loneli_plot, 
                                             ncol = 2) +
  draw_label("Loneliness (standardized)", angle = 90, x = 0.01, y = 0.5, vjust = 0.5, fontfamily = "Times", size = 15.5)

ggsave(filename = "plot_ATE_health_loneli.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 16.5, 
       height = 12,  
       bg="white",
       dpi=1050)