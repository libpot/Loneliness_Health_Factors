### SENSITIVITY ANALYSIS ###

### LIFE STAGE-SPECIFIC (SAMPLE-BASED) CONTRASTS FOR HEALTH FACTORS ###
# e.g., +1 SD and - 1SD compared to mean BMI value #




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
baseline_covariates_t0 <- "age + sex + ethnicity + parental_education + behavioral_monitoring + psych_overcontrol + cold_parenting + parental_alcoholuse + parental_smoking + chronic_condition +"
time_varying_covariates <- c("employment_", "living_situation_", "relationship_", "friends_",
                             "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                             "depression_", "social_support_")
time_varying_covariates_t0 <- paste(paste0(time_varying_covariates, "1"), collapse = " + ")
time_varying_covariates_t1 <- paste(paste0(time_varying_covariates, "2"), collapse = " + ")
time_varying_covariates_t2 <- paste(c( "education_3", paste0(time_varying_covariates, "3")), collapse = " + ")



# descriptive statistics: life stage-specific contrasts for health factors
life_stage_spec_health_contrasts <- UiN_data_analysis %>%
    select(alcohol_use_2:alcohol_use_4, bmi_2:bmi_4, phys_exer_2:phys_exer_4) %>% 
    pivot_longer(cols = everything(), 
                 names_to = "factor",
                 values_to = "value") %>% 
    mutate(wave = str_sub(factor, nchar(factor)-1, nchar(factor)),
           factor = str_remove(factor, "_2|_3|_4"),
           wave = case_when(wave == "_2" ~ "Adolescence (T1)",
                            wave == "_3" ~ "Emerging adulthood (T2)",
                            wave == "_4" ~ "Young adulthood (T3)")) %>% 
    group_by(factor, wave) %>% 
    summarise(minus = mean(value, na.rm = T) - sd(value, na.rm = T),
              mean = mean(value, na.rm = T),
              plus = mean(value, na.rm = T) + sd(value, na.rm = T)) %>% 
    pivot_longer(cols = c(minus, mean, plus),
                   names_to = "level",
                 values_to = "value") %>% 
    pivot_wider(names_from = "wave",
                values_from = value) %>% 
    filter(factor == "alcohol_use" & level %in% c("mean", "plus") |
           factor == "bmi" |
           factor == "phys_exer" & level == "mean") %>% 
  mutate(across(where(is.numeric), ~  formatC(round(.x, 2), format = "f", digits = 1))) 
write.csv(life_stage_spec_health_contrasts, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/life_stage_spec_health_contrasts")



time_varying_covariates_bmi <- c("employment_", "living_situation_", "relationship_", "friends_",
                                 "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                                 "depression_")
time_varying_covariates_t0_bmi <- paste(paste0(time_varying_covariates_bmi, "1"), collapse = " + ")

# Body Mass Index
BMI_std_loneli_models <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, 
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(bmi_1:bmi_5,
                                           ~ (.x - mean(.x))/sd(.x)), 
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("bmi_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "gaussian", weights = weights, data = .x)),
         data_ext = map2(.x = data,
                         .y = exposure_model_t1,
                         ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data_ext = map2(.x = data_ext,
                         .y = exposure_model_t2,
                         ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data_ext = map2(.x = data_ext,
                         .y = exposure_model_t3,
                         ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))),
         outcome_model_l1_b2 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0_bmi, 
                                                           "+ ns(bmi_2, knots = quantile(bmi_2, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_2, probs = c(0.05, 0.95)))")),
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l1_b2_plus = map(.x = outcome_model_l1_b2,
                                 ~ avg_comparisons(.x, variables = list(bmi_2 = c(0,1)), vcov = "HC3")), 
         effect_l1_b2_minus = map(.x = outcome_model_l1_b2,
                                  ~ avg_comparisons(.x, variables = list(bmi_2 = c(0,-1)), vcov = "HC3")),
         outcome_model_l2_b3 = map(.x = data_ext,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(bmi_3, knots = quantile(bmi_3, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_3, probs = c(0.05, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l2_b3_plus = map(.x = outcome_model_l2_b3,
                                 ~ avg_comparisons(.x, variables = list(bmi_3 = c(0,1)), vcov = "HC3")),
         effect_l2_b3_minus = map(.x = outcome_model_l2_b3,
                                  ~ avg_comparisons(.x, variables = list(bmi_3 = c(0,-1)), vcov = "HC3")),
         outcome_model_l3_b4 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(bmi_4, knots = quantile(bmi_4, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bmi_4, probs = c(0.05, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l3_b4_plus = map(.x = outcome_model_l3_b4,
                                 ~ avg_comparisons(.x, variables = list(bmi_4 = c(0,1)), vcov = "HC3")),
         effect_l3_b4_minus = map(.x = outcome_model_l3_b4,
                                  ~ avg_comparisons(.x, variables = list(bmi_4 = c(0,-1)), vcov = "HC3")))

BMI_std_loneli_effects <- select(BMI_std_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_l1_b2_plus:effect_l3_b4_minus, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% 
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`)



### Current smoking ###

smoke_std_loneli_models <- imp_UiN_data %>% 
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
         effect_t1_t2_smoke = map(.x = outcome_model_t1_t2,
                                  ~ avg_comparisons(.x, variables = list(smoking_2 = c(0,1)), vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ smoking_3")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_smoke = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(smoking_3 = c(0,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ smoking_4")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_smoke = map(.x = outcome_model_t3_t4,
                                 ~ avg_comparisons(.x, variables = list(smoking_4 = c(0,1)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
smoke_std_loneli_effects <- select(smoke_std_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_smoke:effect_t3_t4_smoke, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 




### Alcohol Use ###

alcohol_std_loneli_models <- imp_UiN_data %>% 
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
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_2 = c(12.5,38.9)), vcov = "HC3")), # mean = 13.2 drinks, 1 SD above the mean = 45.2 drinks
         effect_t1_t2_absti = map(.x = outcome_model_t1_t2,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_2 = c(12.5,0)), vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(alcohol_use_3, knots = quantile(alcohol_use_3, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(alcohol_use_3, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_heavy = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_3 = c(26.2,64.2)), vcov = "HC3")),
         effect_t2_t3_absti = map(.x = outcome_model_t2_t3,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_3 = c(26.2,0)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(alcohol_use_4, knots = quantile(alcohol_use_4, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(alcohol_use_4, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_heavy = map(.x = outcome_model_t3_t4,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_4 = c(23.4,58.3)), vcov = "HC3")),
         effect_t3_t4_absti = map(.x = outcome_model_t3_t4,
                                  ~ avg_comparisons(.x, variables = list(alcohol_use_4 = c(23.4,0)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
alcohol_std_loneli_effects <- select(alcohol_std_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_heavy:effect_t3_t4_absti, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 






### Physical Activity ###

phys_exer_std_loneli_models <- imp_UiN_data %>% 
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
         effect_t1_t2_moderate = map(.x = outcome_model_t1_t2,
                                     ~ avg_comparisons(.x, variables = list(phys_exer_2 = c(0,176.9)), vcov = "HC3")), 
         outcome_model_t2_t3 = map(.x = data,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(phys_exer_3, knots = quantile(phys_exer_3, probs = 0.50), Boundary.knots = quantile(phys_exer_3, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_moderate = map(.x = outcome_model_t2_t3,
                                     ~ avg_comparisons(.x, variables = list(phys_exer_3 = c(0, 176.7)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(phys_exer_4, knots = quantile(phys_exer_4, probs = 0.50), Boundary.knots = quantile(phys_exer_4, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_moderate = map(.x = outcome_model_t3_t4,
                                     ~ avg_comparisons(.x, variables = list(phys_exer_4 = c(0, 146.5)), vcov = "HC3")))

# extracting and pooling estimates across imputed datasets: marginal (average treatment) effects of 1SD loneliness increase compared to mean and -1SD
phys_exer_std_loneli_effects <- select(phys_exer_std_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_moderate:effect_t3_t4_moderate, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 




### MAIN PLOTS ###

# CAUSAL EFFECT OF HEALTH FACTORS ON LONELINESS #
ATE_health_std_loneli <- bind_rows(BMI_std_loneli_effects, alcohol_std_loneli_effects, 
                                   smoke_std_loneli_effects, phys_exer_std_loneli_effects) %>% 
  ungroup() %>% 
  mutate(path = str_sub(effect,8,12),
         path = case_when(path %in% c("t1_t2", "l1_b2") ~ "adolescence \n (17y)",
                          path %in% c("t2_t3", "l2_b3") ~ "emerg. adulthood \n (22y)",
                          path %in% c("t3_t4", "l3_b4") ~ "young adulthood \n (28y)"),
         contrast = case_when(str_detect(effect, "plus") ~ "overweight",
                              str_detect(effect, "minus") ~ "underweight",
                              str_detect(effect, "heavy") ~ "heavy",
                              str_detect(effect, "absti") ~ "abstinence",
                              str_detect(effect, "smoke") ~ "smoking",
                              str_detect(effect, "moderate") ~ "moderate"),
         across(c(estimate, `2.5 %`, `97.5 %`), 
                ~ case_when(contrast %in% c("abstinence", "underweight") ~ .x * -1,
                            TRUE ~ .x)))

results_ATE_health_std_loneli <- ATE_health_std_loneli %>% 
  mutate(est = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                      " (", formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                      formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")"),
         factor = case_when(contrast %in% c("abstinence", "heavy") ~ "Alcohol use",
                            contrast %in% c("underweight", "overweight") ~ "Body mass index",
                            contrast == "smoking" ~ "Smoking",
                            contrast == "moderate" ~ "Physical activity")) %>% 
  select(factor, contrast, path, est) %>% 
  pivot_wider(names_from = "path",
              values_from = "est") %>% 
  arrange(factor, contrast)
write.csv(results_ATE_health_std_loneli, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_ATE_health_std_loneli")

  


ATE_bmi_loneli_plot <- ATE_health_std_loneli %>% 
  filter(contrast %in% c("overweight", "underweight")) %>% 
  mutate(contrast = factor(contrast, levels = c("underweight", "overweight"))) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Body Mass Index",
       x = "",
       y = "",
       color = "ref.: normative",
       shape = "ref.: normative") +
  theme_minimal() +
  scale_color_manual(values = c("overweight" = "black", "underweight" = "darkgrey")) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(t = 40),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_alcohol_loneli_plot <- ATE_health_std_loneli %>% 
  filter(contrast %in% c("heavy", "abstinence")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Alcohol Use",
       x = "",
       y = "",
       color = "ref.: occasional/low",
       shape = "ref.: occasional/low") +
  theme_minimal() +
  scale_color_manual(values = c("heavy" = "black", "abstinence" = "darkgrey")) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(t = 40, l = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_smoke_loneli_plot <- ATE_health_std_loneli %>% 
  filter(contrast %in% c("smoking")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Current smoking",
       x = "",
       y = "",
       color = "ref.: no smoking",
       shape = "ref.: no smoking") +
  theme_minimal() +
  scale_color_manual(values = c("smoking" = "black")) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.margin = margin(l = 20),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ATE_exer_loneli_plot <- ATE_health_std_loneli %>% 
  filter(contrast %in% c("moderate")) %>% 
  ggplot(aes(x = path, y = estimate, color = contrast, shape = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  labs(title = "Physical Exercise",
       x = "",
       y = "",
       color = "ref.: inactive",
       shape = "ref.: inactive") +
  theme_minimal() +
  scale_color_manual(values = c("moderate" = "black")) +  
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

plot_ATE_health_std_loneli <- cowplot::plot_grid(ATE_alcohol_loneli_plot, ATE_bmi_loneli_plot, 
                                                 ATE_smoke_loneli_plot, ATE_exer_loneli_plot, 
                                                 ncol = 2) +
  draw_label("Causal Effect of Health Factors on Loneliness", fontface = 'bold', fontfamily = "Times", size = 16, hjust = 0.5, y = 0.98) +
  draw_label("Life stage-specific contrasts", fontfamily = "Times", size = 14, hjust = 0.5, y = 0.95) +
  draw_label("Loneliness (standardized)", angle = 90, x = 0.01, y = 0.5, vjust = 0.5, fontfamily = "Times", size = 12)
#draw_label("Prospective path", y = 0.03, x = 0.5, hjust = 0.5, fontfamily = "Times", size = 14)

ggsave(filename = "plot_ATE_health_std_loneli.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 12, 
       height = 9,  
       bg="white",
       dpi=900)
