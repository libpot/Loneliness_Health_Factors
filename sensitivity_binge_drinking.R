### SENSITIVITY ANALYSIS ###

### BINGE DRINKING (HEAVY EPISODIC DRINKING) frequency instead of total alcohol consumption (drinks per 4 weeks) ### 



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
                             "smoking_", "bingedrink_", "phys_exer_", "bmi_", "loneliness_",
                             "depression_", "social_support_")
time_varying_covariates_t0 <- paste(paste0(time_varying_covariates, "1"), collapse = " + ")
time_varying_covariates_t1 <- paste(paste0(time_varying_covariates, "2"), collapse = " + ")
time_varying_covariates_t2 <- paste(c( "education_3", paste0(time_varying_covariates, "3")), collapse = " + ")

# binge drinking -> loneliness
binge_loneli_models <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(bingedrink_1:bingedrink_5, # standardizing loneliness and BMI scores (mean = 0, SD = 1)
                                           ~ .x - 1),
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))), 
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("bingedrink_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "poisson", weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("bingedrink_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "poisson", weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("bingedrink_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "poisson", weights = weights, data = .x)),
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
                                   ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ ns(bingedrink_2, knots = quantile(bingedrink_2, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bingedrink_2, probs = c(0.00, 0.95)))")),
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l1_b2_heavy = map(.x = outcome_model_l1_b2,
                                  ~ avg_comparisons(.x, variables = list(bingedrink_2 = c(2,5)), vcov = "HC3")), # vcov = HC3: robust (sandwich) standard errors (SE)
         effect_l1_b2_absti = map(.x = outcome_model_l1_b2,
                                  ~ avg_comparisons(.x, variables = list(bingedrink_2 = c(2,0)), vcov = "HC3")), # marginal (counter-factual) contrasts
         outcome_model_l2_b3 = map(.x = data_ext,                         
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ ns(bingedrink_3, knots = quantile(bingedrink_3, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bingedrink_3, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l2_b3_heavy = map(.x = outcome_model_l2_b3,
                                  ~ avg_comparisons(.x, variables = list(bingedrink_3 = c(2,5)), vcov = "HC3")),
         effect_l2_b3_absti = map(.x = outcome_model_l2_b3,
                                  ~ avg_comparisons(.x, variables = list(bingedrink_3 = c(2,0)), vcov = "HC3")),
         outcome_model_l3_b4 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ ns(bingedrink_4, knots = quantile(bingedrink_4, probs = c(0.25, 0.50, 0.75)), Boundary.knots = quantile(bingedrink_4, probs = c(0.00, 0.95)))")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l3_b4_heavy = map(.x = outcome_model_l3_b4,
                                  ~ avg_comparisons(.x, variables = list(bingedrink_4 = c(2,5)), vcov = "HC3")),
         effect_l3_b4_absti = map(.x = outcome_model_l3_b4,
                                  ~ avg_comparisons(.x, variables = list(bingedrink_4 = c(2,0)), vcov = "HC3")))

binge_loneli_effects <- select(binge_loneli_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_l1_b2_heavy:effect_l3_b4_absti, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`)


# loneliness -> binge drinking 
loneli_binge_models <- imp_UiN_data %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5,
                                           ~ (.x - mean(.x))/sd(.x)))), # standardizing loneliness scores (mean = 0, SD = 1)
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
                         ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))),
         outcome_model_l1_b2 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("bingedrink_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ loneliness_2")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l1_b2_1sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")), 
         effect_l1_b2_2sd = map(.x = outcome_model_l1_b2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(-1,1)), vcov = "HC3")), 
         outcome_model_l2_b3 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("bingedrink_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ loneliness_3")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l2_b3_1sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         effect_l2_b3_2sd = map(.x = outcome_model_l2_b3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(-1,1)), vcov = "HC3")),
         outcome_model_l3_b4 = map(.x = data_ext,
                                   ~ glm(as.formula(paste0("bingedrink_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ loneliness_4")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_l3_b4_1sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")),
         effect_l3_b4_2sd = map(.x = outcome_model_l3_b4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(-1,1)), vcov = "HC3")))

loneli_binge_effects <- select(loneli_binge_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_l1_b2_1sd:effect_l3_b4_2sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling estimates according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`)



### RESULTS ###
ATE_binge <- rbind(binge_loneli_effects, loneli_binge_effects) %>% 
  ungroup() %>% 
  mutate(contrast = case_when(contrast == "1 - -1" ~ "+1 SD vs. -1 SD",
                              contrast == "1 - 0" ~ "+1 SD vs. Mean",
                              contrast == "5 - 2" ~ "Heavy",
                              contrast == "2 - 0" ~ "None"),
         path = str_sub(effect, 8, 12),
         path = case_when(path == "l1_b2" ~ "adolescence \n (17y)",
                          path == "l2_b3" ~ "emerg. adulthood \n (22y)",
                          path == "l3_b4" ~ "young adulthood \n (28y)"),
         across(c(estimate, `2.5 %`, `97.5 %`), 
                ~ case_when(contrast %in% c("None") ~ .x * -1,
                            TRUE ~ .x)))

results_binge <- ATE_binge %>% 
  mutate(est = paste0(formatC(round(estimate, 2), format = "f", digits = 2), 
                      " (", formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                      formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")")) %>% 
  select(contrast, path, est) %>% 
  pivot_wider(names_from = "path",
              values_from = "est") %>% 
  arrange(contrast)
write.csv(results_binge, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/results_binge")


plot_binge <- ATE_binge %>% 
  filter(contrast != "+1 SD vs. -1 SD") %>% 
  mutate(label = c(rep("Effect of Binge Drinking on Loneliness", 6), rep("Effect of Loneliness on Binge Drinking", 3)),
         contrast = ifelse(contrast == "+1 SD vs. Mean", "+1 SD", contrast)) %>% 
  ggplot(aes(x = path, y = estimate, shape = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.3), color = "black") +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ label, 
             ncol = 2, scales = "free") +
  labs(title = "",
       x = "",
       y = "",
       shape = "Contrast") +
  theme_minimal() +
  theme(strip.text.x = element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_binge.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 9, 
       height = 6,  
       bg="white",
       dpi=700)




