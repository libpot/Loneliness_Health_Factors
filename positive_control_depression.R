##################################################################################
### SENSITIVITY ANALYSIS: INTERPLAY BETWEEN LONELINESS AND DEPRESSIVE SYMPTOMS ### 
###################################################################################

# loading packages
library(tidyverse) # data pre-processing
library(splines) # natural (restricted cubic) splines
library(marginaleffects) # marginal (average) treatment effects
library(sandwich) # robust (sandwich) standard errors
library(ggtext) # fine-tuning the plots
library(cowplot) # arranging multiple plots
library(mice) # pooling estimates across imputations

# loading multiply imputed dataset and the incomplete (original) dataset
UiN_data_weights <- read.csv("N:/durable/Data_analyses/Libor/Social_Gradient_Mental_Health_Lifespan/data/UiN_data_IPAW")
imp_UiN_data_dep5_ext <- read.csv("N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Data/imp_UiN_data_dep5_ext") %>% 
  merge(., UiN_data_weights, by = "id")

# specifying a list of baseline and time-varying covariates
baseline_covariates_t0 <- "age + sex + ethnicity + parental_education + behavioral_monitoring + psych_overcontrol + cold_parenting + parental_alcoholuse + parental_smoking +"
time_varying_covariates <- c("employment_", "living_situation_", "relationship_", "friends_",
                             "smoking_", "alcohol_use_", "phys_exer_", "bmi_", "loneliness_",
                             "depression_", "social_support_")
time_varying_covariates_t0 <- paste(paste0(time_varying_covariates, "1"), collapse = " + ")
time_varying_covariates_t1 <- paste(paste0(time_varying_covariates, "2"), collapse = " + ")
time_varying_covariates_t2 <- paste(c( "education_3", paste0(time_varying_covariates, "3")), collapse = " + ")



### EXPOSURE MODELS ###

exposure_models_loneli <- imp_UiN_data_dep5_ext %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5,
                                           ~ (.x - mean(.x))/sd(.x)),   # standardizing loneliness scores (mean = 0, SD = 1)
                                    across(depression_1:depression_5,
                                           ~ (.x - mean(.x))/sd(.x)),
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

exposure_models_depress <- imp_UiN_data_dep5_ext %>% 
  group_by(.imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    ~ .x %>% mutate(across(loneliness_1:loneliness_5,
                                           ~ (.x - mean(.x))/sd(.x)),   # standardizing loneliness scores (mean = 0, SD = 1)
                                    across(depression_1:depression_5,
                                           ~ (.x - mean(.x))/sd(.x)),
                                    across(smoking_3:smoking_5,
                                           ~ ifelse(.x == "smoking", 1, 0)))),
         exposure_model_t1 = map(.x = data,
                                 ~ glm(as.formula(paste0("depression_2 ~", baseline_covariates_t0, time_varying_covariates_t0)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t2 = map(.x = data,
                                 ~ glm(as.formula(paste0("depression_3 ~", baseline_covariates_t0, time_varying_covariates_t1)), family = "gaussian", weights = weights, data = .x)),
         exposure_model_t3 = map(.x = data,
                                 ~ glm(as.formula(paste0("depression_4 ~", baseline_covariates_t0, time_varying_covariates_t2)), family = "gaussian", weights = weights, data = .x)),
         data = map2(.x = data,
                     .y = exposure_model_t1,
                     ~ .x %>% mutate(conditional_density_t1 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t2,
                     ~ .x %>% mutate(conditional_density_t2 = predict.glm(.y, type = "response"))),
         data = map2(.x = data,
                     .y = exposure_model_t3,
                     ~ .x %>% mutate(conditional_density_t3 = predict.glm(.y, type = "response"))))


### OUTCOME MODELS ###

# loneliness -> depression
depress_outcome_models <- exposure_models_loneli %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("depression_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ loneliness_2")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(loneliness_2 = c(0,1)), vcov = "HC3")),
         outcome_model_t2_t3 = map(.x = data,
                                   ~ glm(as.formula(paste0("depression_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ loneliness_3")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(loneliness_3 = c(0,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("depression_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ loneliness_4")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(loneliness_4 = c(0,1)), vcov = "HC3")))

depress_outcome_effects <- select(depress_outcome_models, starts_with("effect")) %>% 
  pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_1sd, names_to = "effect", values_to = "contrasts") %>% 
  group_by(effect) %>%  nest() %>% 
  mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 

# depression -> loneliness
loneli_outcome_models <- exposure_models_depress %>% 
  mutate(outcome_model_t1_t2 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_3 ~", baseline_covariates_t0, "+ conditional_density_t1 + ", time_varying_covariates_t0, 
                                                           "+ depression_2")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t1_t2_1sd = map(.x = outcome_model_t1_t2,
                                ~ avg_comparisons(.x, variables = list(depression_2 = c(0,1)), vcov = "HC3")),
         outcome_model_t2_t3 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_4 ~", baseline_covariates_t0, "+ conditional_density_t2 + ", time_varying_covariates_t1, 
                                                           "+ depression_3")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t2_t3_1sd = map(.x = outcome_model_t2_t3,
                                ~ avg_comparisons(.x, variables = list(depression_3 = c(0,1)), vcov = "HC3")),
         outcome_model_t3_t4 = map(.x = data,
                                   ~ glm(as.formula(paste0("loneliness_5 ~", baseline_covariates_t0, "+ conditional_density_t3 + ", time_varying_covariates_t2, 
                                                           "+ depression_4")), 
                                         family = "gaussian", weights = weights, data = .x)),
         effect_t3_t4_1sd = map(.x = outcome_model_t3_t4,
                                ~ avg_comparisons(.x, variables = list(depression_4 = c(0,1)), vcov = "HC3")))
         
         # extracting and pooling estimates across imputed datasets: marginal (average) treatment effects of 1SD loneliness increase compared to mean and -1SD
loneli_outcome_effects <- select(loneli_outcome_models, starts_with("effect")) %>% 
           pivot_longer(cols = effect_t1_t2_1sd:effect_t3_t4_1sd, names_to = "effect", values_to = "contrasts") %>% 
           group_by(effect) %>%  nest() %>% 
           mutate(contrasts = map(.x = data, ~ summary(pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
           select(contrasts) %>% unnest() %>% 
           select(effect, contrast, estimate, `2.5 %`, `97.5 %`) 


plot_loneli_dep <- bind_rows(loneli_outcome_effects, depress_outcome_effects) %>% 
  ungroup() %>% 
  mutate(outcome = c(rep("Causal Effect of Depressive Symptoms on Loneliness", 3), rep("Causal Effect of Loneliness on Depressive Symptoms", 3)),
         path = str_sub(effect,8,12),
         path = case_when(path == "t1_t2" ~ "adolescence \n (17y)",
                          path == "t2_t3" ~ "emerg. adulthood \n (22y)",
                          path == "t3_t4" ~ "young adulthood \n (28y)")) %>% 
  ggplot(aes(x = path, y = estimate)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ outcome, 
             ncol = 2, scales = "free") +
  labs(title = "Positive Control Analysis",
       x = "",
       y = "standardized outcome") +
  theme_minimal() +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2)) +
  theme(strip.text.x = ggtext::element_markdown(size = 12),
        text = element_text(family = "Times"),
        plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 10)),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_loneli_dep.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 9, 
       height = 6,  
       bg="white",
       dpi=700)

