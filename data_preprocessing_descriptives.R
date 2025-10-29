######################################################################
### THE DYNAMIC INTERPLAY BETWEEN LONELINESS AND HEALTH BEHAVIORS ###
#####################################################################

           ############################################    
           ### DATA PRE-PROCESSING and DESCRIPTIVES ###    
           ############################################    
 


# loading packages
library(haven) # data loading
library(tidyverse) # data pre-processing
library(mice) # multiple imputation

# loading the data
UiN_data_orig <- read_sav("N:/durable/Data_analyses/Libor/data/UiN_regfile_030124.sav") 



### SELECTING VARIABLES ###

UiN_data_analysis <- UiN_data_orig %>% 
# selecting participants
  filter(!is.na(id), # 0 observations without ID
         Type2_1 %in% c(1:8),   # longitudinal study sub-sample
         partici1 == 1    
         & partici2 == 1 & (partici3 == 1 | partici4 == 1)) %>%  # participated in at least three data collection waves (n = 3072)
  select(id, # identification number
         
# BASELINE COVARIATES
# sociodemographic factors
         Gende1_1, gende, # gender (self-reported and register-based)
         Age1_1, Age2_1, Birth9_3, # age and birth date
         ParEd1_1:ParEd1_9, SOSBAK, # parental education (self-reported and register-based)
         Origi1p1, Origi1p2, # ethnic background
         age_1chl, child_t3,

# early family environment (parenting)
         PM1_1:PM1_6, # behavioral_monitoring
         PBI1_01:PBI1_10, # items 1-5 psychological over-protection/over-control and item 6-10 warmth/care
# parental health behaviors
         Subst1p1, Subst1p2, Subst1p3, Subst3p4, Subst3p5, # parental alcohol use 
         Smok3_02, Smok3_03, # parental smoking
# chronic conditions
        Disa1_03, # asthma
        Disa1_05, # physical disability
        Disa1_02, # writing/reading disorder
# LONELINESS
         starts_with("UCLA") & !ends_with("5") & !ends_with(c("n", "m")), # UCLA loneliness scale
         starts_with("UCLA") & ends_with("5") & !ends_with(c("n", "m")), # direct loneliness measure

# HEALTH BEHAVIORS
         starts_with("Subst") & ends_with(c("a1", "a2")), # alcohol use (frequency and quantity, respectively)
         starts_with("Subst") & ends_with("_1"), # binge drinking
         starts_with("Smoke") & ends_with("_1"), # smoking
         starts_with("Weigh") & ends_with("_1"), # weight
         starts_with("Heigh") & ends_with("_1"), # height
         starts_with("Train") & ends_with(c("_1", "_2", "_3")), # physical activity
         starts_with("LeAc") & ends_with(c("03", "04","06")),
        BFCP2_1a, BFCP2_1b, BFCP3_1a, BFCP3_1b, BFCP4_1a, BFCP4_1b, # friend's smoking
        BFCP2_2a, BFCP2_2b, BFCP3_2a, BFCP3_2b, BFCP4_2a, BFCP4_2b, # friend's alcohol use
        Spor2_3, Spor3_3, Spor4_3, Spor2_6, Spor3_6, Spor4_6, # sports involvement
        Spor2_5a, Spor2_5b, Spor3_5a,Spor3_5b, Spor4_5a,Spor4_5b, # sport activity type
        Spor2_7a, Spor2_7b, Spor3_7a, Spor3_7b, Spor4_7a, Spor4_7b, # sport activity type 
        


# TIME-VARYING COVARIATES
# socio-demographic factors
        eduF1999, eduF2005, # educational attainment
        Educ3_5, Educ3_6, Educ3_8, Educ3_9,
        Educ4_5, Educ4_6, Educ4_8, Educ4_9,
        Occup2_1, Occ2_1a, Occ2_1b, Occ2_1c, Occ2_1d,
ParJo1_1, ParJo1_2, ParJo2_1, ParJo2_2, Occup3_2, Occup4_2, # employment status
# social isolation indicators
        House1_1, House2_1, House3_1, House4_1, # living situation
        Girlf1_1, Partn3_1, Partn3y2, Partn3y3, Marit3_1, Marit4_1, Marit5_1, # relationship status
        SS2b08, SS2c08, SS2d08, SS2e04, 
        SS3b08, SS3c08, SS3d08, SS3e04, 
        SS4b08, SS4c08, SS4d08, SS4e04, 
        SS1b08, SS1c08, SS1d08,  
        SPPA1sa2, SPPA2sa2, SPPA3sa2, SPPA4sa2, # social network size
# depressive symptoms
        starts_with("SCL") & ends_with(c("_07", "_08", "_09", "_10", "_11", "_12")) & !ends_with("tot"), 
# perceived social support (in hypothetical situations)
        SS1b01, SS2b01, SS3b01, SS4b01, # educational choice
        SS1c01, SS2c01, SS3c01, SS4c01, # personal problem
        SS1d01, SS2d01, SS3d01, SS4d01, # done something illegal
                SS2e01, SS3e01, SS4e01) # feeling down
  

         
                
###########################
### RE-CODING VARIABLES ###
###########################

UiN_data_analysis <- UiN_data_analysis %>% 
  
# reverse coding and re-coding some items
  mutate(across(c(starts_with("UCLA") & ends_with(c("1", "2")),
                  PBI1_03, PBI1_04, PBI1_05, PBI1_06, PBI1_09, PBI1_10), # higher values indicate over-controlling and cold parenting
                ~ 5 - .x),
         across(starts_with("PM"), # higher values indicate neglectful parenting (low behavioral monitoring)
                ~ 7 - .x),
         Heigh2_1 = ifelse(Heigh2_1 < 100, Heigh2_1 + 100, Heigh2_1), # coding errors: some heights were missing 1 at the beginning (e.g., 68 instead of 168)
         Heigh3_1 = ifelse(Heigh3_1 < 130, (Heigh2_1+Heigh4_1)/2, Heigh3_1),
         across(c(SS1b01, SS1c01, SS1d01),
                ~ ifelse(.x %in% c(2,3), 0, .x)),
         across(c(Educ3_5:Educ4_9),
                ~ case_when(.x == 3 ~ 1,
                            .x %in% c(1,2,4) ~ 0)),
         across(c(Occup3_2, Occup4_2),
                ~ case_when(.x %in% c(1,2,3) ~ "employed",
                            .x == 4 ~ "unemployed",
                            .x %in% c(5,6,7) ~ "inactive")),
         Girlf2_1 = case_when(Partn3_1 == 1 | Partn3y2 > 94 ~ "single",
                              Partn3y2 <= 94 & Partn3y3 > 94 | Partn3y2 <= 94 & is.na(Partn3y3) ~ "in relationship",
                              TRUE ~ NA),
         Girlf2_supp = rowSums(select(., SS2b08, SS2c08, SS2d08, SS2e04), na.rm = F),
         Girlf1_supp = rowSums(select(., SS1b08, SS1c08, SS1d08), na.rm = F),
         Girlf3_supp = rowSums(select(., SS3b08, SS3c08, SS3d08, SS3e04), na.rm = F),
         Girlf4_supp = rowSums(select(., SS4b08, SS4c08, SS4d08, SS4e04), na.rm = F),
         across(c(SPPA1sa2, SPPA2sa2, SPPA3sa2, SPPA4sa2),
                ~ case_when(.x %in% c(1,2) ~ "little",
                            .x %in% c(3,4) ~ "a lot"))) %>% 
# BASELINE COVARIATES
# sociodemographic factors
  mutate(age_birth = 92 - Birth9_3,
         age = ifelse(is.na(Age1_1), Age2_1 - 2, Age1_1),
         age = ifelse(is.na(age), age_birth, age),
         age = ifelse(age > 25, Age2_1 - 2, age), # re-coding non-sensical values
         sex = ifelse(is.na(Gende1_1), gende, Gende1_1),
         sex = ifelse(sex == 1, "female", "male"),
         ethnicity = case_when(Origi1p1 == 1 | Origi1p2 == 1 ~ "Norway-born",
                               Origi1p1 == 2 & Origi1p2 == 2 ~ "migrant"),
         education_sr = case_when(ParEd1_7 %in% c(1,2,3) | ParEd1_6 %in% c(1,2,3) ~ "college/university",
                                  ParEd1_3 %in% c(1,2,3) | ParEd1_4 %in% c(1,2,3) | ParEd1_5 %in% c(1,2,3) | ParEd1_1 %in% c(1,2,3) | ParEd1_2 %in% c(1,2,3) ~ "lower",
                                  TRUE ~ NA),
         education_reg = case_when(SOSBAK %in% c(4,3) ~ "college/university",
                                   SOSBAK  %in% c(2,1) ~ "lower",
                                   SOSBAK == 9 ~ NA),
         parental_education = ifelse(is.na(education_reg), education_sr, education_reg),
         
# early family environment (parenting)
         behavioral_monitoring = rowMeans(select(., PM1_1:PM1_6), na.rm = T),
         psych_overcontrol = rowMeans(select(., PBI1_01:PBI1_05), na.rm = T), 
         cold_parenting = rowMeans(select(., PBI1_06:PBI1_10), na.rm = T),

# parental health behaviors
        parental_alcoholuse = case_when(Subst1p1 %in% c(4,5) | Subst1p2 %in% c(4,5) | Subst1p3 %in% c(4,5) | Subst3p4 %in% c(4,5) | Subst3p5 %in% c(4,5) ~ "heavy use",
                                         Subst1p1 == 1 & Subst1p2 == 1 & Subst1p3 == 1 ~ "abstinence",
                                         Subst1p1 %in% c(2,3) | Subst1p2 %in% c(2,3) | Subst1p3 %in% c(1,2,3)  ~ "occasional/low use",
                                         TRUE ~ NA),
        parental_smoking = case_when(Smok3_02 == 4 | Smok3_03 == 4 ~ "daily smoking",
                                     Smok3_02 %in% c(1,2,3) | Smok3_03 %in% c(1,2,3) ~ "non-smoking",
                                     TRUE ~ NA),
# chronic conditions
        chronic_condition = case_when(Disa1_05 == 1 | Disa1_03 == 1 | Disa1_02 == 1 ~ 1, 
                                      is.na(Disa1_05) & is.na(Disa1_03) & is.na(Disa1_02) ~ NA,
                                      TRUE ~ 0)) %>% 
  
# TIME-VARYING COVARIATES
# educational attainment
 mutate(education_3 = rowSums(select(., Educ3_5, Educ3_6, Educ3_8, Educ3_9), na.rm = T),
         education_4 = rowSums(select(., Educ4_5, Educ4_6, Educ4_8, Educ4_9), na.rm = T),
         across(c(education_3, education_4), 
                ~ case_when(.x >= 1 ~"college/university",
                            .x == 0 ~ "lower",
                            TRUE ~ NA)),
         across(c(eduF1999, eduF2005),
                ~ case_when(.x %in% c(6,7,8) ~ "college/university",
                            .x %in% c(2,3,4,5) ~ "lower")),
         education_3 = ifelse(is.na(eduF1999), education_3, eduF1999),
         education_4 = ifelse(is.na(eduF2005), education_4, eduF2005),
# employment status 
         work_status_1f = case_when(ParJo1_1 %in% c(1,2) ~ "employed",
                                    ParJo1_1 %in% c(3,4,5,6) ~ "unemployed",
                                    TRUE ~ NA),
          work_status_1m = case_when(ParJo1_2 %in% c(1,2) ~ "employed",
                                     ParJo1_2 %in% c(3,4,5,6) ~ "unemployed",
                                     TRUE ~ NA),
          employment_1 = case_when(work_status_1f == "employed"| work_status_1m == "employed" ~ "employed",
                                    is.na(work_status_1f) ~ work_status_1m,
                                    is.na(work_status_1m) ~ work_status_1f,
                                    work_status_1f == "unemployed" & work_status_1m == "unemployed" ~ "unemployed"),
          work_status_2f = case_when(ParJo2_1 %in% c(1,2) ~ "employed",
                                     ParJo2_1 %in% c(3,4,5,6) ~ "unemployed",
                                     TRUE ~ NA),
          work_status_2m = case_when(ParJo2_2 %in% c(1,2) ~ "employed",
                                     ParJo2_2 %in% c(3,4,5,6) ~ "unemployed",
                                     TRUE ~ NA),
          employment_2 = case_when(work_status_2f == "employed"| work_status_2m == "employed" ~ "employed",
                                    is.na(work_status_2f) ~ work_status_2m,
                                    is.na(work_status_2m) ~ work_status_2f,
                                    work_status_2f == "unemployed" & work_status_2m == "unemployed" ~ "unemployed"),
         employment_3 = Occup3_2,
         employment_4 = Occup4_2,
# social isolation indicators
         across(c(House1_1,House2_1), ~ case_when(.x %in% c(2,3,4,5,6,7,8) ~ "not with both parents",
                                                  .x == 1 ~ "with both parents")),
         across(c(House3_1,House4_1), ~ case_when(.x == 2 ~ "alone",
                                                  .x %in% c(1,3,4,5) ~ "not alone")),
         living_situation_1 = House1_1,
         living_situation_2 = House2_1,
         living_situation_3 = House3_1,
         living_situation_4 = House4_1,
         across(c(Girlf1_supp, Girlf2_supp, Girlf3_supp, Girlf4_supp),
                ~ case_when(.x >= 1 ~ "in relationship",
                            .x == 0 ~ "single")),
         across(c(Marit3_1, Marit4_1),
                ~ case_when(.x == 1 ~ "single",
                            .x %in% c(2,3) ~ "in relationship")),
         relationship_1 = case_when(Girlf1_1 == 1 ~ "in relationship",
                                    Girlf1_1 %in% c(2,3) ~ "single"),
         relationship_1 = ifelse(is.na(relationship_1), Girlf1_supp, relationship_1),
         relationship_2 = ifelse(is.na(Girlf2_supp), Girlf2_1, Girlf2_supp),
         relationship_3 = ifelse(is.na(Girlf3_supp), Marit3_1, Girlf3_supp),
         relationship_4 = ifelse(is.na(Girlf4_supp), Marit4_1, Girlf4_supp),
         friends_1 = SPPA1sa2,
         friends_2 = SPPA2sa2,
         friends_3 = SPPA3sa2,
         friends_4 = SPPA4sa2) %>% 
# social support
  mutate(social_support_1 = rowSums(select(., SS1b01,SS1c01,SS1d01), na.rm = F),
         social_support_2 = rowSums(select(., SS2b01,SS2c01,SS2d01,SS2e01), na.rm = F),
         social_support_3 = rowSums(select(., SS3b01,SS3c01,SS3d01,SS3e01), na.rm = F),
         social_support_4 = rowSums(select(., SS4b01,SS4c01,SS4d01,SS4e01), na.rm = F),
         #across(c(social_support_1:social_support_4),
         #      ~ ifelse(.x == max(.x, na.rm = T), .x - 1, .x)), # too few cases with extremely low social support, merging the two lowest categories
# depressive symptoms
         depression_1 = rowMeans(select(., SCL1_07:SCL1_12), na.rm = T),
         depression_2 = rowMeans(select(., SCL2_07:SCL2_12), na.rm = T),
         depression_3 = rowMeans(select(., SCL3_07:SCL3_12), na.rm = T),
         depression_4 = rowMeans(select(., SCL4_07:SCL4_12), na.rm = T),
         depression_5 = rowMeans(select(., SCL5_07:SCL5_12), na.rm = T)) %>% 

# LONELINESS
  mutate(loneliness_1 = rowMeans(select(., UCLA1_1:UCLA1_4), na.rm = T),
         loneliness_2 = rowMeans(select(., UCLA2_1:UCLA2_4), na.rm = T),
         loneliness_3 = rowMeans(select(., UCLA3_1:UCLA3_4), na.rm = T),
         loneliness_4 = rowMeans(select(., UCLA4_1:UCLA4_4), na.rm = T),
         loneliness_5 = rowMeans(select(., UCLA5_1:UCLA5_4), na.rm = T),
         loneli_direct_1 = UCLA1_5,
         loneli_direct_2 = UCLA2_5,
         loneli_direct_3 = UCLA3_5,
         loneli_direct_4 = UCLA4_5,
         loneli_direct_5 = UCLA5_5) %>% 

# HEALTH BEHAVIORS
  mutate(across(starts_with("Smoke") & ends_with("_1"),
                ~ case_when(.x %in% c(1,2,3) ~ "non-smoking",
                            .x %in% c(4,5) ~ "smoking",
                            TRUE ~ NA)),
         across(c(Train3_2, Train4_2, Train5_2), 
                ~ ifelse(is.na(.x), 0, .x)),
         Subst1a2 = ifelse(Subst1a1 == 0, 0, Subst1a2),
         Subst2a2 = ifelse(Subst2a1 == 0, 0, Subst2a2),
         Subst3a2 = ifelse(Subst3a1 == 0, 0, Subst3a2),
         Subst4a2 = ifelse(Subst4a1 == 0, 0, Subst4a2),
         Subst5a2 = ifelse(Subst5a1 == 0, 0, Subst5a2),
         bmi_1 = Weigh1_1/((Heigh1_1/100)*(Heigh1_1/100)),
         bmi_2 = Weigh2_1/((Heigh2_1/100)*(Heigh2_1/100)),
         bmi_3 = Weigh3_1/((Heigh3_1/100)*(Heigh3_1/100)),
         bmi_4 = Weigh4_1/((Heigh4_1/100)*(Heigh4_1/100)),
         bmi_5 = Weigh5_1/((Heigh5_1/100)*(Heigh5_1/100))) %>% 
  mutate(phys_exer_1 = (LeAc1_03 + LeAc1_04 + LeAc1_06)*60,
         phys_exer_2 = (LeAc2_03 + LeAc2_04 + LeAc2_06)*60,
         phys_exer_3 = case_when(Train3_3 == 1 ~ 0,
                                 Train3_3 == 0 ~ Train3_1*60 + Train3_2,
                                 TRUE ~ NA),
         phys_exer_4 = case_when(Train4_3 == 1 ~ 0,
                                 Train4_3 == 0 ~ Train4_1*60 + Train4_2,
                                 TRUE ~ NA),
         phys_exer_5 = case_when(!is.na(Train5_1) ~ Train5_1*60 + Train5_2,
                                 TRUE ~ NA),
         smoking_1 = Smoke1_1,
         smoking_2 = Smoke2_1,  
         smoking_3 = Smoke3_1,
         smoking_4 = Smoke4_1,
         smoking_5 = Smoke5_1,
         bingedrink_1 = Subst1_1,
         bingedrink_2 = Subst2_1,
         bingedrink_3 = Subst3_1,
         bingedrink_4 = Subst4_1,
         bingedrink_5 = Subst5_1,
         alcohol_use_1 = Subst1a1*Subst1a2,
         alcohol_use_2 = Subst2a1*Subst2a2,
         alcohol_use_3 = Subst3a1*Subst3a2,
         alcohol_use_4 = Subst4a1*Subst4a2,
         alcohol_use_5 = Subst5a1*Subst5a2,
         sport_2 = case_when(Spor2_3 == 1 | Spor2_6 == 1 ~ 1, 
                             Spor2_3 %in% c(2,3) & Spor2_6 %in% c(2,3) ~ 0,
                             TRUE ~ NA),
         sport_3 = case_when(Spor3_3 == 1 | Spor3_6 == 1 ~ 1, 
                             Spor3_3 %in% c(2,3) & Spor3_6 %in% c(2,3) ~ 0,
                             TRUE ~ NA),
         sport_4 = case_when(Spor4_3 == 1 | Spor4_6 == 1 ~ 1, 
                             Spor4_3 %in% c(2,3) & Spor4_6 %in% c(2,3) ~ 0,
                             TRUE ~ NA),
         bfsmoke_2 = case_when(BFCP2_1a == 1 | BFCP2_1b == 1 ~ 1, 
                               BFCP2_1a == 2 & BFCP2_1a == 2 ~ 0,
                               TRUE ~ NA),
         bfsmoke_3 = case_when(BFCP3_1a == 1 | BFCP3_1b == 1 ~ 1, 
                               BFCP3_1a == 2 & BFCP3_1a == 2 ~ 0,
                               TRUE ~ NA),
         bfsmoke_4 = case_when(BFCP4_1a == 1 | BFCP4_1b == 1 ~ 1, 
                               BFCP4_1a == 2 & BFCP4_1a == 2 ~ 0,
                               TRUE ~ NA),
         bfalcohol_2 = case_when(BFCP2_2a == 1 | BFCP2_2b == 1 ~ 1, 
                                 BFCP2_2a == 2 & BFCP2_2a == 2 ~ 0,
                                 TRUE ~ NA),
         bfalcohol_3 = case_when(BFCP3_2a == 1 | BFCP3_2b == 1 ~ 1, 
                                 BFCP3_2a == 2 & BFCP3_2a == 2 ~ 0,
                                 TRUE ~ NA),
         bfalcohol_4 = case_when(BFCP4_2a == 1 | BFCP4_2b == 1 ~ 1, 
                                 BFCP4_2a == 2 & BFCP4_2a == 2 ~ 0,
                                 TRUE ~ NA)) %>%

# FIDELITY CHECK: replacing non-nonsensical values with NA
   mutate(across(starts_with("alcohol_use_"),
                 ~ ifelse(.x > 300, NA, .x)),  # more than 300 drinks per 4 weeks (e.g., 30 times x 10 drinks) = implausible and coded as missing 
          across(starts_with("phys_exer"),
                 ~ ifelse(.x > 2400, NA, .x))) %>%  # more than 2400 min (40h) of exercise per week = implausible and coded as missing 
# coding variable type: binary/ordinal/numeric
  mutate(across(c(age, behavioral_monitoring, psych_overcontrol, cold_parenting,
                  starts_with(c("social_support_", "depression_")),
                  starts_with(c("alcohol_use_", "bingedrink_", "phys_exer_", "bmi_")),
                  starts_with(c("loneliness", "loneli_direct", "UCLA"))),
                ~ as.numeric(.x)),
         across(c(sex, ethnicity, parental_education, 
                  parental_smoking, parental_alcoholuse,
                  education_3, education_4, starts_with("employment"),
                  starts_with(c("living_situation_", "relationship_", "friends_")),
                  starts_with("smoking"),
                  starts_with(c("bfsmoke_", "bfalcohol_", "sport_"))),
                ~ as.factor(.x))) %>% 
         
# CREATING FINAL DATASET: selecting final variables
  select(id, age, sex, ethnicity, parental_education, # baseline socio-demographic characteristics
         behavioral_monitoring, psych_overcontrol, cold_parenting, # early family environment 
         parental_alcoholuse, parental_smoking, # parental health behaviors 
         chronic_condition, # chronic condition
         starts_with(c("social_support_", "depression_")), # social support and depression
         education_3, education_4, starts_with("employment"), # socioeconomic factors
         starts_with(c("living_situation_", "relationship_", "friends_")), # social isolation
         starts_with(c("smoking_", "alcohol_use_", "bingedrink_", "phys_exer_", "bmi_")), # health factors
         starts_with(c("bfsmoke_", "bfalcohol_", "sport_")), # context of health behaviors
         starts_with("loneli"), UCLA1_1:UCLA1_4, UCLA2_1:UCLA2_4, UCLA3_1:UCLA3_4, UCLA4_1:UCLA4_4, UCLA5_1:UCLA5_4) # loneliness
write.csv(UiN_data_analysis, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Data/UiN_data_analysis_hbcontext_ext")



### SAMPLE CHARACTERISTICS ###

sample_characteristics_cat <- UiN_data_analysis %>% 
  select(sex, ethnicity, parental_education, employment_1,
         parental_alcoholuse, parental_smoking,
         chronic_condition,
         living_situation_1, relationship_1, friends_1, 
         smoking_1) %>% 
  mutate(across(everything(), as.factor)) %>% 
  pivot_longer(everything()) %>% 
  count(name, value) %>% 
  group_by(name) %>% 
  mutate(prop = formatC(round(n/sum(n)*100, 2), format = "f", digits = 2),
         stat = paste(n, paste0("(", prop, ")"))) %>% 
  filter(!is.na(value)) %>% 
  select(name, value, stat) 

sample_characteristics_conti <- UiN_data_analysis %>% 
  select(age,  
         behavioral_monitoring, psych_overcontrol, cold_parenting,
         depression_1, social_support_1,
         alcohol_use_1, phys_exer_1, bmi_1,
         loneliness_1) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(stat = paste0(formatC(round(mean(value, na.rm = TRUE), 2), format = "f", digits = 2), " (", 
                          formatC(round(sd(value, na.rm = TRUE), 2), format = "f", digits = 2), ")"),
            median = formatC(round(median(value, na.rm = TRUE), 2), format = "f", digits = 2),
            min_max = paste0(formatC(round(min(value, na.rm = TRUE), 0), format = "f", digits = 0), " to ", 
                             formatC(round(max(value, na.rm = TRUE), 0), format = "f", digits = 0)))
            
sample_characteristics <- rbind(sample_characteristics_cat, sample_characteristics_conti) %>% 
  mutate(name = factor(name, levels = c("sex", "age", "ethnicity", "parental_education", "employment_1",
                                        "behavioral_monitoring", "psych_overcontrol", "cold_parenting",
                                        "parental_alcoholuse", "parental_smoking",
                                        "chronic_condition",
                                        "depression_1", "social_support_1", 
                                        "living_situation_1", "relationship_1", "friends_1", 
                                        "alcohol_use_1", "smoking_1", "phys_exer_1", "bmi_1",
                                        "loneliness_1"), ordered = TRUE)) %>% 
  arrange(name) %>% 
  select(-median)
write.csv(sample_characteristics, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/sample_characteristics")



### DESCRIPTIVES: loneliness and health behaviors across time ###

descriptives_cat <- UiN_data_analysis %>% 
  select(smoking_2:smoking_5) %>% 
  mutate(across(everything(), as.factor)) %>% 
  pivot_longer(everything()) %>% 
  count(name, value) %>% 
  group_by(name) %>% 
  mutate(prop = formatC(round(n/sum(n)*100, 2), format = "f", digits = 2),
         stat = paste(n, paste0("(", prop, ")"))) %>% 
  filter(!is.na(value),
         value == "smoking") %>% 
  select(name, stat) 

descriptives_conti <- UiN_data_analysis %>% 
  select(alcohol_use_2:alcohol_use_5, 
         phys_exer_2:phys_exer_5, 
         bmi_2:bmi_5,
         loneliness_2:loneliness_5) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(stat = paste0(formatC(round(mean(value, na.rm = TRUE), 2), format = "f", digits = 2), " (", 
                          formatC(round(sd(value, na.rm = TRUE), 2), format = "f", digits = 2), ")"))

min_max_across_time <- UiN_data_analysis %>% 
  select(alcohol_use_2:alcohol_use_5, 
         phys_exer_2:phys_exer_5, 
         bmi_2:bmi_5,
         loneliness_2:loneliness_5) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(min_max = paste0(formatC(round(min(value, na.rm = TRUE), 0), format = "f", digits = 0), " to ", 
                             formatC(round(max(value, na.rm = TRUE), 0), format = "f", digits = 0))) %>% 
  separate(name, into = c("name", "wave"), sep = "_(?=\\d+$)") %>% 
  pivot_wider(names_from = "wave",
              values_from = "min_max")
write.csv(min_max_across_time, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/min_max_across_time")

descriptives_across_time <- rbind(descriptives_cat, descriptives_conti) %>% 
  separate(name, into = c("name", "wave"), sep = "_(?=\\d+$)") %>% 
  pivot_wider(names_from = "wave",
              values_from = "stat") 
write.csv(descriptives_across_time, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/descriptives_across_time")


# histogram of loneliness scores 
loneliness_data <- UiN_data_analysis %>% 
  select(loneliness_2:loneliness_4) %>% 
  pivot_longer(names_to = "time_wave",
               values_to = "score",
               cols = loneliness_2:loneliness_4) %>% 
  mutate(time_wave = case_when(time_wave == "loneliness_2" ~ "T1",
                               time_wave == "loneliness_3" ~ "T2",
                               time_wave == "loneliness_4" ~ "T3"))

loneliness_descriptives <- loneliness_data %>% 
  group_by(time_wave) %>% 
  summarise(mean = mean(score, na.rm = T),
            plus = mean+sd(score, na.rm = T),
            minus = mean-sd(score, na.rm = T))
  
loneliness_histograms <- ggplot(loneliness_data, aes(x = score)) +
  geom_histogram(bins = 12) +
  facet_wrap(~ time_wave, scales = "free", ncol = 1) +
  geom_vline(loneliness_descriptives, mapping = aes(xintercept = mean), size = 1, color = "red") +
  geom_vline(loneliness_descriptives, mapping = aes(xintercept = plus), size = 1, color = "blue") +
  geom_vline(loneliness_descriptives, mapping = aes(xintercept = minus), size = 1, color = "blue") +
  labs(title = "Distribution of loneliness scores",
       subtitle = "Mean values (red) and 1 standard deviation below and above the mean (blue)",
       x = "Composite loneliness score ranging from 1 (never) to 4 (often)",
       y = "Count") +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times"),
        strip.text.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) 

ggsave(filename = "loneliness_histograms.jpeg",
       path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
       width = 6, 
       height = 9,  
       bg="white",
       dpi=700)


# histogram of BMI scores 
bmi_data <- UiN_data_analysis %>% 
  select(bmi_2:bmi_4) %>% 
  pivot_longer(names_to = "time_wave",
               values_to = "score",
               cols = bmi_2:bmi_4) %>% 
  mutate(time_wave = case_when(time_wave == "bmi_2" ~ "T1",
                               time_wave == "bmi_3" ~ "T2",
                               time_wave == "bmi_4" ~ "T3"))

bmi_descriptives <- bmi_data %>% 
  group_by(time_wave) %>% 
  summarise(mean = mean(score, na.rm = T),
            plus = mean+sd(score, na.rm = T),
            minus = mean-sd(score, na.rm = T))

bmi_histograms_prespec <- ggplot(bmi_data, aes(x = score)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ time_wave, scales = "free", ncol = 1) +
  geom_vline(bmi_descriptives, mapping = aes(xintercept = mean), size = 1, color = "red") +
  geom_vline(bmi_descriptives, mapping = aes(xintercept = plus), size = 1, color = "blue") +
  geom_vline(bmi_descriptives, mapping = aes(xintercept = minus), size = 1, color = "blue") +
  labs(title = "Distribution of Body Mass Index",
       subtitle = "Mean values (red) and 1 standard deviation below and above the mean (blue)",
       x = "Body Mass Index",
       y = "Count") +
  scale_x_continuous(limits = c(15, 40)) +
  theme_minimal() +
  theme(text = element_text(size = 12, family = "Times"),
        strip.text.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) 

  ggsave(filename = "bmi_histograms_prespef.jpeg",
         path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
         width = 6, 
         height = 9,  
         bg="white",
         dpi=900)
  
  bmi_histograms_sd <- ggplot(bmi_data, aes(x = score)) +
    geom_histogram(bins = 20) +
    facet_wrap(~ time_wave, scales = "free", ncol = 1) +
    geom_vline(xintercept = 17, size = 1, color = "orange") +
    geom_vline(xintercept = 30, size = 1, color = "orange") +
    geom_vline(xintercept = 21.75, size = 1, color = "darkgreen") +
    labs(title = "Distribution of Body Mass Index",
         subtitle = "Normative value (gree) and underweight and obesity (orange)",
         x = "Body Mass Index",
         y = "Count") +
    scale_x_continuous(limits = c(15, 40)) +
    theme_minimal() +
    theme(text = element_text(size = 12, family = "Times"),
          strip.text.x = element_text(size = 14, face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) 
  
  ggsave(filename = "bmi_histograms_sd.jpeg",
         path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
         width = 6, 
         height = 9,  
         bg="white",
         dpi=900)
  
  
alcohol_histograms <- UiN_data_analysis %>% 
  select(alcohol_use_2:alcohol_use_4) %>% 
  pivot_longer(names_to = "time_wave",
               values_to = "score",
               cols = alcohol_use_2:alcohol_use_4) %>% 
  mutate(time_wave = case_when(time_wave == "alcohol_use_2" ~ "T1",
                               time_wave == "alcohol_use_3" ~ "T2",
                               time_wave == "alcohol_use_4" ~ "T3")) %>%
  filter(score <= 80) %>% 
    ggplot(aes(x = score)) +
    geom_histogram() +
    coord_cartesian(xlim = c(0, 80)) +
    facet_wrap(~ time_wave, scales = "free", ncol = 1) +
    geom_vline(xintercept = 0, size = 1, color = "darkgreen") +
    geom_vline(xintercept = 60, size = 1, color = "red") +
    geom_vline(xintercept = 10, size = 1, color = "orange") +
    labs(title = "Distribution of Alcohol Use",
         subtitle = "Abstinence (green), low/occasional use (orange) and heavy use (red)",
         x = "Drinks per 4 weeks",
         y = "Count") +
    theme_minimal() +
    theme(text = element_text(size = 12, family = "Times"),
          strip.text.x = element_text(size = 14, face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) 
  
  ggsave(filename = "alcohol_histograms.jpeg",
         path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
         width = 6, 
         height = 9,  
         bg="white",
         dpi=900)
  
  
exercise_histograms <- UiN_data_analysis %>% 
    select(phys_exer_2:phys_exer_4) %>% 
    pivot_longer(names_to = "time_wave",
                 values_to = "score",
                 cols = phys_exer_2:phys_exer_4) %>% 
    mutate(time_wave = case_when(time_wave == "phys_exer_2" ~ "T1",
                                 time_wave == "phys_exer_3" ~ "T2",
                                 time_wave == "phys_exer_4" ~ "T3")) %>%
    filter(score <= 400) %>% 
    ggplot(aes(x = score)) +
    geom_histogram(bins = 12) +
    coord_cartesian(xlim = c(0, 400)) +
    facet_wrap(~ time_wave, scales = "free", ncol = 1) +
    geom_vline(xintercept = 0, size = 1, color = "red") +
    geom_vline(xintercept = 150, size = 1, color = "darkgreen") +
    labs(title = "Distribution of Physical Exercise",
         subtitle = "Inactivity (red) and moderate exercise (green)",
         x = "Minutes per week/number of activities",
         y = "Count") +
    theme_minimal() +
    theme(text = element_text(size = 12, family = "Times"),
          strip.text.x = element_text(size = 14, face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) 
  
  ggsave(filename = "exercise_histograms.jpeg",
         path = "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Graphs", 
         width = 6, 
         height = 9,  
         bg="white",
         dpi=900)
  
  
  
  


### MISSING DATA EXPLORATION ###

missing_proportions_long <- UiN_data_analysis %>% 
  select(-c(id, starts_with("loneliness_"))) %>% 
  summarise(across(everything(), ~ round(mean(is.na(.))*100))) %>% 
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value")  

missing_proportions_long %>% arrange(desc(value))

missing_proportions_table <- missing_proportions_long %>% 
  mutate(variable = ifelse(grepl("\\d$", variable), variable, paste0(variable, "_1"))) %>%
  separate(variable, into = c("name", "wave"), sep = "_(?=\\d+$)") %>% 
  pivot_wider(names_from = "wave",
              values_from = "value") 
write.csv(missing_proportions_table, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Tables/missing_props")


missing_proportions_long %>% 
  summarise(mean_NA_proportion = mean(value)) 
# on average, 10% of values are missing -> number of recommended imputed datasets will be 10

nrow(na.omit(UiN_data_analysis))






### MULTIPLE IMPUTATION ###

# dry run for specifying predictors and imputation methods
UiN_data_analysis <- select(UiN_data_analysis, -c(loneliness_1:loneliness_5))
init_imp <- mice(UiN_data_analysis, maxit = 0)

# predictor matrix
predictor_matrix <- init_imp$predictorMatrix
predictor_matrix[,"id"] <- 0
predictor_matrix["id",] <- 0
# we are including all variables in imputation models

# imputation methods
imp_method <- init_imp$method
# predictive mean matching (pmm) for continuous data,
# logistic regression (logreg) for binary, and
# polytomous regression (polyreg) for un-ordered categorical, 

# MULTIPLE IMPUTATION
set.seed(12345)
Imp1 <- mice(UiN_data_analysis, 
             m = 10, # number of imputed datasets
             maxit = 5, # number of iterations
             method = imp_method, 
             predictorMatrix = predictor_matrix)

# convergence diagnostics
plot(Imp1) 

imp_UiN_data_hbcontext_ext <-  complete(Imp1, "long") %>% 
  mutate(loneliness_1 = rowMeans(select(., UCLA1_1:UCLA1_4), na.rm = T),
         loneliness_2 = rowMeans(select(., UCLA2_1:UCLA2_4), na.rm = T),
         loneliness_3 = rowMeans(select(., UCLA3_1:UCLA3_4), na.rm = T),
         loneliness_4 = rowMeans(select(., UCLA4_1:UCLA4_4), na.rm = T),
         loneliness_5 = rowMeans(select(., UCLA5_1:UCLA5_4), na.rm = T))

# saving data: long-format
write.csv(imp_UiN_data_hbcontext_ext, "N:/durable/Data_analyses/Libor/Loneliness_Health_Risk_Factors/Data/imp_UiN_data_hbcontext_ext")





# SCALE RELIABILITY

outcome_4i <- ' outcome =~ item_1 + item_2 + item_3 + item_4 '

loneli_reliability <- UiN_data_orig %>% 
  filter(!is.na(id), Type2_1 %in% c(1:8), partici1 == 1 & partici2 == 1 & (partici3 == 1 | partici4 == 1)) %>%  # sample: n = 3072
  select(id, starts_with("UCLA") & !ends_with("5") & !ends_with(c("n", "m"))) %>%  # UCLA loneliness scale
  mutate(across(c(starts_with("UCLA") & ends_with(c("1", "2"))),
                ~ 5 - .x)) %>%    # reverse coding items
  pivot_longer(cols = -id,
               names_to = c("time_wave", "item"),
               names_sep = "_",
               values_to = "symptoms") %>% 
  pivot_wider(names_from = "item",
              values_from = "symptoms") %>% 
  filter(time_wave %in% c("UCLA2", "UCLA3", "UCLA4", "UCLA5")) %>% 
  mutate(item_1 = `1`,
         item_2 = `2`,
         item_3 = `3`,
         item_4 = `4`) %>% 
  group_by(time_wave) %>% 
  nest() %>% 
  mutate(model = map(.x = data,
                     ~ cfa(model = outcome_4i,
                           data = .x,
                           estimator = "MLR", # estimator: robust maximum likelihood
                           missing = "ML")), # full-information likelihood for missing values
         omega = map_dbl(.x = model,
                               ~ semTools::reliability(.x)["omega",])) # McDonaldb

           