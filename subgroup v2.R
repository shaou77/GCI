library(survival)
data5 <- readRDS('data5.rds')
#cabg
group_cabg <- data5[data5$cabg == 1, ]
group_nocabg <- data5[data5$cabg == 0 , ]
table(group_cabg$death_365day)
table(group_nocabg$death_365day)
surv_obj_cabg <- Surv(group_cabg$surtime_from_icu, group_cabg$death_365day)
cabg_cox_model1 <- coxph(surv_obj_cabg ~ group, data = group_cabg)
summary(cabg_cox_model1)
surv_obj_nocabg <- Surv(group_nocabg$surtime_from_icu, group_nocabg$death_365day)
nocabg_cox_model1 <- coxph(surv_obj_nocabg ~ group, data = group_nocabg)
summary(nocabg_cox_model1)
cabg_model_no_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ cabg + group, data = data5)
cabg_model_with_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ cabg * group, data = data5)
cabg_lrt <- anova(model_no_interaction, model_with_interaction, test = "LRT")
#gender
group_male <- data5[data5$gender == 'M', ]
group_female <- data5[data5$gender == 'F' , ]
# Male subgroup
surv_obj_male <- Surv(group_male$surtime_from_icu, group_male$death_365day)
male_cox_model <- coxph(surv_obj_male ~ group, data = group_male)
summary(male_cox_model)
# Female subgroup
surv_obj_female <- Surv(group_female$surtime_from_icu, group_female$death_365day)
female_cox_model <- coxph(surv_obj_female ~ group, data = group_female)
summary(female_cox_model)
# Test for interaction between gender and group
gender_model_no_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ gender + group, data = data5)
gender_model_with_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ gender * group, data = data5)
gender_lrt <- anova(gender_model_no_interaction, gender_model_with_interaction, test = "LRT")
#age
group_younger <- data5[data5$admission_age < 70, ]
surv_obj_younger <- Surv(group_younger$surtime_from_icu, group_younger$death_365day)
younger_cox_model <- coxph(surv_obj_younger ~ group, data = group_younger)
summary(younger_cox_model)
# Subgroup of patients aged 70 or older
group_older <- data5[data5$admission_age >= 70, ]
surv_obj_older <- Surv(group_older$surtime_from_icu, group_older$death_365day)
older_cox_model <- coxph(surv_obj_older ~ group, data = group_older)
summary(older_cox_model)
# Test for interaction between admission_age (as a binary variable) and group
data5$age_category <- ifelse(data5$admission_age < 70, 'younger', 'older')
data5$age_category <- factor(data5$age_category)
age_model_no_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ age_category + group, data = data5)
age_model_with_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ age_category * group, data = data5)
age_lrt <- anova(age_model_no_interaction, age_model_with_interaction, test = "LRT")
# Subgroup analysis for vasoactive
group_vasoactive <- data5[data5$vasoactive == 1, ]
surv_obj_vasoactive <- Surv(group_vasoactive$surtime_from_icu, group_vasoactive$death_365day)
vasoactive_cox_model <- coxph(surv_obj_vasoactive ~ group, data = group_vasoactive)
summary(vasoactive_cox_model)
group_novasoactive <- data5[data5$vasoactive == 0, ]
surv_obj_novasoactive <- Surv(group_novasoactive$surtime_from_icu, group_novasoactive$death_365day)
novasoactive_cox_model <- coxph(surv_obj_novasoactive ~ group, data = group_novasoactive)
summary(novasoactive_cox_model)
vasoactive_model_no_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ vasoactive + group, data = data5)
vasoactive_model_with_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ vasoactive * group, data = data5)
vasoactive_lrt <- anova(vasoactive_model_no_interaction, vasoactive_model_with_interaction, test = "LRT")
# Subgroup analysis for congestive_heart_failure
group_chf <- data5[data5$congestive_heart_failure == 1, ]
surv_obj_chf <- Surv(group_chf$surtime_from_icu, group_chf$death_365day)
chf_cox_model <- coxph(surv_obj_chf ~ group, data = group_chf)
summary(chf_cox_model)
group_nochf <- data5[data5$congestive_heart_failure == 0, ]
surv_obj_nochf <- Surv(group_nochf$surtime_from_icu, group_nochf$death_365day)
nochf_cox_model <- coxph(surv_obj_nochf ~ group, data = group_nochf)
summary(nochf_cox_model)
chf_model_no_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ congestive_heart_failure + group, data = data5)
chf_model_with_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ congestive_heart_failure * group, data = data5)
chf_lrt <- anova(chf_model_no_interaction, chf_model_with_interaction, test = "LRT")
# Subgroup analysis for hba1c
group_hba1c_low <- data5[data5$hba1c < 7, ]
surv_obj_hba1c_low <- Surv(group_hba1c_low$surtime_from_icu, group_hba1c_low$death_365day)
hba1c_low_cox_model <- coxph(surv_obj_hba1c_low ~ group, data = group_hba1c_low)
summary(hba1c_low_cox_model)
group_hba1c_high <- data5[data5$hba1c >= 7, ]
surv_obj_hba1c_high <- Surv(group_hba1c_high$surtime_from_icu, group_hba1c_high$death_365day)
hba1c_high_cox_model <- coxph(surv_obj_hba1c_high ~ group, data = group_hba1c_high)
summary(hba1c_high_cox_model)
data5$hba1c_group <- ifelse(data5$hba1c < 7, 'low', 'high')
data5$hba1c_group <- factor(data5$hba1c_group)
hba1c_model_no_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ hba1c_group + group, data = data5)
hba1c_model_with_interaction <- coxph(Surv(data5$surtime_from_icu, data5$death_365day) ~ hba1c_group * group, data = data5)
hba1c_lrt <- anova(hba1c_model_no_interaction, hba1c_model_with_interaction, test = "LRT")
