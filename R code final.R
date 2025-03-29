library(readr)
library(dplyr)
library(rcssci)
library(naniar)
library(gtsummary)
library(mice)
library(survival)
library(survminer)
library(glmnet)
library(timeROC)
library(stats)
#载入数据，计算比值
data <- read_csv("data.csv", col_types = cols(gender = col_factor(levels = c("F", "M")), 
                                              hospital_expire_flag = col_factor(levels = c("0","1")), 
                                              hypoglycemia = col_factor(levels = c("0", "1"))))
data2 <- data %>%   mutate(dialysis_active = ifelse(is.na(dialysis_active), 0, dialysis_active),
                           surtime_from_admit= ifelse(is.na(surtime_from_admit), 366, surtime_from_admit),
                           surtime_from_icu= ifelse(is.na(surtime_from_icu), 366, surtime_from_icu),
                           death_365day = ifelse(surtime_from_icu<=365,1,0),
                           medium_gs_mmol=(round(medium_gs/18,2)),
                           mean_gs_mmol=(round(avg_glucose/18,2)),
                           index=round(medium_gs/(28.7*hba1c-46.7),4)*100,
                           GCI=round(avg_glucose/(28.7*hba1c-46.7),4)*100) %>% 
  mutate(surtime_from_admit = ifelse(surtime_from_admit > 365, 366, surtime_from_admit),
         surtime_from_icu = ifelse(surtime_from_icu > 365, 366, surtime_from_icu))
rcssci_cox(data=data2, y = "death_365day",x = "GCI",time="surtime_from_icu",prob=0.1,filepath= "output")
##Time-dependent ROC curve estimation
roc2 <- timeROC(T=data2$surtime_from_icu,delta=data2$death_365day,marker=data2$GCI,cause=1,weighting="marginal",times=c(30,180,365),iid=TRUE)
plot(roc2, time=30, lwd=2, col = "#3688f4", title=FALSE,add = FALSE)
plot(roc2, time=180, lwd=2,col = "#f43688", add =TRUE)
plot(roc2, time=365,lwd=2, col = "#88f436", add =TRUE)
title(main="Time Dependent ROC")
legend("bottomright",paste0("AUC of ",c(1,6,12)," month survival:",format(round(roc2$AUC,3),nsmall=2)), col=c("#3688f4","#f43688","#88f436"),lty=1,lwd=2,bty = "n")
# 3分位分组
quantiles <- quantile(data2$GCI, probs = c(1/3,2/3))
data2$group <- cut(data2$GCI,
                   breaks = c(-Inf, quantiles[1],quantiles[2], Inf),
                   labels = c("Low GCI", 'Medium GCI',"High GCI"),
                   include.lowest = TRUE)
data2$group <- factor(data2$group)
data2$group <- relevel(data2$group, ref = "Medium GCI")
cat("Low GCI range: ", "-Inf to", quantiles[1], "\n")
cat("Medium GCI range: ", quantiles[1], "to", quantiles[2], "\n")
cat("High GCI range: ", quantiles[2], "to Inf", "\n")
quantiles2 <- quantile(data2$index, probs = c(1/3,2/3))
data2$group2 <- cut(data2$index,
                   breaks = c(-Inf, quantiles[1],quantiles[2], Inf),
                   labels = c("Low mGCI", 'Medium mGCI',"High mGCI"),
                   include.lowest = TRUE)
data2$group2 <- factor(data2$group2)
data2$group2 <- relevel(data2$group2, ref = "Medium mGCI")
cat("Low mGCI range: ", "-Inf to", quantiles2[1], "\n")
cat("Medium mGCI range: ", quantiles2[1], "to", quantiles2[2], "\n")
cat("High mGCI range: ", quantiles2[2], "to Inf", "\n")
#绘制分布图
library(ggplot2)
# first calculate the median values for each group
medians <- aggregate(cbind(avg_glucose, hba1c) ~ group, data = data2, median)
# Create a scatter plot
ggplot(data2, aes(y = avg_glucose, x = hba1c, color = group)) +
  geom_point(alpha = 0.5) + # Set transparency to see overlapping points
  geom_point(data = medians, aes(y = avg_glucose, x = hba1c), color = "black", size = 3, shape = 17) +
  scale_color_manual(values = c("Low GCI" = "#36f443", "Medium GCI" = "#4336f4", "High GCI" = "#f44336")) +
  theme_minimal() + ylim(50,350)+ xlim(4,16)+
  labs( y = "Average Glucose Level (mg/dL)",
       x = "Glycated Hemoglobin (HbA1c)(%) ",
       color = "GCI Group")
# Create a scatter plot with updated legend labels for 1-year outcome
ggplot(data2, aes(x = hba1c, y = avg_glucose, color = group, shape = factor(death_365day))) +
  geom_point(alpha = 0.7, size = 3) + 
  geom_point(data = medians, aes(x = hba1c, y = avg_glucose), color = "black", size = 4, shape = 17) + # Median points are larger
  scale_color_manual(values = c("Low GCI" = "#8ecfc9", "Medium GCI" = "#ffbe7a", "High GCI" = "#fa7f6f")) +
  scale_shape_manual(values = c("0" = 21, "1" = 16), labels = c("Survived", "Deceased")) +  # Assign labels to shapes
  theme_bw() +
  ylim(50, 350) +
  xlim(4, 14) +
  labs(y = "Average Glucose Level (mg/dL)",
       x = "Glycated Hemoglobin (HbA1c) (%)",
       color = "GCI Group",
       shape = "1-Year Outcome")  # Update legend title
#缺失值处理，10%以上缺失值变量去除。
print(miss_var_summary(data2), n=40)
miss_summary <- miss_var_summary(data2)
cols_to_keep <- miss_summary %>%  filter(pct_miss <= 10) %>%  pull(variable)
data2_filter <- data2 %>% dplyr::select(all_of(cols_to_keep))
# 转换种族，更改因子变量类型
data3 <- data2_filter %>%
  mutate(race = case_when(
    grepl("ASIAN", race) ~ "ASIAN",
    grepl("WHITE", race) ~ "WHITE",
    grepl("BLACK|AFRICAN", race) ~ "BLACK/AFRICAN AMERICAN",
    grepl("HISPANIC|LATINO", race) ~ "HISPANIC OR LATINO",
    TRUE ~ "OTHER"  ),glucose_mean_mmol = round(glucose_mean / 18,2),
    pci = ifelse(is.na(pci), 0, pci),cabg = ifelse(is.na(cabg), 0, cabg))%>% 
  dplyr::select(-stay_id,-heart_rate_mean,-mbp_mean,-resp_rate_mean,-spo2_mean,-weight,-dementia,-aids,-rheumatic_disease,-liver_disease,-gcs_min,-glucose_mean,-chloride, -sofa, -oasis,) %>% 
  rename(median_gs = medium_gs) %>%
  mutate(across(c(ventilation, vasoactive, congestive_heart_failure,
                  cerebrovascular_disease,dialysis_active, chronic_pulmonary_disease,pci,cabg,
                  renal_disease, cancer, death_365day, race), factor))

# 绘制table1
table1 <- data3 %>% tbl_summary(by = group, missing_text = "(Missing)",digits = list(all_continuous() ~ 2)) %>% 
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value=TRUE))
as_gt(table1) %>% gt::gtsave('table1.docx')
mtable1 <- data3 %>% tbl_summary(by = group2, missing_text = "(Missing)",digits = list(all_continuous() ~ 2)) %>% 
  add_p(test.args = all_tests("fisher.test") ~ list(simulate.p.value=TRUE))
as_gt(mtable1) %>% gt::gtsave('mtable1.docx')
#插补缺失值
imputed_data <- mice(data3,seed = 8177)
data5 <- complete(imputed_data, 3)
#开始生存分析
data4 <- data3
data4$death_365day <- as.numeric(data4$death_365day) - 1
surv_obj <- Surv(data4$surtime_from_icu, data4$death_365day)
#K-M曲线
km_fit1 <- survfit(surv_obj ~ group, data = data4)
ggsurvplot(km_fit1, data = data4, pval = TRUE, risk.table = "abs_pct", risk.table.y.text.col = T, risk.table.y.text = FALSE,
           conf.int = F, xlim=c(0,360),break.time.by = 50,
           palette = c('#ffbe7a','#8ecfc9','#fa7f6f'),
           ggtheme = theme_light())
surv_diff <- survdiff(Surv(data4$surtime_from_icu, data4$death_365day) ~ group, data = data4)
surv_diff
#敏感性K-M
km_fit2 <- survfit(surv_obj ~ group2, data = data4)
ggsurvplot(km_fit2, data = data4, pval = TRUE, risk.table = "abs_pct", risk.table.y.text.col = T, risk.table.y.text = FALSE,
           conf.int = F, xlim=c(0,360),break.time.by = 50,
           palette = c('#4336f4','#36f443','#f44336'),
           ggtheme = theme_light())
surv_diff <- survdiff(Surv(data4$surtime_from_icu, data4$death_365day) ~ group2, data = data4)
surv_diff
#cox3个模型
cox_model1 <- coxph(surv_obj ~ group, data = data4)
summary(cox_model1)
cox_model2 <- coxph(surv_obj ~ group+pci+cabg+cancer, data = data4)
summary(cox_model2)
cox_model3 <- coxph(surv_obj ~ apsiii+congestive_heart_failure+group+pci+cabg+cancer, data = data4)
summary(cox_model3)
#敏感性cox3个模型
cox_model1m <- coxph(surv_obj ~ group2, data = data4)
summary(cox_model1m)
cox_model2m <- coxph(surv_obj ~ group2+pci+cabg+cancer, data = data4)
summary(cox_model2m)
cox_model3m <- coxph(surv_obj ~ apsiii+congestive_heart_failure+group2+pci+cabg+cancer, data = data4)
summary(cox_model3m)
#cox第4个模型
data5$death_365day <- as.numeric(data5$death_365day) - 1
surv_obj2 <- Surv(data5$surtime_from_icu, data5$death_365day)
# Define the full Cox model with potential predictors
full_model <- coxph(surv_obj2 ~ pci + cabg + pt + ptt + urineoutput + rdw + potassium + gender + admission_age + race + hematocrit + hemoglobin + platelets + wbc + aniongap + bicarbonate + bun + creatinine + sodium + dialysis_active + ventilation + vasoactive + congestive_heart_failure + cerebrovascular_disease + chronic_pulmonary_disease + renal_disease + cancer + apsiii + group, data = data5)
# Perform stepwise selection based on AIC
stepwise_selected_model <- step(full_model, direction = "both")
summary(stepwise_selected_model)
#另存插补后的data5亚组分析用
saveRDS(data5,file = 'subgroup/data5.rds')
#分段回归
# Assuming data5 is a data frame
group1 <- data5[data5$GCI < 98.54, ]
group2 <- data5[data5$GCI >= 98.54, ]

# Calculate GCI2 for both groups
group1$GCI2 <- group1$GCI / 100
group2$GCI2 <- group2$GCI / 100

# Create data frames explicitly including the desired variables
group1_df <- data.frame(surtime_from_icu = group1$surtime_from_icu, death_365day = group1$death_365day, GCI2 = group1$GCI2)
group2_df <- data.frame(surtime_from_icu = group2$surtime_from_icu, death_365day = group2$death_365day, GCI2 = group2$GCI2)

# Fit the Cox proportional hazards models
cox_model_group1 <- coxph(Surv(group1_df$surtime_from_icu, group1_df$death_365day) ~ GCI2, data = group1_df)
cox_model_group2 <- coxph(Surv(group2_df$surtime_from_icu, group2_df$death_365day) ~ GCI2, data = group2_df)

#次要终点，院内死亡率。
data6 <- complete(imputed_data, 3)
# Fit the initial full model without 'group'
initial_model <- glm(hospital_expire_flag ~ pci + cabg + pt + ptt + urineoutput + rdw + 
                       potassium + gender + admission_age + race + hematocrit + hemoglobin + 
                       platelets + wbc + aniongap + bicarbonate + bun + creatinine + sodium + 
                       dialysis_active + ventilation + vasoactive + congestive_heart_failure + 
                       cerebrovascular_disease + chronic_pulmonary_disease + renal_disease + 
                       cancer + apsiii,
                     family = binomial, data = data6)

# Perform stepwise selection
stepwise_model <- step(initial_model, direction = "both")
# Add 'group' back to the model if it was removed
final_model <- update(stepwise_model, . ~ . + group)
# Summary of the final model
summary(final_model)
#次要终点，低血糖发生。
full_model2 <- glm(hypoglycemia ~ pci + cabg + pt + ptt + urineoutput + rdw + 
                    potassium + gender + admission_age + race + hematocrit + hemoglobin + 
                    platelets + wbc + aniongap + bicarbonate + bun + creatinine + sodium + 
                    dialysis_active + ventilation + vasoactive + congestive_heart_failure + 
                    cerebrovascular_disease + chronic_pulmonary_disease + renal_disease + 
                    cancer + apsiii + group, 
                  family = binomial, data = data6)
# Perform stepwise regression using AIC for variable selection
stepwise_model2 <- step(full_model2, direction = "both")
# Summary of the final stepwise model
summary(stepwise_model2)
#次要终点，ICU free days
data6 <- data6 %>%
  mutate(icufreeday = case_when(
    is.na(surtime_from_icu) ~ 30 - los_icu, 
    (surtime_from_icu <= 30 & hospital_expire_flag == 1) ~ 0,
    (abs(surtime_from_icu - los_icu) <= 1 & hospital_expire_flag == 1 & los_icu <= 28) ~ 0,
    TRUE ~ pmax(0, 30 - los_icu)
  ))
full_model3 <- glm(icufreeday ~ pci + cabg + pt + ptt + urineoutput + rdw + potassium + gender + admission_age + race + hematocrit + hemoglobin + platelets + wbc + aniongap + bicarbonate + bun + creatinine + sodium + dialysis_active + ventilation + vasoactive + congestive_heart_failure + cerebrovascular_disease + chronic_pulmonary_disease + renal_disease + cancer + apsiii + group, data = data6, family = gaussian())
reduced_model <- step(full_model3, direction = "both")
summary(reduced_model)


