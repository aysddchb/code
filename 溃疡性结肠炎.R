library(haven)#读取XPT文件
library(plyr)
library(dplyr)
library(arsenal)#下载用于数据快速预览的包 arsenal，使用 tableby 函数
library(survey)#用于加权情况下的分析
setwd("F://work/101")#设置工作路径，等会可直接读取这个文件夹里的数据


#人口学数据读取，后面的数据相同

demo.f<-read_xpt("DEMO_F.XPT")#读取你的人口学数据，有几个周期就读取几个
colnames(demo.f)#查看变量名
tab1 <- tableby( ~ RIDAGEYR + factor(RIAGENDR) +factor( RIDRETH1) +WTINT2YR +WTMEC2YR+SDMVPSU+SDMVSTRA, 
                 data=demo.f)
summary(tab1, text=TRUE)#查看你提取的变量，并进行简单的描述表，确认提取是否有错误
demo.data<-demo.f[,c('SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH1', 'INDFMPIR', 'DMDEDUC2','DMDMARTL', 'SDMVPSU','SDMVSTRA' ,'WTMEC2YR')]#提取需要的变量


demo.data <- demo.data %>%
  mutate(DMDEDUC2 = case_when(
    DMDEDUC2 %in% c(1, 2) ~ 1,
    DMDEDUC2 == 3 ~ 2,
    DMDEDUC2 %in% c(4, 5) ~ 3,
    TRUE ~ DMDEDUC2  # 对于不符合上述条件的值，保持原值不变
  ))

demo.data <- demo.data %>%
  mutate(DMDMARTL = case_when(
    DMDMARTL %in% c(1, 6) ~ 1,
    DMDMARTL %in% c(2, 3, 4) ~ 2,
    DMDMARTL == 5 ~ 3,
    TRUE ~ DMDMARTL  # 对于不符合上述条件的值，保持原值不变
  ))

#这个是简单的一个周期的提取，这个是两个的，同理
#demo.p<-read_xpt("20/P_DEMO.XPT")
#demo.i<-read_xpt("16/DEMO_I.XPT")
#colnames(demo.p)
#colnames(demo.i)
#tab1 <- tableby( ~ RIDAGEYR + factor(RIAGENDR) +factor( RIDRETH1) +WTINTPRP +WTMECPRP+SDMVPSU+SDMVSTRA, 
#data=demo.p)这一步其实要不要都行，只是确认下提取的没错，数据没有少
#summary(tab1, text=TRUE)
#demo.data.file <- dplyr::bind_rows(list(demo.i, demo.p))按行上下连接
#dim(demo.data.file)
#demo.data <- demo.data.file[,c('SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH1', 'SDMVPSU','SDMVSTRA' )]


#结肠炎提取
arq.f<-read_xpt("ARQ_F.XPT")
arq.data<-arq.f[,c('SEQN','ARQ125C')]

#hypertension
bpq.f<-read_xpt("BPQ_F.XPT")
bpq.data<-bpq.f[,c('SEQN','BPQ020','BPQ030','BPQ057')]
bpq.data <- bpq.data %>%mutate(BPQ020 = ifelse(BPQ020 == 2, 0, BPQ020))
bpq.data <- bpq.data %>%mutate(BPQ030 = ifelse(BPQ030 == 2, 0, BPQ030))
bpq.data <- bpq.data %>%mutate(BPQ057 = ifelse(BPQ057 == 2, 0, BPQ057))


#Diabetes
diq.f<-read_xpt("DIQ_F.XPT")
diq.data<-diq.f[,c('SEQN','DIQ010')]
diq.data <- diq.data %>%mutate(DIQ010 = ifelse(DIQ010 == 2, 0, DIQ010))

#BMI
bmx.f<-read_xpt("BMX_F.XPT")
bmx.data<-bmx.f[,c('SEQN','BMXBMI')]#每一步都要带上SEQN，相当于你的ID CARD
#SMQ
smq.f<-read_xpt("SMQ_F.XPT")
smq.data<-smq.f[,c('SEQN','SMQ040')]

smq.data <- smq.data %>%
  mutate(SMQ040 = case_when(
    SMQ040 %in% c(1, 2) ~ 1,
    SMQ040 == 3 ~ 0,
    TRUE ~ SMQ040  # 对于不符合上述条件的值，保持原值不变
  ))
#ALQ
alq.f<-read_xpt("ALQ_F.XPT")
alq.data<-alq.f[,c('SEQN','ALQ101')]
alq.data <- alq.data %>%mutate(ALQ101 = ifelse(ALQ101 == 2, 0, ALQ101))

#葡萄糖，白蛋白
biopro.f<-read_xpt("BIOPrO_F.XPT")
biopro.data<-biopro.f[,c('SEQN','LBXSGL', 'LBXSAL')]

#CRP
crp.f<-read_xpt("CRP_F.XPT")
crp.data<-crp.f[,c('SEQN','LBXCRP')]
#总胆固醇
tchol.f<-read_xpt("TCHOL_F.XPT")
tchol.data<-tchol.f[,c('SEQN','LBXTC')]



#全连接，根据SEQN进行匹配，左连接，把前面各部分提取的数据综合在一起
data.crude <- plyr::join_all(list(demo.data, bmx.data,bpq.data,diq.data,arq.data,smq.data,alq.data,biopro.data,crp.data,tchol.data),
                          by='SEQN', type='left')

#输出数据用EXCEL清洗数据，R也能做，但是我觉得excel更方便，个人习惯
write.csv(data.crude, "arqdata.csv", row.names = FALSE)



library(tidyverse)
library(gtsummary)
library(tidyr) # drop_na 函数，用于快速去掉 NA
library(survey)
library(haven)

uc.data<-read.csv("arqdatadel.csv")#读取清洗好的数据
#####  生成复杂抽样 NHANES_design ##### 
NHANES_design <- svydesign(
  data = uc.data, 
  ids = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  nest = TRUE, 
  weights = ~WTMEC2YR,
  survey.lonely.psu = "adjust") # 可以加上 survey.lonely.psu = "adjust" 避免1个PSU报错,ids,strata,weights均为权重

summary(NHANES_design)


#高血压1，logistic回归,人口学,加权
m1 <- svyglm(Uc ~ factor(BPQ020)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio, design = NHANES_design,
             family = binomial)

summary(m1)

# 加载所需的包
library(survey)
library(broom)

# 假设 m1 是您的 svyglm 模型

# 计算比值比（Odds Ratios）
or_values <- exp(coef(m1))

# 计算标准误差
std_errors <- summary(m1)$coefficients[, "Std. Error"]

# 计算95%置信区间
conf_int_lower <- exp(coef(m1) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m1) + 1.96 * std_errors)

# 创建一个数据框来展示结果
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

# 输出结果，保留4位小数
print(round(results, 4))


#高血压1，logistic回归,人口学,不加权
m2 <- glm(Uc ~ factor(BPQ020)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio,data = uc.data,
             family = binomial)

summary(m2)
or_values <- exp(coef(m2))
std_errors <- summary(m2)$coefficients[, "Std. Error"]
conf_int_lower <- exp(coef(m2) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m2) + 1.96 * std_errors)
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

print(round(results, 4))

#高血压2，总模型加权
m3 <- svyglm(Uc ~ factor(BPQ020)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio + BMI + factor(Diabetes) + factor(Smoke) + factor(Alcohol) + Glucose + Albumin + Crp, design = NHANES_design,
             family = binomial)

summary(m3)
or_values <- exp(coef(m3))
std_errors <- summary(m3)$coefficients[, "Std. Error"]
conf_int_lower <- exp(coef(m3) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m3) + 1.96 * std_errors)
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

print(round(results, 4))

#高血压，总模型不加
m4 <- glm(Uc ~ factor(BPQ020)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio + BMI + factor(Diabetes) + factor(Smoke) + factor(Alcohol) + Glucose + Albumin + Crp, data = uc.data ,
             family = binomial)

summary(m4)
or_values <- exp(coef(m4))
std_errors <- summary(m4)$coefficients[, "Std. Error"]
conf_int_lower <- exp(coef(m4) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m4) + 1.96 * std_errors)
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

print(round(results, 4))

#临界高血压人口学加权
m5 <- svyglm(Uc ~ factor(BPQ057)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio, design = NHANES_design,
             family = binomial)

summary(m5)
or_values <- exp(coef(m5))
std_errors <- summary(m5)$coefficients[, "Std. Error"]
conf_int_lower <- exp(coef(m5) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m5) + 1.96 * std_errors)
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

print(round(results, 4))
#临界高血压人口学不加权
m6 <- glm(Uc ~ factor(BPQ057)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio,data = uc.data,
          family = binomial)

summary(m6)
or_values <- exp(coef(m6))
std_errors <- summary(m6)$coefficients[, "Std. Error"]
conf_int_lower <- exp(coef(m6) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m6) + 1.96 * std_errors)
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

print(round(results, 4))

#全部加权
m7 <- svyglm(Uc ~ factor(BPQ057)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio + BMI + factor(Diabetes) + factor(Smoke) + factor(Alcohol) + Glucose + Albumin + Crp, design = NHANES_design,
             family = binomial)

summary(m7)
or_values <- exp(coef(m7))
std_errors <- summary(m7)$coefficients[, "Std. Error"]
conf_int_lower <- exp(coef(m7) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m7) + 1.96 * std_errors)
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

print(round(results, 4))

#全部不加权
m8 <- glm(Uc ~ factor(BPQ057)+ Age + factor(Gender) + factor(Ethnicity) + factor(Education) +factor(Marry) + Income.poverty.ratio + BMI + factor(Diabetes) + factor(Smoke) + factor(Alcohol) + Glucose + Albumin + Crp, data = uc.data ,
          family = binomial)

summary(m4)
or_values <- exp(coef(m8))
std_errors <- summary(m8)$coefficients[, "Std. Error"]
conf_int_lower <- exp(coef(m8) - 1.96 * std_errors)
conf_int_upper <- exp(coef(m8) + 1.96 * std_errors)
results <- data.frame(
  Odds_Ratio = or_values,
  Lower_CI_95 = conf_int_lower,
  Upper_CI_95 = conf_int_upper
)

print(round(results, 4))
































data_male <- subset(uc.data, Gender == 1)  # 男性
data_female <- subset(uc.data, Gender == 2)  # 女性

#高血压男女分层回归
m4 <- glm(
  Uc ~factor(Hypertension)+ Year  + factor(Ethnicity) + factor(Diabetes) + BMI, #
  data =data_male,
  family = binomial)
summary(m4)

#年龄分组我用EXCEL分的，年龄亚组

yeard<-read.csv("arqdel.csv")
 year1<-subset(yeard, Yeargrp == 1)  
 year2<-subset(yeard, Yeargrp == 2)  
 year3<-subset(yeard, Yeargrp == 3)#重新提取数据，进行分层  

year3$Hypertension[year3$Hypertension == 2] <- 0#改了高血压编码

m8 <- glm(Uc ~ factor(Hypertension) + factor(Gender)  + factor(Ethnicity) + factor(Diabetes) + BMI, #
          +            data = year3,
          +              family = binomial)