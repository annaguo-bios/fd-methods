############################
# Data Cleaning
############################
library(scales)
library(boot)
library(dplyr)
library(fdcausal)
library(mice)
library(flexCausal)
# Due to data sharing constraints, we are unable to provide direct access to the raw data used for the real data analysis presented in this study. 
# However, the dataset is available for application through the Finnish Social Science Data Archive at 
# https://services.fsd.tuni.fi/catalogue/FSD2076?tab=variables&lang=en&study_language=en. 

dat <- read.csv2("data.csv") # data for analysis

# subset variables of potential interests
dat <- dat[,c("t3","bv4_1","ktu32","t10","t11","t12","t13","t14","l16","l17","k70","k71","k72","k73","k74",
              "k75","k76","k77","koulu7_1","koulu1","koulu4","koulu7_2","koulu7_9","ktu19","ktu31a_4","ktu31a11","ktu31a28","ktu31a29","k102")]

dat <- dat%>% rename(
  sex=t3,
  income2000=bv4_1,
  income1991=ktu32,
  father_income1972=t10,
  father_income1983=t11,
  mother_income1972=t12,
  mother_income1983=t13,
  familySES=t14,
  totalITPA_10=l16,
  totalITPA_12=l17,
  grade_2yr=k70,
  grade_3yr=k71,
  grade_4yr=k72,
  grade_5yr=k73,
  grade_6yr=k74,
  grade_7yr=k75,
  grade_8yr=k76,
  grade_9yr=k77,
  highest_edu=koulu7_1,
  len_of_edu=koulu1,
  num_field=koulu4,
  edu_field=koulu7_2,
  age_start_highedu=koulu7_9,
  qual_job=ktu19,
  len_unemp=ktu31a_4,
  len_gap=ktu31a11,
  age_work=ktu31a28,
  age1991=ktu31a29,
  GPA=k102
) %>% mutate( edu_field = ifelse(edu_field==9, NA, edu_field), edu_field=factor(ifelse(edu_field %in% c(0,1,2),0,1), levels=0:1, labels=c("Art","Science")),
              qual_job=factor(qual_job,levels=0:2,labels=c("No","Somewhat","Yes"),ordered = T),
              familySES=factor(familySES,levels=1:3,labels=c("High","Medium","low"),ordered = T),
              family.income=father_income1972 +father_income1983+ mother_income1972 +mother_income1983,
              sex=factor(sex,levels=c(1,2), labels=c("Male","Female")))

# Selected variable for analysis
dt1 <- dat %>% select(-grade_2yr,-grade_3yr,-grade_4yr,-grade_5yr,-grade_7yr,-grade_8yr,-grade_9yr,
                      -income1991,-father_income1972,-father_income1983,-mother_income1972,-mother_income1983, -totalITPA_12,-GPA,-len_gap,-familySES)

# Write the data frame to a R object
save(list= c("dt1"), file = "data.Rdata")



############################
# Data Analysis
############################

# Load data
load("data.Rdata")

# Estimation
set.seed(7)
dt1 <- mice(dt1, m=1)
dt1 <- complete(dt1,1) %>% mutate(A=ifelse(grade_6yr>median(grade_6yr),1,0), highest_edu=factor(ifelse(highest_edu<6,0,1), levels=0:1,label=c("No","Yes")),
                                  len_unemp=factor(ifelse(len_unemp<=2,0,1),levels=0:1,label=c("No","Yes")))

sixth.est3.SL <- TMLE(a=c(1,0),data=dt1,treatment="A", mediators=c("highest_edu","len_of_edu","age_start_highedu","num_field","edu_field","qual_job","len_unemp","age_work"),
                      outcome="income2000", covariates=c("family.income","totalITPA_10","sex","age1991"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                      lib = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean","SL.xgboost"))

# set.seed(7)
# tmp <- ADMGtmle(a=c(1,0),data=dt1,vertices=c('A','income2000','X','M'),treatment="A", outcome="income2000",
#          multivariate.variables = list(M=c("highest_edu","len_of_edu","age_start_highedu","num_field","edu_field","qual_job","len_unemp","age_work"),
#                                        X=c("family.income","totalITPA_10","sex","age1991")),
#          bi_edges = list(c("A","income2000")),
#          di_edges = list(c("X","A"),c("X","M"),c("X","income2000"), c('A','M'), c('M','income2000')),
#          superlearner.seq = F,
#          superlearner.Y = T,
#          superlearner.A = T,
#          superlearner.M = T,
#          superlearner.L = T,
#          lib.seq = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean","SL.xgboost"),
#          lib.L = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean","SL.xgboost"),
#          lib.M = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean","SL.xgboost"),
#          lib.Y = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean","SL.xgboost"),
#          lib.A = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean","SL.xgboost"))
# 
# save(list = c("tmp"), file="/Users/apple/Desktop/estimation.Rdata")
save(list = c("sixth.est3.SL"), file="estimation.Rdata")
