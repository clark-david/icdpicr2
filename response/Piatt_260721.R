####################################################
#
#    PROGRAM TO REPLICATE
#        PIATT J, INJ EPIDEMIOL 2026; 13:35
#
#    David Clark
#    July 2026
#
####################################################


#CLEAR WORKSPACE
#LOAD REQUIRED PACKAGES, IF NOT LOADED ALREADY
rm(list=ls())
require(tidyverse)
require(janitor)
require(broom)
require(skimr)
require(modelr)
require(icdpicr2)
require(pROC)

#IMPORT TQIP DATA  
#  (You have to get the raw data from the ACS)

d0<-read_csv("~/PUF2022/PUF_TRAUMA.csv")  

#Drop if outcome unspecified or transferred to another acute-care hospital
#Delete unnecessary fields
d1<-filter(d0,HOSPDISCHARGEDISPOSITION>=2 & HOSPDISCHARGEDISPOSITION<=14)
d1<-select(d1,inc_key,SEX,AgeYears,ASIAN,BLACK,ETHNICITY,PRIMARYMETHODPAYMENT,ISS,HOSPDISCHARGEDISPOSITION,
           RACEOTHER,RACE_NA,RACE_UK,WITHDRAWALLST,FINALDISCHARGEDAYS,PRIMARYECODEICD10)

#Get ICD diagnosis data and combine with above
d2<-read_csv("~/PUF2022/PUF_ICDDIAGNOSIS.csv")
d2<-rename(d2,inc_key=Inc_Key,dx=ICDDIAGNOSISCODE)
d3<-left_join(d1,d2,by=("inc_key"))

#Convert to "wide" format and run ICDPIC-R2
d3 <- group_by(d3,inc_key)
d3 <- mutate(d3,seq=row_number())
d3 <- ungroup(d3)
d4 <- pivot_wider(d3,id_cols=c(inc_key,PRIMARYECODEICD10),names_from=seq,names_prefix="I10_DX",values_from=dx)
d4 <- relocate(d4,PRIMARYECODEICD10, .after=I10_DX50)
d4 <- rename(d4,I10_DX51=PRIMARYECODEICD10)
#The next command will take about 15 minutes to run
d5 <- cat_trauma2(d4,"I10_DX",TRUE)

d6 <- left_join(d1,d5,by="inc_key")
write_csv(d6,"~/PUF2022/cat_trauma_result")

########################################################

d6 <- read_csv("~/PUF2022/cat_trauma_result")

#Identify cases with head injury
dh0<-read_csv("~/PUF2022/PUF_AISDIAGNOSIS.csv")  
dh0<-filter(dh0,str_sub(AISPreDot,1,1)=="1")
dh1<-read_csv("~/PUF2022/PUF_AISDIAGNOSIS_LOOKUP.csv")  
dh1<-rename(dh1,AISPreDot=AISPREDOT,AISVersion=AISVERSION)
dh01<-left_join(dh0,dh1,by=c("AISPreDot","AISVersion"))

#Head injuries acknowledged to lack outcome data
(409187-337148)/409187

#Note that not all selected head injuries are classified in ISS region 1
tabyl(dh01,ISSREGION)

#Apparently Dr. Piatt considered AIS=9 maximum rather than missing
#This error will lead to inappropriate models later on
dh01<-group_by(dh01,inc_key)
dh01<-arrange(dh01,-AISSeverity)
dh01<-mutate(dh01,seq=row_number())
dh01<-mutate(dh01,AISmax=max(AISSeverity))
dh01<-ungroup(dh01)
dh01<-filter(dh01,seq==1)
dh01<-arrange(dh01,inc_key,AISSeverity)
tabyl(dh01,AISmax)
#If you want to replicate Dr. Piatt's numbers, use the above
#If you want to correct the error, use the following
#tabyl(dh01,AISSEVERITY,AISSeverity) # Initially identical
#dh01<-mutate(dh01,AISSEVERITY=if_else(AISSeverity==9,0,AISSEVERITY))
#tabyl(dh01,AISSEVERITY,AISSeverity) # Severity 0 replaces 9 (for proper sorting)
#dh01<-group_by(dh01,inc_key)
#dh01<-arrange(dh01,-AISSEVERITY)
#dh01<-mutate(dh01,seq=row_number())
#dh01<-mutate(dh01,AISmax=max(AISSEVERITY))
#dh01<-ungroup(dh01)
#dh01<-filter(dh01,seq==1)
#dh01<-arrange(dh01,inc_key,AISSEVERITY)
#tabyl(dh01,AISmax)

#Merge patient data with diagnoses
#Replicate Table 1
d7<-inner_join(dh01,d6,by="inc_key")
tabyl(d7,AISmax,mxaisbr_H)
#A couple of cells in the published Table 1 are incorrect


#Create other analytic variables

d7<-mutate(d7,died=if_else(HOSPDISCHARGEDISPOSITION==5,1,0))
d7<-mutate(d7,hospice=if_else(HOSPDISCHARGEDISPOSITION==8,1,0)) 
d7<-mutate(d7,diedorhospice=if_else((died==1 | hospice==1),1,0))   
d7<-mutate(d7,dnr=if_else(WITHDRAWALLST==1,1,0))
d7<-mutate(d7,diedorhospiceordnr=if_else(diedorhospice==1 | dnr==1,1,0))
tabyl(d7,dnr,diedorhospice)

d7<-mutate(d7,ISS2=ISS^2)
d7<-mutate(d7,ISS3=ISS^3)
d7<-mutate(d7,riss2=riss^2)
d7<-mutate(d7,riss3=riss^3)
d7<-mutate(d7,head0=if_else(AISmax==0,1,0))
d7<-mutate(d7,head2=if_else(AISmax==2,1,0))
d7<-mutate(d7,head3=if_else(AISmax==3,1,0))
d7<-mutate(d7,head4=if_else(AISmax==4,1,0))
d7<-mutate(d7,head5=if_else(AISmax==5,1,0))
d7<-mutate(d7,head6=if_else(AISmax==6,1,0))
d7<-mutate(d7,icdhead0=if_else(mxaisbr_H==0,1,0))
d7<-mutate(d7,icdhead2=if_else(mxaisbr_H==2,1,0))
d7<-mutate(d7,icdhead3=if_else(mxaisbr_H==3,1,0))
d7<-mutate(d7,icdhead4=if_else(mxaisbr_H==4,1,0))
d7<-mutate(d7,icdhead5=if_else(mxaisbr_H==5,1,0))
d7<-mutate(d7,age0005=if_else(AgeYears>=0 & AgeYears<=5,1,0))
d7<-mutate(d7,age0612=if_else(AgeYears>=6 & AgeYears<=12,1,0))
d7<-mutate(d7,age1317=if_else(AgeYears>=13 & AgeYears<=17,1,0))
d7<-mutate(d7,age1860=if_else(AgeYears>=18 & AgeYears<=60,1,0))
d7<-mutate(d7,age6199=if_else(AgeYears>=61 & AgeYears<=999,1,0))
d7<-mutate(d7,commercial=if_else(PRIMARYMETHODPAYMENT==4,1,0))
d7<-mutate(d7,otherins=if_else(commercial==0,1,0))
d7<-mutate(d7,other0005=otherins*age0005)
d7<-mutate(d7,other0612=otherins*age0612)
d7<-mutate(d7,other1317=otherins*age1317)
d7<-mutate(d7,other6199=otherins*age6199)
d7<-rename(d7,LOS=FINALDISCHARGEDAYS)
d7<-mutate(d7,logLOS=log(LOS))

tabyl(d7,ETHNICITY)
d7<-mutate(d7,HISPANIC=if_else(ETHNICITY==1,1,0))
d7<-mutate(d7,HISPANIC=if_else(is.na(HISPANIC),0,HISPANIC))
d7<-mutate(d7,ETHNICITY_UK=if_else(is.na(ETHNICITY),1,0))

d7<-mutate(d7,ISScat=case_when(
      ISS>=1 & ISS<=8   ~ "ISS0108",
      ISS>=9 & ISS<=15  ~ "ISS0915",
      ISS>=16 & ISS<=24 ~ "ISS1624",
      ISS>=25 & ISS<=40 ~ "ISS2540",
      ISS>=41 & ISS<=49 ~ "ISS4149",
      ISS>=50 & ISS<=75 ~ "ISS5075",
      TRUE ~ "ISSxxxx"
      ))
d7<-mutate(d7,risscat=case_when(
      riss>=0 & riss<=8   ~ "riss0108",
      riss>=9 & riss<=15  ~ "riss0915",
      riss>=16 & riss<=24 ~ "riss1624",
      riss>=25 & riss<=40 ~ "riss2540",
      riss>=41 & riss<=49 ~ "riss4149",
      riss>=50 & riss<=75 ~ "riss5075",
      TRUE ~ "rissxxxx"
      ))
tabyl(d7,ISScat,risscat)
#A couple of very minor errors in Table 2


#ROC AUC
roc(d7$diedorhospice,d7$ISS)
roc(d7$diedorhospice,d7$riss)
roc(d7$diedorhospice,d7$PmortNIS)
roc(d7$diedorhospice,d7$PmortTQP)


#MORTALITY BY MAXIMUM AIS
#Note that mortality for AISmax=2 is less than for AISmax=1
t<-tabyl(d7,AISmax,diedorhospice)
adorn_percentages(t,"row")
t<-tabyl(d7,mxaisbr_H,diedorhospice)
adorn_percentages(t,"row")


#Can't quite replicate Dr. Piatt's models, except the first one below

#However, confirm adjusted odds ratio <1 for AIS=2 in "Ground Truth" model (Table 4)
m3a<-glm(diedorhospice ~ riss + riss2 + riss3
         +icdhead0 + icdhead2 + icdhead3 + icdhead4 + icdhead5 + age0005 + age0612 + age1317 + age6199, 
         data=d7, family="binomial")
exp(coef(m3a))
m3b<-glm(diedorhospice ~ ISS + ISS2 + ISS3 
         + head2 + head3 + head4 + head5 + head6 + age0005 + age0612 + age1317 + age6199,
         data=d7, family="binomial")
exp(coef(m3b))

#Also, confirm adjusted odds ratio <1 for ASIAN in "Ground Truth" model (Table 5)
m4a<-glm(diedorhospice ~ riss + riss2 + riss3
         + icdhead0 + icdhead2 + icdhead3 + icdhead4 + icdhead5 + age0005 + age0612 + age1317 + age6199
         + ASIAN + BLACK + RACEOTHER + RACE_UK + HISPANIC + ETHNICITY_UK 
         + commercial + other0005 + other0612 + other1317 + other6199, 
         data=d7, family="binomial")
exp(coef(m4a))
m4b<-glm(diedorhospice ~ ISS + ISS2 + ISS3
         + head2 + head3 + head4 + head5 + head6 + age0005 + age0612 + age1317 + age6199
         + ASIAN + BLACK + RACEOTHER + RACE_UK + HISPANIC + ETHNICITY_UK 
         + commercial + other0005 + other0612 + other1317 + other6199, 
         data=d7, family="binomial")
exp(coef(m4b))

#And confirm discrepancy for ASIAN category in simple LOS model (Table 6)
m5a<-lm(LOS ~ riss + ASIAN + BLACK + RACEOTHER + RACE_UK, data=d7)
coef(m5a)
m5b<-lm(LOS ~ ISS + ASIAN + BLACK + RACEOTHER + RACE_UK, data=d7)
coef(m5b)

#However, the effect of ASIAN on LOS is more similar after including other confounders
m6a<-lm(LOS ~ riss + +icdhead0 + icdhead2 + icdhead3 + icdhead4 + icdhead5 + age0005 + age0612 + age1317 + age6199
        + ASIAN + BLACK + RACEOTHER + RACE_UK + HISPANIC + ETHNICITY_UK 
        + commercial + other0005 + other0612 + other1317 + other6199, 
        data=d7)
coef(m6a)
m6b<-lm(LOS ~ ISS + head2 + head3 + head4 + head5 + head6 + age0005 + age0612 + age1317 + age6199
        + ASIAN + BLACK + RACEOTHER + RACE_UK + HISPANIC + ETHNICITY_UK 
        + commercial + other0005 + other0612 + other1317 + other6199, 
        data=d7)
coef(m6b)

#EXPLORE THE CONSEQUENCES OF DIFFERENT CATEGORIZATIONS OF MAXIMUM AIS
#Note that cases with AISmax=9 are included in reference category for models

d7<-mutate(d7,hcat=case_when(
  AISmax==1 & mxaisbr_H==0 ~ "h101",
  AISmax==1 & mxaisbr_H==1 ~ "h101",
  AISmax==1 & mxaisbr_H==2 ~ "h12",
  AISmax==1 & mxaisbr_H==3 ~ "h13",
  AISmax==1 & mxaisbr_H==4 ~ "h14",
  AISmax==1 & mxaisbr_H==5 ~ "h15",
  AISmax==2 & mxaisbr_H==0 ~ "h201",
  AISmax==2 & mxaisbr_H==1 ~ "h201",
  AISmax==2 & mxaisbr_H==2 ~ "h22",
  AISmax==2 & mxaisbr_H==3 ~ "h23",
  AISmax==2 & mxaisbr_H==4 ~ "h24",
  AISmax==2 & mxaisbr_H==5 ~ "h25",
  AISmax==3 & mxaisbr_H==0 ~ "h301",
  AISmax==3 & mxaisbr_H==1 ~ "h301",
  AISmax==3 & mxaisbr_H==2 ~ "h32",
  AISmax==3 & mxaisbr_H==3 ~ "h33",
  AISmax==3 & mxaisbr_H==4 ~ "h34",
  AISmax==3 & mxaisbr_H==5 ~ "h35",
  AISmax==4 & mxaisbr_H==0 ~ "h401",
  AISmax==4 & mxaisbr_H==1 ~ "h401",
  AISmax==4 & mxaisbr_H==2 ~ "h42",
  AISmax==4 & mxaisbr_H==3 ~ "h43",
  AISmax==4 & mxaisbr_H==4 ~ "h44",
  AISmax==4 & mxaisbr_H==5 ~ "h45",
  AISmax==5 & mxaisbr_H==0 ~ "h5601",
  AISmax==5 & mxaisbr_H==1 ~ "h5601",
  AISmax==5 & mxaisbr_H==2 ~ "h562",
  AISmax==5 & mxaisbr_H==3 ~ "h563",
  AISmax==5 & mxaisbr_H==4 ~ "h564",
  AISmax==5 & mxaisbr_H==5 ~ "h565",
  AISmax==6 & mxaisbr_H==0 ~ "h5601",
  AISmax==6 & mxaisbr_H==1 ~ "h5601",
  AISmax==6 & mxaisbr_H==2 ~ "h562",
  AISmax==6 & mxaisbr_H==3 ~ "h563",
  AISmax==6 & mxaisbr_H==4 ~ "h564",
  AISmax==6 & mxaisbr_H==5 ~ "h565",
  AISmax==9 & mxaisbr_H==0 ~ "h901",
  AISmax==9 & mxaisbr_H==1 ~ "h901",
  AISmax==9 & mxaisbr_H==2 ~ "h92",
  AISmax==9 & mxaisbr_H==3 ~ "h93",
  AISmax==9 & mxaisbr_H==4 ~ "h94",
  AISmax==9 & mxaisbr_H==5 ~ "h95",
  TRUE ~ "hx"
))

td<-tabyl(d7,hcat,died)
td
adorn_percentages(td,"row")

ta<-tabyl(d7,hcat,ASIAN)
ta
adorn_percentages(ta,"row")

#Convert above table to likelihood ratios
tabyl(d7,ASIAN)
A=10182
N=337148
tlr<-rename(ta,a="1",nota="0")
tlr<-mutate(tlr,n=a+nota)
tlr<-mutate(tlr,lr=(a/n)/( (A-a)/(N-n) ))
tlr

