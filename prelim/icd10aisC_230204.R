

#' Categorize trauma data
#'
#' This function adds Abbreviated Injury Scores (AIS), Injury Severity Scores (ISS), and other descriptors of injury to a dataframe.
#' For each observation this function will
#' \enumerate{
#'    \item assign a severity (AIS) and ISS body region values to each valid ICD-9 or ICD-10 injury diagnosis code,
#'    \item add variables for maximum severity of each body region,
#'    \item calculate ISS, "New ISS", maximum AIS, and a regression-based mortality prediction,
#'    \item select first 4 e-codes/mechanism codes and categorize major mechanism, minor mechanism, and intent
#'}
#'
#'


##########################################################################
#
#      R PROGRAMS TO EXTRACT INJURY DATA
#           FROM TQP PUF (FORMERLY TQIP/NTDB) AND FROM NIS
#           AND ANALYZE ROCMAX OPTIONS FOR ICDPIC-R
#
#      PART C:  USE REGRESSION EFFECTS FROM PART B TO RANK INJURY SEVERITY
#               DETERMINE AIS CUTPOINTS THAT OPTIMIZE ISS
#
#      David Clark, 2022-2023
#
##########################################################################


#1
#CLEAR WORKSPACE, SET WORKING DIRECTORY,
#LOAD REQUIRED PACKAGES, IF NOT LOADED ALREADY

rm(list=ls())
setwd("/Users/davideugeneclark/Documents/icdpicr")
require(tidyverse)
require(skimr)
require(janitor)
require(broom)
require(pROC)


#2
# IMPORT DATA WITH REGRESSION EFFECT OF EACH VALID INJURY DIAGNOSIS

 d0<-read_csv("tqip2020cm_effects.csv")
 d0<-read_csv("nis2020cm_effects.csv")

d1<-rename(d0,icdcm=dx)


#3
# OBTAIN MAXIMUM "EFFECT" FOR EACH PERSON AND BODY REGION

d3<-rename(d1,effect=x)
d3<-group_by(d3,INC_KEY,issbr)
d3<-mutate(d3,a0=ifelse(issbr=="A",max(effect),-99))
d3<-mutate(d3,c0=ifelse(issbr=="C",max(effect),-99))
d3<-mutate(d3,e0=ifelse(issbr=="E",max(effect),-99))
d3<-mutate(d3,f0=ifelse(issbr=="F",max(effect),-99))
d3<-mutate(d3,h0=ifelse(issbr=="H",max(effect),-99))
d3<-mutate(d3,s0=ifelse(issbr=="G",max(effect),-99))
d3<-ungroup(d3)
d3<-ungroup(d3)
d3<-group_by(d3,INC_KEY)
d3<-mutate(d3,idseq=row_number())
d3<-mutate(d3,amax=max(a0))
d3<-mutate(d3,cmax=max(c0))
d3<-mutate(d3,emax=max(e0))
d3<-mutate(d3,fmax=max(f0))
d3<-mutate(d3,hmax=max(h0))
d3<-mutate(d3,smax=max(s0))

d3<-ungroup(d3)

d4<-filter(d3,idseq==1)


########################################################################################

#4
#IMPLEMENT "GREEDY" ALGORITHM TO SEEK OPTIMAL CUTPOINTS FOR SUM OF SQUARES


#Define a function to return a triangular random variate
#   with range a to b and mode m.
trirandom<-function(a,b,m) {
  k=(m-a)/(b-a)
  ru01=runif(1)
  rtri01k=if_else(ru01<=k,sqrt(k*ru01),(1-sqrt((1-k)*(1-ru01))))
  rtriabm=a+(b-a)*rtri01k
  return(rtriabm)
}


#sink("x221114.txt",split=T)
starttime=Sys.time()

cut12=0.1
cut23=0.5
cut34=0.9
cut45=1.6
bestroc=.5
i=0

while(i<1000) {

  i=i+1
  old12=cut12
  old23=cut23
  old34=cut34
  old45=cut45
  if(i%%4==1) {
    cut12=trirandom(-1,cut23,cut12)
  }
  if(i%%4==2) {
    cut23=trirandom(cut12,cut34,cut23)
  }
  if(i%%4==3) {
    cut34=trirandom(cut23,cut45,cut34)
  }
  if(i%%4==0) {
    cut45=trirandom(cut34,2,cut45)
  }

  d4<-mutate(d4,a=case_when(
    amax>=-9 & amax<cut12 ~ 1,
    amax>=cut12 & amax<cut23 ~ 2,
    amax>=cut23 & amax<cut34 ~ 3,
    amax>=cut34 & amax<cut45 ~ 4,
    amax>=cut45 ~ 5,
    TRUE ~ 0
  ) )
  d4<-mutate(d4,c=case_when(
    cmax>=-9 & cmax<cut12 ~ 1,
    cmax>=cut12 & cmax<cut23 ~ 2,
    cmax>=cut23 & cmax<cut34 ~ 3,
    cmax>=cut34 & cmax<cut45 ~ 4,
    cmax>=cut45 ~ 5,
    TRUE ~ 0
  ) )
  d4<-mutate(d4,e=case_when(
    emax>=-9 & emax<cut12 ~ 1,
    emax>=cut12 & emax<cut23 ~ 2,
    emax>=cut23 & emax<cut34 ~ 3,
    emax>=cut34 & emax<cut45 ~ 4,
    emax>=cut45 ~ 5,
    TRUE ~ 0
  ) )
  d4<-mutate(d4,f=case_when(
    fmax>=-9 & fmax<cut12 ~ 1,
    fmax>=cut12 & fmax<cut23 ~ 2,
    fmax>=cut23 & fmax<cut34 ~ 3,
    fmax>=cut34 & fmax<cut45 ~ 4,
    fmax>=cut45 ~ 5,
    TRUE ~ 0
  ) )
  d4<-mutate(d4,h=case_when(
    hmax>=-9 & hmax<cut12 ~ 1,
    hmax>=cut12 & hmax<cut23 ~ 2,
    hmax>=cut23 & hmax<cut34 ~ 3,
    hmax>=cut34 & hmax<cut45 ~ 4,
    hmax>=cut45 ~ 5,
    TRUE ~ 0
  ) )
  d4<-mutate(d4,s=case_when(
    smax>=-9 & smax<cut12 ~ 1,
    smax>=cut12 & smax<cut23 ~ 2,
    smax>=cut23 & smax<cut34 ~ 3,
    smax>=cut34 & smax<cut45 ~ 4,
    smax>=cut45 ~ 5,
    TRUE ~ 0
  ) )
  d4<-mutate(d4,a2=a^2)
  d4<-mutate(d4,c2=c^2)
  d4<-mutate(d4,e2=e^2)
  d4<-mutate(d4,f2=f^2)
  d4<-mutate(d4,h2=h^2)
  d4<-mutate(d4,s2=s^2)
  d4<-mutate(d4,sumsquares=a2+c2+e2+f2+h2+s2)

  testroc<-roc(d4$died,d4$sumsquares,quiet=T)
  newroc=as.numeric(auc(testroc))

  if(newroc<=bestroc) {
    cut12=old12
    cut23=old23
    cut34=old34
    cut45=old45
  }
  if(newroc>bestroc){
    print(c("i:",i,"Cuts:",format(cut12,digits=3),format(cut23,digits=3),format(cut34,digits=3),
            format(cut45,digits=3),"AUC:",format(newroc,digits=4)))
    bestroc=newroc
  }
  if(i%%100==0){
    print(c("i:",i))
  }

} #while


endtime=Sys.time()
endtime-starttime
sink(NULL)

#Best for TQIPcm: .146, .349, .814, 1.67:  R-squared for sum squares = .8408, for ISS = .8425
#Best for NIScm: .0995, .507, .965, 1.63;  R-squared for sum squares = .8065, for ISS = .8074

#################################################################################################


#5
#VERIFY SUM OF SQUARES AND CALCULATE ISS

cut12=.146
cut23=.349
cut34=.814
cut45=1.67

cut12=.0995
cut23=.507
cut34=.965
cut45=1.63

d5<-mutate(d4,a=case_when(
  amax>-9 & amax<cut12 ~ 1,
  amax>=cut12 & amax<cut23 ~ 2,
  amax>=cut23 & amax<cut34 ~ 3,
  amax>=cut34 & amax<cut45 ~ 4,
  amax>=cut45 ~ 5,
  TRUE ~ 0
) )
d5<-mutate(d5,c=case_when(
  cmax>-9 & cmax<cut12 ~ 1,
  cmax>=cut12 & cmax<cut23 ~ 2,
  cmax>=cut23 & cmax<cut34 ~ 3,
  cmax>=cut34 & cmax<cut45 ~ 4,
  cmax>=cut45 ~ 5,
  TRUE ~ 0
) )
d5<-mutate(d5,e=case_when(
  emax>-9 & emax<cut12 ~ 1,
  emax>=cut12 & emax<cut23 ~ 2,
  emax>=cut23 & emax<cut34 ~ 3,
  emax>=cut34 & emax<cut45 ~ 4,
  emax>=cut45 ~ 5,
  TRUE ~ 0
) )
d5<-mutate(d5,f=case_when(
  fmax>-9 & fmax<cut12 ~ 1,
  fmax>=cut12 & fmax<cut23 ~ 2,
  fmax>=cut23 & fmax<cut34 ~ 3,
  fmax>=cut34 & fmax<cut45 ~ 4,
  fmax>=cut45 ~ 5,
  TRUE ~ 0
) )
d5<-mutate(d5,h=case_when(
  hmax>-9 & hmax<cut12 ~ 1,
  hmax>=cut12 & hmax<cut23 ~ 2,
  hmax>=cut23 & hmax<cut34 ~ 3,
  hmax>=cut34 & hmax<cut45 ~ 4,
  hmax>=cut45 ~ 5,
  TRUE ~ 0
) )
d5<-mutate(d5,s=case_when(
  smax>-9 & smax<cut12 ~ 1,
  smax>=cut12 & smax<cut23 ~ 2,
  smax>=cut23 & smax<cut34 ~ 3,
  smax>=cut34 & smax<cut45 ~ 4,
  smax>=cut45 ~ 5,
  TRUE ~ 0
) )
d5<-mutate(d5,a2=a^2)
d5<-mutate(d5,c2=c^2)
d5<-mutate(d5,e2=e^2)
d5<-mutate(d5,f2=f^2)
d5<-mutate(d5,h2=h^2)
d5<-mutate(d5,s2=s^2)
d5<-mutate(d5,sumsquares=a2+c2+e2+f2+h2+s2)
skim(d5,a,c,e,f,h,s,sumsquares)
roc(d5$died,d5$sumsquares)
tabyl(d5,sumsquares) %>%
  adorn_pct_formatting(digits=2)


#Calculate actual ISS
d6<-gather(d5,a2,c2,e2,f2,h2,s2,key="reg",value="regmax")
d6<-group_by(d6,INC_KEY)
d6<-arrange(d6,desc(regmax))
d6<-mutate(d6,regorder=row_number())
d6<-mutate(d6,xiss=ifelse(regorder<=3,cumsum(regmax),0))
d6<-mutate(d6,riss=max(xiss))
d6<-ungroup(d6)
d6<-filter(d6,regorder==1)
d6<-select(d6,-reg,-regmax,-regorder,-xiss)
roc(d6$died,d6$riss)
roc(d6$died,d6$totridge)

#Compare results to ISS in registry data
d7<-mutate(d6,risscat=case_when(
  riss>0 & riss<=8 ~ "01-08",
  riss>8 & riss<=15 ~ "09-15",
  riss>15 & riss<=24 ~ "16-24",
  riss>24 & riss<=40 ~ "25-40",
  riss>40 & riss<=49 ~ "41-49",
  riss>49 & riss<=75 ~ "50-75",
  TRUE ~ "Unk"
) )
d7<-mutate(d7,zisscat=case_when(
  TQPISS>0 & TQPISS<=8 ~ "01-08",
  TQPISS>8 & TQPISS<=15 ~ "09-15",
  TQPISS>15 & TQPISS<=24 ~ "16-24",
  TQPISS>24 & TQPISS<=40 ~ "25-40",
  TQPISS>40 & TQPISS<=49 ~ "41-49",
  TQPISS>49 & TQPISS<=75 ~ "50-75",
  TRUE ~ "Unk"
) )

d7<-mutate(d7,TQPISS=ifelse(TQPISS==-2,NaN,TQPISS))

roc(d7$died,d7$TQPISS)
roc(d7$died,d7$riss)
t<-tabyl(d7,zisscat,died)
t2<-adorn_percentages(t,"row")
t3<-adorn_pct_formatting(t2, digits=2)
t4<-adorn_ns(t3,"front")
t4
t<-tabyl(d7,risscat,died)
t2<-adorn_percentages(t,"row")
t3<-adorn_pct_formatting(t2, digits=2)
t4<-adorn_ns(t3,"front")
t4
t5<-tabyl(d7,zisscat,risscat)
t5
t6<-adorn_percentages(t5,"all")
t7<-adorn_pct_formatting(t6, digits=2)
t8<-adorn_ns(t7,"front")
t8

##############################################################################################



#6
#CREATE TABLE OF INCLUDED DIAGNOSIS CODES WITH EFFECTS AND AIS SCORES, SAVE RESULT
#  Use optimal AIS cutpoints for each source, as determined above

rm(list=ls())
 d1<-read_csv("tqip2020cm_effects.csv")
  d1<-read_csv("nis2020cm_effects.csv")


#Identify and remove duplicates, drop unnecessary variables
d2<-rename(d1,effect=x)
d2<-group_by(d2,dx)
d2<-mutate(d2,dxseq=row_number())
d2<-mutate(d2,dxrep=max(dxseq))
d2<-ungroup(d2)
d3<-filter(d2,dxseq==1)
d3<-arrange(d3,dx)
d3<-select(d3,dx,issbr,effect,ridge_int,dxrep)

####################################################
#  ENTER APPROPRIATE CUTPOINTS AND ASSIGN AIS VALUES
####################################################

#Best for TQIPcm: .146, .349, .814, 1.67:  R-squared for sum squares = .8408, for ISS = .8567
#Best for NIScm: .0995, .507, .965, 1.63;  R-squared for sum squares = .8065, for ISS = .8074

#################################################################################################

cut12=.146
cut23=.349
cut34=.814
cut45=1.67

cut12=.0995
cut23=.507
cut34=.965
cut45=1.63


#Assign AIS values
d4<-mutate(d3,ais=case_when(
  effect>-99 & effect<cut12 ~ 1,
  effect>=cut12 & effect<cut23 ~ 2,
  effect>=cut23 & effect<cut34 ~ 3,
  effect>=cut34 & effect<cut45 ~ 4,
  effect>=cut45 ~ 5,
  TRUE ~ 0
) )


#########################################
#  ONLY RUN THE APPROPRIATE SECTION BELOW
#########################################

#Rename diagnosis code type as appropriate and save results

#####FOR TQIP CM
d5<-rename(d4,icdcm=dx)
d5<-rename(d5,TQIPeffect=effect)
d5<-rename(d5,TQIPint=ridge_int)
d5<-rename(d5,TQIPais=ais)
d5<-rename(d5,TQIPn=dxrep)
d5<-rename(d5,TQIPbr=issbr)
write_csv(d5,"tqip2020cm_ais.csv")

#####FOR NIS CM
d5<-rename(d4,icdcm=dx)
d5<-rename(d5,NISeffect=effect)
d5<-rename(d5,NISint=ridge_int)
d5<-rename(d5,NISais=ais)
d5<-rename(d5,NISn=dxrep)
d5<-rename(d5,NISbr=issbr)
write_csv(d5,"nis2020cm_ais.csv")





