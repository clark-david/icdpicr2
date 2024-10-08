####################################################################
#
#      R PROGRAMS TO EXTRACT INJURY DATA 
#           FROM TQP PUF (FORMERLY TQIP/NTDB) AND FROM NIS
#           AND ANALYZE ROCMAX OPTIONS FOR ICDPIC-R
#
#      PART B: USING RIDGE REGRESSION, 
#                     ESTIMATE INDEPENDENT EFFECT OF EACH DIAGNOSIS
#              ADD ISS BODY REGIONS FOR EACH DIAGNOSIS
#
#      David Clark, 2022-2023
#                     
#####################################################################


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
require(lme4)
require(ggplot2)
require(glmnet)


#2
#OBTAIN DATA FROM PART A

d0<-read_csv("tqip2020cm.csv")
d0<-read_csv("nis2020cm.csv")

d1<-rename(d0,dx=icdcm)

#3
#CONSTRUCT MATRICES FOR RIDGE REGRESSION
# (in pieces due to memory restrictions)


####################################################################
#REPEAT INCREMENTING modno FROM 0 TO 9

modno=9
#Create dummy observation with all diagnoses
# (so all pieces have the same columns)
ddummy<-group_by(d1,dx)
ddummy<-mutate(ddummy,dxseq=row_number())
ddummy<-ungroup(ddummy)
ddummy<-filter(ddummy,dxseq==1)
ddummy<-select(ddummy,INC_KEY,dx,died)
ddummy<-mutate(ddummy,died=0)
ddummy<-mutate(ddummy,INC_KEY=modno)
#Take 10% sample with mod(INC_KEY,10)=modno and append to dummy
dmod<-filter(d1,INC_KEY%%10==modno)
dmod<-bind_rows(ddummy,dmod)
dmod<-arrange(dmod,INC_KEY,dx)

#Convert mortality data to a vector with one observation per person
dmort<-group_by(dmod,INC_KEY)
dmort<-mutate(dmort,idseq=row_number())
dmort<-ungroup(dmort)
dmort<-filter(dmort,idseq==1)
dmort<-select(dmort,died)
matmort<-data.matrix(dmort)

#Convert diagnosis data to wide format and sparse logic (T/F) matrix
#  (This is the part that needs some memory and time)
ddx<-mutate(dmod,x=T)
ddx<-select(ddx,INC_KEY,dx,x)
ddx<-spread(ddx,dx,x,fill=F)
ddx<-select(ddx,-INC_KEY)

smatdx<-Matrix(as.matrix(ddx),sparse=TRUE)
smatdx[1:5,1:5]
object.size(ddx)    #  3.1GB for TQIPcm, 0.6GB for NIScm
object.size(smatdx) #  3.7MB for TQIPcm, 1.0MB for NIScm
rm(ddx,dmort)

#Save matrices named by source and modno
#write(matmort,"tqip2020cm_matmort9",ncolumns=1)
#writeMM(smatdx,"tqip2020cm_smatdx9")
write(matmort,"nis2020cm_matmort9",ncolumns=1)
writeMM(smatdx,"nis2020cm_smatdx9")

##############################################################


#Import partial data matrices
#FOR TQIP:
matmort0<-as.matrix(read.table("tqip2020cm_matmort0"))
smatdx0<-readMM("tqip2020cm_smatdx0")
matmort1<-as.matrix(read.table("tqip2020cm_matmort1"))
smatdx1<-readMM("tqip2020cm_smatdx1")
matmort2<-as.matrix(read.table("tqip2020cm_matmort2"))
smatdx2<-readMM("tqip2020cm_smatdx2")
matmort3<-as.matrix(read.table("tqip2020cm_matmort3"))
smatdx3<-readMM("tqip2020cm_smatdx3")
matmort4<-as.matrix(read.table("tqip2020cm_matmort4"))
smatdx4<-readMM("tqip2020cm_smatdx4")
matmort5<-as.matrix(read.table("tqip2020cm_matmort5"))
smatdx5<-readMM("tqip2020cm_smatdx5")
matmort6<-as.matrix(read.table("tqip2020cm_matmort6"))
smatdx6<-readMM("tqip2020cm_smatdx6")
matmort7<-as.matrix(read.table("tqip2020cm_matmort7"))
smatdx7<-readMM("tqip2020cm_smatdx7")
matmort8<-as.matrix(read.table("tqip2020cm_matmort8"))
smatdx8<-readMM("tqip2020cm_smatdx8")
matmort9<-as.matrix(read.table("tqip2020cm_matmort9"))
smatdx9<-readMM("tqip2020cm_smatdx9")
#FOR NIS
matmort0<-as.matrix(read.table("nis2020cm_matmort0"))
smatdx0<-readMM("nis2020cm_smatdx0")
matmort1<-as.matrix(read.table("nis2020cm_matmort1"))
smatdx1<-readMM("nis2020cm_smatdx1")
matmort2<-as.matrix(read.table("nis2020cm_matmort2"))
smatdx2<-readMM("nis2020cm_smatdx2")
matmort3<-as.matrix(read.table("nis2020cm_matmort3"))
smatdx3<-readMM("nis2020cm_smatdx3")
matmort4<-as.matrix(read.table("nis2020cm_matmort4"))
smatdx4<-readMM("nis2020cm_smatdx4")
matmort5<-as.matrix(read.table("nis2020cm_matmort5"))
smatdx5<-readMM("nis2020cm_smatdx5")
matmort6<-as.matrix(read.table("nis2020cm_matmort6"))
smatdx6<-readMM("nis2020cm_smatdx6")
matmort7<-as.matrix(read.table("nis2020cm_matmort7"))
smatdx7<-readMM("nis2020cm_smatdx7")
matmort8<-as.matrix(read.table("nis2020cm_matmort8"))
smatdx8<-readMM("nis2020cm_smatdx8")
matmort9<-as.matrix(read.table("nis2020cm_matmort9"))
smatdx9<-readMM("nis2020cm_smatdx9")

#Remove initial dummy row from each matrix
matmort0<-as.matrix(matmort0[-1,])
smatdx0<-smatdx0[-1,]
matmort1<-as.matrix(matmort1[-1,])
smatdx1<-smatdx1[-1,]
matmort2<-as.matrix(matmort2[-1,])
smatdx2<-smatdx2[-1,]
matmort3<-as.matrix(matmort3[-1,])
smatdx3<-smatdx3[-1,]
matmort4<-as.matrix(matmort4[-1,])
smatdx4<-smatdx4[-1,]
matmort5<-as.matrix(matmort5[-1,])
smatdx5<-smatdx5[-1,]
matmort6<-as.matrix(matmort6[-1,])
smatdx6<-smatdx6[-1,]
matmort7<-as.matrix(matmort7[-1,])
smatdx7<-smatdx7[-1,]
matmort8<-as.matrix(matmort8[-1,])
smatdx8<-smatdx8[-1,]
matmort9<-as.matrix(matmort9[-1,])
smatdx9<-smatdx9[-1,]

#Combine matrices and save 
matmort<-rbind(matmort0,matmort1,matmort2,matmort3,matmort4,matmort5,matmort6,matmort7,matmort8,matmort9)
smatdx<-rbind(smatdx0,smatdx1,smatdx2,smatdx3,smatdx4,smatdx5,smatdx6,smatdx7,smatdx8,smatdx9)

#FOR TQIP
write(matmort,"tqip2020cm_matmort",ncolumns=1)
writeMM(smatdx,"tqip2020cm_smatdx")
#FOR NIS
write(matmort,"nis2020cm_matmort",ncolumns=1)
writeMM(smatdx,"nis2020cm_smatdx")
####################################################################################


#3
#RIDGE REGRESSION

#Obtain mortality and diagnosis data in matrix format
matmort<-as.matrix(read.table("tqip2020cm_matmort"))
smatdx<-readMM("tqip2020cm_smatdx")
matmort<-as.matrix(read.table("nis2020cm_matmort"))
smatdx<-readMM("nis2020cm_smatdx")

#Convert ngTMatrix (logical) to dgCMatrix (numeric)
#Otherwise cv.glmnet doesn't work
smatdx=smatdx*1

#Determine optimal value of lambda by cross-validation
#Note: Ridge regression if alpha=0 (LASSO if alpha=1)
cvridge<-cv.glmnet(smatdx,matmort,family="binomial",alpha=0) 
#plot(cvridge)
cvridge$lambda.min
log(cvridge$lambda.min)
#Get coefficients for model with lowest binomial deviance
mridge<-coef(cvridge,s="lambda.min")
ridge_eff<-as.data.frame(summary(mridge))
skim(ridge_eff)
intercept=ridge_eff[1,3]


#Obtain corresponding data in original format
d0<-read_csv("tqip2020cm.csv")
d0<-read_csv("nis2020cm.csv")
d0<-rename(d0,dx=icdcm)
d0<-arrange(d0,INC_KEY,dx)


#Extract corresponding list of all diagnoses in the data
d1<-group_by(d0,dx)
d1<-mutate(d1,dxseq=row_number())
d1<-ungroup(d1)
d1<-filter(d1,dxseq==1)
d1<-select(d1,dx)
d1<-arrange(d1,dx)
d1<-mutate(d1,dxseq=row_number())


#Line up rows correctly to merge data
head(mridge)
head(ridge_eff)
head(d1)
d2<-mutate(d1,i=dxseq+1)
head(d2)
ridge_eff2<-slice(ridge_eff,-1)
head(ridge_eff2)
d3<-full_join(d2,ridge_eff2,"i")
head(d3)

d4<-full_join(d0,d3,"dx")
d4<-group_by(d4,INC_KEY)
d4<-mutate(d4,idseq=row_number())
d4<-mutate(d4,x=if_else(is.na(x),0,x))
d4<-mutate(d4,totridge=sum(x))
d4<-ungroup(d4)
d4<-mutate(d4,ridge_int=intercept)
d4<-arrange(d4,INC_KEY,dx)

#Explore results
d4test=filter(d4,idseq==1)
lm(died~totridge,data=d4test)
roc(d4test$died,d4test$totridge,quiet=T)


#4 
# ADD BODY REGION FOR EACH DIAGNOSIS

d5<-mutate(d4,issbr=case_when(
  str_sub(dx,1,3)=="S00" & str_sub(dx,4,4)=="1" ~ "F",
  str_sub(dx,1,3)=="S00" & str_sub(dx,4,4)=="2" ~ "F",
  str_sub(dx,1,3)=="S00" & str_sub(dx,4,4)=="3" ~ "F",
  str_sub(dx,1,3)=="S00" & str_sub(dx,4,4)=="4" ~ "F",
  str_sub(dx,1,3)=="S00" & str_sub(dx,4,4)=="5" ~ "F",
  str_sub(dx,1,3)=="S01" & str_sub(dx,4,4)=="1" ~ "F",
  str_sub(dx,1,3)=="S01" & str_sub(dx,4,4)=="2" ~ "F",
  str_sub(dx,1,3)=="S01" & str_sub(dx,4,4)=="3" ~ "F",
  str_sub(dx,1,3)=="S01" & str_sub(dx,4,4)=="4" ~ "F",
  str_sub(dx,1,3)=="S01" & str_sub(dx,4,4)=="5" ~ "F",
  str_sub(dx,1,3)=="S02" & str_sub(dx,4,4)=="1" ~ "F",
  str_sub(dx,1,3)=="S02" & str_sub(dx,4,4)=="2" ~ "F",
  str_sub(dx,1,3)=="S02" & str_sub(dx,4,4)=="3" ~ "F",
  str_sub(dx,1,3)=="S02" & str_sub(dx,4,4)=="4" ~ "F",
  str_sub(dx,1,3)=="S02" & str_sub(dx,4,4)=="5" ~ "F",
  str_sub(dx,1,3)=="S02" & str_sub(dx,4,4)=="6" ~ "F",
  str_sub(dx,1,3)=="S03" & str_sub(dx,4,4)=="0" ~ "F",
  str_sub(dx,1,3)=="S03" & str_sub(dx,4,4)=="1" ~ "F",
  str_sub(dx,1,3)=="S03" & str_sub(dx,4,4)=="2" ~ "F",
  str_sub(dx,1,3)=="S03" & str_sub(dx,4,4)=="3" ~ "F",
  str_sub(dx,1,3)=="S03" & str_sub(dx,4,4)=="4" ~ "F",
  str_sub(dx,1,3)=="S03" & str_sub(dx,4,4)=="5" ~ "F",
  str_sub(dx,1,3)=="S05"   ~ "F",
  str_sub(dx,1,4)=="S070" ~ "F",
  str_sub(dx,1,4)=="S081" ~ "F",
  str_sub(dx,1,4)=="S092" ~ "F",
  str_sub(dx,1,2)=="S0" ~ "H",
  str_sub(dx,1,2)=="S1" ~ "H",
  str_sub(dx,1,2)=="S2" ~ "C",
  str_sub(dx,1,2)=="S3" ~ "A",
  str_sub(dx,1,2)=="S4" ~ "E",
  str_sub(dx,1,2)=="S5" ~ "E",
  str_sub(dx,1,2)=="S6" ~ "E",
  str_sub(dx,1,2)=="S7" ~ "E",
  str_sub(dx,1,2)=="S8" ~ "E",
  str_sub(dx,1,2)=="S9" ~ "E",
  str_sub(dx,1,4)=="T242" ~ "E",
  str_sub(dx,1,4)=="T79A" ~ "E",
  str_sub(dx,1,1)=="T" ~ "G",
  TRUE ~ "N"  
) )
tabyl(d5,issbr)

#Discard temporary analytic variables and save result
d6<-select(d5,-idseq,-dxseq,-i,-j)
write_csv(d6,"tqip2020cm_effects.csv")
write_csv(d6,"nis2020cm_effects.csv")





