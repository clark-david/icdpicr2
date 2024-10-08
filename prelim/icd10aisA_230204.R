###################################################################################
#
#      R PROGRAMS TO EXTRACT INJURY DATA 
#           FROM TQP PUF (FORMERLY TQIP/NTDB) AND FROM NIS
#           AND ANALYZE ROCMAX OPTIONS FOR ICDPIC-R
#
#      PART A: EXTRACT DATA
#              RESTRICT TO CASES WITH SPECIFIED ICD-10-CM INJURY CODES
#               
#      David Clark, 2022-2023
#                     
####################################################################################


#1
#CLEAR WORKSPACE, SET WORKING DIRECTORY,
#LOAD REQUIRED PACKAGES, IF NOT LOADED ALREADY
rm(list=ls())
setwd("/Users/davideugeneclark/Documents/icdpicr")  
require(tidyverse)
require(janitor)
require(broom)
require(skimr)


#####################################################################################


#2a
#GET TQIP DATA  

#Import raw 2020 TQIP data 
d0<-read_csv("PUF AY 2020/CSV/PUF_TRAUMA.csv")  
d1<-rename(d0,INC_KEY=inc_key,TQPISS=ISS)

#Demonstrate that INC_KEY is unique identifier
d1<-group_by(d1,INC_KEY)
d1<-mutate(d1,seq=row_number())
d1<-ungroup(d1)
tabyl(d1,seq)

#Save recorded outcome data and ISS
d2<-select(d1,INC_KEY,HOSPDISCHARGEDISPOSITION,TQPISS)
write_csv(d2,"tqipmort2020.csv")  

#Obtain ICD-10-CM codes for each TQIP patient
d0<-read_csv("PUF AY 2020/CSV/PUF_ICDDIAGNOSIS.csv")
d1<-rename(d0,INC_KEY=Inc_Key,withdot=ICDDIAGNOSISCODE)

#Identify duplicates
d2<-group_by(d1,INC_KEY,withdot)
d2<-mutate(d2,dxseq=row_number())
d2<-mutate(d2,dxrep=max(dxseq))
d2<-ungroup(d2)
d2test<-filter(d2,dxrep>1)
d2test<-arrange(d2test,INC_KEY,withdot,dxseq,dxrep)
#View(d2test)
#Drop duplicates
d3<-filter(d2,dxseq==1)
d3<-select(d3,INC_KEY,withdot)
d3<-arrange(d3,INC_KEY,withdot)

#Extract parts of each ICD-10-CM code
#Convert to format without dot (e.g., S060X0A not S06.0X0A)
#Restrict to codes with first digit "S" or "T"
#   and seventh digit "A","B","C", or blank  
d4<-arrange(d3,INC_KEY,withdot)
d4<-mutate(d4,predot=str_sub(withdot,1,3))
d4<-mutate(d4,postdot=str_sub(withdot,5,8))
d4<-mutate(d4,icdcm=str_c(predot,postdot))
d4<-mutate(d4,digit1=str_sub(icdcm,1,1))
d4<-mutate(d4,digit7=str_sub(icdcm,7,7))
d4<-filter(d4,digit1=="S" | digit1=="T")
d4<-filter(d4,digit7=="A" | digit7=="B" | digit7=="C" | digit7=="")

d5<-select(d4,INC_KEY,icdcm)
write_csv(d5,"tqipdcode2020.csv")

#Merge TQIP patient diagnoses and outcomes
d1<-read_csv("tqipdcode2020.csv")
d2<-read_csv("tqipmort2020.csv")  
d3<-full_join(d1,d2,by="INC_KEY")  

#Drop if outcome unspecified or if transferred to another acute care hospital
#Generate "died" outcome, and save result
d4<-filter(d3,HOSPDISCHARGEDISPOSITION>=2 & HOSPDISCHARGEDISPOSITION<=14)
d4<-mutate(d4,died=if_else(
  HOSPDISCHARGEDISPOSITION==5,1,0))
tabyl(d4,HOSPDISCHARGEDISPOSITION,died)
write_csv(d4,"tqipmerged2020.csv")


#2b  
#GET NIS DATA

#Import raw 2020 NIS cases 

d0<-read_fwf("NIS_2020/NIS_2020_Core.ASC",
             col_positions=fwf_positions(
               c(1,10,23,35,39,55,62,69,76,83,90,97,104,111,118,125,132,139,146,153,160,167,174,181,188,
                 195,202,209,216,223,230,237,244,251,258,265,272,279,286,293,300,307,314,321,328,335,521,531),
               c(3,11,24,36,41,61,68,75,82,89,96,103,110,117,124,131,138,145,152,159,166,173,180,187,194,
                 201,208,215,222,229,236,243,250,257,264,271,278,285,292,299,306,313,320,327,334,336,530,535),
               c("AGE","DIED","DISPUNIFORM","ELECTIVE","HCUP_ED",
                 "I10_DX1","I10_DX2","I10_DX3","I10_DX4","I10_DX5","I10_DX6","I10_DX7","I10_DX8",
                 "I10_DX9","I10_DX10","I10_DX11","I10_DX12","I10_DX13","I10_DX14","I10_DX15","I10_DX16",
                 "I10_DX17","I10_DX18","I10_DX19","I10_DX20","I10_DX21","I10_DX22","I10_DX23","I10_DX24",
                 "I10_DX25","I10_DX26","I10_DX27","I10_DX28","I10_DX29","I10_DX30","I10_DX31","I10_DX32",
                 "I10_DX33","I10_DX34","I10_DX35","I10_DX36","I10_DX37","I10_DX38","I10_DX39","I10_DX40",
                 "INJURY","INC_KEY","LOS")
             ))
View(d0)

#Include only admissions that were not elective and originated in Emergency Department
tabyl(d0,ELECTIVE,HCUP_ED)
d1<-filter(d0,ELECTIVE==0,HCUP_ED!=0)
tabyl(d1,ELECTIVE,HCUP_ED)
#Include only admissions with a principal AHRQ "injury" diagnosis (S00-T34 and some others)
tabyl(d1,INJURY)
d2<-filter(d1,INJURY==1)
tabyl(d2,INJURY)

#Verify that INC_KEY is a unique identifier
d2<-group_by(d2,INC_KEY)
d2<-mutate(d2,seq=row_number())
d2<-ungroup(d2)
tabyl(d2,seq)

write_csv(d2,"NISraw2020.csv")

#Reshape to "long" format
d3<-gather(d2,I10_DX1,I10_DX2,I10_DX3,I10_DX4,I10_DX5,I10_DX6,I10_DX7,I10_DX8,
           I10_DX9,I10_DX10,I10_DX11,I10_DX12,I10_DX13,I10_DX14,I10_DX15,I10_DX16,
           I10_DX17,I10_DX18,I10_DX19,I10_DX20,I10_DX21,I10_DX22,I10_DX23,I10_DX24,
           I10_DX25,I10_DX26,I10_DX27,I10_DX28,I10_DX29,I10_DX30,I10_DX31,I10_DX32,
           I10_DX33,I10_DX34,I10_DX35,I10_DX36,I10_DX37,I10_DX38,I10_DX39,I10_DX40,
           key="original",value="icdcm")

#Identify duplicates
d3a<-group_by(d3,INC_KEY,icdcm)
d3a<-mutate(d3a,dxseq=row_number())
d3a<-mutate(d3a,dxrep=max(dxseq))
d3a<-ungroup(d3a)
#Drop duplicates (mostly NA)
d3b<-filter(d3a,dxseq==1)
d3b<-select(d3b,INC_KEY,icdcm,DIED,DISPUNIFORM,LOS)
d3b<-arrange(d3b,INC_KEY,icdcm)
View(d3b)

#Extract parts of each ICD-10-CM code
#Restrict to codes with first digit "S" or "T"
#   and seventh digit "A","B","C", or blank  
d4<-arrange(d3b,INC_KEY,icdcm)  
d4<-mutate(d4,digit1=str_sub(icdcm,1,1))
d4<-mutate(d4,digit7=str_sub(icdcm,7,7))
d4<-filter(d4,digit1=="S" | digit1=="T")
d4<-filter(d4,digit7=="A" | digit7=="B" | digit7=="C" | digit7=="")

#Verify that "DIED" corresponds correctly to "DISPUNIFORM==20"
#Rename DIED to died, drop if died not 0 or 1
#Drop if transferred to another acute-care hospital
#Create dummy TQPISS
tabyl(d4,DIED,DISPUNIFORM)
d4<-rename(d4,died=DIED)
d5<-filter(d4,died>=0)
d5<-filter(d5,DISPUNIFORM!=2)
d5<-mutate(d5,longLOS=if_else(died==1 | LOS>=10, 1,0))
d5<-mutate(d5,TQPISS=0)

#Discard temporary analytic variables and save result
d5<-select(d5,INC_KEY,icdcm,died,TQPISS)
write_csv(d5,"nisreshaped2020.csv")  


#######################################################################################  


#3
#PREPARE TQIP OR NIS DATA FOR REGRESSION ANALYSIS    
    
  d0<-read_csv("tqipmerged2020.csv")
  d0<-read_csv("nisreshaped2020.csv")  

#Inspect for otherwise valid codes not ending in "A","B", or "C"  
  d1<-filter(d0,str_sub(icdcm,7,7)=="")
  tabyl(d1,icdcm)
#There are some, especially T31x, T31xx, T32x, T32xx (Burns)  
  
#Restrict to cases with at least one valid icdcm Anatomic Injury Code
#  as defined by National Trauma Data Standard
#  or codes T20-T32 (denoting burns/corrosions)
  d4<-mutate(d0,validcode=case_when(
    str_sub(icdcm,1,1)=="S" & str_sub(icdcm,7,7)=="A" ~ 1,
    str_sub(icdcm,1,1)=="S" & str_sub(icdcm,7,7)=="B" ~ 1,
    str_sub(icdcm,1,1)=="S" & str_sub(icdcm,7,7)=="C" ~ 1,
    str_sub(icdcm,1,1)=="T" & str_sub(icdcm,2,3)=="07" ~ 1,
    str_sub(icdcm,1,1)=="T" & str_sub(icdcm,2,3)=="14" ~ 1,
    str_sub(icdcm,1,1)=="T" & str_sub(icdcm,2,2)=="2"  ~ 1,
    str_sub(icdcm,1,1)=="T" & str_sub(icdcm,2,3)=="30" ~ 1,
    str_sub(icdcm,1,1)=="T" & str_sub(icdcm,2,3)=="31" ~ 1,
    str_sub(icdcm,1,1)=="T" & str_sub(icdcm,2,3)=="32" ~ 1,
    str_sub(icdcm,1,5)=="T79.A" & str_sub(icdcm,7,7)=="A" ~ 1, 
    TRUE ~ 0
    ) )

#Inspect whether valid codes properly identified
  d4test<-group_by(d4,icdcm)
  d4test<-mutate(d4test,dxseq=row_number())
  d4test<-filter(d4test,dxseq==1)
  d4test<-ungroup(d4test)
  d4test<-arrange(d4test,icdcm)
  #View(d4test)

#Identify valid cases (with at least one valid code and died==0 or 1)
#Discard invalid cases and invalid codes  
  d5<-group_by(d4,INC_KEY)
  d5<-mutate(d5,idseq=row_number())
  d5<-mutate(d5,validcase=max(validcode))
  d5<-ungroup(d5)
  tabyl(d5,validcase,validcode)
  d6<-filter(d5,validcase==1,validcode==1)
  d6<-filter(d6,died>=0)
  tabyl(d6,validcase,validcode)
  
#Discard temporary analytic variables, and save results
  d7<-select(d6,INC_KEY,icdcm,died,TQPISS)
  write_csv(d7,"tqip2020cm.csv")
  write_csv(d7,"nis2020cm.csv")
  


  
  
  