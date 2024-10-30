##########################################################################
#
#      R PROGRAMS TO EXTRACT INJURY DATA 
#           FROM TQP PUF (FORMERLY TQIP/NTDB) AND FROM NIS
#           AND ANALYZE ROCMAX OPTIONS FOR ICDPIC-R
#
#      PART D:  MODIFY RESULTS FROM PART C 
#               COMBINE AIS RESULTS FROM TQIP AND NIS
#               EXTEND AIS RESULTS TO ICD CODES NOT IN ORIGINAL DATA
#               SMOOTH AIS RESULTS FOR SIMILAR OR GRADED CODES
#               ALLOW FOR TRUNCATED ICD CODES, INCLUDING BASIC ICD-10
#               PRODUCE FINAL TABLE FOR CALCULATING ISS
#
#      David Clark, 2022-2024
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
require(icdpicr)
require(comorbidity)


#2
#OBTAIN DATA FROM PART C AND COMBINE RESULTS TO GET WEIGHTED AIS

d6<-read_csv("tqip2020cm_ais.csv")
d7<-read_csv("nis2020cm_ais.csv")
d8<-full_join(d6,d7,by="icdcm")
d8<-rename(d8,ICD=icdcm)
d8<-mutate(d8,TQIPn=if_else(is.na(TQIPn),0,TQIPn))
d8<-mutate(d8,NISn=if_else(is.na(NISn),0,NISn))
d8<-mutate(d8,N=TQIPn+NISn)
d8<-mutate(d8,AIS0=round( (TQIPais*TQIPn + NISais*NISn)/N) )
d8<-mutate(d8,AIS0=case_when(
  is.na(AIS0) & is.na(TQIPais) ~ NISais,
  is.na(AIS0) & is.na(NISais) ~ TQIPais,
  TRUE ~ AIS0
))
d8<-mutate(d8,TQIPbr=if_else(is.na(TQIPbr),"X",TQIPbr))
d8<-mutate(d8,NISbr=if_else(is.na(NISbr),"X",NISbr))
d8<-group_by(d8,ICD)
d8<-mutate(d8,BR0=min(TQIPbr,NISbr))
d8<-ungroup(d8)
d8<-arrange(d8,ICD)

d9<-select(d8,-TQIPbr,-NISbr,-TQIPais,-TQIPn,-NISn,-NISais)
write_csv(d9,"ICDAIS_0.csv")


#3
#CREATE LIST OF ALL VALID ICD-10-CM INJURY CODES 
#Including those possibly not encountered in data

#Get descriptions of ICD-10-CM codes 
d0<-read_csv("PUF AY 2020/CSV/PUF_ICDDIAGNOSIS_LOOKUP.csv")
d1<-rename(d0,withdot=ICDDIAGNOSISCODE,description="ICDDiagnosisCode_Desc")
d1<-mutate(d1,predot=str_sub(withdot,1,3))
d1<-mutate(d1,postdot=str_sub(withdot,5,8))
d1<-mutate(d1,ICD=str_c(predot,postdot))
d1<-select(d1,ICD,description)

#Restrict to codes defined by National Trauma Data Standard
#  including fracture codes ending in B or C (denoting open fractures)
#  plus all codes T20-T32 (denoting burns/corrosions)
d1<-mutate(d1,validcode=case_when(
  str_sub(ICD,1,1)=="S" & str_sub(ICD,7,7)=="A" ~ 1,
  str_sub(ICD,1,1)=="S" & str_sub(ICD,7,7)=="B" ~ 1,
  str_sub(ICD,1,1)=="S" & str_sub(ICD,7,7)=="C" ~ 1,
  str_sub(ICD,1,1)=="T" & str_sub(ICD,2,3)=="07" ~ 1,
  str_sub(ICD,1,1)=="T" & str_sub(ICD,2,3)=="14" ~ 1,
  str_sub(ICD,1,1)=="T" & str_sub(ICD,2,2)=="2"  ~ 1,
  str_sub(ICD,1,1)=="T" & str_sub(ICD,2,3)=="30" ~ 1,
  str_sub(ICD,1,1)=="T" & str_sub(ICD,2,3)=="31" ~ 1,
  str_sub(ICD,1,1)=="T" & str_sub(ICD,2,3)=="32" ~ 1,
  str_sub(ICD,1,5)=="T79.A" & str_sub(ICD,7,7)=="A" ~ 1, 
  TRUE ~ 0
) )
d1<-mutate(d1,validcode=if_else(
  str_sub(ICD,7,7)=="D" | str_sub(ICD,7,7)=="S", 0, validcode
  ))
d3<-filter(d1,validcode==1)  
d3<-arrange(d3,ICD)

#Identify any duplicates
d3test<-group_by(d3,ICD)
d3test<-mutate(d3test,dxseq=row_number())
d3test<-mutate(d3test,dxrep=max(dxseq))
d3test<-ungroup(d3test)
tabyl(d3test,dxrep)
#No duplicates

#Extract parts of each code
d3<-mutate(d3,digit1=str_sub(ICD,1,1))
d3<-mutate(d3,digit4=str_sub(ICD,4,4))
d3<-mutate(d3,digit5=str_sub(ICD,5,5))
d3<-mutate(d3,digit6=str_sub(ICD,6,6))
d3<-mutate(d3,digit7=str_sub(ICD,7,7))
d3<-mutate(d3,digits123=str_sub(ICD,1,3))
d3<-mutate(d3,digits1234=str_sub(ICD,1,4))
d3<-mutate(d3,digits12345=str_sub(ICD,1,5))
d3<-mutate(d3,digits123456=str_sub(ICD,1,6))
d3<-mutate(d3,digits123457=str_c(digits12345,digit7,sep=""))
d3<-mutate(d3,digits123467=str_c(digits1234,digit6,digit7,sep=""))

#Add ISS body regions
d4<-mutate(d3,BR1=case_when(
  str_sub(ICD,1,3)=="S00" & str_sub(ICD,4,4)=="1" ~ "F",
  str_sub(ICD,1,3)=="S00" & str_sub(ICD,4,4)=="2" ~ "F",
  str_sub(ICD,1,3)=="S00" & str_sub(ICD,4,4)=="3" ~ "F",
  str_sub(ICD,1,3)=="S00" & str_sub(ICD,4,4)=="4" ~ "F",
  str_sub(ICD,1,3)=="S00" & str_sub(ICD,4,4)=="5" ~ "F",
  str_sub(ICD,1,3)=="S01" & str_sub(ICD,4,4)=="1" ~ "F",
  str_sub(ICD,1,3)=="S01" & str_sub(ICD,4,4)=="2" ~ "F",
  str_sub(ICD,1,3)=="S01" & str_sub(ICD,4,4)=="3" ~ "F",
  str_sub(ICD,1,3)=="S01" & str_sub(ICD,4,4)=="4" ~ "F",
  str_sub(ICD,1,3)=="S01" & str_sub(ICD,4,4)=="5" ~ "F",
  str_sub(ICD,1,3)=="S02" & str_sub(ICD,4,4)=="1" ~ "F",
  str_sub(ICD,1,3)=="S02" & str_sub(ICD,4,4)=="2" ~ "F",
  str_sub(ICD,1,3)=="S02" & str_sub(ICD,4,4)=="3" ~ "F",
  str_sub(ICD,1,3)=="S02" & str_sub(ICD,4,4)=="4" ~ "F",
  str_sub(ICD,1,3)=="S02" & str_sub(ICD,4,4)=="5" ~ "F",
  str_sub(ICD,1,3)=="S02" & str_sub(ICD,4,4)=="6" ~ "F",
  str_sub(ICD,1,3)=="S03" & str_sub(ICD,4,4)=="0" ~ "F",
  str_sub(ICD,1,3)=="S03" & str_sub(ICD,4,4)=="1" ~ "F",
  str_sub(ICD,1,3)=="S03" & str_sub(ICD,4,4)=="2" ~ "F",
  str_sub(ICD,1,3)=="S03" & str_sub(ICD,4,4)=="3" ~ "F",
  str_sub(ICD,1,3)=="S03" & str_sub(ICD,4,4)=="4" ~ "F",
  str_sub(ICD,1,3)=="S03" & str_sub(ICD,4,4)=="5" ~ "F",
  str_sub(ICD,1,3)=="S05"   ~ "F",
  str_sub(ICD,1,4)=="S070" ~ "F",
  str_sub(ICD,1,4)=="S081" ~ "F",
  str_sub(ICD,1,4)=="S092" ~ "F",
  str_sub(ICD,1,2)=="S0" ~ "H",
  str_sub(ICD,1,2)=="S1" ~ "H",
  str_sub(ICD,1,2)=="S2" ~ "C",
  str_sub(ICD,1,2)=="S3" ~ "A",
  str_sub(ICD,1,2)=="S4" ~ "E",
  str_sub(ICD,1,2)=="S5" ~ "E",
  str_sub(ICD,1,2)=="S6" ~ "E",
  str_sub(ICD,1,2)=="S7" ~ "E",
  str_sub(ICD,1,2)=="S8" ~ "E",
  str_sub(ICD,1,2)=="S9" ~ "E",
  str_sub(ICD,1,4)=="T242" ~ "E",
  str_sub(ICD,1,4)=="T79A" ~ "E",
  str_sub(ICD,1,1)=="T" ~ "G",
  TRUE ~ "N"  
) )
tabyl(d4,BR1)

write_csv(d4,"validcodes.csv")


#4
#IDENTIFY GROUPS OF CODES THAT SHOULD HAVE THE SAME SEVERITY
#  Some codes for intracranial injury specify length of coma or vital outcomes (survival/death)
#     This is inappropriate for a predictive anatomic severity scale 
#     Therefore these codes are grouped so they can be assigned the same severity
#  Some codes specify "right", "left", or "unspecified side"
#     These are also grouped so that they can be assigned the same severity
#     (Unfortunately, ICD-10-CM has many ways of specifying laterality)

#Read in list of valid codes created above
d0<-read_csv("validcodes.csv")

#Identify intracranial injuries describing prolonged loss of consciousness (LOC)
#Note: In some cases, these also specify survival/death
d1<-mutate(d0,LOC=if_else(digits123=="S06"
            & (digit6=="2" | digit6=="3" | digit6=="4" | digit6=="5"              
              | digit6=="6" | digit6=="7" | digit6=="8" | digit6=="9") , 1, 0))

#Group LOC injuries with prolonged coma, whether left-, right-, or unspecified-sided
d2<-mutate(d1,LOCgroup=case_when(
   LOC==1 & digits1234=="S060" ~ "S060",
   LOC==1 & digits1234=="S061" ~ "S061",
   LOC==1 & digits1234=="S062" ~ "S062",
   LOC==1 & digits12345=="S0630" ~ "S063_0",
   LOC==1 & digits12345=="S0631" ~ "S063_123",
   LOC==1 & digits12345=="S0632" ~ "S063_123",
   LOC==1 & digits12345=="S0633" ~ "S063_123",
   LOC==1 & digits12345=="S0634" ~ "S063_456",
   LOC==1 & digits12345=="S0635" ~ "S063_456",
   LOC==1 & digits12345=="S0636" ~ "S063_456",
   LOC==1 & digits12345=="S0637" ~ "S063_7",
   LOC==1 & digits12345=="S0638" ~ "S063_8",
   LOC==1 & digits1234=="S064" ~ "S064",
   LOC==1 & digits1234=="S065" ~ "S065",
   LOC==1 & digits1234=="S066" ~ "S066",
   LOC==1 & digits12345=="S0681" ~ "S068_12",
   LOC==1 & digits12345=="S0682" ~ "S068_12",
   LOC==1 & digits12345=="S0689" ~ "S068_9",
   LOC==1 & digits1234=="S069" ~ "S069",
   TRUE ~ "X"
))
tabyl(d2,LOCgroup,LOC)

#Identify other injuries specifying right-, left-, or unspecified-sided
d4<-mutate(d2,hand=case_when(
  str_detect(description,"right")==TRUE ~ "R",
  str_detect(description,"left")==TRUE ~ "L",
  str_detect(description,"unspecified") 
    & str_detect(description,"right")==FALSE & str_detect(description,"left")==FALSE ~ "U",
  TRUE ~ "X"
))
d4<-mutate(d4,RLgroup=if_else(LOC==1,"Z","X"))

d5<-group_by(d4,digits123457)
d5<-mutate(d5,test6_120=case_when(
  hand=="R" & digit6=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit6=="2" & RLgroup=="X" ~ 1,
  hand=="U" & digit6=="0" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d5<-mutate(d5,sum6_120=sum(test6_120))
d5<-mutate(d5,RLgroup=if_else((sum6_120==3)
                         & (digit6=="1" | digit6=="2" | digit6=="0"), "RLU6_120", RLgroup))

d5<-mutate(d5,test6_123=case_when(
  hand=="R" & digit6=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit6=="2" & RLgroup=="X" ~ 1,
  hand=="U" & digit6=="3" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d5<-mutate(d5,sum6_123=sum(test6_123))
d5<-mutate(d5,RLgroup=if_else((sum6_123==3)
                         & (digit6=="1" | digit6=="2" | digit6=="3"), "RLU6_123", RLgroup))

d5<-mutate(d5,test6_456=case_when(
  hand=="R" & digit6=="4" & RLgroup=="X" ~ 1,
  hand=="L" & digit6=="5" & RLgroup=="X" ~ 1,
  hand=="U" & digit6=="6" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d5<-mutate(d5,sum6_456=sum(test6_456))
d5<-mutate(d5,RLgroup=if_else((sum6_456==3)
                         & (digit6=="4" | digit6=="5" | digit6=="6"), "RLU6_456", RLgroup))

d5<-mutate(d5,test6_789=case_when(
  hand=="R" & digit6=="7" & RLgroup=="X" ~ 1,
  hand=="L" & digit6=="8" & RLgroup=="X" ~ 1,
  hand=="U" & digit6=="9" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d5<-mutate(d5,sum6_789=sum(test6_789))
d5<-mutate(d5,RLgroup=if_else((sum6_789==3)
                         & (digit6=="7" | digit6=="8" | digit6=="9"), "RLU6_789", RLgroup))

d5<-mutate(d5,test6_129=case_when(
  hand=="R" & digit6=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit6=="2" & RLgroup=="X" ~ 1,
  hand=="U" & digit6=="9" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d5<-mutate(d5,sum6_129=sum(test6_129))
d5<-mutate(d5,RLgroup=if_else((sum6_129==3)
                         & (digit6=="1" | digit6=="2" | digit6=="9"), "RLU6_129", RLgroup))

d5<-mutate(d5,test6_12=case_when(
  hand=="R" & digit6=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit6=="2" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d5<-mutate(d5,sum6_12=sum(test6_12))
d5<-mutate(d5,RLgroup=if_else((sum6_12==2)
                         & (digit6=="1" | digit6=="2"), "RL6_12", RLgroup))

d5<-ungroup(d5)
tabyl(d5,RLgroup,hand)

d6<-group_by(d5,digits123467)
d6<-mutate(d6,test5_120=case_when(
  hand=="R" & digit5=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit5=="2" & RLgroup=="X" ~ 1,
  hand=="U" & digit5=="0" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d6<-mutate(d6,sum5_120=sum(test5_120))
d6<-mutate(d6,RLgroup=if_else((sum5_120==3)
                         & (digit5=="1" | digit5=="2" | digit5=="0"), "RLU5_120", RLgroup))

d6<-mutate(d6,test5_123=case_when(
  hand=="R" & digit5=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit5=="2" & RLgroup=="X" ~ 1,
  hand=="U" & digit5=="3" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d6<-mutate(d6,sum5_123=sum(test5_123))
d6<-mutate(d6,RLgroup=if_else((sum5_123==3)
                         & (digit5=="1" | digit5=="2" | digit5=="3"), "RLU5_123", RLgroup))

d6<-mutate(d6,test5_456=case_when(
  hand=="R" & digit5=="4" & RLgroup=="X" ~ 1,
  hand=="L" & digit5=="5" & RLgroup=="X" ~ 1,
  hand=="U" & digit5=="6" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d6<-mutate(d6,sum5_456=sum(test5_456))
d6<-mutate(d6,RLgroup=if_else((sum5_456==3)
                         & (digit5=="4" | digit5=="5" | digit5=="6"), "RLU5_456", RLgroup))

d6<-mutate(d6,test5_789=case_when(
  hand=="R" & digit5=="7" & RLgroup=="X" ~ 1,
  hand=="L" & digit5=="8" & RLgroup=="X" ~ 1,
  hand=="U" & digit5=="9" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d6<-mutate(d6,sum5_789=sum(test5_789))
d6<-mutate(d6,RLgroup=if_else((sum5_789==3)
                         & (digit5=="7" | digit5=="8" | digit5=="9"), "RLU5_789", RLgroup))

d6<-mutate(d6,test5_129=case_when(
  hand=="R" & digit5=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit5=="2" & RLgroup=="X" ~ 1,
  hand=="U" & digit5=="9" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d6<-mutate(d6,sum5_129=sum(test5_129))
d6<-mutate(d6,RLgroup=if_else((sum5_129==3)
                         & (digit5=="1" | digit5=="2" | digit5=="9"), "RLU5_129", RLgroup))

d6<-mutate(d6,test5_12=case_when(
  hand=="R" & digit5=="1" & RLgroup=="X" ~ 1,
  hand=="L" & digit5=="2" & RLgroup=="X" ~ 1,
  TRUE ~ 0
))
d6<-mutate(d6,sum5_12=sum(test5_12))
d6<-mutate(d6,RLgroup=if_else((sum5_12==2)
                         & (digit5=="1" | digit5=="2"), "RL5_12", RLgroup))

d6<-ungroup(d6)
tabyl(d6,RLgroup,hand)  
tabyl(d6,RLgroup,LOC)

d7<-select(d6,ICD,description,BR1,LOC,LOCgroup,hand,RLgroup)
           
write_csv(d7,"ICD_categories.csv")


##################################################################################

#5
# ADD DATA FROM PARTS A-C TO LIST OF ALL VALID CODES
# UNIFY AIS BY GROUPS SPECIFIED ABOVE

d1<-read_csv("ICDAIS_0.csv")
d2<-read_csv("ICD_categories.csv")
d3<-full_join(d1,d2,by="ICD")

#Extract parts of each code
d3<-mutate(d3,digit1=str_sub(ICD,1,1))
d3<-mutate(d3,digit4=str_sub(ICD,4,4))
d3<-mutate(d3,digit5=str_sub(ICD,5,5))
d3<-mutate(d3,digit6=str_sub(ICD,6,6))
d3<-mutate(d3,digit7=str_sub(ICD,7,7))
d3<-mutate(d3,digits123=str_sub(ICD,1,3))
d3<-mutate(d3,digits1234=str_sub(ICD,1,4))
d3<-mutate(d3,digits12345=str_sub(ICD,1,5))
d3<-mutate(d3,digits123456=str_sub(ICD,1,6))
d3<-mutate(d3,digits123457=str_c(digits12345,digit7,sep=""))
d3<-mutate(d3,digits123467=str_c(digits1234,digit6,digit7,sep=""))

d3<-arrange(d3,ICD)
d3<-mutate(d3,AIS0=if_else(is.na(AIS0),0,AIS0))
d3<-mutate(d3,N=if_else(is.na(N),0,N))
tabyl(d3,BR0,BR1)
d3<-mutate(d3,BR1=if_else(is.na(BR1),BR0,BR1))
tabyl(d3,BR0,BR1)
d3<-mutate(d3,LOCgroup=if_else(is.na(LOCgroup),"X",LOCgroup))
d3<-mutate(d3,RLgroup=if_else(is.na(RLgroup),"X",RLgroup))

d4<-group_by(d3,digits1234,LOCgroup)
d4<-mutate(d4,NN=sum(N))
d4<-mutate(d4,wAIS=((N/NN)*AIS0))
d4<-mutate(d4,AIS1=if_else(
    LOCgroup=="S060" | LOCgroup=="S061" | LOCgroup=="S062" | LOCgroup=="S064" |
    LOCgroup=="S065" | LOCgroup=="S066" | LOCgroup=="S069", round(sum(wAIS)),AIS0))
d4<-ungroup(d4)

d4<-group_by(d4,digits12345,LOCgroup)
d4<-mutate(d4,NN=sum(N))
d4<-mutate(d4,wAIS=((N/NN)*AIS0))
d4<-mutate(d4,AIS1=if_else(
    LOCgroup=="S063_0" | LOCgroup=="S063_123" | LOCgroup=="S063_456" | LOCgroup=="S063_7" |
    LOCgroup=="S063_8" | LOCgroup=="S068_12" | LOCgroup=="S068_9", round(sum(wAIS)),AIS1))

d5<-group_by(d4,RLgroup,digits123467)
d5<-mutate(d5,NN=sum(N))
d5<-mutate(d5,wAIS=((N/NN)*AIS0))
d5<-mutate(d5,AIS1=if_else(
    substr(RLgroup,1,3)=="RL5" | substr(RLgroup,1,4)=="RLU5", round(sum(wAIS)),AIS1))
d5<-ungroup(d5)

d6<-group_by(d5,RLgroup,digits123457)
d6<-mutate(d6,NN=sum(N))
d6<-mutate(d6,wAIS=((N/NN)*AIS0))
d6<-mutate(d6,AIS1=if_else(
  substr(RLgroup,1,3)=="RL6" | substr(RLgroup,1,4)=="RLU6", round(sum(wAIS)),AIS1))
d6<-ungroup(d6)
d6<-arrange(d6,ICD)

write_csv(d6,"ICDAIS_1.csv")



#6
#IDENTIFY GROUPS OF CODES THAT SHOULD HAVE GRADED SEVERITY

#  Some codes for moderate intracranial injury specify brief coma
#     These are ranked so AIS increases from no coma to brief coma to prolonged coma
#  Burn diagnoses specifying surface area are similaRLgroupy ranked.
#     as are fracture diagnoses specifying severe open vs simple open vs closed 


d0<-read_csv("ICDAIS_1.csv")

d1<-mutate(d0,AIS1=if_else(is.na(AIS1) | AIS1==0, 1, AIS1))
d1<-mutate(d1,AIS2=AIS1)

d1<-mutate(d1,LOCgroup012=if_else(digits123=="S06"
                          & (digit6=="0" | digit6=="1" | digit6=="2") , 1, 0))
d1<-group_by(d1,digits12345,LOCgroup012)
d1<-mutate(d1,minais=min(AIS2))
d1<-mutate(d1,maxais=max(AIS2))
d1<-mutate(d1,AIS2=case_when(
    LOCgroup012==1 & digit6=="0" ~ minais,
    LOCgroup012==1 & digit6=="1" ~ round((minais+maxais)/2),
    TRUE ~ AIS2
    ))

d2<-mutate(d1,Burngroup=case_when(
    digits123=="T31" & digit4=="0" ~ "T31_0",
    digits123=="T31" & digit4=="1" ~ "T31_1",
    digits123=="T31" & digit4=="2" ~ "T31_2",
    digits123=="T31" & digit4=="3" ~ "T31_3",
    digits123=="T31" & digit4=="4" ~ "T31_4",
    digits123=="T31" & digit4=="5" ~ "T31_5",
    digits123=="T31" & digit4=="6" ~ "T31_6",
    digits123=="T31" & digit4=="7" ~ "T31_7",
    digits123=="T31" & digit4=="8" ~ "T31_8",
    digits123=="T31" & digit4=="9" ~ "T31_9",
    digits123=="T32" & digit4=="0" ~ "T32_0",
    digits123=="T32" & digit4=="1" ~ "T32_1",
    digits123=="T32" & digit4=="2" ~ "T32_2",
    digits123=="T32" & digit4=="3" ~ "T32_3",
    digits123=="T32" & digit4=="4" ~ "T32_4",
    digits123=="T32" & digit4=="5" ~ "T32_5",
    digits123=="T32" & digit4=="6" ~ "T32_6",
    digits123=="T32" & digit4=="7" ~ "T32_7",
    digits123=="T32" & digit4=="8" ~ "T32_8",
    digits123=="T32" & digit4=="9" ~ "T32_9",
    TRUE ~ "X"
    ))
d2<-group_by(d2,Burngroup)
d2<-mutate(d2,NN=sum(N))
d2<-mutate(d2,NN=if_else(is.na(NN) | NN==0,1,NN))
d2<-mutate(d2,wAIS=((N/NN)*AIS1))
d2<-mutate(d2,AIS2=if_else(Burngroup!="X",round(sum(wAIS)),AIS2))
d2<-ungroup(d2)

d3<-group_by(d2,digits123456)
d3<-mutate(d3,test7_ABC=case_when(
  digit7=="A" ~ 1,
  digit7=="B" ~ 1,
  digit7=="C" ~ 1,
  TRUE ~ 0
))
d3<-mutate(d3,sum_ABC=sum(test7_ABC))
d3<-mutate(d3,Fxgroup=if_else( (sum_ABC==2)
                & (digit7=="A" | digit7=="B"), "FxAB", "X"))
d3<-mutate(d3,Fxgroup=if_else( (sum_ABC==3)
                & (digit7=="A" | digit7=="B" | digit7=="C"), "FxABC", Fxgroup))
d3<-ungroup(d3)

d4<-group_by(d3,digits123456,Fxgroup)
d4<-mutate(d4,NN=sum(N))
d4<-mutate(d4,NN=if_else(is.na(NN) | NN==0,1,NN))
d4<-mutate(d4,wAIS=((N/NN)*AIS1))
d4<-mutate(d4,maxais=max(AIS1))
d4<-mutate(d4,AIS2=case_when(
    Fxgroup=="FxAB" & digit7=="A" ~ round(sum(wAIS)),
    Fxgroup=="FxAB" & digit7=="B" ~ maxais,
    Fxgroup=="FxABC" & digit7=="A" ~ round(sum(wAIS)),
    Fxgroup=="FxABC" & digit7=="C" ~ maxais,
    Fxgroup=="FxABC" & digit7=="B" ~ round( (sum(wAIS)+maxais)/2 ),
    TRUE ~ AIS2
    ))
d4<-ungroup(d4)

d5<-mutate(d4,description=if_else(is.na(description),"( )",description))
d5<-mutate(d5,AIS2=if_else(
  str_detect(description,"uperficial")==TRUE, 1, AIS2))
d5<-mutate(d5,AIS2=if_else(AIS2==0,1,AIS2))
write_csv(d5,"ICDAIS_2.csv")
  

#7
#ASSIGN AIS TO TRUNCATED CODES (INCLUDING BASIC ICD-10 CODES)
#AIS severity as weighted average of contributing codes 

d1<-read_csv("ICDAIS_2.csv")

d2<-group_by(d1,digits123456)
d2<-mutate(d2,NN=sum(N))
d2<-mutate(d2,NN=if_else(is.na(NN) | NN==0,1,NN))
d2<-mutate(d2,wAIS=((N/NN)*AIS2))
d2<-mutate(d2,AIS3=round(sum(wAIS)))
d2<-mutate(d2,BR3=min(BR1))
d2<-mutate(d2,seq=row_number())
d2<-ungroup(d2)

d123456<-filter(d2,seq==1,Burngroup=="X")
d123456<-select(d123456,digits123456,AIS3,BR3)
d123456<-rename(d123456,ICD=digits123456,AIS=AIS3,BR=BR3)

d3<-group_by(d1,digits12345)
d3<-mutate(d3,NN=sum(N))
d3<-mutate(d3,NN=if_else(is.na(NN) | NN==0,1,NN))
d3<-mutate(d3,wAIS=((N/NN)*AIS2))
d3<-mutate(d3,AIS3=round(sum(wAIS)))
d3<-mutate(d3,BR3=min(BR1))
d3<-mutate(d3,seq=row_number())
d3<-ungroup(d3)

d12345<-filter(d3,seq==1,Burngroup=="X")
d12345<-select(d12345,digits12345,AIS3,BR3)
d12345<-rename(d12345,ICD=digits12345,AIS=AIS3,BR=BR3)

d4<-group_by(d1,digits1234)
d4<-mutate(d4,NN=sum(N))
d4<-mutate(d4,NN=if_else(is.na(NN) | NN==0,1,NN))
d4<-mutate(d4,wAIS=((N/NN)*AIS2))
d4<-mutate(d4,AIS3=round(sum(wAIS)))
d4<-mutate(d4,BR3=min(BR1))
d4<-mutate(d4,seq=row_number())
d4<-ungroup(d4)

d1234<-filter(d4,seq==1,Burngroup=="X")
d1234<-select(d1234,digits1234,AIS3,BR3)
d1234<-rename(d1234,ICD=digits1234,AIS=AIS3,BR=BR3)

d5<-group_by(d1,digits123)
d5<-mutate(d5,NN=sum(N))
d5<-mutate(d5,NN=if_else(is.na(NN) | NN==0,1,NN))
d5<-mutate(d5,wAIS=((N/NN)*AIS2))
d5<-mutate(d5,AIS3=round(sum(wAIS)))
d5<-mutate(d5,BR3=min(BR1))
d5<-mutate(d5,seq=row_number())
d5<-ungroup(d5)

d123<-filter(d5,seq==1,Burngroup=="X")
d123<-select(d123,digits123,AIS3,BR3)
d123<-rename(d123,ICD=digits123,AIS=AIS3,BR=BR3)



####################### FOLLOWING CODE MODIFIED OCTOBER 2024 ###############################

write_csv(d123,"ICDAIS_2A.csv")

d123 <- read_csv("ICDAIS_2A.csv")

d1renamed<-rename(d1,BR=BR1,AIS=AIS2)
d6<-bind_rows(d123,d1234,d12345,d123456,d1renamed)
d6<-mutate(d6,AIS=if_else(AIS==0,1,AIS))
d6<-arrange(d6,ICD)
d6<-select(d6,ICD,AIS,BR,TQIPeffect,TQIPint,NISeffect,NISint)

write_csv(d6,"ICDAIS_4.csv")


#Get list of valid ICD10 (international) codes not obtainable by truncation
#Table modified from Gedeborg et al., Journal of Trauma 2014
#  Excludes diagnoses outside National Trauma Data Standard (see icd10aisA)
#  Adds a few codes (S137,S577,S732,S737,S738) 
#  Indicates codes already in ICDAIS_4.csv
#  Assigns body region

dged <- read_csv("/Users/davideugeneclark/Documents/icdpicr2/gedeborg_modified.csv")
dged <- rename(dged,ICD=icd10)
dged <- rename(dged,BR=br)
dged <- mutate(dged,mort=(totaln-survivors)/totaln)

#Assign AIS, somewhat arbitrarily, based on diagnosis-specific mortality 
#  and compare to ais from truncation already in ICDAIS_4.csv

dged <- mutate(dged,AIS=case_when(
   mort>=0 & mort<.01 ~ 1,
   mort>=.01 & mort<.02 ~ 2,
   mort>=.02 & mort<.1 ~ 3,
   mort>=.1 & mort<.2 ~ 4,
   mort>=.2 & mort<=1 ~ 5,
   TRUE ~ 1
   ) )
tabyl(dged,AIS,ais)
tabyl(dged,AIS,added)

dged <- select(dged,ICD,AIS,BR,already_in,totaln,dsp)

#Add remaining diagnoses to ICDAIS_4.csv
d6 <- read_csv("/Users/davideugeneclark/Documents/icdpicr/ICDAIS_4.csv")
d6 <- bind_rows(d6,dged)
write_csv(d6,"/Users/davideugeneclark/Documents/icdpicr2/ICDAIS_5.csv")


#Unduplicate to create ICD_AIS.csv
#Sort in such a way that information is used from Gedeborg if available
d7<-read_csv("/Users/davideugeneclark/Documents/icdpicr2/ICDAIS_5.csv")
d7<-arrange(d7,ICD,already_in)
d7<-group_by(d7,ICD)
d7<-mutate(d7,seq=row_number())
d7<-ungroup(d7)
d7<-filter(d7,seq==1)
d7<-mutate(d7,version="v241023")
d7<-select(d7,-seq,-already_in,-totaln,-dsp)
write_csv(d7,"/Users/davideugeneclark/Documents/icdpicr2/ICD_AIS_241023.csv")

#Unduplicate to create i10_map_iciss.csv
d8<-read_csv("/Users/davideugeneclark/Documents/icdpicr2/ICDAIS_5.csv")
d8<-arrange(d8,ICD,already_in)
d8<-group_by(d8,ICD)
d8<-select(d8,-TQIPeffect,-TQIPint,-NISeffect,-NISint)
d8<-mutate(d8,seq=row_number())
d8<-ungroup(d8)
d8<-filter(d8,seq==1)
d8<-mutate(d8,base=str_sub(ICD,1,4))
d8<-group_by(d8,base)
d8<-mutate(d8,totaln=max(totaln,na.rm=TRUE))
d8<-mutate(d8,dsp=max(dsp,na.rm=TRUE))
d8<-ungroup(d8)
d8<-mutate(d8,totaln=if_else((totaln<0),NA,totaln))
d8<-mutate(d8,dsp=if_else((dsp<0|dsp>1),NA,dsp))
d8<-mutate(d8,dsp_cons=if_else(totaln<5,1,dsp))
d8<-select(d8,ICD,totaln,dsp,dsp_cons)
d8<-rename(d8,dsp_noncons=dsp)
d8<-rename(d8,dx=ICD)
d8<-mutate(d8,version="v241025")
write_csv(d8,"/Users/davideugeneclark/Documents/icdpicr2/i10_map_iciss_241025.csv")


        