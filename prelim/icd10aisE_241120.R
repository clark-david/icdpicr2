##########################################################################
#
#      R PROGRAMS TO EXTRACT INJURY DATA 
#           FROM TQP PUF (FORMERLY TQIP/NTDB) AND FROM NIS
#           AND ANALYZE ROCMAX OPTIONS FOR ICDPIC-R
#
#      PART E:  MODIFY RESULTS FROM PART D 
#               UPDATE EXTERNAL CAUSE CATEGORIES
#               DERIVE "BARELL MATRIX" CATEGORIES
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





#UPDATE ICD-10 MECHANISM CODE TABLE

#Combine latest CDC External Cause tables for ICD-10-CM

d1a<-read_csv("CDC2021External_Unintentional.csv")
d1b<-select(d1a,ICD10CM,MECHANISM)
d1b<-mutate(d1b,INTENT="Unintentional")
d2a<-read_csv("CDC2021External_Assault.csv")
d2b<-select(d2a,ICD10CM,MECHANISM)
d2b<-mutate(d2b,INTENT="Assault")
d3a<-read_csv("CDC2021External_SelfHarm.csv")
d3b<-select(d3a,ICD10CM,MECHANISM)
d3b<-mutate(d3b,INTENT="Self-inflicted")
d4a<-read_csv("CDC2021External_LegalWar.csv")
d4b<-select(d4a,ICD10CM,MECHANISM)
d4b<-mutate(d4b,INTENT="Legalintervention/War")
d5a<-read_csv("CDC2021External_Undetermined.csv")
d5b<-select(d5a,ICD10CM,MECHANISM)
d5b<-mutate(d5b,INTENT="Undetermined")
d6<-bind_rows(d1b,d2b,d3b,d4b,d5b)
d6<-filter(d6,is.na(ICD10CM)==0)
tabyl(d6,MECHANISM)
d6<-mutate(d6,MECHANISM=case_when(
  MECHANISM=="Bites/Stings, venomous" ~ "Bites and Stings, venomous",
  MECHANISM=="Natural Environmental, Other" ~ "Natural/Environmental, Other",
  TRUE ~ MECHANISM
))
tabyl(d6,MECHANISM)
d6<-rename(d6,ICD10=ICD10CM)
d6<-arrange(d6,ICD10)
write_csv(d6,"mech6.csv")


#ASSIGN "BARELL MATRIX" TYPE CATEGORIES

#Obtain latest CDC matrix and simplify
d1<-read_csv("CDC2021Injury_Matrix.csv")
d2<-rename(d1,BL1=Body_Level_1,BL2=Body_Level_2,BL3=Body_Level_3)
d2<-mutate(d2,BL1=str_sub(BL1,1,3))
d2<-mutate(d2,BL2=str_sub(BL2,1,3))
d2<-mutate(d2,BL=case_when(
  BL1=="Hea" & BL2=="Tra" ~ "H+N/TBI",
  BL1=="Hea" & BL2=="Oth" & BL3=="Face" ~ "H+N/Oth/Face",
  BL1=="Hea" & BL2=="Oth" & BL3=="Eye" ~ "H+N/Oth/Eyes",
  BL1=="Hea" & BL2=="Oth" & BL3=="Other head" ~ "H+N/Oth/Head",
  BL1=="Hea" & BL2=="Oth" & BL3=="Neck" ~ "H+N/Oth/Neck",
  BL1=="Hea" & BL2=="Oth" & BL3=="Head and Neck, Other" ~ "H+N/Oth/HdNk",
  BL1=="Spi" & BL2=="Spi" & BL3=="Cervical SCI" ~ "SpB/SCI/Cerv",
  BL1=="Spi" & BL2=="Spi" & BL3=="Thoracic or dorsal SCI" ~ "SpB/SCI/Thor",
  BL1=="Spi" & BL2=="Spi" & BL3=="Lumbar SCI" ~ "SpB/SCI/Lumb",
  BL1=="Spi" & BL2=="Spi" & BL3=="Sacral/coccygeal SCI" ~ "SpB/SCI/Sacr",
  BL1=="Spi" & BL2=="Ver" & BL3=="Cervical VCI" ~ "SpB/VCI/Cerv",
  BL1=="Spi" & BL2=="Ver" & BL3=="Thoracic or dorsal VCI" ~ "SpB/VCI/Thor",
  BL1=="Spi" & BL2=="Ver" & BL3=="Lumbar VCI" ~ "SpB/VCI/Lumb",
  BL1=="Spi" & BL2=="Ver" & BL3=="Sacral/coccygeal VCI" ~ "SpB/VCI/Sacr",
  BL1=="Tor" & BL2=="Che" ~ "Tor/Thor",
  BL1=="Tor" & BL2=="Abd" ~ "Tor/Abdo",
  BL1=="Tor" & BL2=="Pel" & BL3=="External genitalia" ~ "Tor/Pel/Geni",
  BL1=="Tor" & BL2=="Pel" & BL3=="Pelvic organs" ~ "Tor/Pel/Orgs",
  BL1=="Tor" & BL2=="Pel" & BL3=="Lower back and pelvis" ~ "Tor/Pel/Back",
  BL1=="Tor" & BL2=="Pel" & BL3=="Pelvic girdle" ~ "Tor/Pel/Gird",
  BL1=="Tor" & BL2=="Pel" & BL3=="Buttock" ~ "Tor/Pel/Butt",
  BL1=="Tor" & BL2=="Pel" & BL3=="Other" ~ "Tor/Pel/Othe",
  BL1=="Tor" & BL2=="Oth" ~ "Tor/Othe",
  BL1=="Ext" & BL2=="Upp" & BL3=="Shoulder and upper arm" ~ "Ext/Upp/Shou",
  BL1=="Ext" & BL2=="Upp" & BL3=="Forearm and elbow" ~ "Ext/Upp/Fore",
  BL1=="Ext" & BL2=="Upp" & BL3=="Wrist, hand, and fingers" ~ "Ext/Upp/Hand",
  BL1=="Ext" & BL2=="Upp" & BL3=="Arm, not further specified" ~ "Ext/Upp/Arms",
  BL1=="Ext" & BL2=="Low" & BL3=="Hip" ~ "Ext/Low/Hips",
  BL1=="Ext" & BL2=="Low" & BL3=="Upper leg and thigh" ~ "Ext/Low/Thig",
  BL1=="Ext" & BL2=="Low" & BL3=="Knee" ~ "Ext/Low/Knee",
  BL1=="Ext" & BL2=="Low" & BL3=="Lower leg and ankle" ~ "Ext/Low/Legs",
  BL1=="Ext" & BL2=="Low" & BL3=="Foot and toes" ~ "Ext/Low/Foot ",
  BL1=="Ext" & BL2=="Low" & BL3=="Ankle and foot" ~ "Ext/Low/Ankl",
  BL1=="Ext" & BL2=="Low" & BL3=="Other, multiple, and unspecified" ~ "Ext/Low/Othe",
  BL1=="Unc" & BL2=="Mul" ~ "Unc/Mul",
  BL1=="Unc" & BL2=="Sys" ~ "Unc/Sys",
  BL1=="Uns" ~ "Uns",
  TRUE ~ "X"
))

d2test<-filter(d2,BL=="X")
d2test

d3<-rename(d2,NL1=Nature_Level_1,NL2=Nature_Level_2)
d3<-mutate(d3,NL1=str_sub(NL1,1,3))
d3<-mutate(d3,NL=case_when(
  NL1=="Fra" ~ "Fra",
  NL1=="Dis" ~ "Dis",
  NL1=="Int" ~ "Org",
  NL1=="Ope" ~ "Wnd",
  NL1=="Amp" ~ "Amp",
  NL1=="Blo" ~ "Vas",
  NL1=="Sup" ~ "Sup",
  NL1=="Cru" ~ "Csh",
  NL1=="Bur" & NL2=="Burns" ~ "Bur",
  NL1=="Bur" & NL2=="Corrosions" ~ "Cor",
  NL1=="Eff" ~ "FBo",
  NL1=="Oth" & is.na(NL2) ~ "OtE",
  NL1=="Poi" ~ "Poi",
  NL1=="Tox" ~ "Tox",
  NL1=="Oth" & NL2=="Sprains and strains" ~ "Spr",
  NL1=="Oth" & NL2=="Nerves" ~ "Ner",
  NL1=="Oth" & NL2=="Muscles and tendons" ~ "Mus",
  NL1=="Oth" & NL2=="Other injury" ~ "OtI",
  NL1=="Uns" ~ "Uns",
  TRUE ~ "X"
))

d3test<-filter(d3,NL=="X")
d3test

d4<-filter(d3,NL!="X")
d4<-arrange(d4,ICD_Full)
d4<-mutate(d4,cell=str_c(BL,"//",NL))
d4<-mutate(d4,digits123456=str_c(str_sub(ICD_Full,1,3),str_sub(ICD_Full,5,7)))
d4<-select(d4,digits123456,cell)
write_csv(d4,"barell6digits.csv")
tabyl(d4,cell)


#Determine mortality by cell

d6<-read_csv("tqip2020cm.csv")
#OR
d6<-read_csv("nis2020cm.csv")

d4<-read_csv("barell6digits.csv")

d6<-mutate(d6,digits123456=str_sub(icdcm,1,6)) 

d7<-full_join(d4,d6,by="digits123456")
d7<-group_by(d7,INC_KEY,digits123456)
d7<-mutate(d7,seq=row_number())
d7<-ungroup(d7)
tabyl(d7,seq)

d8<-filter(d7,seq==1,is.na(icdcm)==FALSE,is.na(died)==FALSE)
d8<-group_by(d8,cell)
d8<-mutate(d8,cellseq=row_number())
d8<-mutate(d8,n=max(cellseq))
d8<-mutate(d8,Pm=sum(died)/n)
d8<-ungroup(d8)

d8test=filter(d8,cellseq==1)
d8test=select(d8test,cell,Pm,n)
#write_csv(d8test,"PmortNIScell.csv")
#write_csv(d8test,"PmortTQPcell.csv")

#Obtain weighted mean cell mortality
d11<-read_csv("PmortNIScell.csv")
d12<-read_csv("PmortTQPcell.csv")
d13<-full_join(d11,d12,by="cell")
d13<-mutate(d13,Pm.x=if_else(is.na(Pm.x),0,Pm.x))
d13<-mutate(d13,Pm.y=if_else(is.na(Pm.y),0,Pm.y))
d13<-mutate(d13,n.x=if_else(is.na(n.x),1,n.x))
d13<-mutate(d13,n.y=if_else(is.na(n.y),1,n.y))
d13<-mutate(d13,PmCell=(Pm.x*n.x+Pm.y*n.y)/(n.x+n.y))
d13<-select(d13,cell,PmCell)
write_csv(d13,"PmortCell")

#Allow truncated codes, and then merge with list of all ICD codes
d1<-read_csv("barell6digits.csv")
d2<-read_csv("PmortCell")
d3<-full_join(d1,d2,by="cell")

d45<-mutate(d3,digits12345=str_sub(digits123456,1,5))
d45<-group_by(d45,digits12345)
d45<-mutate(d45,cell=if_else(max(cell)==min(cell),max(cell),"Uns//Uns"))
d45<-mutate(d45,PmCell=if_else(max(PmCell)==min(PmCell),max(PmCell),0.0520))
d45<-mutate(d45,seq=row_number())
d45<-ungroup(d45)
d45<-filter(d45,seq==1)
d45<-mutate(d45,digits123456=digits12345)
d45<-select(d45,digits123456,cell,PmCell)
write_csv(d45,"barell5digits.csv")
d44<-mutate(d3,digits1234=str_sub(digits123456,1,4))
d44<-group_by(d44,digits1234)
d44<-mutate(d44,cell=if_else(max(cell)==min(cell),max(cell),"Uns/Uns"))
d44<-mutate(d44,PmCell=if_else(max(PmCell)==min(PmCell),max(PmCell),0.0520))
d44<-mutate(d44,seq=row_number())
d44<-ungroup(d44)
d44<-filter(d44,seq==1)
d44<-mutate(d44,digits123456=digits1234)
d44<-select(d44,digits123456,cell,PmCell)
write_csv(d44,"barell4digits.csv")
d43<-mutate(d3,digits123=str_sub(digits123456,1,3))
d43<-group_by(d43,digits123)
d43<-mutate(d43,cell=if_else(max(cell)==min(cell),max(cell),"Uns/Uns"))
d43<-mutate(d43,PmCell=if_else(max(PmCell)==min(PmCell),max(PmCell),0.0520))
d43<-mutate(d43,seq=row_number())
d43<-ungroup(d43)
d43<-filter(d43,seq==1)
d43<-mutate(d43,digits123456=digits123)
d43<-select(d43,digits123456,cell,PmCell)
write_csv(d43,"barell3digits.csv")

d3456<-bind_rows(d3,d45,d44,d43)
d3456<-group_by(d3456,digits123456)
d3456<-mutate(d3456,seq=row_number())
d3456<-ungroup(d3456)
d3456<-filter(d3456,seq==1)
d3456<-select(d3456,-seq)
d3456<-arrange(d3456,digits123456)
write_csv(d3456,"barell3456digits.csv")

#Merge with full list of ICD-10 codes

d3456=read_csv("barell3456digits.csv")
d5<-read_csv("ICD_AIS.csv")
d5<-mutate(d5,digits123456=str_sub(ICD,1,6))
d5<-select(d5,ICD,digits123456)

d6<-full_join(d5,d3456,by="digits123456")
d6<-filter(d6,is.na(ICD)==FALSE)
d6<-arrange(d6,ICD)
d6<-mutate(d6,cell=if_else(is.na(cell),"Unk/Unk",cell))
d6<-mutate(d6,PmCell=if_else(is.na(PmCell),0.0520,PmCell))
d6<-select(d6,-digits123456)

write_csv(d6,"ICD_Cells.csv")


#COMORBIDITY SCORES (using package=comorbidity)

d1<-read_csv("NISraw2020.csv")
d2<-gather(d1,I10_DX1,I10_DX2,I10_DX3,I10_DX4,I10_DX5,I10_DX6,I10_DX7,I10_DX8,
           I10_DX9,I10_DX10,I10_DX11,I10_DX12,I10_DX13,I10_DX14,I10_DX15,I10_DX16,
           I10_DX17,I10_DX18,I10_DX19,I10_DX20,I10_DX21,I10_DX22,I10_DX23,I10_DX24,
           I10_DX25,I10_DX26,I10_DX27,I10_DX28,I10_DX29,I10_DX30,I10_DX31,I10_DX32,
           I10_DX33,I10_DX34,I10_DX35,I10_DX36,I10_DX37,I10_DX38,I10_DX39,I10_DX40,
           key="original",value="icdcm")
d3<-comorbidity(d2,id="INC_KEY",code="icdcm",map="charlson_icd10_quan",assign0=FALSE)
d4<-score(d3,weights="quan",assign0=FALSE)
d5<-bind_cols(d1,d4)
d5<-rename(d5,score=...50)


########## ADDED OR MODIFIED 2024 ###############################################



#Derive mechanism and intent for truncated codes, including basic ICD-10
d6<-read_csv("/Users/davideugeneclark/Documents/icdpicr/mech6.csv")

#Add a few common mechanisms not in CDC lists
others<-read.table(header=TRUE, quote="'", text="
   ICD10    MECHANISM                   INTENT
   V00131   'Pedestrian, other'       Unintentional
   V00141   'Pedestrian, other'       Unintentional
   V80010   'Other Land Transport'    Unintentional
   W01198    Fall                     Unintentional
   W01190    Fall                     Unintentional
   X959XX    Firearm                  Assault
   ")
d6a<-bind_rows(d6,others)

d7a<-mutate(d6a,digits12345=str_sub(ICD10,1,5))
d7a<-group_by(d7a,digits12345)
d7a<-mutate(d7a,MECH5=if_else(max(MECHANISM)==min(MECHANISM),max(MECHANISM),
                              "Other Specified"))
d7a<-mutate(d7a,INT5=if_else(max(INTENT)==min(INTENT),max(INTENT),
                             "Undetermined"))
d7a<-mutate(d7a,seq=row_number())
d7a<-ungroup(d7a)
d7b<-filter(d7a,seq==1)
d7b<-select(d7b,digits12345,MECH5,INT5)
d7b<-rename(d7b,ICD10=digits12345,MECHANISM=MECH5,INTENT=INT5)
write_csv(d7b,"mech5.csv")

d8a<-mutate(d6a,digits1234=str_sub(ICD10,1,4))
d8a<-group_by(d8a,digits1234)
d8a<-mutate(d8a,MECH4=if_else(max(MECHANISM)==min(MECHANISM),max(MECHANISM),
                              "Other Specified"))
d8a<-mutate(d8a,INT4=if_else(max(INTENT)==min(INTENT),max(INTENT),
                             "Undetermined"))
d8a<-mutate(d8a,seq=row_number())
d8a<-ungroup(d8a)
d8b<-filter(d8a,seq==1)
d8b<-select(d8b,digits1234,MECH4,INT4)
d8b<-rename(d8b,ICD10=digits1234,MECHANISM=MECH4,INTENT=INT4)
write_csv(d8b,"mech4.csv")

#Drop the following section
#d9a<-mutate(d6,digits123=str_sub(ICD10,1,3))
#d9a<-group_by(d9a,digits123)
#d9a<-mutate(d9a,MECH3=if_else(max(MECHANISM)==min(MECHANISM),max(MECHANISM),
#                             "Other Specified"))
#d9a<-mutate(d9a,INT3=if_else(max(INTENT)==min(INTENT),max(INTENT),
#                             "Undetermined"))
#d9a<-mutate(d9a,seq=row_number())
#d9a<-ungroup(d9a)
#d9b<-filter(d9a,seq==1)
#d9b<-select(d9b,digits123,MECH3,INT3)
#d9b<-rename(d9b,ICD10=digits123,MECHANISM=MECH3,INTENT=INT3)
#write_csv(d9b,"mech3.csv")

#Get all mechanism codes from NIS and TQP (7 digits)
#Assign categories based on truncated 6-digit equivalent, if any
dnis0 <- read_csv("/Users/davideugeneclark/Documents/icdpicr/nisraw2020.csv")
dnis1 <- select(dnis0,-AGE,-DIED,-DISPUNIFORM,-ELECTIVE,-HCUP_ED,-INJURY,-INC_KEY,-LOS,-seq)
dnis2 <- pivot_longer(dnis1,cols=starts_with("I10_DX"),names_to="ecode")
dnis3 <- mutate(dnis2,validmech=case_when(
              str_sub(value,1,1)=="T" & str_sub(value,2,1)=="5" ~ 1,
              str_sub(value,1,1)=="T" & str_sub(value,2,1)=="6" ~ 1,
              str_sub(value,1,1)=="T" & str_sub(value,2,1)=="7" ~ 1,
              str_sub(value,1,1)=="V" ~ 1,
              str_sub(value,1,1)=="W" ~ 1,
              str_sub(value,1,1)=="X" ~ 1,
              str_sub(value,1,1)=="Y" & str_sub(value,2,1)=="0" ~ 1,
              str_sub(value,1,1)=="Y" & str_sub(value,2,1)=="2" ~ 1,
              str_sub(value,1,1)=="Y" & str_sub(value,2,1)=="3" ~ 1,
              TRUE ~ 0
          ))
dnis4 <- filter(dnis3,validmech==1)
dnis4 <- rename(dnis4,ICD10=value)
dnis4 <- select(dnis4,ICD10)

dtqp0 <- read.csv("/Users/davideugeneclark/Documents/icdpicr/PUF AY 2020/CSV/PUF_TRAUMA.csv")  
dtqp1 <- select(dtqp0,PRIMARYECODEICD10)
dtqp2 <- mutate(dtqp1,predot=str_sub(PRIMARYECODEICD10,1,3))
dtqp2 <- mutate(dtqp2,postdot=str_sub(PRIMARYECODEICD10,5,8))
dtqp2 <- mutate(dtqp2,ICD10=str_c(predot,postdot))
dtqp3 <- select(dtqp2,ICD10)

dnistqp0 <- bind_rows(dnis4,dtqp3,d6a)
dnistqp0 <- group_by(dnistqp0,ICD10)
dnistqp0 <- arrange(dnistqp0,MECHANISM)
dnistqp0 <- mutate(dnistqp0,seq=row_number())
dnistqp0 <- ungroup(dnistqp0)
dnistqp0 <- arrange(dnistqp0,ICD10)

dnistqp1 <- filter(dnistqp0,seq==1)
dnistqp1 <- mutate(dnistqp1,digits123456=str_sub(ICD10,1,6))
dnistqp1 <- group_by(dnistqp1,digits123456)
dnistqp1 <- mutate(dnistqp1,MECHANISM=max(MECHANISM,na.rm=TRUE))
dnistqp1 <- mutate(dnistqp1,INTENT=max(INTENT,na.rm=TRUE))
dnistqp1 <- ungroup(dnistqp1)

dnistqp2 <- arrange(dnistqp1,ICD10)
dnistqp2 <- select(dnistqp2,-seq,-digits123456)

d10<-bind_rows(dnistqp2,d7b,d8b)
d10<-arrange(d10,ICD10)
write_csv(d10,"/Users/davideugeneclark/Documents/icdpicr2/ICD_Mech_241120.csv")


#  MAKE LOOKUP TABLES FOR ICDPICR2

etab<-read_csv("/Users/davideugeneclark/Documents/icdpicr2/ICD_Mech_241120.csv")
etab<-rename(etab,dx=ICD10,mechmaj=MECHANISM,intent=INTENT)
etab<-mutate(etab,mechmin="")
i10_map_mech<-distinct(etab)
i10_map_mech<-mutate(i10_map_mech,version="v241120")
write_csv(i10_map_mech,"/Users/davideugeneclark/Documents/icdpicr2/i10_map_mech_241120.csv")

i10_map_sev<-read_csv("/Users/davideugeneclark/Documents/icdpicr2/ICD_AIS_241119.csv")
i10_map_sev<-rename(i10_map_sev,dx=ICD,issbr=BR,severity=AIS)
i10_map_sev<-mutate(i10_map_sev,version="v241120")
write_csv(i10_map_sev,"/Users/davideugeneclark/Documents/icdpicr2/i10_map_sev_241120.csv")

i10_map_frame<-read_csv("/Users/davideugeneclark/Documents/icdpicr/ICD_Cells.csv")
i10_map_frame<-rename(i10_map_frame,dx=ICD)
i10_map_frame<-mutate(i10_map_frame,PsCell=1-PmCell)
i10_map_frame<-mutate(i10_map_frame,version="v241027")
write_csv(i10_map_frame,"/Users/davideugeneclark/Documents/icdpicr2/i10_map_frame_241027.csv")

