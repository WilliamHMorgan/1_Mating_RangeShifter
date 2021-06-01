#
#Mating system, habitat quality and the rate of range expansion 
#
#William H Morgan, Stephen C F Palmer & Xavier Lambin
#
#May 2021



######################################################################################################################################################
#
#Load packages
library(ggplot2)
library(plyr)
library(viridis)
library(ggpubr)
library(lemon)


#define functions
se <- function(x) sqrt(var(x, na.rm = T)/length(x))

#function to return change in occupancy between years of the simulation
occ.chng.func<- function(x) { 
  #diff, lag of 1, gives differences between entries of occupancy vector
  #use c(NA, diff...) to add NA to start of each vector
  c(NA, (diff(x, 1)))
}

#function to return change in range edge between years of the simulation
rng.chng.func<- function(x) { 
  #diff, lag of 1, gives differences between entries of rand extent vector
  #use c(NA, diff...) to add NA to start of each vector
  c(NA, (diff(x, 1)))
}

#define colour pallette
virid.pal<-viridis(9)


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#
#Create files for initilising occupied cells in the landscape with individuals
#
#
#First load the RangeShifter landscape file

setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_HIGHSURVIVAL/Inputs")


#land 0
ls0.07_0<- as.matrix(read.table("Series0discRandom_X20Y1000_p0.07_nr0.txt", skip = 6)) #use the skip argument to only read the matrix, not the meta data

#change row and column names to reflect coordinate system
row.names(ls0.07_0)<- c(999:0)
colnames(ls0.07_0)<- c(0:19)

#return the location of all the habitat cells - code = 2
hab0.07_0<- as.data.frame(which(ls0.07_0 == 2, arr.ind = T, useNames = T))
colnames(hab0.07_0)<- c("y","x") #remember - rows (y) then columns (x)
#now covnert row/column numbers to x and y coordinates
hab0.07_0$x_coord<- hab0.07_0$x - 1
hab0.07_0$y_coord<- 1000-hab0.07_0$y


#sort the data frame by y_coordinate
hab0.07_0<- hab0.07_0[,c("x_coord", "y_coord")]
hab0.07_0<- hab0.07_0[order(hab0.07_0$y_coord),]
hab0.07_0_half<- hab0.07_0
hab0.07_0_half<- hab0.07_0_half[hab0.07_0_half$y_coord < 500,]

#cut this in to slices 100 cells wide
hab0.07_0_half$slice100<- cut(hab0.07_0_half$y_coord, seq(0,500,100), include.lowest = T, labels = F, right =F)

#Select 50% of suitable habitat cells to initialise
#######################################
## For poisson with lamda 2
## 2 sex model
library(dplyr)

half.L2.50.init.cells<- hab0.07_0_half %>% group_by(slice100) %>% sample_frac(.5)

#repeat the data frame twice
half.2sex.L2.50.init.cells<- rbind(half.L2.50.init.cells, half.L2.50.init.cells)


#individuals per patch
#draw values for the number of adults on each patch from zero truncated poisson
library(extraDistr)

half.2sex.L2.50.init.cells$Ninds<- rtpois(nrow(half.2sex.L2.50.init.cells), lambda = 2, a = 0, b = Inf)

#add info on sex - remembering first half is for females (0) and second half is for males (1)
half.2sex.L2.50.init.cells$Sex<- c(rep(0, nrow(half.L2.50.init.cells)), rep(1, nrow(half.L2.50.init.cells)))

#turn in to RangeShifter compatible intitialisation file
half.2sex.L2.50.init.cells$Year<-0
half.2sex.L2.50.init.cells$Species<-0
half.2sex.L2.50.init.cells$Age<-1
half.2sex.L2.50.init.cells$Stage<-1

RngSh_half.2sex_.L2.50.1<- half.2sex.L2.50.init.cells[, which(names(half.2sex.L2.50.init.cells) %in% c("Year", "Species", "x_coord", "y_coord", "Ninds", "Sex", "Age", "Stage"))]
col_order<- c("Year", "Species", "x_coord", "y_coord", "Ninds", "Sex", "Age", "Stage")
RngSh_half.2sex_.L2.50.1<- RngSh_half.2sex_.L2.50.1[, col_order]
colnames(RngSh_half.2sex_.L2.50.1)<- c("Year", "Species", "X", "Y", "Ninds", "Sex", "Age", "Stage")

#"2" refers to the strength of density dependence (1/b = 2) and the lambda of the poisson = 2
write.table(RngSh_half.2sex_.L2.50.1, file = "init_2sex_half_50.L2.1.txt", sep = "\t",row.names = F, quote = F) 


#######################################
## For poisson with lamda 2
##single sex model

half.L2.50.init.cells<- hab0.07_0_half %>% group_by(slice100) %>% sample_frac(.5)


#individuals per patch
#draw values for the number of adults on each patch from zero truncated poisson
library(extraDistr)

half.L2.50.init.cells$fem.per.patch<- rtpois(nrow(half.L2.50.init.cells), lambda = 2, a = 0, b = Inf)


#turn in to RangeShifter compatible intitialisation file
half.L2.50.init.cells$Year<-0
half.L2.50.init.cells$Species<-0
half.L2.50.init.cells$Age<-1
half.L2.50.init.cells$Stage<-1

RngSh_half_.L2.50.1<- half.L2.50.init.cells[, which(names(half.L2.50.init.cells) %in% c("Year", "Species", "x_coord", "y_coord", "fem.per.patch", "Age", "Stage"))]
col_order<- c("Year", "Species", "x_coord", "y_coord", "fem.per.patch", "Age", "Stage")
RngSh_half_.L2.50.1<- RngSh_half_.L2.50.1[, col_order]
colnames(RngSh_half_.L2.50.1)<- c("Year", "Species", "X", "Y", "Ninds", "Age", "Stage")

#"2" refers to the strength of density dependence (1/b = 2) and the lambda of the poisson = 2
write.table(RngSh_half_.L2.50.1, file = "init_fem_half_50.L2.1.txt", sep = "\t",row.names = F, quote = F)




#The above files are used by RangeShifter
#
#
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################




######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#
#Experiment 1####
#
#Proportional occupancy and rate of expansion for unlimited, polygynous and monogamous populations for range of habitat qualities and fecundities
#
######################################################################################################################################################

setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_HIGHSURVIVAL/Outputs")

#Read in population files and add parameter info for habitat quality (dd) and fecundity (f)####


#Unlimited mating####


#Fecund = 2
Asex_2_f2_pop<- read.table("Batch2_Sim6_Land1_Pop.txt", header = T)
Asex_2_f2_pop$dd<-2
Asex_4_f2_pop<- read.table("Batch2_Sim7_Land1_Pop.txt", header = T)
Asex_4_f2_pop$dd<-4
Asex_6_f2_pop<- read.table("Batch2_Sim8_Land1_Pop.txt", header = T)
Asex_6_f2_pop$dd<-6
Asex_8_f2_pop<- read.table("Batch2_Sim9_Land1_Pop.txt", header = T)
Asex_8_f2_pop$dd<-8
Asex_10_f2_pop<- read.table("Batch2_Sim10_Land1_Pop.txt", header = T)
Asex_10_f2_pop$dd<-10
Asex_12_f2_pop<- read.table("Batch11_Sim3_Land1_Pop.txt", header = T)
Asex_12_f2_pop$dd<-12
Asex_14_f2_pop<- read.table("Batch11_Sim4_Land1_Pop.txt", header = T)
Asex_14_f2_pop$dd<-14

Asex_f2_pop<- rbind(Asex_2_f2_pop,Asex_4_f2_pop,Asex_6_f2_pop,Asex_8_f2_pop,Asex_10_f2_pop,Asex_12_f2_pop,Asex_14_f2_pop)

Asex_f2_pop$fecund<-2

#Fecund = 3
Asex_2_f3_pop<- read.table("Batch2_Sim11_Land1_Pop.txt", header = T)
Asex_2_f3_pop$dd<-2
Asex_4_f3_pop<- read.table("Batch2_Sim12_Land1_Pop.txt", header = T)
Asex_4_f3_pop$dd<-4
Asex_6_f3_pop<- read.table("Batch2_Sim13_Land1_Pop.txt", header = T)
Asex_6_f3_pop$dd<-6
Asex_8_f3_pop<- read.table("Batch2_Sim14_Land1_Pop.txt", header = T)
Asex_8_f3_pop$dd<-8
Asex_10_f3_pop<- read.table("Batch2_Sim15_Land1_Pop.txt", header = T)
Asex_10_f3_pop$dd<-10
Asex_12_f3_pop<- read.table("Batch11_Sim5_Land1_Pop.txt", header = T)
Asex_12_f3_pop$dd<-12
Asex_14_f3_pop<- read.table("Batch11_Sim6_Land1_Pop.txt", header = T)
Asex_14_f3_pop$dd<-14

Asex_f3_pop<- rbind(Asex_2_f3_pop,Asex_4_f3_pop,Asex_6_f3_pop,Asex_8_f3_pop,Asex_10_f3_pop,Asex_12_f3_pop,Asex_14_f3_pop)

Asex_f3_pop$fecund<-3

#Fecund = 4
Asex_2_f4_pop<- read.table("Batch2_Sim16_Land1_Pop.txt", header = T)
Asex_2_f4_pop$dd<-2
Asex_4_f4_pop<- read.table("Batch2_Sim17_Land1_Pop.txt", header = T)
Asex_4_f4_pop$dd<-4
Asex_6_f4_pop<- read.table("Batch2_Sim18_Land1_Pop.txt", header = T)
Asex_6_f4_pop$dd<-6
Asex_8_f4_pop<- read.table("Batch2_Sim19_Land1_Pop.txt", header = T)
Asex_8_f4_pop$dd<-8
Asex_10_f4_pop<- read.table("Batch2_Sim20_Land1_Pop.txt", header = T)
Asex_10_f4_pop$dd<-10
Asex_12_f4_pop<- read.table("Batch11_Sim7_Land1_Pop.txt", header = T)
Asex_12_f4_pop$dd<-12
Asex_14_f4_pop<- read.table("Batch11_Sim8_Land1_Pop.txt", header = T)
Asex_14_f4_pop$dd<-14

Asex_f4_pop<- rbind(Asex_2_f4_pop,Asex_4_f4_pop,Asex_6_f4_pop,Asex_8_f4_pop,Asex_10_f4_pop,Asex_12_f4_pop,Asex_14_f4_pop)

Asex_f4_pop$fecund<-4

Asex_ALL_POP<- rbind(Asex_f2_pop,Asex_f3_pop,Asex_f4_pop)



######################################################################################################################################################
######################################################################################################################################################
#Polygynous mating####

#Fecund = 2
SexLim_2_f2_pop<- read.table("Batch1_Sim6_Land1_Pop.txt", header = T)
SexLim_2_f2_pop$dd<-2
SexLim_4_f2_pop<- read.table("Batch1_Sim7_Land1_Pop.txt", header = T)
SexLim_4_f2_pop$dd<-4
SexLim_6_f2_pop<- read.table("Batch1_Sim8_Land1_Pop.txt", header = T)
SexLim_6_f2_pop$dd<-6
SexLim_8_f2_pop<- read.table("Batch1_Sim9_Land1_Pop.txt", header = T)
SexLim_8_f2_pop$dd<-8
SexLim_10_f2_pop<- read.table("Batch1_Sim10_Land1_Pop.txt", header = T)
SexLim_10_f2_pop$dd<-10
SexLim_12_f2_pop<- read.table("Batch10_Sim3_Land1_Pop.txt", header = T)
SexLim_12_f2_pop$dd<-12
SexLim_14_f2_pop<- read.table("Batch10_Sim4_Land1_Pop.txt", header = T)
SexLim_14_f2_pop$dd<-14

SexLim_f2_pop<- rbind(SexLim_2_f2_pop,SexLim_4_f2_pop,SexLim_6_f2_pop,SexLim_8_f2_pop,SexLim_10_f2_pop,SexLim_12_f2_pop,SexLim_14_f2_pop)

SexLim_f2_pop$fecund<-2

#Fecund = 3
SexLim_2_f3_pop<- read.table("Batch1_Sim11_Land1_Pop.txt", header = T)
SexLim_2_f3_pop$dd<-2
SexLim_4_f3_pop<- read.table("Batch1_Sim12_Land1_Pop.txt", header = T)
SexLim_4_f3_pop$dd<-4
SexLim_6_f3_pop<- read.table("Batch1_Sim13_Land1_Pop.txt", header = T)
SexLim_6_f3_pop$dd<-6
SexLim_8_f3_pop<- read.table("Batch1_Sim14_Land1_Pop.txt", header = T)
SexLim_8_f3_pop$dd<-8
SexLim_10_f3_pop<- read.table("Batch1_Sim15_Land1_Pop.txt", header = T)
SexLim_10_f3_pop$dd<-10
SexLim_12_f3_pop<- read.table("Batch10_Sim5_Land1_Pop.txt", header = T)
SexLim_12_f3_pop$dd<-12
SexLim_14_f3_pop<- read.table("Batch10_Sim6_Land1_Pop.txt", header = T)
SexLim_14_f3_pop$dd<-14

SexLim_f3_pop<- rbind(SexLim_2_f3_pop,SexLim_4_f3_pop,SexLim_6_f3_pop,SexLim_8_f3_pop,SexLim_10_f3_pop,SexLim_12_f3_pop,SexLim_14_f3_pop)

SexLim_f3_pop$fecund<-3

#Fecund = 4
SexLim_2_f4_pop<- read.table("Batch1_Sim16_Land1_Pop.txt", header = T)
SexLim_2_f4_pop$dd<-2
SexLim_4_f4_pop<- read.table("Batch1_Sim17_Land1_Pop.txt", header = T)
SexLim_4_f4_pop$dd<-4
SexLim_6_f4_pop<- read.table("Batch1_Sim18_Land1_Pop.txt", header = T)
SexLim_6_f4_pop$dd<-6
SexLim_8_f4_pop<- read.table("Batch1_Sim19_Land1_Pop.txt", header = T)
SexLim_8_f4_pop$dd<-8
SexLim_10_f4_pop<- read.table("Batch1_Sim20_Land1_Pop.txt", header = T)
SexLim_10_f4_pop$dd<-10
SexLim_12_f4_pop<- read.table("Batch10_Sim7_Land1_Pop.txt", header = T)
SexLim_12_f4_pop$dd<-12
SexLim_14_f4_pop<- read.table("Batch10_Sim8_Land1_Pop.txt", header = T)
SexLim_14_f4_pop$dd<-14

SexLim_f4_pop<- rbind(SexLim_2_f4_pop,SexLim_4_f4_pop,SexLim_6_f4_pop,SexLim_8_f4_pop,SexLim_10_f4_pop,SexLim_12_f4_pop,SexLim_14_f4_pop)

SexLim_f4_pop$fecund<-4

SexLim_ALL_POP<- rbind(SexLim_f2_pop,SexLim_f3_pop,SexLim_f4_pop)

######################################################################################################################################################
######################################################################################################################################################
#Monogamous mating####

#Fecund = 2
SexAllee_2_f2_pop<- read.table("Batch1_Sim46_Land1_Pop.txt", header = T)
SexAllee_2_f2_pop$dd<-2
SexAllee_4_f2_pop<- read.table("Batch1_Sim47_Land1_Pop.txt", header = T)
SexAllee_4_f2_pop$dd<-4
SexAllee_6_f2_pop<- read.table("Batch1_Sim48_Land1_Pop.txt", header = T)
SexAllee_6_f2_pop$dd<-6
SexAllee_8_f2_pop<- read.table("Batch1_Sim49_Land1_Pop.txt", header = T)
SexAllee_8_f2_pop$dd<-8
SexAllee_10_f2_pop<- read.table("Batch1_Sim50_Land1_Pop.txt", header = T)
SexAllee_10_f2_pop$dd<-10
SexAllee_12_f2_pop<- read.table("Batch10_Sim19_Land1_Pop.txt", header = T)
SexAllee_12_f2_pop$dd<-12
SexAllee_14_f2_pop<- read.table("Batch10_Sim20_Land1_Pop.txt", header = T)
SexAllee_14_f2_pop$dd<-14

SexAllee_f2_pop<- rbind(SexAllee_2_f2_pop,SexAllee_4_f2_pop,SexAllee_6_f2_pop,SexAllee_8_f2_pop,SexAllee_10_f2_pop,SexAllee_12_f2_pop,SexAllee_14_f2_pop)

SexAllee_f2_pop$fecund<-2

#Fecund = 3
SexAllee_2_f3_pop<- read.table("Batch1_Sim51_Land1_Pop.txt", header = T)
SexAllee_2_f3_pop$dd<-2
SexAllee_4_f3_pop<- read.table("Batch1_Sim52_Land1_Pop.txt", header = T)
SexAllee_4_f3_pop$dd<-4
SexAllee_6_f3_pop<- read.table("Batch1_Sim53_Land1_Pop.txt", header = T)
SexAllee_6_f3_pop$dd<-6
SexAllee_8_f3_pop<- read.table("Batch1_Sim54_Land1_Pop.txt", header = T)
SexAllee_8_f3_pop$dd<-8
SexAllee_10_f3_pop<- read.table("Batch1_Sim55_Land1_Pop.txt", header = T)
SexAllee_10_f3_pop$dd<-10
SexAllee_12_f3_pop<- read.table("Batch10_Sim21_Land1_Pop.txt", header = T)
SexAllee_12_f3_pop$dd<-12
SexAllee_14_f3_pop<- read.table("Batch10_Sim22_Land1_Pop.txt", header = T)
SexAllee_14_f3_pop$dd<-14

SexAllee_f3_pop<- rbind(SexAllee_2_f3_pop,SexAllee_4_f3_pop,SexAllee_6_f3_pop,SexAllee_8_f3_pop,SexAllee_10_f3_pop,SexAllee_12_f3_pop,SexAllee_14_f3_pop)

SexAllee_f3_pop$fecund<-3

#Fecund = 4
SexAllee_2_f4_pop<- read.table("Batch1_Sim56_Land1_Pop.txt", header = T)
SexAllee_2_f4_pop$dd<-2
SexAllee_4_f4_pop<- read.table("Batch1_Sim57_Land1_Pop.txt", header = T)
SexAllee_4_f4_pop$dd<-4
SexAllee_6_f4_pop<- read.table("Batch1_Sim58_Land1_Pop.txt", header = T)
SexAllee_6_f4_pop$dd<-6
SexAllee_8_f4_pop<- read.table("Batch1_Sim59_Land1_Pop.txt", header = T)
SexAllee_8_f4_pop$dd<-8
SexAllee_10_f4_pop<- read.table("Batch1_Sim60_Land1_Pop.txt", header = T)
SexAllee_10_f4_pop$dd<-10
SexAllee_12_f4_pop<- read.table("Batch10_Sim23_Land1_Pop.txt", header = T)
SexAllee_12_f4_pop$dd<-12
SexAllee_14_f4_pop<- read.table("Batch10_Sim24_Land1_Pop.txt", header = T)
SexAllee_14_f4_pop$dd<-14

SexAllee_f4_pop<- rbind(SexAllee_2_f4_pop,SexAllee_4_f4_pop,SexAllee_6_f4_pop,SexAllee_8_f4_pop,SexAllee_10_f4_pop,SexAllee_12_f4_pop,SexAllee_14_f4_pop)

SexAllee_f4_pop$fecund<-4

SexAllee_ALL_POP<- rbind(SexAllee_f2_pop,SexAllee_f3_pop,SexAllee_f4_pop)


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#Poportional occupancy in core####
#
#Core defined as 1st half of landscape (500 cells high)
#
#There are 664 cells in the 1st half of the landscape

init.Asex_ALL_POP<- Asex_ALL_POP[Asex_ALL_POP$y<500,]
init.Asex_ALL_POP$occ.breeding<- ifelse(init.Asex_ALL_POP$NInd_stage1 >0, 1, 0) 

#sum the number of occupied patches per year
init.Asex_ALL_POP.summary<- ddply(init.Asex_ALL_POP, c("Rep", "Year", "dd", "fecund"), summarise, total.occ.breeding=sum(occ.breeding))
init.Asex_ALL_POP.summary$prop.occ<- init.Asex_ALL_POP.summary$total.occ.breeding / 664 #664 is number of cells in intitialised zone
init.Asex_ALL_POP.summary$mating<- factor("Unlimited")
init.Asex_ALL_POP.summary$un.occ<- 664 - init.Asex_ALL_POP.summary$total.occ.breeding


init.SexLim_ALL_POP<- SexLim_ALL_POP[SexLim_ALL_POP$y<500,]
init.SexLim_ALL_POP$occ.breeding<- ifelse(init.SexLim_ALL_POP$Nfemales_stage1 >0 & init.SexLim_ALL_POP$Nmales_stage1 >0, 1, 0) 

init.SexLim_ALL_POP.summary<- ddply(init.SexLim_ALL_POP, c("Rep", "Year", "dd", "fecund"), summarise, total.occ.breeding=sum(occ.breeding))
init.SexLim_ALL_POP.summary$prop.occ<- init.SexLim_ALL_POP.summary$total.occ.breeding / 664 #664 is number of cells in intitialised zone
init.SexLim_ALL_POP.summary$mating<- factor("Polygyny")
init.SexLim_ALL_POP.summary$un.occ<- 664 - init.SexLim_ALL_POP.summary$total.occ.breeding


init.SexAllee_ALL_POP<- SexAllee_ALL_POP[SexAllee_ALL_POP$y<500,]
init.SexAllee_ALL_POP$occ.breeding<- ifelse(init.SexAllee_ALL_POP$Nfemales_stage1 >0 & init.SexAllee_ALL_POP$Nmales_stage1 >0, 1, 0) 

init.SexAllee_ALL_POP.summary<- ddply(init.SexAllee_ALL_POP, c("Rep", "Year", "dd", "fecund"), summarise, total.occ.breeding=sum(occ.breeding))
init.SexAllee_ALL_POP.summary$prop.occ<- init.SexAllee_ALL_POP.summary$total.occ.breeding / 664 #664 is number of cells in intitialised zone
init.SexAllee_ALL_POP.summary$mating<- factor("Monogamy")
init.SexAllee_ALL_POP.summary$un.occ<- 664 - init.SexAllee_ALL_POP.summary$total.occ.breeding

occupancy.inits.zone<- rbind(init.Asex_ALL_POP.summary, init.SexLim_ALL_POP.summary, init.SexAllee_ALL_POP.summary)


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#
#
#Equilibrium occupancy####
#
#Which populations could maintain equilibrium proportion occupancy in the core of the range?
#Only those simulations where proportional occupancy is not declining are included in further analyses
#

#First remove those that went extinct by year 50
table(occupancy.inits.zone[occupancy.inits.zone$Year==50,]$dd,occupancy.inits.zone[occupancy.inits.zone$Year==50,]$fecund,occupancy.inits.zone[occupancy.inits.zone$Year==50,]$mating)

occupancy.inits.zone.0<- occupancy.inits.zone[occupancy.inits.zone$fecund !=2 | occupancy.inits.zone$dd >12 | occupancy.inits.zone$mating !='Monogamy',]
occupancy.inits.zone.0<- occupancy.inits.zone.0[occupancy.inits.zone.0$fecund !=3 | occupancy.inits.zone.0$dd >6 | occupancy.inits.zone.0$mating !='Monogamy',]
occupancy.inits.zone.0<- occupancy.inits.zone.0[occupancy.inits.zone.0$fecund !=4 | occupancy.inits.zone.0$dd >4 | occupancy.inits.zone.0$mating !='Monogamy',]

occupancy.inits.zone.0<- occupancy.inits.zone.0[occupancy.inits.zone.0$fecund !=2 | occupancy.inits.zone.0$dd >8 |occupancy.inits.zone.0$mating !='Polygyny',]
occupancy.inits.zone.0<- occupancy.inits.zone.0[occupancy.inits.zone.0$fecund !=3 | occupancy.inits.zone.0$dd >4 |occupancy.inits.zone.0$mating !='Polygyny',]
occupancy.inits.zone.0<- occupancy.inits.zone.0[occupancy.inits.zone.0$fecund !=4 | occupancy.inits.zone.0$dd >4 |occupancy.inits.zone.0$mating !='Polygyny',]

occupancy.inits.zone.0<- occupancy.inits.zone.0[occupancy.inits.zone.0$fecund !=2 | occupancy.inits.zone.0$dd !=2 |occupancy.inits.zone.0$mating !='Unlimited',]

table(occupancy.inits.zone.0$dd,occupancy.inits.zone.0$fecund,occupancy.inits.zone.0$mating)


#derive change in occupany between 20 - 50 years

occupancy.inits.zone1<- ddply(occupancy.inits.zone.0, c("Rep","dd","fecund","mating"), mutate, occ.chng=occ.chng.func(prop.occ))
occupancy.inits.zone2<- na.omit(occupancy.inits.zone1)


delta.occ<- ddply(occupancy.inits.zone2[occupancy.inits.zone2$Year>=20,], c("dd","fecund","mating"),summarise, mean.delta=mean(occ.chng), se.delta=se(occ.chng))
delta.occ$mean.plus.error<- delta.occ$mean.delta+delta.occ$se.delta
delta.occ$persist<- ifelse(delta.occ$mean.plus.error<0,0,1)

# New facet label names for fecund variable
fecund.labs <- c("Mean fecundity = 2.0", "Mean fecundity = 3.0", "Mean fecundity = 4.0")
names(fecund.labs)<- c("2","3","4")
panel.labels<- data.frame(label=c("a","b","c"), fecund=c(2,3,4))

ggplot()+
  geom_point(data=delta.occ, aes(x=factor(dd), y=mean.delta, colour=mating), position=position_dodge(width=0.3))+
  geom_errorbar(data=delta.occ, aes(x=factor(dd), ymin=mean.delta-se.delta, ymax=mean.delta+se.delta, colour=mating), width=0.2, position=position_dodge(width=0.3))+
  scale_colour_manual(values = c(virid.pal[2], virid.pal[4], virid.pal[8]),labels=c("Not mate limited","Polygynous","Monogamous"))+
  facet_rep_wrap(~fecund,nrow=3,ncol=1,as.table = F, labeller = labeller(fecund = fecund.labs),strip.position='right')+
  geom_hline(yintercept =  0, linetype=2)+
  ylab("Change in occupancy")+
  xlab("Habitat quality 1/b")+
  theme_classic()

#save as high resolution jpeg
#ggsave("Delta_occupancy_quasi_equilibrium.jpg",dpi=300,device="jpeg",width=18,height=18, units="cm")



######################################################################################################################################################
######################################################################################################################################################
# Proportional occupancy after 50 years for simulations with quasi equilibrium occupancy####

#Add rows to the data frame for simulations that went extinct to give full factorial design i.e. they have occupancy of 0 after 50 years

#Unlmited
unlim.to.add<- data.frame(dd=rep(2,5),fecund=rep(2,5), mating=c("Unlimited","Unlimited","Unlimited","Unlimited","Unlimited"),prop.occ=rep(0,5))
#Polygyny
poly.to.add<- data.frame(dd=c(rep(seq(2,8,by=2),each=5),rep(seq(2,4,by=2),each=5),rep(seq(2,4,by=2),each=5)),
                         fecund=c(rep(2,20),rep(3,10),rep(4,10)), mating="Polygyny",prop.occ=rep(0,40))
#Monogamy
mono.to.add.1<- data.frame(dd=c(rep(seq(2,10,by=2),each=5),rep(seq(2,6,by=2),each=5),rep(seq(2,4,by=2),each=5)),
                           fecund=c(rep(2,25),rep(3,15),rep(4,10)), mating="Monogamy",prop.occ=rep(0,50))
mono.to.add.2<- data.frame(dd=c(rep(12,2)),
                           fecund=c(rep(2,2)), mating="Monogamy",prop.occ=rep(0,2))

mono.to.add<- rbind(mono.to.add.1,mono.to.add.2)

#Retrieve occupancy after 50 years from simulations that did not go extinct and bind with objects just created to give complete set of parameter combinations
occupancy.inits.zone50<- occupancy.inits.zone[occupancy.inits.zone$Year==50,]

occupancy.inits.zone50.sub<- occupancy.inits.zone50[,c("dd","fecund","mating","prop.occ")]


occupancy.inits.zone50.complete<- rbind(occupancy.inits.zone50.sub,unlim.to.add,poly.to.add,mono.to.add)

#Check we have a complete set of simulations for the statistical analysis
table(occupancy.inits.zone50.complete$dd,occupancy.inits.zone50.complete$fecund,occupancy.inits.zone50.complete$mating)


#Linear model of proportional occupancy####

occupancy.inits.zone50.complete$dd.fac<- factor(occupancy.inits.zone50.complete$dd)
occupancy.inits.zone50.complete$fecund.fac<- factor(occupancy.inits.zone50.complete$fecund)

prev1<- lm(prop.occ~dd.fac*fecund.fac*mating, data=occupancy.inits.zone50.complete)

summary(prev1)$adj.r.squared


########
# See variance explained by each main effect and interaction
prev1.fit<- anova(prev1)
prev1.fit.ssq<- prev1.fit$`Sum Sq`
var.expl<- prev1.fit.ssq/sum(prev1.fit.ssq)*100
round(var.expl,1)
cbind(row.names(prev1.fit),round(var.expl,1))

#return mean and standar error of proportional occupancy for each parameter combination
plot.occ<- ddply(occupancy.inits.zone50.complete, c("fecund", "dd", "mating"), summarise, mean.occ=mean(prop.occ), se.occ=se(prop.occ))


# New facet label names for fecund variable
fecund.labs <- c("Mean fecundity = 2.0", "Mean fecundity = 3.0", "Mean fecundity = 4.0")
names(fecund.labs)<- c("2","3","4")
panel.labels<- data.frame(label=c("a","b","c"), fecund=c(2,3,4))

#create plot ocject
occ.plot.object<-ggplot()+
  geom_point(data=plot.occ, aes(x=dd, y=mean.occ, colour=mating),size=0.85)+
  geom_line(data=plot.occ, aes(x=dd, y=mean.occ, colour=mating,group=mating))+
  geom_errorbar(data=plot.occ, aes(x=dd, ymin=mean.occ-se.occ, ymax=mean.occ+se.occ, colour=mating),width=0)+
  scale_colour_manual(values = c(virid.pal[2], virid.pal[4], virid.pal[8]),labels=c("Not mate limited","Polygynous","Monogamous"))+
  geom_text(data=panel.labels, aes(x=3,y=0.9,label=label))+
  ylab("Metapopulation size")+
  xlab("")+
  theme_classic()+
  scale_x_continuous(breaks=seq(2,14,by=2))+
  facet_rep_wrap(~fecund,nrow=1,ncol=3,as.table = F, labeller = labeller(fecund = fecund.labs))+
  theme( strip.background = element_blank(), axis.text.x = element_blank(),legend.position = 'top',legend.title = element_blank())

#save as high resolution jpeg
#ggsave("occupancyVsHabitat_facet_by_fecund.jpeg",dpi=300,device="jpeg",width=18,height=7, units="cm")


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#Rate of range expansion####
#
######################################################################################################################################################

setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_HIGHSURVIVAL/Outputs")

#Read in population files and add parameter info for habitat quality (dd) and fecundity (f)####
######################################################################################################################################################

#Unlimited mating####


#Fecund = 2
Asex_2_f2_Range<- read.table("Batch2_Sim6_Land1_Range.txt", header = T)
Asex_2_f2_Range$dd<-2
Asex_4_f2_Range<- read.table("Batch2_Sim7_Land1_Range.txt", header = T)
Asex_4_f2_Range$dd<-4
Asex_6_f2_Range<- read.table("Batch2_Sim8_Land1_Range.txt", header = T)
Asex_6_f2_Range$dd<-6
Asex_8_f2_Range<- read.table("Batch2_Sim9_Land1_Range.txt", header = T)
Asex_8_f2_Range$dd<-8
Asex_10_f2_Range<- read.table("Batch2_Sim10_Land1_Range.txt", header = T)
Asex_10_f2_Range$dd<-10
Asex_12_f2_Range<- read.table("Batch11_Sim3_Land1_Range.txt", header = T)
Asex_12_f2_Range$dd<-12
Asex_14_f2_Range<- read.table("Batch11_Sim4_Land1_Range.txt", header = T)
Asex_14_f2_Range$dd<-14

Asex_f2_Range<- rbind(Asex_2_f2_Range,Asex_4_f2_Range,Asex_6_f2_Range,Asex_8_f2_Range,Asex_10_f2_Range,Asex_12_f2_Range,Asex_14_f2_Range)

Asex_f2_Range$fecund<-2

#Fecund = 3
Asex_2_f3_Range<- read.table("Batch2_Sim11_Land1_Range.txt", header = T)
Asex_2_f3_Range$dd<-2
Asex_4_f3_Range<- read.table("Batch2_Sim12_Land1_Range.txt", header = T)
Asex_4_f3_Range$dd<-4
Asex_6_f3_Range<- read.table("Batch2_Sim13_Land1_Range.txt", header = T)
Asex_6_f3_Range$dd<-6
Asex_8_f3_Range<- read.table("Batch2_Sim14_Land1_Range.txt", header = T)
Asex_8_f3_Range$dd<-8
Asex_10_f3_Range<- read.table("Batch2_Sim15_Land1_Range.txt", header = T)
Asex_10_f3_Range$dd<-10
Asex_12_f3_Range<- read.table("Batch11_Sim5_Land1_Range.txt", header = T)
Asex_12_f3_Range$dd<-12
Asex_14_f3_Range<- read.table("Batch11_Sim6_Land1_Range.txt", header = T)
Asex_14_f3_Range$dd<-14

Asex_f3_Range<- rbind(Asex_2_f3_Range,Asex_4_f3_Range,Asex_6_f3_Range,Asex_8_f3_Range,Asex_10_f3_Range,Asex_12_f3_Range,Asex_14_f3_Range)

Asex_f3_Range$fecund<-3

#Fecund = 4
Asex_2_f4_Range<- read.table("Batch2_Sim16_Land1_Range.txt", header = T)
Asex_2_f4_Range$dd<-2
Asex_4_f4_Range<- read.table("Batch2_Sim17_Land1_Range.txt", header = T)
Asex_4_f4_Range$dd<-4
Asex_6_f4_Range<- read.table("Batch2_Sim18_Land1_Range.txt", header = T)
Asex_6_f4_Range$dd<-6
Asex_8_f4_Range<- read.table("Batch2_Sim19_Land1_Range.txt", header = T)
Asex_8_f4_Range$dd<-8
Asex_10_f4_Range<- read.table("Batch2_Sim20_Land1_Range.txt", header = T)
Asex_10_f4_Range$dd<-10
Asex_12_f4_Range<- read.table("Batch11_Sim7_Land1_Range.txt", header = T)
Asex_12_f4_Range$dd<-12
Asex_14_f4_Range<- read.table("Batch11_Sim8_Land1_Range.txt", header = T)
Asex_14_f4_Range$dd<-14

Asex_f4_Range<- rbind(Asex_2_f4_Range,Asex_4_f4_Range,Asex_6_f4_Range,Asex_8_f4_Range,Asex_10_f4_Range,Asex_12_f4_Range,Asex_14_f4_Range)

Asex_f4_Range$fecund<-4

Asex_ALL_Range<- rbind(Asex_f2_Range,Asex_f3_Range,Asex_f4_Range)
Asex_ALL_Range$mating<- "Unlimited"



######################################################################################################################################################
######################################################################################################################################################
#Polygynous mating####


#Fecund = 2
SexLim_2_f2_Range<- read.table("Batch1_Sim6_Land1_Range.txt", header = T)
SexLim_2_f2_Range$dd<-2
SexLim_4_f2_Range<- read.table("Batch1_Sim7_Land1_Range.txt", header = T)
SexLim_4_f2_Range$dd<-4
SexLim_6_f2_Range<- read.table("Batch1_Sim8_Land1_Range.txt", header = T)
SexLim_6_f2_Range$dd<-6
SexLim_8_f2_Range<- read.table("Batch1_Sim9_Land1_Range.txt", header = T)
SexLim_8_f2_Range$dd<-8
SexLim_10_f2_Range<- read.table("Batch1_Sim10_Land1_Range.txt", header = T)
SexLim_10_f2_Range$dd<-10
SexLim_12_f2_Range<- read.table("Batch10_Sim3_Land1_Range.txt", header = T)
SexLim_12_f2_Range$dd<-12
SexLim_14_f2_Range<- read.table("Batch10_Sim4_Land1_Range.txt", header = T)
SexLim_14_f2_Range$dd<-14

SexLim_f2_Range<- rbind(SexLim_2_f2_Range,SexLim_4_f2_Range,SexLim_6_f2_Range,SexLim_8_f2_Range,SexLim_10_f2_Range,SexLim_12_f2_Range,SexLim_14_f2_Range)

SexLim_f2_Range$fecund<-2

#Fecund = 3
SexLim_2_f3_Range<- read.table("Batch1_Sim11_Land1_Range.txt", header = T)
SexLim_2_f3_Range$dd<-2
SexLim_4_f3_Range<- read.table("Batch1_Sim12_Land1_Range.txt", header = T)
SexLim_4_f3_Range$dd<-4
SexLim_6_f3_Range<- read.table("Batch1_Sim13_Land1_Range.txt", header = T)
SexLim_6_f3_Range$dd<-6
SexLim_8_f3_Range<- read.table("Batch1_Sim14_Land1_Range.txt", header = T)
SexLim_8_f3_Range$dd<-8
SexLim_10_f3_Range<- read.table("Batch1_Sim15_Land1_Range.txt", header = T)
SexLim_10_f3_Range$dd<-10
SexLim_12_f3_Range<- read.table("Batch10_Sim5_Land1_Range.txt", header = T)
SexLim_12_f3_Range$dd<-12
SexLim_14_f3_Range<- read.table("Batch10_Sim6_Land1_Range.txt", header = T)
SexLim_14_f3_Range$dd<-14

SexLim_f3_Range<- rbind(SexLim_2_f3_Range,SexLim_4_f3_Range,SexLim_6_f3_Range,SexLim_8_f3_Range,SexLim_10_f3_Range,SexLim_12_f3_Range,SexLim_14_f3_Range)

SexLim_f3_Range$fecund<-3

#Fecund = 4
SexLim_2_f4_Range<- read.table("Batch1_Sim16_Land1_Range.txt", header = T)
SexLim_2_f4_Range$dd<-2
SexLim_4_f4_Range<- read.table("Batch1_Sim17_Land1_Range.txt", header = T)
SexLim_4_f4_Range$dd<-4
SexLim_6_f4_Range<- read.table("Batch1_Sim18_Land1_Range.txt", header = T)
SexLim_6_f4_Range$dd<-6
SexLim_8_f4_Range<- read.table("Batch1_Sim19_Land1_Range.txt", header = T)
SexLim_8_f4_Range$dd<-8
SexLim_10_f4_Range<- read.table("Batch1_Sim20_Land1_Range.txt", header = T)
SexLim_10_f4_Range$dd<-10
SexLim_12_f4_Range<- read.table("Batch10_Sim7_Land1_Range.txt", header = T)
SexLim_12_f4_Range$dd<-12
SexLim_14_f4_Range<- read.table("Batch10_Sim8_Land1_Range.txt", header = T)
SexLim_14_f4_Range$dd<-14

SexLim_f4_Range<- rbind(SexLim_2_f4_Range,SexLim_4_f4_Range,SexLim_6_f4_Range,SexLim_8_f4_Range,SexLim_10_f4_Range,SexLim_12_f4_Range,SexLim_14_f4_Range)

SexLim_f4_Range$fecund<-4

SexLim_ALL_Range<- rbind(SexLim_f2_Range,SexLim_f3_Range,SexLim_f4_Range)
SexLim_ALL_Range$mating<- "Polygyny"


######################################################################################################################################################
######################################################################################################################################################
#Monogamous mating####


#Fecund = 2
SexAllee_2_f2_Range<- read.table("Batch1_Sim46_Land1_Range.txt", header = T)
SexAllee_2_f2_Range$dd<-2
SexAllee_4_f2_Range<- read.table("Batch1_Sim47_Land1_Range.txt", header = T)
SexAllee_4_f2_Range$dd<-4
SexAllee_6_f2_Range<- read.table("Batch1_Sim48_Land1_Range.txt", header = T)
SexAllee_6_f2_Range$dd<-6
SexAllee_8_f2_Range<- read.table("Batch1_Sim49_Land1_Range.txt", header = T)
SexAllee_8_f2_Range$dd<-8
SexAllee_10_f2_Range<- read.table("Batch1_Sim50_Land1_Range.txt", header = T)
SexAllee_10_f2_Range$dd<-10
SexAllee_12_f2_Range<- read.table("Batch10_Sim19_Land1_Range.txt", header = T)
SexAllee_12_f2_Range$dd<-12
SexAllee_14_f2_Range<- read.table("Batch10_Sim20_Land1_Range.txt", header = T)
SexAllee_14_f2_Range$dd<-14

SexAllee_f2_Range<- rbind(SexAllee_2_f2_Range,SexAllee_4_f2_Range,SexAllee_6_f2_Range,SexAllee_8_f2_Range,SexAllee_10_f2_Range,SexAllee_12_f2_Range,SexAllee_14_f2_Range)

SexAllee_f2_Range$fecund<-2

#Fecund = 3
SexAllee_2_f3_Range<- read.table("Batch1_Sim51_Land1_Range.txt", header = T)
SexAllee_2_f3_Range$dd<-2
SexAllee_4_f3_Range<- read.table("Batch1_Sim52_Land1_Range.txt", header = T)
SexAllee_4_f3_Range$dd<-4
SexAllee_6_f3_Range<- read.table("Batch1_Sim53_Land1_Range.txt", header = T)
SexAllee_6_f3_Range$dd<-6
SexAllee_8_f3_Range<- read.table("Batch1_Sim54_Land1_Range.txt", header = T)
SexAllee_8_f3_Range$dd<-8
SexAllee_10_f3_Range<- read.table("Batch1_Sim55_Land1_Range.txt", header = T)
SexAllee_10_f3_Range$dd<-10
SexAllee_12_f3_Range<- read.table("Batch10_Sim21_Land1_Range.txt", header = T)
SexAllee_12_f3_Range$dd<-12
SexAllee_14_f3_Range<- read.table("Batch10_Sim22_Land1_Range.txt", header = T)
SexAllee_14_f3_Range$dd<-14

SexAllee_f3_Range<- rbind(SexAllee_2_f3_Range,SexAllee_4_f3_Range,SexAllee_6_f3_Range,SexAllee_8_f3_Range,SexAllee_10_f3_Range,SexAllee_12_f3_Range,SexAllee_14_f3_Range)

SexAllee_f3_Range$fecund<-3

#Fecund = 4
SexAllee_2_f4_Range<- read.table("Batch1_Sim56_Land1_Range.txt", header = T)
SexAllee_2_f4_Range$dd<-2
SexAllee_4_f4_Range<- read.table("Batch1_Sim57_Land1_Range.txt", header = T)
SexAllee_4_f4_Range$dd<-4
SexAllee_6_f4_Range<- read.table("Batch1_Sim58_Land1_Range.txt", header = T)
SexAllee_6_f4_Range$dd<-6
SexAllee_8_f4_Range<- read.table("Batch1_Sim59_Land1_Range.txt", header = T)
SexAllee_8_f4_Range$dd<-8
SexAllee_10_f4_Range<- read.table("Batch1_Sim60_Land1_Range.txt", header = T)
SexAllee_10_f4_Range$dd<-10
SexAllee_12_f4_Range<- read.table("Batch10_Sim23_Land1_Range.txt", header = T)
SexAllee_12_f4_Range$dd<-12
SexAllee_14_f4_Range<- read.table("Batch10_Sim24_Land1_Range.txt", header = T)
SexAllee_14_f4_Range$dd<-14

SexAllee_f4_Range<- rbind(SexAllee_2_f4_Range,SexAllee_4_f4_Range,SexAllee_6_f4_Range,SexAllee_8_f4_Range,SexAllee_10_f4_Range,SexAllee_12_f4_Range,SexAllee_14_f4_Range)

SexAllee_f4_Range$fecund<-4

SexAllee_ALL_Range<- rbind(SexAllee_f2_Range,SexAllee_f3_Range,SexAllee_f4_Range)
SexAllee_ALL_Range$mating<- "Monogamy"


######################################################################################################################################################
#
#
#Bind all together

All.range<- rbind(Asex_ALL_Range[,c("Rep","Year","max_Y","dd","fecund","mating")],
                  SexLim_ALL_Range[,c("Rep","Year","max_Y","dd","fecund","mating")],
                  SexAllee_ALL_Range[,c("Rep","Year","max_Y","dd","fecund","mating")])

All.range$edge<- All.range$max_Y/100


######################################################################################################################################################
######################################################################################################################################################
#Plot rate of expansion####
#
#
#Remove those simulations where quasi-equilibrium was not maintained


complete.range<- All.range[All.range$fecund !=2 | All.range$dd >12 | All.range$mating !='Monogamy',]
complete.range<- complete.range[complete.range$fecund !=3 | complete.range$dd >6 | complete.range$mating !='Monogamy',]
complete.range<- complete.range[complete.range$fecund !=4 | complete.range$dd >4 | complete.range$mating !='Monogamy',]

complete.range<- complete.range[complete.range$fecund !=2 | complete.range$dd >8 |complete.range$mating !='Polygyny',]
complete.range<- complete.range[complete.range$fecund !=3 | complete.range$dd >6 |complete.range$mating !='Polygyny',]
complete.range<- complete.range[complete.range$fecund !=4 | complete.range$dd >4 |complete.range$mating !='Polygyny',]

complete.range<- complete.range[complete.range$fecund !=2 | complete.range$dd !=2 |complete.range$mating !='Unlimited',]
complete.range<- complete.range[complete.range$fecund !=3 | complete.range$dd !=2 |complete.range$mating !='Unlimited',]

#Extract data for every fifth year
complete.range<- complete.range[order(complete.range$Year), ]
complete.range<- complete.range[complete.range$Year==0|
                                  complete.range$Year==5|
                                  complete.range$Year==10|
                                  complete.range$Year==15|
                                  complete.range$Year==20|
                                  complete.range$Year==25|
                                  complete.range$Year==30|
                                  complete.range$Year==35|
                                  complete.range$Year==40|
                                  complete.range$Year==45|
                                  complete.range$Year==50,]

#use range change function defined above to derive change in range edge then divide by 5
complete.range<- ddply(complete.range, c("Rep", "dd", "fecund", "mating"), mutate, rng.chng = rng.chng.func(edge)) 
complete.range$rows.per.year<- complete.range$rng.chng/5

#Calculate average and SE of rate of expansion
#By mating, habitat qauity and fecundity
mean.rate.mating<- ddply(complete.range, c("mating", "dd", "fecund"), summarise, mean_rate=mean(rows.per.year, na.rm = T), se_rate=se(rows.per.year))
#By mating system only
means.by.mating<- ddply(complete.range, "mating", summarise, overall.mean=mean(rows.per.year, na.rm = T), se_rate=se(rows.per.year))

panel.labels2<- data.frame(label=c("d","e","f"), fecund=c(2,3,4))

#Create plot object
rate.plot.object<- ggplot()+
  geom_point(data=mean.rate.mating, aes(x=dd, y=mean_rate, colour=mating),size=1)+#
  geom_line(data=mean.rate.mating, aes(x=dd, y=mean_rate, colour=mating, group=mating))+
  geom_errorbar(data=mean.rate.mating, aes(x=dd, ymin=mean_rate-se_rate,ymax=mean_rate+se_rate,colour=mating), width=0)+
  scale_colour_manual(values = c(virid.pal[8], virid.pal[4], virid.pal[2]),guide = guide_legend(override.aes = list(colour='white')))+
  geom_text(data=panel.labels2, aes(x=3,y=8,label=label))+
  theme_classic()+
  scale_x_continuous(breaks=seq(2,14,by=2))+
  ylab(expression(Rate~of~spread~(rows~generation^-1)))+
  xlab("Habitat quality 1/b (females/patch)")+
  facet_rep_wrap(~fecund,nrow=1,ncol=3,as.table = F, labeller = labeller(fecund = fecund.labs))+
  theme(strip.background = element_blank(),strip.text.x = element_blank(), 
        legend.text = element_blank(),legend.title = element_blank(), legend.position = 'top')


#Plot occupancy and Rate of expansion together

ggarrange(occ.plot.object,rate.plot.object,nrow=2,ncol = 1)

#ggsave("BIG_rate_occ_combo_plot.jpeg",dpi=400,device="jpeg",width=18,height=16, units="cm")


######################################################################################################################################################
######################################################################################################################################################
#Linear model of rate of expansion####
rate.models<-complete.range

rate.models$dd.fac<- factor(rate.models$dd)
rate.models$fecund.fac<- factor(rate.models$fecund)


rate1<- lm(rows.per.year~dd.fac*fecund.fac*mating,  data=rate.models)

summary(rate1)$adj.r.squared
anova(rate1)

########
# See variance explained by each main effect and interaction
rate1.fit<- anova(rate1)
rate1.fit.ssq<- rate1.fit$`Sum Sq`
var.expl.rate<- rate1.fit.ssq/sum(rate1.fit.ssq)*100
round(var.expl.rate,1)
cbind(row.names(rate1.fit),round(var.expl.rate,1))



######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#
#
#Female mating rate####


#Return the number of females on each patch that remained unmated

#Polygyny:

#Add column denoting if only females are on the patch i.e. females remain unmated
SexLim_ALL_POP$single.fem<- ifelse(SexLim_ALL_POP$Nmales_stage1 == 0 & SexLim_ALL_POP$Nfemales_stage1 > 0 , 1, 0)

#multiply number of adult females by "single.sex" column - those patches that are single sex will have nr females x 1, both sexes nr females x 0
SexLim_ALL_POP$unmated.fem<- SexLim_ALL_POP$single.fem * SexLim_ALL_POP$Nfemales_stage1

SexLim_ALL_POP$mating<- factor("Polygyny")

SexLim_ALL_POP.to.bind<- SexLim_ALL_POP[,c("Rep","Year","x","y","unmated.fem","Nfemales_stage1","mating","NInd","fecund")]

#Monogamy:

#Males mate with only 1 females so just subtract number of males from number of females
SexAllee_ALL_POP$unmated.fem<- SexAllee_ALL_POP$Nfemales_stage1 - SexAllee_ALL_POP$Nmales_stage1
#set any negative values to zero
SexAllee_ALL_POP$unmated.fem[SexAllee_ALL_POP$unmated.fem < 0]<- 0

SexAllee_ALL_POP$mating<- factor("Monogamy")

SexAllee_ALL_POP.to.bind<- SexAllee_ALL_POP[,c("Rep","Year","x","y","unmated.fem","Nfemales_stage1","mating","NInd","fecund")]


#bind polygyny and monogamy together
Mate.Lim_Vs_Allee_MateRate<- rbind(SexLim_ALL_POP.to.bind,SexAllee_ALL_POP.to.bind)


#Model of female mating rate - local patch scale
#
#testing effect of local population size and mating system on female mating rate while controlling for fecundity

Mate.Lim_Vs_Allee_mod1<- glm(cbind(unmated.fem, Nfemales_stage1) ~ mating * NInd + fecund , family = 'binomial', data = Mate.Lim_Vs_Allee_MateRate)


###plot the model predictions

#create "fake" data set with all parameter combinations of interest
Mate.Fake.Allee <- data.frame(expand.grid(NInd = seq(min(Mate.Lim_Vs_Allee_MateRate$NInd),max(Mate.Lim_Vs_Allee_MateRate$NInd),0.1),
                                          mating = c('Monogamy', 'Polygyny'),
                                          fecund=3))

#predict response and bind with fake data
Mate.p1.Fake.Allee<- predict(Mate.Lim_Vs_Allee_mod1, newdata = Mate.Fake.Allee, type = 'response') 
Mate.new.Fake.Allee<- cbind(Mate.Fake.Allee, Mate.p1.Fake.Allee)

#add "Not mate limited" to the data frame for plotting
nml<- data.frame(NInd=c(0,24),mating="Not mate limited", fecund=3,Mate.p1.Fake.Allee=rep(0,2))

Mate.new.Fake.Allee.plot<- rbind(Mate.new.Fake.Allee,nml)

#Plot object
ggplot()+
  geom_line(data = Mate.new.Fake.Allee.plot, aes(x = NInd, y = Mate.p1.Fake.Allee, colour = mating,
                                                 group = mating), size =1.5)+
  scale_colour_manual(values = c(virid.pal[8], virid.pal[4],virid.pal[2]))+
  ylab("Proportion unmated females")+
  xlab("Local pop size")+
  theme_classic()+
  theme(legend.position="top",legend.title = element_blank())

#save as high resolution jpeg
#ggsave("mate_rate.jpeg",dpi=300,device="jpeg",width=4,height=4, units="in")





######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#
#
#Experiment 2####
#
#Rate of expansion and female mating under different dispersal and settlement rules

#No density dependence in settlement rule
setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_HIGHSURVIVAL/Outputs")

######################################################################################################################################################
#Poly - NUll

#Create "empty" data frame
SexLim0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Polygyny'), 3), settle=rep(factor('DI'),3), finding=rep(factor('N'),3), dd=c(6,8,10))
#Read in the range files
SexLim_6_f6_0.005<- read.table("Batch1_Sim13_Land1_Range.txt", header=T)
SexLim_6_f6_0.005$dd<-6
SexLim_8_f6_0.005<- read.table("Batch1_Sim14_Land1_Range.txt", header=T)
SexLim_8_f6_0.005$dd<-8
SexLim_10_f6_0.005<- read.table("Batch1_Sim15_Land1_Range.txt", header=T)
SexLim_10_f6_0.005$dd<-10

#bind together
SexLim_all<- rbind(SexLim_6_f6_0.005,SexLim_8_f6_0.005,SexLim_10_f6_0.005)

#convert max_y from metres to rows
SexLim_all$edge<- SexLim_all$max_Y/100

#select every fifth year
SexLim_all_5yr<- SexLim_all[SexLim_all$Year==0|SexLim_all$Year==5|SexLim_all$Year==10|SexLim_all$Year==15|SexLim_all$Year==20|SexLim_all$Year==25|SexLim_all$Year==30|SexLim_all$Year==35|SexLim_all$Year==40|SexLim_all$Year==45|SexLim_all$Year==50,]
#calculate change in range using range change function defined above and divide by 5
SexLim_all_5yr<- ddply(SexLim_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
SexLim_all_5yr$rng.chng<- SexLim_all_5yr$rng.chng/5
#remove NA values
SexLim_all_5yr<- SexLim_all_5yr[!is.na(SexLim_all_5yr$rng.chng),]
#calculate mean and SE of rate of range expansion
SexLim_all_5yr_summary<- ddply(SexLim_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

#add these values to the "empty" data frame created above
SexLim0.005$mean.rate<- SexLim_all_5yr_summary$mean.rate
SexLim0.005$se.rate<- SexLim_all_5yr_summary$se.rate


#repeat this process for all mating and dispersal combinations...

######################################################################################################################################################
#Monog - Null

#Create empty data frame
SexAllee0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Monogamy'), 3), settle=rep(factor('DI'),3), finding=rep(factor('N'),3), dd=c(6,8,10))
#Read in the range files
SexAllee_6_f6_0.005<- read.table("Batch1_Sim53_Land1_Range.txt", header=T)
SexAllee_6_f6_0.005$dd<- 6
SexAllee_8_f6_0.005<- read.table("Batch1_Sim54_Land1_Range.txt", header=T)
SexAllee_8_f6_0.005$dd<- 8
SexAllee_10_f6_0.005<- read.table("Batch1_Sim55_Land1_Range.txt", header=T)
SexAllee_10_f6_0.005$dd<- 10

#bind together
SexAllee_all<- rbind(SexAllee_6_f6_0.005,SexAllee_8_f6_0.005,SexAllee_10_f6_0.005)

#convert max_y from metres to rows
SexAllee_all$edge<- SexAllee_all$max_Y/100

SexAllee_all_5yr<- SexAllee_all[SexAllee_all$Year==0|SexAllee_all$Year==5|SexAllee_all$Year==10|SexAllee_all$Year==15|SexAllee_all$Year==20|SexAllee_all$Year==25|SexAllee_all$Year==30|SexAllee_all$Year==35|SexAllee_all$Year==40|SexAllee_all$Year==45|SexAllee_all$Year==50,]
SexAllee_all_5yr<- ddply(SexAllee_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
SexAllee_all_5yr$rng.chng<- SexAllee_all_5yr$rng.chng/5
SexAllee_all_5yr<- SexAllee_all_5yr[!is.na(SexAllee_all_5yr$rng.chng),]
SexAllee_all_5yr_summary<- ddply(SexAllee_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

SexAllee0.005$mean.rate<- SexAllee_all_5yr_summary$mean.rate
SexAllee0.005$se.rate<- SexAllee_all_5yr_summary$se.rate

#WHEN F=6 AND 1/B=6 NO POPS LAST FOR 50 YEARS


######################################################################################################################################################
#Poly - Mate search

#Create empty data frame
LimMate0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Polygyny'), 3), settle=rep(factor('DI'),3), finding=rep(factor('Y'),3), dd=c(6,8,10))
#Read in the range files
LimMate_6_f6_0.005<- read.table("Batch1_Sim33_Land1_Range.txt", header=T)
LimMate_6_f6_0.005$dd<-6
LimMate_8_f6_0.005<- read.table("Batch1_Sim34_Land1_Range.txt", header=T)
LimMate_8_f6_0.005$dd<-8
LimMate_10_f6_0.005<- read.table("Batch1_Sim35_Land1_Range.txt", header=T)
LimMate_10_f6_0.005$dd<-10

#bind together
LimMate_all<- rbind(LimMate_6_f6_0.005,LimMate_8_f6_0.005,LimMate_10_f6_0.005)

#convert max_y from metres to rows
LimMate_all$edge<- LimMate_all$max_Y/100

LimMate_all_5yr<- LimMate_all[LimMate_all$Year==0|LimMate_all$Year==5|LimMate_all$Year==10|LimMate_all$Year==15|LimMate_all$Year==20|LimMate_all$Year==25|LimMate_all$Year==30|LimMate_all$Year==35|LimMate_all$Year==40|LimMate_all$Year==45|LimMate_all$Year==50,]
LimMate_all_5yr<- ddply(LimMate_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
LimMate_all_5yr$rng.chng<- LimMate_all_5yr$rng.chng/5
LimMate_all_5yr<- LimMate_all_5yr[!is.na(LimMate_all_5yr$rng.chng),]
LimMate_all_5yr_summary<- ddply(LimMate_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

LimMate0.005$mean.rate<- LimMate_all_5yr_summary$mean.rate
LimMate0.005$se.rate<- LimMate_all_5yr_summary$se.rate


######################################################################################################################################################
#Mono - Mate search

#Create empty data frame
AlleeMate0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Monogamy'), 3), settle=rep(factor('DI'),3), finding=rep(factor('Y'),3), dd=rep(c(6,8,10)))
#Read in the range files
AlleeMate_6_f6_0.005<- read.table("Batch1_Sim73_Land1_Range.txt", header=T)
AlleeMate_6_f6_0.005$dd<-6
AlleeMate_8_f6_0.005<- read.table("Batch1_Sim74_Land1_Range.txt", header=T)
AlleeMate_8_f6_0.005$dd<-8
AlleeMate_10_f6_0.005<- read.table("Batch1_Sim75_Land1_Range.txt", header=T)
AlleeMate_10_f6_0.005$dd<-10

#bind together
AlleeMate_all<- rbind(AlleeMate_6_f6_0.005,AlleeMate_8_f6_0.005,AlleeMate_10_f6_0.005)

#convert max_y from metres to rows
AlleeMate_all$edge<- AlleeMate_all$max_Y/100

AlleeMate_all_5yr<- AlleeMate_all[AlleeMate_all$Year==0|AlleeMate_all$Year==5|AlleeMate_all$Year==10|AlleeMate_all$Year==15|AlleeMate_all$Year==20|AlleeMate_all$Year==25|AlleeMate_all$Year==30|AlleeMate_all$Year==35|AlleeMate_all$Year==40|AlleeMate_all$Year==45|AlleeMate_all$Year==50,]
AlleeMate_all_5yr<- ddply(AlleeMate_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
AlleeMate_all_5yr$rng.chng<- AlleeMate_all_5yr$rng.chng/5
AlleeMate_all_5yr<- AlleeMate_all_5yr[!is.na(AlleeMate_all_5yr$rng.chng),]
AlleeMate_all_5yr_summary<- ddply(AlleeMate_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

AlleeMate0.005$mean.rate<- AlleeMate_all_5yr_summary$mean.rate
AlleeMate0.005$se.rate<- AlleeMate_all_5yr_summary$se.rate


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

#Dendity dependent settlement rule

#****IMPORTANT NOTE***
#RangeShifter files have the same names for density dependent settlement and density independent settlement
#Ensure they are saved in a different folder such that the correct fil is called

setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_DDS_HIGHSURVIVAL/Outputs")

######################################################################################################################################################
#Poly - no mate finding

#Create empty data frame
DDS_SexLim0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Polygyny'), 3), settle=rep(factor('DD'),3), finding=rep(factor('N'),3), dd=c(6,8,10))
#Read in the range files
DDS_SexLim_6_f6_0.005<- read.table("Batch1_Sim13_Land1_Range.txt", header=T)
DDS_SexLim_6_f6_0.005$dd<-6
DDS_SexLim_8_f6_0.005<- read.table("Batch1_Sim14_Land1_Range.txt", header=T)
DDS_SexLim_8_f6_0.005$dd<-8
DDS_SexLim_10_f6_0.005<- read.table("Batch1_Sim15_Land1_Range.txt", header=T)
DDS_SexLim_10_f6_0.005$dd<-10

#bind together
DDS_SexLim_all<- rbind(DDS_SexLim_6_f6_0.005,DDS_SexLim_8_f6_0.005,DDS_SexLim_10_f6_0.005)

#convert max_y from metres to rows
DDS_SexLim_all$edge<- DDS_SexLim_all$max_Y/100

DDS_SexLim_all_5yr<- DDS_SexLim_all[DDS_SexLim_all$Year==0|DDS_SexLim_all$Year==5|DDS_SexLim_all$Year==10|DDS_SexLim_all$Year==15|DDS_SexLim_all$Year==20|DDS_SexLim_all$Year==25|DDS_SexLim_all$Year==30|DDS_SexLim_all$Year==35|DDS_SexLim_all$Year==40|DDS_SexLim_all$Year==45|DDS_SexLim_all$Year==50,]
DDS_SexLim_all_5yr<- ddply(DDS_SexLim_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
DDS_SexLim_all_5yr$rng.chng<- DDS_SexLim_all_5yr$rng.chng/5
DDS_SexLim_all_5yr<- DDS_SexLim_all_5yr[!is.na(DDS_SexLim_all_5yr$rng.chng),]
DDS_SexLim_all_5yr_summary<- ddply(DDS_SexLim_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

DDS_SexLim0.005$mean.rate<- DDS_SexLim_all_5yr_summary$mean.rate
DDS_SexLim0.005$se.rate<- DDS_SexLim_all_5yr_summary$se.rate

######################################################################################################################################################
#Mono - no mate finding

#Create empty data frame
DDS_SexAllee0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Monogamy'), 3), settle=rep(factor('DD'),3), finding=rep(factor('N'),3), dd=rep(c(6,8,10)))
#Read in the range files
DDS_SexAllee_6_f6_0.005<- read.table("Batch1_Sim53_Land1_Range.txt", header=T)
DDS_SexAllee_6_f6_0.005$dd<-6
DDS_SexAllee_8_f6_0.005<- read.table("Batch1_Sim54_Land1_Range.txt", header=T)
DDS_SexAllee_8_f6_0.005$dd<-8
DDS_SexAllee_10_f6_0.005<- read.table("Batch1_Sim55_Land1_Range.txt", header=T)
DDS_SexAllee_10_f6_0.005$dd<-10

#bind together
DDS_SexAllee_all<- rbind(DDS_SexAllee_6_f6_0.005,DDS_SexAllee_8_f6_0.005,DDS_SexAllee_10_f6_0.005)

#convert max_y from metres to rows
DDS_SexAllee_all$edge<- DDS_SexAllee_all$max_Y/100

DDS_SexAllee_all_5yr<- DDS_SexAllee_all[DDS_SexAllee_all$Year==0|DDS_SexAllee_all$Year==5|DDS_SexAllee_all$Year==10|DDS_SexAllee_all$Year==15|DDS_SexAllee_all$Year==20|DDS_SexAllee_all$Year==25|DDS_SexAllee_all$Year==30|DDS_SexAllee_all$Year==35|DDS_SexAllee_all$Year==40|DDS_SexAllee_all$Year==45|DDS_SexAllee_all$Year==50,]
DDS_SexAllee_all_5yr<- ddply(DDS_SexAllee_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
DDS_SexAllee_all_5yr$rng.chng<- DDS_SexAllee_all_5yr$rng.chng/5
DDS_SexAllee_all_5yr<- DDS_SexAllee_all_5yr[!is.na(DDS_SexAllee_all_5yr$rng.chng),]
DDS_SexAllee_all_5yr_summary<- ddply(DDS_SexAllee_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

DDS_SexAllee0.005$mean.rate<- DDS_SexAllee_all_5yr_summary$mean.rate
DDS_SexAllee0.005$se.rate<- DDS_SexAllee_all_5yr_summary$se.rate


######################################################################################################################################################
#Poly - mate finding


#Create empty data frame
DDS_LimMate0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Polygyny'), 3), settle=rep(factor('DD'),3), finding=rep(factor('Y'),3), dd=c(6,8,10))
#Read in the range files
DDS_LimMate_6_f6_0.005<- read.table("Batch1_Sim33_Land1_Range.txt", header=T)
DDS_LimMate_6_f6_0.005$dd<-6
DDS_LimMate_8_f6_0.005<- read.table("Batch1_Sim34_Land1_Range.txt", header=T)
DDS_LimMate_8_f6_0.005$dd<-8
DDS_LimMate_10_f6_0.005<- read.table("Batch1_Sim35_Land1_Range.txt", header=T)
DDS_LimMate_10_f6_0.005$dd<-10

#bind together
DDS_LimMate_all<- rbind(DDS_LimMate_6_f6_0.005,DDS_LimMate_8_f6_0.005,DDS_LimMate_10_f6_0.005)

#convert max_y from metres to rows
DDS_LimMate_all$edge<- DDS_LimMate_all$max_Y/100

DDS_LimMate_all_5yr<- DDS_LimMate_all[DDS_LimMate_all$Year==0|DDS_LimMate_all$Year==5|DDS_LimMate_all$Year==10|DDS_LimMate_all$Year==15|DDS_LimMate_all$Year==20|DDS_LimMate_all$Year==25|DDS_LimMate_all$Year==30|DDS_LimMate_all$Year==35|DDS_LimMate_all$Year==40|DDS_LimMate_all$Year==45|DDS_LimMate_all$Year==50,]
DDS_LimMate_all_5yr<- ddply(DDS_LimMate_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
DDS_LimMate_all_5yr$rng.chng<- DDS_LimMate_all_5yr$rng.chng/5
DDS_LimMate_all_5yr<- DDS_LimMate_all_5yr[!is.na(DDS_LimMate_all_5yr$rng.chng),]
DDS_LimMate_all_5yr_summary<- ddply(DDS_LimMate_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

DDS_LimMate0.005$mean.rate<- DDS_LimMate_all_5yr_summary$mean.rate
DDS_LimMate0.005$se.rate<- DDS_LimMate_all_5yr_summary$se.rate


######################################################################################################################################################
#Mono - mate finding

#Create empty data frame
DDS_AlleeMate0.005<- data.frame(PSM=rep(0.005,3),mating=rep(factor('Monogamy'), 3), settle=rep(factor('DD'),3), finding=rep(factor('Y'),3), dd=rep(c(6,8,10)))
#Read in the range files
DDS_AlleeMate_6_f6_0.005<- read.table("Batch1_Sim73_Land1_Range.txt", header=T)
DDS_AlleeMate_6_f6_0.005$dd<-6
DDS_AlleeMate_8_f6_0.005<- read.table("Batch1_Sim74_Land1_Range.txt", header=T)
DDS_AlleeMate_8_f6_0.005$dd<-8
DDS_AlleeMate_10_f6_0.005<- read.table("Batch1_Sim75_Land1_Range.txt", header=T)
DDS_AlleeMate_10_f6_0.005$dd<-10

#bind together
DDS_AlleeMate_all<- rbind(DDS_AlleeMate_6_f6_0.005,DDS_AlleeMate_8_f6_0.005,DDS_AlleeMate_10_f6_0.005)

#convert max_y from metres to rows
DDS_AlleeMate_all$edge<- DDS_AlleeMate_all$max_Y/100

DDS_AlleeMate_all_5yr<- DDS_AlleeMate_all[DDS_AlleeMate_all$Year==0|DDS_AlleeMate_all$Year==5|DDS_AlleeMate_all$Year==10|DDS_AlleeMate_all$Year==15|DDS_AlleeMate_all$Year==20|DDS_AlleeMate_all$Year==25|DDS_AlleeMate_all$Year==30|DDS_AlleeMate_all$Year==35|DDS_AlleeMate_all$Year==40|DDS_AlleeMate_all$Year==45|DDS_AlleeMate_all$Year==50,]
DDS_AlleeMate_all_5yr<- ddply(DDS_AlleeMate_all_5yr, c("Rep", "dd"), mutate, rng.chng=rng.chng.func(edge))
DDS_AlleeMate_all_5yr$rng.chng<- DDS_AlleeMate_all_5yr$rng.chng/5
DDS_AlleeMate_all_5yr<- DDS_AlleeMate_all_5yr[!is.na(DDS_AlleeMate_all_5yr$rng.chng),]
DDS_AlleeMate_all_5yr_summary<- ddply(DDS_AlleeMate_all_5yr, "dd", summarise, mean.rate=mean(rng.chng), se.rate=se(rng.chng))

DDS_AlleeMate0.005$mean.rate<- DDS_AlleeMate_all_5yr_summary$mean.rate
DDS_AlleeMate0.005$se.rate<- DDS_AlleeMate_all_5yr_summary$se.rate


######################################################################################################################################################
######################################################################################################################################################
#
#Bind all together

range.ext.0.005<- rbind(SexLim0.005, LimMate0.005, SexAllee0.005, AlleeMate0.005, 
                        DDS_SexLim0.005, DDS_LimMate0.005, DDS_SexAllee0.005, DDS_AlleeMate0.005)


#grab the values from the one sex simulations for comparision
comp_unlimited<- mean.rate.mating[mean.rate.mating$mating=="Unlimited"&mean.rate.mating$dd<12&mean.rate.mating$dd>4&mean.rate.mating$fecund==3,
                                  c("dd","fecund","mean_rate","se_rate")]

#colour palettes for figure
virid.pal2<-viridis(option='plasma',end=0.8,9)

#dodge points so error bars do not overlap
dodge1<-c(rep(c(-.02,.02),c(2,3)),rep(c(-.02,.02),c(2,2)),rep(c(-.02,.02),c(3,3)),rep(c(-.02,.02),c(2,3)))
M.S.labs<- c("No mate-search\nDI settlement","Mate-search\nDI settlement",
             "No mate-search\nDD settlement","Mate-search\nDD settlement")

panel.labs.ab<- data.frame(label=c("a","b"),mating=factor(c("Polygyny","Monogamy")))
panel.labs.ab$mating <- ordered(panel.labs.ab$mating , levels = c("Polygyny", "Monogamy"))

#plot object
ggplot()+
  geom_point(data=range.ext.0.005[range.ext.0.005$mean.rate>0,], aes(x=(dd+dodge1), y=mean.rate, colour=interaction(finding,settle)), size=1.5)+
  geom_line(data=range.ext.0.005[range.ext.0.005$mean.rate>0,], aes(x=(dd+dodge1), y=mean.rate, group=interaction(finding,settle), colour=interaction(finding,settle)), size=.5)+
  geom_errorbar(data=range.ext.0.005[range.ext.0.005$mean.rate>0,], aes(x=(dd+dodge1), ymin=mean.rate-se.rate,ymax=mean.rate+se.rate, group=interaction(finding,settle), colour=interaction(finding,settle)), width=0)+
  scale_colour_manual(values=c(virid.pal2[1],virid.pal2[6],virid.pal2[3],virid.pal2[9]), labels=M.S.labs)+
  theme_classic()+
  ylab(expression(Rate~of~spread~(rows~generation^-1)))+
  #xlim(0.3,1)+
  xlab("Habitat quality 1/b (females/patch)")+
  facet_rep_grid(.~mating)+
  geom_point(data=comp_unlimited, aes(x=dd, y=mean_rate), size=1.5)+
  geom_line(data=comp_unlimited, aes(x=dd, y=mean_rate), linetype=2)+
  geom_errorbar(data=comp_unlimited, aes(x=dd, ymin=mean_rate-se_rate,ymax=mean_rate+se_rate), width=0)+
  geom_text(data=panel.labs.ab, aes(x=6.5,y=8,label=label))+
  theme( strip.background = element_blank(),legend.position="top",legend.title = element_blank())


#save as high resolution jpeg
#ggsave("Rate_informed_dispersal_VS_unlimited.jpeg",dpi=400,device="jpeg",width=15,height=9, units="cm")





######################################################################################################################################################
######################################################################################################################################################
#
#Visulasing simulations - the effect of mate finding####

setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_DDS_HIGHSURVIVAL/Outputs")

#For dd = 10 and fecund = 3
#
#show breeding vs non-breeding cells

#no mate searching
DDS_SexLim_10_f3_pop<- read.table("Batch1_Sim15_Land1_Pop.txt", header = T)
DDS_SexLim_10_f3_pop$finding<- 0
#mate searching
DDS_LimMate_10_f6_pop<- read.table("Batch1_Sim35_Land1_Pop.txt", header = T)
DDS_LimMate_10_f6_pop$finding<- 1

#bind together
vis1<- rbind(DDS_SexLim_10_f3_pop,DDS_LimMate_10_f6_pop)

vis2<- vis1[vis1$NInd>0,]# & vis1$Rep==0

#Add column denoting if it was a breeding patch (males and females) or not (females only)

vis2$Breeding<- ifelse(vis2$Nmales_stage1 == 0 & vis2$Nfemales_stage1 > 0 , 'N', 'Y')



#Create two plot objects and combine with ggarrange()

#First plot object: No mate search
vis.finding.0<- vis2[vis2$finding==0,]


vis.plot.no.finding<- ggplot()+
  geom_jitter(data = vis.finding.0[vis.finding.0$y>400 & vis.finding.0$Year==10|
                                     vis.finding.0$y>400 & vis.finding.0$Year==20|
                                     vis.finding.0$y>400 & vis.finding.0$Year==30|
                                     vis.finding.0$y>400 & vis.finding.0$Year==40|
                                     vis.finding.0$y>400 & vis.finding.0$Year==50,], aes(x=x, y=y, colour=Breeding), height = 3, width = 1, size=0.2)+
  scale_colour_manual(values=c("white","black"),guide=F)+
  facet_grid(~ Year, as.table = F)+
  theme_void()+
  ylim(400,1000)+
  labs(title = "No mate-search")+
  theme(panel.background=element_rect(fill='grey50'), strip.text.x = element_blank(),plot.margin=margin(0.5,0.5,0.5,0.5,unit = 'cm'),
        title = element_text(size=10),panel.spacing = unit(0.1, "lines"))

#Second plot object: Mate search
vis.finding.1<- vis2[vis2$finding==1,]


vis.plot.mate.finding<- ggplot()+
  geom_jitter(data = vis.finding.1[vis.finding.1$y>400 & vis.finding.1$Year==10|
                                     vis.finding.1$y>400 & vis.finding.1$Year==20|
                                     vis.finding.1$y>400 & vis.finding.1$Year==30|
                                     vis.finding.1$y>400 & vis.finding.1$Year==40|
                                     vis.finding.1$y>400 & vis.finding.1$Year==50,], aes(x=x, y=y, colour=Breeding), height = 3, width = 1, size=0.2)+
  scale_colour_manual(values=c("white","black"),guide=F)+
  facet_grid(~ Year, as.table = F)+
  theme_void()+
  ylim(400,1000)+
  labs(title = "Mate-search")+
  theme(panel.background=element_rect(fill='grey50'), strip.text.x = element_blank(),plot.margin=margin(0.5,0.5,0.5,0.5,unit = 'cm'),
        title = element_text(size=10),panel.spacing = unit(0.1, "lines"))

#Plot
ggarrange(vis.plot.no.finding,vis.plot.mate.finding,nrow = 1)

#save as high resolution jpeg
#ggsave("visualsing simulation.jpeg",dpi=300,device="jpeg",width=11,height=11, units="cm")




######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
#
#
#Mating rate throughout the range####


#No density dependence in settlement rule
#****IMPORTANT NOTE***
#RangeShifter files have the same names for density dependent settlement and density independent settlement
#Ensure they are saved in a different folder such that the correct fil is called

setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_HIGHSURVIVAL/Outputs")

######################################################################################################################################################
#Polygyny - no mate finding
SexLim_10_f6_POP<- read.table("Batch1_Sim15_Land1_Pop.txt", header=T)
#Mate-search: (confusingly i called this "finding")
SexLim_10_f6_POP$finding<- 'N'
#Settlement is density independent
SexLim_10_f6_POP$settle<- 'DI'
SexLim_10_f6_POP$single.fem<- ifelse(SexLim_10_f6_POP$Nmales_stage1 == 0 & SexLim_10_f6_POP$Nfemales_stage1 > 0 , 1, 0)
SexLim_10_f6_POP$unmated.fem<- SexLim_10_f6_POP$single.fem * SexLim_10_f6_POP$Nfemales_stage1
#Whether a patch is occupied by a breeding populations
SexLim_10_f6_POP$occ<- ifelse(SexLim_10_f6_POP$Nfemales_stage1>0 & SexLim_10_f6_POP$Nmales_stage1>0, 1, 0)
#whether a patch is occupied by any population
SexLim_10_f6_POP$occ.any<- ifelse(SexLim_10_f6_POP$NInd>0,1,0)

####merge with range info
SexLim_ALL_POP_RNG<- merge(SexLim_10_f6_POP, SexLim_10_f6_0.005, all = T)
SexLim_ALL_POP_RNG$edge<- SexLim_ALL_POP_RNG$max_Y / 100

####create landscape slices - cut occupied range in to slices that are 10 rows wide

#First define the location of occupied cells relative to the range edge
SexLim_ALL_POP_RNG$rel.edge<- SexLim_ALL_POP_RNG$edge - SexLim_ALL_POP_RNG$y

#Then cut in to slices
SexLim_ALL_POP_RNG$edge.10<- cut(SexLim_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)
#
#Note - slice 101 is the range edge (as defined as further forward breeding patch), including any patch that is 0-9 rows behind the edge
# Slices 100 and less are those in front of the range edge - eight no longer occupied or occupied by a non-breeding population

#female mating rate at the patch level
SexLim_ALL_POP_RNG$mate.rate<- (SexLim_ALL_POP_RNG$Nfemales_stage1 - SexLim_ALL_POP_RNG$unmated.fem) / SexLim_ALL_POP_RNG$Nfemales_stage1

#Derive the number of suitable patches (i.e. Cells) per landscape slice

#1. turn "edge.10" in to the number of rows behind the range edge the FURTHEST cell could be and still belong to that slice
SexLim_ALL_POP_RNG$furthest.y<- (SexLim_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range edge the CLOSEST cell could be and still belong to that slice
SexLim_ALL_POP_RNG$closest.y<- (SexLim_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0


#Now count the number of suitable patches in each of the landscape slices each year of the simulation:

cellsINslice<- rep(NA, nrow(SexLim_ALL_POP_RNG))

for(i in 1:nrow(SexLim_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (SexLim_ALL_POP_RNG$edge[i] - (SexLim_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (SexLim_ALL_POP_RNG$edge[i] - SexLim_ALL_POP_RNG$closest.y[i])] )
}

SexLim_ALL_POP_RNG$cellsINslice<- cellsINslice

#This can then be used later to calculate the mating rate and proportional occupancy in each "neighbourhood" (i.e. landscape slice)



#Repeat the same process for all mating system and dispersal strategy combinations
######################################################################################################################################################
#Monogamy - no mate finding
SexAllee_10_f6_POP<- read.table("Batch1_Sim55_Land1_Pop.txt", header=T)
SexAllee_10_f6_POP$finding<- 'N'
SexAllee_10_f6_POP$settle<- 'DI'
SexAllee_10_f6_POP$unmated.fem<- SexAllee_10_f6_POP$Nfemales_stage1 - SexAllee_10_f6_POP$Nmales_stage1
SexAllee_10_f6_POP$unmated.fem[SexAllee_10_f6_POP$unmated.fem < 0]<- 0
SexAllee_10_f6_POP$occ<- ifelse(SexAllee_10_f6_POP$Nfemales_stage1>0 & SexAllee_10_f6_POP$Nmales_stage1>0, 1, 0)
SexAllee_10_f6_POP$occ.any<- ifelse(SexAllee_10_f6_POP$NInd>0,1,0)

####merge with range info
SexAllee_ALL_POP_RNG<- merge(SexAllee_10_f6_POP, SexAllee_10_f6_0.005, all = T)
SexAllee_ALL_POP_RNG$edge<- SexAllee_ALL_POP_RNG$max_Y / 100

####create landscape slices
SexAllee_ALL_POP_RNG$rel.edge<- SexAllee_ALL_POP_RNG$edge - SexAllee_ALL_POP_RNG$y
SexAllee_ALL_POP_RNG$edge.10<- cut(SexAllee_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)

#feamle mating rate
SexAllee_ALL_POP_RNG$mate.rate<- (SexAllee_ALL_POP_RNG$Nfemales_stage1 - SexAllee_ALL_POP_RNG$unmated.fem) / SexAllee_ALL_POP_RNG$Nfemales_stage1

#Cells per slice
#1. turn "edge.10" in to the number of rows behind the range the FURTHEST cell could be to belong to that slice
SexAllee_ALL_POP_RNG$furthest.y<- (SexAllee_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range the CLOSEST cell could be to belong to that slice
SexAllee_ALL_POP_RNG$closest.y<- (SexAllee_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0


cellsINslice<- rep(NA, nrow(SexAllee_ALL_POP_RNG))

for(i in 1:nrow(SexAllee_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (SexAllee_ALL_POP_RNG$edge[i] - (SexAllee_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (SexAllee_ALL_POP_RNG$edge[i] - SexAllee_ALL_POP_RNG$closest.y[i])] )
}

SexAllee_ALL_POP_RNG$cellsINslice<- cellsINslice



#summarise by landscape slice####
SexLim_ALL_POP_RNG_summary<- ddply(SexLim_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(SexLim_ALL_POP_RNG_summary)<- c("Rep", "Year","edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
SexLim_ALL_POP_RNG_summary$mating<- as.factor('Polygyny');SexLim_ALL_POP_RNG_summary$finding<- as.factor('N');SexLim_ALL_POP_RNG_summary$settle<- as.factor('DI')


SexAllee_ALL_POP_RNG_summary<- ddply(SexAllee_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(SexAllee_ALL_POP_RNG_summary)<- c("Rep", "Year", "edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
SexAllee_ALL_POP_RNG_summary$mating<- as.factor('Monogamy');SexAllee_ALL_POP_RNG_summary$finding<- as.factor('N');SexAllee_ALL_POP_RNG_summary$settle<- as.factor('DI')



######################################################################################################################################################
#Polygyny - mate finding
LimMate_10_f6_POP<- read.table("Batch1_Sim35_Land1_Pop.txt", header=T)
LimMate_10_f6_POP$finding<-'Y'
LimMate_10_f6_POP$settle<- 'DI'
LimMate_10_f6_POP$single.fem<- ifelse(LimMate_10_f6_POP$Nmales_stage1 == 0 & LimMate_10_f6_POP$Nfemales_stage1 > 0 , 1, 0)
LimMate_10_f6_POP$unmated.fem<- LimMate_10_f6_POP$single.fem * LimMate_10_f6_POP$Nfemales_stage1
LimMate_10_f6_POP$occ<- ifelse(LimMate_10_f6_POP$Nfemales_stage1>0 & LimMate_10_f6_POP$Nmales_stage1>0, 1, 0)
LimMate_10_f6_POP$occ<- ifelse(LimMate_10_f6_POP$Nfemales_stage1>0 & LimMate_10_f6_POP$Nmales_stage1>0, 1, 0)
LimMate_10_f6_POP$occ.any<- ifelse(LimMate_10_f6_POP$NInd>0,1,0)

####merge with range info
LimMate_ALL_POP_RNG<- merge(LimMate_10_f6_POP, LimMate_10_f6_0.005, all = T)
LimMate_ALL_POP_RNG$edge<- LimMate_ALL_POP_RNG$max_Y / 100

####create landscape slices
LimMate_ALL_POP_RNG$rel.edge<- LimMate_ALL_POP_RNG$edge - LimMate_ALL_POP_RNG$y
LimMate_ALL_POP_RNG$edge.10<- cut(LimMate_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)

#feamle mating rate
LimMate_ALL_POP_RNG$mate.rate<- (LimMate_ALL_POP_RNG$Nfemales_stage1 - LimMate_ALL_POP_RNG$unmated.fem) / LimMate_ALL_POP_RNG$Nfemales_stage1

#Cells per slice
#1. turn "edge.10" in to the number of rows behind the range the FURTHEST cell could be to belong to that slice
LimMate_ALL_POP_RNG$furthest.y<- (LimMate_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range the CLOSEST cell could be to belong to that slice
LimMate_ALL_POP_RNG$closest.y<- (LimMate_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0


cellsINslice<- rep(NA, nrow(LimMate_ALL_POP_RNG))

for(i in 1:nrow(LimMate_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (LimMate_ALL_POP_RNG$edge[i] - (LimMate_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (LimMate_ALL_POP_RNG$edge[i] - LimMate_ALL_POP_RNG$closest.y[i])] )
}

LimMate_ALL_POP_RNG$cellsINslice<- cellsINslice


######################################################################################################################################################
#Monogamy - mate finding
AlleeMate_10_f6_POP<- read.table("Batch1_Sim75_Land1_Pop.txt", header=T)
AlleeMate_10_f6_POP$finding<-'Y'
AlleeMate_10_f6_POP$settle<- 'DI'
AlleeMate_10_f6_POP$unmated.fem<- AlleeMate_10_f6_POP$Nfemales_stage1 - AlleeMate_10_f6_POP$Nmales_stage1
AlleeMate_10_f6_POP$unmated.fem[AlleeMate_10_f6_POP$unmated.fem < 0]<- 0
AlleeMate_10_f6_POP$occ<- ifelse(AlleeMate_10_f6_POP$Nfemales_stage1>0 & AlleeMate_10_f6_POP$Nmales_stage1>0, 1, 0)
AlleeMate_10_f6_POP$occ.any<- ifelse(AlleeMate_10_f6_POP$NInd>0,1,0)

####merge with range info
AlleeMate_ALL_POP_RNG<- merge(AlleeMate_10_f6_POP, AlleeMate_10_f6_0.005, all = T)
AlleeMate_ALL_POP_RNG$edge<- AlleeMate_ALL_POP_RNG$max_Y / 100

####create landscape slices
AlleeMate_ALL_POP_RNG$rel.edge<- AlleeMate_ALL_POP_RNG$edge - AlleeMate_ALL_POP_RNG$y
AlleeMate_ALL_POP_RNG$edge.10<- cut(AlleeMate_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)

#feamle mating rate
AlleeMate_ALL_POP_RNG$mate.rate<- (AlleeMate_ALL_POP_RNG$Nfemales_stage1 - AlleeMate_ALL_POP_RNG$unmated.fem) / AlleeMate_ALL_POP_RNG$Nfemales_stage1

#Cells per slice
#1. turn "edge.10" in to the number of rows behind the range the FURTHEST cell could be to belong to that slice
AlleeMate_ALL_POP_RNG$furthest.y<- (AlleeMate_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range the CLOSEST cell could be to belong to that slice
AlleeMate_ALL_POP_RNG$closest.y<- (AlleeMate_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0


cellsINslice<- rep(NA, nrow(AlleeMate_ALL_POP_RNG))

for(i in 1:nrow(AlleeMate_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (AlleeMate_ALL_POP_RNG$edge[i] - (AlleeMate_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (AlleeMate_ALL_POP_RNG$edge[i] - AlleeMate_ALL_POP_RNG$closest.y[i])] )
}

AlleeMate_ALL_POP_RNG$cellsINslice<- cellsINslice




#summarise by landscape slice####
LimMate_ALL_POP_RNG_summary<- ddply(LimMate_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(LimMate_ALL_POP_RNG_summary)<- c("Rep", "Year", "edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
LimMate_ALL_POP_RNG_summary$mating<- as.factor('Polygyny');LimMate_ALL_POP_RNG_summary$finding<- as.factor('Y');LimMate_ALL_POP_RNG_summary$settle<- as.factor('DI')


AlleeMate_ALL_POP_RNG_summary<- ddply(AlleeMate_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(AlleeMate_ALL_POP_RNG_summary)<- c("Rep", "Year", "edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
AlleeMate_ALL_POP_RNG_summary$mating<- as.factor('Monogamy');AlleeMate_ALL_POP_RNG_summary$finding<- as.factor('Y');AlleeMate_ALL_POP_RNG_summary$settle<- as.factor('DI')





######################################################################################################################################################
######################################################################################################################################################
#Dendity dependent settlement rule

#****IMPORTANT NOTE***
#RangeShifter files have the same names for density dependent settlement and density independent settlement
#Ensure they are saved in a different folder such that the correct fil is called

setwd("C:/Users/whm20/OneDrive/Documents/Documents/PhD/University PC back-up/Wills_stuff/RangeShifter_chptr/Experiment_5a_DDS_HIGHSURVIVAL/Outputs")

######################################################################################################################################################
#Polygyny - no mate finding
DDS_SexLim_10_f6_POP<- read.table("Batch1_Sim15_Land1_Pop.txt", header=T)
DDS_SexLim_10_f6_POP$finding<-'N'
DDS_SexLim_10_f6_POP$settle<-'DD'
DDS_SexLim_10_f6_POP$single.fem<- ifelse(DDS_SexLim_10_f6_POP$Nmales_stage1 == 0 & DDS_SexLim_10_f6_POP$Nfemales_stage1 > 0 , 1, 0)
DDS_SexLim_10_f6_POP$unmated.fem<- DDS_SexLim_10_f6_POP$single.fem * DDS_SexLim_10_f6_POP$Nfemales_stage1
DDS_SexLim_10_f6_POP$occ<- ifelse(DDS_SexLim_10_f6_POP$Nfemales_stage1>0 & DDS_SexLim_10_f6_POP$Nmales_stage1>0, 1, 0)
DDS_SexLim_10_f6_POP$occ.any<- ifelse(DDS_SexLim_10_f6_POP$NInd>0,1,0)

####merge with range info
DDS_SexLim_ALL_POP_RNG<- merge(DDS_SexLim_10_f6_POP, DDS_SexLim_10_f6_0.005, all = T)
DDS_SexLim_ALL_POP_RNG$edge<- DDS_SexLim_ALL_POP_RNG$max_Y / 100

####create landscape slices
DDS_SexLim_ALL_POP_RNG$rel.edge<- DDS_SexLim_ALL_POP_RNG$edge - DDS_SexLim_ALL_POP_RNG$y
DDS_SexLim_ALL_POP_RNG$edge.10<- cut(DDS_SexLim_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)

#feamle mating rate
DDS_SexLim_ALL_POP_RNG$mate.rate<- (DDS_SexLim_ALL_POP_RNG$Nfemales_stage1 - DDS_SexLim_ALL_POP_RNG$unmated.fem) / DDS_SexLim_ALL_POP_RNG$Nfemales_stage1

#Cells per slice

#1. turn "edge.10" in to the number of rows behind the range the FURTHEST cell could be to belong to that slice
DDS_SexLim_ALL_POP_RNG$furthest.y<- (DDS_SexLim_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range the CLOSEST cell could be to belong to that slice
DDS_SexLim_ALL_POP_RNG$closest.y<- (DDS_SexLim_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0

cellsINslice<- rep(NA, nrow(DDS_SexLim_ALL_POP_RNG))

for(i in 1:nrow(DDS_SexLim_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (DDS_SexLim_ALL_POP_RNG$edge[i] - (DDS_SexLim_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (DDS_SexLim_ALL_POP_RNG$edge[i] - DDS_SexLim_ALL_POP_RNG$closest.y[i])] )
}

DDS_SexLim_ALL_POP_RNG$cellsINslice<- cellsINslice

######################################################################################################################################################
#Monogamy - no mate finding
DDS_SexAllee_10_f6_POP<- read.table("Batch1_Sim55_Land1_Pop.txt", header=T)
DDS_SexAllee_10_f6_POP$finding<-'N'
DDS_SexAllee_10_f6_POP$settle<-'DD'
DDS_SexAllee_10_f6_POP$unmated.fem<- DDS_SexAllee_10_f6_POP$Nfemales_stage1 - DDS_SexAllee_10_f6_POP$Nmales_stage1
DDS_SexAllee_10_f6_POP$unmated.fem[DDS_SexAllee_10_f6_POP$unmated.fem < 0]<- 0
DDS_SexAllee_10_f6_POP$occ<- ifelse(DDS_SexAllee_10_f6_POP$Nfemales_stage1>0 & DDS_SexAllee_10_f6_POP$Nmales_stage1>0, 1, 0)
DDS_SexAllee_10_f6_POP$occ.any<- ifelse(DDS_SexAllee_10_f6_POP$NInd>0,1,0)

####merge with range info
DDS_SexAllee_ALL_POP_RNG<- merge(DDS_SexAllee_10_f6_POP, DDS_SexAllee_10_f6_0.005, all = T)
DDS_SexAllee_ALL_POP_RNG$edge<- DDS_SexAllee_ALL_POP_RNG$max_Y / 100

####create landscape slices
DDS_SexAllee_ALL_POP_RNG$rel.edge<- DDS_SexAllee_ALL_POP_RNG$edge - DDS_SexAllee_ALL_POP_RNG$y
DDS_SexAllee_ALL_POP_RNG$edge.10<- cut(DDS_SexAllee_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)

#feamle mating rate
DDS_SexAllee_ALL_POP_RNG$mate.rate<- (DDS_SexAllee_ALL_POP_RNG$Nfemales_stage1 - DDS_SexAllee_ALL_POP_RNG$unmated.fem) / DDS_SexAllee_ALL_POP_RNG$Nfemales_stage1

#Cells per slice

#1. turn "edge.10" in to the number of rows behind the range the FURTHEST cell could be to belong to that slice
DDS_SexAllee_ALL_POP_RNG$furthest.y<- (DDS_SexAllee_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range the CLOSEST cell could be to belong to that slice
DDS_SexAllee_ALL_POP_RNG$closest.y<- (DDS_SexAllee_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0

cellsINslice<- rep(NA, nrow(DDS_SexAllee_ALL_POP_RNG))

for(i in 1:nrow(DDS_SexAllee_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (DDS_SexAllee_ALL_POP_RNG$edge[i] - (DDS_SexAllee_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (DDS_SexAllee_ALL_POP_RNG$edge[i] - DDS_SexAllee_ALL_POP_RNG$closest.y[i])] )
}

DDS_SexAllee_ALL_POP_RNG$cellsINslice<- cellsINslice




#summarise by landscape slice####
DDS_SexLim_ALL_POP_RNG_summary<- ddply(DDS_SexLim_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(DDS_SexLim_ALL_POP_RNG_summary)<- c("Rep", "Year", "edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
DDS_SexLim_ALL_POP_RNG_summary$mating<- as.factor('Polygyny');DDS_SexLim_ALL_POP_RNG_summary$finding<- as.factor('N');DDS_SexLim_ALL_POP_RNG_summary$settle<- as.factor('DD')


DDS_SexAllee_ALL_POP_RNG_summary<- ddply(DDS_SexAllee_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(DDS_SexAllee_ALL_POP_RNG_summary)<- c("Rep", "Year", "edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
DDS_SexAllee_ALL_POP_RNG_summary$mating<- as.factor('Monogamy');DDS_SexAllee_ALL_POP_RNG_summary$finding<- as.factor('N');DDS_SexAllee_ALL_POP_RNG_summary$settle<- as.factor('DD')


######################################################################################################################################################
#Polygyny - mate finding
DDS_LimMate_10_f6_POP<- read.table("Batch1_Sim35_Land1_Pop.txt", header=T)
DDS_LimMate_10_f6_POP$finding<-'Y'
DDS_LimMate_10_f6_POP$setle<-'DD'
DDS_LimMate_10_f6_POP$single.fem<- ifelse(DDS_LimMate_10_f6_POP$Nmales_stage1 == 0 & DDS_LimMate_10_f6_POP$Nfemales_stage1 > 0 , 1, 0)
DDS_LimMate_10_f6_POP$unmated.fem<- DDS_LimMate_10_f6_POP$single.fem * DDS_LimMate_10_f6_POP$Nfemales_stage1
DDS_LimMate_10_f6_POP$occ<- ifelse(DDS_LimMate_10_f6_POP$Nfemales_stage1>0 & DDS_LimMate_10_f6_POP$Nmales_stage1>0, 1, 0)
DDS_LimMate_10_f6_POP$occ.any<- ifelse(DDS_LimMate_10_f6_POP$NInd>0,1,0)

####merge with range info
DDS_LimMate_ALL_POP_RNG<- merge(DDS_LimMate_10_f6_POP, DDS_LimMate_10_f6_0.005, all = T)
DDS_LimMate_ALL_POP_RNG$edge<- DDS_LimMate_ALL_POP_RNG$max_Y / 100

####create landscape slices
DDS_LimMate_ALL_POP_RNG$rel.edge<- DDS_LimMate_ALL_POP_RNG$edge - DDS_LimMate_ALL_POP_RNG$y
DDS_LimMate_ALL_POP_RNG$edge.10<- cut(DDS_LimMate_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)

#feamle mating rate
DDS_LimMate_ALL_POP_RNG$mate.rate<- (DDS_LimMate_ALL_POP_RNG$Nfemales_stage1 - DDS_LimMate_ALL_POP_RNG$unmated.fem) / DDS_LimMate_ALL_POP_RNG$Nfemales_stage1

#Cells per slice

#1. turn "edge.10" in to the number of rows behind the range the FURTHEST cell could be to belong to that slice
DDS_LimMate_ALL_POP_RNG$furthest.y<- (DDS_LimMate_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range the CLOSEST cell could be to belong to that slice
DDS_LimMate_ALL_POP_RNG$closest.y<- (DDS_LimMate_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0

cellsINslice<- rep(NA, nrow(DDS_LimMate_ALL_POP_RNG))

for(i in 1:nrow(DDS_LimMate_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (DDS_LimMate_ALL_POP_RNG$edge[i] - (DDS_LimMate_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (DDS_LimMate_ALL_POP_RNG$edge[i] - DDS_LimMate_ALL_POP_RNG$closest.y[i])] )
}

DDS_LimMate_ALL_POP_RNG$cellsINslice<- cellsINslice


######################################################################################################################################################
#Monogamy - mate finding
DDS_AlleeMate_10_f6_POP<- read.table("Batch1_Sim75_Land1_Pop.txt", header=T)
DDS_AlleeMate_10_f6_POP$finding<-'Y'
DDS_AlleeMate_10_f6_POP$settle<-'DD'
DDS_AlleeMate_10_f6_POP$unmated.fem<- DDS_AlleeMate_10_f6_POP$Nfemales_stage1 - DDS_AlleeMate_10_f6_POP$Nmales_stage1
DDS_AlleeMate_10_f6_POP$unmated.fem[DDS_AlleeMate_10_f6_POP$unmated.fem < 0]<- 0
DDS_AlleeMate_10_f6_POP$occ<- ifelse(DDS_AlleeMate_10_f6_POP$Nfemales_stage1>0 & DDS_AlleeMate_10_f6_POP$Nmales_stage1>0, 1, 0)
DDS_AlleeMate_10_f6_POP$occ.any<- ifelse(DDS_AlleeMate_10_f6_POP$NInd>0,1,0)

####merge with range info
DDS_AlleeMate_ALL_POP_RNG<- merge(DDS_AlleeMate_10_f6_POP, DDS_AlleeMate_10_f6_0.005, all = T)
DDS_AlleeMate_ALL_POP_RNG$edge<- DDS_AlleeMate_ALL_POP_RNG$max_Y / 100

####create landscape slices
DDS_AlleeMate_ALL_POP_RNG$rel.edge<- DDS_AlleeMate_ALL_POP_RNG$edge - DDS_AlleeMate_ALL_POP_RNG$y
DDS_AlleeMate_ALL_POP_RNG$edge.10<- cut(DDS_AlleeMate_ALL_POP_RNG$rel.edge,seq(-1000,1000,10), include.lowest = T, labels = F, right =F)

#feamle mating rate
DDS_AlleeMate_ALL_POP_RNG$mate.rate<- (DDS_AlleeMate_ALL_POP_RNG$Nfemales_stage1 - DDS_AlleeMate_ALL_POP_RNG$unmated.fem) / DDS_AlleeMate_ALL_POP_RNG$Nfemales_stage1

#Cells per slice

#1. turn "edge.10" in to the number of rows behind the range the FURTHEST cell could be to belong to that slice
DDS_AlleeMate_ALL_POP_RNG$furthest.y<- (DDS_AlleeMate_ALL_POP_RNG$edge.10 - 100)*10
#1a turn "edge.10" in to the number of rows behind the range the CLOSEST cell could be to belong to that slice
DDS_AlleeMate_ALL_POP_RNG$closest.y<- (DDS_AlleeMate_ALL_POP_RNG$furthest.y - 10)
#2.load the landscape - hab0.07_0

cellsINslice<- rep(NA, nrow(DDS_AlleeMate_ALL_POP_RNG))

for(i in 1:nrow(DDS_AlleeMate_ALL_POP_RNG)){
  cellsINslice[i]<- length(hab0.07_0$y_coord[hab0.07_0$y_coord >=  (DDS_AlleeMate_ALL_POP_RNG$edge[i] - (DDS_AlleeMate_ALL_POP_RNG$furthest.y[i]-1)) & 
                                               hab0.07_0$y_coord <= (DDS_AlleeMate_ALL_POP_RNG$edge[i] - DDS_AlleeMate_ALL_POP_RNG$closest.y[i])] )
}

DDS_AlleeMate_ALL_POP_RNG$cellsINslice<- cellsINslice



#summarise by landscape slice####
DDS_LimMate_ALL_POP_RNG_summary<- ddply(DDS_LimMate_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(DDS_LimMate_ALL_POP_RNG_summary)<- c("Rep", "Year", "edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
DDS_LimMate_ALL_POP_RNG_summary$mating<- as.factor('Polygyny');DDS_LimMate_ALL_POP_RNG_summary$finding<- as.factor('Y');DDS_LimMate_ALL_POP_RNG_summary$settle<- as.factor('DD')


DDS_AlleeMate_ALL_POP_RNG_summary<- ddply(DDS_AlleeMate_ALL_POP_RNG, c("Rep", "Year", "edge.10", "cellsINslice"), function(df)
  c(sum(df$occ.any), sum(df$occ), mean(df$mate.rate, na.rm = T))) 
colnames(DDS_AlleeMate_ALL_POP_RNG_summary)<- c("Rep", "Year", "edge.10", "cellsINslice", "No.occ.any", "No.occ", "mean.mate.rate")
DDS_AlleeMate_ALL_POP_RNG_summary$mating<- as.factor('Monogamy');DDS_AlleeMate_ALL_POP_RNG_summary$finding<- as.factor('Y');DDS_AlleeMate_ALL_POP_RNG_summary$settle<- as.factor('DD')


######################################################################################################################################################
######################################################################################################################################################

#Bind all together
All_mate_rate<- rbind(SexLim_ALL_POP_RNG_summary, SexAllee_ALL_POP_RNG_summary, LimMate_ALL_POP_RNG_summary, AlleeMate_ALL_POP_RNG_summary,
                      DDS_SexLim_ALL_POP_RNG_summary, DDS_SexAllee_ALL_POP_RNG_summary, DDS_LimMate_ALL_POP_RNG_summary, DDS_AlleeMate_ALL_POP_RNG_summary)

#Calculate proportional occupancy (breeding populations) in each neighborhood (i.e. landscape slice)
All_mate_rate$prop.breed.occ<- All_mate_rate$No.occ / All_mate_rate$cellsINslice

#Calculate proportional occupancy (all breeding and non-breeding populations) in each neighborhood (i.e. landscape slice)
All_mate_rate$prop.any.occ<- All_mate_rate$No.occ.any / All_mate_rate$cellsINslice

#Split prop.breed.occ in to 10 equally sized bins, each 0.1 wide
All_mate_rate$prop.breed.occ.fac<-  cut(All_mate_rate$prop.breed.occ, seq(0,1,0.1), include.lowest = T, labels = F, right =F)


#Subset of the part of the range we want to see i.e. the edge and a few rows back
All_mate.rate.edge2core<- All_mate_rate[All_mate_rate$edge.10>100 & All_mate_rate$edge.10 < 107,]

#create label objects
mating.labs <- c("Polygyny", "Monogamy")
names(mating.labs)<- c("Lim","Allee")
finding.labs <- c("No mate finding","Mate finding")
names(finding.labs)<- c("N","Y")
settle.labs <- c("Density independent settlement","Density dependent settlement")
names(settle.labs)<- c("DI","DD")


#calculate overall mean and SE for mating rate in each landscape slice
tester1<- ddply(All_mate_rate,c("mating","finding","settle","prop.breed.occ.fac"),summarise, mean.of.mean=mean(mean.mate.rate,na.rm=T), se.of.mean=se(mean.mate.rate))

tester1$prop.breed.occ.fac<- tester1$prop.breed.occ.fac/10 # for plotting purposes


######################################################################################################################################################
#
#Create plot insets that show how occupancy changes from edge to core - these use the overall mean occupancy in landscape slices

All_mean_occ.summary<- ddply(All_mate.rate.edge2core, c("edge.10","mating","settle","finding"), summarise,
                             mean.any.occ=mean(prop.any.occ,na.rm = T),mean.breed.occ=mean(prop.breed.occ,na.rm = T), all.mean.mate.rate=mean(mean.mate.rate,na.rm = T))


p1<- ggplot()+
  geom_line(data=All_mean_occ.summary[All_mean_occ.summary$mating=='Polygyny' & All_mean_occ.summary$settle=='DI',], 
            aes(x=factor(edge.10),y=mean.breed.occ, colour=finding, group=finding))+
  scale_colour_viridis(option='magma', discrete=T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  scale_linetype(guide=F)+
  ylab("P")+
  xlab(expression(Edge %<->% Core))+
  theme_classic()+
  theme( strip.background = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1))+
  scale_x_discrete(labels = c("Edge", "", "", "", "", ""))+
  geom_hline(yintercept = 0.8, linetype=2)

p2<- ggplot()+
  geom_line(data=All_mean_occ.summary[All_mean_occ.summary$mating=='Polygyny' & All_mean_occ.summary$settle=='DD',], 
            aes(x=factor(edge.10),y=mean.breed.occ, colour=finding, group=finding))+
  scale_colour_viridis(option='magma', discrete=T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  scale_linetype(guide=F)+
  ylab("P")+
  xlab(expression(Edge %<->% Core))+
  theme_classic()+
  theme( strip.background = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1))+
  scale_x_discrete(labels = c("Edge", "", "", "", "", ""))+
  geom_hline(yintercept = 0.8, linetype=2)

p3<- ggplot()+
  geom_line(data=All_mean_occ.summary[All_mean_occ.summary$mating=='Monogamy' & All_mean_occ.summary$settle=='DI',], 
            aes(x=factor(edge.10),y=mean.breed.occ, colour=finding, group=finding))+
  scale_colour_viridis(option='magma', discrete=T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  scale_linetype(guide=F)+
  ylab("P")+
  xlab(expression(Edge %<->% Core))+
  theme_classic()+
  theme( strip.background = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1))+
  geom_hline(yintercept = 0.8, linetype=2)

p4<- ggplot()+
  geom_line(data=All_mean_occ.summary[All_mean_occ.summary$mating=='Monogamy' & All_mean_occ.summary$settle=='DD',], 
            aes(x=factor(edge.10),y=mean.breed.occ, colour=finding, group=finding))+
  scale_colour_viridis(option='magma', discrete=T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  scale_linetype(guide=F)+
  ylab("P")+
  xlab(expression(Edge %<->% Core))+
  theme_classic()+
  theme( strip.background = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1))+
  geom_hline(yintercept = 0.8, linetype=2)

######################################################################################################################################################
#
#Create plots of proportion unmated females as a function of Neighborhood occupancy
#Including inserts created above

just.a.line<- data.frame(x=0:1,y=1:0)
just.another.line<- data.frame(x=c(0.8,0.8),y=c(0,0.6))

#create standard jitter to avoid error bar overlap
jitter2<- c(rep(c(0.005,-0.005),each=20),rep(c(0.01,-0.01),each=20))
tester1$jitter.x<-tester1$prop.breed.occ.fac+jitter2

m1<- ggplot()+
  annotation_custom(ggplotGrob(p1), xmin = 0.4, xmax = 1, 
                    ymin = 0.6, ymax = 1.1)+
  geom_point(data=tester1[tester1$mating=='Polygyny' & tester1$settle=='DI',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding),size=1.25,shape=16)+
  geom_errorbar(data=tester1[tester1$mating=='Polygyny' & tester1$settle=='DI',], aes(x=jitter.x, ymax=1-mean.of.mean+se.of.mean, ymin=1-mean.of.mean-se.of.mean, colour=finding),width=0)+
  geom_line(data=tester1[tester1$mating=='Polygyny' & tester1$settle=='DI',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding, group=finding))+
  ylab("Proportion unmated females")+
  xlab("")+
  theme_classic()+
  theme( strip.background = element_blank())+
  scale_colour_viridis(option='magma', discrete = T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  ylim(0,1)+
  ggtitle("Polygyny", "Density independent settlement")+
  coord_fixed()+
  geom_line(data=just.another.line,aes(x=x,y=y), linetype=2,size=1)+
  annotate("text", x = 0.2, y = 1, label = "a")


m2<- ggplot()+
  annotation_custom(ggplotGrob(p2), xmin = 0.4, xmax = 1, 
                    ymin = 0.6, ymax = 1.1)+
  geom_point(data=tester1[tester1$mating=='Polygyny' & tester1$settle=='DD',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding),size=1.25,shape=16)+
  geom_errorbar(data=tester1[tester1$mating=='Polygyny' & tester1$settle=='DD',], aes(x=jitter.x, ymax=1-mean.of.mean+se.of.mean, ymin=1-mean.of.mean-se.of.mean, colour=finding),width=0)+
  geom_line(data=tester1[tester1$mating=='Polygyny' & tester1$settle=='DD',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding,group=finding))+
  ylab("Proportion unmated females")+
  xlab("")+
  theme_classic()+
  theme( strip.background = element_blank())+
  scale_colour_viridis(option='magma', discrete = T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  ylim(0,1)+
  ggtitle("", "Density dependent settlement")+
  coord_fixed()+
  geom_line(data=just.another.line,aes(x=x,y=y), linetype=2,size=1)+
  annotate("text", x = 0.2, y = 1, label = "b")

m3<- ggplot()+
  annotation_custom(ggplotGrob(p3), xmin = 0.4, xmax = 1, 
                    ymin = 0.6, ymax = 1.1)+
  geom_point(data=tester1[tester1$mating=='Monogamy' & tester1$settle=='DI',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding),size=1.25,shape=16)+
  geom_errorbar(data=tester1[tester1$mating=='Monogamy' & tester1$settle=='DI',], aes(x=jitter.x, ymax=1-mean.of.mean+se.of.mean, ymin=1-mean.of.mean-se.of.mean, colour=finding),width=0)+
  geom_line(data=tester1[tester1$mating=='Monogamy' & tester1$settle=='DI',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding,group=finding))+
  ylab("Proportion unmated females")+
  xlab("Neighbourhood occupancy (P)")+
  theme_classic()+
  theme( strip.background = element_blank())+
  scale_colour_viridis(option='magma', discrete = T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  ylim(0,1)+
  ggtitle("Monogamy", "Density independent settlement")+
  coord_fixed()+
  geom_line(data=just.another.line,aes(x=x,y=y), linetype=2,size=1)+
  annotate("text", x = 0.2, y = 1, label = "c")


m4<- ggplot()+
  annotation_custom(ggplotGrob(p4), xmin = 0.4, xmax = 1, 
                    ymin = 0.6, ymax = 1.1)+
  geom_point(data=tester1[tester1$mating=='Monogamy' & tester1$settle=='DD',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding),size=1.25,shape=16)+
  geom_errorbar(data=tester1[tester1$mating=='Monogamy' & tester1$settle=='DD',], aes(x=jitter.x, ymax=1-mean.of.mean+se.of.mean, ymin=1-mean.of.mean-se.of.mean, colour=finding),width=0)+
  geom_line(data=tester1[tester1$mating=='Monogamy' & tester1$settle=='DD',], aes(x=jitter.x, y=1-mean.of.mean, colour=finding,group=finding))+
  ylab("Proportion unmated females")+
  xlab("Neighbourhood occupancy (P)")+
  theme_classic()+
  theme( strip.background = element_blank())+
  scale_colour_viridis(option='magma', discrete = T, end=0.8, guide=F)+
  scale_shape(guide=F)+
  ylim(0,1)+
  ggtitle("", "Density dependent settlement")+
  coord_fixed()+
  geom_line(data=just.another.line,aes(x=x,y=y), linetype=2,size=1)+
  annotate("text", x = 0.2, y = 1, label = "d")


#Plot the whole lot
ggarrange(m1,m2,m3,m4,nrow=2,ncol=2)



#save as high resolution jpeg
#ggsave("neighbourhood occ vs mate rate_WITH INSETS.jpeg",dpi=300,device="jpeg",width=18,height=18, units="cm")






