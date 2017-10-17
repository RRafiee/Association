##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, October 2017
# Research Fellow, Queen's University Belfast
# Box plot and t-test 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
############ Box plot ############# 17/10/17
setwd("/home/reza/Documents/Live/")
library(reshape2)
library(lattice)
#dataset <- read.table("sample.txt",  header = TRUE,na.strings="NA", dec=".", strip.white=TRUE)
dataset <- read.csv("LUSC_DDRDSubgroup4_StageII_III_Treated_EMTScore_FABRCAAlterations__DFS_Boxplot_171017.csv")
dataset <- read.csv("LUSC_1DDRDSubgroups_StageII_III_Treated_EMTOrdered_171017_boxplot.csv")
dataset <- read.csv("LUAD_1DDRDSubgroups_StageII_III_OrderbyDDRDScore_withEMT_206_171017_boxplot.csv")
dat.m <- melt(dataset,id.vars='Label')
#bwplot(value~Label | paste0(variable), data=dat.m,   main="Mine vs Other", layout=c(2,1), par.settings = list(box.rectangle = list(fill= rep(c('blue','red'),2))))
maintext <- paste(levels(dataset$Label)[1], " vs. ",levels(dataset$Label)[2],sep = "")
bwplot(value~Label | paste0(variable), data=dat.m,   main=maintext, layout=c(1,1),aspect = "fill", par.settings = list(box.rectangle = list(fill= rep(c('white','grey'),2))))
t.test(value ~ Label, data = dat.m)
#as.character(dataset$Label) #levels of the factor 

# ##################################################################################################
# #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # End
# ##################################################################################################
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##################################################################################################
