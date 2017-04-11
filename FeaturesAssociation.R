##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Finding association between clinical, pathological and molecular features of an infant cohort
# Written by Dr Reza Rafiee
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script gets a csv file including all variables for finding association using Fisher or chi-squared test

# Install all required packages as well as dependencies
# install_all_packages_automatic(ellipse)
# install_all_packages_automatic(corrgram)

# Loading libraries
library(corrgram)
library(ellipse)
library(corrplot)


setwd("~/My Projects at NICR/2014/Infant/Our New Material")
## loading main data (whole cohort: all subgroups)
data_main <- read.csv("InfantCohort14Features.csv", header=T) # updated 2017, 
data_main <- data_main[,-15] # removing column 15 ("CDKN0AB")



# Choose one of the following options and cohort:
#++++++++++++++++++++++++++++++++++++++++++++++++++++
#1)++++++++++++++++++++++++++++++++++++++++++++++++++++
## loading data, only 14 first features (without CNVs)
DataFeturesMat <- data_main[,1:14]
filename_association <- paste("FeatureAssociationAnalysisInfant.csv")
filename_association_pvalues <- paste("FeatureAssociationAnalysisInfant_pvalues.csv")
Flag_subgroup <- FALSE

#++++++++++++++++++++++++++++++++++++++++++++++++++++
#2)++++++++++++++++++++++++++++++++++++++++++++++++++++
# loading data,choosing samples with SHH molecular subgroup
DataFeturesMat <- data_main[which(data_main$SHH == 1),]
# DataFeturesMat$SHH  # in column 7
# DataFeturesMat$Grp3 # in column 8
# DataFeturesMat$Grp4 # in column 9
# DataFeturesMat$Grp3Grp4 # in column 10
DataFeturesMat <- DataFeturesMat[,-c(7:10)] # exclude subgroup features
#choosing SHH samples with having pathology data
#DataFeturesMat <- Data_main[c(which(Data_main$DN_MBEN == 0), which(Data_mainMat$DN_MBEN == 1)),]
#DataFeturesMat$Loss.9q # in column 53
#DataFeturesMat$Gain.9p # in column 50
CNVs_index <- c(which(colnames(DataFeturesMat) == c("Loss.9q","Gain.9p")))
DataFeturesMat <- DataFeturesMat[,c(1:10,CNVs_index)]
filename_association <- paste("FeatureAssociationAnalysisInfant_SHH.csv")
filename_association_pvalues <- paste("FeatureAssociationAnalysisInfant_SHH_pvalues.csv")
Flag_subgroup <- TRUE 

#++++++++++++++++++++++++++++++++++++++++++++++++++++
#3)++++++++++++++++++++++++++++++++++++++++++++++++++++
# loading data,choosing samples with Grp3 molecular subgroup
DataFeturesMat <- data_main[which(data_main$Grp3 == 1),]
# choosing Grp3 samples with having pathology data
# DataFeturesMat <- Data_main[c(which(Data_main$DN_MBEN == 0), which(Data_mainMat$DN_MBEN == 1)),]
DataFeturesMat <- DataFeturesMat[,-c(7:10)] # exclude subgroup features
# DataFeturesMat$Gain.1q 
# DataFeturesMat$Loss.3 
# DataFeturesMat$Loss.4q 
# DataFeturesMat$Gain.5
# DataFeturesMat$Gain.7 
# DataFeturesMat$Loss.8 
# DataFeturesMat$Loss.10 
# DataFeturesMat$Loss.10q
# DataFeturesMat$Loss.11
# DataFeturesMat$Loss.11q
# DataFeturesMat$Loss.13q
# DataFeturesMat$Loss.15q
# DataFeturesMat$Loss.16
# DataFeturesMat$Loss.16q
# DataFeturesMat$Loss.17p
# DataFeturesMat$Gain.17
# DataFeturesMat$iso17q
# DataFeturesMat$Gain.18
# DataFeturesMat$Loss.20

CNVs_index <- c(which(colnames(DataFeturesMat) == c("Gain.1q","Loss.3")),which(colnames(DataFeturesMat) == "Loss.4q"), which(colnames(DataFeturesMat) == "Gain.5"),
     which(colnames(DataFeturesMat) == "Gain.7"), which(colnames(DataFeturesMat) == "Loss.8"), which(colnames(DataFeturesMat) == "Loss.10"),
     which(colnames(DataFeturesMat) == "Loss.10q"), which(colnames(DataFeturesMat) == "Loss.11"), which(colnames(DataFeturesMat) == "Loss.11q"), 
     which(colnames(DataFeturesMat) == "Loss.13q"), which(colnames(DataFeturesMat) == "Loss.15q"), which(colnames(DataFeturesMat) == "Loss.16"),
     which(colnames(DataFeturesMat) == "Loss.16q"), which(colnames(DataFeturesMat) == "Loss.17p"),which(colnames(DataFeturesMat) == "Gain.17q"),
     which(colnames(DataFeturesMat) == "Gain.18"), which(colnames(DataFeturesMat) == "Loss.20"))

DataFeturesMat <- DataFeturesMat[,c(1:10,CNVs_index)]

filename_association <- paste("FeatureAssociationAnalysisInfant_Grp3.csv")
filename_association_pvalues <- paste("FeatureAssociationAnalysisInfant_Grp3_pvalues.csv")
Flag_subgroup <- TRUE
#++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++


colnames(DataFeturesMat)
dim(DataFeturesMat)
#[1] 174  14 (All subgroups)
#[1] 61 16 (SHH)
#[1] 64 28 (Grp3)

## Categories to test for the selected features 
categories <- c(1:dim(DataFeturesMat)[2]) 

data1 <- DataFeturesMat
tests <- vector()
category1 <- vector()
category2 <- vector()
pvalues <- vector()
  
Flag1 <- TRUE

for(i in 1:length(categories))
  {
    j <- i+1
    while (j <= length(categories))
     {
      
      ## Check that i doesn't equal j and Removing non-Meaningfull comparsion (By Reza, 14 May 2014)
      
      
      if (Flag_subgroup)
       {
        Flag1 <-  ifelse((i==4&&j==5),F,ifelse((i==5&&j==4),F,      # DN_MBEN(col. 4), CLA (col. 5), LCA(col. 6)
                                               ifelse((i==4&&j==6),F,ifelse((i==6&&j==4),F,
                                                                            ifelse((i==5&&j==6),F,ifelse((i==6&&j==5),F,T))))))
                                                                                                         
       } 
      else {
        
        Flag1 <-  ifelse((i==4&&j==5),F,ifelse((i==5&&j==4),F,      # DN_MBEN(col. 1), CLA (col. 2), LCA(in col.3)
                  ifelse((i==4&&j==6),F,ifelse((i==6&&j==4),F,
                  ifelse((i==5&&j==6),F,ifelse((i==6&&j==5),F,

                  ifelse((i==7&&j==8),F,ifelse((i==8&&j==7),F,      #SHH(7), Grp3(8), Grp4(9) and Grp3Grp4(10)
                  ifelse((i==7&&j==9),F,ifelse((i==9&&j==7),F,
                  ifelse((i==7&&j==10),F,ifelse((i==10&&j==7),F,
                  ifelse((i==8&&j==9),F,ifelse((i==9&&j==8),F,
                  ifelse((i==8&&j==10),F,ifelse((i==10&&j==8),F,
                  ifelse((i==9&&j==10),F,ifelse((i==10&&j==9),F,T))))))))))))))))))
       }
      
      
      if((i!=j) && (Flag1))
      {  
        ## Select categories and get rid of missing data
        test <- cbind(data1[,categories[i]],data1[,categories[j]])   # Original is data, 13th October 2016
        ## Get rid of missing data
        
        #test <- subset(test,test[,1]!=9) #Original
        #test <- subset(test,test[,2]!=9) #Original
        
        test <- subset(test,test[,1]!="NA")
        test <- subset(test,test[,2]!="NA")
        
        
        ## Check if either column has a single category and ignore if
        ## that is the case
        
        if(dim(table(test[,1])) > 1 & dim(table(test[,2])) > 1)  #&
        {
          
          ## If table == 2X2, do Fisher, else do Chi squared
          
          if(dim(table(test[,1],test[,2]))[1] == 2 &       #&
             (dim(table(test[,1],test[,2]))[2] == 2))
          {
            # do fisher test
            
            pval <- fisher.test(table(test[,1],test[,2]))$p.value
            pvalTest <- "Fisher"
          } else { # do chi squared
            pval <- chisq.test(table(test[,1],test[,2]))$p.value
            pvalTest <- "Chi"
          }
          
          tests <- c(tests,pvalTest)
          pvalues <- c(pvalues,pval)
          category1 <- c(category1,colnames(data1)[categories[i]])
          category2 <- c(category2,colnames(data1)[categories[j]])
        } # end if greater than one category in each column
        
      } # end if(i!=j)
      j <- j+1
    } # end j while loop
} # end i loop
  
results <- data.frame(category1,category2,tests,pvalues)
  
## Adjusting p-values for multiple comparsions using benferroni
results_adjusted_1 <- p.adjust(as.numeric(results[,4]), method = "bonferroni", n = nrow(results))
  
## Adjusting p-values for multiple comparsions using Benjamini-Hochberg or "fdr"  (1995)
results_adjusted_2 <- p.adjust(as.numeric(results[,4]), method = "BH", n = nrow(results))
  
# in Benjamini-Hochberg, it replaces all the missing p-values with 1. This would be conservative.
#Use multiple imputation where you repeatedly impute the missing p-values with values between 0.05 and 1.
#This would also be acceptable, since conditional on the p-value being larger than 0.05 and the null hypothesis being true
#for all of them (conservative assumption), this would be the actual distribution. 

## Adjusting p-values for multiple comparsions using Benjamini-Yekutieli (2001)
results_adjusted_3 <- p.adjust(as.numeric(results[,4]), method = "BY", n = nrow(results))

results_Total_1 <- cbind(results,results_adjusted_1,results_adjusted_2,results_adjusted_3)
  
colnames (results_Total_1) <- c("category1","category2","tests","Uncorrected_pvalues","Bonferroni","Benjamini-Hochberg","Benjamini-Yekutieli")

write.csv(results_Total_1, file=filename_association)


data_sample_indexed <- as.matrix(results_Total_1)
data_sample_indexed <- data_sample_indexed[,-3]  # remove column test (Fisher)

# Replacing variable names with numbers
for(k in 1:ncol(data1))
{
  data_sample_indexed[which(results[,1] == colnames(data1)[k]),1] <- as.integer(c(rep.int(k,length(which(results[,1] == colnames(data1)[k])))))
  if (length(which(results[,2] == colnames(data1)[k])) != 0)
  {
    data_sample_indexed[which(results[,2] == colnames(data1)[k]),2] <- as.integer(c(rep.int(k,length(which(results[,2] == colnames(data1)[k])))))
  }
}

data_sample_indexed <- as.data.frame(data_sample_indexed) # from the previous level

data_sample_indexed <- t(apply(data_sample_indexed,1,as.numeric))
colnames(data_sample_indexed) <- c("category1","category2","Uncorrected_pvalues","Bonferroni","Benjamini-Hochberg","Benjamini-Yekutieli")


Corr_Mat_correctedBH <- matrix(nrow=length(categories), ncol=length(categories),"")
colnames(Corr_Mat_correctedBH) <- colnames(data1) 
rownames(Corr_Mat_correctedBH) <- colnames(data1) 


for (ix in 1:nrow(data_sample_indexed))
{
  #ix <- 1
  ic1 <- as.integer(data_sample_indexed[ix,1])
  ic2 <- as.integer(data_sample_indexed[ix,2])
  Corr_Mat_correctedBH[ic1,ic2] <- data_sample_indexed[ix,5] # corrrected p-values: Benjamini.Hochberg
  Corr_Mat_correctedBH[ic2,ic1] <- data_sample_indexed[ix,5] # corrrected p-values: Benjamini.Hochberg
}

Corr_Mat_correctedBH <- apply(Corr_Mat_correctedBH,1,as.numeric)

colnames(Corr_Mat_correctedBH) <- colnames(data1) 
rownames(Corr_Mat_correctedBH) <- colnames(data1) 


write.csv(Corr_Mat_correctedBH, file=filename_association_pvalues)

sample_data1 <- Corr_Mat_correctedBH

colnames(Corr_Mat_correctedBH) <- colnames(data1)
rownames(Corr_Mat_correctedBH) <- colnames(data1)

#sample_data1 <- replace(sample_data1, )
sample_data1[sample_data1 > 0.05 ] <- 1

par(mfrow=c(1,1))
corrgram(sample_data1,lower.panel=panel.shade, upper.panel=NULL,col.regions=colorRampPalette(c("black","darkgreen","grey")))
plotcorr(sample_data1, outline = TRUE, col = 'red', numbers = FALSE,type = "full",
       xlab = "Association between molecular, clinical and pathology features", ylab = "",
       cex = 0.75*par("cex"), mar = 0.1 + c(2,2,4,2), bty = "n")


#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End - Part One
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
