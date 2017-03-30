##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Finding association between clinical, pathological and molecular features of an infant cohort
# Written by Dr Reza Rafiee
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script gets a csv file including all variables for finding association 

# Install all required packages as well as dependencies
install_all_packages_automatic(ellipse)
install_all_packages_automatic(corrgram)

# Loading libraries
library(corrgram)
library(ellipse)
library(corrplot)

## Load data
data1 <- read.csv("InfantCohort14Features.csv", header=T) # updated 2017, 
DataFeturesMat <- data1[,1:14]

colnames(DataFeturesMat)
dim(DataFeturesMat)
#[1] 174  14

## Categories to test
categories <- c(1:dim(DataFeturesMat)[2])  # for the whole cohort


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
      
      
      Flag1 <-  ifelse((i==4&&j==5),F,ifelse((i==5&&j==4),F,      # DN_MBEN(col. 1), CLA (col. 2), LCA(in col.3)
                ifelse((i==4&&j==6),F,ifelse((i==6&&j==4),F,
                ifelse((i==5&&j==6),F,ifelse((i==6&&j==5),F,

                ifelse((i==7&&j==8),F,ifelse((i==8&&j==7),F,      #SHH(7), Grp3(8), Grp4(9) and Grp3Grp4(10)
                ifelse((i==7&&j==9),F,ifelse((i==9&&j==7),F,
                ifelse((i==7&&j==10),F,ifelse((i==10&&j==7),F,
                ifelse((i==8&&j==9),F,ifelse((i==9&&j==8),F,
                ifelse((i==8&&j==10),F,ifelse((i==10&&j==8),F,
                ifelse((i==9&&j==10),F,ifelse((i==10&&j==9),F,T))))))))))))))))))
      

      if((i!=j) && (Flag1))
      {  
        ## Select categories and get rid of missing data
        test <- cbind(data1[,categories[i]],data1[,categories[j]])   
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

filename_association <- paste("FeatureAssociationAnalysisInfant.csv")
write.csv(results_Total_1, file=filename_association)


data_sample_indexed <- as.matrix(results_Total_1)
data_sample_indexed <- data_sample_indexed[,-3]  # remove column test (Fisher)

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


Corr_Mat_correctedBH <- matrix(nrow=14, ncol=14,"")
colnames(Corr_Mat_correctedBH) <- colnames(data1) # c("Male","TRvsNotTR","MStage2plus","DN_MBEN","Classic","LCA","SHH","G3","G4","G3/G4","P53.IHC","MYCamp","MYCNamp","iso17q")
rownames(Corr_Mat_correctedBH) <- colnames(data1) 


for (ix in 1:nrow(data_sample_indexed))
{
  ic1 <- as.integer(data_sample_indexed[ix,1])
  ic2 <- as.integer(data_sample_indexed[ix,2])
  Corr_Mat_correctedBH[ic1,ic2] <- data_sample_indexed[ix,5] # corrrected p-values: Benjamini.Hochberg
  Corr_Mat_correctedBH[ic2,ic1] <- data_sample_indexed[ix,5]  
}

Corr_Mat_correctedBH <- apply(Corr_Mat_correctedBH,1,as.numeric)

colnames(Corr_Mat_correctedBH) <- colnames(data1) #c("Male","TRvsNotTR","MStage2plus","DN_MBEN","Classic","LCA","SHH","G3","G4","G3/G4","P53.IHC","MYCamp","MYCNamp","iso17q")
rownames(Corr_Mat_correctedBH) <- colnames(data1) #c("Male","TRvsNotTR","MStage2plus","DN_MBEN","Classic","LCA","SHH","G3","G4","G3/G4","P53.IHC","MYCamp","MYCNamp","iso17q")

sample_data1 <- Corr_Mat_correctedBH

colnames(Corr_Mat_correctedBH) <- colnames(data1)
rownames(Corr_Mat_correctedBH) <- colnames(data1)

#sample_data1 <- replace(sample_data1, )
sample_data1[sample_data1 > 0.05 ] <- 1

par(mfrow=c(1,1))
corrgram(sample_data1,lower.panel=panel.shade, upper.panel=NULL,col.regions=colorRampPalette(c("black","darkgreen","grey")))
plotcorr(sample_data1, col = colorRampPalette(c("white", "navy"))(100))

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
