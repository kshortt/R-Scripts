---
title: "Plot_ReadCTs_KS"
author: "KAShortt"
date: "2017"
output: html_document
---
  
## The purpose of this program is to graph read count data.

R.version.string
Sys.Date()

#set the working path
setwd("U:\APAP screen_crispr liver\Analysis\Test")

library("ggplot2")

#import data 
MyData <- read.table(file="Test_All_andPlasmid_rep1_rep2.normalized.txt",sep="\t", header=T,stringsAsFactor=F)


data2 <- cbind(MyData[,3:4],MyData[,7:14],MyData[,5:6])

#log2 (readcounts) for each sample
data <- log(data2,base=2)



#boxplot for each sample
par(font.axis = 2,font.lab=2,family="sans",ps=12,las=2)
boxplot(data, frame = FALSE,border = "steelblue",ylim=c(-5,20),ylab = "log2 sgRNA counts")



data3 <- cbind(data[,11:12],data[,1:10])

data4<-array(dim=c(12,12))
for (i in 1:ncol(data3))
{
  for (j in 1:ncol(data3))
  {
    pv<- wilcox.test(data3[,i], data3[,j])$p.value
    data4[i,j] <- pv    
  }
}



colnames(data4) <- c("Plasmid_rep1","Plasmid_rep2",colnames(data[,1:10]))
rownames(data4) <- c("Plasmid_rep1","Plasmid_rep2",colnames(data[,1:10]))


data4


library(gplots)

col<- colorRampPalette(c("steelblue","white"))(50)
par(font.axis = 2,font.lab=2,family="sans",ps=10)

heatmap.2(data4, col=col, symm=TRUE,scale ="none",Rowv=NA,Colv=NA,key=TRUE,density.info="none", trace="none", key.title ="p value", key.xlab = "", key.ylab = "",keysize = 1.5)

## End of file