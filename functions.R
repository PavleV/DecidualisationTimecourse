# load libraries

library(tidyverse)
library(ggplot2)


# load data

AH_EL_RNA_ALLREPS <- read.csv("./Data/AH_EL_RNA_ALLREPS.csv")


# Functions

# subset data according to gene name or ensembl ID
subset.FUN <- function(mydata=AH_EL_RNA_ALLREPS, geneName=NULL, ensemblID=NULL){

  if(!is.null(geneName) & !is.null(ensemblID)){
    return(subset(mydata, GeneName==geneName & EnsemblID ==ensemblID))
  }
  if(!is.null(geneName) & is.null(ensemblID)){
    return(subset(mydata, GeneName==geneName))
  }
  if(is.null(geneName) & !is.null(ensemblID)){
    return(subset(mydata, EnsemblID ==ensemblID))
  }

}

#arrange subsetted data for plotting

arrange.FUN <- function(mydata, geneName=NULL, ensemblID=NULL){

  mydata.subset <- subset.FUN(mydata=mydata, geneName=geneName, ensemblID=ensemblID) %>%
    pivot_longer(cols = colnames(mydata)[2:40])

  mydata.subset <- cbind(mydata.subset,str_split(mydata.subset$name,pattern="_", simplify = T)) %>%
    colnames()[3:6] <- c("Sample","TPM","Biopsy","Hours")

  return(mydata.subset)

}




plotRNA.FUN <- function(mydata, geneName=NULL, ensemblID=NULL, biopsies = c("S169","S170","S506","S507"), times = c(0,3,6,9,12,18,24,36,48,96)){

  data.sub <- arrange.FUN(mydata=mydata, geneName=geneName, ensemblID=ensemblID)

  data.sub <- subset(data.sub, Biopsy %in% biopsies & Hours %in% times)

  p1 <- ggplot(data.sub, aes(x=as.numeric(Hours),y=TPM))+
    stat_summary(geom = "line", fun.y = mean)+
    theme_bw()+theme(legend.position="none")+
    scale_x_continuous(breaks=times, minor_breaks = NULL)+
    labs(y="Transcripts per Million", x="Hours")


  print(p1)
}

