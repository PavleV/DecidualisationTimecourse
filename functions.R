# load libraries

library(tidyverse)
library(ggplot2)
library(GenomicRanges)


# load data

AH_EL_RNA_ALLREPS <- read.csv("./Data/AH_EL_RNA_ALLREPS.csv")


# Functions

# subset data according to gene name or ensembl ID
subset.RNA.FUN <- function(mydata=AH_EL_RNA_ALLREPS, geneName=NULL, ensemblID=NULL){

  if(!is.null(geneName) & !is.null(ensemblID)){
    return(subset(mydata, GeneName %in% geneName & EnsemblID == ensemblID))
  }
  if(!is.null(geneName) & is.null(ensemblID)){
    return(subset(mydata, GeneName %in% geneName))
  }
  if(is.null(geneName) & !is.null(ensemblID)){
    return(subset(mydata, EnsemblID == ensemblID))
  }

}

#arrange subsetted data for plotting

arrange.FUN <- function(mydata, geneName=NULL, ensemblID=NULL){

  mydata.subset <- subset.RNA.FUN(mydata=mydata, geneName = geneName, ensemblID=ensemblID)
  mydata.subset <-  pivot_longer(mydata.subset, cols = colnames(mydata)[2:40])

  mydata.subset <- cbind(mydata.subset,str_split(mydata.subset$name,pattern="_", simplify = T))
  colnames(mydata.subset)[3:6] <- c("Sample","TPM","Biopsy","Hours")

  return(mydata.subset)

}




plotRNA.FUN <- function(mydata, geneName=NULL, ensemblID=NULL, biopsies = c("S169","S170","S506","S508"), times = c(0,3,6,9,12,18,24,36,48,96)){

  if(is.null(geneName) & is.null(ensemblID)){
    return(
      print(
        ggplot()+theme_bw()
      )
    )
  }

  data.sub <- arrange.FUN(mydata=mydata, geneName=geneName, ensemblID=ensemblID)

  data.sub <- subset(data.sub, Biopsy %in% biopsies & Hours %in% times)

  if(length(geneName) == 1){

  p1 <- ggplot(subset(data.sub))+
    geom_line(aes(x=as.numeric(Hours),y=TPM,colour=Biopsy), size = 2)+
    stat_summary(aes(x=as.numeric(Hours),y=TPM), geom = "line", fun.y = mean, size = 2)+
    theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=10), title= element_text(size=20) )+
    scale_x_continuous(breaks=times, minor_breaks = NULL)+
    labs(y="Transcripts per Million", x="Hours")+ggtitle(paste0(geneName,ensemblID))

  }

  if(length(geneName) >= 2){

    p1 <- ggplot(subset(data.sub))+
      stat_summary(aes(x=as.numeric(Hours),y=TPM,colour=GeneName), geom = "line", fun.y = mean, size = 2)+
      theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=10), title= element_text(size=20) )+
      scale_x_continuous(breaks=times, minor_breaks = NULL)+
      labs(y="Transcripts per Million", x="Hours")#+ggtitle(paste0(geneName,ensemblID))

  }


  print(p1)
}



#load ATAC count matrix and genomic coordinates

ATAC_countsmatrix_cleaned <- read.csv("./Data/ATAC_countsmatrix_cleaned_230616.csv", row.names=1)

Alex_ATAC_peaks_hg19_201110 <- read.delim("./Data/Alex_ATAC_peaks_hg19_201110.bed", header=FALSE)

AllPeaks.granges <- GRanges(seqnames = Alex_ATAC_peaks_hg19_201110[,1], ranges = IRanges(start=Alex_ATAC_peaks_hg19_201110[,2], end=Alex_ATAC_peaks_hg19_201110[,3]), mcols = data.frame(PeakID = Alex_ATAC_peaks_hg19_201110[,4]))


# function for extracting genomic coordinates from text input and generating genomic ranges object

extractCoord <- function(string.input){

  coordinates <- str_split_1(string.input, pattern=regex("[:,_-[:space:]]"))

  return(GRanges(seqnames = coordinates[1], ranges = IRanges(start=as.numeric(coordinates[2]), end=as.numeric(coordinates[3]))))


}

# function for subsetting data matrix given a set of coordinates

subset.ATAC.FUN <- function(mydata = ATAC_countsmatrix_cleaned, coordinate.key = AllPeaks.granges, coordinates){

  mydata.subset <- mydata[subjectHits(findOverlaps(extractCoord(coordinates), coordinate.key)),]

  mydata.subset$PeakID <- row.names(mydata.subset)

  mydata.subset <-  pivot_longer(mydata.subset, cols = colnames(mydata))

  mydata.subset <- cbind(mydata.subset, str_split(mydata.subset$name,pattern="_", simplify = T))
  colnames(mydata.subset)[4:5] <- c("Biopsy","Hours")

  return(mydata.subset)

}


# function for plotting ATAC peaks based on set of valid genomic coordinates


plotATAC.FUN <- function(mydata = ATAC_countsmatrix_cleaned, coordinate.key = AllPeaks.granges, coordinates){

  plot.data <- subset.ATAC.FUN(mydata = mydata, coordinate.key = coordinate.key, coordinates = coordinates)

  p1 <- ggplot(plot.data)+
    stat_summary(aes(x=as.numeric(Hours),y=value,colour=PeakID), geom = "line", fun.y = mean, size = 2)+
    theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=10), title= element_text(size=20) )+
    labs(y="Counts", x="Hours")

  print(p1)

}


