# load libraries

library(tidyverse)
library(ggplot2)
library(GenomicRanges)


# load data

AH_EL_RNA_ALLREPS <- readRDS("./Data/RNA_timecourse_TPM.rds")

Gene_key_hg38 <- read.delim("./Data/GRCh38_key_230713.txt")
Gene_key_hg38 <- subset(Gene_key_hg38, !is.na(GeneName))

geneKey.ranges <- GRanges(seq=Gene_key_hg38$Chromosome,IRanges(start=Gene_key_hg38$Start, end=Gene_key_hg38$End), strand=Gene_key_hg38$Strand,mcols=Gene_key_hg38[,c("EnsemblID","GeneName")])


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

  mydata.subset <-  pivot_longer(mydata.subset, cols = colnames(mydata)[2:41])

  mydata.subset <- cbind(mydata.subset,str_split(mydata.subset$name,pattern="_", simplify = T))
  colnames(mydata.subset)[3:6] <- c("Sample","TPM","Biopsy","Hours")

  return(mydata.subset)

}

# function for calculating zscores

zscore.FUN <- function(mydata){

  biopsies <- unique(mydata$Biopsy)
  genes <- unique(mydata$GeneName)

  mydata$zscore <- NA

  for(i in 1:length(biopsies)){
    for(j in 1:length(genes)){

      mydata$zscore[mydata$GeneName == genes[j] & mydata$Biopsy == biopsies[i]] <- scale(mydata$TPM[mydata$GeneName == genes[j] & mydata$Biopsy == biopsies[i]])

    }
  }

  mydata

}


plotRNA.FUN <- function(mydata, geneName=NULL, ensemblID=NULL, biopsies = c("S169","S170","S506","S508"), times = c(0,3,6,9,12,18,24,36,48,96), yaxis= "TPM"){

  if(is.null(geneName) & is.null(ensemblID)){
    return(
      print(
        ggplot()+theme_bw()
      )
    )
  }

  data.sub <- arrange.FUN(mydata=mydata, geneName=geneName, ensemblID=ensemblID)

  data.sub <- subset(data.sub, Biopsy %in% biopsies & Hours %in% times)

  # TPM as y-axis

  if(yaxis == "TPM" ){

    if(length(geneName) == 1){

      p1 <- ggplot(subset(data.sub))+
        geom_line(aes(x=as.numeric(Hours),y=TPM,colour=Biopsy), size = 2)+
        stat_summary(aes(x=as.numeric(Hours),y=TPM), geom = "line", fun.y = mean, size = 2)+
        theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
        scale_x_continuous(breaks=times, minor_breaks = NULL)+
        labs(y="Transcripts per Million", x="Hours")+ggtitle("RNA-seq")

    }

    if(length(geneName) >= 2){

      p1 <- ggplot(subset(data.sub))+
        stat_summary(aes(x=as.numeric(Hours),y=TPM,colour=GeneName), geom = "line", fun.y = mean, size = 2)+
        theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
        scale_x_continuous(breaks=times, minor_breaks = NULL)+
        labs(y="Transcripts per Million", x="Hours")+ggtitle("RNA-seq")

    }
  }

  # log10 as y-axis

  if(yaxis == "log"){

    if(length(geneName) == 1){

      p1 <- ggplot(subset(data.sub))+
        geom_line(aes(x=as.numeric(Hours),y=TPM,colour=Biopsy), size = 2)+
        stat_summary(aes(x=as.numeric(Hours),y=TPM), geom = "line", fun.y = mean, size = 2)+
        theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
        scale_x_continuous(breaks=times, minor_breaks = NULL)+
        scale_y_continuous(trans='log10')+
        labs(y="Transcripts per Million", x="Hours")+ggtitle("RNA-seq")

    }

    if(length(geneName) >= 2){

      p1 <- ggplot(subset(data.sub))+
        stat_summary(aes(x=as.numeric(Hours),y=TPM,colour=GeneName), geom = "line", fun.y = mean, size = 2)+
        theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
        scale_x_continuous(breaks=times, minor_breaks = NULL)+
        scale_y_continuous(trans='log10')+
        labs(y="Transcripts per Million", x="Hours")+ggtitle("RNA-seq")

    }
  }

  # zscore as y-axis

  if(yaxis == "zscore"){

    data.sub <- zscore.FUN(data.sub)

    if(length(geneName) == 1){

      p1 <- ggplot(subset(data.sub))+
        geom_line(aes(x=as.numeric(Hours),y=zscore,colour=Biopsy), size = 2)+
        stat_summary(aes(x=as.numeric(Hours),y=zscore), geom = "line", fun.y = mean, size = 2)+
        theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
        scale_x_continuous(breaks=times, minor_breaks = NULL)+
        labs(y="Transcripts per Million (normalised)", x="Hours")+ggtitle("RNA-seq")

    }

    if(length(geneName) >= 2){

      p1 <- ggplot(subset(data.sub))+
        stat_summary(aes(x=as.numeric(Hours),y=zscore,colour=GeneName), geom = "line", fun.y = mean, size = 2)+
        theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
        scale_x_continuous(breaks=times, minor_breaks = NULL)+
        labs(y="Transcripts per Million (normalised)", x="Hours")+ggtitle("RNA-seq")

    }
  }


  print(p1)
}



#load ATAC count matrix and genomic coordinates

ATAC_peaks <- read.delim("./Data/ATAC_peaks.bed", header=FALSE)

AllPeaks.granges <- GRanges(seqnames = ATAC_peaks[,1], ranges = IRanges(start=ATAC_peaks[,2], end=ATAC_peaks[,3]), mcols = data.frame(PeakID = ATAC_peaks[,4]))

ATAC_TPMmatrix <- readRDS("./Data/ATAC_timecourse_TPM.rds")

# function for extracting genomic coordinates from text input and generating genomic ranges object

extractCoord <- function(string.input){

  coordinates <- str_split_1(string.input, pattern=regex("[:,_-[:space:]]"))

  if (length(coordinates) != 3 || any(is.na(as.numeric(coordinates[2:3])))) {
    stop("Invalid genomic region provided. Please provide a valid genomic region in the format 'chr:start-end'.")
  }

  return(GRanges(seqnames = coordinates[1], ranges = IRanges(start=as.numeric(coordinates[2]), end=as.numeric(coordinates[3]))))


}

# extend coordinates

grexpand <- function(x, upstream=0, downstream=0){
  if (any(strand(x) == "*")){
    warning("'*' ranges were treated as '+'")}
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

grcontract <- function(x, gene.region){
  if (any(strand(x) == "*")){
    warning("'*' ranges were treated as '+'")}
  on_plus <- strand(x) == "+" | strand(x) == "*"

  if(gene.region == "TSS"){
    new_start <- ifelse(on_plus, start(x), end(x)-1)
    new_end <- ifelse(on_plus,start(x)+1, end(x))
  }
  if(gene.region == "TTS"){

    new_start <- ifelse(on_plus, end(x)-1, start(x))
    new_end <- ifelse(on_plus, end(x), start(x)+1)

  }
  if(gene.region == "Whole"){
    new_start <- start(x)
    new_end <- end(x)
  }

  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}


# find coordinates based on gene names

extractCoordfromGene <- function(genes=NULL, geneCoordKey=geneKey.ranges, add.upstream = 0, add.downstream = 0, gene.region = "Whole"){

  coordinates <- subset(geneKey.ranges, mcols.GeneName %in% genes)
  coordinates <- grcontract(coordinates, gene.region = gene.region)

  coordinates <- grexpand(coordinates, upstream = add.upstream, downstream= add.downstream)

  return(coordinates)

}


# function for subsetting data matrix given a set of coordinates

subset.ATAC.FUN <- function(mydata = ATAC_TPMmatrix, coordinate.key = AllPeaks.granges, coordinates){

  if (length(coordinates) == 0 || any(is.na(coordinates$start)) || any(is.na(coordinates$end))) {
    stop("Invalid gene name or genomic location provided. Please provide a valid gene name or genomic location.")
  }

  mydata.subset <- mydata[subjectHits(findOverlaps(coordinates, coordinate.key)),]

  mydata.subset$PeakID <- row.names(mydata.subset)

  mydata.subset <-  pivot_longer(mydata.subset, cols = colnames(mydata))

  mydata.subset <- cbind(mydata.subset, str_split(mydata.subset$name,pattern="_", simplify = T))
  colnames(mydata.subset)[4:5] <- c("Biopsy","Hours")

  return(mydata.subset)

}


# function for plotting ATAC peaks based on set of valid genomic coordinates


plotATAC.FUN <- function(mydata = ATAC_TPMmatrix, coordinate.key = AllPeaks.granges, coordinates, times = c(0,3,6,9,12,18,24,36,48,96)){

  if (length(coordinates) == 0 || any(is.na(coordinates))) {
    stop("Invalid gene name or genomic location provided. Please provide a valid gene name or genomic location.")
  }

  plot.data <- subset.ATAC.FUN(mydata = mydata, coordinate.key = coordinate.key, coordinates = coordinates)

  p1 <- ggplot(plot.data)+
    stat_summary(aes(x=as.numeric(Hours),y=value,colour=PeakID), geom = "line", fun.y = mean, size = 2)+
    scale_x_continuous(breaks=times, minor_breaks = NULL)+
    theme_bw()+theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text= element_text(size=15), title= element_text(size=20) )+
    labs(y="RPM", x="Hours")+ggtitle("ATAC-seq")

  print(p1)

}


plot_error_message <- function(message) {
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", main = "")
  text(x = 0.5, y = 0.5, label = message, cex = 1.5, col = "red", font = 2)
}




