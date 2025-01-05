#!/usr/bin/env Rscript
#06.05.2024
## script to create Manhattan plots for REGENIE results
## Aakash Nepal
rm(list=ls())

## get the arguments from the command line
#args <- commandArgs(trailingOnly=T)

## little options
#options(stringsAsFactors = F)

setwd("/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs")

## --> import parameters <-- ##
library(data.table)
#outt<- "Lymphocytes"
# Read REGENIE files and filter
#files <- c(list.files(pattern = paste0(outt,"\\.regenie.gz")))
#files<- "GWAS_chr1_height.regenie.gz"

## import labels to add names to the plot 
lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
## add identifier to match with data set
outt<- "Height"
for(outt in lab.phe$V1){
  lab.phe$id <- paste0("manplot_", outt) 
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs__ukb_md.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_ukb.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_ukb.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  ## Plotting the zoomed-in Manhattan plot for chromosome 9
  zoomed_res <- res.man[CHROM == 10 & GENPOS >= 101136163 & GENPOS <= 106111147]
  #top_variant <- zoomed_res[which.max(LOG10P)]
  top_variants <- zoomed_res[order(-LOG10P)][1:10]
  head(res.man)
  png(paste0("./00_manhattan_plots/", outt,"_ukb.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(zoomed_res)
  ## add title
  mtext(paste(outt,"_ukb"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
  
}

for(outt in lab.phe$V1){
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs__ukb_md.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_ukb.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_ukb.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0("./00_manhattan_plots/", outt,"_ukb.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_ukb"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}



#####
outt<- "Lymphocytes"
setwd("/Aakash/11_GWAS_EHR_W_H/output_step2")

## --> import parameters <-- ##
#library(data.table)
#outt<- "MCV"
# Read REGENIE files and filter
#files <- c(list.files(pattern = paste0(outt,"\\.regenie.gz")))
#files<- "GWAS_chr1_height.regenie.gz"

for(outt in lab.phe$V1){
  ## import labels to add names to the plot 
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_md.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_md.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_md.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0("./00_manhattan_plots/", outt,"_md.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_md"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()

}


for(outt in lab.phe$V1){
  ##################################################################
  lab.phe    <- data.table::fread("Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_max.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_max.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_max.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0("./00_manhattan_plots/", outt,"_max.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_max"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}

for(outt in lab.phe$V1){
  #################################
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_min.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_min.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_min.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0("./00_manhattan_plots/", outt,"_min.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_min"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}
#######################################################

for(outt in lab.phe$V1){
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_mean.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_mean.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_mean.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0( "./00_manhattan_plots/",outt,"_mean.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_mean"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}
##################################
for(outt in lab.phe$V1){
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_median.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_median.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_median.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0( "./00_manhattan_plots/",outt,"_median.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_median"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}
##############################
#07.05.2024

setwd("/Aakash/12_GWAS_EHR_UKB/output/20_combined_allchrs_UKB_only")

## import labels to add names to the plot 
lab.phe    <- data.table::fread("Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
outt<- "Height"
## add identifier to match with data set
for(outt in lab.phe$V1){
  lab.phe$id <- paste0("manplot_", outt) 
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0("./00_manhattan_plots/", outt,"_only_ukb.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  ## Plotting the zoomed-in Manhattan plot for chromosome 9
  zoomed_res <- res.man[CHROM == 10 & GENPOS >= 101136163 & GENPOS <= 106111147]
  top_variants_o_ukb <- zoomed_res[order(-LOG10P)][1:10]
  ## Manhattan plot
  plot.manhattan(zoomed_res)
  
  ## add title
  mtext(paste(outt,"_only_ukb"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
  
}

########################
#08.05.2024

setwd("Aakash/12_GWAS_EHR_UKB/output/23_combined_allchrs_EHR_only_qt")

## import labels to add names to the plot 
lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)

for(outt in lab.phe$V1){
  ##################################################################
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_max.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_max.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_max.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0("./00_manhattan_plots/", outt,"_max_only_EHR_qt.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_max_only_EHR"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}

for(outt in lab.phe$V1){
  #################################
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_min.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_min.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_min.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0("./00_manhattan_plots/", outt,"_min_only_EHR.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_min_only_EHR"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}
#######################################################

for(outt in lab.phe$V1){
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_mean.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_mean.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_mean.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0( "./00_manhattan_plots/",outt,"_mean_only_EHR.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_mean_only_EHR"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}
##################################
for(outt in lab.phe$V1){
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_median.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_median.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_median.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0( "./00_manhattan_plots/",outt,"_median_only_EHR.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_median_only_EHR"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}


####################################################################################
setwd("/Aakash/12_GWAS_EHR_UKB/output/19_combined_allchrs_all_pheno")
outt<- lab.phe$V1[1]
for(outt in lab.phe$V1){
  lab.phe    <- data.table::fread("/Aakash/11_GWAS_EHR_W_H/scripts/20_names.txt",header= F)
  ## add identifier to match with data set
  lab.phe$id <- paste0("manplot_", outt) 
  
  #---------------------------#
  ##--    Manhattan plot   --##
  #---------------------------#
  
  ## import results
  res.man    <- paste0(outt,"_allchrs_.tsv.gz")               
  res.man    <- data.table::fread(res.man, fill=T)
  
  ## load regional results to indicate number of hits
  if(file.size(paste0( outt, "_regional_sentinels_ukb.txt")) > 0){
    res.region <- read.table(paste0( outt, "_regional_sentinels_ukb.txt"), header=F, sep="")
  }else{
    res.region <- array(data=NA, dim=c(0,5))
  }
  
  ## --> Manhattan plot <-- ##
  source("/Aakash/11_GWAS_EHR_W_H/scripts/plot_manhattan.R")
  
  png(paste0( "./00_manhattan_plots/",outt,"_ind_all.png") , width = 16, height = 6, res=300, units = "cm")
  par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
  
  ## Manhattan plot
  plot.manhattan(res.man)
  
  ## add title
  mtext(paste(outt,"_all_EHR"),
        side = 3, cex=.5, font=1)
  
  ## add number of hits
  legend("topright", bty="n", cex=.5, legend = paste0("Regional sentinels: ", nrow(res.region)))
  dev.off()
}
