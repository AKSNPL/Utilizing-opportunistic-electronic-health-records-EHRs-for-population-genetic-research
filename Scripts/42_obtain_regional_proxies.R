#!/usr/bin/env Rscript

## script to extract regional signals in LD with given lead variants

rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

## set path for computing
setwd("/Aakash/12_GWAS_EHR_UKB")

## --> import parameters <-- ##

## get regional coordinates
chr.s <- as.numeric(args[1])
pos.s <- as.numeric(args[2])
pos.e <- as.numeric(args[3])
var.i <- strsplit(args[4], "\\|")[[1]]

## get regional coordinates
#tmp   <- read.table("input/region.query.txt")
#j     <- 8
#chr.s <- as.numeric(tmp$V1[j])
#pos.s <- as.numeric(tmp$V2[j])
#pos.e <- as.numeric(tmp$V3[j])
#var.i <- strsplit(tmp$V4[j], "\\|")[[1]]

#----------------------------#
##--  obtain SNP dosages  --##
#----------------------------#

## package for data loading
require(data.table)

## obtain SNP list from GWAS stats
snp.list <- fread(cmd=paste0("zcat input/04_others/allchrs_RSIDS.tsv.gz |", 
                             " awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, " '{if(($1 == chr && $2 >= low && $2 <= upp) || NR == 1) print $0}'"))

## write list of SNPs to be queried to file 
write.table(snp.list$ID, paste("tmpdir/snplist", "bla_ehr_only_qt", chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)

## import function to do so
source("scripts/14_get_LD.R")

## import 
foo              <- get.LD(chr.s, pos.s, pos.e, "bla_ehr_only_qt")
## separate out into two data sets to ease downstream computation
snp.dat          <- foo[[1]]
snp.info         <- foo[[2]]
## delete and clear up
rm(foo); gc()

#----------------------------#
##--      compute LD      --##
#----------------------------#

## set the names of variants to be queried
jj            <- snp.info$id[which(snp.info$MarkerName %in% var.i)]
## compute LD for selected markers
#gg<- t(cor(snp.dat[, "X1"], snp.dat[, "X6"])^2)
tmp           <- t(cor(snp.dat[, jj], snp.dat[, snp.info$id])^2)
colnames(tmp) <- jj
## convert to data frame
tmp           <- data.frame(id=rownames(tmp), tmp)
## write to file apply threshold
tmp$max.ld    <- apply(tmp[,-1, drop=F], 1, max)
## add MarkerName
tmp           <- merge(snp.info[, c("id", "MarkerName")], tmp)
## subset and store
write.table(subset(tmp, max.ld >= .1), paste("output/09_proxies/bla_ehr_only_qt_LD.proxies", chr.s, pos.s, pos.e, sep="."), row.names=F, sep="\t")

#----------------------------#
##--       cleaning       --##
#----------------------------#

system(paste("rm tmpdir/tmp.bla_ehr_only_qt", chr.s, pos.s, pos.e, "dosage", sep="."))

