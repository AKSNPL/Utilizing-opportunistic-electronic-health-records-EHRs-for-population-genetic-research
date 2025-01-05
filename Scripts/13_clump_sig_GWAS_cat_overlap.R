rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("/Aakash/12_GWAS_EHR_UKB")
for (i in c(
  "data.table",
  "lubridate",
  "tidyr",
  "ggplot2",
  "dplyr",
  "stringr",
  "arrow",
  "readxl",
  "Hmisc",
  "tidyverse",
  "here"
)
) {
  suppressPackageStartupMessages(
    library(i, character.only = TRUE
    ))
}

## set root and parallel settings ####
# Knitr should use the project root and not the script location as root
# base.dir refers to the plot location, that should remain with the script
knitr::opts_knit$set(
  root.dir = here()
)

# Give data.table enough threads
writeLines(paste0("Threads available: ", parallel::detectCores()))
writeLines(paste0("Threads given to data.table: ", parallel::detectCores() / 2))
setDTthreads(parallel::detectCores() / 2)

# load the required data ####

#############################
#### import lead signals ####
#############################

## import regional lead signals
res.region <- fread("output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_md.txt")
res.region_ukb <- fread("output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_ukb.txt")
res.region_min <- fread("output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_min.txt")
res.region_max <- fread("output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_max.txt")
res.region_mean <- fread("output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_mean.txt")
res.region_median <- fread("output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_median.txt")
res.region_ukb$id <- paste(res.region_ukb$id,"ukb",sep= "_")
res.region_min$id <- paste(res.region_min$id,"min",sep= "_")
res.region_mean$id <- paste(res.region_mean$id,"mean",sep= "_")
res.region_max$id <- paste(res.region_max$id,"max",sep= "_")
res.region_median$id <- paste(res.region_median$id,"median",sep= "_")
res.region$id <- paste(res.region$id,"md",sep= "_")
res.region_merged<- rbind(res.region,res.region_ukb,res.region_min,res.region_mean,res.region_max,res.region_median)
#res.region<- res.region_merged%>%
#  filter(id=="Weight_md" | id == "Weight_ukb")
####set to clump from
#test<- fread("../input/Results.retinal.risk.states.GWAS.20240105_md.txt")
res.region<- res.region_merged
## create MarkerName
res.region[, MarkerName := paste0("chr", CHROM, ":", GENPOS, "_", pmin(as.character(ALLELE0), as.character(ALLELE1)), "_", pmax(as.character(ALLELE0), as.character(ALLELE1)))]

## sort
res.region <- res.region[order(CHROM, region_start)]

## from https://stackoverflow.com/questions/28938147/how-to-flatten-merge-overlapping-time-periods

tmp        <- res.region %>%
  ## sort
  arrange(CHROM, region_start) %>% 
  group_by(CHROM) %>%
  ## create index for regions
  mutate(indx = c(0, cumsum(as.numeric(lead(region_start)) >
                              cummax(as.numeric(region_end)))[-n()])) %>%
  group_by(CHROM, indx) %>%
  summarise(region_start = data.table::first(region_start), region_end = data.table::last(region_end)) %>% as.data.table()
## 

### the code down should pull all markers (SNPs with chr_pos_allele notations) per respective new non-overlapping regions

## set keys to indicate on what to define new regions on
setkey(tmp, CHROM, region_start, region_end)

## add information to regional results-> Identifying Non-overlapping Regions
foo        <- foverlaps(res.region, tmp, by.x=c("CHROM", "region_start", "region_end"))
## the prefix 'i' denotes the former regional boundaries before collapsing

## create new index for regions
foo$group  <- paste(foo$CHROM, foo$indx, sep="_")
## assign
res.region <- foo

## create data frame to query each region
region     <- lapply(unique(res.region$group), function(x){
  
  ## get the group of interest
  tmp <- res.region[ group == x, c("MarkerName", "CHROM", "region_start", "region_end")]
  ## return subset
  return(data.table(tmp[1, c("CHROM", "region_start", "region_end")], MarkerName = paste(unique(tmp$MarkerName), collapse = "|")))
  
})
## combine
region     <- do.call(rbind, region)


# DONE
## write to file to extract proxies
#write.table(region, "output/08_res_sentinels_all/region.query_all.txt", sep="\t", col.names = F, row.names = F, quote = F)


#Load the saved file
#region <- read.table("region.query.txt", header=T)
#
######################################
####       import proxy SNPs      ####
######################################

## get all proxy files (MHC might be missing)
jj <- dir(pattern="*all_LD.pr*","./output/09_proxies/")


## collate
ld.proxies <- lapply(jj, function(x){
  
  ## import LD-matrix
  tmp           <- fread(paste0("./output/09_proxies/", x), data.table = F)
  ## ensure no duplications
  ii            <- table(tmp$MarkerName)
  tmp           <- tmp %>% filter(MarkerName %in% names(ii[ii == 1]))
  ## assign rownames
  rownames(tmp) <- tmp$MarkerName
  
  ## get all SNPs
  jj  <- names(tmp) # column names of the ld.proxie files
  jj  <- jj[-which(jj %in% c("id", "MarkerName", "max.ld"))] #do not select those that are other than  c() 
  
  ## convert to table for each SNP of interest (MakerName by MarkerName)
  foo <- lapply(jj, function(c){
    
    ## get MarkerName of interest
    kk              <- tmp[which(tmp[, c] >= .1), c, drop=F]
    ## get SNPs of interest
    kk              <- data.frame(MarkerName.1=tmp$MarkerName[which(tmp$id == c)], 
                                  MarkerName.2=rownames(kk),
                                  r2=kk[,1])
    return(kk)
    
  })
  return(do.call(rbind, foo))
  
})
## combine into one large table
ld.proxies <- do.call(rbind, ld.proxies)
## convert
ld.proxies <- as.data.table(ld.proxies)

######################################
####     LD-clumping of signals   ####
######################################

## do as a graph
require(igraph)

cx<- ld.proxies[ r2 >= .7]
## convert to graph (use max to allow for edge weights)
ld.sub     <- graph_from_data_frame(ld.proxies[ r2 >= .7])

ld.proxie<- ld.proxies %>%
  filter(MarkerName.1== "rs78339356")
head(ld.proxies)
## get all separate components
ld.sub     <- components(ld.sub)$membership
## convert to data frame

ld.sub     <- data.table(ID=names(ld.sub), R2.group=ld.sub)

## add to results 
res.region <- merge(res.region, ld.sub, by.x = "MarkerName", by.y = "ID", all.x=T)

table(res.region$R2.group, useNA = "always")
length(unique(res.region$R2.group))
length(unique(res.region$MarkerName))
summary(res.region$R2.group)

##  
summary(res.region$id)

aa<- res.region %>%
  filter(id=="Weight_ukb") %>%
  drop_na()

ab<- res.region %>%
  filter(id=="Weight_md") %>%
  drop_na()
#summary(res.region$R2.group)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#      1    1042    2162    2142    3197    4267     236 

###############Visualize the signals########



#remove NAs
#res.region <- res.region[complete.cases(res.region[, "R2.group"]), ]

#===============================SAVING
#SAVE RES.RESGION AT THIS STAGE - 
#saveRDS(res.region, "res.region_for_matrix.rds")
#save res.region .. 
#saveRDS(res.region, "res.region_ind_for_matrix.rds")





###############################################
####          map to build 38              ####
###############################################

#-------------------------------------#
##--       get all proxies         --##
#-------------------------------------#

## get unique list of variants
snp.lift   <- data.table(MarkerName=unique(c(ld.proxies$MarkerName.1, ld.proxies$MarkerName.2, res.region$MarkerName)))
dim(snp.lift) #

## add rsIDs #155k snps
ukbb.snps  <- fread("/14_phecode_GWAS/04_V2G/01_VEP_annotation/data/QCed.UKBB.variants.Berlin.20220404.txt")
dim(ukbb.snps)

## I will need to create a new ukbb.snps_new?
#using this script
#/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/04_V2G/01_VEP_annotation/scripts

ukbb.snps  <- ukbb.snps[MarkerName %in% snp.lift$MarkerName]
dim(ukbb.snps)#1 023 314 seem to have rsid and the rest not. Is it because 
# not all snps are addressed in that 04_V2G table or because simply theres not rsid annotation available
#for 2k SNPs ?

## I need to see this for my data wrt old comment "some are missing, and some are twice (multi-allelic)"


#table(duplicated(ukbb.snps$rsid))

# FALSE    TRUE 
# 1022830     484

#-------------------------------------#
##--      lift over to build 38    --##
#-------------------------------------#

## R package to do so
require(liftOver)

## edit chromosome
ukbb.snps$chr.hg39 <- ifelse(ukbb.snps$chromosome == 23, "X", ukbb.snps$chromosome)

## import chain
chain              <- import.chain("/14_phecode_GWAS/04_V2G/03_eQTL/input/hg19ToHg38.over.chain")

## create GRanges object
grObject           <- GRanges(seqnames = paste0("chr", ukbb.snps$chr.hg39), ranges=IRanges(start = ukbb.snps$position, end = ukbb.snps$position, names=ukbb.snps$MarkerName))

## now do the mapping
tmp                <- as.data.frame(liftOver(grObject, chain))
## rename
names(tmp)         <- c("group", "MarkerName", "seqnames", "pos.hg38", "end", "width", "strand")
## add to the data
ukbb.snps          <- merge(ukbb.snps, tmp[, c("MarkerName", "pos.hg38")], all.x=T)
## some will otherwise be lost

######################################
####     overlap GWAS catalog     ####
######################################

## import download from 31/10/2023
gwas.catalogue           <- fread("/23_GWAS_retinal_risk_states/02_clump_regional_sentinels/input/alternative", sep="\t", header=T)

## rename some columns
names(gwas.catalogue)[8] <- "TRAIT"

## prune GWAS catalogue data
gwas.catalogue           <- subset(gwas.catalogue, !is.na(`OR or BETA`) & is.finite(`OR or BETA`))
## generate risk allele and drop everything w/o this information
gwas.catalogue$riskA     <- sapply(gwas.catalogue$`STRONGEST SNP-RISK ALLELE`, function(x) strsplit(x,"-")[[1]][2])
## delete spaces at the end
gwas.catalogue$riskA     <- trimws(gwas.catalogue$riskA, which = "b")
## drop interaction entries
ii                       <- grep("[0-9]", gwas.catalogue$riskA)
gwas.catalogue           <- gwas.catalogue[-ii,]
## only genome-wide significant ones
gwas.catalogue           <- subset(gwas.catalogue, PVALUE_MLOG > 7.3)
## N = 381,006 entries

## create another entry to possible merge on (careful, genome build 38 mapping)
gwas.catalogue[, snp.id := paste0(ifelse(CHR_ID == "X", 23, CHR_ID), ":", CHR_POS)]

###############################################
####  assign each variant to consequences  ####
###############################################

#-------------------------------------#
##--        create mapping         --##
#-------------------------------------#

## get all unique variants from regional results
res.var    <- unique(res.region[, c("MarkerName", "ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1")])
## n 

## load function to do so
source("/23_GWAS_retinal_risk_states/02_clump_regional_sentinels/scripts/map_gwas_catalog.R")

## strong LD
res.gwas.8 <- map.gwas.catalog(gwas.catalogue, res.var, ld.proxies, ld.thr = .8)
## moderate LD
res.gwas.1 <- map.gwas.catalog(gwas.catalogue, res.var, ld.proxies, ld.thr = .1)

## combine with regional results
res.region <- merge(res.region, res.gwas.8, by = c("MarkerName", "ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1")) 
res.region <- merge(res.region, res.gwas.1, by = c("MarkerName", "ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1"), suffixes=c(".r2.8", ".r2.1")) 

#--------------------------------#
##--      gene annotation     --##
#--------------------------------#

## get list of variants to query

## import VEP assignment:
vep.snps   <- fread("/01_VEP/data/VEP_QCed.UKBB.variants.Berlin.20222307_final.txt")
## create MarkerName
vep.snps[, MarkerName := paste0("chr", CHROM, ":", POS, "_", pmin(REF, ALT), "_", pmax(REF, ALT))]
## subset to ld.proxies of interest
jj         <- ld.proxies[ r2 >= .8]
vep.snps   <- vep.snps[ MarkerName %in% unique(c(res.region$MarkerName, jj$MarkerName.1, jj$MarkerName.2))]
## some SNPs are missing...

# just above 126K? how so? bc of filtering or so many missing snps?
length(unique(vep.snps$MarkerName)) #126 692

## add annotation of lead variants
jj.anno    <- merge(jj, vep.snps, by.x="MarkerName.1", by.y="MarkerName", all.x=T)
## add annotation for proxies
jj.anno    <- merge(jj.anno, vep.snps, by.x="MarkerName.2", by.y="MarkerName", all.x=T, suffixes = c(".lead", ".proxy"))

## import consequence table
rank.cons  <- data.frame(consequence=c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion",
                                       "missense_variant", "protein_altering_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
                                       "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
                                       "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
                                       "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation", "regulatory_region_variant", "feature_truncation"), 
                         rank=1:35)

## assign to each proxy variant the most consequential score
jj.anno[, rank.proxy := sapply(Consequence.proxy, function(x){
  ## split multiple assignments
  x <- strsplit(x, "&")[[1]]
  # report most severe
  x <- subset(rank.cons, consequence %in% x)
  if(nrow(x) > 0){
    return(min(x$rank))
  }else{
    return(36)
  }
})]

## order
jj.anno    <- jj.anno[ order(MarkerName.1, rank.proxy)]
## keep only the most severe
jj.anno[, ind := 1:.N, by="MarkerName.1"]
jj.anno    <- jj.anno[ ind == 1]
## keep only columns needed
jj.anno    <- jj.anno[, c("MarkerName.1","Consequence.lead", "IMPACT.lead", "SYMBOL.lead", "Gene.lead", "CADD_PHRED.lead",
                          "MarkerName.2", "r2", "POS.proxy", "ID.proxy", "Consequence.proxy", "IMPACT.proxy", "Gene.proxy", 
                          "SYMBOL.proxy", "CADD_PHRED.proxy")]
## add to the results
res.region <- merge(res.region, jj.anno, by.x="MarkerName", by.y="MarkerName.1")
#saveRDS(res.region, "res.region_gene_annot.RDS")
#you can load this and proceed with the downstream steps

#--------------------------------#
##--       closest gene       --##
#--------------------------------#

## add closest gene
require(biomaRt)

## get data on build 37
gene.ensembl                       <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) 

## obtain a list of all protein coding genes
tmp.genes                          <- getBM(attributes = c('chromosome_name', 'start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'transcription_start_site'),
                                            filters = c('chromosome_name'),
                                            values = list(c(1:22, "X")),
                                            mart = gene.ensembl)
## convert to data table to ease query
tmp.genes                          <- as.data.table(tmp.genes)
## convert some entries to numeric
tmp.genes$start_position           <- as.numeric(tmp.genes$start_position)
tmp.genes$end_position             <- as.numeric(tmp.genes$end_position)
tmp.genes$transcription_start_site <- as.numeric(tmp.genes$transcription_start_site)

## loop through each entry and define closest gene (protein coding)
tmp          <- lapply(1:nrow(res.region), function(x){
  
  print(x)
  
  ## get the gene of potential interest (2MB window)
  tmp           <- tmp.genes[ chromosome_name == ifelse(res.region$CHROM[x] == 23, "X", res.region$CHROM[x]) & start_position >= res.region$GENPOS[x]-1e6 & end_position <= res.region$GENPOS[x]+1e6 ]
  
  ## test whether enough protein coding gens are in the region
  if(nrow(subset(tmp, gene_biotype %in% c("protein_coding", "processed_transcript"))) > 0){
    ## restrict to protein encoding genes for now
    tmp <- subset(tmp, gene_biotype %in% c("protein_coding", "processed_transcript"))
  }
  
  ## compute distance to any gene (gene body not TSS)
  tmp$dist.body <- apply(tmp[, c("start_position", "end_position")], 1, function(k) min(abs(k-res.region$GENPOS[x])))
  ## compute distance to TSS
  tmp$dist.tss  <- abs(tmp$transcription_start_site - res.region$GENPOS[x])
  ## sort by distance to gene body
  tmp           <- tmp[order(tmp$dist.body), ]
  ## get position for TSS
  jj            <- which.min(tmp$dist.tss)
  
  ## return three different entries
  return(data.table(res.region[x, ], 
                    closest.gene.body=paste0(tmp$external_gene_name[1], " (", tmp$dist.body[1], ") - ", tmp$ensembl_gene_id[1]),
                    closest.gene.tss=paste0(tmp$external_gene_name[jj], " (", tmp$dist.tss[jj], ") - ", tmp$ensembl_gene_id[jj]),
                    closest.genes=paste(unique(subset(tmp, dist.body <= 5e5)$external_gene_name), collapse = "|")))
  
  
  
})
## combine everything
tmp     <- do.call(rbind, tmp)
#save.image()

## reassign
res.region <- tmp

## write to file
#write.table(res.region, "Results.GWAS_ukb_md_genes20240318.txt", row.names = F, sep="\t")
## write to file
write.table(res.region, "/Aakash/12_GWAS_EHR_UKB/output/10_result_regions/Results.GWAS_all_genes20240407.txt", row.names = F, sep="\t")

#####################################################################################################################
#####################################################################################################################
####                                              END OF SCRIPT                                                  ####
#####################################################################################################################
#####################################################################################################################
