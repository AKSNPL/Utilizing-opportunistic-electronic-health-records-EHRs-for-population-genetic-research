################################################
## function to compile overlap with GWAS catalog
## at different ld thresholds

map.gwas.catalog <- function(gwas.catalogue, res.var, ld.proxies, ld.thr=.8){
  
  ## 'gwas.catalogue' -- data from the GWAS catalog
  ## 'res.var'        -- variants to annotate
  ## 'ld.proxies'     -- table of LD proxies for variants to annotate
  ## 'ld.thr'         -- LD threshold to be used
  
  ## subset proxies to required threshold
  ld.proxies <- subset(ld.proxies, r2 >= ld.thr)
  
  ## do in parallel
  require(doMC)
  registerDoMC(12)
  
  ## go through each variant and 1) identify all proxies, 2) map to GWAS catalog findings, and 3) reduce redundancy
  res.gwas       <- mclapply(1:nrow(res.var), function(x){
    
    ## get the variant of interest
    tmp              <- unlist(res.var[x, "MarkerName"])
    
    ## get all possible proxies
    tmp              <- ld.proxies[MarkerName.1 == tmp | MarkerName.1 == tmp]
    ## proceed only, if any proxies
    if(nrow(tmp) > 0){
      snp              <- unique(unlist(tmp[, c("MarkerName.1", "MarkerName.2"), with=T]))
    }else{
      snp              <- unlist(res.var[x, "MarkerName"])
      tmp              <- data.frame(MarkerName.1=unlist(res.var[x, "MarkerName"]), MarkerName.2=unlist(res.var[x, "MarkerName"]), r2=1)
    }
    ## get mapping with GChr38 coordinates
    snp              <- ukbb.snps[MarkerName %in% snp]
    ## combine with LD information
    snp$MarkerName.1 <- unlist(res.var[x, "MarkerName"])
    ## adapt names
    names(snp)       <- gsub("MarkerName$", "MarkerName.2", names(snp))
    ## add LD
    snp              <- merge(snp, tmp, by=c("MarkerName.1", "MarkerName.2"))
    ## create snp id to optimize merging multiple mappings
    snp$snp.id       <- paste0(snp$chromosome, ":", snp$pos.hg38)
    
    ## create two versions of mapping
    snp.rsid         <- merge(snp, gwas.catalogue, by.x="rsid", by.y="SNPS")
    snp.pos          <- merge(snp, gwas.catalogue, by = "snp.id")
    ## edit
    snp.rsid$snp.id  <- snp.rsid$snp.id.x
    snp.rsid$snp.id.x <- snp.rsid$snp.id.y <- NULL
    snp.pos$SNPS     <- snp.pos$rsid
    snp.rsid$SNPS    <- snp.rsid$rsid
    ## combine
    snp              <- unique(rbind(snp.rsid, snp.pos))
    
    ## prepare return
    if(nrow(snp) > 0){
      ## sort 
      snp              <- snp[order(TRAIT, MarkerName.1, -r2)]
      ## create indicator
      snp[, ind := 1:.N, by=c("TRAIT", "MarkerName.1")]
      ## keep only one finding per trait
      snp              <- snp[ind == 1]
      
      ## report summary back
      snp              <- data.table(rsid.gwas=paste(sort(unique(snp$rsid)), collapse = "||"),
                                     trait_reported=paste(sort(unique(snp$TRAIT)), collapse = "||"),
                                     mapped_trait=paste(sort(unique(snp$MAPPED_TRAIT)), collapse = "||"),
                                     study_id=paste(sort(unique(snp$`STUDY ACCESSION`)), collapse = "||"),
                                     source_gwas=paste(sort(unique(snp$PUBMEDID)), collapse = "||"),
                                     num_reported=nrow(snp))
      
    }else{
      
      ## report summary back
      snp              <- data.table(rsid.gwas="",
                                     trait_reported="",
                                     mapped_trait="",
                                     study_id="",
                                     source_gwas="",
                                     num_reported=0)
    }
    
    ## return data set
    return(data.table(res.var[x,], snp))
  }, mc.cores=12)
  ## combine 
  res.gwas   <- do.call(rbind, res.gwas)
  
  ## return results
  return(res.gwas)
  
}