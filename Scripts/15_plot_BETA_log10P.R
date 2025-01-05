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


# list of phenotypes:
phenotypes <- c("ALP", "ALT", "Albumin", "Basophills", "CRP", "Calcium", "Cholesterol", 
               "Creatinine", "DBP", "Eosinophills", "Glucose", "HDL", 
               "Haematocritperc", "Haemoglobinconc", "HbA1c", "Lymphocytes", "MCHbconc", 
               "MCV", "Monocytes", "Neutrophills", "Platelets", "RBC", "SBP", 
               "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides")

#phenotypes<-c("FVC","FEV1")
#name<- "ALT"
## 1)EHR_md_UKB #### ############################################################################
tryCatch(
  expr = { 
    skip_to_next <- FALSE
    for (name in phenotypes){
    cat(name,"\n")
    # load the required data ####
    r2.sentinels<- fread("/Aakash/12_GWAS_EHR_UKB/output/10_result_regions/Results.GWAS_all_genes20240407.txt")
  
    # filter the data for required 
    r2.sentinels<- r2.sentinels %>%
      filter(id==paste0(name,"_md")| id== paste0(name,"_ukb"))
    
    # convert it to a dataframe
    filtered_aa <- as.data.frame(r2.sentinels)
    
    # select the required columns:
    r2.sentinels<- r2.sentinels %>% 
      dplyr::select(BETA,LOG10P,id,ID,R2.group,SE)
    
    # Find duplicated R2 groups
    duplicated_groups <- duplicated(r2.sentinels$R2.group) | duplicated(r2.sentinels$R2.group, fromLast = TRUE)
    
    # Subset the data to keep only rows with duplicated R2 groups
    filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
      filter(id==paste0(name,"_md"))
    
    ## Algorithm ####
    # load the data for the phenotypes
    data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs__ukb_md.tsv.gz"))
    
    # pick up required rows for the selected ids
    ll<- c()
    for(ids in filtered_r2.sentinels$ID){
      data_filt<- c()
      #filter the data:
      data_filt <- data %>%
        filter(ID==ids)
      
      #bind the required data
      ll<- rbind(ll,data_filt)
      
    }
    
    # assign the name
    ll$id<- paste0(name,"_ukb")
    
    # select the required columns
    ll<- ll%>%
      dplyr::select(BETA,LOG10P,id,ID,SE)
    
    # select the required columns
    filtered_r2.sentinels<- filtered_r2.sentinels%>%
      dplyr::select(BETA,LOG10P,id,ID,SE)
    
    # filter the required column
    r2.sentinels1 <- r2.sentinels %>%
      filter(id==paste0(name,"_md"))
    
    mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
    
    mainl1<- rbind(ll,filtered_r2.sentinels)
    mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
    
    mainl0<- mainl2 %>%
      distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
    
    mainl0<- mainl0 %>% group_by(ID) %>%
      mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
      distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
    
    #data1<- data%>%
    #  filter(ID == "rs41280328")
    
    ######################################
    
    filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
      filter(id==paste0(name,"_ukb"))
    
    data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs_md.tsv.gz"))
    ll<- c()
    for(ids in filtered_r2.sentinels$ID){
      #data_filt<- c()
      #filter the data:
      data_filt <- data %>%
        filter(ID==ids)
      ll<- rbind(ll,data_filt)
      
    }
    
    ll$id<- paste0(name,"_md")
    
    ll<- ll%>%
      dplyr::select(BETA,LOG10P,id,ID,SE)
    
    filtered_r2.sentinels<- filtered_r2.sentinels%>%
      dplyr::select(BETA,LOG10P,id,ID,SE)
    
    r2.sentinels1 <- r2.sentinels %>%
      filter(id==paste0(name,"_ukb"))
    
    mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
    
    mainl1<- rbind(ll,filtered_r2.sentinels)
    mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
    
    mainl3<- mainl2 %>%
      distinct(ID,id,LOG10P,BETA,R2.group, .keep_all = TRUE)
    
    mainl3<- mainl3 %>% group_by(ID) %>%
      mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
      distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
    
    main00<- rbind(mainl0, mainl3)
    
    # write the table
    write.table(main00, paste("output/11_BETA_LOG10P_plots/", paste0(name,"_pdata.txt"), sep="\t"))
    
    ## plot ####
    # Convert filtered_aa to data.frame
    library(reshape2)
    library(tidyr)
    library(dplyr)
    
    filt_ukb<- main00%>%
      filter(id==paste0(name,"_ukb"))
    
    filt_md <- main00%>%
      filter(id==paste0(name,"_md"))
    
    
    merged_t<- merge(filt_md,filt_ukb,by="R2.group")
    
    merged_t<- merged_t %>%
      drop_na()
    
    merged_t<- merged_t %>%
      dplyr::select(BETA.x,BETA.y,LOG10P.x,LOG10P.y,ID.x,ID.y,id.x,R2.group)%>%
      mutate(BETA.x = abs(BETA.x))%>%
      mutate(BETA.y= abs(BETA.y))
    
    # Plot BETA.x against BETA.y
    #p1<- plot(merged_t$LOG10P.x, merged_t$LOG10P.y, 
    #     xlab = "LOG10P_EHR", ylab = "LOG10P_UKB")
    
    # Optionally, you can add a regression line
    #p1<- abline(lm(BETA.y ~ BETA.x, data = merged_t), col = "red")
    
    p1 <- ggplot(merged_t, aes(x = LOG10P.x, y = LOG10P.y)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(x = "LOG10P_EHR_md", y = "LOG10P_UKB") +
      facet_wrap(~id.x, scales = "free")
    
    #p1
    # Plot BETA.x against BETA.y
    #p2<- plot(merged_t$BETA.x, merged_t$BETA.y, 
    #     xlab = "BETA_EHR", ylab = "BETA_UKB")
    ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_LOG10P_md_ukb.png"), sep = ""), plot = p1, width = 10, height = 6, dpi = 500)
    
    p2 <- ggplot(merged_t, aes(x = BETA.x, y = BETA.y)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(x = "BETA_EHR_md", y = "BETA_UKB")+
      facet_wrap(~id.x, scales = "free")
    #p2
    ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_BETA_md_ukb.png"), sep = ""), plot = p2, width = 10, height = 6, dpi = 500)

  }
  },
  error = function(e){
    message('Caught an error!')
    skip_to_next <<- TRUE
    #print(e)
  }
)



## 2)EHR_mean_UKB #### ############################################################################
tryCatch(
  expr = { 
    skip_to_next <- FALSE
    for (name in phenotypes){
      cat(name,"\n")
      # load the required data ####
      r2.sentinels<- fread("/Aakash/12_GWAS_EHR_UKB/output/10_result_regions/Results.GWAS_all_genes20240407.txt")
      
      # filter the data for required 
      r2.sentinels<- r2.sentinels %>%
        filter(id==paste0(name,"_mean")| id== paste0(name,"_ukb"))
      
      # convert it to a dataframe
      filtered_aa <- as.data.frame(r2.sentinels)
      
      # select the required columns:
      r2.sentinels<- r2.sentinels %>% 
        dplyr::select(BETA,LOG10P,id,ID,R2.group,SE)
      
      # Find duplicated R2 groups
      duplicated_groups <- duplicated(r2.sentinels$R2.group) | duplicated(r2.sentinels$R2.group, fromLast = TRUE)
      
      # Subset the data to keep only rows with duplicated R2 groups
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_mean"))
      
      ## Algorithm ####
      # load the data for the phenotypes
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs__ukb_md.tsv.gz"))
      
      # pick up required rows for the selected ids
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        
        #bind the required data
        ll<- rbind(ll,data_filt)
        
      }
      
      # assign the name
      ll$id<- paste0(name,"_ukb")
      
      # select the required columns
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # select the required columns
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # filter the required column
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_mean"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl0<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      mainl0<- mainl0 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      #data1<- data%>%
      #  filter(ID == "rs41280328")
      
      ######################################
      
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_ukb"))
      
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs_mean.tsv.gz"))
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        #data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        ll<- rbind(ll,data_filt)
        
      }
      
      ll$id<- paste0(name,"_mean")
      
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_ukb"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl3<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group, .keep_all = TRUE)
      
      mainl3<- mainl3 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      main00<- rbind(mainl0, mainl3)
      
      # write the table
      write.table(main00, paste("output/11_BETA_LOG10P_plots/", paste0(name,"_pdata_mean.txt"), sep="\t"))
      
      ## plot ####
      # Convert filtered_aa to data.frame
      library(reshape2)
      library(tidyr)
      library(dplyr)
      
      filt_ukb<- main00%>%
        filter(id==paste0(name,"_ukb"))
      
      filt_md <- main00%>%
        filter(id==paste0(name,"_mean"))
      
      
      merged_t<- merge(filt_md,filt_ukb,by="R2.group")
      
      merged_t<- merged_t %>%
        drop_na()
      
      merged_t<- merged_t %>%
        dplyr::select(BETA.x,BETA.y,LOG10P.x,LOG10P.y,ID.x,ID.y,id.x,R2.group)%>%
        mutate(BETA.x = abs(BETA.x))%>%
        mutate(BETA.y= abs(BETA.y))
      
      # Plot BETA.x against BETA.y
      #p1<- plot(merged_t$LOG10P.x, merged_t$LOG10P.y, 
      #     xlab = "LOG10P_EHR", ylab = "LOG10P_UKB")
      
      # Optionally, you can add a regression line
      #p1<- abline(lm(BETA.y ~ BETA.x, data = merged_t), col = "red")
      
      p1 <- ggplot(merged_t, aes(x = LOG10P.x, y = LOG10P.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "LOG10P_EHR_mean", y = "LOG10P_UKB")+
        facet_wrap(~id.x, scales = "free")
      
      #p1
      # Plot BETA.x against BETA.y
      #p2<- plot(merged_t$BETA.x, merged_t$BETA.y, 
      #     xlab = "BETA_EHR", ylab = "BETA_UKB")
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_LOG10P_mean_ukb.png"), sep = ""), plot = p1, width = 10, height = 6, dpi = 500)
      
      p2 <- ggplot(merged_t, aes(x = BETA.x, y = BETA.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "BETA_EHR_mean", y = "BETA_UKB")+
        facet_wrap(~id.x, scales = "free")
      #p2
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_BETA_mean_ukb.png"), sep = ""), plot = p2, width = 10, height = 6, dpi = 500)
      
    }
  },
  error = function(e){
    message('Caught an error!')
    skip_to_next <<- TRUE
    #print(e)
  }
)


## 3)EHR_median_UKB #### ############################################################################
tryCatch(
  expr = { 
    skip_to_next <- FALSE
    for (name in phenotypes){
      cat(name,"\n")
      # load the required data ####
      r2.sentinels<- fread("/Aakash/12_GWAS_EHR_UKB/output/10_result_regions/Results.GWAS_all_genes20240407.txt")
      
      # filter the data for required 
      r2.sentinels<- r2.sentinels %>%
        filter(id==paste0(name,"_median")| id== paste0(name,"_ukb"))
      
      # convert it to a dataframe
      filtered_aa <- as.data.frame(r2.sentinels)
      
      # select the required columns:
      r2.sentinels<- r2.sentinels %>% 
        dplyr::select(BETA,LOG10P,id,ID,R2.group,SE)
      
      # Find duplicated R2 groups
      duplicated_groups <- duplicated(r2.sentinels$R2.group) | duplicated(r2.sentinels$R2.group, fromLast = TRUE)
      
      # Subset the data to keep only rows with duplicated R2 groups
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_median"))
      
      ## Algorithm ####
      # load the data for the phenotypes
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs__ukb_md.tsv.gz"))
      
      # pick up required rows for the selected ids
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        
        #bind the required data
        ll<- rbind(ll,data_filt)
        
      }
      
      # assign the name
      ll$id<- paste0(name,"_ukb")
      
      # select the required columns
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # select the required columns
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # filter the required column
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_median"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl0<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      mainl0<- mainl0 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      #data1<- data%>%
      #  filter(ID == "rs41280328")
      
      ######################################
      
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_ukb"))
      
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs_median.tsv.gz"))
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        #data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        ll<- rbind(ll,data_filt)
        
      }
      
      ll$id<- paste0(name,"_median")
      
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_ukb"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl3<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group, .keep_all = TRUE)
      
      mainl3<- mainl3 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      main00<- rbind(mainl0, mainl3)
      
      # write the table
      write.table(main00, paste("output/11_BETA_LOG10P_plots/", paste0(name,"_pdata_median.txt"), sep="\t"))
      
      ## plot ####
      # Convert filtered_aa to data.frame
      library(reshape2)
      library(tidyr)
      library(dplyr)
      
      filt_ukb<- main00%>%
        filter(id==paste0(name,"_ukb"))
      
      filt_md <- main00%>%
        filter(id==paste0(name,"_median"))
      
      
      merged_t<- merge(filt_md,filt_ukb,by="R2.group")
      
      merged_t<- merged_t %>%
        drop_na()
      
      merged_t<- merged_t %>%
        dplyr::select(BETA.x,BETA.y,LOG10P.x,LOG10P.y,ID.x,ID.y,id.x,R2.group)%>%
        mutate(BETA.x = abs(BETA.x))%>%
        mutate(BETA.y= abs(BETA.y))
      
      # Plot BETA.x against BETA.y
      #p1<- plot(merged_t$LOG10P.x, merged_t$LOG10P.y, 
      #     xlab = "LOG10P_EHR", ylab = "LOG10P_UKB")
      
      # Optionally, you can add a regression line
      #p1<- abline(lm(BETA.y ~ BETA.x, data = merged_t), col = "red")
      
      p1 <- ggplot(merged_t, aes(x = LOG10P.x, y = LOG10P.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "LOG10P_EHR_median", y = "LOG10P_UKB")+
        facet_wrap(~id.x, scales = "free")
      
      #p1
      # Plot BETA.x against BETA.y
      #p2<- plot(merged_t$BETA.x, merged_t$BETA.y, 
      #     xlab = "BETA_EHR", ylab = "BETA_UKB")
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_LOG10P_median_ukb.png"), sep = ""), plot = p1, width = 10, height = 6, dpi = 500)
      
      p2 <- ggplot(merged_t, aes(x = BETA.x, y = BETA.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "BETA_EHR_median", y = "BETA_UKB")+
        facet_wrap(~id.x, scales = "free")
      #p2
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_BETA_median_ukb.png"), sep = ""), plot = p2, width = 10, height = 6, dpi = 500)
      
    }
  },
  error = function(e){
    message('Caught an error!')
    skip_to_next <<- TRUE
    #print(e)
  }
)

## 4)EHR_minimum_UKB #### ############################################################################
tryCatch(
  expr = { 
    skip_to_next <- FALSE
    for (name in phenotypes){
      cat(name,"\n")
      # load the required data ####
      r2.sentinels<- fread("/Aakash/12_GWAS_EHR_UKB/output/10_result_regions/Results.GWAS_all_genes20240407.txt")
      
      # filter the data for required 
      r2.sentinels<- r2.sentinels %>%
        filter(id==paste0(name,"_min")| id== paste0(name,"_ukb"))
      
      # convert it to a dataframe
      filtered_aa <- as.data.frame(r2.sentinels)
      
      # select the required columns:
      r2.sentinels<- r2.sentinels %>% 
        dplyr::select(BETA,LOG10P,id,ID,R2.group,SE)
      
      # Find duplicated R2 groups
      duplicated_groups <- duplicated(r2.sentinels$R2.group) | duplicated(r2.sentinels$R2.group, fromLast = TRUE)
      
      # Subset the data to keep only rows with duplicated R2 groups
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_min"))
      
      ## Algorithm ####
      # load the data for the phenotypes
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs__ukb_md.tsv.gz"))
      
      # pick up required rows for the selected ids
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        
        #bind the required data
        ll<- rbind(ll,data_filt)
        
      }
      
      # assign the name
      ll$id<- paste0(name,"_ukb")
      
      # select the required columns
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # select the required columns
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # filter the required column
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_min"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl0<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      mainl0<- mainl0 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      #data1<- data%>%
      #  filter(ID == "rs41280328")
      
      ######################################
      
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_ukb"))
      
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs_min.tsv.gz"))
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        #data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        ll<- rbind(ll,data_filt)
        
      }
      
      ll$id<- paste0(name,"_min")
      
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_ukb"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl3<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group, .keep_all = TRUE)
      
      mainl3<- mainl3 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      main00<- rbind(mainl0, mainl3)
      
      # write the table
      write.table(main00, paste("output/11_BETA_LOG10P_plots/", paste0(name,"_pdata_min.txt"), sep="\t"))
      
      ## plot ####
      # Convert filtered_aa to data.frame
      library(reshape2)
      library(tidyr)
      library(dplyr)
      
      filt_ukb<- main00%>%
        filter(id==paste0(name,"_ukb"))
      
      filt_md <- main00%>%
        filter(id==paste0(name,"_min"))
      
      
      merged_t<- merge(filt_md,filt_ukb,by="R2.group")
      
      merged_t<- merged_t %>%
        drop_na()
      
      merged_t<- merged_t %>%
        dplyr::select(BETA.x,BETA.y,LOG10P.x,LOG10P.y,ID.x,ID.y,id.x,R2.group)%>%
        mutate(BETA.x = abs(BETA.x))%>%
        mutate(BETA.y= abs(BETA.y))
      
      # Plot BETA.x against BETA.y
      #p1<- plot(merged_t$LOG10P.x, merged_t$LOG10P.y, 
      #     xlab = "LOG10P_EHR", ylab = "LOG10P_UKB")
      
      # Optionally, you can add a regression line
      #p1<- abline(lm(BETA.y ~ BETA.x, data = merged_t), col = "red")
      
      p1 <- ggplot(merged_t, aes(x = LOG10P.x, y = LOG10P.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "LOG10P_EHR_min", y = "LOG10P_UKB")+
        facet_wrap(~id.x, scales = "free")
      
      #p1
      # Plot BETA.x against BETA.y
      #p2<- plot(merged_t$BETA.x, merged_t$BETA.y, 
      #     xlab = "BETA_EHR", ylab = "BETA_UKB")
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_LOG10P_min_ukb.png"), sep = ""), plot = p1, width = 10, height = 6, dpi = 500)
      
      p2 <- ggplot(merged_t, aes(x = BETA.x, y = BETA.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "BETA_EHR_min", y = "BETA_UKB")+
        facet_wrap(~id.x, scales = "free")
      #p2
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_BETA_min_ukb.png"), sep = ""), plot = p2, width = 10, height = 6, dpi = 500)
      
    }
  },
  error = function(e){
    message('Caught an error!')
    skip_to_next <<- TRUE
    #print(e)
  }
)


## 5)EHR_maximum_UKB #### ############################################################################
tryCatch(
  expr = { 
    skip_to_next <- FALSE
    for (name in phenotypes){
      cat(name,"\n")
      # load the required data ####
      r2.sentinels<- fread("/Aakash/12_GWAS_EHR_UKB/output/10_result_regions/Results.GWAS_all_genes20240407.txt")
      
      # filter the data for required 
      r2.sentinels<- r2.sentinels %>%
        filter(id==paste0(name,"_max")| id== paste0(name,"_ukb"))
      
      # convert it to a dataframe
      filtered_aa <- as.data.frame(r2.sentinels)
      
      # select the required columns:
      r2.sentinels<- r2.sentinels %>% 
        dplyr::select(BETA,LOG10P,id,ID,R2.group,SE)
      
      # Find duplicated R2 groups
      duplicated_groups <- duplicated(r2.sentinels$R2.group) | duplicated(r2.sentinels$R2.group, fromLast = TRUE)
      
      # Subset the data to keep only rows with duplicated R2 groups
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_max"))
      
      ## Algorithm ####
      # load the data for the phenotypes
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs__ukb_md.tsv.gz"))
      
      # pick up required rows for the selected ids
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        
        #bind the required data
        ll<- rbind(ll,data_filt)
        
      }
      
      # assign the name
      ll$id<- paste0(name,"_ukb")
      
      # select the required columns
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # select the required columns
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      # filter the required column
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_max"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl0<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      mainl0<- mainl0 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      #data1<- data%>%
      #  filter(ID == "rs41280328")
      
      ######################################
      
      filtered_r2.sentinels <- r2.sentinels[!duplicated_groups, ] %>%
        filter(id==paste0(name,"_ukb"))
      
      data<- fread(paste0("output/07_combined_allchrs/",name,"_allchrs_max.tsv.gz"))
      ll<- c()
      for(ids in filtered_r2.sentinels$ID){
        #data_filt<- c()
        #filter the data:
        data_filt <- data %>%
          filter(ID==ids)
        ll<- rbind(ll,data_filt)
        
      }
      
      ll$id<- paste0(name,"_max")
      
      ll<- ll%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      filtered_r2.sentinels<- filtered_r2.sentinels%>%
        dplyr::select(BETA,LOG10P,id,ID,SE)
      
      r2.sentinels1 <- r2.sentinels %>%
        filter(id==paste0(name,"_ukb"))
      
      mainl<- merge(ll,r2.sentinels1, by=c("ID"), all.x=TRUE)
      
      mainl1<- rbind(ll,filtered_r2.sentinels)
      mainl2<- rbind(mainl1,r2.sentinels1, fill=TRUE) 
      
      mainl3<- mainl2 %>%
        distinct(ID,id,LOG10P,BETA,R2.group, .keep_all = TRUE)
      
      mainl3<- mainl3 %>% group_by(ID) %>%
        mutate(R2.group=ifelse(is.na(R2.group),mean(R2.group,na.rm=TRUE),R2.group)) %>%
        distinct(ID,id,LOG10P,BETA,R2.group,.keep_all = TRUE)
      
      main00<- rbind(mainl0, mainl3)
      
      # write the table
      write.table(main00, paste("output/11_BETA_LOG10P_plots/", paste0(name,"_pdata_max.txt"), sep="\t"))
      
      ## plot ####
      # Convert filtered_aa to data.frame
      library(reshape2)
      library(tidyr)
      library(dplyr)
      
      filt_ukb<- main00%>%
        filter(id==paste0(name,"_ukb"))
      
      filt_md <- main00%>%
        filter(id==paste0(name,"_max"))
      
      
      merged_t<- merge(filt_md,filt_ukb,by="R2.group")
      
      merged_t<- merged_t %>%
        drop_na()
      
      merged_t<- merged_t %>%
        dplyr::select(BETA.x,BETA.y,LOG10P.x,LOG10P.y,ID.x,ID.y,id.x,R2.group)%>%
        mutate(BETA.x = abs(BETA.x))%>%
        mutate(BETA.y= abs(BETA.y))
      
      # Plot BETA.x against BETA.y
      #p1<- plot(merged_t$LOG10P.x, merged_t$LOG10P.y, 
      #     xlab = "LOG10P_EHR", ylab = "LOG10P_UKB")
      
      # Optionally, you can add a regression line
      #p1<- abline(lm(BETA.y ~ BETA.x, data = merged_t), col = "red")
      
      p1 <- ggplot(merged_t, aes(x = LOG10P.x, y = LOG10P.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "LOG10P_EHR_max", y = "LOG10P_UKB")+
        facet_wrap(~id.x, scales = "free")
      
      #p1
      # Plot BETA.x against BETA.y
      #p2<- plot(merged_t$BETA.x, merged_t$BETA.y, 
      #     xlab = "BETA_EHR", ylab = "BETA_UKB")
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_LOG10P_max_ukb.png"), sep = ""), plot = p1, width = 10, height = 6, dpi = 500)
      
      p2 <- ggplot(merged_t, aes(x = BETA.x, y = BETA.y)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(x = "BETA_EHR_max", y = "BETA_UKB")+
        facet_wrap(~id.x, scales = "free")
      #p2
      ggsave(paste("output/11_BETA_LOG10P_plots/", paste0(name,"_BETA_max_ukb.png"), sep = ""), plot = p2, width = 10, height = 6, dpi = 500)
      
    }
  },
  error = function(e){
    message('Caught an error!')
    skip_to_next <<- TRUE
    #print(e)
  }
)




