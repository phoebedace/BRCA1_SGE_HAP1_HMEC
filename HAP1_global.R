#### Setup ####

#load packages
library(readxl)
library(mclust)
library(tidyverse)


#setup info
date <- "260304"
setwd("/Users/dacep/Downloads/Run_SGE_scripts/")
options(scipen=999) 

#### Import HAP1 data ####
##### SNVs #####
BRCA1u1_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/220421_BRCA1u1_tHDR_pos_merge_df.csv")
BRCA1x5_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x5_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x10a_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x10a_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x10h2_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/250320_BRCA1x10h2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x10h_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x10h_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x12a_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x12a_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x12b_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x12b_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x11_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x11alt_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x6a_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x6a_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x1_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/220421_BRCA1x1_tHDR_pos_merge_df.csv")
BRCA1i11_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/220421_BRCA1i11_tHDR_pos_merge_df.csv")
BRCA1x10d_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x10d_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")
BRCA1x17q_tHDR_pos_merge_df <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/dfs/240913_BRCA1x17_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv")

#add column to each with name of SGE region
BRCA1u1_tHDR_pos_merge_df$sge_region <- "u1"
BRCA1x5_tHDR_pos_merge_df$sge_region <- "x5"
BRCA1x10a_tHDR_pos_merge_df$sge_region <- "x10a"
BRCA1x10h2_tHDR_pos_merge_df$sge_region <- "x10h2"
BRCA1x10h_tHDR_pos_merge_df$sge_region <- "x10h"
BRCA1x12a_tHDR_pos_merge_df$sge_region <- "x12a"
BRCA1x12b_tHDR_pos_merge_df$sge_region <- "x12b"
BRCA1x11_tHDR_pos_merge_df$sge_region <- "x11"
BRCA1x6a_tHDR_pos_merge_df$sge_region <- "x6a"
BRCA1x1_tHDR_pos_merge_df$sge_region <- "x1"
BRCA1i11_tHDR_pos_merge_df$sge_region <- "i11"
BRCA1x10d_tHDR_pos_merge_df$sge_region <- "x10d"
BRCA1x17q_tHDR_pos_merge_df$sge_region <- "x17q" #q includes p.R1699Q

#bind together all the dfs
BRCA1_all_regions_tHDR_pos_merge_df <- rbind(
  BRCA1u1_tHDR_pos_merge_df,
  BRCA1x5_tHDR_pos_merge_df,
  BRCA1x10a_tHDR_pos_merge_df,
  BRCA1x10h2_tHDR_pos_merge_df,
  BRCA1x10h_tHDR_pos_merge_df,
  BRCA1x12a_tHDR_pos_merge_df,
  BRCA1x12b_tHDR_pos_merge_df,
  BRCA1x11_tHDR_pos_merge_df,
  BRCA1x6a_tHDR_pos_merge_df,
  BRCA1x1_tHDR_pos_merge_df,
  BRCA1i11_tHDR_pos_merge_df,
  BRCA1x10d_tHDR_pos_merge_df,
  BRCA1x17q_tHDR_pos_merge_df
)

##### Indels #####
#Promoter
df_u1_indels <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/250109_BRCA1u1_indels.csv") %>% 
  mutate(sge_region = "u1")

#5' UTR
df_x1_indels <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/250109_BRCA1x1_indels.csv") %>% 
  mutate(sge_region = "x1")

df_u1_x1_indels <- rbind(df_u1_indels, df_x1_indels)


#### Adding additional information ####
#hg38, cHGVS, pHGVS, calculating function scores and RNA scores

#load hg38 conversion and rename columns ready to join dfs
hg38_conversion <- read.csv("230420_BRCA1_hg19_hg38_conversion.csv")
hg38_conversion <- hg38_conversion %>% 
  dplyr::rename("Chrom" = "chr",
         "pos" = "hg19")

#add hg38 coords to dataframe
BRCA1_all_regions_tHDR_pos_merge_df <-
  left_join(x = BRCA1_all_regions_tHDR_pos_merge_df,
            y = hg38_conversion,
            by = c("Chrom", "pos"))

#add pHGVS
BRCA1_all_regions_tHDR_pos_merge_df <- BRCA1_all_regions_tHDR_pos_merge_df %>% 
  mutate(protPos = ifelse(sge_region == "x17" | sge_region == "x17q", protPos - 21, protPos),
         CDSpos = ifelse(sge_region == "x17" | sge_region == "x17q", CDSpos - 63, CDSpos)) %>% 
  mutate(pHGVS = case_when(
  !is.na(Intron) ~ NA_character_,
  conseq == "5PRIME_UTR" ~ NA_character_,
  #include p.R1699Q variant
  !is.na(Exon) | edit_string == "104-X-A,105-X-A" ~ paste0("p.", oAA, protPos, nAA)
)) %>% 
 #add reverse complement for Ref allele
   mutate(Ref_RC = case_when(
    Ref == "A" ~ "T",
    Ref == "C" ~ "G",
    Ref == "T" ~ "A",
    Ref == "G" ~ "C",
    Ref == "CC" ~ "GG"
  )) %>% 
  #add function score (mean and individual for each rep, plus difference between reps) + RNA score
  mutate(rL41_post_pre_function_score = log2(tHDR_post_pre_ratio_synnorm),
         rL42_post_pre_function_score = log2(rL42_tHDR_post_pre_ratio_synnorm),
         diff_post_pre_function_scores = abs(rL41_post_pre_function_score - rL42_post_pre_function_score),
         function_score = (rL41_post_pre_function_score + rL42_post_pre_function_score) / 2,
         RNA_score = case_when(
           !is.na(CDSpos) ~ log2(tHDR_rna_pre_ratio_rL41rL42_mean_synnorm))) %>% 
  mutate(SpliceAI_max = pmax(SpliceAI.acc.gain, SpliceAI.acc.loss, SpliceAI.don.gain, SpliceAI.don.loss))
  

#Adding cHGVS
#distance to splice only goes to +/-20 into the intron - add in for introns that need longer
BRCA1_all_regions_tHDR_pos_merge_df <- BRCA1_all_regions_tHDR_pos_merge_df %>% 
  mutate(Dst2Splice = case_when(is.na(Dst2Splice) & sge_region == "x10d" ~ as.integer(pos-41243452),
                                is.na(Dst2Splice) & sge_region == "x3" & Intron == "3/23" ~ as.integer(pos-41267743),
                                is.na(Dst2Splice) & sge_region == "x3" & Intron == "2/23" ~ as.integer(41267796-pos), 
                                is.na(Dst2Splice) & sge_region == "x12a" & Intron == "11/23" ~ as.integer(41234592-pos),   
                                TRUE ~ Dst2Splice))
#R1699Q information
BRCA1_all_regions_tHDR_pos_merge_df <- BRCA1_all_regions_tHDR_pos_merge_df %>% 
  mutate(pHGVS = case_when(oAA == "R" & protPos == "1699" & nAA == "Q" ~ "p.R1699Q",
                           TRUE ~ pHGVS),
         hg38 = case_when(pos == 41215947 & protPos == "1699" & nAA == "Q" ~ 43063930,
                          TRUE ~ hg38),
  )

#adding cHGVS
BRCA1_all_regions_tHDR_pos_merge_df <- BRCA1_all_regions_tHDR_pos_merge_df %>% 
  mutate(cHGVS = case_when(
    (is.na(CDSpos)&(sge_region == "x1" | sge_region == "u1")) ~ paste0("c.",41277268-as.numeric(as.character(pos)), Ref_RC,">",rev_comp),
    (is.na(CDSpos)&(Intron == "3/23")&(sge_region == "x3")) ~ paste0("c.134+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "2/23")&(sge_region == "x3")) ~ paste0("c.81-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "5/23")&(sge_region == "x5")) ~ paste0("c.301+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "4/23")&(sge_region == "x5")) ~ paste0("c.213-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "5/23")&(sge_region == "x6a")) ~ paste0("c.302-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "9/23")&(sge_region == "x10a")) ~ paste0("c.671-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "10/23")&(sge_region == "x10d")) ~ paste0("c.4096+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "10/23")&(sge_region == "x11")) ~ paste0("c.4097-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "11/23")&(sge_region == "x11")) ~ paste0("c.4185+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "11/23")&(sge_region == "i11")) ~ paste0("c.4185+", 41242961-pos, Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "11/23")&(sge_region == "x12a")) ~ paste0("c.4186-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "12/23")&(sge_region == "x12b")) ~ paste0("c.4357+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "17/23")&(sge_region == "x17q")) ~ paste0("c.5075-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "18/23")&(sge_region == "x17q")) ~ paste0("c.5152+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    !is.na(CDSpos) ~ paste0("c.", CDSpos, Ref_RC, ">", rev_comp)
  ))


#### Filtering ####
#filtering on variant frequency and differences between reps
BRCA1_all_regions_tHDR_pos_merge_df_filtered <- BRCA1_all_regions_tHDR_pos_merge_df %>% 
   filter(lib_pseudo_freq > 0.0001) %>%  
   filter(tHDR_pre_pseudo_freq >= 0.00001 & rL42_tHDR_pre_pseudo_freq >= 0.00001) %>% 
   filter((diff_post_pre_function_scores < 1.5) | (diff_post_pre_function_scores > 1.5 & (rL41_post_pre_function_score < -1 & rL42_post_pre_function_score < -1)))

#filtering out specific variants
 BRCA1_all_regions_tHDR_pos_merge_df_filtered <- BRCA1_all_regions_tHDR_pos_merge_df_filtered %>% 
   filter(pos != 41243496) %>%  #filter out x10d variants between 2 PAM edits: (labelled as missense but without PAM edits would be nonsense)
   filter(!(pos == 41256962 & alt == "A")) %>%  #x5 variant that in presence of PAM site introduces splice site
   filter(!(pos == 41246779 & alt == "A")) #x10a variant that in absence of PAM would introduce splice site with high SpliceAI score.
 
 #filtering out variants without cHGVS (variants affected by PAM edits)
  BRCA1_all_regions_tHDR_pos_merge_df_filtered <- BRCA1_all_regions_tHDR_pos_merge_df_filtered %>% 
   filter(!is.na(cHGVS))
 

  #### Are variants significantly depleted? ####
 # global scaling of filtered data by SYN and NS variants
 # general formula:  b/(a-c)*(A-C) where b is snv's function score, a is sge_region med syn, c is sge_region med nonsense, A is global med syn, C is global med nonsense
BRCA1_function_score_median_syn <- median(BRCA1_all_regions_tHDR_pos_merge_df_filtered[which(BRCA1_all_regions_tHDR_pos_merge_df_filtered$conseq == 'SYNONYMOUS'),c('function_score')])
BRCA1_function_score_median_ns <- median(BRCA1_all_regions_tHDR_pos_merge_df_filtered[which(BRCA1_all_regions_tHDR_pos_merge_df_filtered$conseq == 'STOP_GAINED') ,c('function_score')])
 
##### Coding regions #####
regions <- c("x12a", "x12b", "x11", "x10a", "x10h", "x10h2", "x10d", "x6a", "x5", "x17q") #coding only
  
BRCA1region_final_df = ""
 all_regions_df = data.frame()
 
 add_fdr <- function(region) {
   BRCA1region_final_df <- BRCA1_all_regions_tHDR_pos_merge_df_filtered[which(BRCA1_all_regions_tHDR_pos_merge_df_filtered$sge_region == region),] 
   BRCA1region_function_score_median_syn <- median(BRCA1region_final_df[which(BRCA1region_final_df$conseq == 'SYNONYMOUS'),c('function_score')])
   BRCA1region_function_score_median_ns <- median(BRCA1region_final_df[which(BRCA1region_final_df$conseq == 'STOP_GAINED'),c('function_score')])
   BRCA1region_final_df$function_score_sns <- (BRCA1region_final_df$function_score/(BRCA1region_function_score_median_syn-BRCA1region_function_score_median_ns))*(BRCA1_function_score_median_syn-BRCA1_function_score_median_ns)
   #have added in a function score limit for null because x17 has no RNA scores and some low synonymous
   #function score >-0.8 is equivalent to filtering x17 scores with low RNA in Findlay 2018.
   BRCA1region_null_fss <- BRCA1region_final_df[which(BRCA1region_final_df$conseq == "SYNONYMOUS" & BRCA1region_final_df$RNA_score > -1 & BRCA1region_final_df$function_score > -0.8),]$function_score_sns
   BRCA1region_null_mean <- mean(BRCA1region_null_fss)
   BRCA1region_null_sd <- sd(BRCA1region_null_fss)
   BRCA1region_final_df$pvalues <- pnorm(BRCA1region_final_df$function_score_sns, BRCA1region_null_mean, sd=BRCA1region_null_sd) 
   BRCA1region_final_df$fdr <- p.adjust(BRCA1region_final_df$pvalues, method = "BH")
   BRCA1region_final_df$fs_sig <- "FDR NA"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr > 0.5),]$fs_sig <- "FDR > 0.5"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.5),]$fs_sig <- "FDR < 0.5"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.2),]$fs_sig <- "FDR < 0.2"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.05),]$fs_sig <- "FDR < 0.05"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.01),]$fs_sig <- "FDR < 0.01"
   
   #SAME BUT FOR NONSENSE VARIANTS - to create LoF/intermediate threshold
   #NB this is not the final method for determining LoF vs intermediate - see LoF vs intermediate section below
   BRCA1_nonsense <- BRCA1region_final_df[which(BRCA1region_final_df$conseq == "STOP_GAINED" & BRCA1region_final_df$function_score_sns < -0.9),]$function_score_sns 
   BRCA1_nonsense_mean <- mean(BRCA1_nonsense)   
   BRCA1_nonsense_sd <- sd(BRCA1_nonsense)
   BRCA1region_final_df$pvalues_n <- pnorm(BRCA1region_final_df$function_score_sns, BRCA1_nonsense_mean, BRCA1_nonsense_sd, lower.tail = FALSE)
   BRCA1region_final_df$fdr_n <- p.adjust(BRCA1region_final_df$pvalues_n, method = "BH")
   BRCA1region_final_df <- BRCA1region_final_df %>%
     mutate(
       fs_sig_n = case_when(
         fdr_n > 0.5 ~ "FDR > 0.5",
         fdr_n < 0.5 & fdr_n >= 0.2 ~ "FDR < 0.5",
         fdr_n < 0.2 & fdr_n >= 0.05 ~ "FDR < 0.2",
         fdr_n < 0.05 & fdr_n >= 0.01 ~ "FDR < 0.05",
         fdr_n < 0.01 ~ "FDR < 0.01",
         TRUE ~ "FDR NA"
       ))
   if (any(BRCA1region_final_df$fs_sig == "FDR < 0.01")) {
     # Calculate max0.01 and min0.05
     upper_max0.01 <- max(BRCA1region_final_df[BRCA1region_final_df$fs_sig == "FDR < 0.01", ]$function_score_sns)
     upper_min0.05 <- min(BRCA1region_final_df[BRCA1region_final_df$fs_sig != "FDR < 0.01", ]$function_score_sns)
     
     # Calculate fs_threshold (based on syn vars)
     BRCA1region_final_df$fs_threshold_upper <- (upper_max0.01 + upper_min0.05) / 2
   } else {
     # If no values where fs_sig == "FDR < 0.01", set max0.01, min0.05, and fs_threshold to the minimum function score
     upper_max0.01 <- NA
     upper_min0.05 <- NA
     BRCA1region_final_df$fs_threshold_upper <- NA
   }
   if (any(BRCA1region_final_df$fs_sig_n == "FDR < 0.01")) {
     # Calculate max0.01 and min0.05
     lower_min0.01 <- min(BRCA1region_final_df[BRCA1region_final_df$fs_sig_n == "FDR < 0.01", ]$function_score_sns)
     lower_max0.05 <- max(BRCA1region_final_df[BRCA1region_final_df$fs_sig_n != "FDR < 0.01", ]$function_score_sns)
     
     # Calculate fs_threshold (based on nonsense vars)
     BRCA1region_final_df$fs_threshold_lower <- (lower_min0.01 + lower_max0.05) / 2
   } else {
     # If no values where fs_sig == "FDR < 0.01", set max0.01, min0.05, and fs_threshold to the minimum function score
     lower_min0.01 <- NA
     lower_max0.05 <- NA
     BRCA1region_final_df$fs_threshold_lower <- NA
   }
   
   #write.csv(x = BRCA1region_final_df, file = paste0("dfs/", date, "_", region, "_filtered_df.csv"), row.names = FALSE)
   return(BRCA1region_final_df)
 }
 
 for (region in regions) {
   output <- add_fdr(region)
   all_regions_df <- rbind(all_regions_df, output) #NB this is coding regions only - can't normalise using syn/nonsense vars in non coding expts
 }
 
 
 ##### Non-coding regions #####
 regions <- c("i11", "u1", "x1") #non-coding only
 
 BRCA1region_final_df = ""
 nc_regions_df = data.frame()
 
 add_fdr_nc <- function(region) {
   BRCA1region_final_df <- BRCA1_all_regions_tHDR_pos_merge_df_filtered[which(BRCA1_all_regions_tHDR_pos_merge_df_filtered$sge_region == region),] 
   BRCA1region_final_df$function_score_sns <- BRCA1region_final_df$function_score
   BRCA1region_null_fss <- BRCA1region_final_df$function_score_sns
   BRCA1region_null_mean <- mean(BRCA1region_null_fss)
   BRCA1region_null_sd <- sd(BRCA1region_null_fss)
   BRCA1region_final_df$pvalues <- pnorm(BRCA1region_final_df$function_score_sns, BRCA1region_null_mean, sd=BRCA1region_null_sd) 
   BRCA1region_final_df$fdr <- p.adjust(BRCA1region_final_df$pvalues, method = "BH")
   BRCA1region_final_df$fs_sig <- "FDR NA"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr > 0.5),]$fs_sig <- "FDR > 0.5"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.5),]$fs_sig <- "FDR < 0.5"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.2),]$fs_sig <- "FDR < 0.2"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.05),]$fs_sig <- "FDR < 0.05"
   BRCA1region_final_df[which(BRCA1region_final_df$fdr < 0.01),]$fs_sig <- "FDR < 0.01"
   #write.csv(x = BRCA1region_final_df, file = paste0("dfs/", date, "_", region, "_filtered_df.csv"), row.names = FALSE)
   if (any(BRCA1region_final_df$fs_sig == "FDR < 0.01")) {
     # Calculate max0.01 and min0.05
     max0.01 <- max(BRCA1region_final_df[BRCA1region_final_df$fs_sig == "FDR < 0.01", ]$function_score_sns)
     min0.05 <- min(BRCA1region_final_df[BRCA1region_final_df$fs_sig != "FDR < 0.01", ]$function_score_sns)
     
     # Calculate fs_threshold
     BRCA1region_final_df$fs_threshold_upper <- (max0.01 + min0.05) / 2
   } else {
     # If no values where fs_sig == "FDR < 0.01", set max0.01, min0.05, and fs_threshold to the minimum function score
     max0.01 <- NA
     min0.05 <- NA
     BRCA1region_final_df$fs_threshold_upper <- NA
   }
   BRCA1region_final_df$pvalues_n <- NA
   BRCA1region_final_df$fdr_n <- NA
   BRCA1region_final_df$fs_sig_n <- NA
   BRCA1region_final_df$fs_threshold_lower <- BRCA1region_final_df$fs_threshold_upper
   
   return(BRCA1region_final_df)
 }
 
 
 for (region in regions) {
   output <- add_fdr_nc(region)
   nc_regions_df <- rbind(nc_regions_df, output) #NB this is NON-coding regions only - can't normalise using syn/nonsense vars in non coding expts
 }
 
 
 #merging coding regions with non-coding df
 coding_nc_regions_df <- rbind(all_regions_df, nc_regions_df)
 
 #adding classifications based on fs_thresholds
 coding_nc_regions_df <- coding_nc_regions_df %>% 
   mutate(HAP1_func_class = case_when(function_score_sns > fs_threshold_upper ~ "Neutral",
                                      function_score_sns < fs_threshold_lower ~ "LoF",
                                      (function_score_sns < fs_threshold_upper) & (function_score_sns > fs_threshold_lower) ~ "Intermediate")) 
 

 ##### LoF vs intermediate #####
 #Alternative method to determine LoF vs intermediate - this was the method that was used for final function classes.
 #uses nonsense that are significantly depleted for null distribution
 
 coding_nc_regions_df_thresh2 <- coding_nc_regions_df
 
 BRCA1_nonsense <- coding_nc_regions_df_thresh2[which(coding_nc_regions_df_thresh2$conseq == "STOP_GAINED" & coding_nc_regions_df_thresh2$fdr < 0.01),]$function_score_sns 
 BRCA1_nonsense_mean <- mean(BRCA1_nonsense)   
 BRCA1_nonsense_sd <- sd(BRCA1_nonsense)
 coding_nc_regions_df_thresh2$pvalues_n2 <- pnorm(coding_nc_regions_df_thresh2$function_score_sns, BRCA1_nonsense_mean, BRCA1_nonsense_sd, lower.tail = FALSE)
 coding_nc_regions_df_thresh2$fdr_n2 <- p.adjust(coding_nc_regions_df_thresh2$pvalues_n2, method = "BH")
 coding_nc_regions_df_thresh2 <- coding_nc_regions_df_thresh2 %>%
   mutate(
     fs_sig_n2 = case_when(
       fdr_n2 > 0.5 ~ "FDR > 0.5",
       fdr_n2 < 0.5 & fdr_n2 >= 0.2 ~ "FDR < 0.5",
       fdr_n2 < 0.2 & fdr_n2 >= 0.05 ~ "FDR < 0.2",
       fdr_n2 < 0.05 & fdr_n2 >= 0.01 ~ "FDR < 0.05",
       fdr_n2 < 0.01 ~ "FDR < 0.01",
       TRUE ~ "FDR NA"
     ))
 
 if (any(coding_nc_regions_df_thresh2$fs_sig_n2 == "FDR < 0.01")) {
   # Calculate max0.01 and min0.05
   lower_min0.01 <- min(coding_nc_regions_df_thresh2[coding_nc_regions_df_thresh2$fs_sig_n2 == "FDR < 0.01", ]$function_score_sns)
   lower_max0.05 <- max(coding_nc_regions_df_thresh2[coding_nc_regions_df_thresh2$fs_sig_n2 != "FDR < 0.01", ]$function_score_sns)
   
   # Calculate fs_threshold (based on nonsense vars)
   coding_nc_regions_df_thresh2$fs_threshold_lower2 <- (lower_min0.01 + lower_max0.05) / 2
 } else {
   # If no values where fs_sig == "FDR < 0.01", set max0.01, min0.05, and fs_threshold to the minimum function score
   lower_min0.01 <- NA
   lower_max0.05 <- NA
   coding_nc_regions_df_thresh2$fs_threshold_lower2 <- NA
 }
 
 coding_nc_regions_df_thresh2 <- coding_nc_regions_df_thresh2 %>% 
   mutate(HAP1_func_class_v2 = case_when(function_score_sns > fs_threshold_upper ~ "Neutral",
                                      function_score_sns < fs_threshold_lower2 ~ "LoF",
                                      (function_score_sns < fs_threshold_upper) & (function_score_sns > fs_threshold_lower2) ~ "Intermediate")) 
 
 
 coding_nc_regions_df <- coding_nc_regions_df_thresh2
 #HAP1_func_class_v2 should be used instead of HAP1_func_class
 
 #### Add additional data #### 
 ##### ClinVar #####
 #Data downloaded from https://www.ncbi.nlm.nih.gov/clinvar
 clinvar_download <- read.table(file = "clinvar/240913_clinvar.txt", sep="\t", header = TRUE, fill = TRUE, na.strings = "")
 
 clinvar_download <- clinvar_download %>% 
   dplyr::rename("Review.status" = "Germline.review.status",
          "Clinvar_interpretation_Sep24" = "Germline.classification",
          "Last_reviewed" = "Germline.date.last.evaluated")
 
 clinvar_download[, c("Ref", "alt")] <- stringr::str_match(
   string = clinvar_download$Canonical.SPDI,
   pattern = ".+([ACGT]):([ATCG])"
 )[,2:3]
 
 clinvar_download2 <- clinvar_download %>% 
   select(GRCh37Location, Ref, alt, GRCh38Location, Clinvar_interpretation_Sep24, Review.status, Last_reviewed, Canonical.SPDI, VariationID) %>% 
   mutate(Clinvar_Sep24_histo = case_when(grepl(pattern = "ncertain", x = Clinvar_interpretation_Sep24) ~ "VUS/conflicting interpretations",
                                          grepl(pattern = "interpretations", x = Clinvar_interpretation_Sep24) ~ "VUS/conflicting interpretations",
                                          grepl(pattern = "athogenic", x = Clinvar_interpretation_Sep24) ~ "Pathogenic/Likely pathogenic",
                                          grepl(pattern = "enign", x = Clinvar_interpretation_Sep24) ~ "Benign/Likely benign",
                                          #  grepl(pattern = "ikely", x = Clinvar_interpretation_Sep24) & grepl(pattern = "athogenic", x = Clinvar_interpretation_Sep24) ~ "Likely pathogenic",
                                          TRUE ~ Clinvar_interpretation_Sep24),
          Clinvar_link = paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", VariationID))
 
coding_nc_regions_df_newclinvar <- left_join(y = clinvar_download2,
                                              x = coding_nc_regions_df,
                                              by = c("pos" = "GRCh37Location",
                                                     "Ref" = "Ref",
                                                     "alt" = "alt"))
 
 coding_nc_regions_df_newclinvar$Clinvar_interpretation_Sep24 <- ifelse(is.na(coding_nc_regions_df_newclinvar$Clinvar_interpretation_Sep24), "Absent", coding_nc_regions_df_newclinvar$Clinvar_interpretation_Sep24)
 coding_nc_regions_df_newclinvar$Clinvar_Sep24_histo <- ifelse(is.na(coding_nc_regions_df_newclinvar$Clinvar_Sep24_histo), "Absent", coding_nc_regions_df_newclinvar$Clinvar_Sep24_histo)
 coding_nc_regions_df_newclinvar$Clinvar_Sep24_histo <- factor(coding_nc_regions_df_newclinvar$Clinvar_Sep24_histo, levels = c("Absent", "VUS/conflicting interpretations", "Benign/Likely benign", "Pathogenic/Likely pathogenic"))
 

 ##### SpliceAI #####
 #SpliceAI run with 500bp distance
 df_spliceAI_scores_new <- read.csv(file = "~/OneDrive - The Francis Crick Institute/R/250127_new_spliceAI_scores/250703_SpliceAI_scores_500bp_all_vars.csv")   
 
 df_spliceAI_scores_new2 <- df_spliceAI_scores_new %>% 
   select(cHGVS, hg38, Transcript, AG, AL, DG, DL) %>% 
   dplyr::rename("SpliceAI_transcript" = "Transcript")
 
 #joining new spliceAI with SGE data
 coding_nc_regions_df2 <- left_join(x = coding_nc_regions_df_newclinvar,
                                   y = df_spliceAI_scores_new2,
                                   by = c("cHGVS", "hg38")) %>% 
   mutate(SAI_lookup = paste0("chr17-", hg38, "-", Ref, "-", alt))
 

 #### Removing duplicated scores ####
 #keep x12b where there is overlap with x12a
   coding_nc_regions_df_nodups <- coding_nc_regions_df2 %>% 
   mutate(x12_duplicate_score = case_when(sge_region == "x12a" & pos > 41234500 & pos < 41234526 ~ "Yes",
                                          TRUE ~ "No")) %>% 
   filter(x12_duplicate_score == "No") %>% 
   select(-x12_duplicate_score)
   
 #also remove duplicate scores in x1/u1 - keep u1 scores as assay much cleaner
 coding_nc_regions_df_nodups2 <- coding_nc_regions_df_nodups %>% 
   mutate(x1_duplicate_score = case_when(sge_region == "x1" & pos >= 41277354 & pos <= 41277397 ~ "Yes",
                                         TRUE ~ "No")) %>% 
   filter(x1_duplicate_score == "No") %>% 
   select(-x1_duplicate_score)
 

#### Defining truth set and GMM ####
df_truthset <- coding_nc_regions_df_nodups2 %>% 
  filter(Clinvar_Sep24_histo != "VUS/conflicting interpretations" & Clinvar_Sep24_histo != "Absent") %>% 
  filter(Review.status != "criteria provided, single submitter" & Review.status != "no assertion criteria provided") %>% 
  #Remove benign variants with a high spliceAI score
  filter(!(Clinvar_Sep24_histo == "Benign/Likely benign" & SpliceAI_max > 0.2))

#Mixture modelling
  df_truthset_gmm <- df_truthset %>% 
    mutate(numeric_label = case_when(Clinvar_Sep24_histo == "Benign/Likely benign" ~ 1,
                                     Clinvar_Sep24_histo == "Pathogenic/Likely pathogenic" ~ 0
                                     )) %>% 
    select(cHGVS, numeric_label, function_score_sns)
  
  gmm_model <- Mclust(df_truthset_gmm$function_score_sns, G = 2)
  summary(gmm_model)
  
  
  df_vars_for_pvalues <- coding_nc_regions_df_nodups2
  
  pathogenic_index <- which.min(gmm_model$parameters$mean)  # Lower mean = pathogenic
  df_truthset_gmm$prob_pathogenic <- gmm_model$z[, pathogenic_index]
  df_vars_for_pvalues$prob_pathogenic <- predict(gmm_model, newdata = df_vars_for_pvalues$function_score_sns)$z[, pathogenic_index]
  
   df_vars_for_pvalues$prob_pathogenic <- predict(gmm_model, newdata = df_vars_for_pvalues$function_score_sns)$z[, which.min(gmm_model$parameters$mean)]
  
  prior_prob <- 0.1
  Threshold <- c(-Inf, 0.053, 0.23, 0.48, 2.1, 4.3, 18.7, Inf)
  Labels <- c("BS3", "BS3_mod", "BS3_supp", "indeterminate", "PS3_supp", "PS3_mod", "PS3")
  
  df_vars_for_pvalues2 <- df_vars_for_pvalues %>% 
    select(cHGVS, pHGVS, pos, hg38, Ref, alt, sge_region, conseq, Clinvar_interpretation_Sep24, Clinvar_Sep24_histo, Clinvar_link, Review.status, function_score_sns, prob_pathogenic, HAP1_func_class_v2) %>% 
    mutate(OddsPath = (prob_pathogenic*(1-prior_prob))/((1-prob_pathogenic)*prior_prob)) %>% 
    mutate(evidence_code = cut(OddsPath, breaks = Threshold, labels = Labels))
 
 #importing evidence points conversion 
  evidence_code_points <- read.csv("250207_evidence_code_points.csv")
  
  #variants with high function scores outside of distribution assigned BS3
  df_vars_for_pvalues2_corrected <- df_vars_for_pvalues2 %>% 
    mutate(OddsPath_corrected = case_when(function_score_sns > 0.35 ~ 0.05,
                                     TRUE ~ OddsPath)) %>% 
    mutate(evidence_code = case_when(function_score_sns > 0.35 ~ "BS3",
                                     TRUE ~ evidence_code))
  
  
  df_vars_for_pvalues2_corrected_with_points <- left_join(x = df_vars_for_pvalues2_corrected,
                                                  y = evidence_code_points,
                                                  by = "evidence_code")
  
  df_vars_for_pvalues2_corrected_with_points$evidence_code <- factor(df_vars_for_pvalues2_corrected_with_points$evidence_code, levels = c("BS3", "BS3_mod", "BS3_supp", "indeterminate", "PS3_supp", "PS3_mod", "PS3"))
  
  ##### Modifying evidence codes for neutral/intermediate pathogenic variants #####
  
  df_vars_for_pvalues2_corrected_with_points_modified <- df_vars_for_pvalues2_corrected_with_points %>% 
    dplyr::rename("HAP1_evidence_code_original" = "evidence_code",
                  "HAP1_points_original" = "points") %>% 
    #Adjust HAP1 evidence codes: neutral path become indeterminate. Int path become P_supporting.
    mutate(HAP1_ev_code_adj = case_when(HAP1_func_class_v2 == "Neutral" & (HAP1_evidence_code_original == "PS3" | HAP1_evidence_code_original == "PS3_supp" | HAP1_evidence_code_original == "PS3_mod") ~ "indeterminate",
                                        HAP1_func_class_v2 == "Intermediate" & (HAP1_evidence_code_original == "PS3" | HAP1_evidence_code_original == "PS3_mod") ~ "PS3_supp",
                                        TRUE ~ HAP1_evidence_code_original),
           HAP1_points_adj = case_when(HAP1_func_class_v2 == "Neutral" & (HAP1_evidence_code_original == "PS3" | HAP1_evidence_code_original == "PS3_supp" | HAP1_evidence_code_original == "PS3_mod") ~ 0,
                                       HAP1_func_class_v2 == "Intermediate" & (HAP1_evidence_code_original == "PS3" | HAP1_evidence_code_original == "PS3_mod") ~ 1,
                                       TRUE ~ HAP1_points_original)
           )


  HAP1_points_mclust_with_rules2 <- df_vars_for_pvalues2_corrected_with_points_modified %>% 
     select(cHGVS, prob_pathogenic, HAP1_evidence_code_original, HAP1_points_original, OddsPath_corrected, function_score_sns, HAP1_ev_code_adj, HAP1_points_adj) %>% 
    dplyr::rename("HAP1_prob_pathogenic" = "prob_pathogenic")



#### DATA EXPORT FOR SUPPLEMENTARY TABLES 1 AND 2 ####
##### HAP1 data - Supplementary Table 1 #####
#final version of adjusted evidence codes (22/7/25)
HAP1_ev_codes_for_supp <- df_vars_for_pvalues2_corrected_with_points_modified %>% 
  select(cHGVS, sge_region, prob_pathogenic, OddsPath, OddsPath_corrected, HAP1_evidence_code_original, HAP1_points_original, HAP1_ev_code_adj, HAP1_points_adj)
  
truthset_list_df <- df_truthset %>% 
  mutate(cHGVS_region = paste0(cHGVS, sge_region))

truthset_list <- as.character(truthset_list_df$cHGVS_region)

coding_nc_regions_df_export <- left_join(x = coding_nc_regions_df2,
          y = HAP1_ev_codes_for_supp, #df_vars_for_pvalues2_corrected_supp,
          by = c("cHGVS", "sge_region"))


#label which variants are in the truthset
coding_nc_regions_df_export2 <- coding_nc_regions_df_export %>% 
  mutate(cHGVS_region = paste0(cHGVS, sge_region)) %>% 
  mutate(in_truthset = case_when(cHGVS_region %in% truthset_list ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  select(-cHGVS_region)

#add JASPAR scores
jaspar_data <- read.csv("jaspar_data.csv")

coding_nc_regions_df_export3 <- left_join(x = coding_nc_regions_df_export2,
                                          y = jaspar_data,
                                          by = "cHGVS")

#Label variants used as final score (i.e. where a variant is assayed in 2 experiments, which one is used?)
df_final_score <- coding_nc_regions_df_nodups2 %>%
  mutate(cHGVS_region = paste0(cHGVS, sge_region))
final_scores <- as.character(df_final_score$cHGVS_region)

coding_nc_regions_df_export4 <- coding_nc_regions_df_export3 %>% 
  mutate(cHGVS_region = paste0(cHGVS, sge_region)) %>% 
  mutate(final_score = case_when(cHGVS_region %in% final_scores ~ TRUE,
                                 TRUE ~ FALSE))

coding_nc_regions_df_export5 <- coding_nc_regions_df_export4 %>% 
  select(hg38, sge_region, conseq, Ref_RC, rev_comp, Dst2Splice, cHGVS, pHGVS, Ref, alt, protPos, oAA, nAA, CDSpos, Intron, Exon, function_score_sns, rL41_post_pre_function_score, rL42_post_pre_function_score, RNA_score, fdr, fs_threshold_upper, fs_threshold_lower2, HAP1_func_class_v2, Clinvar_interpretation_Sep24, Clinvar_Sep24_histo, Review.status, CADD.phred, variant, edit_string, tHDR_pre, tHDR_post, tHDR_rna, rL42_tHDR_pre, rL42_tHDR_post, rL42_tHDR_rna, tHDR_lib, tHDR_neg, prob_pathogenic, OddsPath, OddsPath_corrected, HAP1_evidence_code_original, HAP1_points_original, HAP1_ev_code_adj, HAP1_points_adj, in_truthset, Relative.score, AL, AG, DL, DG, SpliceAI_transcript, final_score) %>% 
  mutate(cHGVS = case_when(hg38 > 43125364 ~ NA_character_,
                           pHGVS == "p.R1699Q" ~ "c.5096_5097delinsAA",
                           TRUE ~ cHGVS),
         Exon = case_when(hg38 > 43125364 ~ NA_character_,
                           TRUE ~ Exon),
         conseq = case_when(hg38 > 43125364 ~ "Promoter",
                          TRUE ~ conseq),
         #No RNA scores for x17 - change RNA score to NA and raw counts to NA
         RNA_score = case_when(sge_region == "x17q" ~ NA_real_,
                               TRUE ~ RNA_score),
         #RNA scores for intronic regions are 1 - make NA
         tHDR_rna = case_when(sge_region == "x17q" ~ NA_integer_,
                              is.na(RNA_score) ~ NA_integer_,
                              TRUE ~ tHDR_rna),
         rL42_tHDR_rna = case_when(sge_region == "x17q" ~ NA_integer_,
                                   is.na(RNA_score) ~ NA_integer_,
                                   TRUE ~ rL42_tHDR_rna),
         conseq = case_when(conseq == "NON_SYNONYMOUS" ~ "Missense",
                            conseq == "STOP_GAINED" ~ "Nonsense",
                            conseq == "SYNONYMOUS" ~ "Synonymous",
                            conseq == "INTRONIC" ~ "Intronic",
                            conseq == "SPLICE_SITE" ~ "Splice region",
                            conseq == "CANONICAL_SPLICE" ~ "Canonical splice",
                            conseq == "5PRIME_UTR" ~ "5' UTR",
                            TRUE ~ conseq),
         pHGVS = case_when(sge_region == "u1" | sge_region == "x1" ~ NA_character_,
                           TRUE ~ pHGVS),
         #need to correct position for R1699Q - Should be 17_43063929_CC_TT hg38
         #cHGVS also needs correcting: c.5096_5097delinsAA
         hg38 = case_when(pHGVS == "p.R1699Q" ~ 43063929,
                          TRUE ~ hg38),
         fs_threshold_lower2 = case_when(fs_threshold_upper < fs_threshold_lower2 ~ NA_real_,
                                         TRUE ~ fs_threshold_lower2),
         sge_region = case_when(sge_region == "u1" ~ "promoter",
                                sge_region == "x1" ~ "5' UTR",
                                sge_region == "x5" ~ "exon 5",
                                sge_region == "x6a" ~ "exon 6",
                                sge_region == "x17q" ~ "exon 17",
                                sge_region == "x12a" ~ "exon 12 (5')",
                                sge_region == "x12b" ~ "exon 12 (3')",
                                sge_region == "x10a" ~ "exon 10 (5')",
                                sge_region == "x10d" ~ "exon 10 (3')",
                                sge_region == "x10h2" ~ "exon 10 (mid 1)",
                                sge_region == "x10h" ~ "exon 10 (mid 2)",
                                sge_region == "x11" ~ "exon 11",
                                sge_region == "i11" ~ "intron 11"
                                ),
         evidence_code_original = case_when(HAP1_evidence_code_original == "PS3" ~ "P_strong",
                                   HAP1_evidence_code_original == "PS3_mod" ~ "P_moderate",
                                   HAP1_evidence_code_original == "PS3_supp" ~ "P_supporting",
                                   HAP1_evidence_code_original == "BS3" ~ "B_strong",
                                   HAP1_evidence_code_original == "BS3_mod" ~ "B_moderate",
                                   HAP1_evidence_code_original == "BS3_supp" ~ "B_supporting",
                                   TRUE ~ HAP1_evidence_code_original),
         evidence_code_adjusted = case_when(HAP1_ev_code_adj == "PS3" ~ "P_strong",
                                            HAP1_ev_code_adj == "PS3_mod" ~ "P_moderate",
                                            HAP1_ev_code_adj == "PS3_supp" ~ "P_supporting",
                                            HAP1_ev_code_adj == "BS3" ~ "B_strong",
                                            HAP1_ev_code_adj == "BS3_mod" ~ "B_moderate",
                                            HAP1_ev_code_adj == "BS3_supp" ~ "B_supporting",
                                            TRUE ~ HAP1_ev_code_adj),
         ) %>%
  mutate(SpliceAI_max = pmax(AG, AL, DG, DL)) %>% 
  #Get rid of spliceAI scores for variants not in transcript
  mutate(SpliceAI_max = case_when(SpliceAI_transcript != "ENST00000357654.9" ~ NA_real_,
                                  TRUE ~ SpliceAI_max)) %>% 
  dplyr::rename("ClinVar" = "Clinvar_interpretation_Sep24",
                "ClinVar_simple" = "Clinvar_Sep24_histo",
                "Consequence" = "conseq",
                "SGE_region" = "sge_region",
                "Alt" = "alt",
                "final_function_score" = "function_score_sns",
                "r1_function_score" = "rL41_post_pre_function_score",
                "r2_function_score" = "rL42_post_pre_function_score",
                "RNA_score" = "RNA_score",
                "q_value" = "fdr",
                "function_class" = "HAP1_func_class_v2",
                "r1_D5" = "tHDR_pre",
                "r2_D5" = "rL42_tHDR_pre",
                "r1_D14" = "tHDR_post",
                "r2_D14" = "rL42_tHDR_post",
                "r1_D5_RNA" = "tHDR_rna",
                "r2_D5_RNA" = "rL42_tHDR_rna",
                "library" = "tHDR_lib",
                "negative" = "tHDR_neg",
                "threshold_upper" = "fs_threshold_upper",
                "threshold_lower" = "fs_threshold_lower2",
                "E2F_binding_score" = "Relative.score",
                "ClinVar_review_status" = "Review.status",
                "points_original" = "HAP1_points_original",
                "points_adjusted" = "HAP1_points_adj"
                )


#Removing intron 1 variants from main df
coding_nc_regions_df_export_correct_cHGVS <- coding_nc_regions_df_export5 %>% 
  filter(Intron != "1/23" | is.na(Intron))

#Correcting cHGVS for intron 1 variants
coding_nc_regions_df_export_cHGVS_corrected <- coding_nc_regions_df_export5 %>% 
  filter(Intron == "1/23") %>% 
  mutate(cHGVS = paste0("c.-20+", -Dst2Splice, Ref_RC, ">", rev_comp))

#binding back together
coding_nc_regions_df_export6 <- 
rbind(coding_nc_regions_df_export_correct_cHGVS,
      coding_nc_regions_df_export_cHGVS_corrected
      )

#select columns
coding_nc_regions_df_export7 <- coding_nc_regions_df_export6 %>%
  select(hg38, Alt, Ref, cHGVS, pHGVS, protPos, oAA, nAA, Consequence, CDSpos, SGE_region, final_function_score, q_value, r1_function_score, r2_function_score, RNA_score, threshold_upper, threshold_lower, function_class, final_score, variant, edit_string, ClinVar, ClinVar_simple, ClinVar_review_status, AG, AL, DG, DL, SpliceAI_max, E2F_binding_score, library, negative, r1_D5, r1_D14, r2_D5, r2_D14, r1_D5_RNA, r2_D5_RNA, in_truthset, prob_pathogenic, OddsPath, OddsPath_corrected, evidence_code_original, points_original, evidence_code_adjusted, points_adjusted) %>% 
  arrange(desc(hg38))

#write.csv(x = coding_nc_regions_df_export7, file = paste0("~/Dropbox (The Francis Crick)/BRCA1_HAP1_HMEC_manuscript/SUPPLEMENTARY_TABLES/", date, "_ST1.csv"), row.names = FALSE) 

dim(coding_nc_regions_df_export7)
coding_nc_regions_df_export7 %>% 
  group_by(final_score) %>% 
  summarise(n())
#4,657 scores for 4,470 unique variants

##### HAP1 indel scores - Supplementary Table 2 #####

df_x1_indels_v2 <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/250711_BRCA1x1_indels.csv") %>% 
  mutate(sge_region = "x1")
df_u1_indels_v2 <- read.csv(file = "/Users/dacep/Downloads/Run_SGE_scripts/250711_BRCA1u1_indels.csv") %>% 
  mutate(sge_region = "u1")

df_u1_x1_indels_v2 <- rbind(df_x1_indels_v2,
                            df_u1_indels_v2)

  df_indels <- df_u1_x1_indels_v2 %>% 
    #filter variants with difference >1.5 between replicates - removes 4 variants
    mutate(diff = abs(r1_function_score - r2_function_score)) %>% 
    filter(diff < 1.5) %>% 
    #filter out any reads with non-HDR mutations
    filter(!grepl(pattern = "73-X-C", x = edit_string)) %>%
    filter(!grepl(pattern = "158-X-A", x = edit_string)) %>% 
    filter(!grepl(pattern = "155-X-A", x = edit_string)) %>% 
    filter(!grepl(pattern = "145-X-G", x = edit_string)) %>% 
    filter(!grepl(pattern = "95-X-T", x = edit_string)) %>% 
    filter(!grepl(pattern = "155-X-C", x = edit_string)) %>%
    filter(!grepl(pattern = "68-X-G", x = edit_string)) %>%
    filter(!grepl(pattern = "64-X-G", x = edit_string)) %>%
    filter(!grepl(pattern = "155-X-C", x = edit_string)) %>%
    #filter vars with one PAM edit, when score for both is available: 
    filter(edit_string != "88-D5,144-X-G") %>%
    filter(edit_string != "93-D1,144-X-G") %>%
    filter(edit_string != "94-X-T,143-D1") %>%
    filter(edit_string != "63-X-G,152-D1") %>%
  select(variant, cHGVS, cigar, edit_string, indel_size, pos_first_deleted, post_last_deleted, sge_region, r1_function_score, r2_function_score, function_score, pre, post, rL42_pre, rL42_post, lib, neg) %>% 
  dplyr::rename("pos_last_deleted" = "post_last_deleted") %>% 
  mutate(sge_region = case_when(sge_region == "u1" ~ "promoter",
                         sge_region == "x1" ~ "5' UTR"
                         ))

  hg38_conversion_indels <- hg38_conversion %>% 
  select(pos, hg38)

df_indels2 <- left_join(x = df_indels, 
          y = hg38_conversion_indels,
          by = c("pos_first_deleted" = "pos")) %>% 
  dplyr::rename("pos_first_deleted_hg38" = "hg38")

df_indels3 <- left_join(x = df_indels2, 
                        y = hg38_conversion_indels,
                        by = c("pos_last_deleted" = "pos")) %>% 
  dplyr::rename("pos_last_deleted_hg38" = "hg38")

df_indels4 <- df_indels3 %>% 
  select(variant, cigar, edit_string, indel_size, pos_first_deleted_hg38, pos_last_deleted_hg38, sge_region, r1_function_score, r2_function_score, function_score, pre, post, rL42_pre, rL42_post, lib, neg) %>% 
  dplyr::rename( "r1_day_5" = "pre",
                 "r2_day_5" = "rL42_pre",
                 "r1_day_14" = "post",
                 "r2_day_14" = "rL42_post",
                 "library" = "lib",
                 "negative" = "neg",
                 "mean_function_score" = "function_score"
                 )

#write.csv(x = df_indels4, file = paste0("~/Dropbox (The Francis Crick)/BRCA1_HAP1_HMEC_manuscript/SUPPLEMENTARY_TABLES/", date, "_ST2.csv"), row.names = FALSE) 

