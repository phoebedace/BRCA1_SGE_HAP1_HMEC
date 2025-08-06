#### Setup ####

#load packages
library(readxl)
library(patchwork)
library(scales)
library(plotROC)
library(PRROC)
library(gridExtra)
library(broom)
library(tidyverse)
library(ggpubr)
library(corrplot)
library(networkD3)
library(mclust)
library(ggstar)
library(htmlwidgets)
library(pROC)
library(networkD3)
library(htmlwidgets)

#setup info
date <- "250605"
setwd("/Users/dacep/Library/CloudStorage/OneDrive-TheFrancisCrickInstitute/R/HMEC_ALL_V2/")
options(scipen=999)
options(digits = 15)

##### General functions #####

#function to convert a list of columns to numeric
convert_to_numeric <- function(data, columns) {
  for (col in columns) {
    data[[col]] <- as.numeric(data[[col]])
  }
  return(data)
}

#### Importing HAP1 data ####
df_HAP1_new <- read.csv("~/Downloads/Run_SGE_scripts/250605_Phoebe_BRCA1_SGE_data_all_regions_newclinvar_nodups.csv")

#get rid of the old function classes and replace column name with new so I don't have to change all the names to "v2"
df_HAP1_new <- df_HAP1_new %>% 
  select(-HAP1_func_class) %>% 
  dplyr::rename("HAP1_func_class" = "HAP1_func_class_v2")

HAP1_neut <- df_HAP1_new %>% 
  filter(HAP1_func_class == "Neutral")

#Making list of neutral variants - only use variants from this list for null distribution for fdr calc
HAP1_neut_list <- as.character(HAP1_neut$cHGVS)

#### Import HMEC data ####
#Outputs from 'Process_pipeline_output_HMEC.R' script
df_x3_D7_D14 <- read.csv("R_script_output/240918_BRCA1x3_D7D14_rHMEC3rHMEC4_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x3_D7_D14o <- read.csv("R_script_output/240918_BRCA1x3_D7D14o_rHMEC3rHMEC4_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x3_D7_D21 <- read.csv("R_script_output/240918_BRCA1x3_D7D21_rHMEC3rHMEC4_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x3_D7_D21o <- read.csv("R_script_output/240918_BRCA1x3_D7D21o_rHMEC3rHMEC4_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17_D7_D14 <- read.csv("R_script_output/240920_BRCA1x17_D7D14_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17_D7_D14o <- read.csv("R_script_output/240920_BRCA1x17_D7D14o_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17_D7_D21 <- read.csv("R_script_output/240920_BRCA1x17_D7D21_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17_D7_D21o <- read.csv("R_script_output/240920_BRCA1x17_D7D21o_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x10d_D7_D14 <- read.csv("R_script_output/240918_BRCA1x10d_D7D14_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x10d_D7_D14o <- read.csv("R_script_output/240918_BRCA1x10d_D7D14o_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x10d_D7_D21 <- read.csv("R_script_output/240918_BRCA1x10d_D7D21_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x10d_D7_D21o <- read.csv("R_script_output/240918_BRCA1x10d_D7D21o_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x5_D7_D14 <- read.csv("R_script_output/240918_BRCA1x5_D7D14_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x5_D7_D21 <- read.csv("R_script_output/240918_BRCA1x5_D7D21_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x10a_D7_D14 <- read.csv("R_script_output/240918_BRCA1x10a_D7D14_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x10a_D7_D21 <- read.csv("R_script_output/240918_BRCA1x10a_D7D21_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x12a_D7_D14 <- read.csv("R_script_output/240918_BRCA1x12a_D7D14_rHMEC3rHMEC4_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x12a_D7_D14o <- read.csv("R_script_output/240918_BRCA1x12a_D7D14_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x12a_D7_D21 <- read.csv("R_script_output/240918_BRCA1x12a_D7D21_rHMEC3rHMEC4_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x12a_D7_D21o <- read.csv("R_script_output/240918_BRCA1x12a_D7D21_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x11_D7_D14 <- read.csv("R_script_output/240918_BRCA1x11alt_D7D14_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x11_D7_D14o <- read.csv("R_script_output/240918_BRCA1x11alt_D7D14o_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x11_D7_D21 <- read.csv("R_script_output/240918_BRCA1x11alt_D7D21_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x11_D7_D21o <- read.csv("R_script_output/240918_BRCA1x11alt_D7D21o_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17q_D7_D14 <- read.csv("R_script_output/240920_BRCA1x17q_D7D14_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17q_D7_D14o <- read.csv("R_script_output/240920_BRCA1x17q_D7D14o_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17q_D7_D21 <- read.csv("R_script_output/240920_BRCA1x17q_D7D21_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x17q_D7_D21o <- read.csv("R_script_output/240920_BRCA1x17q_D7D21o_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x6a_D7_D14 <- read.csv("R_script_output/240926_BRCA1x6a_D7D14_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x6a_D7_D14o <- read.csv("R_script_output/240926_BRCA1x6a_D7D14o_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x6a_D7_D21 <- read.csv("R_script_output/240926_BRCA1x6a_D7D21_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_x6a_D7_D21o <- read.csv("R_script_output/240926_BRCA1x6a_D7D21o_rHMEC1rHMEC2_tHDR_pos_merge_df_SYNNORM_CORRECTED.csv", header = TRUE)
df_u1_D7_D14 <- read.csv("R_script_output/241008_BRCA1u1_D7D14_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_u1_custom.csv", header = TRUE)
df_u1_D7_D14o <- read.csv("R_script_output/241008_BRCA1u1_D7D14o_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_u1_custom.csv", header = TRUE)
df_u1_D7_D21 <- read.csv("R_script_output/241008_BRCA1u1_D7D21_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_u1_custom.csv", header = TRUE)
df_u1_D7_D21o <- read.csv("R_script_output/241008_BRCA1u1_D7D21o_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_u1_custom.csv", header = TRUE)
df_i11_D7_D14 <- read.csv("R_script_output/241008_BRCA1i11_D7D14_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_i11_custom.csv", header = TRUE)
df_i11_D7_D14o <- read.csv("R_script_output/241008_BRCA1i11_D7D14o_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_i11_custom.csv", header = TRUE)
df_i11_D7_D21 <- read.csv("R_script_output/241008_BRCA1i11_D7D21_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_i11_custom.csv", header = TRUE)
df_i11_D7_D21o <- read.csv("R_script_output/241008_BRCA1i11_D7D21o_rHMEC5rHMEC6_tHDR_pos_merge_df_SYNNORM_CORRECTED_i11_custom.csv", header = TRUE)

#adding on identifiers so all dfs can be merged
df_x3_D7_D14_2 <- df_x3_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x3")
df_x3_D7_D14o_2 <- df_x3_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "x3")
df_x3_D7_D21_2 <- df_x3_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x3")
df_x3_D7_D21o_2 <- df_x3_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "x3")
df_x17_D7_D14_2 <- df_x17_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x17")
df_x17_D7_D14o_2 <- df_x17_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "x17")
df_x17_D7_D21_2 <- df_x17_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x17")
df_x17_D7_D21o_2 <- df_x17_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "x17")
df_x10d_D7_D14_2 <- df_x10d_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x10d")
df_x10d_D7_D14o_2 <- df_x10d_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "x10d")
df_x10d_D7_D21_2 <- df_x10d_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x10d")
df_x10d_D7_D21o_2 <- df_x10d_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "x10d")
df_x5_D7_D14_2 <- df_x5_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x5")
df_x5_D7_D21_2 <- df_x5_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x5")
df_x10a_D7_D14_2 <- df_x10a_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x10a")
df_x10a_D7_D21_2 <- df_x10a_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x10a")
df_x12a_D7_D14_2 <- df_x12a_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x12a")
df_x12a_D7_D14o_2 <- df_x12a_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "x12a")
df_x12a_D7_D21_2 <- df_x12a_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x12a")
df_x12a_D7_D21o_2 <- df_x12a_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "x12a")
df_x11_D7_D14_2 <- df_x11_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x11")
df_x11_D7_D14o_2 <- df_x11_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "x11")
df_x11_D7_D21_2 <- df_x11_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x11")
df_x11_D7_D21o_2 <- df_x11_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "x11")
df_x17q_D7_D14_2 <- df_x17q_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x17q")
df_x17q_D7_D14o_2 <- df_x17q_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "x17q")
df_x17q_D7_D21_2 <- df_x17q_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x17q")
df_x17q_D7_D21o_2 <- df_x17q_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "x17q")
df_x6a_D7_D14_2 <- df_x6a_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "x6a")
df_x6a_D7_D14o_2 <- df_x6a_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "x6a")
df_x6a_D7_D21_2 <- df_x6a_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "x6a")
df_x6a_D7_D21o_2 <- df_x6a_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "x6a")
df_u1_D7_D14_2 <- df_u1_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "u1")
df_u1_D7_D14o_2 <- df_u1_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "u1")
df_u1_D7_D21_2 <- df_u1_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "u1")
df_u1_D7_D21o_2 <- df_u1_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "u1")
df_i11_D7_D14_2 <- df_i11_D7_D14 %>% 
  mutate(expt = "D7D14",
         sge_region = "i11")
df_i11_D7_D14o_2 <- df_i11_D7_D14o %>% 
  mutate(expt = "D7D14o",
         sge_region = "i11")
df_i11_D7_D21_2 <- df_i11_D7_D21 %>% 
  mutate(expt = "D7D21",
         sge_region = "i11")
df_i11_D7_D21o_2 <- df_i11_D7_D21o %>% 
  mutate(expt = "D7D21o",
         sge_region = "i11")


#making all columns in all dfs identical
dfs <- list(df_x3_D7_D14_2, df_x3_D7_D14o_2, df_x3_D7_D21_2, df_x3_D7_D21o_2,
            df_x17_D7_D14_2, df_x17_D7_D14o_2, df_x17_D7_D21_2, df_x17_D7_D21o_2,
            df_x10d_D7_D14_2, df_x10d_D7_D14o_2, df_x10d_D7_D21_2, df_x10d_D7_D21o_2,
            df_x5_D7_D14_2, df_x5_D7_D21_2,
            df_x10a_D7_D14_2, df_x10a_D7_D21_2,
            df_x12a_D7_D14_2, df_x12a_D7_D14o_2, df_x12a_D7_D21_2, df_x12a_D7_D21o_2,
            df_x11_D7_D14_2, df_x11_D7_D14o_2, df_x11_D7_D21_2, df_x11_D7_D21o_2,
            df_x17q_D7_D14_2, df_x17q_D7_D14o_2, df_x17q_D7_D21_2, df_x17q_D7_D21o_2,
            df_x6a_D7_D14_2, df_x6a_D7_D14o_2, df_x6a_D7_D21_2, df_x6a_D7_D21o_2,
            df_u1_D7_D14_2, df_u1_D7_D14o_2, df_u1_D7_D21_2, df_u1_D7_D21o_2,
            df_i11_D7_D14_2, df_i11_D7_D14o_2, df_i11_D7_D21_2, df_i11_D7_D21o_2
            )

colClean <- function(x) {
  if (is.data.frame(x)) {
    # Check if 'sge_region' column exists in the data frame
    if ("sge_region" %in% colnames(x)) {
      pattern <- x[["sge_region"]]  # Accessing 'sge_region' as a vector
      colnames(x) <- gsub(pattern = pattern, replacement = "region", colnames(x))
      colnames(x) <- gsub(pattern = "rHMEC3", replacement = "rHMEC1", colnames(x))
      colnames(x) <- gsub(pattern = "rHMEC4", replacement = "rHMEC2", colnames(x))
      colnames(x) <- gsub(pattern = "rHMEC5", replacement = "rHMEC1", colnames(x))
      colnames(x) <- gsub(pattern = "rHMEC6", replacement = "rHMEC2", colnames(x))
    }
  }
  return(x)
}

cleaned_dfs <- lapply(dfs, colClean)

#merge all dfs
all_dfs <- bind_rows(cleaned_dfs)


#load hg38 conversion and rename columns ready to join dfs
hg38_conversion <- read.csv("databases/230420_BRCA1_hg19_hg38_conversion.csv")
hg38_conversion <- hg38_conversion %>% 
  dplyr::rename("Chrom" = "chr",
         "pos" = "hg19")


#add hg38 coords to dataframe
all_dfs <-
  left_join(x = all_dfs,
            y = hg38_conversion,
            by = c("Chrom", "pos"))

#add columns
all_dfs <- all_dfs %>% 
  #correct mis-numbering of exon 17 experiments
  mutate(protPos = ifelse(sge_region == "x17" | sge_region == "x17q", protPos - 21, protPos),
         CDSpos = ifelse(sge_region == "x17" | sge_region == "x17q", CDSpos - 63, CDSpos)) %>% 
  mutate(pHGVS = case_when(
    !is.na(Intron) ~ NA_character_,
    conseq == "5PRIME_UTR" ~ NA_character_,
    !is.na(Exon) ~ paste0("p.", oAA, protPos, nAA)
  )) %>% 
  #add reverse complement for Ref allele
  mutate(Ref_RC = case_when(
    Ref == "A" ~ "T",
    Ref == "C" ~ "G",
    Ref == "T" ~ "A",
    Ref == "G" ~ "C",
    Ref == "CC" ~ "GG"
  )) %>% 
mutate(rHMEC1_post_pre_function_score = log2(tHDR_post_pre_ratio_synnorm),
       rHMEC2_post_pre_function_score = log2(rL42_tHDR_post_pre_ratio_synnorm),
       diff_post_pre_function_scores = abs(rHMEC1_post_pre_function_score - rHMEC2_post_pre_function_score),
       post_pre_function_score = (rHMEC1_post_pre_function_score + rHMEC2_post_pre_function_score) / 2,
       SpliceAI_max = pmax(SpliceAI.acc.gain, SpliceAI.acc.loss, SpliceAI.don.gain, SpliceAI.don.loss)
)


#### Add ClinVar data ####
clinvar_download <- read.table(file = "databases/240913_clinvar.txt", sep="\t", header = TRUE, fill = TRUE, na.strings = "")

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
  mutate(Clinvar_Sep24_simple = case_when(grepl(pattern = "onflicting", x = Clinvar_interpretation_Sep24) ~ "Conflicting interpretations",
                                          grepl(pattern = "athogenic", x = Clinvar_interpretation_Sep24) ~ "Pathogenic/Likely pathogenic",
                                          grepl(pattern = "enign", x = Clinvar_interpretation_Sep24) ~ "Benign/Likely benign",
                                          grepl(pattern = "bsent", x = Clinvar_interpretation_Sep24) | grepl(pattern = "not provided", x = Clinvar_interpretation_Sep24) ~ "Absent",
                                          
                                          #  grepl(pattern = "ikely", x = Clinvar_interpretation_Sep24) & grepl(pattern = "athogenic", x = Clinvar_interpretation_Sep24) ~ "Likely pathogenic",
                                          TRUE ~ Clinvar_interpretation_Sep24),
         Clinvar_Sep24_histo = case_when(grepl(pattern = "ncertain", x = Clinvar_Sep24_simple) ~ "VUS/conflicting interpretations",
                                         grepl(pattern = "interpretations", x = Clinvar_Sep24_simple) ~ "VUS/conflicting interpretations",
                                         grepl(pattern = "athogenic", x = Clinvar_Sep24_simple) ~ "Pathogenic/Likely pathogenic",
                                         grepl(pattern = "enign", x = Clinvar_Sep24_simple) ~ "Benign/Likely benign",
                                         #  grepl(pattern = "ikely", x = Clinvar_interpretation_Sep24) & grepl(pattern = "athogenic", x = Clinvar_interpretation_Sep24) ~ "Likely pathogenic",
                                         TRUE ~ Clinvar_Sep24_simple),
         Clinvar_link = paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", VariationID))

#Add ClinVar data to HMEC data
all_dfs <- left_join(y = clinvar_download2,
                     x = all_dfs,
                     by = c("pos" = "GRCh37Location",
                            "Ref" = "Ref",
                            "alt" = "alt"))

all_dfs$Clinvar_interpretation_Sep24 <- ifelse(is.na(all_dfs$Clinvar_interpretation_Sep24), "Absent", all_dfs$Clinvar_interpretation_Sep24)
all_dfs$Clinvar_Sep24_histo <- ifelse(is.na(all_dfs$Clinvar_Sep24_histo), "Absent", all_dfs$Clinvar_Sep24_histo)
all_dfs$Clinvar_Sep24_histo <- factor(all_dfs$Clinvar_Sep24_histo, levels = c("Absent", "VUS/conflicting interpretations", "Benign/Likely benign", "Pathogenic/Likely pathogenic"))

#### Adding additional information ####
#distance to splice only goes to +/-20 into the intron - add in for introns that need longer
all_dfs <- all_dfs %>% 
  mutate(Dst2Splice = case_when(is.na(Dst2Splice) & sge_region == "x10d" ~ as.integer(pos-41243452),
                                is.na(Dst2Splice) & sge_region == "x3" & Intron == "3/23" ~ as.integer(pos-41267743),
                                is.na(Dst2Splice) & sge_region == "x3" & Intron == "2/23" ~ as.integer(41267796-pos), 
                                is.na(Dst2Splice) & sge_region == "x12a" & Intron == "11/23" ~ as.integer(41234592-pos),   
                                           TRUE ~ Dst2Splice))

#correcting information for x17 p.R1699Q
all_dfs <- all_dfs %>% 
  mutate(pHGVS = case_when(oAA == "R" & protPos == "1699" & nAA == "Q" ~ "p.R1699Q",
                           TRUE ~ pHGVS),
         hg38 = case_when(pos == 41215947 & protPos == "1699" & nAA == "Q" ~ 43063930,
                            TRUE ~ hg38),
  )


#adding cHGVS
all_dfs <- all_dfs %>% 
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
    (is.na(CDSpos)&(Intron == "11/23")&(sge_region == "x12a2")) ~ paste0("c.4186-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "12/23")&(sge_region == "x12b")) ~ paste0("c.4357+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "17/23")&(sge_region == "x17")) ~ paste0("c.5075-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "18/23")&(sge_region == "x17")) ~ paste0("c.5152+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "17/23")&(sge_region == "x17q")) ~ paste0("c.5075-", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    (is.na(CDSpos)&(Intron == "18/23")&(sge_region == "x17q")) ~ paste0("c.5152+", abs(as.numeric(as.character(Dst2Splice))), Ref_RC,">", rev_comp),
    !is.na(CDSpos) ~ paste0("c.", CDSpos, Ref_RC, ">", rev_comp)),
    canvar_link = paste0("https://canvaruk.org/result/BRCA1/", cHGVS)
    )


#### Calculating function scores ####

D14_all_dfs <- all_dfs %>%
  filter(expt == "D7D14") %>%
  filter(!is.na(cHGVS)) %>%
  mutate(D14_D7_PPR_r1 = tHDR_post_pre_ratio_synnorm,
         D14_D7_PPR_r2 = rL42_tHDR_post_pre_ratio_synnorm)

D21_all_dfs <- all_dfs %>%
  filter(expt == "D7D21") %>%
  filter(!is.na(cHGVS)) %>%
  mutate(D21_D7_PPR_r1 = tHDR_post_pre_ratio_synnorm,
         D21_D7_PPR_r2 = rL42_tHDR_post_pre_ratio_synnorm)

D14o_all_dfs <- all_dfs %>%
  filter(expt == "D7D14o") %>%
  filter(!is.na(cHGVS)) %>%
  mutate(D14_D7_PPR_r1 = tHDR_post_pre_ratio_synnorm,
         D14_D7_PPR_r2 = rL42_tHDR_post_pre_ratio_synnorm)

D21o_all_dfs <- all_dfs %>%
  filter(expt == "D7D21o") %>%
  filter(!is.na(cHGVS)) %>%
  mutate(D21_D7_PPR_r1 = tHDR_post_pre_ratio_synnorm,
         D21_D7_PPR_r2 = rL42_tHDR_post_pre_ratio_synnorm)

  
  ut_all_dfs_merge <- full_join(x = D14_all_dfs,
                                y = D21_all_dfs,
                                #by = "cHGVS")
                                by = c("cHGVS", "pHGVS", "sge_region", "pos", "hg38", "variant", "cigar", "pre", "rL42_pre", "lib", "Ref", "alt", "conseq", "oAA", "protPos", "nAA", "Intron", "Exon", "Clinvar_interpretation_Sep24", "Clinvar_Sep24_histo", "Clinvar_link", "canvar_link", "CDSpos", "Dst2Splice", "CADD.raw", "CADD.phred", "rev_comp", "Ref_RC", "Review.status", "Last_reviewed", "VariationID", "lib_pseudo_freq", "tHDR_pre_pseudo_freq", "rL42_tHDR_pre_pseudo_freq", "within_2bp_of_pam_edit"),
                                suffix = c(".D14", ".D21"))
  
  olap_all_dfs_merge <- full_join(x = D14o_all_dfs,
                                y = D21o_all_dfs,
                                #by = "cHGVS")
                                by = c("cHGVS", "pHGVS", "sge_region", "pos", "hg38", "variant", "cigar", "pre", "rL42_pre", "lib", "Ref", "alt", "conseq", "oAA", "protPos", "nAA", "Intron", "Exon", "Clinvar_interpretation_Sep24", "Clinvar_Sep24_histo", "Clinvar_link", "canvar_link", "CDSpos", "Dst2Splice", "CADD.raw", "CADD.phred", "rev_comp", "Ref_RC", "Review.status", "Last_reviewed", "VariationID", "lib_pseudo_freq", "tHDR_pre_pseudo_freq", "rL42_tHDR_pre_pseudo_freq", "within_2bp_of_pam_edit"),
                                suffix = c(".D14", ".D21"))
  
  #label dfs with treatment so they can be merged
  ut_all_dfs_merge$expt <- "ut"
  olap_all_dfs_merge$expt <- "olaparib"
  
#calculate function scores
  ut_all_dfs_merge2 <- ut_all_dfs_merge %>% 
    mutate(D21_D14_PPR_r1 = D21_D7_PPR_r1 / D14_D7_PPR_r1,
           D21_D14_PPR_r2 = D21_D7_PPR_r2 / D14_D7_PPR_r2,
           combined_FS_r1 = (log2(D21_D14_PPR_r1) + log2(D14_D7_PPR_r1))/2, #PPR = post/pre ratio
           combined_FS_r2 = (log2(D21_D14_PPR_r2) + log2(D14_D7_PPR_r2))/2,
           diff_r1_r2 = abs(combined_FS_r1-combined_FS_r2),
           mean_combined_fs = (combined_FS_r1+combined_FS_r2)/2)
  
  olap_all_dfs_merge2 <- olap_all_dfs_merge %>% 
    mutate(D21_D14_PPR_r1 = D21_D7_PPR_r1 / D14_D7_PPR_r1,
           D21_D14_PPR_r2 = D21_D7_PPR_r2 / D14_D7_PPR_r2,
           combined_FS_r1 = (log2(D21_D14_PPR_r1) + log2(D14_D7_PPR_r1))/2, #PPR = post/pre ratio
           combined_FS_r2 = (log2(D21_D14_PPR_r2) + log2(D14_D7_PPR_r2))/2,
           diff_r1_r2 = abs(combined_FS_r1-combined_FS_r2),
           mean_combined_fs = (combined_FS_r1+combined_FS_r2)/2)
  
  #changing name of final function score to avoid needing to change later 
  ut_all_dfs_merge2$post_pre_function_score <- ut_all_dfs_merge2$mean_combined_fs
  olap_all_dfs_merge2$post_pre_function_score <- olap_all_dfs_merge2$mean_combined_fs
  
  #join dfs for filtering
  all_all_dfs_merge2 <- rbind(ut_all_dfs_merge2, olap_all_dfs_merge2)

#### Filtering ####
##### filter variants affected by PAM edits ##### 
all_dfs_filtered_1 <- all_all_dfs_merge2 %>% 
  filter(pos != 41243496) %>%  #filter out x10d variants between 2 PAM edits: (labelled as missense but without PAM edits would be nonsense)
  filter(!(pos == 41256962 & alt == "A")) %>%  #x5 variant that in presence of PAM site introduces splice site
  filter(!(pos == 41246779 & alt == "A")) %>% #x10a variant that in absence of PAM would introduce splice site with high spliceAI score.
  filter(!(pos == 41215948 & alt == "A")) %>% #x17
  filter(!(pos == 41215947 & alt == "T")) %>%  #x17
  filter(pos != 41267770) %>% #x3 PAM
  filter(!is.na(cHGVS)) #filtering out variants that don't actually exist because PAM site messes them up - cHGVS didn't get added as incorrect


##### filter based on lib + day 7 freq #####
all_dfs_filtered_freq <- all_dfs_filtered_1 %>% 
  filter(lib_pseudo_freq > 0.0001) %>%  #filters 31/36 variants?
  filter(tHDR_pre_pseudo_freq >= 0.00001 & rL42_tHDR_pre_pseudo_freq >= 0.00001) #%>%  #filters 1354 variants
  
##### filter variants scoring discordantly between reps #####
all_dfs_filtered_2 <- all_dfs_filtered_freq %>% 
  filter((diff_r1_r2 < 2) | (diff_r1_r2 >= 2 & ((combined_FS_r1 < -1 & combined_FS_r2 < -1) | (combined_FS_r1 > -1 & combined_FS_r2 > -1)))) 

      
#### Are variants significantly depleted? #### 
add_fdr_info <- function(region) {
  global_median_syn <- median(df_input[which(df_input$conseq == 'SYNONYMOUS'),c('post_pre_function_score')])
  global_median_ns <- median(df_input[which(df_input$conseq == 'STOP_GAINED'),c('post_pre_function_score')])
  BRCA1_region_final_df <- df_input[which(df_input$sge_region == region),] 
  BRCA1_region_function_score_median_syn <- median(BRCA1_region_final_df[which(BRCA1_region_final_df$conseq == 'SYNONYMOUS'),c('post_pre_function_score')])
  BRCA1_region_function_score_median_ns <- median(BRCA1_region_final_df[which(BRCA1_region_final_df$conseq == 'STOP_GAINED'),c('post_pre_function_score')])
  BRCA1_region_final_df$function_score_sns <- (BRCA1_region_final_df$post_pre_function_score/(BRCA1_region_function_score_median_syn-BRCA1_region_function_score_median_ns))*(global_median_syn-global_median_ns)
  BRCA1_region_null_fss <- BRCA1_region_final_df[which(BRCA1_region_final_df$conseq == "SYNONYMOUS" & (BRCA1_region_final_df$cHGVS %in% HAP1_neut_list | BRCA1_region_final_df$sge_region == "x3")),]$function_score_sns
  BRCA1_region_null_mean <- mean(BRCA1_region_null_fss)
  BRCA1_region_null_sd <- sd(BRCA1_region_null_fss)
  BRCA1_region_final_df$pvalues <- pnorm(BRCA1_region_final_df$function_score_sns, BRCA1_region_null_mean, sd=BRCA1_region_null_sd) 
  BRCA1_region_final_df$fdr <- p.adjust(BRCA1_region_final_df$pvalues, method = "BH") # Calculate FDR for each p-value
  BRCA1_region_final_df <- BRCA1_region_final_df %>%
    mutate(
      fs_sig = case_when(
        fdr > 0.5 ~ "FDR > 0.5",
        fdr < 0.5 & fdr >= 0.2 ~ "FDR < 0.5",
        fdr < 0.2 & fdr >= 0.05 ~ "FDR < 0.2",
        fdr < 0.05 & fdr >= 0.01 ~ "FDR < 0.05",
        fdr < 0.01 ~ "FDR < 0.01",
        TRUE ~ "FDR NA"
      )
    )
  
  if (any(BRCA1_region_final_df$fs_sig == "FDR < 0.01")) {
    # Calculate max0.01 and min0.05
    max0.01 <- max(BRCA1_region_final_df[BRCA1_region_final_df$fs_sig == "FDR < 0.01", ]$function_score_sns)
    min0.05 <- min(BRCA1_region_final_df[BRCA1_region_final_df$fs_sig != "FDR < 0.01", ]$function_score_sns)
    
    # Calculate fs_threshold (based on syn vars)
    BRCA1_region_final_df$fs_threshold <- (max0.01 + min0.05) / 2
  } else {
    # If no values where fs_sig == "FDR < 0.01", set max0.01, min0.05, and fs_threshold to the minimum function score
    max0.01 <- NA
    min0.05 <- NA
    BRCA1_region_final_df$fs_threshold <- NA
  }
  return(BRCA1_region_final_df)
}

regions <- unique(all_dfs_filtered_2$sge_region[which(all_dfs_filtered_2$sge_region != "u1" & all_dfs_filtered_2$sge_region != "i11" & all_dfs_filtered_2$expt == "ut")])
df_input <- all_dfs_filtered_2 %>%  filter(expt == "ut")
BRCA1_region_final_df <- data.frame()  # Initialize as an empty data frame
df_ut <- data.frame()

for (region in regions) {
  output <- add_fdr_info(region)
  df_ut <- rbind(df_ut, output) #NB this is coding regions only
}


#do same for olap data
regions <- unique(all_dfs_filtered_2$sge_region[which(all_dfs_filtered_2$sge_region != "u1" & all_dfs_filtered_2$sge_region != "i11" & all_dfs_filtered_2$expt == "olaparib")])
df_input <- all_dfs_filtered_2 %>%  filter(expt == "olaparib")
BRCA1_region_final_df <- data.frame()  # Initialize as an empty data frame
df_olap <- data.frame()

for (region in regions) {
  output <- add_fdr_info(region)
  df_olap <- rbind(df_olap, output) #NB this is coding regions only
}

add_fdr_info_nc <- function(region) {
  BRCA1_region_final_df <- df_input[which(df_input$sge_region == region),] 
  BRCA1_region_final_df$function_score_sns <- BRCA1_region_final_df$post_pre_function_score #(BRCA1_region_final_df$post_pre_function_score/(BRCA1_region_function_score_median_syn-BRCA1_region_function_score_median_ns))*(median_syn_global-median_non_global)
  BRCA1_region_null_fss <- BRCA1_region_final_df[which(BRCA1_region_final_df$cHGVS %in% HAP1_neut_list),]$function_score_sns
  
  BRCA1_region_null_mean <- mean(BRCA1_region_null_fss)
  BRCA1_region_null_sd <- sd(BRCA1_region_null_fss)
  BRCA1_region_final_df$pvalues <- pnorm(BRCA1_region_final_df$function_score_sns, BRCA1_region_null_mean, sd=BRCA1_region_null_sd) 
  BRCA1_region_final_df$fdr <- p.adjust(BRCA1_region_final_df$pvalues, method = "BH") # Calculate FDR for each p-value
  BRCA1_region_final_df <- BRCA1_region_final_df %>%
    mutate(
      fs_sig = case_when(
        fdr > 0.5 ~ "FDR > 0.5",
        fdr < 0.5 & fdr >= 0.2 ~ "FDR < 0.5",
        fdr < 0.2 & fdr >= 0.05 ~ "FDR < 0.2",
        fdr < 0.05 & fdr >= 0.01 ~ "FDR < 0.05",
        fdr < 0.01 ~ "FDR < 0.01",
        TRUE ~ "FDR NA"
      )
    )
  
  if (any(BRCA1_region_final_df$fs_sig == "FDR < 0.01")) {
    # Calculate max0.01 and min0.05
    max0.01 <- max(BRCA1_region_final_df[BRCA1_region_final_df$fs_sig == "FDR < 0.01", ]$function_score_sns)
    min0.05 <- min(BRCA1_region_final_df[BRCA1_region_final_df$fs_sig != "FDR < 0.01", ]$function_score_sns)
    
    # Calculate fs_threshold (based on syn vars)
    BRCA1_region_final_df$fs_threshold <- (max0.01 + min0.05) / 2
  } else {
    # If no values where fs_sig == "FDR < 0.01", set max0.01, min0.05, and fs_threshold to the minimum function score
    max0.01 <- NA
    min0.05 <- NA
    BRCA1_region_final_df$fs_threshold <- NA
  }
  
  return(BRCA1_region_final_df)
}

#non-coding, untreated
regions <- c("i11", "u1")
df_input <- all_dfs_filtered_2 %>%  filter(expt == "ut")
BRCA1_region_final_df <- data.frame()  # Initialize as an empty data frame
df_nc_ut <- data.frame()

for (region in regions) {
  output <- add_fdr_info_nc(region)
  df_nc_ut <- rbind(df_nc_ut, output) #NB this is NON-coding regions only
}

#non-coding, olaparib treated
df_input <- all_dfs_filtered_2 %>%  filter(expt == "olaparib")
BRCA1_region_final_df <- data.frame()  # Initialize as an empty data frame
df_nc_olap <- data.frame()

for (region in regions) {
  output <- add_fdr_info_nc(region)
  df_nc_olap <- rbind(df_nc_olap, output) #NB this is NON-coding regions only
}


###### Binding dfs back together ###### 
all_dfs_final <- rbind(df_ut, df_nc_ut, df_olap, df_nc_olap)
all_dfs_final <- all_dfs_final %>% 
  mutate(HMEC_class = case_when(fdr < 0.01 ~ "LoF",
                                TRUE ~ "Neutral"))


#### Merge x3 HMEC with 2018 HAP1 data ####

#The only SGE region I have HMEC but not new HAP1 data for is x3 - get from Findlay et al. 2018
df_HAP1_old <- read_excel("HAP1_data/Supplementary_Table_1_revised_w_WT.xlsx",
                          sheet = "Sheet1",
                          range = cell_limits(c(4, 1), c(NA, NA)),
                          na = c("", "NA", "n/a"))
df_HAP1_old_x3 <- df_HAP1_old %>% 
  filter(experiment == "X3") %>% 
  dplyr::rename("pos" = "position (hg19)",
                "Chrom" = "chromosome",
                "Ref" = "reference")


# modifying old HAP1 df to match new
df_HAP1_old_x3_renamed <- df_HAP1_old_x3 %>% 
  select(pos, Ref, alt, experiment, consequence, transcript_variant, aa_pos, aa_ref, aa_alt, function.score.mean, function.score.r1, function.score.r2, func.class, mean.rna.score) %>% 
  mutate(pHGVS = case_when(!is.na(aa_pos) ~ paste0("p.", aa_ref, aa_pos, aa_alt),
                           TRUE ~ NA_character_),
         func.class = case_when(func.class == "LOF" ~ "LoF",
                                func.class == "INT" ~ "Intermediate",
                                func.class == "FUNC" ~ "Neutral"),
         experiment = case_when(experiment == "X3" ~ "x3"),
         HAP1_FDR = NA_character_) %>% 
  select(pos, Ref, alt, experiment, consequence, transcript_variant, pHGVS, function.score.mean, function.score.r1, function.score.r2, func.class, mean.rna.score, HAP1_FDR) %>% 
  mutate(conseq = case_when(consequence == "Intronic" ~ "INTRONIC",
                            consequence == "Synonymous" ~ "SYNONYMOUS",
                            consequence == "Missense" ~ "NON_SYNONYMOUS",
                            consequence == "Splice region" ~ "SPLICE_SITE",
                            consequence == "Canonical splice" ~ "CANONICAL_SPLICE",
                            consequence == "Nonsense" ~ "STOP_GAINED")) %>% 
  mutate(fs_threshold_upper = -0.748,
         fs_threshold_lower2 = -1.328) %>% 
  select(-consequence)
names(df_HAP1_old_x3_renamed) <- c("pos", "Ref", "alt", "sge_region", "cHGVS", "pHGVS", "HAP1_function_score_mean", "HAP1_function_score_r1", "HAP1_function_score_r2", "HAP1_func_class", "HAP1_mean_RNA_score", "HAP1_FDR", "conseq", "fs_threshold_upper", "fs_threshold_lower2")
  df_HAP1_old_x3_renamed$HAP1_data_source <- "old"

#merge HMEC df with Greg's HAP1 scores - x3 only
df_HMEC_old_HAP1_merge_x3 <- inner_join(x = all_dfs_final,
                                     y = df_HAP1_old_x3_renamed,
                                     by = c("pos", "alt"),
                                     suffix = c(".HMEC", ".HAP1"))

#modifying new HAP1 df
df_HAP1_new_renamed <- df_HAP1_new %>%
  select(pos, Ref, alt, sge_region, cHGVS, pHGVS, function_score_sns, rL41_post_pre_function_score, rL42_post_pre_function_score, HAP1_func_class, RNA_score, fdr, conseq, fs_threshold_upper, fs_threshold_lower2)
names(df_HAP1_new_renamed) <- c("pos", "Ref", "alt", "sge_region", "cHGVS", "pHGVS", "HAP1_function_score_mean", "HAP1_function_score_r1", "HAP1_function_score_r2", "HAP1_func_class", "HAP1_mean_RNA_score", "HAP1_FDR", "conseq", "fs_threshold_upper", "fs_threshold_lower2")
df_HAP1_new_renamed$HAP1_data_source <- "new"

#binding the new and old dataframes together
df_HAP1_old_and_new <- rbind(df_HAP1_new_renamed, df_HAP1_old_x3_renamed)


#### Calculating mean scores for exon 17 (4 reps) ####
# 4 reps for exon 17 in 2 experiments (x17 and x17q)
#x17q also contains p.R1699Q variant

df_x17x17q <- all_dfs_final %>% 
  filter(sge_region == "x17" | sge_region == "x17q")

#this function takes the x17 and x17q data and creates a mean score for variants where I have 2 scores
#otherwise uses the one score I already have
#re-calculates fdr based on the new scores
x17_expt_merge <- function(experiment, df) {
  df_expt <- df[which(df$expt == experiment),]
  df_expt <- df_expt %>% 
    #select which columns I want - can only be columns that have same value regardless of whether x17/x17q
    select(pos, hg38, alt, Ref, oAA, nAA, protPos, pHGVS, cHGVS, expt, sge_region, conseq, Clinvar_interpretation_Sep24, Clinvar_Sep24_histo, Clinvar_link, canvar_link, Intron, Exon, CDSpos, Dst2Splice, CADD.raw, CADD.phred, rev_comp, Ref_RC, within_2bp_of_pam_edit, Review.status, Last_reviewed, VariationID, function_score_sns, fdr, fs_sig, fs_threshold, combined_FS_r1, combined_FS_r2, HMEC_class, SpliceAI_max.D14)
  variables = c("function_score_sns", "combined_FS_r1", "combined_FS_r2", "fdr", "fs_sig", "fs_threshold", "HMEC_class")
  df_x17x17q_spread <- df_expt %>%
    gather("property", "value", dplyr::all_of(variables)) %>%
    mutate(sge_region_property = paste(sge_region, property, sep = "_")) %>%
    select(-property, -sge_region) %>%
    spread(sge_region_property, value)
  
  #convert specific columns to numeric
  columns_to_convert <- c("x17q_function_score_sns", "x17_function_score_sns", "x17q_combined_FS_r1", "x17q_combined_FS_r2", "x17_combined_FS_r1", "x17_combined_FS_r2")
  df_x17x17q_spread <- convert_to_numeric(df_x17x17q_spread, columns_to_convert)
  
  #calculate mean function score
  df_x17x17q_spread_mean <- df_x17x17q_spread %>% 
    dplyr::mutate(mean_function_score_sns = case_when(!is.na(x17q_function_score_sns) & !is.na(x17_function_score_sns) ~ (x17q_function_score_sns+x17_function_score_sns)/2,
                                                      !is.na(x17q_function_score_sns) & is.na(x17_function_score_sns) ~ x17q_function_score_sns,
                                                      is.na(x17q_function_score_sns) & !is.na(x17_function_score_sns) ~ x17_function_score_sns)) 

  #add new fdrs
  BRCA1_region_final_df <- df_x17x17q_spread_mean
  BRCA1_region_null_fss <- BRCA1_region_final_df[which(BRCA1_region_final_df$conseq == "SYNONYMOUS" & BRCA1_region_final_df$cHGVS != "c.5130A>G" & BRCA1_region_final_df$cHGVS != "c.5127A>G"),]$mean_function_score_sns
  BRCA1_region_null_mean <- mean(BRCA1_region_null_fss)
  BRCA1_region_null_sd <- sd(BRCA1_region_null_fss)
  BRCA1_region_final_df$pvalues <- pnorm(BRCA1_region_final_df$mean_function_score_sns, BRCA1_region_null_mean, sd=BRCA1_region_null_sd) 
  BRCA1_region_final_df$fdr <- p.adjust(BRCA1_region_final_df$pvalues, method = "BH") # Calculate FDR for each p-value
  BRCA1_region_final_df <- BRCA1_region_final_df %>%
    mutate(
      fs_sig = case_when(
        fdr > 0.5 ~ "FDR > 0.5",
        fdr < 0.5 & fdr >= 0.2 ~ "FDR < 0.5",
        fdr < 0.2 & fdr >= 0.05 ~ "FDR < 0.2",
        fdr < 0.05 & fdr >= 0.01 ~ "FDR < 0.05",
        fdr < 0.01 ~ "FDR < 0.01",
        TRUE ~ "FDR NA"
      )
    )
  
  if (any(BRCA1_region_final_df$fs_sig == "FDR < 0.01")) {
    # Calculate max0.01 and min0.05
    max0.01 <- max(BRCA1_region_final_df[BRCA1_region_final_df$fs_sig == "FDR < 0.01", ]$mean_function_score_sns)
    min0.05 <- min(BRCA1_region_final_df[BRCA1_region_final_df$fs_sig != "FDR < 0.01", ]$mean_function_score_sns)
    
    # Calculate fs_threshold (based on syn vars)
    BRCA1_region_final_df$fs_threshold <- (max0.01 + min0.05) / 2
  } else {
    # If no values where fs_sig == "FDR < 0.01", set max0.01, min0.05, and fs_threshold to the minimum function score
    max0.01 <- NA
    min0.05 <- NA
    BRCA1_region_final_df$fs_threshold <- NA
  }
  
  #add new columns
  BRCA1_region_final_df <- BRCA1_region_final_df %>% 
    mutate(mean_HMEC_class = case_when(fdr < 0.01 ~ "LoF",
                                       TRUE ~ "Neutral"),
           mismatch = case_when(x17_HMEC_class != x17q_HMEC_class ~ TRUE,
                                TRUE ~ FALSE),
           diff_x17_x17q = abs(x17_function_score_sns - x17q_function_score_sns))
  
  
  return(BRCA1_region_final_df)
}

df_x17x17q_output <- data.frame()
expts <- c("ut", "olaparib")

#run function for each experiment
for (expt in expts) {
  output <- x17_expt_merge(expt, df_x17x17q)
  df_x17x17q_output <- rbind(df_x17x17q_output, output)
}


#Apply filters
df_x17x17q_output_filtered_original <- df_x17x17q_output %>%
  filter(((diff_x17_x17q < 2) | (diff_x17_x17q >= 2 & ((x17_function_score_sns < -1 & x17q_function_score_sns < -1) | (x17_function_score_sns > -1 & x17q_function_score_sns > -1)))) | is.na(diff_x17_x17q))
 
df_x17x17q_output_filtered <- df_x17x17q_output_filtered_original


##### Combine mean x17 data (17m) with rest of data #####

#Get x17 df into same format as the all_dfs_final df below:
df_x17x17q_output_filtered_reformatted <- df_x17x17q_output_filtered %>% 
  dplyr::rename(HMEC_class = mean_HMEC_class,
                function_score_sns = mean_function_score_sns,
                combined_FS_r1 = x17_function_score_sns,
                combined_FS_r2 = x17q_function_score_sns,
  ) %>% 
  mutate(sge_region = "x17m") %>% #m = mean
  select(pos, hg38, alt, Ref, oAA, nAA, protPos, pHGVS, cHGVS, expt, sge_region, conseq, Clinvar_interpretation_Sep24, Clinvar_Sep24_histo, Clinvar_link, canvar_link, Intron, Exon, CDSpos, Dst2Splice, CADD.raw, CADD.phred, rev_comp, Ref_RC, within_2bp_of_pam_edit, Review.status, Last_reviewed, VariationID, function_score_sns, fdr, fs_sig, fs_threshold, combined_FS_r1, combined_FS_r2, HMEC_class, SpliceAI_max.D14) %>% 
  mutate(LoF = case_when(fs_sig == "FDR < 0.01" ~ TRUE,
                         TRUE ~ FALSE))

all_dfs_final_select <- all_dfs_final %>% 
  #get rid of non-meaned x17 data
  filter(sge_region != "x17" & sge_region != "x17q") %>% 
  select(pos, hg38, alt, Ref, oAA, nAA, protPos, pHGVS, cHGVS, expt, sge_region, conseq, Clinvar_interpretation_Sep24, Clinvar_Sep24_histo, Clinvar_link, canvar_link, Intron, Exon, CDSpos, Dst2Splice, CADD.raw, CADD.phred, rev_comp, Ref_RC, within_2bp_of_pam_edit, Review.status, Last_reviewed, VariationID, function_score_sns, fdr, fs_sig, fs_threshold, combined_FS_r1, combined_FS_r2, HMEC_class, SpliceAI_max.D14) %>% 
  mutate(LoF = case_when(fs_sig == "FDR < 0.01" ~ TRUE,
                         TRUE ~ FALSE))

  
#join with mean x17 data (x17m)
all_dfs_final_select_incl_x17 <- rbind(all_dfs_final_select, df_x17x17q_output_filtered_reformatted)

#renaming to avoid having to rename throughout script
df_HMEC_spread <- all_dfs_final_select_incl_x17


#### Modifying HAP1 dataframe to allow merging with HMEC data ####
#rename x17q to x17m in HAP1 df so that I can merge with HMEC
df_HAP1_old_and_new2 <- df_HAP1_old_and_new %>% 
  mutate(sge_region = case_when(sge_region == "x17q" ~ "x17m",
                                TRUE ~ sge_region))

#Normalising Findlay et al. exon 3 data to new HAP1 data
#New HAP1 data median syn/non
new_HAP1_med_score_NON <- df_HAP1_old_and_new2 %>% 
  filter(sge_region != "x3") %>% 
  filter(conseq == "STOP_GAINED") %>% 
  group_by(conseq) %>% 
  summarise(median_fs = median(HAP1_function_score_mean)) %>% 
  pull()

new_HAP1_med_score_SYN <- df_HAP1_old_and_new2 %>% 
  filter(sge_region != "x3") %>% 
  filter(conseq == "SYNONYMOUS") %>% 
  group_by(conseq) %>% 
  summarise(median_fs = median(HAP1_function_score_mean)) %>% 
  pull()

#x3 2018 median syn/non
old_x3_HAP1_med_score_NON <- df_HAP1_old_and_new2 %>% 
  filter(sge_region == "x3") %>% 
  filter(conseq == "STOP_GAINED") %>% 
  group_by(conseq) %>% 
  summarise(median_fs = median(HAP1_function_score_mean)) %>% 
  pull()

old_x3_HAP1_med_score_SYN <- df_HAP1_old_and_new2 %>% 
  filter(sge_region == "x3") %>% 
  filter(conseq == "SYNONYMOUS") %>% 
  group_by(conseq) %>% 
  summarise(median_fs = median(HAP1_function_score_mean)) %>% 
  pull()

x3_normalised <- df_HAP1_old_and_new2 %>% 
  filter(sge_region == "x3") %>% 
  mutate(HAP1_function_score_mean = HAP1_function_score_mean/(old_x3_HAP1_med_score_SYN-old_x3_HAP1_med_score_NON)*(new_HAP1_med_score_SYN-new_HAP1_med_score_NON))

df_HAP1_old_and_new2_x3norm <- rbind(df_HAP1_old_and_new2 %>% 
                                       filter(sge_region != "x3"),
                                     x3_normalised)


#keep vars with scores in both HAP1 and HMEC (at least one HMEC score)
df_HMEC_HAP1_all_spread_ut <- inner_join(x = df_HMEC_spread %>% filter(expt == "ut"),
                                     y = df_HAP1_old_and_new2_x3norm,
                                     by = c("pos", "Ref", "alt", "cHGVS", "pHGVS", "conseq")) 


df_HMEC_HAP1_all_spread_olap <- inner_join(x = df_HMEC_spread %>% filter(expt == "olaparib"),
                                         y = df_HAP1_old_and_new2_x3norm,
                                         by = c("pos", "Ref", "alt",  "cHGVS", "pHGVS", "conseq")) 



#combine ut and olaparib dataframes
df_HMEC_HAP1_all_spread <- rbind(df_HMEC_HAP1_all_spread_ut, df_HMEC_HAP1_all_spread_olap)


  #### Assigning combined function classes ####
df_HMEC_HAP1_all_spread2 <- df_HMEC_HAP1_all_spread %>% 
  mutate(sge_region = sge_region.x) %>% 
  mutate(conseq = case_when(conseq == "5PRIME_UTR" & hg38 > 43125364 ~ "PROMOTER",
                            TRUE ~ conseq)) %>% 
  mutate(final_func_class = case_when(HAP1_func_class == "LoF" & HMEC_class == "LoF" ~ "LoF",
                                      HAP1_func_class == "Neutral" & HMEC_class == "Neutral" ~ "Neutral",
                                      HAP1_func_class == "LoF" & HMEC_class == "Neutral" ~ "HAP1 LoF, HMEC neutral",
                                      HAP1_func_class == "Intermediate" & HMEC_class == "LoF" ~ "HAP1 int, HMEC LoF",
                                      HAP1_func_class == "Intermediate" & HMEC_class == "Neutral" ~ "HAP1 int, HMEC neutral",
                                      HAP1_func_class == "Neutral" & HMEC_class == "LoF" ~ "HAP1 neutral, HMEC LoF"
                                      ))
df_HMEC_HAP1_all_spread2$final_func_class <- factor(df_HMEC_HAP1_all_spread2$final_func_class, levels = c("Neutral", "HAP1 int, HMEC neutral", "HAP1 LoF, HMEC neutral", "LoF", "HAP1 neutral, HMEC LoF", "HAP1 int, HMEC LoF"))

df_HMEC_HAP1_all_spread2$conseq <- factor(df_HMEC_HAP1_all_spread2$conseq, levels = c("INTRONIC", "5PRIME_UTR", "PROMOTER", "NON_SYNONYMOUS", "SPLICE_SITE", "CANONICAL_SPLICE", "STOP_GAINED", "SYNONYMOUS"))


#### Adding additional data ####
##### gnomAD #####

gnomadv4 <- read.csv("databases/gnomAD_v4.1.0_17-43044295-43170245_2025_02_12_16_51_25.csv")

df_HMEC_HAP1_all_spread3 <- left_join(df_HMEC_HAP1_all_spread2,
                                  gnomadv4,
                                  by = c("hg38" = "Position",
                                         "Ref" = "Reference",
                                         "alt" = "Alternate")
)


##### FoldX #####
foldX_RING <- read.csv("foldx/240124_FoldX_output_RING_domain.csv")
foldX_CC <- read.csv("foldx/240126_FoldX_output_CC_domain.csv")
heptad_positions <- read.csv("databases/cc_domain_heptad_repeats.csv")
foldX_BRCT <- read.csv("foldx/240126_FoldX_output_BRCT_domain.csv")

foldX_BRCT2 <- foldX_BRCT %>% 
  mutate(domain = "BRCT")
foldX_RING2 <- foldX_RING %>% 
  mutate(domain = "RING")
foldX_CC2 <- foldX_CC %>% 
  mutate(domain = "Coiled-coil")

foldX_all <- rbind(foldX_BRCT2, foldX_CC2, foldX_RING2)

foldX_all_join <- left_join(df_HMEC_HAP1_all_spread2,
          foldX_all,
          by = "pHGVS") %>% 
  filter(!is.na(FoldX_score))


##### Alpha Missense #####
alphamissense <- read_tsv("databases/AlphaMissense-Search-P38398_2.tsv")
alphamissense <- alphamissense %>% 
  mutate(pHGVS = paste0("p.", a.a.1, position, a.a.2)) %>% 
  select(pHGVS, position, a.a.2, `pathogenicity score`, `pathogenicity class`)


df_AM_merge_all <- inner_join(x = df_HMEC_HAP1_all_spread2,
                             y = alphamissense,
                             by = c("pHGVS"))

  
##### EVE scores #####
  eve_scores <- read.csv("databases/241024_EVE_scores_BRCA1_HUMAN.csv", na.strings = c(""))
  
  df_HMEC_HAP1_all_spread2_protein <- df_HMEC_HAP1_all_spread2 %>% 
    filter(expt == "ut")
  
  dim(df_HMEC_HAP1_all_spread2_protein) 
  
eve_merge <- left_join(x = df_HMEC_HAP1_all_spread2_protein,
                       y = eve_scores,
                       by = c("protPos" = "position",
                              "oAA" = "wt_aa",
                              "nAA" = "mt_aa"))
  
##### All of Us #####
AOU_data <- read.csv("~/Downloads/Run_SGE_scripts/240626_AllOfUs_BRCA1.csv", na.strings = "")

df_AOU_ut <- df_HMEC_HAP1_all_spread2 %>% 
  filter(expt == "ut") %>% 
  mutate(variantId = paste("17", hg38, Ref, alt, sep = "-"))

df_AOU2 <- left_join(x = df_AOU_ut,
                     y = AOU_data,
                     by = "variantId") %>% 
  mutate(present_in_AOU = case_when(is.na(alleleCount) ~ FALSE,
                                    TRUE ~ TRUE))

df_present_in_AOU <- df_AOU2 %>% 
  group_by(final_func_class, present_in_AOU) %>% 
  summarise(n = n()) %>% 
  mutate(pct = n/sum(n)*100)
df_present_in_AOU$final_func_class <- factor(df_present_in_AOU$final_func_class, levels = c("Neutral", "HAP1 int, HMEC neutral", "HAP1 LoF, HMEC neutral", "LoF", "HAP1 neutral, HMEC LoF"))


##### Comparison to other functional assays #####
#From Lyra et al.
df_other_assays <- read_excel("~/OneDrive - The Francis Crick Institute/BRCA1_assays_Lyra_et_al.xlsx",
                          sheet = "STable 1",
                          #range = cell_limits(c(4, 1), c(NA, NA)),
                          na = c("", "NA", "n/a")) %>% 
  dplyr::rename("pHGVS" = "T5",
                "exon" = "T1",
                "oAA" = "T2",
                "protPos" = "T3",
                "nAA" = "T4",
                )
#making sure columns get read in as 1/0, not T/F
df_other_assays[, 10:140] <- lapply(df_other_assays[, 10:140], as.numeric) 

#exclude T131 (Findlay et al. 2018 assay)
df_other_assays2_T131 <- df_other_assays %>% 
  select(-T131) %>% 
  mutate(
    count_zeros = rowSums(select(., 10:139) == 0, na.rm = TRUE),
    count_ones = rowSums(select(., 10:139) == 1, na.rm = TRUE))


df_sge_data_and_other_assays_T131 <- inner_join(x = df_HMEC_HAP1_all_spread2 %>% filter(expt == "ut"),
                                                y = df_other_assays2_T131,
                                                by = c("pHGVS", "oAA", "nAA", "protPos")) %>% 
  mutate(diff_ones_zeros = count_ones - count_zeros,
         total_assays = count_ones + count_zeros,
         pct_path = count_ones/total_assays*100,
         pct_benign = count_zeros/total_assays*100,
         final_func_class_labs = str_replace_all(final_func_class, ",", ",\n")
  ) %>% 
  filter(final_func_class != "HAP1 neutral, HMEC LoF")
df_sge_data_and_other_assays_T131$final_func_class_labs <- factor(df_sge_data_and_other_assays_T131$final_func_class_labs, levels = c("Neutral", "HAP1 int,\n HMEC neutral", "HAP1 LoF,\n HMEC neutral",  "LoF")) 
df_sge_data_and_other_assays_T131_2 <- df_sge_data_and_other_assays_T131 %>% 
  filter(!(count_zeros == 0 & count_ones == 0)) #get rid of vars where no other assay data exists


#### Defining truthset and GMM ####
df_truthset <- df_HMEC_HAP1_all_spread2 %>% 
  filter(expt == "ut") %>% 
  filter(Clinvar_Sep24_histo != "VUS/conflicting interpretations" & Clinvar_Sep24_histo != "Absent") %>% 
  filter(Review.status != "criteria provided, single submitter" & Review.status != "no assertion criteria provided") %>% 
  #Remove benign variants with a high spliceAI score
  filter(!(Clinvar_Sep24_histo == "Benign/Likely benign" & SpliceAI_max.D14 >0.2))


#Mixture modelling
df_truthset_gmm <- df_truthset %>% 
  mutate(numeric_label = case_when(Clinvar_Sep24_histo == "Benign/Likely benign" ~ 1,
                                   Clinvar_Sep24_histo == "Pathogenic/Likely pathogenic" ~ 0
  )) %>% 
  select(cHGVS, numeric_label, function_score_sns)


gmm_model <- Mclust(df_truthset_gmm$function_score_sns, G=2, modelNames=c("V"))
summary(gmm_model)

df_vars_for_pvalues <- df_HMEC_HAP1_all_spread2 %>% 
  filter(expt == "ut")

pathogenic_index <- which.min(gmm_model$parameters$mean)  # Lower mean = pathogenic
df_truthset_gmm$prob_pathogenic <- gmm_model$z[, pathogenic_index]
df_vars_for_pvalues$prob_pathogenic <- predict(gmm_model, newdata = df_vars_for_pvalues$function_score_sns)$z[, pathogenic_index]

df_vars_for_pvalues$prob_pathogenic <- predict(gmm_model, newdata = df_vars_for_pvalues$function_score_sns)$z[, which.min(gmm_model$parameters$mean)]

prior_prob <- 0.1
Threshold <- c(-Inf, 0.053, 0.23, 0.48, 2.1, 4.3, 18.7, Inf)
Labels <- c("BS3", "BS3_mod", "BS3_supp", "indeterminate", "PS3_supp", "PS3_mod", "PS3")

df_vars_for_pvalues2 <- df_vars_for_pvalues %>% 
  select(cHGVS, pHGVS, pos, hg38, Ref, alt, sge_region, conseq, HAP1_function_score_mean, HAP1_func_class, final_func_class, Clinvar_interpretation_Sep24, Clinvar_Sep24_histo, Clinvar_link, Review.status, function_score_sns, HMEC_class, prob_pathogenic, SpliceAI_max.D14) %>% 
  mutate(OddsPath = (prob_pathogenic*(1-prior_prob))/((1-prob_pathogenic)*prior_prob)) %>% 
  mutate(evidence_code = cut(OddsPath, breaks = Threshold, labels = Labels))

#adding points for evidence codes
evidence_code_points <- read.csv("databases/250207_evidence_code_points.csv")

new_variants_with_points <- left_join(x = df_vars_for_pvalues2,
                                                y = evidence_code_points,
                                                by = "evidence_code")

##### Adjusting evidence codes #####
  new_variants_with_points_adjusted <- new_variants_with_points %>% 
    dplyr::rename("HMEC_ev_code_original" = "evidence_code",
                  "HMEC_points_original" = "points") %>% 
    mutate(HMEC_ev_code_adj = case_when(HMEC_class == "Neutral" & (HMEC_ev_code_original == "PS3" | HMEC_ev_code_original == "PS3_supp" | HMEC_ev_code_original == "PS3_mod") ~ "indeterminate",
                                        TRUE ~ HMEC_ev_code_original),
           HMEC_points_adj = case_when(HMEC_class == "Neutral" & (HMEC_ev_code_original == "PS3" | HMEC_ev_code_original == "PS3_supp" | HMEC_ev_code_original == "PS3_mod") ~ 0,
                                       TRUE ~ HMEC_points_original)
    )
  
#Combining with HAP1 points
  HAP1_points_final <- read.csv("~/Dropbox (The Francis Crick)/BRCA1_HAP1_HMEC_manuscript/Data/250722_HAP1_ev_codes_and_points_mclust_corrected_with_final_rules_applied.csv")
  
  
  df_HAP1_HMEC_points_adj_final <- left_join(x = new_variants_with_points_adjusted,
                                             y = HAP1_points_final,
                                             by = "cHGVS",
                                             suffix = c(".HMEC", ".HAP1")) %>% 
    filter(!is.na(HAP1_ev_code_adj)) %>% 
    mutate(
      total_points_original = HAP1_points_original + HMEC_points_original,
      total_points_both_adj = HAP1_points_adj + HMEC_points_adj
    )
  

  
  ##### Sankey for HAP+HMEC adjustment (final adjustments applied) #####
  df_sankey_HAP1_to_HMEC_points_both_adj <- df_HAP1_HMEC_points_adj_final %>% 
    select(HAP1_ev_code_adj, HMEC_ev_code_adj) %>% 
    mutate(source = case_when(HAP1_ev_code_adj == "BS3" ~ 0,
                              HAP1_ev_code_adj == "BS3_mod" ~ 1,
                              HAP1_ev_code_adj == "BS3_supp" ~ 2,
                              HAP1_ev_code_adj == "indeterminate" ~ 3,
                              HAP1_ev_code_adj == "PS3_supp" ~ 4,
                              HAP1_ev_code_adj == "PS3_mod" ~ 5,
                              HAP1_ev_code_adj == "PS3" ~ 6,
    )) %>% 
    mutate(target = case_when(HMEC_ev_code_adj == "BS3" ~ 7,
                              HMEC_ev_code_adj == "BS3_mod" ~ 8,
                              HMEC_ev_code_adj == "BS3_supp" ~ 9,
                              HMEC_ev_code_adj == "indeterminate" ~ 10,
                              HMEC_ev_code_adj == "PS3_supp" ~ 11,
                              HMEC_ev_code_adj == "PS3_mod" ~ 12,
                              HMEC_ev_code_adj == "PS3" ~ 13,
    )) %>% 
    group_by(source, target) %>% 
    summarise(value = n())
  
  
  df_sankey_HMEC_to_points_both_adj <- df_HAP1_HMEC_points_adj_final %>% 
    select(total_points_both_adj, HMEC_ev_code_adj) %>% 
    mutate(source = case_when(HMEC_ev_code_adj == "BS3" ~ 7,
                              HMEC_ev_code_adj == "BS3_mod" ~ 8,
                              HMEC_ev_code_adj == "BS3_supp" ~ 9,
                              HMEC_ev_code_adj == "indeterminate" ~ 10,
                              HMEC_ev_code_adj == "PS3_supp" ~ 11,
                              HMEC_ev_code_adj == "PS3_mod" ~ 12,
                              HMEC_ev_code_adj == "PS3" ~ 13,
    )) %>% 
    mutate(target = case_when(total_points_both_adj == -8 ~ 14,
                              total_points_both_adj == -6 ~ 15,
                              total_points_both_adj == -5 ~ 16,
                              total_points_both_adj == -4 ~ 17,
                              total_points_both_adj == -3 ~ 18,
                              total_points_both_adj == -2 ~ 19,
                              total_points_both_adj == -1 ~ 20,
                              total_points_both_adj == 0 ~ 21,
                              total_points_both_adj == 1 ~ 22,
                              total_points_both_adj == 2 ~ 23,
                              total_points_both_adj == 3 ~ 24,
                              total_points_both_adj == 4 ~ 25,
                              total_points_both_adj == 5 ~ 26,
                              total_points_both_adj == 6 ~ 27,
                              total_points_both_adj == 8 ~ 28
                              
    )) %>% 
    group_by(source, target) %>% 
    summarise(value = n())
  
   df_sankey_3_column_points_both_adj <- rbind(df_sankey_HAP1_to_HMEC_points_both_adj, df_sankey_HMEC_to_points_both_adj)
  
  color_scale_points <- 'd3.scaleOrdinal()
                .domain(["BS3", "BS3_mod", "BS3_supp", "indeterminate", "PS3_supp", "PS3_mod", "PS3", "-8", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4", "5", "6", "8"])
                .range(["#36648B", "#4F94CD", "#63B8FF", "#EEE9E9", "#FFC1C1", "#EE6363", "#CD0000", "#122a47", "#2b446b", "#3e5e84", "#36648B", "#2f77a8", "#4F94CD", "#63B8FF", "#EEE9E9", "#FFC1C1", "#EE6363", "#cc4a4a", "#CD0000", "#9b1e1e", "#770505", "#770505"])'
  
  
 
  plot_sankey_3_column_HAP1_HMEC_points_both_adj <- sankeyNetwork(Links = df_sankey_3_column_points_both_adj, Nodes = nodes_points3, Source = "source",
                                                                  Target = "target", Value = "value", NodeID = "name",
                                                                  units = "Overlap",
                                                                  fontSize = 0,
                                                                  nodeWidth = 60,
                                                                  fontFamily = "Helvetica",
                                                                  colourScale = color_scale_points)  
  
  plot_sankey_final_render <- onRender(plot_sankey_3_column_HAP1_HMEC_points_both_adj, '
  function(el) {
    var colors = {};

    // Save the node colours by name
    d3.select(el).selectAll(".node rect").each(function(d) {
      colors[d.name] = d3.select(this).style("fill");
    });

    // Style the links to match the source node colour
    d3.select(el).selectAll(".link")
      .style("stroke", function(d) {
        return colors[d.source.name];
      })
      .style("stroke-opacity", 0.6)
      .style("fill", "none")
      .style("stroke-width", function(d) { return Math.max(1, d.dy); })
      .style("pointer-events", "none"); // disables hover changes entirely
  }
')
  
  #save
  saveWidget(plot_sankey_final_render, "250724_sankey_plot_final.html", selfcontained = TRUE)
  
  #modify the html - put columns in correct order etc
  #save with new names
  
  #FINAL PLOT
  webshot("250724_sankeyNetwork_edited.html", "250724_sankeyNetwork_for_illustrator.pdf", vwidth = 1000, vheight = 800)
  
  
 
#### EXPORT FOR SUPPLEMENTARY TABLE 3 ####
##### Adding data from external sources #####

df_export_ut <- df_HMEC_HAP1_all_spread2 %>% 
  filter(expt == "ut") %>% 
  dplyr::rename("q_value_ut" = "fdr",
                "final_function_score_ut" = "function_score_sns",
                "threshold_ut" = "fs_threshold",
                "r1_function_score_ut" = "combined_FS_r1",
                "r2_function_score_ut" = "combined_FS_r2",
                "HMEC_class_ut" = "HMEC_class"
  ) %>% 
  select(-expt)


df_export_olap <- df_HMEC_HAP1_all_spread2 %>% 
  filter(expt == "olaparib") %>% 
  select(cHGVS, function_score_sns, fdr, fs_threshold, combined_FS_r1, combined_FS_r2, HMEC_class) %>% 
  dplyr::rename("q_value_olaparib" = "fdr",
                "final_function_score_olaparib" = "function_score_sns",
                "threshold_olaparib" = "fs_threshold",
                "r1_function_score_olaparib" = "combined_FS_r1",
                "r2_function_score_olaparib" = "combined_FS_r2",
                "HMEC_class_olaparib" = "HMEC_class"
                )

df_export_all <- left_join(x = df_export_ut,
          y = df_export_olap,
          by = "cHGVS") %>% 
  select(-sge_region.x, -Clinvar_link, -canvar_link, -Dst2Splice, -CADD.raw, -rev_comp, -Ref_RC, -within_2bp_of_pam_edit, -VariationID, -LoF) %>% 
  dplyr::rename("sge_region.HAP1" = "sge_region.y")

#Add AOU data
df_AOU_export <- df_AOU2 %>% 
  select(cHGVS, alleleCount, alleleFrequency, present_in_AOU) %>% 
  dplyr::rename("alleleCount_AOU" = "alleleCount",
                "alleleFrequency_AOU" = "alleleFrequency"
                )

df_export_all2 <- left_join(x = df_export_all,
          y = df_AOU_export,
          by = "cHGVS")

#Add gnomad
gnomadv4_export <- gnomadv4 %>% 
  select(Position, Reference, Alternate, Allele.Count, Allele.Frequency) %>% 
  dplyr::rename("alleleCount_gnomad" = "Allele.Count",
                "alleleFrequency_gnomad" = "Allele.Frequency")
  

df_export_all3 <- left_join(df_export_all2,
                            gnomadv4_export,
                            by = c("hg38" = "Position",
                                   "Ref" = "Reference",
                                   "alt" = "Alternate")) %>% 
  mutate(present_in_gnomad = case_when(alleleCount_gnomad > 0 ~ TRUE,
                                       is.na(alleleCount_gnomad) ~ FALSE))

#add AM scores
alphamissense_export <- alphamissense %>% 
  select(-position, -a.a.2) %>% 
  dplyr::rename("AlphaMissense_score" = `pathogenicity score`,
                "AlphaMissense_class" = `pathogenicity class`
                )

df_export_all4 <- left_join(x = df_export_all3,
           y = alphamissense_export,
           by = c("pHGVS"))

#add EVE scores
eve_export <- eve_scores %>% 
  select(position, wt_aa, mt_aa, EVE_scores_ASM) %>% 
  dplyr::rename("EVE_score" = "EVE_scores_ASM")

df_export_all5 <- left_join(x = df_export_all4,
                       y = eve_export,
                       by = c("protPos" = "position",
                              "oAA" = "wt_aa",
                              "nAA" = "mt_aa"))

assay_data <- df_sge_data_and_other_assays_T131_2 %>% 
  select(cHGVS, count_ones, count_zeros, total_assays, pct_path, pct_benign) %>% 
  dplyr::rename("count_abnormal" = "count_ones",
                "count_normal" = "count_zeros",
                "pct_abnormal" = "pct_path",
                "pct_normal" = "pct_benign")

df_export_all6 <- left_join(x = df_export_all5,
                            y = assay_data,
                            by = "cHGVS")

#Adding Oddspath  
#HMEC first to make sure I keep HMEC x3 info:
HMEC_truthset <- as.character(df_truthset$cHGVS)

df_HMEC_points_export <- new_variants_with_points_adjusted %>% 
  select(cHGVS, prob_pathogenic, OddsPath, HMEC_ev_code_original, HMEC_points_original, HMEC_ev_code_adj, HMEC_points_adj) %>% 
  dplyr::rename("HMEC_OddsPath" = "OddsPath",
                "HMEC_evidence_code_original" = "HMEC_ev_code_original",
                "HMEC_evidence_code_adjusted" = "HMEC_ev_code_adj",
                "HMEC_prob_pathogenic" = "prob_pathogenic",
                "HMEC_points_adjusted" = "HMEC_points_adj") %>% 
  mutate(in_HMEC_truthset = case_when(cHGVS %in% HMEC_truthset ~ TRUE,
                                      TRUE ~ FALSE))

df_export_all7 <- left_join(x = df_export_all6,
                            y = df_HMEC_points_export,
                            by = "cHGVS")


#Now add the HAP1 OddsPath info
df_HAP1_HMEC_points_export <- df_HAP1_HMEC_points_adj_final %>% 
  select(cHGVS, HAP1_prob_pathogenic, OddsPath_corrected, HAP1_ev_code_adj, HAP1_points_adj, total_points_both_adj) %>% 
  dplyr::rename("HAP1_OddsPath" = "OddsPath_corrected",
                "total_points" = "total_points_both_adj",
                "HAP1_points_adjusted" = "HAP1_points_adj",
                "HAP1_evidence_code_adjusted" = "HAP1_ev_code_adj"
  )

df_export_all8 <- left_join(x = df_export_all7,
                            y = df_HAP1_HMEC_points_export,
                            by = "cHGVS")

#Adding new SpliceAI scores
df_new_SAI_scores <- read.csv("~/Downloads/SpliceAI_new_scores.csv")

df_export_all9 <- left_join(x = df_export_all8,
                            y = df_new_SAI_scores,
                            by = "cHGVS")

df_export_all10 <- left_join(x = df_export_all9,
                             y = foldX_all,
                             by = "pHGVS") %>% 
  select(-domain)


##### getting read counts from previous df ##### 
  
  all_dfs_ut <-
  all_dfs_final %>% 
    filter(expt == "ut") %>% 
    filter(sge_region != "x17" & sge_region != "x17q") %>% 
    select(hg38, Ref, alt, sge_region, conseq, cHGVS, pHGVS, tHDR_neg.D14, tHDR_lib.D14, tHDR_pre.D14, tHDR_post.D14, tHDR_post.D21, rL42_tHDR_pre.D14, rL42_tHDR_post.D14, rL42_tHDR_post.D21, combined_FS_r1, combined_FS_r2, mean_combined_fs, function_score_sns) %>% 
    dplyr::rename("neg_ut.r1r2" = "tHDR_neg.D14",
                  "lib_ut.r1r2" = "tHDR_lib.D14",
                  "r1_D7_ut" = "tHDR_pre.D14",
                  "r2_D7_ut" = "rL42_tHDR_pre.D14",
                  "r1_D14_ut" = "tHDR_post.D14",
                  "r2_D14_ut" = "rL42_tHDR_post.D14",
                  "r1_D21_ut" = "tHDR_post.D21",
                  "r2_D21_ut" = "rL42_tHDR_post.D21",
                  "r1_function_score_ut" = "combined_FS_r1",
                  "r2_function_score_ut" = "combined_FS_r2",
                  "mean_function_score_ut" = "mean_combined_fs",
                  "final_function_score_ut" = "function_score_sns",
                  )
    
    all_dfs_olap <- all_dfs_final %>% 
    filter(expt == "olaparib") %>%
    filter(sge_region != "x17" & sge_region != "x17q") %>% 
    select(cHGVS, sge_region, tHDR_pre.D14, tHDR_neg.D14, tHDR_lib.D14, tHDR_post.D14, tHDR_post.D21, rL42_tHDR_pre.D14, rL42_tHDR_post.D14, rL42_tHDR_post.D21, combined_FS_r1, combined_FS_r2, mean_combined_fs, function_score_sns) %>% 
    dplyr::rename("neg_olaparib.r1r2" = "tHDR_neg.D14",
                  "lib_olaparib.r1r2" = "tHDR_lib.D14",
                  "r1_D7_olaparib" = "tHDR_pre.D14",
                  "r2_D7_olaparib" = "rL42_tHDR_pre.D14",
                  "r1_D14_olaparib" = "tHDR_post.D14",
                  "r2_D14_olaparib" = "rL42_tHDR_post.D14",
                  "r1_D21_olaparib" = "tHDR_post.D21",
                  "r2_D21_olaparib" = "rL42_tHDR_post.D21",
                  "r1_function_score_olaparib" = "combined_FS_r1",
                  "r2_function_score_olaparib" = "combined_FS_r2",
                  "mean_function_score_olaparib" = "mean_combined_fs",
                  "final_function_score_olaparib" = "function_score_sns",
    )
    
    #Join ut and olap back together - this is all data with read count information
    all_dfs_ut_olap_joined <- left_join(x = all_dfs_ut, 
               y = all_dfs_olap,
               by = c("cHGVS", "sge_region"))# %>% 
    

df_export_all10_no_x17 <- df_export_all10 %>% 
  filter(sge_region != "x17m")

df_export_all10_x17 <- df_export_all10 %>% 
  filter(sge_region == "x17m")


##### exon 17 scores #####
#getting raw reads from earlier dataframe
#x17 - need to get raw reads from this (4 reps +/-olap)   
x17_df_ut <- all_dfs_final %>% 
  filter(expt == "ut") %>%
  filter(sge_region == "x17") %>% 
  select(cHGVS, sge_region, tHDR_pre.D14, tHDR_neg.D14, tHDR_lib.D14, tHDR_post.D14, tHDR_post.D21, rL42_tHDR_pre.D14, rL42_tHDR_post.D14, rL42_tHDR_post.D21, combined_FS_r1, combined_FS_r2, function_score_sns) %>% 
  dplyr::rename("neg_ut" = "tHDR_neg.D14",
                "lib_ut" = "tHDR_lib.D14",
                "r1_D7_ut" = "tHDR_pre.D14",
                "r2_D7_ut" = "rL42_tHDR_pre.D14",
                "r1_D14_ut" = "tHDR_post.D14",
                "r2_D14_ut" = "rL42_tHDR_post.D14",
                "r1_D21_ut" = "tHDR_post.D21",
                "r2_D21_ut" = "rL42_tHDR_post.D21",
                "r1_function_score_ut" = "combined_FS_r1",
                "r2_function_score_ut" = "combined_FS_r2",
                "r1r2function_score_ut" = "function_score_sns",
  )

x17q_df_ut <- all_dfs_final %>% 
  filter(expt == "ut") %>%
  filter(sge_region == "x17q") %>% 
  select(cHGVS, sge_region, tHDR_pre.D14, tHDR_neg.D14, tHDR_lib.D14, tHDR_post.D14, tHDR_post.D21, rL42_tHDR_pre.D14, rL42_tHDR_post.D14, rL42_tHDR_post.D21, combined_FS_r1, combined_FS_r2, function_score_sns) %>% 
  dplyr::rename("neg_ut" = "tHDR_neg.D14",
                "lib_ut" = "tHDR_lib.D14",
                "r3_D7_ut" = "tHDR_pre.D14",
                "r4_D7_ut" = "rL42_tHDR_pre.D14",
                "r3_D14_ut" = "tHDR_post.D14",
                "r4_D14_ut" = "rL42_tHDR_post.D14",
                "r3_D21_ut" = "tHDR_post.D21",
                "r4_D21_ut" = "rL42_tHDR_post.D21",
                "r3_function_score_ut" = "combined_FS_r1",
                "r4_function_score_ut" = "combined_FS_r2",
                "r3r4function_score_ut" = "function_score_sns",
  )


x17_df_olap <- all_dfs_final %>% 
  filter(expt == "olaparib") %>%
  filter(sge_region == "x17") %>% 
  select(cHGVS, sge_region, tHDR_pre.D14, tHDR_neg.D14, tHDR_lib.D14, tHDR_post.D14, tHDR_post.D21, rL42_tHDR_pre.D14, rL42_tHDR_post.D14, rL42_tHDR_post.D21, combined_FS_r1, combined_FS_r2, function_score_sns) %>% 
  dplyr::rename("neg_olaparib" = "tHDR_neg.D14",
                "lib_olaparib" = "tHDR_lib.D14",
                "r1_D7_olaparib" = "tHDR_pre.D14",
                "r2_D7_olaparib" = "rL42_tHDR_pre.D14",
                "r1_D14_olaparib" = "tHDR_post.D14",
                "r2_D14_olaparib" = "rL42_tHDR_post.D14",
                "r1_D21_olaparib" = "tHDR_post.D21",
                "r2_D21_olaparib" = "rL42_tHDR_post.D21",
                "r1_function_score_olaparib" = "combined_FS_r1",
                "r2_function_score_olaparib" = "combined_FS_r2",
                "r1r2function_score_olaparib" = "function_score_sns",
  )

x17q_df_olap <- all_dfs_final %>% 
  filter(expt == "olaparib") %>%
  filter(sge_region == "x17q") %>% 
  select(cHGVS, sge_region, tHDR_pre.D14, tHDR_neg.D14, tHDR_lib.D14, tHDR_post.D14, tHDR_post.D21, rL42_tHDR_pre.D14, rL42_tHDR_post.D14, rL42_tHDR_post.D21, combined_FS_r1, combined_FS_r2, function_score_sns) %>% 
  dplyr::rename("neg_olaparib" = "tHDR_neg.D14",
                "lib_olaparib" = "tHDR_lib.D14",
                "r3_D7_olaparib" = "tHDR_pre.D14",
                "r4_D7_olaparib" = "rL42_tHDR_pre.D14",
                "r3_D14_olaparib" = "tHDR_post.D14",
                "r4_D14_olaparib" = "rL42_tHDR_post.D14",
                "r3_D21_olaparib" = "tHDR_post.D21",
                "r4_D21_olaparib" = "rL42_tHDR_post.D21",
                "r3_function_score_olaparib" = "combined_FS_r1",
                "r4_function_score_olaparib" = "combined_FS_r2",
                "r3r4function_score_olaparib" = "function_score_sns",
  )


x17m_df_ut <- full_join(x = x17_df_ut,
                        y = x17q_df_ut,
                        by = "cHGVS",
                        suffix = c(".r1r2", ".r3r4")) %>%
  select(-sge_region.r1r2, -sge_region.r3r4) %>% 
  mutate(sge_region = "x17m")


x17m_df_olap <- full_join(x = x17_df_olap,
                          y = x17q_df_olap,
                          by = "cHGVS",
                          suffix = c(".r1r2", ".r3r4")) %>%
  select(-sge_region.r1r2, -sge_region.r3r4) %>% 
  mutate(sge_region = "x17m")

#This data frame contains read counts for all x17 reps
df_x17m_with_read_counts <- full_join(x = x17m_df_ut,
                                      y = x17m_df_olap,
                                      by = c("cHGVS", "sge_region"))

#keep only columns I need - cHGVS, replicate function scores, final score.
df_x17_final_scores <- df_x17x17q_output_filtered %>% 
  select(cHGVS, expt, x17_combined_FS_r1, x17_combined_FS_r2, x17q_combined_FS_r1, x17q_combined_FS_r2, mean_function_score_sns, fdr)

df_x17_final_scores_ut <- df_x17_final_scores %>% 
  filter(expt == "ut") %>% 
  select(-expt)
df_x17_final_scores_olap <- df_x17_final_scores %>% 
  filter(expt == "olaparib") %>% 
  select(-expt)

df_x17_final_scores2 <- left_join(x = df_x17_final_scores_ut,
          y = df_x17_final_scores_olap,
          by = "cHGVS",
          suffix = c(".ut", ".olaparib"))

#joining df with reads and df with individual function scores for each rep 
#must make sure the columns match up here

df_x17_final_scores3 <- df_x17_final_scores2 %>% 
  dplyr::rename("r1_function_score_olaparib" = "x17_combined_FS_r1.olaparib",
                "r2_function_score_olaparib" = "x17_combined_FS_r2.olaparib",
                "r3_function_score_olaparib" = "x17q_combined_FS_r1.olaparib",
                "r4_function_score_olaparib" = "x17q_combined_FS_r2.olaparib",
                "r1_function_score_ut" = "x17_combined_FS_r1.ut",
                "r2_function_score_ut" = "x17_combined_FS_r2.ut",
                "r3_function_score_ut" = "x17q_combined_FS_r1.ut",
                "r4_function_score_ut" = "x17q_combined_FS_r2.ut"
                )


df_x17_scores_and_counts <- 
  left_join(x = df_x17_final_scores3,
          y = df_x17m_with_read_counts,
          by = c("cHGVS", "r1_function_score_olaparib", "r2_function_score_olaparib", "r3_function_score_olaparib", "r4_function_score_olaparib", "r1_function_score_ut", "r2_function_score_ut", "r3_function_score_ut", "r4_function_score_ut"))  %>% 
  mutate(mean_function_score_olaparib = (r1_function_score_olaparib + r2_function_score_olaparib + r3_function_score_olaparib + r4_function_score_olaparib)/4,
         mean_function_score_ut = (r1_function_score_ut + r2_function_score_ut + r3_function_score_ut + r4_function_score_ut)/4) %>% 
    dplyr::rename("final_function_score_ut" = "mean_function_score_sns.ut",
                   "final_function_score_olaparib" = "mean_function_score_sns.olaparib",
                  "q_value_ut" = "fdr.ut",
                  "q_value_olaparib" = "fdr.olaparib"
                   ) %>% 
    select(-r1r2function_score_ut, -r1r2function_score_olaparib, -r3r4function_score_ut, -r3r4function_score_olaparib) #don't need these intermediate scores
  

  ##### Joining x17 with rest of data ##### 
  all_dfs_ut_olap_joined2 <- all_dfs_ut_olap_joined %>% 
    select(-Ref, -alt, -hg38, -conseq, -pHGVS)

df_x17_scores_and_counts2 <- df_x17_scores_and_counts %>% 
  select(-q_value_ut, -q_value_olaparib) #get rid of these columns as not present in all_dfs_ut_olap_joined2
 
#Bind together scores and counts for x17 and rest of data
df_all_scores_and_counts <- bind_rows(all_dfs_ut_olap_joined2, df_x17_scores_and_counts2)


##### Joining counts and scores with df_export_all10 #####  

  df_export_all11 <- left_join(x = df_export_all10,
          y = df_all_scores_and_counts,
          by = c("cHGVS", "sge_region", "final_function_score_ut", "final_function_score_olaparib")) %>%   #checked that all these columns match
          select(-r1_function_score_ut.x, -r1_function_score_olaparib.x, -r2_function_score_ut.x, -r2_function_score_olaparib.x) #get rid of r1 and r2 values that are incorrect for x17
   
##### Adding new SpliceAI scores #####

   df_spliceAI_scores_new <- read.csv(file = "~/OneDrive - The Francis Crick Institute/R/250127_new_spliceAI_scores/250703_SpliceAI_scores_500bp_all_vars.csv")   
   
   df_spliceAI_scores_new2 <- df_spliceAI_scores_new %>% 
     select(-Ref, -alt)
   
   df_export_all12 <- left_join(x = df_export_all11,
                                y = df_spliceAI_scores_new2,
                                by = c("cHGVS", "hg38"))
   
   df_export_all13 <- df_export_all12 %>% 
     select(-SpliceAI_max.D14, -AG.x, -AL.x, -DG.x, -DL.x, -SpliceAI_transcript) %>% 
     dplyr::rename("AG" = "AG.y",
                   "DG" = "DG.y",
                   "AL" = "AL.y",
                   "DL" = "DL.y",
                   "SpliceAI_transcript" = "Transcript"
                   ) %>% 
     mutate(SpliceAI_max = case_when(SpliceAI_transcript != "ENST00000357654.9" ~ NA_real_,
                                     TRUE ~ SpliceAI_max)) %>% 
     select(-SpliceAI_transcript)
  
   ##### Cleaning up final columns and export #####
   df_export_all14 <- df_export_all13 %>% 
     mutate(cHGVS = case_when(hg38 > 43125364 ~ NA_character_,
                              pHGVS == "p.R1699Q" ~ "c.5096_5097delinsAA",
                              TRUE ~ cHGVS),
            Exon = case_when(hg38 > 43125364 ~ NA_character_,
                             TRUE ~ Exon),
            conseq = case_when(hg38 > 43125364 ~ "Promoter",
                               TRUE ~ conseq),
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
            hg38 = case_when(pHGVS == "p.R1699Q" ~ 43063929,
                             TRUE ~ hg38),
            sge_region = case_when(sge_region == "u1" ~ "promoter",
                                   sge_region == "x1" ~ "5' UTR",
                                   sge_region == "x3" ~ "exon 3",
                                   sge_region == "x5" ~ "exon 5",
                                   sge_region == "x6a" ~ "exon 6",
                                   sge_region == "x17q" ~ "exon 17",
                                   sge_region.HAP1 == "x17" ~ "exon 17",
                                   sge_region.HAP1 == "x17m" ~ "exon 17",
                                   sge_region == "x12a" ~ "exon 12 (5')",
                                   sge_region == "x12b" ~ "exon 12 (3')",
                                   sge_region == "x10a" ~ "exon 10 (5')",
                                   sge_region == "x10d" ~ "exon 10 (3')",
                                   sge_region == "x10h2" ~ "exon 10 (mid 1)",
                                   sge_region == "x10h" ~ "exon 10 (mid 2)",
                                   sge_region == "x11" ~ "exon 11",
                                   sge_region == "i11" ~ "intron 11"),
            sge_region.HAP1 = case_when(sge_region.HAP1 == "u1" ~ "promoter",
                                        sge_region.HAP1 == "x1" ~ "5' UTR",
                                        sge_region.HAP1 == "x3" ~ "exon 3",
                                        sge_region.HAP1 == "x5" ~ "exon 5",
                                        sge_region.HAP1 == "x6a" ~ "exon 6",
                                        sge_region.HAP1 == "x17q" ~ "exon 17",
                                        sge_region.HAP1 == "x17" ~ "exon 17",
                                        sge_region.HAP1 == "x17m" ~ "exon 17",
                                        sge_region.HAP1 == "x12a" ~ "exon 12 (5')",
                                        sge_region.HAP1 == "x12b" ~ "exon 12 (3')",
                                        sge_region.HAP1 == "x10a" ~ "exon 10 (5')",
                                        sge_region.HAP1 == "x10d" ~ "exon 10 (3')",
                                        sge_region.HAP1 == "x10h2" ~ "exon 10 (mid 1)",
                                        sge_region.HAP1 == "x10h" ~ "exon 10 (mid 2)",
                                        sge_region.HAP1 == "x11" ~ "exon 11",
                                        sge_region.HAP1 == "i11" ~ "intron 11"),
            HAP1_evidence_code_adjusted = case_when(HAP1_evidence_code_adjusted == "PS3" ~ "P_strong",
                                                    HAP1_evidence_code_adjusted == "PS3_mod" ~ "P_moderate",
                                                    HAP1_evidence_code_adjusted == "PS3_supp" ~ "P_supporting",
                                                    HAP1_evidence_code_adjusted == "BS3" ~ "B_strong",
                                                    HAP1_evidence_code_adjusted == "BS3_mod" ~ "B_moderate",
                                                    HAP1_evidence_code_adjusted == "BS3_supp" ~ "B_supporting",
                                                    TRUE ~ HAP1_evidence_code_adjusted),
            HMEC_evidence_code_original = case_when(HMEC_evidence_code_original == "PS3" ~ "P_strong",
                                                    HMEC_evidence_code_original == "PS3_mod" ~ "P_moderate",
                                                    HMEC_evidence_code_original == "PS3_supp" ~ "P_supporting",
                                                    HMEC_evidence_code_original == "BS3" ~ "B_strong",
                                                    HMEC_evidence_code_original == "BS3_mod" ~ "B_moderate",
                                                    HMEC_evidence_code_original == "BS3_supp" ~ "B_supporting",
                                                    TRUE ~ HMEC_evidence_code_original),
            HMEC_evidence_code_adjusted = case_when(HMEC_evidence_code_adjusted == "PS3" ~ "P_strong",
                                                    HMEC_evidence_code_adjusted == "PS3_mod" ~ "P_moderate",
                                                    HMEC_evidence_code_adjusted == "PS3_supp" ~ "P_supporting",
                                                    HMEC_evidence_code_adjusted == "BS3" ~ "B_strong",
                                                    HMEC_evidence_code_adjusted == "BS3_mod" ~ "B_moderate",
                                                    HMEC_evidence_code_adjusted == "BS3_supp" ~ "B_supporting",
                                                    TRUE ~ HMEC_evidence_code_adjusted),
            HMEC_class_ut = case_when(HMEC_class_ut == "LoF" ~ "Depleted",
                                      TRUE ~ HMEC_class_ut),
            HMEC_class_olaparib = case_when(HMEC_class_olaparib == "LoF" ~ "Depleted",
                                            TRUE ~ HMEC_class_olaparib),
            final_func_class = case_when(final_func_class == "HAP1 neutral, HMEC LoF" ~ "HAP1 neutral, HMEC depleted",
                                         final_func_class == "HAP1 int, HMEC LoF" ~ "HAP1 int, HMEC depleted",
                                         final_func_class == "Neutral" ~ "Concordant neutral",
                                         final_func_class == "LoF" ~ "Concordant LoF",
                                         TRUE ~ final_func_class
                                         )
            ) %>%
     dplyr::rename("ClinVar" = "Clinvar_interpretation_Sep24",
                   "ClinVar_simple" = "Clinvar_Sep24_histo",
                   "Consequence" = "conseq",
                   "SGE_region" = "sge_region",
                   "Alt" = "alt",
                   "ClinVar_review_status" = "Review.status",
                   "HAP1_SGE_region" = "sge_region.HAP1",
                   "r1_function_score_olaparib" = "r1_function_score_olaparib.y",
                   "r2_function_score_olaparib" = "r2_function_score_olaparib.y",
                   "r1_function_score_ut" = "r1_function_score_ut.y",
                   "r2_function_score_ut" = "r2_function_score_ut.y"
     ) %>%
     mutate(HAP1_data_source = case_when(HAP1_data_source == "old" ~ "Findlay et al. 2018",
                                         TRUE ~ HAP1_data_source)) %>% 
     select(-pos, -Intron, -Exon, -Last_reviewed, -fs_sig, -fs_threshold_upper, -fs_threshold_lower2, -HAP1_function_score_r1, -HAP1_function_score_r2, -HAP1_mean_RNA_score, -HAP1_FDR, -threshold_ut, -threshold_olaparib, -mean_function_score_ut, -mean_function_score_olaparib, -HAP1_prob_pathogenic, -HAP1_OddsPath) %>% 
     #
     select(hg38, Alt, Ref, cHGVS, pHGVS, protPos, oAA, nAA, Consequence, CDSpos, SGE_region, HAP1_SGE_region, HAP1_function_score_mean, HAP1_func_class, HAP1_data_source, final_function_score_ut, q_value_ut, HMEC_class_ut, final_func_class, r1_function_score_ut, r2_function_score_ut, r3_function_score_ut, r4_function_score_ut, r1_function_score_olaparib, r2_function_score_olaparib, r3_function_score_olaparib, r4_function_score_olaparib, final_function_score_olaparib, q_value_olaparib, HMEC_class_olaparib, ClinVar, ClinVar_simple, ClinVar_review_status, CADD.phred, AlphaMissense_score, AlphaMissense_class, EVE_score, FoldX_score, AG, AL, DG, DL, SpliceAI_max, present_in_gnomad, alleleCount_gnomad, alleleFrequency_gnomad, present_in_AOU, alleleCount_AOU, alleleFrequency_AOU, total_assays, count_normal, count_abnormal, pct_normal, pct_abnormal, HAP1_evidence_code_adjusted, HAP1_points_adjusted, in_HMEC_truthset, HMEC_prob_pathogenic, HMEC_OddsPath, HMEC_evidence_code_original, HMEC_points_original, HMEC_evidence_code_adjusted, HMEC_points_adjusted, total_points, lib_ut.r1r2, neg_ut.r1r2, lib_olaparib.r1r2, neg_olaparib.r1r2, r1_D7_ut, r1_D14_ut, r1_D21_ut, r1_D7_olaparib, r1_D14_olaparib, r1_D21_olaparib, r2_D7_ut, r2_D14_ut, r2_D21_ut, r2_D7_olaparib, r2_D14_olaparib, r2_D21_olaparib, lib_ut.r3r4, neg_ut.r3r4, lib_olaparib.r3r4, neg_olaparib.r3r4, r3_D7_ut, r3_D14_ut, r3_D21_ut, r3_D7_olaparib, r3_D14_olaparib, r3_D21_olaparib, r4_D7_ut, r4_D14_ut, r4_D21_ut, r4_D7_olaparib, r4_D14_olaparib, r4_D21_olaparib) %>% 
     arrange(desc(hg38))

   #write.csv(x = df_export_all14, file = paste0("~/Dropbox (The Francis Crick)/BRCA1_HAP1_HMEC_manuscript/SUPPLEMENTARY_TABLES/", date, "_HAP1-HMEC_combined_data.csv"), row.names = FALSE) 
   
##### RNF168 dataframe #####
df_RNF168_SGE <- read.csv(file = "output_dataframes/250317_HMEC_RNF168_R-C2_x12a_SNS_COMBINED_FSs_incl_olap.csv")
   
df_RNF168_export <- df_RNF168_SGE %>% 
  select(sge_region, expt, hg38, Ref, alt, cHGVS, pHGVS, conseq, tHDR_pre.D14, rL42_tHDR_pre.D14, tHDR_post.D14, rL42_tHDR_post.D14, tHDR_post.D21, rL42_tHDR_post.D21, tHDR_lib.D14, tHDR_neg.D14, combined_FS_r1, combined_FS_r2, mean_combined_fs, function_score_sns, fdr, fs_threshold, HMEC_class) %>% 
  mutate(sge_region = case_when(sge_region == "Rx12a" ~ "HMEC RNF168-KO",
                                sge_region == "x12a" ~ "HMEC RNF168-WT"),
         expt = case_when(expt == "ut" ~ "untreated",
                               TRUE ~ expt)) %>% 
  dplyr::rename("r1_day_7" = "tHDR_pre.D14",
                "r2_day_7" = "rL42_tHDR_pre.D14",
                "r1_day_14" = "tHDR_post.D14",
                "r2_day_14" = "rL42_tHDR_post.D14",
                "r1_day_21" = "tHDR_post.D21",
                "r2_day_21" = "rL42_tHDR_post.D21",
                "library" = "tHDR_lib.D14",
                "negative" = "tHDR_neg.D14",
                "r1_function_score" = "combined_FS_r1",
                "r2_function_score" = "combined_FS_r2",
                "mean_function_score" = "mean_combined_fs",
                "final_function_score" = "function_score_sns",
                "q_value" = "fdr",
                "threshold" = "fs_threshold",
                "cell_line" = "sge_region",
                "treatment" = "expt",
                "function_class" = "HMEC_class",
                "ref" = "Ref"
                ) %>% 
  mutate(sge_region = "exon 12 (5')",
         conseq = case_when(conseq == "NON_SYNONYMOUS" ~ "Missense",
                            conseq == "STOP_GAINED" ~ "Nonsense",
                            conseq == "SYNONYMOUS" ~ "Synonymous",
                            conseq == "INTRONIC" ~ "Intronic",
                            conseq == "SPLICE_SITE" ~ "Splice region",
                            conseq == "CANONICAL_SPLICE" ~ "Canonical splice",
                            conseq == "5PRIME_UTR" ~ "5' UTR",
                            TRUE ~ conseq)) %>% 
  select(1, 2, 24, 3:23) %>% 
  dplyr::rename("consequence" = "conseq") %>% 
  select(-mean_function_score)

#write.csv(x = df_RNF168_export, file = paste0("~/Dropbox (The Francis Crick)/BRCA1_HAP1_HMEC_manuscript/SUPPLEMENTARY_TABLES/", date, "_RNF168_lines_SGE_scores.csv"), row.names = FALSE) 


