#Produces plots used in main and supplementary figures
#Requires files in data_for_plots

#### Setup ####
library(tidyverse)
library(readxl)
library(plotROC)

#Set working directory
#Ensure data_for_plots is saved in the working directory
setwd("path/to/working/directory/")
supplementary_tables <- "data_for_plots/SUPPLEMENTARY_TABLES_1-8.xlsx"

options(scipen=999)
threshold_lower <- -0.7991285
xy_label_size <- 6
text_size <- 6

##### Read in additional data #####
HAP1_thresholds <- read.csv("data_for_plots/HAP1_thresholds.csv")
HAP1_variants_unfiltered <- read.csv("data_for_plots/HAP1_variants_unfiltered.csv")
promoter_phyloP_scores <- read.csv("data_for_plots/promoter_phyloP_scores.csv")
NU7441_editing_rates <- read.csv("data_for_plots/NU7441_HDR_rates.csv")
CC_domain_foldX <- read.csv("data_for_plots/CC_domain_FoldX.csv")
HMEC_thresholds <- read.csv("data_for_plots/HMEC_thresholds.csv")
HMEC_variants_unfiltered <- read.csv("data_for_plots/HMEC_variants_unfiltered.csv")
BRCA1_essentiality <- read.csv("data_for_plots/BRCA1_essentiality.csv")
PALB2_essentiality <- read.csv("data_for_plots/PALB2_essentiality.csv")
RNF168_essentiality <- read.csv("data_for_plots/RNF168_essentiality.csv")
sge_data_2018 <- read.csv("data_for_plots/Findlay_et_al_2018_data.csv")

##### Read in supplementary tables #####
#Supplementary Table 1
ST1 <- read_excel(supplementary_tables,
                  sheet = "ST1",
                  range = cell_limits(c(4, 1), c(NA, NA)),
                  na = c("", "NA", "n/a"))

ST1$OddsPath <- as.numeric(ST1$OddsPath)
ST1$OddsPath_corrected <- as.numeric(ST1$OddsPath_corrected)

ST1 <- left_join(x = ST1,
                 y = HAP1_thresholds,
                 by = "SGE_region")

#Supplementary Table 2
ST2 <- read_excel(supplementary_tables,
                  sheet = "ST2",
                  range = cell_limits(c(4, 1), c(NA, NA)),
                  na = c("", "NA", "n/a"))

#Supplementary Table 3
ST3 <- read_excel(supplementary_tables,
                  sheet = "ST3",
                  range = cell_limits(c(4, 1), c(NA, NA)),
                  na = c("", "NA", "n/a"),
                  guess_max = 2157)

ST3$final_func_class <- factor(ST3$final_func_class, levels = c("Concordant neutral", "HAP1 intermediate, HMEC neutral", "HAP1 LoF, HMEC neutral", "Concordant LoF", "HAP1 neutral, HMEC depleted", "HAP1 intermediate, HMEC depleted"))

ST3 <- left_join(x = ST3,
                 y = HMEC_thresholds,
                 by = "SGE_region")


#Supplementary Table 4
ST4 <- read_excel(supplementary_tables,
                  sheet = "ST4",
                  range = cell_limits(c(4, 1), c(NA, NA)),
                  na = c("", "NA", "n/a"))

#Supplementary Table 5
ST5_HAP1 <- read_excel(supplementary_tables,
                  sheet = "ST5",
                  range = cell_limits(c(4, 2), c(13, 11)),
                  na = c("", "NA", "n/a"))

ST5_HAP1_meta <- ST5_HAP1[7:9,]

ST5_HAP1_HMEC <- read_excel(supplementary_tables,
                       sheet = "ST5",
                       range = cell_limits(c(30, 2), c(48, 11)),
                       na = c("", "NA", "n/a"))

ST5_HAP1_HMEC_meta <- ST5_HAP1_HMEC[13:18,]
ST5_HAP1_HMEC_meta$OR <- as.numeric(ST5_HAP1_HMEC_meta$OR)
ST5_HAP1_HMEC_meta$LCI <- as.numeric(ST5_HAP1_HMEC_meta$LCI)
ST5_HAP1_HMEC_meta$UCI <- as.numeric(ST5_HAP1_HMEC_meta$UCI)

#Supplementary Table 7
ST7 <- read_excel(supplementary_tables,
                       sheet = "ST7",
                       range = cell_limits(c(3, 1), c(NA, NA)),
                       na = c("", "NA", "n/a"))

ST7$`-LOG10(FDR)` <- as.numeric(ST7$`-LOG10(FDR)`)

#Supplementary Table 8
ST8 <- read_excel(supplementary_tables,
                  sheet = "ST8",
                  range = cell_limits(c(3, 1), c(NA, NA)),
                  na = c("", "NA", "n/a"))


##### Colour palettes #####

consequence_colours = c("Nonsense" = "#CC6677",
                       "Canonical splice" = "#DDCC77",
                       "Intronic" = "#646CFF",
                       "Missense" = "#326C37",
                       "Splice region" = "#66B1A3",
                       "Synonymous" = "#88CCEE",
                       "5' UTR" = "#929000",
                       "Promoter" = "gray80",
                       "Missense_grey" = "gray80")

indel_colours = c("5' UTR" = "#dd9f31",
                  "Promoter" = "#dd9f31",
                  "1" = "#2E3393",
                  "5" = "#AECDCE")

combined_function_class_colours <- c(
  "Concordant LoF" = "#E73577",
  "Concordant neutral" = "#4C77B9",
  "HAP1 intermediate, HMEC neutral" = "#88CCEE",
  "HAP1 LoF, HMEC neutral" = "#F7A73D",
  "HAP1 neutral, HMEC depleted" = "#0F0A3F",
  "HAP1 intermediate, HMEC depleted" = "#BABABA")

combined_function_class_colours2 <- c(
  "Concordant-LoF" = "#E73577",
  "Concordant-neutral" = "#4C77B9",
  "HAP1-int, HMEC-neutral" = "#88CCEE",
  "HAP1-LoF, HMEC-neutral" = "#F7A73D")

domain_colours <- c(
  "Coiled-coil" = "#F0E442",
  "RING" = "#CC79A7",
  "BRCT" = "#0F0A3F")

clinvar_colours = c("Absent" = "gray80",
                    "VUS/conflicting interpretations" = "#DCA237",
                    "Benign/Likely benign" = "#0433FF",
                    "Pathogenic/Likely pathogenic" = "#882255") 

HAP1_func_class_colours <- c(
  "Neutral" = "#4C77B9",
  "Intermediate" = "#E2C64D",
  "LoF" = "#E73577"
)

func_class_colours_HAP1_2018 <- c(
  "LOF" = "#E73577",
  "FUNC" = "#4C77B9",
  "INT" = "#E2C64D"
)

HMEC_func_class_colours <- c(
  "Neutral" = "#4C77B9",
  "Depleted" = "#E73577"
)

evidence_code_colours <- c("B_strong" = "steelblue4",
                           "B_moderate" = "steelblue3",
                           "B_supporting" = "steelblue1",
                           "indeterminate" = "snow2",
                           "P_supporting" = "rosybrown1",
                           "P_moderate" = "indianred2",
                           "P_strong" = "red3")

E2F_binding_colours <- c("TRUE" = "#C14023",
                         "FALSE" = "gray80")

NU7441_region_colours <- c(
  "promoter" = "#cc79a7",
  "exon 10 (3')" = "#f0e442",
  "exon 11" = "#009e73"
)

frame_colours <- c("In-frame" = "#83B25D",
                   "Frameshift" = "#305133")

significant_colours = c(
  "yes" = "#dd7676",
  "no" = "gray80"
)


#### MAIN FIGURES ####
##### Figure 1 #####
###### 1b ######
df_fig1b <- ST1 %>% 
  filter(SGE_region != "exon 5" & SGE_region != "exon 17")

df_fig1b$Consequence <- factor(df_fig1b$Consequence, levels = c("Synonymous", "Nonsense", "Missense", "Canonical splice", "Splice region", "Intronic", "Promoter", "5' UTR"))

ggplot(data = df_fig1b) +
  geom_jitter(mapping = aes(x = Consequence,
                            y = final_function_score,
                            colour = Consequence),
              size = 0.5) +
  ylab('Function score') +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = 'None'
  ) +
  scale_colour_manual(values = consequence_colours) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_discrete(labels = c("Synonymous", "Nonsense", "Missense", "Canonical\nsplice", "Splice region", "Intronic", "Promoter", "5' UTR"))


###### 1c ######
df_fig1c <- ST1 %>% 
  filter(SGE_region != "exon 5" & SGE_region != "exon 17") %>% 
  mutate(threshold_lower = case_when(threshold_upper < threshold_lower ~ NA_real_,
                                         TRUE ~ threshold_lower)) %>% 
  mutate(region = factor(SGE_region, levels = c("promoter", "5' UTR", "exon 5", "exon 6", "exon 10 (5')", "exon 10 (mid 1)", "exon 10 (mid 2)", "exon 10 (3')", "exon 11", "intron 11", "exon 12 (5')", "exon 12 (3')", "exon 17"))) 

custom_labels <- c(
  "exon 5" = "Exon 5",
  "intron 11" = "Intron 11",
  "exon 10 (5')" = "Exon 10 (5')",
  "exon 10 (3')" = "Exon 10 (3')",
  "exon 11" = "Exon 11",
  "exon 12 (5')" = "Exon 12 (5')",
  "exon 12 (3')" = "Exon 12 (3')",
  "exon 17" = "Exon 17",
  "exon 6" = "Exon 6 (5')",
  "promoter" = "Promoter",
  "5' UTR" = "5' UTR",
  "exon 10 (mid 2)" = "Exon 10 (mid 2)",
  "exon 10 (mid 1)" = "Exon 10 (mid 1)"
)

ggplot(df_fig1c, aes(y = final_function_score,
                                 x = hg38)) +
  geom_point(aes(colour = Consequence),
             size = 0.5) + 
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  ylab('Function score') +
  xlab('Chr. 17 position (hg38)') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  theme(strip.text = element_text(face = "bold",
                                  size = 7)) +
  scale_color_manual(values = consequence_colours) +
  geom_hline(aes(yintercept = threshold_upper), lty = 2) +
  geom_hline(aes(yintercept = threshold_lower), lty = 2, color = "gray70") +
  scale_x_reverse(breaks = function(x) round(seq(min(x)+30, max(x)-30, length.out = 2), -1)) + # Fixed 3 breaks
  facet_wrap(~region, ncol = 4, scales = "free_x", labeller = labeller(region = custom_labels)) # 3x4 layout (adjust `ncol` for different layouts)


##### Figure 2 #####
###### 2a ######
df_fig2a <- ST1 %>% 
  filter(SGE_region == "exon 11" | SGE_region == "exon 12 (5')" | SGE_region == "exon 12 (3')") %>% 
  filter(!is.na(pHGVS))


ggplot(df_fig2a,
       aes(y = final_function_score,
           x = as.numeric(CDSpos))) + 
  geom_rect(xmin = 4266.5, ymax = 1, xmax = 4176.5, ymin = -6, fill = "gray90") + #Coiled-coil domain
  geom_point(aes(colour = Consequence),
             size = 0.5) + 
  coord_cartesian(ylim = c(-4, 1)) +
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1)) +
  scale_x_continuous(breaks = c(4120, 4170, 4220, 4270, 4320),
  ) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.position = "None") + 
  ylab('Function score') + 
  xlab('CDS position') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) + 
  geom_segment(x = 4357, y = -0.332941287358147, xend = 4253, yend = -0.332941287358147, lty = 2) + #exon 12 (3') threshold
  geom_segment(x = 4252, y = -0.270354885547281, xend = 4186, yend = -0.270354885547281, lty = 2) + #exon 12 (5') threshold
  geom_segment(x = 4185, y = -0.261369546949105, xend = 4097, yend = -0.261369546949105, lty = 2) + #exon 11 threshold
  geom_segment(x = 4357, y = -0.798994339561149, xend = 4097, yend = -0.798994339561149, lty = 2, color = "gray70") + #lower threshold for all
  scale_color_manual(values = consequence_colours)

###### 2c ######
df_fig2c <- left_join(x = ST1,
                      y = CC_domain_foldX,
                      by = "pHGVS") %>% 
  filter(!is.na(FoldX_score),
         final_score == TRUE)

ggplot(data = df_fig2c) +
  geom_point(aes(x = final_function_score,
                 y = FoldX_score,
                 color = Consequence),
             size = 0.7) +
  ylab('FoldX score') +
  xlab('Function score') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_color_manual(values =  consequence_colours)

cor(df_fig2c$final_function_score, df_fig2c$FoldX_score, use = "pairwise.complete.obs", method = "spearman") 

###### 2d ######
df_fig2d <- ST1 %>% 
  filter(SGE_region != "exon 5" & SGE_region != "exon 17")

ggplot(data = df_fig2d,
       mapping = aes(x = final_function_score,
                     y = RNA_score,
                     color = Consequence)) +
  geom_point(size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.title=element_blank(),
        legend.key = element_blank(),
        legend.position = 'None'
  )+
  ylab('RNA score') +
  xlab('Function score') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  coord_cartesian(ylim = c(-6, 2.2)) +
  scale_color_manual(values = consequence_colours)


###### 2e ######
df_fig2e <- ST1 %>% 
  filter(SGE_region != "exon 5" & SGE_region != "exon 17")

ggplot(data = df_fig2e,
       mapping = aes(x = SpliceAI_max,
                     y = RNA_score,
                     color = Consequence)) +
  geom_point(size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.title=element_blank(),
        legend.key = element_blank(),
        legend.position='None'
  ) +
  ylab('RNA score') +
  xlab('Max SpliceAI score') +
  coord_cartesian(ylim = c(-6, 2.2)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_color_manual(values = consequence_colours)


###### 2f ######
df_fig2f <- ST1 %>% 
  filter(SGE_region == "exon 10 (5')") %>% 
  filter(hg38 < 43094781)

ggplot(df_fig2f,
       aes(y = RNA_score,
           x = hg38)) +
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  ylab('RNA score') +
  xlab('Chr. 17 position (hg38)') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(values = consequence_colours) +
  scale_x_reverse(breaks = c(43094775, 43094750))


###### 2g ######
df_fig2g <- ST1 %>% 
  filter(SGE_region == "promoter")

ggplot(df_fig2g,
       aes(y = final_function_score,
           x = hg38)) + 
  geom_rect(xmin = -43125412, xmax = -43125404, ymin = -3, ymax = 1, fill="gray90") +
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  coord_cartesian(ylim = c(-3, 1),
                  xlim = c(43125459, 43125337)) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0, 1)) +
  scale_x_reverse(breaks = c(43125440, 43125400, 43125360),
                  labels = c("43125440", "43125400", "43125360")) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('Function score') + 
  xlab('Chr. 17 position (hg38)') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_color_manual(values = indel_colours)


###### 2h ######

df_fig2h <- ST2 %>% 
  mutate(mid_point = (pos_first_deleted_hg38 + pos_last_deleted_hg38)/2)
df_fig2h$indel_size <- as.character(df_fig2h$indel_size)


ggplot(df_fig2h,
       aes(y = mean_function_score,
           x = mid_point)) + 
  geom_rect(xmax = -43125412, xmin = -43125404, ymin = -3, ymax = 1, fill="gray90") +
  geom_point(aes(colour = indel_size),
             size = 0.7) +
  coord_cartesian(ylim = c(-3, 1),
                  xlim = c(43125459, 43125337)) +
  scale_y_continuous(breaks = c(-3, -2, -1, 0, 1)) +
  scale_x_reverse(breaks = c(43125440, 43125400, 43125360),
                  labels = c("43125440", "43125400", "43125360")) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('Function score') + 
  xlab('Chr. 17 position (hg38)') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_color_manual(values = indel_colours)

###### 2i ######
df_fig2i <- ST1 %>% 
  filter(!is.na(E2F_binding_score))

ggplot(data = df_fig2i,
       mapping = aes(x = final_function_score,
                     y = E2F_binding_score)) +
  geom_point(size = 0.7,
             colour = "#dd9f31") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position='None'
  )+
  ylab('E2F binding score') +
  xlab('Function score')+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) 

##### Figure 3 #####
###### 3b ######
df_fig3b <- ST3

ggplot(data = df_fig3b,
       aes(x = final_function_score_ut,
           fill = Consequence)) +
  geom_histogram(alpha = 1, binwidth = 0.1) +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.position = 'None',
        legend.key = element_blank()
  ) +
  coord_cartesian(xlim = c(-4, 2)) +
  ylab('SNVs') +
  xlab('Function score') +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_fill_manual(values = consequence_colours) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2))

###### 3c ######
df_fig3c <- ST3

ggplot(data = df_fig3c,
       aes(y = final_function_score_ut,
           x = HAP1_function_score_mean)) + 
  geom_point(aes(colour = final_func_class),
             size = 0.5) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.position = "None") + 
  ylab('HMEC function score') + 
  xlab('HAP1 function score') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0)) +
  scale_color_manual(values = combined_function_class_colours)

###### 3d ######
df_fig3d <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted") %>% 
  filter(SGE_region != "promoter") %>%  #don't use 5'UTR for this as not many variant selected for HMEC assay, proportion is a bit irrelevant
  dplyr::group_by(Consequence, final_func_class) %>% 
  dplyr::summarise(n = n()) %>%
  mutate(pct = n/sum(n)*100)
df_fig3d$Consequence <- factor(df_fig3d$Consequence, levels = c("Synonymous", "Nonsense",  "Missense", "Canonical splice", "Splice region", "Intronic"))

df_fig3d_totals <- df_fig3d %>% 
  dplyr::group_by(Consequence) %>% 
  dplyr::summarise(total = sum(n))

ggplot(df_fig3d,
       aes(x = Consequence, y = pct, fill = final_func_class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  coord_cartesian(ylim = c(0, 105)) +
  scale_fill_manual(values = combined_function_class_colours) +
  scale_x_discrete(expand = c(0.1, 0.1), labels = c("Synonymous", "Nonsense",  "Missense", "Canonical splice", "Splice region", "Intronic")) +
  geom_text(data = df_fig3d_totals, aes(x = Consequence, y = 100, label = total),
            vjust = -0.5, size = 5, inherit.aes = FALSE, size.unit = "pt")

###### 3e ######
df_fig3e <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted") %>%
  mutate(category = case_when(SGE_region == "exon 3" | SGE_region == "exon 5" ~ "RING",
                              SGE_region == "exon 17" ~ "BRCT",
                              protPos >=1393 & protPos <= 1422 ~ "CC",
                              SGE_region == "exon 10 (5')" | SGE_region == "exon 10 (3')" | SGE_region == "exon 10 (mid 1)" | SGE_region == "exon 10 (mid 2)" ~ "Exon 10",
                              TRUE ~ "Other"
  )) %>% 
  filter(Consequence == "Missense") %>% 
  filter(category != "Other") %>%
  dplyr::group_by(category, final_func_class) %>% 
  dplyr::summarise(n = n()) %>%
  mutate(pct = n/sum(n)*100)
df_fig3e$category <- factor(df_fig3e$category, levels = c("BRCT", "RING", "CC", "Exon 10"))

df_fig3e_totals <- df_fig3e %>% 
  group_by(category) %>% 
  summarise(total = sum(n)) %>% 
  filter(!is.na(category))

ggplot(df_fig3e,
  aes(x = category, y = pct, fill = final_func_class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none" 
  ) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  coord_cartesian(ylim = c(0, 105)) +
  scale_fill_manual(values = combined_function_class_colours) +
  geom_text(data = df_fig3e_totals, aes(x = category, y = 100, label = total),
            vjust = -0.5, size = 5, inherit.aes = FALSE, size.unit = "pt")


###### 3f ######
df_fig3f <- ST3 %>% 
  filter(SGE_region == "exon 11" | SGE_region == "exon 12 (5')")

ggplot(data = df_fig3f,
       aes(y = final_function_score_ut,
           x = HAP1_function_score_mean)) + 
  geom_point(aes(colour = Consequence),
             size = 0.7) + 
 theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.position = "None") + 
  labs(title = "Exons 11 and 12, n = 534") +
  ylab('HMEC function score') + 
  xlab('HAP1 function score') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = consequence_colours)


###### 3g ######
df_fig3g <- ST3 %>% 
  filter(SGE_region == "intron 11")

ggplot(df_fig3g,
       aes(y = final_function_score_ut,
           x = HAP1_function_score_mean)) + 
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.position = "None") + 
  labs(title = "Intron 11, n = 20") +
  ylab('HMEC function score') + 
  xlab('HAP1 function score') + 
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = consequence_colours) + 
  geom_text(x = -2.31990869679722 + 0.8, y = -1.64223597670913, label = "c.4185+4108C>G", size = text_size, size.unit = "pt") +
  geom_text(x = -1.37295907396138 + 0.8, y = -0.912227959661949, label = "c.4185+4105C>T", size = text_size, size.unit = "pt")

##### Figure 4 #####
###### 4a ######
df_fig4a <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted")

ggplot(data = df_fig4a,
       aes(y = CADD.phred,
           x = final_func_class)) +
  geom_jitter(aes(colour = final_func_class),
              size = 0.7) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#000000") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  ylab('CADD score') +
  xlab('') +
  ggtitle('CADD, n = 2,146') +
  theme(plot.title = element_text(size = 7, hjust = 0.5, face = "bold")) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  theme(legend.position = "None") +
  scale_colour_manual(values = combined_function_class_colours) + 
  scale_x_discrete(labels = c("Concordant neutral" = "c-neutral", "HAP1 intermediate, HMEC neutral" = "HAP1 int,\nHMEC \nneutral", "HAP1 LoF, HMEC neutral" = "HAP1 LoF,\nHMEC \nneutral", "Concordant LoF" = "c-LoF"))


###### 4b ######
df_fig4b <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted") %>% 
  filter(!is.na(AlphaMissense_score))

ggplot(data = df_fig4b,
       aes(y = AlphaMissense_score,
           x = final_func_class)) +
  geom_jitter(aes(colour = final_func_class),
              size = 0.7) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#000000") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  ylab('AlphaMissense pathogenicity score') +
  xlab('') +
  ggtitle('AlphaMissense, n = 1,200') +
  theme(plot.title = element_text(size = 7, hjust = 0.5, face = "bold")) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None"
  ) +
  scale_x_discrete(labels = c("Concordant neutral" = "c-neutral", "HAP1 intermediate, HMEC neutral" = "HAP1 int,\nHMEC \nneutral", "HAP1 LoF, HMEC neutral" = "HAP1 LoF,\nHMEC \nneutral", "Concordant LoF" = "c-LoF")) +
  scale_colour_manual(values = combined_function_class_colours)


###### 4c ######
df_fig4c <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted") %>% 
  filter(!is.na(EVE_score))

ggplot(data = df_fig4c,
      aes(y = EVE_score,
          x = final_func_class)) +
  geom_jitter(aes(colour = final_func_class),
              size = 0.7) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#000000") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  ylab('EVE score') +
  xlab('') +
  ggtitle('EVE, n = 1,200') +
  theme(plot.title = element_text(size = 7, hjust = 0.5, face = "bold")) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank(),
        legend.position = "none") +
  theme(legend.title = element_blank()) + 
  scale_x_discrete(labels = c("Concordant neutral" = "c-neutral", "HAP1 intermediate, HMEC neutral" = "HAP1 int,\nHMEC \nneutral", "HAP1 LoF, HMEC neutral" = "HAP1 LoF,\nHMEC \nneutral", "Concordant LoF" = "c-LoF")) +
  scale_colour_manual(values = combined_function_class_colours)

###### 4d ######
df_fig4d <- ST3 %>% 
  mutate(domain = case_when(SGE_region == "exon 3" | SGE_region == "exon 5" ~ "RING",
                            SGE_region == "exon 17" ~ "BRCT",
                            protPos >=1393 & protPos <= 1422 ~ "Coiled-coil",
                            TRUE ~ "Other")) %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted",
         SpliceAI_max < 0.2,
         domain == "RING" | domain == "Coiled-coil" | domain == "BRCT",
         !is.na(FoldX_score))
  
ggplot(data =  df_fig4d,
       aes(x = final_func_class, y = FoldX_score)
) +
  geom_jitter(aes(colour = domain),
              size = 0.7) + 
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#000000") +
  ylab('FoldX score') +
  xlab('') +
  ggtitle('FoldX, n = 589') +
  theme(plot.title = element_text(size = 7, hjust = 0.5, face = "bold")) +
  coord_cartesian(ylim = c(-3, 25)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.title = element_blank()) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank(),
        legend.position = "none"
  ) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
  scale_x_discrete(labels = c("Concordant neutral" = "c-neutral", "HAP1 intermediate, HMEC neutral" = "HAP1 int,\nHMEC \nneutral", "HAP1 LoF, HMEC neutral" = "HAP1 LoF,\nHMEC \nneutral", "Concordant LoF" = "c-LoF")) +
  scale_colour_manual(values = domain_colours)

###### 4e ######
df_fig4e <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted",
         Consequence == "Synonymous" | Consequence == "Canonical splice" | Consequence == "Splice region" | Consequence == "Intronic")

ggplot(data = df_fig4e,
       aes(y = SpliceAI_max,
           x = final_func_class)) +
  geom_jitter(aes(color = Consequence),
              size = 0.7) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#000000") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  ylab('SpliceAI max score') +
  ggtitle('SpliceAI, n = 827') +
  theme(plot.title = element_text(size = 7, hjust = 0.5, face = "bold")) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        axis.title.x = element_blank(),
        strip.background = element_blank()) +
  scale_x_discrete(labels = c("Concordant neutral" = "c-neutral", "HAP1 intermediate, HMEC neutral" = "HAP1 int,\nHMEC \nneutral", "HAP1 LoF, HMEC neutral" = "HAP1 LoF,\nHMEC \nneutral", "Concordant LoF" = "c-LoF")) +
  scale_color_manual(values = consequence_colours) +
  theme(legend.position = "None")        

###### 4f ######
df_fig4f <- ST3 %>% 
  mutate(final_func_class_combined = case_when(final_func_class == "HAP1 LoF, HMEC neutral" | final_func_class == "HAP1 intermediate, HMEC neutral" ~ "HAP1-HMEC discordant",
                                               TRUE ~ final_func_class)) %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted",
         total_assays > 2)
df_fig4f$final_func_class_combined <- factor(df_fig4f$final_func_class_combined, levels = c("Concordant neutral", "HAP1-HMEC discordant", "Concordant LoF"))

ggplot(data = df_fig4f,
       aes(y = percent_abnormal,
           x = final_func_class_combined)) +
  geom_jitter(aes(colour = final_func_class),
              size = 0.7) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#000000") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  ylab('Percentage of assays for\nwhich function abnormal') +
  xlab('') +
  ggtitle('Function in other BRCA1 assays\nn = 98') +
  theme(plot.title = element_text(size = 7, hjust = 0.5, face = "bold")) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  theme(legend.position = "None") +
  scale_color_manual(values = combined_function_class_colours) +
  scale_x_discrete(labels = c("Concordant neutral" = "c-neutral", "HAP1-HMEC, discordant" = "HAP1-HMEC,\ndiscordant", "Concordant LoF" = "c-LoF"))

##### Figure 5 #####
###### 5a ######
df_fig5a <- ST3 %>% 
  mutate(plot_alpha = case_when(HMEC_class_ut == "Neutral" & HMEC_class_olaparib == "Depleted" ~ 1,
                                TRUE ~ 0.2)) %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted",
         final_func_class != "HAP1 intermediate, HMEC depleted",
         !is.na(final_function_score_ut) & !is.na(final_function_score_olaparib)) 

ggplot() + 
  geom_point(data = df_fig5a,
             aes(x = final_function_score_ut,
                 y = final_function_score_olaparib, 
                 colour = final_func_class,
                 alpha = plot_alpha),
             size = 0.5
  ) + 
  scale_alpha_identity() +
  coord_cartesian(xlim = c(-3.5, 1.7),
                  ylim = c(-3.5, 1.7)
  ) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.position = "None") + 
  ylab('HMEC function score (+ olaparib)') + 
  xlab('HMEC function score (untreated)') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  geom_abline(slope = 1, lty = 2) +
  scale_color_manual(values = combined_function_class_colours)

cor(df_fig5a$final_function_score_ut, df_fig5a$final_function_score_olaparib, use = "pairwise.complete.obs", method = "pearson") 

###### 5b ######
df_fig5bc <- ST3 %>% 
  filter(SGE_region == "exon 10 (3')")

ggplot(data = df_fig5bc,
       aes(y = final_function_score_ut, 
           x = hg38)) + 
  geom_point(aes(colour = Consequence),
             size = 0.5) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('HMEC function score') + 
  xlab('Chr. 17 position (hg38)') + 
  labs(title = "Untreated") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  geom_hline(yintercept = df_fig5bc$HMEC_threshold_ut, lty = 2) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
        strip.background = element_blank()) + 
  scale_color_manual(values = consequence_colours) + 
  scale_x_reverse(breaks = c(43091500, 43091420))

###### 5c ######
ggplot(data = df_fig5bc,
       aes(y = final_function_score_olaparib, 
           x = hg38)) + 
  geom_point(aes(colour = Consequence),
             size = 0.5) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('HMEC function score') + 
  xlab('Chr. 17 position (hg38)') + 
  labs(title = "Untreated") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  geom_hline(yintercept = df_fig5bc$HMEC_threshold_olaparib, lty = 2) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
        strip.background = element_blank()) + 
  scale_color_manual(values = consequence_colours) + 
  scale_x_reverse(breaks = c(43091500, 43091420))


###### 5e ######
#join HAP1 function scores
join_to_ST4 <- ST1 %>% 
  filter(final_score == TRUE) %>% 
  select(cHGVS, function_class) %>% 
  dplyr::rename("HAP1_function_class" = "function_class")

df_fig5ef <- left_join(x = ST4,
                       y = join_to_ST4,
                       by = "cHGVS") %>%
  #colour missense variants that were not LoF in HAP1 grey
  mutate(consequence_for_plot = case_when(consequence == "Missense" & HAP1_function_class != "LoF" ~ "Missense_grey",
                                          TRUE ~ consequence))

df_fig5e <- df_fig5ef %>% 
  filter(treatment == "untreated") %>% 
  pivot_wider(names_from = "cell_line",
              values_from = c("final_function_score", "consequence_for_plot"),
              id_cols = "cHGVS") %>% 
  select(1:4) %>% 
  dplyr::rename("consequence_for_plot" = "consequence_for_plot_HMEC RNF168-KO")

ggplot() + 
  geom_point(data = df_fig5e,
             aes(y = `final_function_score_HMEC RNF168-KO`, 
                 x = `final_function_score_HMEC RNF168-WT`,
                 colour = consequence_for_plot),
             size = 0.5) + 
  labs(title = "Untreated") +
  coord_cartesian(ylim = c(-2.1, 0.7),
                  xlim = c(-2.1, 0.7)) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('HMEC function score (RNF168-/-)') + 
  xlab('HMEC function score (RNF168+/+)') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
        strip.background = element_blank()) + 
  geom_abline(slope = 1, lty = 2) +
  scale_color_manual(values = consequence_colours)


###### 5f ######
df_fig5f <- df_fig5ef %>% 
  filter(treatment == "olaparib") %>% 
  pivot_wider(names_from = "cell_line",
              values_from = c("final_function_score", "consequence_for_plot"),
              id_cols = "cHGVS") %>% 
  select(1:4) %>% 
  dplyr::rename("consequence_for_plot" = "consequence_for_plot_HMEC RNF168-KO")

ggplot() + 
  geom_point(data = df_fig5f,
             aes(y = `final_function_score_HMEC RNF168-KO`, 
                 x = `final_function_score_HMEC RNF168-WT`,
                 colour = consequence_for_plot),
             size = 0.5) + 
  labs(title = "+ olaparib") +
  coord_cartesian(ylim = c(-2.1, 0.7),
                  xlim = c(-2.1, 0.7)) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('HMEC function score (RNF168-/-)') + 
  xlab('HMEC function score (RNF168+/+)') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
        strip.background = element_blank()) + 
  geom_abline(slope = 1, lty = 2) +
  scale_color_manual(values = consequence_colours)


##### Figure 6 #####
###### 6a ######

#variants with B/P ClinVar classification with minimum one star review status and scored in both cell lines
df_fig6ab <- ST3 %>% 
  filter(ClinVar_simple != "Absent") %>%
  select(cHGVS, Consequence, SGE_region, ClinVar_simple, ClinVar_review_status, final_function_score_ut, HAP1_function_score_mean) %>% 
  filter(ClinVar_review_status != "no assertion criteria provided" & ClinVar_review_status != "no classification provided" & !is.na(ClinVar_review_status)) %>% 
  mutate(classification = case_when(grepl(pattern = "enign", x = ClinVar_simple) ~ 0,
                                    grepl(pattern = "athogenic", x = ClinVar_simple) ~ 1)) %>% 
  filter(!is.na(classification))

fig6a <- ggplot(df_fig6ab,
                aes(d = classification,
                    m = -HAP1_function_score_mean)) +
  geom_roc(n.cuts=0,cutoffs.at = NULL) +
  style_roc(xlab = "1 - Specificity",
            major.breaks = c(0, 0.25, 0.5, 0.75, 1),
            minor.breaks = NULL) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank())

fig6a_annotated <- fig6a + annotate("text", x = .5, y = .5, label = paste("AUC =", round(calc_auc(fig6a)$AUC, 3)), size = 6 * 0.35)
fig6a_annotated


###### 6b ######
fig6b <- ggplot(df_fig6ab,
                aes(d = classification,
                    m = -final_function_score_ut)) +
  geom_roc(n.cuts=0,cutoffs.at = NULL) +
  style_roc(xlab = "1 - Specificity",
            major.breaks = c(0, 0.25, 0.5, 0.75, 1),
            minor.breaks = NULL) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank())

fig6b_annotated <- fig6b + annotate("text", x = .5, y = .5, label = paste("AUC =", round(calc_auc(fig6b)$AUC, 3)), size = 6 * 0.35)
fig6b_annotated

###### 6c ######
df_fig6c <- ST3 %>% 
  #add ClinVar info for p.R1699Q
  mutate(ClinVar_simple = case_when(pHGVS == "p.R1699Q" ~ "Pathogenic/Likely pathogenic",
                                    TRUE ~ ClinVar_simple))

ggplot(data = df_fig6c,
       aes(y = final_function_score_ut,
           x = HAP1_function_score_mean)) + 
  geom_point(aes(colour = ClinVar_simple),
             size = 0.5) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None"
  ) + 
  ylab('HMEC function score') + 
  xlab('HAP1 function score') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0)) +
  scale_color_manual(values = clinvar_colours)

###### 6d ######
df_fig6d <- df_fig6c %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted" & final_func_class != "HAP1 intermediate, HMEC depleted") %>% 
  group_by(final_func_class, ClinVar_simple) %>% 
  summarise(n = n()) %>% 
  group_by(final_func_class) %>% 
  mutate(pct = n/sum(n)*100) %>% 
  mutate(final_func_class_renamed = case_when(final_func_class == "Concordant neutral" ~ "Concordant-neutral",
                                              final_func_class == "Concordant LoF" ~ "Concordant-LoF",
                                              final_func_class == "HAP1 LoF, HMEC neutral" ~ "HAP1-LoF, HMEC-neutral",
                                              final_func_class == "HAP1 intermediate, HMEC neutral" ~ "HAP1-int, HMEC-neutral"
                                              ))

df_fig6d_totals <- df_fig6d %>% 
  group_by(final_func_class_renamed) %>% 
  summarise(total = sum(n))
  
df_fig6d$final_func_class_renamed <- factor(df_fig6d$final_func_class_renamed, levels = c("Concordant-neutral", "HAP1-int, HMEC-neutral", "HAP1-LoF, HMEC-neutral", "Concordant-LoF"))
df_fig6d$ClinVar_simple <- factor(df_fig6d$ClinVar_simple, levels = c("Absent", "VUS/conflicting interpretations", "Benign/Likely benign", "Pathogenic/Likely pathogenic"))

ggplot(df_fig6d,
       aes(x = final_func_class_renamed,
           y = pct,
           fill = ClinVar_simple)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),                   
    axis.title.x = element_blank(),                 
    axis.title.y = element_blank(),                
    legend.title = element_blank(),                    
    legend.position = "none"
  ) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black")) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  coord_cartesian(ylim = c(0, 104)) +
  scale_fill_manual(values = clinvar_colours) + 
  geom_text(data = df_fig6d_totals, aes(x = final_func_class_renamed, y = 100, label = total),
            vjust = -0.5, size = 5*0.35, inherit.aes = FALSE) 


###### 6e ######
df_fig6e <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted" & final_func_class != "HAP1 intermediate, HMEC depleted") %>% 
  group_by(final_func_class, present_in_gnomad) %>% 
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)*100) %>% 
  filter(present_in_gnomad == TRUE)

ggplot(df_fig6e,
       aes(x = final_func_class, y = pct,
           fill = final_func_class)) +
  geom_bar(stat = "identity") +
  labs(title = "gnomAD") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),                  
    axis.title.y = element_blank(),                   
    plot.title = element_text(hjust = 0.5, size = text_size),
    legend.position = "none"
  ) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black")) +
  geom_text(data = df_fig6e, aes(x = final_func_class, y = pct + 1, label = n),
            vjust = -0.5, size = 5*0.35, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 23)) +
  scale_fill_manual(values = combined_function_class_colours)

###### 6f ######
df_fig6f <- ST3 %>% 
  filter(final_func_class != "HAP1 neutral, HMEC depleted" & final_func_class != "HAP1 intermediate, HMEC depleted") %>% 
  group_by(final_func_class, present_in_AOU) %>% 
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)*100) %>% 
  filter(present_in_AOU == TRUE)

ggplot(df_fig6f,
       aes(x = final_func_class, y = pct,
           fill = final_func_class)) +
  geom_bar(stat = "identity") +
  labs(title = "All of Us") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),                  
    axis.title.y = element_blank(),                   
    plot.title = element_text(hjust = 0.5, size = text_size),
    legend.position = "none"
  ) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black")) +
  geom_text(data = df_fig6f, aes(x = final_func_class, y = pct + 1, label = n),
            vjust = -0.5, size = 5*0.35, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 14)) +
  scale_fill_manual(values = combined_function_class_colours)

###### 6g ######
ST5_HAP1_meta$`HAP1 function class` <- factor(ST5_HAP1_meta$`HAP1 function class`, levels = c("LoF", "Intermediate", "Neutral"))

ggplot(ST5_HAP1_meta) +
  geom_pointrange(aes(x = `HAP1 function class`,
                      y = OR,
                      ymin = LCI,
                      ymax = UCI,
                      colour = `HAP1 function class`)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  coord_flip(ylim = c(0, 12)) + 
  labs(x = "",
       y = "Odds ratio (95% CI)",
  ) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 6),
        strip.background = element_blank()) +
  scale_color_manual(values = HAP1_func_class_colours,
                     breaks = names(HAP1_func_class_colours)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12))


###### 6h ######
df_fig6h <- ST5_HAP1_HMEC_meta %>% 
  mutate(final_func_class = case_when(`Combined function class` == "Concordant neutral" ~ "Concordant-neutral",
                                      `Combined function class` == "HAP1 int, HMEC neutral" ~ "HAP1-int, HMEC-neutral",
                                      `Combined function class` == "HAP1 LoF, HMEC neutral" ~ "HAP1-LoF, HMEC-neutral",
                                      `Combined function class` == "Concordant LoF" ~ "Concordant-LoF",
                                      TRUE ~ NA_character_)) %>% 
  filter(!is.na(final_func_class))

df_fig6h$final_func_class <- factor(df_fig6h$final_func_class, levels = c("Concordant-LoF", "HAP1-LoF, HMEC-neutral", "HAP1-int, HMEC-neutral", "Concordant-neutral"))

ggplot(df_fig6h) +
  geom_pointrange(aes(x = final_func_class,
                      y = OR,
                      ymin = LCI,
                      ymax = UCI,
                      colour = final_func_class)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  coord_flip() + 
  labs(x = "",
       y = "Odds ratio (95% CI)",
  ) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 6),
        strip.background = element_blank()) +
  scale_color_manual(values = combined_function_class_colours2,
                     breaks = names(combined_function_class_colours2))


###### 6i ######
df_fig6i <- ST1 %>% 
  filter(final_score == TRUE)

ggplot(data = df_fig6i) +
  geom_histogram(aes(x = final_function_score,
                     fill = evidence_code_adjusted),
                 binwidth = 0.05) +
  xlab('HAP1 function score') +
  ylab('') +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None"
  ) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_fill_manual(values = evidence_code_colours)

###### 6j ######
df_fig6j <- ST3

ggplot(data = df_fig6j) +
  geom_histogram(aes(x = final_function_score_ut,
                     fill = HMEC_evidence_code_adjusted),
                 binwidth = 0.05) +
  xlab('HMEC function score') +
  ylab('') +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None"
  ) +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_fill_manual(values = evidence_code_colours)


#### SUPPLEMENTARY FIGURES ####
##### Figure S1 #####
###### S1a ######

ggplot(data = HAP1_variants_unfiltered,
       mapping = aes(x = r1_function_score,
                     y = r2_function_score,
                     color = Consequence)) +
  geom_point(size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = 'None') +
  ylab('Rep. 2 function score') +
  xlab('Rep. 1 function score') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size),
        axis.text.x = element_text(size = xy_label_size),
        strip.background = element_blank()) +
  geom_abline(slope = 1,
              lty = 2) +
  coord_cartesian(xlim = c(-6, 2),
                  ylim = c(-6, 2)) +
  scale_color_manual(values = consequence_colours)


###### S1b ######
df_figs1b <- ST1 %>% 
  filter(SGE_region == "exon 17")

df_figs1b <- inner_join(x = df_figs1b,
                        y = sge_data_2018,
                        by = "cHGVS")

ggplot(df_figs1b,
       aes(y = function.score.mean, 
           x = final_function_score)) + 
  geom_point(aes(colour = consequence),
             size = 0.7) + 
  coord_cartesian(ylim = c(-3.5, 0.5),
                  xlim = c(-3.5, 0.5)) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('Exon 17 function scores (Findlay et al. 2018)') + 
  xlab('New exon 17 function score') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) + 
  scale_color_manual(values = consequence_colours) + 
  geom_vline(xintercept = df_figs1b$threshold_upper, lty = 2) + 
  geom_vline(xintercept = threshold_lower, lty = 2, colour = "gray70") +
  geom_hline(yintercept = -0.748, lty = 2) +
  geom_hline(yintercept = -1.328, lty = 2, colour = "gray70")


###### S1c ######
df_figs1c <- ST1 %>% 
  filter(SGE_region != "exon 5" & SGE_region != "exon 17") %>% 
  filter(final_score == TRUE) %>% 
  filter(Consequence == "Missense")

ggplot(df_figs1c,
       aes(x = final_function_score, 
           fill = function_class)) + 
  geom_histogram(alpha = 1, binwidth = 0.1) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.position = 'None',
        legend.key = element_blank()
  ) +
  ylab('SNVs') + 
  xlab('Function score') + 
  ggtitle('New regions') +
  theme(plot.title = element_text(size = 6, hjust = 0.5, face = "bold")) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) + 
  scale_fill_manual(values = HAP1_func_class_colours) + 
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1)) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600)) +
  coord_cartesian(xlim = c(-4.5, 1))

###### S1d ######
df_figs1d <- sge_data_2018 %>% 
  filter(consequence == "Missense")

ggplot(df_figs1d,
       aes(x = function.score.mean, 
           fill = func.class)) +
  geom_histogram(alpha = 1, binwidth = 0.1) +  
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.position = 'None',
        legend.key = element_blank()) +
  ylab('SNVs') + 
  xlab('Function score') + 
  ggtitle('RING/BRCT domains\n(Findlay et al. 2018)') +
  theme(plot.title = element_text(size = 6, hjust = 0.5, face = "bold")) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) + 
  scale_fill_manual(values = func_class_colours_HAP1_2018) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1)) +
  coord_cartesian(xlim = c(-4.5, 1))

##### Figure S2 #####
sge_regions <- unique(as.character(ST1$SGE_region))

df_figs2 <- ST1
df_figs2$rep_correlation <- 0

calc_cor_per_region <- function(region) {
  df <- ST1 %>%
    filter(SGE_region == region) %>% 
    filter(!is.na(r1_function_score) & !is.na(r2_function_score))
  correlation <- cor(df$r1_function_score, df$r2_function_score, use = "pairwise.complete.obs", method = "pearson")
  return(correlation)
}

for (reg in sge_regions) {
  output <- calc_cor_per_region(reg)
  print(output)
  df_figs2$rep_correlation[which(df_figs2$SGE_region == reg)] <- output
}

df_figs2 <- df_figs2 %>% 
  mutate(SGE_region = factor(SGE_region, levels = c("promoter", "5' UTR", "exon 5", "exon 6", "exon 10 (5')", "exon 10 (mid 1)", "exon 10 (mid 2)", "exon 10 (3')", "exon 11", "intron 11", "exon 12 (5')", "exon 12 (3')", "exon 17"))) %>% 
  mutate(label = sprintf("italic(r) == '%.2f'", rep_correlation))

ggplot(df_figs2,
       aes(y = r1_function_score, 
           x = r2_function_score)) +
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  ylab('Rep. 1 function score') +
  xlab('Rep. 2 function score') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  scale_color_manual(values = consequence_colours) +
  scale_y_continuous(breaks = c(-6, -5, -4, -3, -2, -1, 0, 1)) +
  scale_x_continuous(breaks = c(-6, -5, -4, -3, -2, -1, 0, 1)) +
  geom_abline(slope = 1, lty = 2) +
  theme(strip.text = element_text(face = "bold",
                                  size = 7)) +
  facet_wrap(~SGE_region, ncol = 4, scales = "fixed", labeller = labeller(SGE_region = custom_labels)) +
  geom_text(
    aes(label = label),
    x = -6.5, y = 0, 
    parse = TRUE, 
    hjust = -0.1, vjust = 1.1,
    size = 7*0.35
  )


##### Figure S3 #####
###### S3b ######
df_figs3b <- sge_data_2018 %>% 
  filter(experiment == "X5")


ggplot(df_figs3b,
       aes(y = function.score.mean, 
           x = hg38)) + 
  geom_point(aes(colour = consequence),
             size = 0.7) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('Function score') + 
  xlab('Chr. 17 position (hg38)') + 
  labs(title = "Exon 5 (Findlay et al. 2018)") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 0.5,
                                  face = "bold",
                                  size = text_size),
        strip.background = element_blank()) + 
  scale_color_manual(values = consequence_colours) + 
  geom_hline(yintercept = -0.748, lty = 2) +
  geom_hline(yintercept = -1.328, lty = 2, colour = "gray70") +
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0)) +
  scale_x_reverse(breaks = c(43104960, 43104930, 43104900, 43104870))


###### S3c ######
df_figs3c <- ST1 %>% 
  filter(SGE_region == "exon 5")

df_figs3c <- inner_join(x = df_figs3c,
                        y = sge_data_2018,
                        by = "cHGVS")

ggplot(df_figs3c,
       aes(y = function.score.mean, 
           x = final_function_score)) + 
  geom_point(aes(colour = consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('Exon 5 function scores (Findlay et al. 2018)') + 
  xlab('New exon 5 function scores') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) + 
  scale_color_manual(values = consequence_colours) + 
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0)) +
  geom_vline(xintercept = threshold_lower, lty = 2, colour = "gray70") +
  geom_vline(xintercept = df_figs3c$threshold_upper, lty = 2) +
  geom_hline(yintercept = -0.748, lty = 2) +
  geom_hline(yintercept = -1.328, lty = 2, colour = "gray70")


###### S3d ######
df_figs3d <- ST1 %>% 
  filter(SGE_region == "exon 5")

ggplot(df_figs3d,
       aes(y = final_function_score, 
           x = hg38)) + 
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('Function score') + 
  xlab('Chr. 17 position (hg38)') + 
  labs(title = "Exon 5 (new)") +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 0.5,
                                  face = "bold",
                                  size = text_size),
        strip.background = element_blank()) + 
  scale_color_manual(values = consequence_colours) + 
  geom_hline(yintercept = df_figs3d$threshold_upper, lty = 2) +
  geom_hline(yintercept = threshold_lower, lty = 2, colour = "gray70") +
  scale_x_reverse(breaks = c(43104960, 43104930, 43104900, 43104870))


##### Figure S4 #####
###### S4b ######
ggplot() +
  geom_point(data = promoter_phyloP_scores,
             aes(x = final_function_score,
                 y = verPhyloP,
                 colour = binding_region),
             size = 0.7) +
  geom_point(data = promoter_phyloP_scores %>% filter(binding_region == "TRUE"),
             aes(x = final_function_score,
                 y = verPhyloP,
                 colour = binding_region),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  xlab('Function score') +
  ylab('Vertebrate PhyloP score') +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_color_manual(values = E2F_binding_colours)


###### S4c ###### 
df_figs4c <- ST1 %>% 
  filter(SGE_region == "intron 11")

ggplot(data = df_figs4c,
       aes(y = r1_function_score,
           x = r2_function_score)) + 
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('Rep. 1 function score') + 
  xlab('Rep. 2 function score') + 
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_y_continuous(breaks = c(-2, -1, 0)) +
  scale_x_continuous(breaks = c(-2, -1, 0)) +
  geom_abline(slope = 1, lty = 2) +
  scale_color_manual(values = consequence_colours)

###### S4d ######
df_figs4d <- ST1 %>% 
  filter(SGE_region == "intron 11")

ggplot(data = df_figs4d,
       mapping = aes(x = AG, 
                     y = DG)) +
  geom_point(aes(colour = function_class),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.title=element_blank(),
        legend.key = element_blank(),
        legend.position = 'None') +
  coord_cartesian(ylim = c(0, 0.3),
                  xlim = c(0, 0.3)) +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3)) +
  scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3)) +
  ylab('SpliceAI donor gain score') +
  xlab('SpliceAI acceptor gain score') +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_color_manual(values = HAP1_func_class_colours)


##### Figure S5 #####
###### S5b ######

ggplot(BRCA1_essentiality, 
       aes(x = frame, y = log2(ratio_normalised))) +
  geom_jitter(aes(colour = frame),
              size = 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black") +
  coord_cartesian(ylim = c(-9, 3.3)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.title.x = element_blank(),     
        legend.title = element_blank(),
        legend.position = "None") +
  theme(text = element_text(size = 6),
        axis.text.y = element_text(size = 6, colour = "#000000"),
        axis.text.x = element_text(size = 6, colour = "#000000"),
        strip.background = element_blank()) + 
  scale_colour_manual(values = frame_colours) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Frameshift", "In-frame")), 
                     label = "p.signif",  
                     bracket.size = 0.5,  
                     tip.length = 0.02,   # 
                     label.y = 2.2) + #
  scale_y_continuous(breaks = c(-7.5, -5, -2.5, 0, 2.5)) +
  ylab("Indel function score")


###### S5d ###### 
ggplot(NU7441_editing_rates, 
       aes(x = factor(Concentration), y = HDR_rate_pct)) +
  geom_line(aes(group = SGE_region,
                colour = SGE_region)) +
  geom_point(aes(colour = SGE_region),
             size = 1) + 
  labs(x = "NU7441 concentration (nM)", y = "HDR editing rate (%)") +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  scale_color_manual(values = NU7441_region_colours)

###### S5f ###### 
genes_to_label = c("BRCA1", "BRCA2", "RAD51", "PALB2", "SLFN11", "RNF168")

df_s5f <- ST7 %>% 
  mutate(significant = case_when(FDR < 0.05 ~ "yes",
                                 TRUE ~ "no"),
         to_label = case_when(Symbol %in% genes_to_label ~ Symbol,
                              TRUE ~ ""))

labelled_points <- df_s5f %>% 
  filter(Symbol %in% genes_to_label)

ggplot(df_s5f, aes(x = logFC_shrunk, y=`-LOG10(FDR)`)) + 
  geom_point(aes(colour = significant),
             size = 0.5) + 
  geom_point(data = labelled_points,
             size = 0.5,
             colour = "black") +
  geom_point(data = labelled_points,
             size = 0.2, 
             colour = "#dd7676") + 
  geom_text(aes(label = to_label), size = xy_label_size*0.35, vjust = 1.5) +
  labs(x = "Enrichment in HAP1 versus HMEC (logFC)") +
  scale_colour_manual(values = significant_colours) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = 'None'
  ) +
  theme(text = element_text(size = xy_label_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white'))


###### S5g ###### 
df_figs5g <- ST8 %>% 
  group_by(Collection) %>%
  mutate(
    ID = fct_reorder(ID, NES)
  ) %>%
  ungroup()

ggplot(df_figs5g,
       aes(
         x = NES,
         y = ID,
         colour = p.adjust,
         size = setSize
       )) +
  geom_point(aes(fill = p.adjust),
             shape = 21,
             colour = "black",
             stroke = 0.3) +
  facet_wrap(~ Collection, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(
               GOBP = "Gene ontology biological process",
               Reactome = "Reactome"
             ))) +
  scale_fill_gradient(
    low = "#de6764",
    high = "#3b7eb8",
    name = "FDR"
  ) +
  scale_size_continuous(name = "Gene set size") +
  scale_y_discrete(
    labels = function(x) {
      x %>%
        gsub("^(REACTOME_|GOBP_)", "", .) %>%
        gsub("_", " ", .) %>%
        str_wrap(width = 30)
    }
  ) +
  labs(
    x = "Normalised enrichment score (NES)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = xy_label_size),
    text = element_text(size = xy_label_size),
    axis.text.y = element_text(size = xy_label_size, colour = "black"),
    axis.text.x = element_text(size = xy_label_size, colour = "black"),
    strip.background = element_blank(),
    legend.text = element_text(size = xy_label_size),
    strip.text.x = element_text(face = "bold")) 


###### S6a ######
ggplot(HMEC_variants_unfiltered,
       aes(y = r1_function_score,
           x = r2_function_score)) + 
  geom_point(aes(colour = Consequence),
             size = 0.4) + 
  coord_cartesian(ylim = c(-5, 2.6),
                  xlim = c(-5, 2.6)) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None"
  ) + 
  ylab('HMEC rep. 1 function score') + 
  xlab('HMEC rep. 2 function score') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) + 
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2)) +
  geom_abline(slope = 1, lty = 2) +
  scale_color_manual(values = consequence_colours)

cor(HMEC_variants_unfiltered$r1_function_score, HMEC_variants_unfiltered$r2_function_score, use = "pairwise.complete.obs", method = "pearson")


###### S6b ######
ggplot(ST3,
       aes(y = r1_function_score_ut,
           x = r2_function_score_ut)) + 
  geom_point(aes(colour = Consequence),
             size = 0.4) + 
  coord_cartesian(ylim = c(-5, 2.6),
                  xlim = c(-5, 2.6)) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None"
  ) + 
  ylab('HMEC rep. 1 function score') + 
  xlab('HMEC rep. 2 function score') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) + 
  scale_y_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2)) +
  geom_abline(slope = 1, lty = 2) +
  scale_color_manual(values = consequence_colours)

cor(ST3$r1_function_score_ut, ST3$r2_function_score_ut, use = "pairwise.complete.obs", method = "pearson")


###### S6c ######
df_figs6c <- ST3 %>% 
  filter(Consequence == "Synonymous" | Consequence == "Nonsense")

ggplot(data = df_figs6c,
       aes(x = Consequence,
           y = final_function_score_ut)) +
  geom_jitter(aes(colour = Consequence),
              size = 0.5) +
  ylab('Function score') +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "#000000") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) + 
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = 'None'
  ) +
  coord_cartesian(ylim = c(-3.2, 2.1)) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Nonsense", "Synonymous")), 
                     label = "p.signif",  
                     bracket.size = 0.5,  
                     tip.length = 0.02,  
                     label.y = 1.5) + 
  scale_y_continuous(breaks = c(-3, -2, -1, 0, 1)) +
  scale_colour_manual(values = consequence_colours)

wilcox.test(df_figs6c$final_function_score_ut[which(df_figs6c$Consequence == "Synonymous")], df_figs6c$final_function_score_ut[which(df_figs6c$Consequence == "Nonsense")])$p.value


###### S6e ######
df_figs6e <- ST3
df_figs6e$HMEC_class_ut <- factor(df_figs6e$HMEC_class_ut, levels = c("Neutral", "Depleted"))

ggplot(data = df_figs6e,
       aes(x = final_function_score_ut,
           fill = HMEC_class_ut)) +
  geom_histogram(alpha = 1, binwidth = 0.1) +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.position = 'None',
        legend.key = element_blank()
  ) +
  coord_cartesian(xlim = c(-4, 2)) +
  ylab('SNVs') +
  xlab('Function score') +
  theme(panel.background = element_rect(fill = 'white',
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  scale_fill_manual(values = HMEC_func_class_colours) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2))


##### Figure S7 #####
HMEC_sge_regions <- unique(as.character(ST3$SGE_region))

df_figs7 <- ST3
df_figs7$rep_correlation <- 0

calc_cor_per_region_HMEC <- function(region) {
  df <- ST3 %>%
    filter(SGE_region == region) %>% 
    filter(!is.na(r1_function_score_ut) & !is.na(r2_function_score_ut))
  correlation <- cor(df$r1_function_score_ut, df$r2_function_score_ut, use = "pairwise.complete.obs", method = "pearson")
  return(correlation)
}

for (reg in HMEC_sge_regions) {
  output <- calc_cor_per_region_HMEC(reg)
  print(output)
  df_figs7$rep_correlation[which(df_figs7$SGE_region == reg)] <- output
}

df_figs7 <- df_figs7 %>% 
  mutate(SGE_region = factor(SGE_region, levels = c("promoter", "exon 3", "exon 5", "exon 6", "exon 10 (5')", "exon 10 (3')", "exon 11", "intron 11", "exon 12 (5')", "exon 17"))) %>% 
  mutate(label = sprintf("italic(r) == '%.2f'", rep_correlation))

ggplot(df_figs7,
       aes(y = r1_function_score_ut, 
           x = r2_function_score_ut)) +
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  ylab('Rep. 1 function score') +
  xlab('Rep. 2 function score') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  scale_color_manual(values = consequence_colours) +
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1)) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1)) +
  geom_abline(slope = 1, lty = 2) +
  theme(strip.text = element_text(face = "bold",
                                  size = 7)) +
  facet_wrap(~SGE_region, ncol = 4, scales = "fixed", labeller = labeller(SGE_region = custom_labels)) +
  geom_text(
    aes(label = label),
    x = -4, y = 1, 
    parse = TRUE, 
    hjust = -0.1, vjust = 1.1,
    size = 7*0.35
  )


##### Figure S8 #####
df_figs8 <- ST3 %>% 
  dplyr::mutate(region = factor(SGE_region, levels = c("promoter", "exon 3", "exon 5", "exon 6", "exon 10 (5')", "exon 10 (3')", "exon 11", "intron 11", "exon 12 (5')", "exon 17"))) 

custom_labels <- c(
  "exon 5" = "Exon 5",
  "intron 11" = "Intron 11",
  "exon 10 (5')" = "Exon 10 (5')",
  "exon 10 (3')" = "Exon 10 (3')",
  "exon 10 (mid 1)" = "Exon 10 (mid 1)",
  "exon 10 (mid 2)" = "Exon 10 (mid 2)",
  "exon 11" = "Exon 11",
  "exon 12 (5')" = "Exon 12 (5')",
  "exon 12 (3')" = "Exon 12 (3')",
  "exon 17" = "Exon 17",
  "exon 3" = "Exon 3",
  "exon 6" = "Exon 6 (5')",
  "promoter" = "Promoter"
)

ggplot(df_figs8,
       aes(y = final_function_score_ut,
           x = hg38)) +
  geom_point(aes(colour = Consequence),
             size = 0.7) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  ylab('Function score') +
  xlab('Chr. 17 position (hg38)') +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  scale_color_manual(values = consequence_colours) +
  geom_hline(aes(yintercept = HMEC_threshold_ut), lty = 2) +
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0, 1)) +
  theme(strip.text = element_text(face = "bold",
                                  size = 7)) +
  scale_x_reverse(breaks = function(x) round(seq(min(x)+30, max(x)-30, length.out = 2), -1)) + # Fixed 3 breaks
  facet_wrap(~region, ncol = 4, scales = "free_x", labeller = labeller(region = custom_labels)) # 3x4 layout (adjust `ncol` for different layouts)


##### Figure S9 #####
###### S9a ######
df_figs9a <- ST3 %>% 
  filter(HAP1_SGE_region == "exon 12 (5')")

ggplot(data = df_figs9a,
       aes(y = final_function_score_olaparib,
           x = HAP1_function_score_mean)) + 
  geom_point(aes(colour = Consequence),
             size = 0.7) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None"
  ) + 
  labs(title = "Exon 12 (5')") +
  ylab('HMEC function score (+ olaparib)') + 
  xlab('HAP1 function score') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) + 
  geom_hline(yintercept = df_figs9a$HMEC_threshold_olaparib, lty = 2) +
  geom_vline(xintercept = HAP1_thresholds$threshold_upper[which(HAP1_thresholds$SGE_region == "exon 12 (5')")], lty = 2) +
  geom_vline(xintercept = threshold_lower, lty = 2, colour = "gray70") +
  scale_color_manual(values = consequence_colours)

###### S9c ######
PALB2_essentiality$frame <- factor(PALB2_essentiality$frame, levels = c("In-frame", "Frameshift"))

ggplot(PALB2_essentiality,
      aes(x = frame, y = log2_ratio_normalised)) + 
  geom_jitter(aes(x = frame,
                  y = log2_ratio_normalised, 
                  colour = frame),
              size = 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.title.x = element_blank(),
        legend.position = "none"
  ) +
  coord_cartesian(ylim = c(-3, 1.4)) +
  ylab("Indel function score") +
  ggtitle("PALB2") +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
        strip.background = element_blank()) + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Frameshift", "In-frame")), 
                     label = "p.signif",
                     tip.length = 0.02, 
                     label.y = 1.1,
                     size = 0.25,
                     textsize = 5*0.35) +
  scale_color_manual(values = frame_colours)

wilcox.test(PALB2_essentiality$log2_ratio_normalised[which(PALB2_essentiality$frame == "Frameshift")], PALB2_essentiality$log2_ratio_normalised[which(PALB2_essentiality$frame == "In-frame")])


###### S9d ######
RNF168_essentiality$frame <- factor(RNF168_essentiality$frame, levels = c("In-frame", "Frameshift"))

ggplot(RNF168_essentiality,
       aes(x = frame, y = log2_ratio_normalised)) + 
  geom_jitter(aes(x = frame,
                  y = log2_ratio_normalised, 
                  colour = frame),
              size = 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.title.x = element_blank(),
        legend.position = "none"
  ) +
  coord_cartesian(ylim = c(-3, 1.4)) +
  ylab("Indel function score") +
  ggtitle("RNF168") +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
        strip.background = element_blank()) + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Frameshift", "In-frame")), 
                     label = "p.signif",
                     bracket.size = 0.5,  
                     tip.length = 0.02,   
                     label.y = 1.1,
                     size = 0.25,
                     textsize = 5*0.35) +
  scale_color_manual(values = frame_colours)

wilcox.test(RNF168_essentiality$log2_ratio_normalised[which(RNF168_essentiality$frame == "Frameshift")], RNF168_essentiality$log2_ratio_normalised[which(RNF168_essentiality$frame == "In-frame")])


##### Figure S10 #####
###### S10a ######
df_figs10a <- ST1 %>% 
  filter(final_score == TRUE) %>% 
  filter(ClinVar_simple != "Absent") %>%
  select(cHGVS, Consequence, SGE_region, ClinVar_simple, ClinVar_review_status, final_function_score) %>% 
  filter(ClinVar_review_status != "no assertion criteria provided" & ClinVar_review_status != "no classification provided" & !is.na(ClinVar_review_status)) %>% 
  mutate(classification = case_when(grepl(pattern = "enign", x = ClinVar_simple) ~ 0,
                                    grepl(pattern = "athogenic", x = ClinVar_simple) ~ 1)) %>% 
  filter(!is.na(classification))

figs10a <- ggplot(df_figs10a,
                aes(d = classification,
                    m = -final_function_score)) +
  geom_roc(n.cuts=0,cutoffs.at = NULL) +
  style_roc(xlab = "1 - Specificity",
            major.breaks = c(0, 0.25, 0.5, 0.75, 1),
            minor.breaks = NULL) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank())

figs10a.annotated <- figs10a + annotate("text", x = .5, y = .5, label = paste("AUC =", round(calc_auc(figs10a)$AUC, 3)), size = 6 * 0.35)
figs10a.annotated

###### S10b ######
df_fig6ab <- ST3 %>% 
  filter(ClinVar_simple != "Absent") %>%
  select(cHGVS, Consequence, SGE_region, ClinVar_simple, ClinVar_review_status, final_function_score_ut, HAP1_function_score_mean) %>% 
  filter(ClinVar_review_status != "no assertion criteria provided" & ClinVar_review_status != "no classification provided" & !is.na(ClinVar_review_status)) %>% 
  mutate(classification = case_when(grepl(pattern = "enign", x = ClinVar_simple) ~ 0,
                                    grepl(pattern = "athogenic", x = ClinVar_simple) ~ 1)) %>% 
  filter(!is.na(classification))

df_figs10bc <- df_fig6ab %>% 
  filter(ClinVar_review_status != "criteria provided, single submitter") %>% 
  filter(!is.na(classification))

figs10b <- ggplot(df_figs10bc,
                  aes(d = classification,
                      m = -HAP1_function_score_mean)) +
  geom_roc(n.cuts=0,cutoffs.at = NULL) +
  style_roc(xlab = "1 - Specificity",
            major.breaks = c(0, 0.25, 0.5, 0.75, 1),
            minor.breaks = NULL) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank())

figs10b.annotated <- figs10b + annotate("text", x = .5, y = .5, label = paste("AUC =", round(calc_auc(figs10b)$AUC, 3)), size = 6 * 0.35)
figs10b.annotated


###### S10c ######
figs10c <- ggplot(df_figs10bc,
                  aes(d = classification,
                      m = -final_function_score_ut)) +
  geom_roc(n.cuts=0,cutoffs.at = NULL) +
  style_roc(xlab = "1 - Specificity",
            major.breaks = c(0, 0.25, 0.5, 0.75, 1),
            minor.breaks = NULL) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank())

figs10c.annotated <- figs10c + annotate("text", x = .5, y = .5, label = paste("AUC =", round(calc_auc(figs10c)$AUC, 3)), size = 6 * 0.35)
figs10c.annotated

###### S10d ######
df_figs10d <- ST1 %>% 
  filter(SGE_region == "exon 17") 

ggplot(df_figs10d,
       aes(y = final_function_score, 
           x = hg38)) + 
  geom_point(aes(colour = Consequence),
             size = 0.5) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('HAP1 function score') + 
  xlab('Chr. 17 position (hg38)') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  geom_hline(yintercept = df_figs10d$threshold_upper, lty = 2) +
  geom_hline(yintercept = threshold_lower, lty = 2, colour = "grey") +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  scale_color_manual(values = consequence_colours) + 
  scale_x_reverse(breaks = c(43063960, 43063920, 43063880))


###### S10e ######
df_figs10e <- ST3 %>% 
  filter(SGE_region == "exon 17") 

ggplot(data = df_figs10e,
       aes(y = final_function_score_ut, 
           x = hg38)) + 
  geom_point(aes(colour = Consequence),
             size = 0.5) + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") + 
  ylab('HMEC function score') + 
  xlab('Chr. 17 position (hg38)') + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  geom_hline(yintercept = df_figs10e$HMEC_threshold_ut, lty = 2) +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size, colour = "#000000"),
        axis.text.x = element_text(size = xy_label_size, colour = "#000000"),
        strip.background = element_blank()) +
  scale_color_manual(values = consequence_colours) + 
  scale_x_reverse(breaks = c(43063960, 43063920, 43063880))


###### S10f ######
df_figs10f <- ST1 %>%
  filter(in_truthset == TRUE)

ggplot(data = df_figs10f) +
  geom_histogram(aes(x = final_function_score,
                     fill = ClinVar_simple)) +
  xlab('HAP1 function score') +
  ylab('') +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size),
        axis.text.x = element_text(size = xy_label_size),
        strip.background = element_blank()) +
  scale_fill_manual(values = clinvar_colours)


###### S10g ######
df_figs10g <- ST3 %>%
  filter(in_HMEC_truthset == TRUE)

ggplot(data = df_figs10g) +
  geom_histogram(aes(x = final_function_score_ut,
                     fill = ClinVar_simple)) +
  xlab('HMEC function score') +
  ylab('') +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  theme(text = element_text(size = text_size),
        axis.text.y = element_text(size = xy_label_size),
        axis.text.x = element_text(size = xy_label_size),
        strip.background = element_blank()) +
  scale_fill_manual(values = clinvar_colours)


###### S10h ######
df_figs10h <- ST1 %>% 
  filter(!is.na(evidence_code_original))

ggplot(data = df_figs10h) +
  geom_histogram(aes(x = final_function_score,
                     fill = evidence_code_original),
                 binwidth = 0.05) +
  xlab('HAP1 function score') +
  ylab('') +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_fill_manual(values = evidence_code_colours)

###### S10i ######
df_figs10i <- ST3

ggplot(data = df_figs10i) + 
  geom_histogram(aes(x = final_function_score_ut,
                     fill = HMEC_evidence_code_original),
                 binwidth = 0.05) +
  xlab('HMEC function score') +
  ylab('') +
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'white')) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "None") +
  theme(text = element_text(size = text_size, colour = "black"),
        axis.text.y = element_text(size = xy_label_size, colour = "black"),
        axis.text.x = element_text(size = xy_label_size, colour = "black"),
        strip.background = element_blank()) +
  scale_fill_manual(values = evidence_code_colours)
