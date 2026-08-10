# ################################################################################################################################################## #
# ## MGB Young Bladder Cancer Patients ## ---------------------------------------------------------------------------------------------------------- #
# ## Master Script ## ------------------------------------------------------------------------------------------------------------------------------ #
# ################################################################################################################################################## #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
###  Last Edit 2026-07


rm(list = ls())  # clear memory
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Global set up --------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

# File dates --------------------------------------------------------------------------------------------------------------------------------------- #
#datafiles
date.data_rpdr <- "20230607" #raw data from RPDR
date.data_clin <- "20230922" #clinical data, manual review from EPIC
#analyzed
date.cleaned  <- "20231020"  # file generated after cleaning data
date.analyzed <- "20260703"  # most recent analysis (today if updating)

# File locations ----------------------------------------------------------------------------------------------------------------------------------- #
local.data.flag <- 1 #use 1 for local, 0 for remote drive D
# local file paths removed for publishing purposes
if(local.data.flag){
  drive.datasuppository <- "C:~"
  drive.dataMGBonco     <- "C:~"
  drive.dataMSKgene     <- "C:~"
  } else {
    drive.datasuppository <- "D:~"
    drive.dataMGBonco     <- "D:~"
    drive.dataMSKgene     <- ""
  }
folder.RPDR.lastpull  <- paste(drive.datasuppository, "20230707_RPDR_RC/", sep="") 
drive.workspace       <- "C:~"
folder.results        <- "Results/"
folder.figures        <- "Figures/"

# Load packages ------------------------------------------------------------------------------------------------------------------------------------ #
# list.packages <- c("maftools","tidyverse","readxl","arsenal","tableone","ggpubr","ggsci","survival","survminer")
# install.packages(list.packages)
# library(BiocManager)
library(maftools)
library(tidyverse)
library(readxl)
library(arsenal)
library(tableone)
library(ggpubr)
library(ggsci)
library(survival)
library(survminer)
library(stringr)
library(readr)

# load data MGB Institutional Data ----------------------------------------------------------------------------------------------------------------- #

#initial available oncopanels used
df.gene <- read.csv( file = paste(drive.dataMGBonco, "20230712_DFCI_BladderOncopanel_cleaned_gene.csv", sep=""), stringsAsFactors = FALSE) %>%
  mutate(TERT   = if_else(HARMONIZED_HUGO_GENE_NAME == "TERT", 1,0,0),
         TP53   = if_else(HARMONIZED_HUGO_GENE_NAME == "TP53", 1,0,0),
         FGFR3  = if_else(HARMONIZED_HUGO_GENE_NAME == "FGFR3", 1,0,0),
         KDM6A  = if_else(HARMONIZED_HUGO_GENE_NAME == "KDM6A", 1,0,0),
         PIK3CA = if_else(HARMONIZED_HUGO_GENE_NAME == "PIK3CA", 1,0,0),
         CDKN2A = if_else(HARMONIZED_HUGO_GENE_NAME == "CDKN2A", 1,0,0),
         ARID1A = if_else(HARMONIZED_HUGO_GENE_NAME == "ARID1A", 1,0,0),
         CDKN2B = if_else(HARMONIZED_HUGO_GENE_NAME == "CDKN2B", 1,0,0),
         KMT2D  = if_else(HARMONIZED_HUGO_GENE_NAME == "KMT2D", 1,0,0),
         CDKN1A = if_else(HARMONIZED_HUGO_GENE_NAME == "CDKN1A", 1,0,0),
         RB1    = if_else(HARMONIZED_HUGO_GENE_NAME == "RB1", 1,0,0),
         CCND1  = if_else(HARMONIZED_HUGO_GENE_NAME == "CCND1", 1,0,0),
         ERBB2  = if_else(HARMONIZED_HUGO_GENE_NAME == "ERBB2", 1,0,0),
         STAG2  = if_else(HARMONIZED_HUGO_GENE_NAME == "STAG2", 1,0,0),
         EP300  = if_else(HARMONIZED_HUGO_GENE_NAME == "EP300", 1,0,0),
         
         ARID2  = if_else(HARMONIZED_HUGO_GENE_NAME == "ARID2", 1,0,0),
         ATM    = if_else(HARMONIZED_HUGO_GENE_NAME == "ATM", 1,0,0),
         ATRX   = if_else(HARMONIZED_HUGO_GENE_NAME == "ATRX", 1,0,0),
         BRCA1  = if_else(HARMONIZED_HUGO_GENE_NAME == "BRCA1", 1,0,0),
         CHEK2  = if_else(HARMONIZED_HUGO_GENE_NAME == "CHEK2", 1,0,0),
         CIITA  = if_else(HARMONIZED_HUGO_GENE_NAME == "CIITA", 1,0,0),
         CREBBP = if_else(HARMONIZED_HUGO_GENE_NAME == "CREBBP", 1,0,0),
         CUX1   = if_else(HARMONIZED_HUGO_GENE_NAME == "CUX1", 1,0,0),
         NOTCH1 = if_else(HARMONIZED_HUGO_GENE_NAME == "NOTCH1", 1,0,0),
         NOTCH2 = if_else(HARMONIZED_HUGO_GENE_NAME == "NOTCH2", 1,0,0),
         SLX4   = if_else(HARMONIZED_HUGO_GENE_NAME == "SLX4", 1,0,0),
         APC    = if_else(HARMONIZED_HUGO_GENE_NAME == "APC", 1,0,0),
         ERCC6  = if_else(HARMONIZED_HUGO_GENE_NAME == "ERCC6", 1,0,0),
         FANCM  = if_else(HARMONIZED_HUGO_GENE_NAME == "FANCM", 1,0,0),
         FH     = if_else(HARMONIZED_HUGO_GENE_NAME == "FH", 1,0,0),
         KMT2A  = if_else(HARMONIZED_HUGO_GENE_NAME == "KNT2A", 1,0,0),
         NF2    = if_else(HARMONIZED_HUGO_GENE_NAME == "NF2", 1,0,0),
         PRKDC  = if_else(HARMONIZED_HUGO_GENE_NAME == "PRKDC", 1,0,0))
# nrow(df.gene); n_distinct(df.gene$DFCI_MRN); n_distinct(df.gene$UNIQUE_SAMPLE_ID)
# 13833 row; 810 MRNs; 832 samples 

#below line will re-run data cleaning script and generate df.clin and df.onc
#generates table of excluded patients
flag.save.tabs <- F
source("BladderYoung_02C_datacleaning_202607.R")
rm(flag.save.tabs)
# create dataset for clinical outcome analysis
df.clin_analysis <- df.clin %>% filter(Exclude_clin == 0)

df.MGB_bladder <- df.clin_analysis %>%
  mutate(sample = EPIC_PMRN,
         NAC = if_else(NAC==1, 1,0,0)) %>%
  mutate(Age_Dx.bin = factor(case_when(Age_Dx < 45 ~ "1.<45",
                                       Age_Dx < 50 ~ "2.45-49",
                                       Age_Dx < 55 ~ "3.50-54")),
         Gender = factor(if_else(Gender == "Male", "1.Male","0.Female"), levels = c("0.Female","1.Male")),
         Race = factor(if_else(Race == "1.White", "1.White","0.Non-White"), levels = c("0.Non-White","1.White")),
         Tobacco = factor(if_else(Tobacco == 0, "0.No","1.Yes",NA)),
         Stage.NMIBC = factor(case_when(RC.T_stage == "0.T0"  ~ "0.<T2",
                                        RC.T_stage == "1.Ta"  ~ "0.<T2",
                                        RC.T_stage == "2.Tis" ~ "0.<T2",
                                        RC.T_stage == "3.T1"  ~ "0.<T2",
                                        RC.T_stage == "4.T2"  ~ "1.>=T2",
                                        RC.T_stage == "5.T3"  ~ "1.>=T2",
                                        RC.T_stage == "6.T4"  ~ "1.>=T2")),
         Systemic_Treatment = factor(case_when((Adjuvant_Immuno == 1) & (NAC == 1)              ~ "3.Both",
                                               (Adjuvant_Immuno == 1) & (Adjuvant_Chemo != "0") ~ "3.Both",
                                               (Adjuvant_Immuno == 1)                           ~ "2.Immunotherapy",
                                               (NAC == 1)                                       ~ "1.Chemo",
                                               (Adjuvant_Chemo != "0")                          ~ "1.Chemo",
                                               .default = "0.None")),
         Intravesical_Treatment = factor(case_when((BCG == 1) & (intravesical_chemo == 1) ~ "3.Both",
                                                   (BCG == 1)                             ~ "1.BCG",
                                                   (intravesical_chemo == 1) ~ "2.Chemo",
                                                   .default = "0.None")))
MGB.colnames <- c("sample", "EPIC_PMRN", "DFCI_MRN",
                  "Age_Dx.bin", "Age_Dx","t_RC_FU", "BMI", "Gender", "Race", "Tobacco",
                  "Stage.NMIBC", "Intravesical_Treatment", "Systemic_Treatment", "NAC", "NAC_response",
                  "RC.Diversion_Conduit", "RC.Robotic", "RC.LOS", "RC.Op_Time",
                  "Comp30d.Readmission", "Comp30d.Ileus", "Comp30d.transfusion", "Comp30d.DVT_PE", "Comp30d.RTOR", "RC.T_stage",
                  "Adjuvant_Immuno", "Adjuvant_Chemo", "BCG", "intravesical_chemo"); length(MGB.colnames)
df.MGB_bladder <- df.MGB_bladder[,MGB.colnames]

#obtain gene data
df.MGB_onc <- left_join(df.MGB_bladder, df.gene, by = "DFCI_MRN"); n_distinct(df.MGB_onc$EPIC_PMRN)
df.MGB_onc$Oncopanel_match <- if_else(!is.na(df.MGB_onc$UNIQUE_SAMPLE_ID), 1,0,0)
# #names of genes
# colnames(df.MGB_onc[,c(119:121,156:(length(df.onc)))])
temp.onc <- df.MGB_onc %>%
  select(EPIC_PMRN, Oncopanel_match,
         # top 20
         TP53, KMT2D, KDM6A, RB1, ARID1A,
         ARID2, ATM, ATRX, BRCA1, CHEK2,
         CIITA, CREBBP, CUX1, NOTCH1, NOTCH2,
         PIK3CA, SLX4, APC, CDKN1A, ERCC6,
         FANCM, FH, KMT2A, NF2, PRKDC,
         # also in MSK top 15
         TERT, FGFR3, CDKN2A, CDKN2B, CCND1, ERBB2, STAG2, EP300) %>%
  group_by(EPIC_PMRN) %>%
  summarize(Oncopanel_match = max(Oncopanel_match),
            TP53 = max(TP53), KMT2D = max(KMT2D), KDM6A = max(KDM6A), RB1 = max(RB1), ARID1A = max(ARID1A),
            
            ARID2 = max(ARID2), ATM = max(ATM), ATRX = max(ATRX), BRCA1 = max(BRCA1), CHEK2 = max(CHEK2),
            CIITA = max(CIITA), CREBBP = max(CREBBP), CUX1 = max(CUX1), NOTCH1 = max(NOTCH1), NOTCH2 = max(NOTCH2),
            PIK3CA = max(PIK3CA), SLX4 = max(SLX4), APC = max(APC), CDKN1A = max(CDKN1A), ERCC6 = max(ERCC6),
            FANCM = max(FANCM), FH = max(FH), KMT2A = max(KMT2A), NF2 = max(NF2), PRKDC = max(PRKDC),
            # also in MSK top 15
            TERT = max(TERT), FGFR3 = max(FGFR3), CDKN2A = max(CDKN2A), CDKN2B = max(CDKN2B),
            CCND1 = max(CCND1), ERBB2 = max(ERBB2), STAG2 = max(STAG2), EP300 = max(EP300)) %>%
  distinct()

df.MGB_bladder <- left_join(df.MGB_bladder, temp.onc, by = "EPIC_PMRN")
rm(temp.onc)

# load data MSK IMPACT Data ------------------------------------------------------------------------------------------------------------------------ #

#load MSK IMPACT data
df.MSK_bladder <- read_excel(paste(drive.dataMSKgene, "IMPACTBladderStats51822RevisionSupplement.xlsx", sep=""), sheet = "Sheet1", trim_ws = TRUE) %>%
  rename("sample"                 = "IMPACT Sample ID",
         "Smoking_status"         = "Smoking status",
         "Age_Dx"                 = "Age at Dx",
         "Specimen_Type"          = "Specimen Type",
         "Intravesical_Treatment" = "Intravesical Treatment",
         "Systemic_Treatment"     = "Systemic Treatment") %>%
  mutate(Gender = factor(if_else(Sex == "0", "1.Male","0.Female"), levels = c("0.Female","1.Male")),
         Stage = as.numeric(Stage),
         Age_Dx = round(as.numeric(if_else(Age_Dx == ".", NA, Age_Dx)),1),
         Tobacco = factor(if_else(Smoking_status == 0, "0.No","1.Yes",NA)),
         Stage.NMIBC = factor(case_when(Stage == 1 ~ "0.<T2",
                                        Stage == 2 ~ "0.<T2",
                                        Stage == 3 ~ "0.<T2",
                                        Stage == 4 ~ "1.>=T2",
                                        Stage == 5 ~ "1.>=T2")),
         Systemic.Therapy = factor(case_when(Systemic_Treatment == 0 ~ "0.No",
                                             Systemic_Treatment == 1 ~ "1.Yes",
                                             Systemic_Treatment == 2 ~ "1.Yes",
                                             Systemic_Treatment == 3 ~ "1.Yes"))) %>%
  mutate(flag_age55 = case_when(Age_Dx < 0 ~ 2,
                                is.na(Age_Dx) ~ 3,
                                Age_Dx <55 ~ 1,
                                Age_Dx >= 55 ~ 0),
         Age_Dx.bin = factor(case_when(Age_Dx < 45 ~ "1.<45",
                                       Age_Dx < 50 ~ "2.45-49",
                                       Age_Dx < 55 ~ "3.50-54")),
         Specimen_Type = if_else(Specimen_Type == "Partial Cx", "Partial cystectomy", Specimen_Type),
         StageSpecimen = factor(case_when(StageSpecimen == 1 ~ "1.Low.Grade",
                                          StageSpecimen == 2 ~ "2.High.Grade.Non-Invasive",
                                          StageSpecimen == 3 ~ "3.High.Grade.Invasive",
                                          StageSpecimen == 4 ~ "4.Metastatic")),
         Stage = factor(case_when(Stage == 1 ~ "1.Low.Grade",
                                  Stage == 2 ~ "2.High.Grade",
                                  Stage == 3 ~ "3.CIS",
                                  Stage == 4 ~ "4.Muscle.Invasive",
                                  Stage == 5 ~ "5.Metastatic")),
         Histology = factor(case_when(Histology == 0 ~ "0.Urothelial.Carcinoma",
                                      Histology == 1 ~ "1.Squamous.Cell",
                                      Histology == 2 ~ "2.Plasmacytoid",
                                      Histology == 3 ~ "3.Small.Cell.Neuroendocrine",
                                      Histology == 4 ~ "4.Micropapillary",
                                      Histology == 5 ~ "5.Glandular",
                                      Histology == 6 ~ "6.Other")),
         Intravesical_Treatment = factor(case_when(Intravesical_Treatment == 0 ~ "0.None",
                                                   Intravesical_Treatment == 1 ~ "1.BCG",
                                                   Intravesical_Treatment == 2 ~ "2.Chemo",
                                                   Intravesical_Treatment == 3 ~ "3.Both")),
         Systemic_Treatment = factor(case_when(Systemic_Treatment == 0 ~ "0.None",
                                               Systemic_Treatment == 1 ~ "1.Chemo",
                                               Systemic_Treatment == 2 ~ "2.Immunotherapy",
                                               Systemic_Treatment == 3 ~ "3.Both")))
MSK.colnames <- colnames(df.MSK_bladder)[c(1,73,74,6,69,70,71,10,11,3,7,2,9,8,4,5,72,12:68)]; length(MSK.colnames); print(MSK.colnames)

summary(freqlist(~flag_age55, data=df.MSK_bladder, addNA = T))
# remove those w flag_age55 == 2 (n=1) and == 3 (n=41); 1271 total remain
df.MSK_bladder <- df.MSK_bladder[,MSK.colnames] %>%
  filter(flag_age55 <= 1)

# add filter for RC and TUR patients
df.MSK_bladder <- df.MSK_bladder %>%
  mutate(RC_TUR_filter = factor(case_when(Specimen_Type == "Cystectomy" ~ "0.RC",
                                          Specimen_Type == "TUR"        ~ "1.TUR")))

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Analysis -------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

# MSK gene mutations ------------------------------------------------------------------------------------------------------------------------------- #

flag.save.figs <- F
flag.save.tabs <- F
source("BladderYoung_03B_genemutationanalysis_202607.R")
rm(flag.save.figs, flag.save.tabs)

# Data characterization ---------------------------------------------------------------------------------------------------------------------------- #

#legacy scripts included below but use this one:
flag.save.figs <- F
flag.save.tabs <- F
source("BladderYoung_03_analysis_characterization_202607.R")
rm(flag.save.figs, flag.save.tabs)

# Oncopanel ---------------------------------------------------------------------------------------------------------------------------------------- #
# oncoplot
df.clin <- df.clin_analysis  #uses same cohort as clinical analysis

flag.save.figs <- F
flag.save.tabs <- F
source("BladderYoung_04_analysis_oncopanel_202607.R") 

# Gene Bar Plot ------------------------------------------------------------------------------------------------------------------------------------ #

df.clin <- df.clin_analysis  #uses same cohort as clinical analysis

flag.save.figs <- F
source("BladderYoung_05_analysis_genebarplot_202607.R")

# MSK-IMPACT comutation plot ----------------------------------------------------------------------------------------------------------------------- #

flag.save.figs <- F
flag.save.tabs <- F
source("BladderYoung_06_analysis_oncopanel_MSK_202607.R") 

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- End ------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
sessionInfo()