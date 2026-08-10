# ################################################################################################################################################## #
# ## MGB Young Bladder Cancer Patients ## ---------------------------------------------------------------------------------------------------------- #
# ## Data Cleaning ## ------------------------------------------------------------------------------------------------------------------------------ #
# ################################################################################################################################################## #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
###  Last Edit 2026-07


# must load packages and file names from main script first

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Load data ------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

#prelim clinical annotation file (from manual review)
df.clin <- read.csv(file = paste(drive.datasuppository, date.data_clin, "_MGB_BladderCa_YoungPts_clinical.csv", sep=""), stringsAsFactors = FALSE)
# nrow(df.clin); n_distinct(df.clin$EPIC_PMRN)  # 182 MRNs


# clean -------------------------------------------------------------------------------------------------------------------------------------------- #

# dates
df.clin <- df.clin %>% mutate_at(vars(c("Date_of_Birth", "Date_bladderdx", "Date_RC", "Date_Discharge_RC", "Date_Last_FU_Encounter", "Date_Of_Death",
                                        "Date_Recurrence_Local", "Date_Recurrence_Distant", "Date_Prog", "Date_TURBT", "Comp.Readmission",
                                        "Comp.RTOR", "Comp.transfusion", "Comp.DVT_PE", "Comp.NGT", "Comp.Ileus", "Comp.SBO",
                                        "Comp.hernia")),
                                 ~ as.Date(., format = "%m/%d/%Y"))
# clean variables
df.clin <- df.clin %>%
  # recalculate to ensure rigor (for t_* use months and restore decimal)
  mutate(t_Dx_RC  = as.numeric(difftime(Date_RC,                Date_bladderdx, units = "days")) /(365.25/12),
         t_Dx_FU  = as.numeric(difftime(Date_Last_FU_Encounter, Date_bladderdx, units = "days")) /(365.25/12),
         t_RC_FU  = as.numeric(difftime(Date_Last_FU_Encounter, Date_RC, units = "days")) /(365.25/12),
         t_OS     = as.numeric(difftime(Date_Last_FU_Encounter, Date_RC, units = "days")) /(365.25/12),  #same as t_RC_FU
         t_CSS    = as.numeric(difftime(Date_Last_FU_Encounter, Date_RC, units = "days")) /(365.25/12),  #same as t_RC_FU
         min30dFU = if_else(as.numeric(difftime(Date_Last_FU_Encounter, Date_RC, units = "days")) < 30, 0,1,0),  #for calculating clinical outcomes
         
         Event_OS  = if_else((Dead == 1), 1,0,0),
         Event_CSS = if_else((Dead == 1) &
                               (Death_from_Disease == "1" | Death_from_Disease == "mesenteric ischemia 2/2 SMA thrombus now sp SBR" ), 1,0,0),
         Event_PFS = if_else(!is.na(Date_Prog), 1,0,0),
         
         Age_RC    = trunc(as.numeric(difftime(Date_RC,        Date_of_Birth, units = "days")) /365.25),
         Age_Dx    = trunc(as.numeric(difftime(Date_bladderdx, Date_of_Birth, units = "days")) /365.25),
         Age_Death = trunc(as.numeric(difftime(Date_Of_Death,  Date_of_Birth, units = "days")) /365.25)) %>%
  mutate(t_PFS = if_else((Event_PFS == 0), t_OS, as.numeric(difftime(Date_Prog, Date_RC, units = "days")) /(365.25/12),NA)) %>%
  # other variables
  mutate(Age_RC.bin = factor(case_when(Age_RC <  45 ~ "1.<45",
                                       Age_RC <  50 ~ "2.45-49",
                                       Age_RC >= 45 ~ "3.50-54")),
         Race = factor(case_when(Race == "White"    ~ "1.White",
                                 Race == "Asian"    ~ "2.Asian",
                                 Race == "Black"    ~ "3.Black",
                                 Race == "Hispanic" ~ "4.Hispanic",
                                 Race == "Hispanic" ~ "4.Hispanic",
                                 Race %in% c("American Indian or Alaska Native","Other","Unknown/Missing") ~ "5.Other")),
         RC.Transfusion       = if_else(RC.PRBC >0, 1,0,0),  #flag if any transfusion
         RC.Op_Time           = RC.Op_Time/60,  #convert to hours
         Gender               = factor(Gender, levels = c("Male", "Female")),
         RC.Diversion_Type    = factor(RC.Diversion_Type),
         RC.Diversion_Conduit = if_else((RC.Diversion_Type == "1.Conduit") | (RC.Diversion_Type == "1.Conduit.aborted"), 1,0),
         RC.NED               = if_else((RC.pT == "T0") & (RC.pN == "0") & (RC.Histology=="negative"), 1,0),
         RC.T_stage = factor(case_when(RC.pT == "T0"  ~ "0.T0",
                                       RC.pT == "Ta"  ~ "1.Ta",
                                       RC.pT == "Tis" ~ "2.Tis",
                                       RC.pT == "T1"  ~ "3.T1",
                                       (RC.pT == "T2") | (RC.pT == "T2a") | (RC.pT == "T2b") ~ "4.T2",
                                       (RC.pT == "T3") | (RC.pT == "T3a") | (RC.pT == "T3b") ~ "5.T3",
                                       (RC.pT == "T4") | (RC.pT == "T4a") | (RC.pT == "T4b") ~ "6.T4")),
         RC.pT = factor(RC.pT),
         RC.pN = factor(RC.pN),
         RC.Variant    = if_else((RC.Variant.Histology    == "0"), 0,1,0),
         TURBT.Variant = if_else((TURBT.Variant.Histology == "0" | TURBT.Variant.Histology == ""), 0,1,0),
         TURBT.pT = factor(TURBT.pT),
         TURBT.T_stage = factor(case_when(TURBT.pT == "T0"  ~ "0.T0",
                                          TURBT.pT == "Ta"  ~ "1.Ta",
                                          TURBT.pT == "Tis" ~ "2.Tis",
                                          TURBT.pT == "T1"  ~ "3.T1",
                                          TURBT.pT == "T2"  ~ "4.T2",
                                          (TURBT.pT == "T3") | (TURBT.pT == "T3a") ~ "5.T3")),
         # place holders
         Exclude_clin = if_else((Exclude == 1) | (Exclude.Description %in% c("0.Aborted.RC.mets","0.Salvage.RC")), 1,0,0)) %>%
  # complications variables: readmission based on dc date, others based on RC date
  mutate( t_Comp.Readmission = as.numeric(Comp.Readmission - Date_Discharge_RC),
          Comp30d.Readmission = if_else( (Comp.Readmission - Date_Discharge_RC) <=30    ,1,0,0),
          Comp6mo.Readmission = if_else( (Comp.Readmission - Date_Discharge_RC) <=182.5 ,1,0,0),
          Comp1yr.Readmission = if_else( (Comp.Readmission - Date_Discharge_RC) <=365   ,1,0,0),
          
          t_Comp.RTOR = as.numeric(Comp.RTOR - Date_RC),
          Comp30d.RTOR = if_else( (Comp.RTOR - Date_RC) <=30    ,1,0,0),
          Comp6mo.RTOR = if_else( (Comp.RTOR - Date_RC) <=182.5 ,1,0,0),
          Comp1yr.RTOR = if_else( (Comp.RTOR - Date_RC) <=365   ,1,0,0),
          
          t_Comp.transfusion = as.numeric(Comp.transfusion - Date_RC),
          Comp30d.transfusion = if_else( (Comp.transfusion - Date_RC) <=30    ,1,0,0),
          Comp6mo.transfusion = if_else( (Comp.transfusion - Date_RC) <=182.5 ,1,0,0),
          Comp1yr.transfusion = if_else( (Comp.transfusion - Date_RC) <=365   ,1,0,0),
          
          t_Comp.DVT_PE = as.numeric(Comp.DVT_PE - Date_RC),
          Comp30d.DVT_PE = if_else( (Comp.DVT_PE - Date_RC) <=30    ,1,0,0),
          Comp6mo.DVT_PE = if_else( (Comp.DVT_PE - Date_RC) <=182.5 ,1,0,0),
          Comp1yr.DVT_PE = if_else( (Comp.DVT_PE - Date_RC) <=365   ,1,0,0),
          
          t_Comp.Ileus = as.numeric(Comp.Ileus - Date_RC),
          Comp30d.Ileus = if_else( (Comp.Ileus - Date_RC) <=30    ,1,0,0),
          Comp6mo.Ileus = if_else( (Comp.Ileus - Date_RC) <=182.5 ,1,0,0),
          Comp1yr.Ileus = if_else( (Comp.Ileus - Date_RC) <=365   ,1,0,0),
          
          t_Comp.SBO = as.numeric(Comp.SBO - Date_RC),
          Comp30d.SBO = if_else( (Comp.SBO - Date_RC) <=30    ,1,0,0),
          Comp6mo.SBO = if_else( (Comp.SBO - Date_RC) <=182.5 ,1,0,0),
          Comp1yr.SBO = if_else( (Comp.SBO - Date_RC) <=365   ,1,0,0),
          
          t_Comp.hernia = as.numeric(Comp.hernia - Date_RC),
          Comp30d.hernia = if_else( (Comp.hernia - Date_RC) <=30    ,1,0,0),
          Comp6mo.hernia = if_else( (Comp.hernia - Date_RC) <=182.5 ,1,0,0),
          Comp1yr.hernia = if_else( (Comp.hernia - Date_RC) <=365   ,1,0,0) ) %>%
  select(
    EPIC_PMRN, DFCI_MRN, DFCI, Oncopanel_match,
    Exclude_gene = Exclude, Exclude.Description, min30dFU, Exclude_clin,
    
    t_Dx_RC, t_Dx_FU, t_RC_FU, t_OS, t_CSS, Dead, Event_OS, Death_from_Disease, Event_CSS, t_PFS, Event_PFS,
    
    Age_RC, Age_RC.bin, Age_Dx, Age_Death, BMI, Gender, Race, Tobacco,
    
    RC.Diversion_Conduit, RC.Diversion_Type, RC.Robotic, RC.LOS, RC.Op_Time, RC.EBL, RC.Fluids, RC.Transfusion, RC.PRBC,
    TURBT.Num, BCG, intravesical_chemo, NAC, NAC_response, Adjuvant_Chemo, Adjuvant_Immuno, Adjuvant_XRT,
    
    RC.ACCESSION.NUM, RC.NED, RC.T_stage, RC.pT, RC.pN, RC.Histology, RC.Variant, RC.Variant.Histology,
    Date_TURBT, TURBT.ACCESSION.NUM,
    TURBT.T_stage, TURBT.pT, TURBT.Histology, TURBT.Grade, TURBT.Variant, TURBT.Variant.Histology, TURBT.pLVI, TURBT.MIBC, TURBT.CIS, TURBT.Size3cm,
    
    Comp.Readm.Reason,
    Comp.Readmission, t_Comp.Readmission, Comp30d.Readmission, Comp6mo.Readmission, Comp1yr.Readmission,
    Comp.RTOR,        t_Comp.RTOR,        Comp30d.RTOR,        Comp6mo.RTOR,        Comp1yr.RTOR,
    Comp.transfusion, t_Comp.transfusion, Comp30d.transfusion, Comp6mo.transfusion, Comp1yr.transfusion,
    Comp.DVT_PE,      t_Comp.DVT_PE,      Comp30d.DVT_PE,      Comp6mo.DVT_PE,      Comp1yr.DVT_PE,
    Comp.NGT,
    Comp.Ileus,       t_Comp.Ileus,       Comp30d.Ileus,       Comp6mo.Ileus,       Comp1yr.Ileus,
    Comp.SBO,         t_Comp.SBO,         Comp30d.SBO,         Comp6mo.SBO,         Comp1yr.SBO,
    Comp.hernia,      t_Comp.hernia,      Comp30d.hernia,      Comp6mo.hernia,      Comp1yr.hernia,
    
    Date_of_Birth, Date_bladderdx, Date_RC, Date_Discharge_RC, Date_Last_FU_Encounter, Date_Of_Death,
    Date_Recurrence_Local, Date_Recurrence_Distant, Date_Prog, Recurrence_Location, Recurrence_Tx_Type,
    State, Zip, Marital_status, Insurance, Notes)


# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Join Oncopanel -------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

df.onc <- left_join(df.clin, df.gene, by = "DFCI_MRN")
n_distinct(df.onc$EPIC_PMRN)

df.onc$Oncopanel_match <- if_else(!is.na(df.onc$UNIQUE_SAMPLE_ID), 1,0,0)
#names of genes
colnames(df.onc[,c(119:121,156:(length(df.onc)))])
#"ARID1A" "ATM"    "ERCC2"  "TP53"   "KMT2D"  "KDM6A"  "RB1"
temp.onc <- df.onc %>%
  select(EPIC_PMRN, Oncopanel_match, ARID1A, ATM, ERCC2, TP53, KMT2D, KDM6A, RB1) %>%
  group_by(EPIC_PMRN) %>%
  summarize(Oncopanel_match = max(Oncopanel_match),
            ARID1A = max(ARID1A),
            ATM    = max(ATM),
            ERCC2  = max(ERCC2),
            KMT2D  = max(KMT2D),
            KDM6A  = max(KDM6A),
            RB1    = max(RB1),
            TP53   = max(TP53)) %>%
  distinct()

df.clin <- left_join(df.clin, temp.onc, by = "EPIC_PMRN") %>%
  mutate(Oncopanel_match.x = Oncopanel_match.y, Oncopanel_match.y = NULL) %>%
  rename(Oncopanel_match = Oncopanel_match.x)
rm(temp.onc)


# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Exclusions ------------------------------------------------------------------------------------------------------------------------------------ #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

# remove those not bladder cancer
summary(freqlist(~Exclude_gene, data=df.clin, addNA = T))
summary(freqlist(~Exclude_gene+Oncopanel_match, data=df.clin, addNA = T))
summary(freqlist(~Exclude_gene+Exclude.Description, data=df.clin %>% filter(Exclude_gene==1), addNA = T))

df.clin <- df.clin %>% filter(Exclude_gene == 0) %>% mutate(Exclude_gene = NULL)  #remove those not bladder cancer, remove variable
df.onc  <- df.onc  %>% filter(Exclude_gene == 0) %>% mutate(Exclude_gene = NULL)  #remove those not bladder cancer, remove variable

# now for clinical exclusions (don't take out of dataframe)
summary(freqlist(~Exclude_clin, data=df.clin, addNA = T))
summary(freqlist(~Exclude_clin+Exclude.Description, data=df.clin, addNA = T))
summary(freqlist(~min30dFU, data=df.clin, addNA = T))
summary(freqlist(~Exclude_clin+min30dFU, data=df.clin, addNA = T))
n_distinct(df.clin %>% filter(min30dFU==1 & Exclude_clin==0) %>% select(EPIC_PMRN))

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- save ------------------------------------------------------------------------------------------------------------------------------------------ #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# export as csv
if(flag.save.tabs){
  print(paste("Tables updated ",date.cleaned, sep=""))
  write.csv(df.clin, file = paste(drive.datasuppository, date.cleaned, "_MGB_BladderCa_YoungPts_clinical_clean.csv",sep="") )
  write.csv(df.onc,  file = paste(drive.datasuppository, date.cleaned, "_MGB_BladderCa_YoungPts_w_oncopanel_clean.csv",sep="") )
}