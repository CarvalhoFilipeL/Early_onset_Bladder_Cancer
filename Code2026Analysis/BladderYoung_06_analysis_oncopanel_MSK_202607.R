# ################################################################################################################################################## #
# ## MGB Young Bladder Cancer Patients ## ---------------------------------------------------------------------------------------------------------- #
# ## Analysis: MSK Comutation Plot ## -------------------------------------------------------------------------------------------------------------- #
# ################################################################################################################################################## #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
###  Last Edit 2026-07


# requires main script

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Prep work ------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

#clinical
df.clin.msk     <- df.MSK_bladder %>%  #all MSK IMPACT patients
  mutate( Tumor_Sample_Barcode = sample )
df.clin.msk.RC  <- df.clin.msk    %>%  #RC only
  filter(RC_TUR_filter == "0.RC")
df.clin.msk.TUR <- df.clin.msk    %>%  #TUR only
  filter(RC_TUR_filter == "1.TUR")

#mutations
df.mutations.MSK <- read_table("~[REPLACE]~") %>%  # need to put local filename in for placeholder "~[REPLACE]~"
  mutate( Tumor_Sample_Barcode = dbSNP_Val_Status )

df.mutations.MSK <- df.mutations.MSK %>%
  left_join( df.clin.msk[,c(76,75,2)], by = "Tumor_Sample_Barcode" ) %>%
  mutate( Tumor_Sample_Barcode = Tumor_Sample_Barcode,
          Start_Position = Start_Position) %>%
  rename( Hugo_Symbol = Hugo_Symbol,
          Chromosome = Chromosome,
          End_Position = End_Position,
          Reference_Allele = Reference_Allele,
          Tumor_Seq_Allele2 = Tumor_Seq_Allele2,
          Variant_Classification = Variant_Classification,
          Variant_Type = Variant_Type)
df.mutations.MSK.RC <- df.mutations.MSK %>%
  filter(RC_TUR_filter == "0.RC")
df.mutations.MSK.TUR <- df.mutations.MSK %>%
  filter(RC_TUR_filter == "1.TUR")

maf.msk.all <- read.maf(maf = df.mutations.MSK, clinicalData = df.clin.msk)
maf.msk.RC  <- read.maf(maf = df.mutations.MSK.RC, clinicalData = df.clin.msk.RC)
maf.msk.TUR <- read.maf(maf = df.mutations.MSK.TUR, clinicalData = df.clin.msk.TUR)

#early-onset
maf.msk.young.all <- read.maf(maf = df.mutations.MSK %>% filter(flag_age55 == 1),
                              clinicalData = df.clin.msk %>% filter(flag_age55 == 1))
maf.msk.young.RC  <- read.maf(maf = df.mutations.MSK.RC %>% filter(flag_age55 == 1),
                              clinicalData = df.clin.msk.RC %>% filter(flag_age55 == 1))
maf.msk.young.TUR <- read.maf(maf = df.mutations.MSK.TUR %>% filter(flag_age55 == 1),
                              clinicalData = df.clin.msk.TUR %>% filter(flag_age55 == 1))
#late-onset
maf.msk.old.all <- read.maf(maf = df.mutations.MSK %>% filter(flag_age55 == 0),
                              clinicalData = df.clin.msk %>% filter(flag_age55 == 0))
maf.msk.old.RC  <- read.maf(maf = df.mutations.MSK.RC %>% filter(flag_age55 == 0),
                              clinicalData = df.clin.msk.RC %>% filter(flag_age55 == 0))
maf.msk.old.TUR <- read.maf(maf = df.mutations.MSK.TUR %>% filter(flag_age55 == 0),
                              clinicalData = df.clin.msk.TUR %>% filter(flag_age55 == 0))

#all ages
if(flag.save.tabs){ # export gene mutations as csv
  write.csv(maf.msk.all@gene.summary,
            file = paste(folder.results, "Table_S3_1a_all_MSK_mutations_", date.analyzed, ".csv",sep="") )
  write.csv(maf.msk.RC@gene.summary,
            file = paste(folder.results, "Table_S3_1b_RC_MSK_mutations_", date.analyzed, ".csv",sep="") )
  write.csv(maf.msk.TUR@gene.summary,
            file = paste(folder.results, "Table_S3_1c_TUR_MSK_mutations_", date.analyzed, ".csv",sep="") )
  #whole file
  write.csv(maf.msk.all@data,
            file = paste(folder.results, "Table_S3_MSK_mutations_raw_", date.analyzed, ".csv",sep="") )
}
#early-onset
if(flag.save.tabs){ # export gene mutations as csv
  write.csv(maf.msk.young.all@gene.summary,
            file = paste(folder.results, "Table_S3_2a_young_all_MSK_mutations_", date.analyzed, ".csv",sep="") )
  write.csv(maf.msk.young.RC@gene.summary,
            file = paste(folder.results, "Table_S3_2b_young_RC_MSK_mutations_", date.analyzed, ".csv",sep="") )
  write.csv(maf.msk.young.TUR@gene.summary,
            file = paste(folder.results, "Table_S3_2c_young_TUR_MSK_mutations_", date.analyzed, ".csv",sep="") )
}
#late-onset
if(flag.save.tabs){ # export gene mutations as csv
  write.csv(maf.msk.old.all@gene.summary,
            file = paste(folder.results, "Table_S3_3a_old_all_MSK_mutations_", date.analyzed, ".csv",sep="") )
  write.csv(maf.msk.old.RC@gene.summary,
            file = paste(folder.results, "Table_S3_3b_old_RC_MSK_mutations_", date.analyzed, ".csv",sep="") )
  write.csv(maf.msk.old.TUR@gene.summary,
            file = paste(folder.results, "Table_S3_3c_old_TUR_MSK_mutations_", date.analyzed, ".csv",sep="") )
}

# Visualization ------------------------------------------------------------------------------------------------------------------------------------ #

#MSK RC
if(flag.save.figs){
  #early-onset
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_RC/", "maf_summary_MSK_RC_earlyonset_",
            date.analyzed, ".png", sep=""))
  plotmafSummary(maf = maf.msk.young.RC, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_RC/", "maf_comutationplot_100_MSK_RC_earlyonset_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.msk.young.RC, top = 100)
  dev.off()
  
  #late-onset
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_RC/", "maf_summary_MSK_RC_lateonset_",
            date.analyzed, ".png", sep=""))
  plotmafSummary(maf = maf.msk.old.RC, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_RC/", "maf_comutationplot_25_MSK_RC_lateonset_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.msk.old.RC, top = 25)
  dev.off()
}

#MSK TUR
if(flag.save.figs){
  #early-onset
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_TUR/", "maf_summary_MSK_TUR_earlyonset_",
            date.analyzed, ".png", sep=""))
  plotmafSummary(maf = maf.msk.young.TUR, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_TUR/", "maf_comutationplot_25_MSK_TUR_earlyonset_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.msk.young.TUR, top = 25)
  dev.off()
  
  #late-onset
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_TUR/", "maf_summary_MSK_TUR_lateonset_",
            date.analyzed, ".png", sep=""))
  plotmafSummary(maf = maf.msk.old.TUR, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/MSK_TUR/", "maf_comutationplot_25_MSK_TUR_lateonset_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.msk.old.TUR, top = 25)
  dev.off()
}
