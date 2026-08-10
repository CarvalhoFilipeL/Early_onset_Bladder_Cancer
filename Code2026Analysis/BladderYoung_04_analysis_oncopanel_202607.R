# ################################################################################################################################################## #
# ## MGB Young Bladder Cancer Patients ## ---------------------------------------------------------------------------------------------------------- #
# ## Analysis: Oncopanel ## ------------------------------------------------------------------------------------------------------------------------ #
# ################################################################################################################################################## #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
###  Last Edit 2026-07


# requires main script

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Prep work ------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# raw files
#clinical: df.clin
df.clinical <- df.clin %>% filter( Oncopanel_match == 1 ) %>%
  mutate( Tumor_Sample_Barcode = EPIC_PMRN )
#mutations: df.onc
df.mutations <- df.onc %>%
  filter( Oncopanel_match == 1 ) %>%
  mutate( Tumor_Sample_Barcode = EPIC_PMRN,
          Start_Position = POSITION) %>%
  rename( Hugo_Symbol = HARMONIZED_HUGO_GENE_NAME,
          Chromosome = CHROMOSOME,
          End_Position = POSITION,
          Reference_Allele = REF_ALLELE,
          Tumor_Seq_Allele2 = ALT_ALLELE,
          Variant_Classification = HARMONIZED_VARIANT_CLASS,
          Variant_Type = VARIANT_TYPE)


maf.youngpts <- read.maf(maf = df.mutations, clinicalData = df.clinical)

if(flag.save.tabs){ # export gene mutations as csv
  write.csv(maf.youngpts@gene.summary,
            file = paste(folder.results, "Table_S3_MGB_mutations_", date.analyzed, ".csv",sep="") )
}
# Visualization ------------------------------------------------------------------------------------------------------------------------------------ #
maf.youngpts
getGeneSummary(maf.youngpts)

if(flag.save.figs){
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/", "maf_summary_",
            date.analyzed, ".png", sep=""))
  plotmafSummary(maf = maf.youngpts, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/", "maf_oncoplot_100_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.youngpts, top = 100)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/", "maf_oncoplot_50_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.youngpts, top = 50)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/", "maf_oncoplot_25_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.youngpts, top = 25)
  dev.off()
  
  png(paste(folder.figures, "MafTools_Visualization/UpdatedRevision202607/", "maf_oncoplot_15_",
            date.analyzed, ".png", sep=""))
  oncoplot(maf = maf.youngpts, top = 15)
  dev.off()
}