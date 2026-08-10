# ################################################################################################################################################## #
# ## MGB Young Bladder Cancer Patients ## ---------------------------------------------------------------------------------------------------------- #
# ## Analysis: characterization ## ----------------------------------------------------------------------------------------------------------------- #
# ################################################################################################################################################## #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
###  Last Edit 2026-07


# requires main script

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- MGB: Institutional Cohort Characteristics ----------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #

table_one_vars <- colnames(df.MGB_bladder)[-c(1:3)]; table_one_vars

table_one_factors <- table_one_vars[-c(2:4,15:16)]; table_one_factors

#overall
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MGB_bladder,
                                   includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MGB_bladder,
                                     includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tbl <- left_join(temp_tab_overall, temp_tab_strat_age %>%
                        rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54"),
                      by = c("Variable"))

# first row N
temp.tbl[1,2] <- paste(temp.tbl[1,2], " (100.0)", sep="")
temp.tbl[1,3:5] <- temp.tbl[3:5,2]
temp.tbl <- temp.tbl[-c(2:5),]

rm(temp_tab_overall, temp_tab_strat_age)


# Alternate format (condensed for 2 levels, median IQRs for continuous) ---------------------------------------------------------------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MGB_bladder,
                                   includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MGB_bladder,
                                     includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tblB <- left_join(temp_tab_overall, temp_tab_strat_age %>%
                         rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54"),
                       by = c("Variable"))
rm(temp_tab_overall, temp_tab_strat_age)

temp.tblB <- temp.tblB[c(6:8,27:28),]

# Combine ------------------------------------------------------------------------------------------------------------------------------------------ #

temp.tbl <- bind_rows(temp.tbl, temp.tblB) %>%
  rownames_to_column("original_row")

# rearrange rows
temp.tbl <- temp.tbl[c(1:2,82,3,83,4,84,5:23,85,24,86,25:81),] %>%
  remove_rownames()

# replace these gene rows w zeros
temp.tbl[c(76, 82, 83, 86),c(3:6)] <- "0 (0.0)"  # KMT2A, CDKN2B, CCND1, EP300 (manual bug fix)

temp.tbl$original_row <- temp.tbl$Variable
temp.tbl$Variable <- c(
  "N (%)", "Age, mean (SD)", "Age, median [IQR]",
  "Follow up - months, mean (SD)", "Follow up - months, median [IQR]", "BMI, mean (SD)", "BMI, median [IQR]",
  "Male (%)", "Race = White (%)", "Smoking (%)", "Stage ≥T2, Muscle Invasive (%)",
  "Intravesical Treatment (%)", "None", "BCG", "Chemo", "Both",
  "Systemic Treatment (%)", "None", "Chemo", "Immunotherapy", "Both",
  "Neoadjuvant Chemo (%)", "Neoadjuvant Chemo Response (%)",
  "Ileal Conduit (%)","Robotic (%)",
  "Length of stay - days, mean (SD)", "Length of stay - days, median [IQR]",
  "Operative time - hours, mean (SD)", "Operative time - hours, median [IQR]",
  "Readmission (%)", "Ileus (%)", "Post-operative Transfusion (%)", "DVT/PE (%)", "Reoperation (%)",
  "Pathology, T Stage (%)", "T0", "Ta", "Tis", "T1", "T2", "T3", "T4",
  "Adjuvant Immunotherapy (%)",
  "Adjuvant Chemotherapy (%)", "None", "Cisplatin / Gemcitabine", "Paclitaxel / Gemcitabine", "MVAC",
  "Fluorouracil / Mitomycin", "Cisplatin / Etoposide",
  "BCG (%)", "Intravesical Chemo (%)", "Matched Oncopanel (%)",
  table_one_vars[c(28:60)])

first.term <- function(my.string){
  as.numeric(unlist(strsplit(my.string, " "))[1])
}
temp.tbl$Rank <- sapply(temp.tbl$Overall, first.term)
temp.tbl$Rank[1:53] <- 100
temp.tbl <- temp.tbl %>% arrange(desc(Rank))
temp.tbl$Rank <- NULL

tbl.MGB_full_characteristics_s1 <- temp.tbl
rm(temp.tblB, temp.tbl, first.term)

if(flag.save.tabs){ # export as csv
  write.csv(tbl.MGB_full_characteristics_s1,
            file = paste(folder.results, "Table_S1_MGB_full_characteristics_", date.analyzed, ".csv",sep="") )
}


# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- MSK IMPACT: Cohort Characteristics ------------------------------------------------------------------------------------------------------------ #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

# ------------------- #
# entire cohort first #
# ------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #

table_one_vars <- colnames(df.MSK_bladder %>% select(-RC_TUR_filter))[-c(1,14:17)]; table_one_vars

table_one_factors <- table_one_vars[-c(3)]; table_one_factors

#overall
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MSK_bladder,
                                   includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MSK_bladder,
                                     includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# split at age 55 ---------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_55 <- CreateTableOne(vars = table_one_vars,
                                    factorVars = table_one_factors,
                                    strata = "flag_age55",
                                    data = df.MSK_bladder,
                                    includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_55 <- temp_tab_strat_55[,c(1,3,2,4)]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tbl <- left_join(temp_tab_overall, temp_tab_strat_55 %>%
                        rename("Under 55" = "X1", "Over 55" = "X0", "p.55" = "p"),
                      by = c("Variable"))
temp.tbl <- left_join(temp.tbl, temp_tab_strat_age %>%
                        rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54", "p.age" = "p"),
                      by = c("Variable"))

# first row N
temp.tbl[1,3:4] <- c(paste(temp.tbl[1,3], " (",round(100*as.numeric(temp.tbl[1,3])/as.numeric(temp.tbl[1,2]),1),")", sep=""),
                     paste(temp.tbl[1,4], " (",round(100*as.numeric(temp.tbl[1,4])/as.numeric(temp.tbl[1,2]),1),")", sep=""))
temp.tbl[1,2] <- paste(temp.tbl[1,2], " (100.0)", sep="")
temp.tbl[1,6:8] <- temp.tbl[4:6,2]

temp.tbl <- temp.tbl[-c(2:6),]

rm(temp_tab_overall, temp_tab_strat_age, temp_tab_strat_55)


# Alternate format (condensed for 2 levels, median IQRs for continuous) ---------------------------------------------------------------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MSK_bladder,
                                   includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MSK_bladder,
                                     includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# split at age 55 ---------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_55 <- CreateTableOne(vars = table_one_vars,
                                    factorVars = table_one_factors,
                                    strata = "flag_age55",
                                    data = df.MSK_bladder,
                                    includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_55 <- temp_tab_strat_55[,c(1,3,2,4)]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tblB <- left_join(temp_tab_overall, temp_tab_strat_55 %>%
                        rename("Under 55" = "X1", "Over 55" = "X0", "p.55" = "p"),
                      by = c("Variable"))
temp.tblB <- left_join(temp.tblB, temp_tab_strat_age %>%
                        rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54", "p.age" = "p"),
                      by = c("Variable"))
rm(temp_tab_overall, temp_tab_strat_age, temp_tab_strat_55)

# Combine ------------------------------------------------------------------------------------------------------------------------------------------ #

temp.tbl <- bind_rows(temp.tbl, temp.tblB[7,]) %>%
  rownames_to_column("original_row")

# rearrange rows
temp.tbl <- temp.tbl[c(1:2,99,3:98),] %>%
  remove_rownames()

temp.tbl$original_row <- temp.tbl$Variable
temp.tbl$Variable <- c(
  "N (%)", "Age, mean (SD)", "Age, median [IQR]",
  "Male (%)", "Smoking (%)", "Stage ≥T2, Muscle Invasive (%)",
  "Intravesical Treatment (%)", "None", "BCG", "Chemo", "Both",
  "Systemic Treatment (%)", "None", "Chemo", "Immunotherapy", "Both",
  
  "Specimen Stage (%)", "Low Grade", "High Grade, Non-Invasive", "High Grade, Invasive", "Metastatic",
  "Specimen Type (%)","Cystectomy","Metastasis", "Partial Cystectomy", "TUR", "Urethral biopsy", "Urethrectomy",
  "Stage (%)", "Low Grade", "High Grade", "CIS", "Muscle-Invasive", "Metastatic",
  "Histology (%)", "Urothelial Carcinoma", "Squamous Cell", "Plasmacytoid", "Small Cell, Neuroendocrine", "Micropapillary","Glandular","Other",
  "Altered", str_sub(temp.tbl$original_row[c(44:99)], end = -9))

tbl.MSK_full_characteristics_s2 <- temp.tbl
rm(temp.tblB, temp.tbl)

# --------------------- #
# RC only cohort second #
# --------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #

table_one_vars <- colnames(df.MSK_bladder %>% select(-RC_TUR_filter))[-c(1,14:17)]; table_one_vars

table_one_factors <- table_one_vars[-c(3)]; table_one_factors

#overall
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MSK_bladder %>% filter(RC_TUR_filter == "0.RC"),
                                   includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MSK_bladder %>% filter(RC_TUR_filter == "0.RC"),
                                     includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# split at age 55 ---------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_55 <- CreateTableOne(vars = table_one_vars,
                                    factorVars = table_one_factors,
                                    strata = "flag_age55",
                                    data = df.MSK_bladder %>% filter(RC_TUR_filter == "0.RC"),
                                    includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_55 <- temp_tab_strat_55[,c(1,3,2,4)]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tbl <- left_join(temp_tab_overall, temp_tab_strat_55 %>%
                        rename("Under 55" = "X1", "Over 55" = "X0", "p.55" = "p"),
                      by = c("Variable"))
temp.tbl <- left_join(temp.tbl, temp_tab_strat_age %>%
                        rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54", "p.age" = "p"),
                      by = c("Variable"))

# first row N
temp.tbl[1,3:4] <- c(paste(temp.tbl[1,3], " (",round(100*as.numeric(temp.tbl[1,3])/as.numeric(temp.tbl[1,2]),1),")", sep=""),
                     paste(temp.tbl[1,4], " (",round(100*as.numeric(temp.tbl[1,4])/as.numeric(temp.tbl[1,2]),1),")", sep=""))
temp.tbl[1,2] <- paste(temp.tbl[1,2], " (100.0)", sep="")
temp.tbl[1,6:8] <- temp.tbl[4:6,2]

temp.tbl <- temp.tbl[-c(2:6),]

rm(temp_tab_overall, temp_tab_strat_age, temp_tab_strat_55)


# Alternate format (condensed for 2 levels, median IQRs for continuous) ---------------------------------------------------------------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MSK_bladder %>% filter(RC_TUR_filter == "0.RC"),
                                   includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MSK_bladder %>% filter(RC_TUR_filter == "0.RC"),
                                     includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# split at age 55 ---------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_55 <- CreateTableOne(vars = table_one_vars,
                                    factorVars = table_one_factors,
                                    strata = "flag_age55",
                                    data = df.MSK_bladder %>% filter(RC_TUR_filter == "0.RC"),
                                    includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_55 <- temp_tab_strat_55[,c(1,3,2,4)]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tblB <- left_join(temp_tab_overall, temp_tab_strat_55 %>%
                         rename("Under 55" = "X1", "Over 55" = "X0", "p.55" = "p"),
                       by = c("Variable"))
temp.tblB <- left_join(temp.tblB, temp_tab_strat_age %>%
                         rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54", "p.age" = "p"),
                       by = c("Variable"))
rm(temp_tab_overall, temp_tab_strat_age, temp_tab_strat_55)

# Combine ------------------------------------------------------------------------------------------------------------------------------------------ #

temp.tbl <- bind_rows(temp.tbl, temp.tblB[7,]) %>%
  rownames_to_column("original_row")

# rearrange rows
temp.tbl <- temp.tbl[c(1:2,91,3:90),] %>%
  remove_rownames()

# replace these gene rows w zeros
temp.tbl[c(75, 85, 88),c(3:5,7:9)] <- "0 (0.0)"  # ARID2, POLE, MTOR (manual bug fix)

temp.tbl$original_row <- temp.tbl$Variable
temp.tbl$Variable <- c(
  "N (%)", "Age, mean (SD)", "Age, median [IQR]",
  "Male (%)", "Smoking (%)", "Stage ≥T2, Muscle Invasive (%)",
  "Intravesical Treatment (%)", "None", "BCG", "Chemo", "Both",
  "Systemic Treatment (%)", "None", "Chemo", "Immunotherapy", "Both",
  
  "Specimen Stage (%)", "Low Grade", "High Grade, Non-Invasive", "High Grade, Invasive",
  "Cystectomy",
  "Stage (%)", "Low Grade", "High Grade", "CIS", "Muscle-Invasive",
  "Histology (%)", "Urothelial Carcinoma", "Squamous Cell", "Plasmacytoid", "Small Cell, Neuroendocrine", "Micropapillary","Glandular","Other",
  str_sub(temp.tbl$original_row[c(35:91)], end = -9))

tbl.MSK_full_characteristics_s2a_RC_only <- temp.tbl
rm(temp.tblB, temp.tbl)

# --------------------- #
# TUR only cohort third #
# --------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #

table_one_vars <- colnames(df.MSK_bladder %>% select(-RC_TUR_filter))[-c(1,14:17)]; table_one_vars

table_one_factors <- table_one_vars[-c(3)]; table_one_factors

#overall
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MSK_bladder %>% filter(RC_TUR_filter == "1.TUR"),
                                   includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MSK_bladder %>% filter(RC_TUR_filter == "1.TUR"),
                                     includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# split at age 55 ---------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_55 <- CreateTableOne(vars = table_one_vars,
                                    factorVars = table_one_factors,
                                    strata = "flag_age55",
                                    data = df.MSK_bladder %>% filter(RC_TUR_filter == "1.TUR"),
                                    includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = NULL, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = NULL) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_55 <- temp_tab_strat_55[,c(1,3,2,4)]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tbl <- left_join(temp_tab_overall, temp_tab_strat_55 %>%
                        rename("Under 55" = "X1", "Over 55" = "X0", "p.55" = "p"),
                      by = c("Variable"))
temp.tbl <- left_join(temp.tbl, temp_tab_strat_age %>%
                        rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54", "p.age" = "p"),
                      by = c("Variable"))

# first row N
temp.tbl[1,3:4] <- c(paste(temp.tbl[1,3], " (",round(100*as.numeric(temp.tbl[1,3])/as.numeric(temp.tbl[1,2]),1),")", sep=""),
                     paste(temp.tbl[1,4], " (",round(100*as.numeric(temp.tbl[1,4])/as.numeric(temp.tbl[1,2]),1),")", sep=""))
temp.tbl[1,2] <- paste(temp.tbl[1,2], " (100.0)", sep="")
temp.tbl[1,6:8] <- temp.tbl[4:6,2]

temp.tbl <- temp.tbl[-c(2:6),]

rm(temp_tab_overall, temp_tab_strat_age, temp_tab_strat_55)


# Alternate format (condensed for 2 levels, median IQRs for continuous) ---------------------------------------------------------------------------- #

# Overall ------------------------------------------------------------------------------------------------------------------------------------------ #
temp_tab_overall <- CreateTableOne(vars = table_one_vars,
                                   factorVars = table_one_factors,
                                   data = df.MSK_bladder %>% filter(RC_TUR_filter == "1.TUR"),
                                   includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")

# Age bins ----------------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_age <- CreateTableOne(vars = table_one_vars,
                                     factorVars = table_one_factors,
                                     strata = "Age_Dx.bin",
                                     data = df.MSK_bladder %>% filter(RC_TUR_filter == "1.TUR"),
                                     includeNA = FALSE) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_age <- temp_tab_strat_age[,1:5]

# split at age 55 ---------------------------------------------------------------------------------------------------------------------------------- #
temp_tab_strat_55 <- CreateTableOne(vars = table_one_vars,
                                    factorVars = table_one_factors,
                                    strata = "flag_age55",
                                    data = df.MSK_bladder %>% filter(RC_TUR_filter == "1.TUR"),
                                    includeNA = F) %>%
  print(printToggle = FALSE,
        catDigits = 1, contDigits = 2, pDigits = 3, format = c("fp", "f", "p", "pf")[1],
        quote = FALSE, explain = TRUE, test = TRUE, noSpaces = TRUE, padColnames = FALSE, varLabels = FALSE,
        cramVars = table_one_factors, showAllLevels = F, dropEqual = FALSE,
        exact = NULL, nonnormal = setdiff(table_one_vars, table_one_factors)) %>%
  data.frame() %>% rownames_to_column("Variable")
temp_tab_strat_55 <- temp_tab_strat_55[,c(1,3,2,4)]

# Formatting --------------------------------------------------------------------------------------------------------------------------------------- #
temp.tblB <- left_join(temp_tab_overall, temp_tab_strat_55 %>%
                         rename("Under 55" = "X1", "Over 55" = "X0", "p.55" = "p"),
                       by = c("Variable"))
temp.tblB <- left_join(temp.tblB, temp_tab_strat_age %>%
                         rename("Under.45" = "X1..45", "45.to.49" = "X2.45.49", "50.to.54" = "X3.50.54", "p.age" = "p"),
                       by = c("Variable"))
rm(temp_tab_overall, temp_tab_strat_age, temp_tab_strat_55)

# Combine ------------------------------------------------------------------------------------------------------------------------------------------ #

temp.tbl <- bind_rows(temp.tbl, temp.tblB[7,]) %>%
  rownames_to_column("original_row")

# rearrange rows
temp.tbl <- temp.tbl[c(1:2,91,3:90),] %>%
  remove_rownames()

# replace these gene rows w zeros
temp.tbl[c(75),c(3:5,7:9)] <- "0 (0.0)"  # ARID2 (manual bug fix)

temp.tbl$original_row <- temp.tbl$Variable
temp.tbl$Variable <- c(
  "N (%)", "Age, mean (SD)", "Age, median [IQR]",
  "Male (%)", "Smoking (%)", "Stage ≥T2, Muscle Invasive (%)",
  "Intravesical Treatment (%)", "None", "BCG", "Chemo", "Both",
  "Systemic Treatment (%)", "None", "Chemo", "Immunotherapy", "Both",
  
  "Specimen Stage (%)", "Low Grade", "High Grade, Non-Invasive", "High Grade, Invasive",
  "TUR",
  "Stage (%)", "Low Grade", "High Grade", "CIS", "Muscle-Invasive",
  "Histology (%)", "Urothelial Carcinoma", "Squamous Cell", "Plasmacytoid", "Small Cell, Neuroendocrine", "Micropapillary","Glandular","Other",
  str_sub(temp.tbl$original_row[c(35:91)], end = -9))

tbl.MSK_full_characteristics_s2b_TUR_only <- temp.tbl
rm(temp.tblB, temp.tbl)

# tbl.MSK_full_characteristics_s2 is entire cohort
# tbl.MSK_full_characteristics_s2a_RC_only is only RC patients in cohort
# tbl.MSK_full_characteristics_s2b_TUR_only is only TUR patients in cohort

if(flag.save.tabs){ # export as csv
  write.csv(tbl.MSK_full_characteristics_s2,
            file = paste(folder.results, "Table_S2_MSK_full_characteristics_", date.analyzed, ".csv",sep="") )
  
  write.csv(tbl.MSK_full_characteristics_s2a_RC_only,
            file = paste(folder.results, "Table_S2a_MSK_full_characteristics_RC_", date.analyzed, ".csv",sep="") )
  
  write.csv(tbl.MSK_full_characteristics_s2b_TUR_only,
            file = paste(folder.results, "Table_S2b_MSK_full_characteristics_TUR_", date.analyzed, ".csv",sep="") )
}
