# ################################################################################################################################################## #
# ## MGB Young Bladder Cancer Patients ## ---------------------------------------------------------------------------------------------------------- #
# ## Analysis: Bar Plot ## ------------------------------------------------------------------------------------------------------------------------- #
# ################################################################################################################################################## #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
###  Last Edit 2026-07


# requires main script

# -------------------------------------------------------------------------------------------------------------------------------------------------- #
# -- Bar Plot -------------------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #

list.gene.filter <- c("TP53", "KMT2D", "ARID1A", "RB1", "TERT", "FGFR3")

temp.1.mgb <- tbl.MGB_full_characteristics_s1[54:86,c(2:3)] %>%
  rename("Gene" = Variable, "MGB < 55" = Overall)

first.term <- function(my.string){
  as.numeric(unlist(strsplit(my.string, " "))[1])
}
second.term <- function(my.string){
  unlist(strsplit(my.string, " "))[2]
}

# temp.1.mgb$N.MGB    <- sapply(temp.1.mgb$"MGB Overall", first.term)
temp.1.mgb$"MGB < 55" <- as.numeric(gsub("[()]", "", sapply(temp.1.mgb$"MGB < 55", second.term)))
temp.1.mgb <- temp.1.mgb %>%
  filter(Gene %in% list.gene.filter)

rm(first.term, second.term)

#MSK whole cohort
temp.2.msk <- tbl.MSK_genes[,c(11,1,14,16,18,22:24)] %>%
  rename("Gene" = Variable,
         "MSK Overall" = Percent.Overall,
         "MSK < 55" = Percent.Under.55, "MSK >= 55" = Percent.Over.55,
         "MSK < 45" = Percent.Under.45, "MSK 45-49" = Percent.45.to.49, "MSK 50-54" = Percent.50.to.54) %>%
  filter(Gene %in% list.gene.filter)
#MSK RC cohort
temp.3.msk.rc <- tbl.MSK_genes_RC[,c(11,1,14,16,18,22:24)] %>%
  rename("Gene" = Variable,
         "MSK RC Overall" = Percent.Overall,
         "MSK RC < 55" = Percent.Under.55, "MSK RC >= 55" = Percent.Over.55) %>%
  filter(Gene %in% list.gene.filter)
#MSK TUR cohort
temp.3.msk.tur <- tbl.MSK_genes_TUR[,c(11,1,14,16,18,22:24)] %>%
  rename("Gene" = Variable,
         "MSK TUR Overall" = Percent.Overall,
         "MSK TUR < 55" = Percent.Under.55, "MSK TUR >= 55" = Percent.Over.55) %>%
  filter(Gene %in% list.gene.filter)

#combine tables
temp.tbl <- left_join(temp.1.mgb, temp.2.msk, by = "Gene") %>%
  left_join(temp.3.msk.rc, by = "Gene") %>%
  left_join(temp.3.msk.tur, by = "Gene")

temp.tbl <- temp.tbl[,c(1,2,12,13,19,20)] %>% pivot_longer(cols = c(2:6), names_to = 'Cohort', values_to = 'Percent') %>%
  mutate(Gene = factor(Gene, levels = list.gene.filter),
         Cohort = factor(Cohort, levels = c("MGB < 55",
                                            "MSK RC < 55",
                                            "MSK RC >= 55",
                                            "MSK TUR < 55",
                                            "MSK TUR >= 55")))

temp.fig <- ggplot(data = temp.tbl,
                   aes(x = Gene, y = Percent, fill = Cohort)) +
  # ggtitle("Figure 1: Gene Mutation Rates") +
  geom_bar(position = 'dodge',
           stat = 'identity') +
  # geom_text(data=df.US.census, aes(x = race, y = census, label = format(census, nsmall = 1)), size = 4.5, vjust = -0.5) +
  scale_fill_manual(values = c('azure4','darkseagreen','darkolivegreen','lightskyblue','steelblue')) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25), limits = c(0,105)) + 
  labs(x = 'Gene', y = '\nMutation Rate, %') + 
  theme_bw() +
  theme(plot.title = element_text(family = "serif", face = "bold", size = 20),
        panel.grid.major.x = element_blank(), # remove vertical gridlines
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = c(0.9,0.75), legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  )
temp.fig

if(flag.save.figs){
  png(paste(folder.figures, date.analyzed, '_fig_bargraph_genemutations.png', sep=''),
      width=7, height=7, units = 'in', res=300)
  temp.fig
  dev.off()
}
rm(temp.tbl, temp.1.mgb, temp.2.msk, temp.3.msk.rc, temp.3.msk.tur, temp.fig)

