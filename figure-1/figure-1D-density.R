#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-1")
densityPath <- file.path(dataPath, "figure-1/density/PD")

#load PSMs
filePath <-
  file.path(densityPath, c("20240517_lfq_dda_z1to4_True_noNormalization.pdResult",
                          "LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01_noNorm.pdResult"))
sampleNamesPath <- writeSampleNames(filePath, outputPath = path, outputName = "density-sample-names.csv")
data <- readPtmGroups(dataPath = filePath, sampleNamesPath = sampleNamesPath)
#data[, .N, keyby=condition]

#conditions
data[, condition_software := factor(gsub("(.*)-.*", "\\1", condition), c("CHIMERYS", "Sequest"))]

#filter organism
organismLevels <-
  c("Escherichia coli (strain K12)", "Homo sapiens", "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)")
organismLabels <- c("E. coli", "Human", "Yeast")
organismRatios <- setNames(log2(c(0.25, 1, 2)), organismLabels)
data <- data[organism %in% organismLevels]
data[, organism := factor(droplevels(organism), organismLevels, organismLabels)]
data[, .N, keyby=.(condition_software, organism)]

#calculate ratios
processIntensity(data, "intensity")
contrasts <- c("CHIMERYS-CondA_vs_CHIMERYS-CondB", "Sequest-CondA_vs_Sequest-CondB")
contrastLabels <- c("CHIMERYS (Minora MS1 Quan)", "Sequest-HT (Minora MS1 Quan)")
lsRatios <- testContrasts(dataTable = data, observationColumn = "ptm_group_J", contrasts, minQuanReplicates = 2L)
lsRatios[, contrastLabel := factor(contrast, contrasts, contrastLabels)]
lsRatios[, .N, keyby=contrast]
dtMaLines <- data.table(YINTERCEPT = organismRatios, organism = factor(organismLabels))
msaid_organism <- c("Human" = msaid_blue, "Yeast" = msaid_orange, "E. coli" = msaid_darkgray)
lsRatios[, hasMethionine := ifelse(grepl("M", ptm_group_J), "has Met", "no Met")]

fwrite(lsRatios, file.path(file.path(dataPath, "figure-1"), "density.csv"))

ggplot(lsRatios, aes(x=ratio, color=organism)) +
  geom_density(linewidth=0.25) +
  geom_vline(data=dtMaLines, aes(xintercept=YINTERCEPT, color=organism),
             linetype = "dashed", linewidth = 0.25, show.legend = F) +
  scale_color_manual(NULL, values = msaid_organism) +
  scale_x_continuous(breaks = pretty_breaks(), limits = c(-4, 3)) +
  guides(fill = guide_legend(override.aes = list(color = NA, size = 2))) +
  facet_grid(rows = vars(contrastLabel), cols = vars(hasMethionine)) +
  xlab("Log2 fold change") + ylab("Density") +
  theme(legend.position = "right", plot.background = element_rect(fill = "transparent", colour = NA))
