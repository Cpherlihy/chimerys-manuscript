#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-E9")
figurePath <- file.path(dataPath, "data/figure-E9")

filePaths <- file.path(dataPath, c("LFQ_Bench_multispecies/DDA/Chimerys/20240517_lfq_dda_z1to4_True_noNormalization.pdResult",
                                   "LFQ_Bench_multispecies/DIA/Chimerys/20240517_lfq_dia_z1to4_v2x7x9_apex_True.pdResult",
                                   "LFQ_Bench_multispecies/DIA/Chimerys/20240522_LFQ3-DIA-Minora_noNorm.pdResult"))
sampleNames <- writeSampleNames(filePaths, outputPath = path, outputName = "sample_names_LFQ3_density.csv")
dtPtmLfq <- readPtmGroups(filePaths, sampleNames, loadBackup = F, writeBackup = F)
dtPtmLfq[, condition_MS := factor(gsub("^(.*)-.*-.*$", "\\1", condition), c("DDA", "DIA"))]
dtPtmLfq[, condition_Quan := factor(gsub("^.*-(.*)-.*$", "\\1", condition), c("CHIMERYS", "Minora"))]

#filter organism
organismLevels <-
  c("Escherichia coli (strain K12)", "Homo sapiens", "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)")
organismLabels <- c("E. coli", "Human", "Yeast")
organismRatios <- setNames(log2(c(0.25, 1, 2)), organismLabels)
dtPtmLfq <- dtPtmLfq[organism %in% organismLevels]
dtPtmLfq[, organism := factor(droplevels(organism), organismLevels, organismLabels)]
dtPtmLfq[, .N, keyby=.(condition_MS, condition_Quan, organism)]

#calculate ratios
processIntensity(dtPtmLfq, "intensity")
contrasts <- c("DDA-Minora-CondA_vs_DDA-Minora-CondB",
               "DIA-Minora-CondA_vs_DIA-Minora-CondB",
               "DIA-CHIMERYS-CondA_vs_DIA-CHIMERYS-CondB")
contrastLabels <- c("DDA (Minora MS1 Quan)",
                    "DIA (Minora MS1 Quan)",
                    "DIA (CHIMERYS MS2 Quan)")
lsRatios <- testContrasts(dataTable = dtPtmLfq, observationColumn = "ptm_group_J", contrasts, minQuanReplicates = 2L)
lsRatios[, contrastLabel := factor(contrast, contrasts, contrastLabels)]
lsRatios[, .N, keyby=contrast]
dtMaLines <- data.table(YINTERCEPT = organismRatios, organism = factor(organismLabels))
msaid_organism <- c("Human" = msaidTheme$blue, "Yeast" = msaidTheme$orange, "E. coli" = msaidTheme$grayDark)
lsRatios[, hasMethionine := ifelse(grepl("M", ptm_group_J), "has Met", "no Met")]

fwrite(lsRatios, file.path(figurePath, "figure-E9C-density.csv"))
