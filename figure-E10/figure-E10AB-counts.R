#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-E10")
figurePath <- file.path(dataPath, "data/figure-E10")


# LFQ3 ====
filePaths <- file.path(dataPath, c("LFQ_Bench_multispecies/DDA/Chimerys/20240517_lfq_dda_z1to4_True_noNormalization.pdResult",
                                  "LFQ_Bench_multispecies/DIA/Chimerys/20240517_lfq_dia_z1to4_v2x7x9_apex_True.pdResult"))
sampleNames <- writeSampleNames(filePaths, outputPath = path, outputName = "sample_names_LFQ3.csv")
dtPsmLfq <- readPsms(filePaths, sampleNames, loadBackup = F, writeBackup = F)
dtProtGrpLfq <- readProtGroups(filePaths, sampleNamesPath = sampleNames, loadBackup = F, writeBackup = F)
dtProtGrpLfq[, condition_MS := factor(gsub("^(.*)-.*$", "\\1", condition), c("DDA", "DIA"))]
dtProtGrpLfq[, condition_org := factor(gsub("^.*-(.*)$", "\\1", condition), c("CondA", "CondB"))]

#process intensities and calculate min 2 data completeness in each condition
processIntensity(dtProtGrpLfq, "intensity")
contrasts <- c("DDA-CondA_vs_DDA-CondB", "DIA-CondA_vs_DIA-CondB")
contrastLabels <- c("DDA (Minora MS1 Quan)", "DIA (CHIMERYS MS2 Quan)")

ratioProtGrpLfq <- testContrasts(dtProtGrpLfq, "protein", contrasts, minQuanReplicates = 0L)
ratioProtGrpLfq[, isQuanMin2Each := nQuan1st>=2 & nQuan2nd>=2]
countProtGrpLfq <- ratioProtGrpLfq[, .N, keyby = .(contrast, isQuanMin2Each)]
fwrite(countProtGrpLfq, file.path(figurePath, "figure-E10B-LFQ3.csv"))

countPsmLfq <- dtPsmLfq[, .N, keyby=.(condition, sample)]
fwrite(countPsmLfq, file.path(figurePath, "figure-E10A-LFQ3.csv"))


# Astral 14min ====
filePaths <- file.path(dataPath, c("Astral_14min/DDA/20240517_astral_dda_z1to4_True_noNormalization.pdResult",
                                   "Astral_14min/DIA/20240517_astral_dia_z1to4_True.pdResult"))
sampleNames <- writeSampleNames(filePaths, outputPath = path, outputName = "sample_names_Ast14.csv")
dtPsmAst14 <- readPsms(filePaths, sampleNames, loadBackup = F, writeBackup = F)
dtProtGrpAst14 <- readProtGroups(filePaths, sampleNamesPath = sampleNames, loadBackup = F, writeBackup = F)

#process intensities and calculate min 2 data completeness in each condition
processIntensity(dtProtGrpAst14, "intensity")

condProtGrpAst14 <- testConditions(dtProtGrpAst14, "protein", "condition", minQuanReplicates = 0L)
condProtGrpAst14 <- condProtGrpAst14[, .N, keyby=.(condition, protein, nQuan)]
condProtGrpAst14[, isQuanMin2 := nQuan>=2]
countProtGrpAst14 <- condProtGrpAst14[, .N, keyby = .(condition, isQuanMin2)]
fwrite(countProtGrpAst14, file.path(figurePath, "figure-E10B-Astral14.csv"))

countPsmAst14 <- dtPsmAst14[, .N, keyby=.(condition, sample)]
fwrite(countPsmAst14, file.path(figurePath, "figure-E10A-Astral14.csv"))


# Astral 30min ====
filePaths <- file.path(dataPath, c("Astral_30min/DDA/20240517_astral_dda_30min_z1to4_v2x7x9_apex_noNormalization_True.pdResult",
                                   "Astral_30min/DIA/20240517_astral_dia_30min_z1to4_v2x7x9_apex_True.pdResult"))
sampleNames <- writeSampleNames(filePaths, outputPath = path, outputName = "sample_names_Ast30.csv")
dtPsmAst30 <- readPsms(filePaths, sampleNames, loadBackup = F, writeBackup = F)
dtProtGrpAst30 <- readProtGroups(filePaths, sampleNamesPath = sampleNames, loadBackup = F, writeBackup = F)

#process intensities and calculate min 2 data completeness in each condition
processIntensity(dtProtGrpAst30, "intensity")

condProtGrpAst30 <- testConditions(dtProtGrpAst30, "protein", "condition", minQuanReplicates = 0L)
condProtGrpAst30 <- condProtGrpAst30[, .N, keyby=.(condition, protein, nQuan)]
condProtGrpAst30[, isQuanMin2 := nQuan>=2]
countProtGrpAst30 <- condProtGrpAst30[, .N, keyby = .(condition, isQuanMin2)]
fwrite(countProtGrpAst30, file.path(figurePath, "figure-E10B-Astral30.csv"))

countPsmAst30 <- dtPsmAst30[, .N, keyby=.(condition, sample)]
fwrite(countPsmAst30, file.path(figurePath, "figure-E10A-Astral30.csv"))


