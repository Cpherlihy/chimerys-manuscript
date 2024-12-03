suppressMessages(source(here::here("scripts/load-dependencies.R")))
path <- file.path(here::here(), "figure-6")
figurePath <- file.path(dataPath, "data/figure-6")


##LFQ3_IDs_CD.csv ====
#GOAL: global peptide group IDs at FDR, filled by "Quantified in >=2 replicates" (same as 3A)
filePath <- file.path(dataPath, c("LFQ_Bench_multispecies/DDA/Chimerys/20240517_lfq_dda_z1to4_True_noNormalization.pdResult",
                                  "LFQ_Bench_multispecies/DIA/Chimerys/20240517_lfq_dia_z1to4_v2x7x9_apex_True.pdResult"))
sampleNamesPath <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names_LFQ3.csv")
data_lfq3 <- readPtmGroups(filePath, sampleNamesPath, loadBackup=F, writeBackup=F)

#conditions
data_lfq3[, condition_MS := factor(gsub("(.*)-.*", "\\1", condition), c("DDA", "DIA"))]

#calculate ratios
processIntensity(data_lfq3, "intensity")
contrasts <- c("DDA-CondA_vs_DDA-CondB", "DIA-CondA_vs_DIA-CondB")
contrastLabels <- c("DDA (Minora MS1 Quan)", "DIA (CHIMERYS MS2 Quan)")
dtRatiosMin2_lfq3 <- testContrasts(data_lfq3, "ptm_group_J", contrasts, minQuanReplicates = 2L)
dtRatiosMin0_lfq3 <- testContrasts(data_lfq3, "ptm_group_J", contrasts, minQuanReplicates = 0L)
dtRatiosMin0_lfq3[, isQuanMin2Each := nQuan1st>=2 & nQuan2nd>=2]
dtObsMin3_lfq3 <- testContrasts(data_lfq3, "ptm_group_J", contrasts, minQuanReplicates = 3L, output = "observations")
dtIdPtmGroup_lfq3 <- dtRatiosMin0_lfq3[, .N, keyby = .(contrast, isQuanMin2Each)]
dtIdPtmGroup_lfq3[, contrastLabel := factor(contrast, contrasts, contrastLabels)]

fwrite(dtIdPtmGroup_lfq3, file.path(figurePath, "LFQ3_IDs_CD.csv"))


##Astral_14min_IDs_CD.csv ====
filePath <- file.path(dataPath, c("Astral_14min/DDA/20240517_astral_dda_z1to4_True_noNormalization.pdResult",
                                  "Astral_14min/DIA/20240517_astral_dia_z1to4_True.pdResult"))
sampleNamesPath <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names_Astral_14min.csv")
data_Ast14 <- readPtmGroups(filePath, sampleNamesPath)

#check MBR peptide groups
data_Ast14[, .N, keyby=.(condition, qvalue_psm_label)]

#conditions
data_Ast14[, condition_MS := factor(gsub("(.*)-.*", "\\1", condition), c("DDA", "DIA"))]
data_Ast14[, condition_method := factor(gsub(".*-(.*)", "\\1", condition), c("cycle", "top100", "top75", "top50"))]
data_Ast14[condition_MS=="DIA" & is.na(condition_method), condition_method := "cycle"]

#filter conditions needed
data_Ast14 <- data_Ast14[condition %in% c("DIA-cycle", "DDA-cycle")]

#calculate CVs and mean intensities
processIntensity(data_Ast14, "intensity")
dtCvsMin2_Ast14 <- testConditions(data_Ast14, "ptm_group_J", minQuanReplicates = 2L)
dtCvsMin0_Ast14 <- testConditions(data_Ast14, "ptm_group_J", minQuanReplicates = 0L, removeNaIntensity = F)
dtCvsMin0_Ast14[, condition_MS := factor(gsub("(.*)-.*", "\\1", condition), c("DDA", "DIA"))]
dtCvsMin0Unique_Ast14 <- unique(dtCvsMin0_Ast14[, .SD, .SDcols=c("condition_MS", "ptm_group_J", "nQuan")])
dtCvsMin0Unique_Ast14[, isQuanMin2 := nQuan>=2]
dtIdPtmGroup_Ast14 <- dtCvsMin0Unique_Ast14[, .N, keyby = .(condition_MS, isQuanMin2)]

fwrite(dtIdPtmGroup_Ast14, file.path(figurePath, "Astral_14min_IDs_CD.csv"))


##Astral_30min_IDs_CD.csv ====
filePath <- file.path(dataPath, c("Astral_30min/DDA/20240517_astral_dda_30min_z1to4_v2x7x9_apex_noNormalization_True.pdResult",
                                  "Astral_30min/DIA/20240517_astral_dia_30min_z1to4_v2x7x9_apex_True.pdResult"))
sampleNamesPath <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names_Astral_30min.csv")
data_Ast30 <- readPtmGroups(filePath, sampleNamesPath, loadBackup = T)

#check MBR peptide groups
data_Ast30[, .N, keyby=.(condition, qvalue_psm_label)]
#OPTIONAL: exclude MBR
#data_Ast30[qvalue_psm_label!="<=0.01", intensity := NA]

#conditions
data_Ast30[, condition_MS := factor(gsub("(.*)-.*", "\\1", condition), c("DDA", "DIA"))]
data_Ast30[, condition_method := factor(gsub(".*-(.*)", "\\1", condition), c("cycle", "2p5", "3p5"))]

#filter conditions needed
data_Ast30 <- data_Ast30[condition %in% c("DIA-2p5", "DDA-cycle")]

#calculate CVs and mean intensities
processIntensity(data_Ast30, "intensity")
dtCvsMin2_Ast30 <- testConditions(data_Ast30, "ptm_group_J", minQuanReplicates = 2L)
dtCvsMin0_Ast30 <- testConditions(data_Ast30, "ptm_group_J", minQuanReplicates = 0L, removeNaIntensity = F)
dtCvsMin0_Ast30[, condition_MS := factor(gsub("(.*)-.*", "\\1", condition), c("DDA", "DIA"))]
dtCvsMin0Unique_Ast30 <- unique(dtCvsMin0_Ast30[, .SD, .SDcols=c("condition_MS", "ptm_group_J", "nQuan")])
dtCvsMin0Unique_Ast30[, isQuanMin2 := nQuan>=2]
dtIdPtmGroup_Ast30 <- dtCvsMin0Unique_Ast30[, .N, keyby = .(condition_MS, isQuanMin2)]

fwrite(dtIdPtmGroup_Ast30, file.path(figurePath, "Astral_30min_IDs_CD.csv"))


##LFQ3_CVs.csv ====
#melt CVs and filter on same peptides
dtCv_lfq3 <- melt(dtRatiosMin2_lfq3, measure.vars = c("cv1st", "cv2nd"),
                  id.vars = c("contrast", "ptm_group_J"),
                  variable.name = "condition", value.name = "cv")
dtCv_lfq3[, condition := factor(condition, c("cv1st", "cv2nd"), c("CondB", "CondA"))]
dtCv_lfq3[, contrastLabel := factor(contrast, contrasts, contrastLabels)]
setkey(dtCv_lfq3, contrast, ptm_group_J)

dtCvShared_lfq3 <- dtCv_lfq3[dtCv_lfq3[, .N, by=ptm_group_J], on="ptm_group_J"][N==4 & !grepl("M", ptm_group_J)]
dtCvShared_lfq3[, .(mean_CV = mean(cv), median_CV = median(cv), .N), keyby = .(contrastLabel)]

fwrite(dtCvShared_lfq3, file.path(figurePath, "LFQ3_CVs.csv"))


##Astral_14min_CVs.csv ====
#make unique filter on same peptides
dtCvsMin2_Ast14[, condition_MS := factor(gsub("(.*)-.*", "\\1", condition), c("DDA", "DIA"))]
dtCv_Ast14 <- unique(dtCvsMin2_Ast14[, .SD, .SDcols=c("condition_MS", "ptm_group_J", "cv")])
setkey(dtCv_Ast14, condition_MS, ptm_group_J)
dtCvShared_Ast14 <- dtCv_Ast14[dtCv_Ast14[, .N, by=ptm_group_J], on="ptm_group_J"][N==2 & !grepl("M", ptm_group_J)]

fwrite(dtCvShared_Ast14, file.path(figurePath, "Astral_14min_CVs.csv"))


##Astral_30min_CVs.csv ====
#make unique filter on same peptides
dtCvsMin2_Ast30[, condition_MS := factor(gsub("(.*)-.*", "\\1", condition), c("DDA", "DIA"))]
dtCv_Ast30 <- unique(dtCvsMin2_Ast30[, .SD, .SDcols=c("condition_MS", "ptm_group_J", "cv")])
setkey(dtCv_Ast30, condition_MS, ptm_group_J)

dtCvShared_Ast30 <- dtCv_Ast30[dtCv_Ast30[, .N, by=ptm_group_J], on="ptm_group_J"][N==2 & !grepl("M", ptm_group_J)]

fwrite(dtCvShared_Ast30, file.path(figurePath, "Astral_30min_CVs.csv"))



##LFQ3_MA.csv ====
#filter on same peptides
dtMaShared_lfq3 <- dtRatiosMin2_lfq3[dtRatiosMin2_lfq3[, .N, by=ptm_group_J],
                                     on="ptm_group_J"][N==2 & !grepl("M", ptm_group_J)]
setkey(dtMaShared_lfq3, contrast, ptm_group_J)

#plot MA
dtOrganism_lfq3 <- data_lfq3[, .(organism = organism[1]), keyby=ptm_group_J]
organismLevels <-
  c("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", "Homo sapiens", "Escherichia coli (strain K12)")
organismLabels <- c("Yeast", "Human", "E. coli")

dtOrganism_lfq3 <- dtOrganism_lfq3[organism %in% organismLevels]
dtOrganism_lfq3[, organism := factor(droplevels(organism), organismLevels, organismLabels)]
dtOrganism_lfq3[, .N, keyby=organism]

dtMa_lfq3 <- dtOrganism_lfq3[dtMaShared_lfq3[!is.na(ratio) & !grepl("M", ptm_group_J)], on="ptm_group_J"]
dtMa_lfq3 <- dtMa_lfq3[!is.na(organism)]
dtMa_lfq3[, contrastLabel := factor(contrast, contrasts, contrastLabels)]
dtMa_lfq3[, .(.N, mean(ratio)), keyby=.(contrastLabel, organism)]

fwrite(dtMa_lfq3, file.path(figurePath, "LFQ3_MA.csv"))


##LFQ3_cor_DIA.csv ====
#GOAL: global peptide group IDs at FDR, filled by "Quantified in >=2 replicates" (same as 3A)
filePath <- file.path(dataPath, c("LFQ_Bench_multispecies/DIA/Chimerys/20240522_LFQ3-DIA-Minora_noNorm.pdResult",
                                  "LFQ_Bench_multispecies/DIA/Chimerys/20240517_lfq_dia_z1to4_v2x7x9_apex_True.pdResult"))
sampleNamesPath <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names_LFQ3_DIA.csv")
data_lfq3 <- readPtmGroups(filePath, sampleNamesPath, loadBackup=F, writeBackup=F)

#calculate ratios
processIntensity(data_lfq3, "intensity")
contrasts <- c("DIA-Minora-CondA_vs_DIA-Minora-CondB", "DIA-CHIMERYS-CondA_vs_DIA-CHIMERYS-CondB")
contrastLabels <- c("DDA (Minora MS1 Quan)", "DIA (CHIMERYS MS2 Quan)")
dtObsMin3_lfq3 <- testContrasts(data_lfq3, "ptm_group_J", contrasts, minQuanReplicates = 3L, output = "observations")

#GOAL: mean intensity of quantified in all groups
dtObsMin3Mean_lfq3 <- dtObsMin3_lfq3[, .(intensityMean = mean(intensity)), by=.(condition, ptm_group_J)]
dtObsMin3Mean_lfq3[, condition_MS := factor(gsub("^DIA-(.*)-.*", "\\1", condition), c("CHIMERYS", "Minora"))]
dtObsMin3Mean_lfq3[, id_cond := paste0(ptm_group_J, "_", gsub("DIA-.*-(.*)", "\\1", condition))]
dtObsMin3Mean_lfq3[, nConditionMS := .N, by=id_cond]
dtObsMin3Mean_lfq3 <- dtObsMin3Mean_lfq3[nConditionMS==2 & !grepl("M", ptm_group_J)]
dtObsMin3Mean_lfq3 <- dcast(dtObsMin3Mean_lfq3, id_cond~condition_MS, value.var = "intensityMean")

fwrite(dtObsMin3Mean_lfq3, file.path(figurePath, "LFQ3_cor_DIA.csv"))
