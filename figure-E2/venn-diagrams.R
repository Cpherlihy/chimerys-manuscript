#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-S3-overlap")
figurePath <- file.path(dataPath, "data/figure-S3")
searchEnginePath <- file.path(dataPath, "LFQ_Bench_human")

#Load data
sampleNamesPath <- writeSampleNames(searchEnginePath, outputPath = path, outputName = "sample_names.csv")
data_psm <- readPsms(searchEnginePath, sampleNamesPath, loadBackup = F, writeBackup = F)
data_ptm <- readPtmGroups(searchEnginePath, sampleNamesPath, loadBackup = F, writeBackup = F)
data_pcm <- readPcms(searchEnginePath, sampleNamesPath, loadBackup = F, writeBackup = F)

#converge to PTM groups
parseModificationsColumn(dataTable = data_psm)
data_psm_ptm <- data_psm[condition %in% c("Metamorpheus", "MSGFplus"),
                         .N, keyby = c("condition", "sample", "ptm_group_J")]

parseModificationsColumn(dataTable = data_pcm)
data_pcm_ptm <- data_pcm[condition %in% c("MSFragger", "MSFragger-DDAplus"),
                         .N, keyby = c("condition", "sample", "ptm_group_J")]

#PSM level ====
psm_se <- c("MSF-200", "MS-GF-plus", "MaxQuant", "Metamorpheus", "MS-Amanda", "Comet", "Sequest-HT")
#data_psm[condition %in% psm_se, .N, keyby = c("condition")]
data_venn_psm_duell <-
  rbind(data_psm[condition %in% psm_se, .(condition = factor("7 other search engines"), .N), by = c("psm_J")],
        data_psm[condition %in% "CHIMERYS", .(condition = factor("CHIMERYS"), .N), by = c("psm_J")])
fwrite(data_venn_psm_duell, file.path(figurePath, "figure-E2A-venn-psms.csv"), quote = T)


#PTM group level - native only ====
ptm_native_se <- c("MaxQuant", "Comet", "SequestHT", "MSAmanda")
#data_ptm[condition %in% ptm_native_se, .N, keyby = c("condition")]
data_venn_ptm_native_duell <-
  rbind(data_psm[condition %in% ptm_native_se, .(condition = factor("4 other search engines"), .N), by = c("ptm_group_J")],
        data_ptm[condition %in% "CHIMERYS", .(condition = factor("CHIMERYS"), .N), by = c("ptm_group_J")])
fwrite(data_venn_ptm_native_duell, file.path(figurePath, "figure-E2B-venn-peptides-native.csv"), quote = T)


#PTM group level - all ====
ptm_native_se <- c("MaxQuant", "Comet", "SequestHT", "MSAmanda")
data_ptm_joined <-
  rbind(data_ptm[condition %in% ptm_native_se, .(condition = factor("7 other search engines"), .N), by = c("ptm_group_J")],
        data_pcm_ptm[condition %in% "MSF-200", .(condition = factor("7 other search engines"), .N), by = c("ptm_group_J")],
        data_psm_ptm[condition %in% c("Metamorpheus", "MS-GF-plus"),
                     .(condition = factor("7 other search engines"), .N), by = c("ptm_group_J")])
data_venn_ptm_all_duell <- rbind(data_ptm_joined,
                                 data_ptm[condition %in% "CHIMERYS", .(condition = factor("CHIMERYS"), .N), by = c("ptm_group_J")])
fwrite(data_venn_ptm_all_duell, file.path(figurePath, "figure-E2C-venn-peptides-all.csv"))
