#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-1")
deconPath <- file.path(dataPath, "figure-1/decon-mirror-upset")

##load PSMs, precursors and peptide groups
#`dataPath` needs to point to the folder containing search engine research result folders
sampleNamesPath <- writeSampleNames(deconPath, outputPath = path, outputName = "decon-mirror-upset-sample-names2.csv")
data_psm <- readPsms(deconPath, sampleNamesPath)
data_pcm <- readPcms(deconPath, sampleNamesPath)
data_ptm <- readPtmGroups(deconPath, sampleNamesPath)

##deconvolution mirror plot
rawScans <- 115649L
data_sub <- data_psm[sample %in% "CHIMERYS_1" & scan_ms2 %in% rawScans, .SD,
                     .SDcols = c("sample", "scan_ms2", "mz_ratio", "ptm", "charge", "nce", "score_coefficient_normalized")]
#export file for plotting
fwrite(data_sub, file.path(dirname(deconPath), "decon-mirror.csv"))


##upset plot
#converge PSMs to PTM groups
parseModificationsColumn(dataTable = data_psm)
data_psm_ptm <- data_psm[condition %in% c("Metamorpheus", "MSGFplus"),
                         .N, keyby = c("condition", "sample", "ptm_group_J")]
count_psm_ptm <- data_psm_ptm[, .(level = "PSM-level FDR", .N), keyby = c("condition")]

#converge PCMs to PTM groups
parseModificationsColumn(dataTable = data_pcm)
data_pcm_ptm <- data_pcm[condition %in% c("MSFragger", "MSFragger-DDAplus"),
                         .N, keyby = c("condition", "sample", "ptm_group_J")]
count_pcm_ptm <- data_pcm_ptm[, .(level = "PCM-level FDR", .N), keyby = c("condition")]


#subset data for upset plot
data_upset <- rbind(cbind(data_psm_ptm, level = "PSM-level FDR"),
                    cbind(data_pcm_ptm, level = "PCM-level FDR"),
                    data_ptm[, .(.N, level = "PTM group-level FDR"),
                             by = c("condition", "sample", "ptm_group_J")])
data_upset <- data_upset[!condition %in% c("Metamorpheus-Percolator", "MSFragger-full-isolation")]

fwrite(data_upset, file.path(dirname(deconPath), "upset.csv"))
