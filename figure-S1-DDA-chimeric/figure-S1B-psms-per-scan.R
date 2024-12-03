#setup
source(here::here("scripts/load-dependencies.R"))
path <- paste0(here::here(), "/figure-S1-DDA-chimeric")
figurePath <- file.path(dataPath, "data/figure-S1")

filePath <- file.path(dataPath, "LFQ_Bench_human/Chimerys/LFQ_01_CHIMERYS_v2x7x9_apex_True.pdResult")
sampleNamesPath <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names.csv")

#PSMs
data_psm <- readPsms(filePath, sampleNamesPath)
scan_counts <- data_psm[sample=="CHIMERYS_1", .(ms2_scans = uniqueN(psm_J)), by = scan_ms2][, .N, keyby = ms2_scans]
scan_counts[, N_label := format(N, big.mark=",")]

fwrite(scan_counts, file.path(figurePath, "psms-per-scan.csv"))
