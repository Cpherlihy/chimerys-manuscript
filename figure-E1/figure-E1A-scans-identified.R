#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-E1")
figurePath <- file.path(dataPath, "data/figure-E1")

filePath <- file.path(dataPath, "LFQ_Bench_human/Chimerys/LFQ_01_CHIMERYS_v2x7x9_apex_True.pdResult")
sampleNamesPath <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names.csv")

#PSMs
data_psm <- readPsms(filePath, sampleNamesPath, loadBackup = F, writeBackup = F)
data_scans <- data_psm[, .(file = factor("DDA"), scan_ms2)][, .N, keyby=.(file, scan_ms2)]
data_scans[, is_identified := T]

#raw file
rawPath <- file.path(dataPath, "_external-raw-files/LFQ_Bench_human/LFQ_Orbitrap_DDA_Human_01.raw")
if(file.exists(rawPath)) {
  raw_scans <- as.data.table(readIndex(rawPath))
  # fwrite(raw_scans, file.path(figurePath, "intermediate/raw_scans.csv"))
} else {
  warning("'LFQ_Orbitrap_DDA_Human_01.raw' from external PXD028735 not found, loading backup instead")
  raw_scans <- fread(file.path(figurePath, "intermediate/raw_scans.csv"))
}
raw_scans[, file := factor("DDA")]
raw_scans <- raw_scans[MSOrder=="Ms2", ]

setkey(raw_scans, file, scan)
raw_scans <- data_scans[raw_scans]
raw_scans[is.na(is_identified), is_identified := F]
raw_scans[is.na(N), N := 0]

counts_scans <- raw_scans[, .N, keyby = is_identified]
counts_scans[, N_rel := N/sum(N)*100]
counts_scans[, N_label := paste0(format(N, big.mark = ","), "\n(", round(N_rel, 1), "%)")]
counts_scans[order(-N_label), ypos := cumsum(N_rel)-0.5*N_rel]

fwrite(counts_scans, file.path(figurePath, "figure-E1A-scans.csv"))
