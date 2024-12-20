#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-S2")
figurePath <- file.path(dataPath, "data/figure-S2")


#load CHIMERYS results
filePath <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys/20240517_lfq_dia_z1to4_v2x7x9_apex_True.pdResult")
sampleNames <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names.csv")
data_psm <- readPsms(filePath, sampleNames, loadBackup = F, writeBackup = F)

data_study <- readSampleNames(sampleNames)
data_psm[, file := factor(sample, data_study$sample_name, paste0(data_study$raw_file, ".raw"))]

data_scans <- data_psm[, .N, keyby=.(file, scan_ms2)]
data_scans[, is_identified := T]


#load raw file scans
rawFiles <- c("LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.raw",
              "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_02.raw",
              "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_02.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_03.raw")
rawPaths <- file.path(dataPath, "_external-raw-files/LFQ_Bench_multispecies", rawFiles)

if(any(!file.exists(rawPaths))) {
  cl <- makeForkCluster(length(rawPaths))
  registerDoParallel(cl)
  raw_scans <- foreach(rawPath = rawPaths, .combine = rbind, .packages = "rawrr") %dopar% {
    temp <- as.data.table(readIndex(rawPath))
    return(cbind(file = basename(rawPath), temp))
  }
  stopCluster(cl)
  fwrite(raw_scans, file = file.path(figurePath, "intermediate/raw_scans.csv"))
} else {
  files <- paste(basename(rawPaths)[!file.exists(rawPaths)], collapse = ", ")
  warning(paste(files, "from external PXD028735 not found, loading backup instead"))
  raw_scans <- fread(file.path(figurePath, "intermediate/raw_scans.csv"))
}

raw_scans[, file := factor(file, rawFiles)]
raw_scans <- raw_scans[MSOrder=="Ms2", ]


#scans identified
setkey(raw_scans, file, scan)
raw_scans <- data_scans[raw_scans]
raw_scans[is.na(is_identified), is_identified := F]
raw_scans[is.na(N), N := 0]

counts_scans <- raw_scans[, .N, keyby = is_identified]
counts_scans[, N_rel := N/sum(N)*100]
counts_scans[, N_label := paste0(format(N, big.mark = ","), "\n(", round(N_rel, 1), "%)")]
counts_scans[order(-N_label), ypos := cumsum(N_rel)-0.5*N_rel]

fwrite(counts_scans, file.path(figurePath, "figure-S2A-scans.csv"))


#PSMs per scan
scan_counts <- data_psm[, .(ms2_scans = uniqueN(psm_J)), by = .(sample, scan_ms2)][, .N, keyby = ms2_scans]
scan_counts[, N_label := format(N, big.mark=",")]

fwrite(scan_counts, file.path(figurePath, "figure-S2B-psms.csv"))
