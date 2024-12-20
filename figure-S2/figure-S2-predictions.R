#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-S2")
figurePath <- file.path(dataPath, "data/figure-S2")

filePath <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys/20240517_lfq_dia_z1to4_v2x7x9_apex_True.pdResult")
sampleNames <- writeSampleNames(filePath, outputPath = path, outputName = "sample_names.csv")
data_psm <- readPsms(filePath, sampleNames, loadBackup = F, writeBackup = F)

ls_psms <- split(data_psm, by="sample")

# ---- merge function ----
rawFiles <- c("LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.raw",
              "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_02.raw",
              "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_02.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_03.raw")
rawPaths <- file.path(dataPath, "_external-raw-files/LFQ_Bench_multispecies", rawFiles)


merged <- foreach(ls_psm = ls_psms, rawPath = rawPaths, iRaw = seq_along(rawFiles), .combine = rbind) %do% {
  print(paste("Starting", iRaw))

  #peptide predictions
  count_psm <- ls_psm[, .(n_PSM = uniqueN(psm_J)), keyby = scan_ms2]
  peptides <- ls_psm[scan_ms2 %in% count_psm$scan_ms2,
                     .(scan_ms2, mz_ratio, ptm, charge, nce,
                       score_coefficient_normalized, score_coefficient_lasso)]
  peptides[, ptm := gsub("[unimod:", "[UNIMOD:", ptm, fixed = T)]
  setkey(peptides, scan_ms2, ptm)
  #remove IL ambiguous peptides and re-normalize normalized CHIMERYS coefficient
  peptides[, ptm_J := gsub("(?<!\\[..)(I|L)", "J", ptm, perl=T)]
  peptides <- peptides[, .SD[1], by=.(scan_ms2, ptm_J, charge)]
  peptides[, score_coefficient_normalized := score_coefficient_normalized/sum(score_coefficient_normalized), by=scan_ms2]
  peptides[, id := 1L:.N]
  setkey(peptides, id)

  #INFERYS predictions (requires local INFERYS API and internal package)
  # rawActivation <- c("HCD", "CID")[1]
  # inferysApi <- "localhost:8081"
  # inferysModel <- "inferys_3.0.0_fragmentation"
  # predictions <- predictSpectrum(sequences = peptides$ptm,
  #                                charges = peptides$charge,
  #                                collisionEnergies = peptides$nce,
  #                                activationType = rawActivation,
  #                                model = inferysModel,
  #                                spectrumApiIP = inferysApi,
  #                                simplify = T)
  # predictions <- rbindlist(predictions, idcol = "id")
  # setkey(predictions, id)
  # setnames(predictions, c("id", "annotation", "charge_ms2", "intensity_ms2", "mz_ms2"))
  # predictions <- peptides[predictions]
  # predictions[, ion := gsub("^(.).*$", "\\1", annotation)]
  # predictions[, position := as.integer(gsub("^.(\\d*)($|-.*$)", "\\1", annotation))]
  # setkey(predictions, scan_ms2)
  # #renormalize fragment ion intensities to sum 1 per id
  # predictions[, intensity_normalized_ms2 := intensity_ms2/sum(intensity_ms2), by=id]
  # write_fst(predictions, file.path(figurePath, "intermediate", paste0(iRaw, "-predictions.fst")), compress = 100)
  predictions <- read_fst(file.path(figurePath, "intermediate", paste0(iRaw, "-predictions.fst")), as.data.table = T)

  #raw data
  #extract spectra (adjust number of cores to those available on your machine)
  # rawScans <- peptides[, unique(scan_ms2)]
  # peaks <- foreach(scans = split(rawScans, ceiling(seq_along(rawScans)/500)), .combine = rbind) %do% {
  #   spectra <- readSpectrum(rawPath, scans)
  #   peaks <- rbindlist(lapply(spectra, function(x) x[c("scan", "centroid.mZ", "centroid.intensity")] ))
  #   setnames(peaks, c("scan", "mz", "intensity"))
  #   return(peaks)
  # }
  # write_fst(peaks, file.path(figurePath, "intermediate", paste0(iRaw, "-peaks.fst")), compress = 100)
  peaks <- read_fst(file.path(figurePath, "intermediate", paste0(iRaw, "-peaks.fst")), as.data.table = T)

  #merge
  # merged <- mergeSpectra(predictions[, .(scan_ms2, mz_ms2)],
  #                        peaks[, .(scan, mz, intensity)],
  #                        mzTruthName = "mz_ms2",
  #                        scanObservedName = "scan")
  # setnames(merged, c("scan", "mz"), c("scan_ms2", "mz_ms2"))
  # #expand
  # setkey(predictions, scan_ms2, mz_ms2)
  # merged <- merged[predictions]
  # merged[, intensity_lasso_ms2 := intensity_normalized_ms2 * score_coefficient_lasso]
  # merged[, intensity_lasso_scaled_ms2 := score_coefficient_lasso/sum(score_coefficient_lasso), by=position]
  # merged[, intensity_scaled_ms2 := intensity_normalized_ms2 * score_coefficient_normalized / uniqueN(scan_ms2)]
  # merged[, n_pos := .N, by = position]
  # merged[!is.na(mzMatch), n_pos_matched := .N, by = position]
  # merged[, n_ptm_pred := uniqueN(ptm), by=.(scan_ms2, mz_ms2)]
  # merged[, is_ptm_pred_shared := n_ptm_pred>1]
  # merged[, mz_ms2_200 := factor(ceiling(mz_ms2/200)*200)]
  # merged[, mz_ms2_300 := factor(ceiling(mz_ms2/300)*300)]
  # merged[!is.na(mzMatch), n_ptm_match := uniqueN(ptm), by=.(scan_ms2, mzMatch)]
  # merged[!is.na(mzMatch), is_ptm_match_shared := n_ptm_match>1]
  # merged[!is.na(mzMatch), mzMatch_10 := factor(ceiling(mzMatch/10)*10)]
  # merged[!is.na(mzMatch), mzMatch_200 := factor(ceiling(mzMatch/200)*200)]
  # merged[!is.na(mzMatch), mzMatch_300 := factor(ceiling(mzMatch/300)*300)]
  # merged[!is.na(mzMatch), .(mzMatch, mzMatch_10, mzMatch_200, mzMatch_300)]
  # write_fst(merged, file.path(figurePath, "intermediate", paste0(iRaw, "-merged.fst")), compress = 100)
  merged <- read_fst(file.path(figurePath, "intermediate", paste0(iRaw, "-merged.fst")), as.data.table = T)
  return(cbind(file = factor(basename(rawPath)), merged))
}
write_fst(merged, file.path(figurePath, "intermediate/merged.fst"), compress = 100)


# ---- reload data ----
# merged <- read_fst(file.path(figurePath, "intermediate/merged.fst"), as.data.table = T)

# ---- predicted fragments ----
merged[n_ptm_pred>1, .(.N,
                       N_rel = round(.N/n_pos[1]*100, 2),
                       int_lasso = scientific(sum(intensity_lasso_ms2), 1),
                       int_lasso_rel = round(sum(intensity_lasso_scaled_ms2)*100, 2)),
       keyby=position]

#plot positions
count_pos <- merged[, .N, keyby=.(position, is_ptm_pred_shared)]
count_pos[, N_rel := N/sum(N), by=position]
count_pos[, N_rel_label := paste0(round(N_rel*100, 0), "%")]
fwrite(count_pos, file.path(figurePath, "figure-S2CE-predicted-position.csv"))


#plot 200 mz-bins
count_200 <- merged[, .N, keyby=.(mz_ms2_200, is_ptm_pred_shared)]
mz_lab <- count_200[, paste0(c(0, 0, paste0(">", as.character(mz_ms2_200)[1:(.N-2)])), "-\n", mz_ms2_200)]
count_200[, mzMatch_label := factor(mz_lab, unique(mz_lab))]
count_200[, N_rel := N/sum(N), by=mz_ms2_200]
count_200[, N_rel_label := paste0(round(N_rel*100, 2), "%")]
fwrite(count_200, file.path(figurePath, "figure-S2DF-predicted-mz.csv"))


# ---- matched fragments ----
merged[n_ptm_match>1, .(.N,
                        N_rel = round(.N/n_pos_matched[1]*100, 2),
                        int_lasso = scientific(sum(intensity_lasso_ms2), 1),
                        int_lasso_rel = round(sum(intensity_lasso_scaled_ms2)*100, 2)),
       keyby=position]


#plot positions
count_pos <- merged[mzMatched==T, .N, keyby=.(position, is_ptm_match_shared)]
count_pos[, N_rel := N/sum(N), by=position]
count_pos[, N_rel_label := paste0(round(N_rel*100, 0), "%")]
fwrite(count_pos, file.path(figurePath, "figure-S1GI-matched-position.csv"))


#plot 200 mz-bins
count_200 <- merged[mzMatched==T, .N, keyby=.(mzMatch_200, is_ptm_match_shared)]
mz_lab <- count_200[, paste0(c(0, 0, paste0(">", as.character(mzMatch_200)[1:(.N-2)])), "-\n", mzMatch_200)]
count_200[, mzMatch_label := factor(mz_lab, unique(mz_lab))]
count_200[, N_rel := N/sum(N), by=mzMatch_200]
count_200[, N_rel_label := paste0(round(N_rel*100, 2), "%")]
fwrite(count_200, file.path(figurePath, "figure-S1HJ-matched-mz.csv"))
