#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-E9")
figurePath <- file.path(dataPath, "data/figure-E9")

nPtm <- function(sequence, ptm = "+15.9949") {
  nSeq <- nchar(sequence)
  nSeqNoPtm <- nchar(gsub(ptm, "", sequence, fixed = T))
  as.integer((nSeq - nSeqNoPtm) / nchar(ptm))
}

read_zodiaq <- function(data_path, unimods, sample_names_path) {
  stopifnot(all(c("1", "4", "5", "35") %in% names(unimods)))
  filePaths <- list.files(data_path, "_fullOutput\\.csv$", full.names = T)
  zodiaq <- foreach(filePath = filePaths, .combine = rbind) %do% {
    #print(paste("Loading", basename(filePath)))
    tmp <- fread(filePath, fill=T)
    tmp[, file := tools::file_path_sans_ext(basename(fileName))]
    tmp[, is_decoy := as.logical(isDecoy)]
    tmp[, MaCC_Score := cosine * shared**0.2]
    #create peptide groups
    #unimod 5 is likely misannotated in library, see section above, count as Ac
    tmp[, peptide := gsub(paste0("C", unimods["4"]), "C", peptide, fixed = T)]
    tmp[, nAc := nPtm(peptide, unimods["1"])]
    tmp[, nCa := nPtm(peptide, unimods["5"])]
    tmp[, nOx := nPtm(peptide, unimods["35"])]
    tmp[, peptide := gsub(unimods["1"], "", peptide, fixed = T)]
    tmp[, peptide := gsub(unimods["5"], "", peptide, fixed = T)]
    tmp[, peptide := gsub(unimods["35"], "", peptide, fixed = T)]
    tmp[, ptm_group := paste0(peptide, "_", nAc, "x1_", nCa, "x5_0x21_",
                              nOx, "x35_0x765_0x766")]
    #collapse onto ptm_group level
    tmp <- tmp[, .(score = max(MaCC_Score), peptide = peptide[1]),
               by=.(file, ptm_group, is_decoy)]
    return(tmp)
  }
  setkey(zodiaq, file)
  #map sample names
  dtSamples <- readSampleNames(sample_names_path)
  dtSamples <- dtSamples[, .(file = raw_file, sample = sample_name, condition = condition_name)]
  setkey(dtSamples, file)
  zodiaq <- dtSamples[zodiaq]

  #calculate file-local FDR
  setorder(zodiaq, sample, -score)
  zodiaq[, q_value := cumsum(is_decoy) / cumsum(!is_decoy), by=sample]
  setorder(zodiaq, sample, score)
  zodiaq[, q_value := cummin(q_value), by = sample]
  zodiaq[is_decoy==FALSE & q_value<=0.01, .N]
  zodiaq[is_decoy==FALSE & q_value<=0.01, uniqueN(ptm_group)]

  #calculate dataset-global FDR (via best global candidate)
  fdr_global <- zodiaq[, .(score = max(score)), by=.(ptm_group, is_decoy)]
  setorder(fdr_global, -score)
  fdr_global[, q_value_global := cumsum(is_decoy) / cumsum(!is_decoy)]
  setorder(fdr_global, score)
  fdr_global[, q_value_global := cummin(q_value_global)]
  fdr_global[, c("is_decoy", "score") := NULL]
  setkey(fdr_global, ptm_group)
  setkey(zodiaq, ptm_group)
  zodiaq <- fdr_global[zodiaq]
  setkey(zodiaq, sample, ptm_group)
  zodiaq[is_decoy==FALSE & q_value_global<=0.01, .N]
  zodiaq[is_decoy==FALSE & q_value_global<=0.01, uniqueN(ptm_group)]

  #define conditions
  zodiaq[, replicate := as.integer(gsub("^.*-.*-.*-.*-.*-.*_(.)$", "\\1", sample))]
  zodiaq[, condition_faims := factor(gsub("^(.*)-.*-.*-.*-.*-.*$", "\\1", condition),
                                     c("faims", "noFaims"))]
  zodiaq[, condition_study := factor(gsub("^.*-(.*)-.*-.*-.*-.*$", "\\1", condition),
                                     c("res", "it", "iw", "ov"))]
  zodiaq[, condition_resolution := as.integer(gsub("^.*-.*-(.*)k-.*-.*-.*$", "\\1", condition))]
  zodiaq[, condition_it := as.integer(gsub("^.*-.*-.*-(.*)ms-.*-.*$", "\\1", condition))]
  zodiaq[, condition_iw := as.integer(gsub("^.*-.*-.*-.*-(.*)th-.*$", "\\1", condition))]
  zodiaq[, condition_overlap := factor(gsub("^.*-.*-.*-.*-.*-(.*)$", "\\1", condition),
                                       c("0.5ov", "1ov", "3ov", "4ov", "6ov", "8ov"))]
  return(zodiaq[])
}

#Load libraries for peptide length distributions
lib_zodiaq <- fread(file.path(dataPath, "Direct_infusion/human_lib.tsv")) #3108165
#fix decoy column (UniprotID is sorted, if first entry target then not a decoy)
lib_zodiaq[, isDecoy := grepl("^\\d{1,2}/DECOY", UniprotID)]
lib_zodiaq[, isDecoyLabel := factor(ifelse(isDecoy, "Decoy", "Target"), c("Target", "Decoy"))]
pep_zodiaq <- lib_zodiaq[, .(n_aa = nchar(PeptideSequence[1])),
                         by=.(isDecoyLabel, PeptideGroupLabel)] #273373
fwrite(pep_zodiaq, file.path(figurePath, "figure-E9A-csodiaq.csv"))

lib_inferys <- fread(file.path(dataPath, "Direct_infusion/human_lib_inferys_filtered.tsv")) #2599104
lib_inferys[, isDecoyLabel := factor(ifelse(isDecoy, "Decoy", "Target"), c("Target", "Decoy"))]
pep_inferys <- lib_inferys[, .(n_aa = nchar(PeptideSequence[1])),
                           by=.(isDecoyLabel, PeptideGroupLabel)] #259916
fwrite(pep_inferys, file.path(figurePath, "figure-E9B-chimerys.csv"))


#Load search engine counts
search_results_folder <- file.path(dataPath, "Direct_infusion/CsoDIAq")
#search_results_folder <- "/mnt/storage1/0_Users/Alex/Projects/direct-infusion/2024-09_i81_Meyer-2020-DISPA"

#complete csoDIAq
unimods <- c("1" = "(UniMod:1)", "4" = "(UniMod:4)", "5" = "(UniMod:5)", "35" = "(UniMod:35)")
study_path <- "ZO_2024-09-27_zodiaq-mcf7-tsv-library"
zodiaq_original <-
  read_zodiaq(data_path = file.path(search_results_folder, study_path),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))

#Filtered to overlap
unimods <- c("1" = "(UniMod:1)", "4" = "(UniMod:4)", "5" = "(UniMod:5)", "35" = "(UniMod:35)")
study_path <- "ZO_2024-10-14_zodiaq-mcf7-tsv-library-filtered-to-inferys"
zodiaq_overlap <-
  read_zodiaq(data_path = file.path(search_results_folder, study_path),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))

#Targets INFERYS
unimods <- c("1" = "[UNIMOD:1]", "4" = "[UNIMOD:4]", "5" = "[UNIMOD:5]", "35" = "[UNIMOD:35]")
study_path <- "ZO_2024-10-17_inferys_targets_zodiaq_decoys_noNL_morethan400mz"
zodiaq_T_inferys <-
  read_zodiaq(data_path = file.path(search_results_folder, study_path),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))

#Decoys INFERYS
unimods <- c("1" = "[UNIMOD:1]", "4" = "[UNIMOD:4]", "5" = "[UNIMOD:5]", "35" = "[UNIMOD:35]")
study_path <- "ZO_2024-10-17_zodiaq_targets_inferys_decoys_noNL_morethan400mz"
zodiaq_D_inferys <-
  read_zodiaq(data_path = file.path(search_results_folder, study_path),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))

#Both INFERYS
unimods <- c("1" = "[UNIMOD:1]", "4" = "[UNIMOD:4]", "5" = "[UNIMOD:5]", "35" = "[UNIMOD:35]")
study_path <- "ZO_2024-10-17_inferys_targets_inferys_decoys_noNL_morethan400mz"
zodiaq_both_inferys <-
  read_zodiaq(data_path = file.path(search_results_folder, study_path),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))

#CHIMERYS digest
unimods <- c("1" = "[UNIMOD:1]", "4" = "[UNIMOD:4]", "5" = "[UNIMOD:5]", "35" = "[UNIMOD:35]")
study_path <- "ZO_2024-10-17_inferys_filtered_noNL_morethan400mz"
chimerys_digest <-
  read_zodiaq(data_path = file.path(search_results_folder, study_path),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))

types <- c("CsoDIAq original",
           "CsoDIAq overlap",
           "CsoDIAq overlap\nTargets INFERYS",
           "CsoDIAq overlap\nDecoys INFERYS",
           "CsoDIAq overlap\nBoth INFERYS",
           "CHIMERYS digest\nBoth INFERYS")
results <-
  rbind(cbind(type = types[1L], zodiaq_original),
        cbind(type = types[2L], zodiaq_overlap),
        cbind(type = types[3L], zodiaq_T_inferys),
        cbind(type = types[4L], zodiaq_D_inferys),
        cbind(type = types[5L], zodiaq_both_inferys),
        cbind(type = types[6L], chimerys_digest))
results[, type := factor(type, types)]
results_local <- results[is_decoy==FALSE & q_value<=0.01]
count_local <- results_local[, .N, keyby=.(type, sample)]

fwrite(count_local, file.path(figurePath, "figure-E9C-counts.csv"))
