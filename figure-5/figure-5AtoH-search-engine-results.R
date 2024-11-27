#setup
path <- file.path(here::here(), "figure-5")
source(here::here("scripts/load-dependencies.R"))

nPtm <- function(sequence, ptm = "+15.9949") {
  nSeq <- nchar(sequence)
  nSeqNoPtm <- nchar(gsub(ptm, "", sequence, fixed = T))
  as.integer((nSeq - nSeqNoPtm) / nchar(ptm))
}

read_zodiaq_protein <- function(data_path, sample_names_path, fasta_file_path) {
  fasta <- fread(fasta_file_path)
  fasta[, proteins := paste0("PROT", protNr)]
  fasta <- fasta[, .(proteins, protein)]
  setkey(fasta, proteins)

  filePaths <- list.files(data_path, "_fullOutput\\.csv$", full.names = T)
  zodiaq <- foreach(filePath = filePaths, .combine = rbind) %do% {
    #print(paste("Loading", basename(filePath)))
    tmp <- fread(filePath, fill=T)
    tmp[, file := tools::file_path_sans_ext(basename(fileName))]
    tmp[, is_decoy := as.logical(isDecoy)]
    tmp[, MaCC_Score := cosine * shared**0.2]
    #map proteins and delete shared peptides
    tmp[, proteins := gsub("^.*\\|(.*)\\|.*$", "\\1", protein)]
    tmp[, protein := NULL]
    setkey(tmp, proteins)
    tmp <- fasta[tmp]
    tmp <- tmp[!grepl(";", protein)]
    tmp[, proteins := protein]
    tmp[, protein := NULL]

    #kick out shared peptides, then collapse PSMs to proteins
    cols <- c("file", "proteins", "is_decoy", "MaCC_Score")
    tmp <- tmp[!grepl(";", proteins), .SD, .SDcols = cols]
    tmp <- tmp[tmp[order(-MaCC_Score),
                            .I[1], by = c("file", "proteins", "is_decoy")]$V1]
    tmp[is_decoy==TRUE, proteins := paste0("DECOY_", proteins)]
    setnames(tmp, "MaCC_Score", "score")

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
  #zodiaq[is_decoy==FALSE & q_value<=0.01, .N]
  #zodiaq[is_decoy==FALSE & q_value<=0.01, uniqueN(proteins)]

  #calculate dataset-global FDR (via best global candidate)
  fdr_global <- zodiaq[, .(score = max(score)), by=.(proteins, is_decoy)]
  setorder(fdr_global, -score)
  fdr_global[, q_value_global := cumsum(is_decoy) / cumsum(!is_decoy)]
  setorder(fdr_global, score)
  fdr_global[, q_value_global := cummin(q_value_global)]
  fdr_global[, c("is_decoy", "score") := NULL]
  setkey(fdr_global, proteins)
  setkey(zodiaq, proteins)
  zodiaq <- fdr_global[zodiaq]
  setkey(zodiaq, sample, proteins)
  #zodiaq[is_decoy==FALSE & q_value_global<=0.01, .N]
  #zodiaq[is_decoy==FALSE & q_value_global<=0.01, uniqueN(proteins)]

  #define conditions
  zodiaq[, replicate := as.integer(gsub("^.*-.*-.*-.*-.*-.*_(.)$", "\\1", sample))]
  zodiaq[, condition_faims := factor(gsub("^(.*)-.*-.*-.*-.*-.*$", "\\1", condition),
                                     c("faims", "noFaims"), c("FAIMS", "no FAIMS"))]
  zodiaq[, condition_study := factor(gsub("^.*-(.*)-.*-.*-.*-.*$", "\\1", condition),
                                     c("res", "it", "iw", "ov"))]
  zodiaq[, condition_resolution := as.integer(gsub("^.*-.*-(.*)k-.*-.*-.*$", "\\1", condition))]
  zodiaq[, condition_it := as.integer(gsub("^.*-.*-.*-(.*)ms-.*-.*$", "\\1", condition))]
  zodiaq[, condition_iw := as.integer(gsub("^.*-.*-.*-.*-(.*)th-.*$", "\\1", condition))]
  zodiaq[, condition_overlap := factor(gsub("^.*-.*-.*-.*-.*-(.*)$", "\\1", condition),
                                       c("0.5ov", "1ov", "3ov", "4ov", "6ov", "8ov"))]
  return(zodiaq[])
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
    tmp[, peptide := gsub(paste0("C", unimods["4"]), "C", peptide, fixed = T)]
    tmp[, nAc := nPtm(peptide, unimods["1"])]
    tmp[, nCa := nPtm(peptide, unimods["5"])]
    tmp[, nOx := nPtm(peptide, unimods["35"])]
    tmp[, peptide := gsub(unimods["1"], "", peptide, fixed = T)]
    tmp[, peptide := gsub(unimods["5"], "", peptide, fixed = T)]
    tmp[, peptide := gsub(unimods["35"], "", peptide, fixed = T)]
    tmp[, ptm_group := paste0(peptide, "_", nAc, "x1_", nCa, "x5_0x21_",
                              nOx, "x35_0x765_0x766")]
    tmp[, protein := gsub("^.*\\|(.*)\\|.*$", "\\1", protein)]
    #collapse onto ptm_group level
    tmp <- tmp[, .(score = max(MaCC_Score), peptide = peptide[1],
                   proteins = paste(sort(unique(protein)), collapse = ";")),
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
  #zodiaq[is_decoy==FALSE & q_value<=0.01, .N]
  #zodiaq[is_decoy==FALSE & q_value<=0.01, uniqueN(ptm_group)]

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
  #zodiaq[is_decoy==FALSE & q_value_global<=0.01, .N]
  #zodiaq[is_decoy==FALSE & q_value_global<=0.01, uniqueN(ptm_group)]

  #define conditions
  zodiaq[, replicate := as.integer(gsub("^.*-.*-.*-.*-.*-.*_(.)$", "\\1", sample))]
  zodiaq[, condition_faims := factor(gsub("^(.*)-.*-.*-.*-.*-.*$", "\\1", condition),
                                     c("faims", "noFaims"), c("FAIMS", "no FAIMS"))]
  zodiaq[, condition_study := factor(gsub("^.*-(.*)-.*-.*-.*-.*$", "\\1", condition),
                                     c("res", "it", "iw", "ov"))]
  zodiaq[, condition_resolution := as.integer(gsub("^.*-.*-(.*)k-.*-.*-.*$", "\\1", condition))]
  zodiaq[, condition_it := as.integer(gsub("^.*-.*-.*-(.*)ms-.*-.*$", "\\1", condition))]
  zodiaq[, condition_iw := as.integer(gsub("^.*-.*-.*-.*-(.*)th-.*$", "\\1", condition))]
  zodiaq[, condition_overlap := factor(gsub("^.*-.*-.*-.*-.*-(.*)$", "\\1", condition),
                                       c("0.5ov", "1ov", "3ov", "4ov", "6ov", "8ov"))]
  return(zodiaq[])
}

#Load CHIMERYS search results
filePath <- file.path(dataPath, "figure-5/CHIMERYS-DISPA-mcf7_MS2_2024-10-14.pdResult")
sampleNamesPath <- suppressMessages(writeSampleNames(filePath, outputPath = path))
#load PSMs
diPsms <- suppressMessages(readPsms(filePath, sampleNamesPath, loadBackup = T,
                                    loadDecoys = T, filterFdr = F))
diPsms[, score_coefficient_lasso := log2(score_coefficient_lasso)]
setnames(diPsms, "score_svm", "score")
parseModificationsColumn(dataTable = diPsms)

#collapse PSMs to peptide groups
diPtm <- diPsms[diPsms[order(-score, -score_coefficient_lasso),
                       .I[1], by = c("sample", "ptm_group", "is_decoy")]$V1]
cols <- c("sample", "condition", "replicate", "ptm_group", "is_decoy", "score", "score_coefficient_lasso")
diPtm <- diPtm[, .SD, .SDcols = cols]

#file-local FDR (monotonized)
setorder(diPtm, sample, -score)
diPtm[, q_value := cumsum(is_decoy) / cumsum(!is_decoy), by = sample]
setorder(diPtm, sample, score)
diPtm[, q_value := cummin(q_value), by = sample]
diPtm[is_decoy==FALSE & q_value<=0.01, .N] #70204
diPtm[is_decoy==FALSE & q_value<=0.01, uniqueN(ptm_group)] #4335

#data set-global FDR (via best global candidate, monotonized)
fdrGlobal <- diPtm[, .(score = max(score)), by=.(ptm_group, is_decoy)]
setorder(fdrGlobal, -score)
fdrGlobal[, q_value_global := cumsum(is_decoy) / cumsum(!is_decoy)]
setorder(fdrGlobal, score)
fdrGlobal[, q_value_global := cummin(q_value_global)]
fdrGlobal[, c("is_decoy", "score") := NULL]
setkey(fdrGlobal, ptm_group)
setkey(diPtm, ptm_group)
diPtm <- fdrGlobal[diPtm]
setkey(diPtm, sample, ptm_group)
diPtm[is_decoy==FALSE & q_value_global<=0.01, .N] #100992
diPtm[is_decoy==FALSE & q_value_global<=0.01, uniqueN(ptm_group)] #2733

#load peptide groups
diPtm[, condition_faims := factor(gsub("^(.*)-.*-.*-.*-.*-.*$", "\\1", condition),
                                  c("faims", "noFaims"), c("FAIMS", "no FAIMS"))]
diPtm[, condition_study := factor(gsub("^.*-(.*)-.*-.*-.*-.*$", "\\1", condition),
                                  c("res", "it", "iw", "ov"))]
diPtm[, condition_resolution := as.integer(gsub("^.*-.*-(.*)k-.*-.*-.*$", "\\1", condition))]
diPtm[, condition_it := as.integer(gsub("^.*-.*-.*-(.*)ms-.*-.*$", "\\1", condition))]
diPtm[, condition_iw := as.integer(gsub("^.*-.*-.*-.*-(.*)th-.*$", "\\1", condition))]
diPtm[, condition_overlap := factor(gsub("^.*-.*-.*-.*-.*-(.*)$", "\\1", condition),
                                    c("0.5ov", "1ov", "3ov", "4ov", "6ov", "8ov"))]
write_fst(diPtm, file.path(dataPath, "figure-5/diPtm.fst"))

#kick out shared peptides, then collapse PSMs to proteins
#diMap <- diPsms[, .(protein = unlist(strsplit(proteins, ";"))), keyby=.(sample, psm)]
cols <- c("sample", "condition", "replicate", "proteins", "is_decoy", "score", "score_coefficient_lasso")
diProt <- diPsms[!grepl(";", proteins), .SD, .SDcols = cols]
diProt <- diProt[diProt[order(-score, -score_coefficient_lasso),
                        .I[1], by = c("sample", "proteins", "is_decoy")]$V1]
diProt[is_decoy==TRUE, proteins := paste0("DECOY_", proteins)]

#file-local FDR (monotonized)
setorder(diProt, sample, -score)
diProt[, q_value := cumsum(is_decoy) / cumsum(!is_decoy), by = sample]
setorder(diProt, sample, score)
diProt[, q_value := cummin(q_value), by = sample]
diProt[is_decoy==FALSE & q_value<=0.01, .N] #15457
diProt[is_decoy==FALSE & q_value<=0.01, uniqueN(proteins)] #978

#data set-global FDR (via best global candidate, monotonized)
fdrGlobal <- diProt[, .(score = max(score)), by=.(proteins, is_decoy)]
setorder(fdrGlobal, -score)
fdrGlobal[, q_value_global := cumsum(is_decoy) / cumsum(!is_decoy)]
setorder(fdrGlobal, score)
fdrGlobal[, q_value_global := cummin(q_value_global)]
fdrGlobal[, c("is_decoy", "score") := NULL]
setkey(fdrGlobal, proteins)
setkey(diProt, proteins)
diProt <- fdrGlobal[diProt]
setkey(diProt, sample, proteins)
diProt[is_decoy==FALSE & q_value_global<=0.01, .N] #26361
diProt[is_decoy==FALSE & q_value_global<=0.01, uniqueN(proteins)] #468

#load peptide groups
diProt[, condition_faims := factor(gsub("^(.*)-.*-.*-.*-.*-.*$", "\\1", condition),
                                   c("faims", "noFaims"), c("FAIMS", "no FAIMS"))]
diProt[, condition_study := factor(gsub("^.*-(.*)-.*-.*-.*-.*$", "\\1", condition),
                                   c("res", "it", "iw", "ov"))]
diProt[, condition_resolution := as.integer(gsub("^.*-.*-(.*)k-.*-.*-.*$", "\\1", condition))]
diProt[, condition_it := as.integer(gsub("^.*-.*-.*-(.*)ms-.*-.*$", "\\1", condition))]
diProt[, condition_iw := as.integer(gsub("^.*-.*-.*-.*-(.*)th-.*$", "\\1", condition))]
diProt[, condition_overlap := factor(gsub("^.*-.*-.*-.*-.*-(.*)$", "\\1", condition),
                                     c("0.5ov", "1ov", "3ov", "4ov", "6ov", "8ov"))]
write_fst(diProt, file.path(dataPath, "figure-5/diProt.fst"))


#Load csoDIAq results CHIMERYS digest
unimods <- c("1" = "[UNIMOD:1]", "4" = "[UNIMOD:4]", "5" = "[UNIMOD:5]", "35" = "[UNIMOD:35]")
csodiaq_digest <-
  read_zodiaq(data_path = file.path(dataPath, "figure-5/ZO_2024-10-17_inferys_filtered_noNL_morethan400mz"),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))
write_fst(csodiaq_digest, file.path(dataPath, "figure-5/csodiaq_digest.fst"))

#Load csoDIAq results CHIMERYS digest
csodiaq_digest_prot <-
  read_zodiaq_protein(data_path = file.path(dataPath, "figure-5/ZO_2024-10-17_inferys_filtered_noNL_morethan400mz"),
                      sample_names_path = file.path(path, "sample_names.csv"),
                      fasta_file_path = file.path(dataPath, "figure-5/fasta_mapping.csv"))
write_fst(csodiaq_digest_prot, file.path(dataPath, "figure-5/csodiaq_digest_prot.fst"))
