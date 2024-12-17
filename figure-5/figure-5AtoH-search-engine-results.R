#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-5")
figurePath <- file.path(dataPath, "data/figure-5")

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
filePath <- file.path(dataPath, "Direct_infusion/CHIMERYS-DISPA-mcf7_MS2_2024-10-14.pdResult")
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
write_fst(diPtm, file.path(figurePath, "intermediate/diPtm.fst"))

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
write_fst(diProt, file.path(figurePath, "intermediate/diProt.fst"))


#Load csoDIAq results CHIMERYS digest
unimods <- c("1" = "[UNIMOD:1]", "4" = "[UNIMOD:4]", "5" = "[UNIMOD:5]", "35" = "[UNIMOD:35]")
csodiaq_digest <-
  read_zodiaq(data_path = file.path(dataPath, "Direct_infusion/ZO_2024-10-17_inferys_filtered_noNL_morethan400mz"),
              unimods = unimods, sample_names_path = file.path(path, "sample_names.csv"))
write_fst(csodiaq_digest, file.path(figurePath, "intermediate/csodiaq_digest.fst"))

#Load csoDIAq results CHIMERYS digest
csodiaq_digest_prot <-
  read_zodiaq_protein(data_path = file.path(dataPath, "Direct_infusion/ZO_2024-10-17_inferys_filtered_noNL_morethan400mz"),
                      sample_names_path = file.path(path, "sample_names.csv"),
                      fasta_file_path = file.path(figurePath, "intermediate/fasta_mapping.csv"))
write_fst(csodiaq_digest_prot, file.path(figurePath, "intermediate/csodiaq_digest_prot.fst"))


## Data processing ====
#peptide groups
diPtm <- read_fst(file.path(figurePath, "intermediate/diPtm.fst"), as.data.table = T)
csodiaq_digest <- read_fst(file.path(figurePath, "intermediate/csodiaq_digest.fst"), as.data.table = T)

namesBoth <- intersect(names(diPtm), names(csodiaq_digest))
namesChim <- c(namesBoth, "score_coefficient_lasso")
results <- rbind(cbind(type = "CHIMERYS", diPtm[, .SD, .SDcols = namesChim]),
                 cbind(type = "CsoDIAq", csodiaq_digest[, .SD, .SDcols = namesBoth]),
                 fill = TRUE)
shapeBoth <- c("no FAIMS" = "dotted", "FAIMS" = "solid")
colBoth <- c("CHIMERYS" = msaid_blue, "CsoDIAq" = msaid_green)
results[, type := factor(type, names(colBoth))]
results_local <- results[is_decoy==FALSE & q_value<=0.01]
results_local[, nType := .N, by=.(sample, ptm_group)]
results_local[, nTypeLabel := factor(nType, 1L:2L, c("CHIMERYS", "Shared"))]

byCol <- c("type", "condition", "replicate", "sample", "condition_faims", "condition_iw",
           "condition_study", "condition_it", "condition_resolution")
counts_local <- results_local[, .N, keyby = byCol]
fwrite(counts_local, file.path(figurePath, "figure-5ABC-counts.csv"))

rel_local <- counts_local[, .("CHIMERYS vs CsoDIAq" = N[type=="CHIMERYS"]/N[type=="CsoDIAq"]), by=sample]
rel_local <- melt(rel_local, id.vars = "sample", variable.name = "type", value.name = "rel")
fwrite(rel_local, file.path(figurePath, "figure-5D-relative.csv"))

results_sub <- results_local[type=="CHIMERYS" & sample=="faims-it-60k-20ms-1th-0.5ov_1" &
                               !is.na(score_coefficient_lasso)]
fwrite(results_sub, file.path(figurePath, "figure-5F-sensitivity.csv"))

#proteins
diProt <- read_fst(file.path(figurePath, "intermediate/diProt.fst"), as.data.table = T)
csodiaq_digest_prot <- read_fst(file.path(figurePath, "intermediate/csodiaq_digest_prot.fst"), as.data.table = T)

namesBoth <- intersect(names(diProt), names(csodiaq_digest_prot))
namesChim <- c(namesBoth, "score_coefficient_lasso")
proteins <- rbind(cbind(type = "CHIMERYS", diProt[, .SD, .SDcols = namesChim]),
                  cbind(type = "CsoDIAq", csodiaq_digest_prot[, .SD, .SDcols = namesBoth]),
                  fill = TRUE)
proteins[, proteins := gsub("Cont_", "", proteins, fixed = T)]
proteins[, type := factor(type, names(colBoth))]
proteins_local <- proteins[is_decoy==FALSE & q_value<=0.01]
proteins_local[, nType := .N, by=.(sample, proteins)]
proteins_local[nType==2, nTypeLabel := "Shared"]
proteins_local[nType==1 & type=="CHIMERYS", nTypeLabel := "CHIMERYS"]
proteins_local[nType==1 & type=="CsoDIAq", nTypeLabel := "CsoDIAq"]
proteins_local[, nTypeLabel := factor(nTypeLabel, c("CHIMERYS", "CsoDIAq", "Shared"))]

#rank
proteins_sub <- proteins_local[type=="CHIMERYS" & sample=="faims-it-60k-20ms-1th-0.5ov_1" & !is.na(score_coefficient_lasso)]
setorder(proteins_sub, -score_coefficient_lasso)
proteins_sub[, rank := 1L:.N]
fwrite(proteins_sub, file.path(figurePath, "figure-5G-rank.csv"))

#proteins sub
proteins_sub2 <- proteins_local[sample=="faims-it-60k-20ms-1th-0.5ov_1"]
#proteins_sub2[, nTypeLabel := factor(nTypeLabel, c("CHIMERYS", "CsoDIAq", "Shared"))]
proteins_sub2_count <- proteins_sub2[, .N, keyby=.(type, nTypeLabel)]
fwrite(proteins_sub2_count, file.path(figurePath, "figure-5E-protein.csv"))

colShared <- c("CHIMERYS" = msaid_blue, "CsoDIAq" = msaid_green, "Shared" = msaid_orange)
colPathway <- c("CHIMERYS" = msaid_blue, "CsoDIAq" = msaid_green,
                "Shared" = msaid_orange, "Pathway" = msaid_gray)


## Cytoscape ====
#export protein list
protein_export <- proteins_sub2[, .(protein = sort(unique(proteins))), keyby = nTypeLabel]
#fwrite(protein_export, file.path(figurePath, "cytoscape/cytoscape_proteins.csv"))


## import into Cytoscape and perform functional enrichment, then reimport ==


#map gene names onto protein list
gene_map <- fread(file.path(figurePath, "cytoscape/nodes-all.csv"),
                  select = c("display name", "query term"))
#fwrite(gene_map, file.path(figurePath, "cytoscape/nodes-all.csv"))
setnames(gene_map, c("display name", "query term"), c("gene", "protein"))
setkey(gene_map, protein)
setkey(protein_export, protein)
protein_export <- gene_map[protein_export]
#protein_export[is.na(gene)] #6 proteins were not mapped by STRING
protein_export[, protein := NULL]
protein_export <- protein_export[!is.na(gene)]
setkey(protein_export, gene)

#read back enrichment of all proteins
cytoscape_chimerys_path <- file.path(figurePath, "cytoscape/enrichment-chimerys.csv")
enrich_chimerys <- cbind(type="CHIMERYS", fread(cytoscape_chimerys_path))
cytoscape_csodiaq_path <- file.path(figurePath, "cytoscape/enrichment-csodiaq.csv")
enrich_csodiaq <- cbind(type="CsoDIAq", fread(cytoscape_csodiaq_path))
enrich <- rbind(enrich_chimerys, enrich_csodiaq)[category=="KEGG Pathways"]
names(enrich) <- make.names(names(enrich))
setkey(enrich, type, category, FDR.value)
enrich[, FDR.value.neg.log10 := -log10(FDR.value)]

#map gene counts
enrich <- enrich[, .(X..background.genes = X..background.genes[1],
                     genes = unique(unlist(strsplit(genes, "|", fixed = T))),
                     FDR.value.neg.log10=max(0, FDR.value.neg.log10[type=="CHIMERYS"])),
                 by=description]
setkey(enrich, genes)
enrich <- protein_export[enrich]
pie <- enrich[, .(variable = c("Shared", "CHIMERYS", "CsoDIAq", "Pathway"),
                  value = c(sum(nTypeLabel=="Shared"),
                            sum(nTypeLabel=="CHIMERYS"),
                            sum(nTypeLabel=="CsoDIAq"),
                            X..background.genes[1]-.N),
                  FDR.value.neg.log10 = round(FDR.value.neg.log10[1], 2)),
              by=description]
pie <- pie[value!=0]
pie[, variable := factor(variable, c("Pathway", "Shared", "CsoDIAq", "CHIMERYS"))]
setorder(pie, -FDR.value.neg.log10, variable)
pie[, description := paste0(description, "\n(-log10 q-value ", FDR.value.neg.log10, ")")]
pie[, description := factor(description, unique(pie$description))]
setkey(pie, description, variable)

term <- unique(pie$description)[c(1, 2, 3, 6)]
pie[, csum := rev(cumsum(rev(value))), by = description]
pie[, pos := value/2 + data.table::shift(csum, type = "lead"), by = description]
pie[, pos := ifelse(is.na(pos), value/2, pos), by = description]

fwrite(pie, file.path(figurePath, "figure-5H-cytoscape.csv"))
