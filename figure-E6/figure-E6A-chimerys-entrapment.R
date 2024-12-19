# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-E6")
figurePath <- file.path(dataPath, "data/figure-E6")

efdrLevels <- c("REGULAR", "PEPTIDE", "CONCATENATED")
efdrLabels <- c("Classic eFDR", "Peptide eFDR", "Concatenated eFDR")

# ---- pathsToData ----
## chimerys
pathToPdResult <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/20240503_lfq_dia_entrapment_z1to4_localPcmFdr_True.pdResult")
pathToMs8 <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/main_search-2.ms8")

pathToPdResult_pepEntr <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/20240503_lfq_dia_entrapment_peptides_z1to4_localPcmFdr_True.pdResult")
pathToMs8_pepEntr <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/main_search-3.ms8")

pathToPdResult_concatEntr <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/20240524_lfq_dia_entrapment_concat_z1to4_localPcmFdr_True.pdResult")
pathToMs8_concatEntr <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/main_search-54.ms8")

## fastas
pathToFasta <- file.path(dataPath, "FASTA/CHIMERYS_Benchmark_human-canonical_yeast_ecoli_AND_jpr_2022_contaminants_mimic.fasta")
pathToFasta_peptides <- file.path(dataPath, "FASTA/CHIMERYS_Benchmark_human-canonical_yeast_ecoli_AND_jpr_2022_contaminants_mimic_peptides.fasta")
pathToFasta_concat <- file.path(dataPath, "FASTA/CHIMERYS_Benchmark_human-canonical_yeast_ecoli_AND_jpr_2022_contaminants_mimic_concat_limit1e5.fasta")

## digest
pathToDigest <- file.path(dataPath, "FASTA/CHIMERYS_Benchmark_human-canonical_yeast_ecoli_AND_jpr_2022_contaminants_mimic_digested_Mass0to10000.txt")


# ---- read data including short sanity checks ----
proteinsFasta <- readProteinsFromFastas(pathToFasta)
proteinsPeptideFasta <- readProteinsFromFastas(pathToFasta_peptides)
proteinsConcatFasta <- readProteinsFromFastas(pathToFasta_concat)
digestConcatFasta <- readDigest(pathToDigest)

pdresult_entr <- readPdResult_localPcmGrouper(pathToPdResult, pathToMs8, proteinsFasta)
any(is.na(pdresult_entr$ORGANISM)) # FALSE
table(pdresult_entr$ENTRAPMENT)
table(pdresult_entr$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_entr[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 63048

pdresult_pepEntr <- readPdResult_localPcmGrouper(pathToPdResult_pepEntr, pathToMs8_pepEntr, proteinsPeptideFasta)
any(is.na(pdresult_pepEntr$ORGANISM)) # FALSE
table(pdresult_pepEntr$ENTRAPMENT)
table(pdresult_pepEntr$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_pepEntr[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 62481

pdresult_concatEntr <- readPdResult_localPcmGrouper(pathToPdResult_concatEntr,
                                                    pathToMs8_concatEntr,
                                                    proteinsConcatFasta,
                                                    digestConcatFasta)
any(is.na(pdresult_concatEntr$ORGANISM)) # FALSE
table(pdresult_concatEntr$ENTRAPMENT)
table(pdresult_concatEntr$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_concatEntr[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 62183

# ---- check entrapment decoys ----
any(is.na(pdresult_concatEntr$DECOY))
any(is.na(pdresult_concatEntr$ENTRAPMENT))

# ---- combine results from different entrapment analyses ----
pdresult_entr[,SOFTWARE:='REGULAR']
pdresult_pepEntr[,SOFTWARE:='PEPTIDE']
pdresult_concatEntr[,SOFTWARE:='CONCATENATED']
combined <- rbindlist(list(pdresult_entr,
                           pdresult_pepEntr,
                           pdresult_concatEntr),
                      fill = T)

table(combined$ORGANISM)
combined[,ORGANISM:=gsub('(^|;)+ENTRAPMENT', '', ORGANISM)]
combined[,ORGANISM:=gsub('^;', '', ORGANISM)]
combined[ORGANISM=='',ORGANISM:='ENTRAPMENT']
combined[DECOY==1,ORGANISM:='DECOY']
table(combined$ORGANISM)
table(combined$SOFTWARE)

combined[,ORGANISM_FASTA_J:=ORGANISM]
combined[,ENTRAPMENT_FASTA_J:=ENTRAPMENT]

# ---- calculate fdr and entrapment fdr ----
combined <- setEntrapmentFDR_local(combined[!is.na(SVMSCORE)])

# ---- remove decoys ----
combined <- combined[DECOY == 0]

# ---- add condition ----
table(combined$SAMPLE)
combined[,CONDITION:=ifelse(grepl("Condition_A", SAMPLE),
                            "A",
                            "B")]

# ---- intermediate save ----
write.fst(combined, file.path(dataPath, 'data/figure-S8/fst-backup/20241111_figureS8a_pdresult_combined_pcms_localPcmEfdr_apexQuan_entr.fst'), compress = 100)

data_efdr <- combined[Q_VALUE<=0.055]
data_efdr[, SOFTWARE := factor(SOFTWARE, efdrLevels, efdrLabels)]
data_efdr[, Q_VALUE_BIN := ceiling(Q_VALUE*10000)/10000]
mean_efdrChimerys <- data_efdr[, .(ENTRAPMENT_Q_VALUE = mean(ENTRAPMENT_Q_VALUE),
                                   ENTRAPMENT_Q_VALUE_1 = mean(ENTRAPMENT_Q_VALUE_1)),
                               by=.(SOFTWARE, CONDITION, SAMPLE, Q_VALUE_BIN)]
setnames(mean_efdrChimerys, "Q_VALUE_BIN", "Q_VALUE")
mean_efdrChimerys[, SAMPLE := gsub("^LFQ_Orbitrap_AIF_(Condition_._Sample_Alpha_0.)$", "\\1", SAMPLE)]

mean_efdrChimerys <-
  rbind(mean_efdrChimerys[, .(Q_VALUE = 0, ENTRAPMENT_Q_VALUE = 0, ENTRAPMENT_Q_VALUE_1 = 0),
                          by=.(SOFTWARE, CONDITION, SAMPLE)], mean_efdrChimerys)

mean_efdrChimerys[, ENTRAPMENT_Q_VALUE_MIXED := ifelse(SOFTWARE=="Concatenated eFDR",
                                                       ENTRAPMENT_Q_VALUE_1, ENTRAPMENT_Q_VALUE)]
fwrite(mean_efdrChimerys, file.path(figurePath, "figure-E6A-chimerys.csv"))

# # ---- re-read intermediately saved data ----
# combined <- read.fst(file.path(dataPath, 'data/figure-S8/fst-backup/20241111_figureS8a_pdresult_combined_pcms_localPcmEfdr_apexQuan_entr.fst'), as.data.table = T)
# combined[SOFTWARE == 'CONCATENATED', ENTRAPMENT_Q_VALUE:=ENTRAPMENT_Q_VALUE_1]
#
# # ---- plot fdr vs efdr; preliminary plot for testing ----
# ggplot(data = combined,
#        mapping = aes(x = Q_VALUE, y = ENTRAPMENT_Q_VALUE,
#                      color = CONDITION,
#                      group = paste0(SOFTWARE, '_', SAMPLE))) +
#   geom_abline(intercept = 0, slope = 1, col = 1, lty = 2) +
#   geom_abline(intercept = 0, slope = 1.5, col = 'lightgray', lty = 2) +
#   geom_abline(intercept = 0, slope = 2/3, col = 'lightgray', lty = 2) +
#   geom_vline(xintercept = 0.01, col = 2, lty = 2) +
#   geom_line() +
#   facet_wrap(~SOFTWARE) +
#   theme(text = element_text(size = 12),
#         panel.grid = element_blank(),
#         # legend.position = 'top',
#         legend.position = 'top',
#         legend.text = element_text( size = 12),
#         axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         axis.ticks = element_line(linewidth = .5),
#         legend.margin = margin(c(0,0,0,0)),
#         panel.border = element_blank(),
#         axis.line = element_line(size = 0.5),
#         strip.background = element_blank()) +
#   labs(y = 'File-local precursor entrapment FDR',
#        x = 'File-local precursor FDR') +
#   xlim(c(0,0.05)) +
#   ylim(c(0,0.05))
