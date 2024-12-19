# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-S8-entrapment")
efdrLevels <- c("REGULAR", "PEPTIDE", "CONCATENATED")
efdrLabels <- c("Classic eFDR", "Peptide eFDR", "Concatenated eFDR")

# ---- pathsToData ----
## diann
pathToTsv <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN_timsTof/20240924_lfq_timstof_entrapment_height_report.tsv")
pathToTsv_pepEntr <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN_timsTof/20240924_lfq_timstof_entrapment_peptides_height_report.tsv")
pathToTsv_concatEntr <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN_timsTof/20240924_lfq_timstof_entrapment_concat_height_report.tsv")

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

diann_entr <- readDiann(pathToTsv, proteinsFasta)
any(is.na(diann_entr$ORGANISM_FASTA_J)) # FALSE
table(diann_entr$ENTRAPMENT)
table(diann_entr$ENTRAPMENT_FASTA_J)
table(diann_entr$ORGANISM_FASTA_J) # ENTRAPMENT EXISTS
uniqueN(diann_entr[Q_VALUE <= 0.01, PCM_J_ID]) # 60282

diann_pepEntr <- readDiann(pathToTsv_pepEntr, proteinsPeptideFasta)
any(is.na(diann_pepEntr$ORGANISM_FASTA_J)) # FALSE
table(diann_pepEntr$ENTRAPMENT)
table(diann_pepEntr$ENTRAPMENT_FASTA_J)
table(diann_pepEntr$ORGANISM_FASTA_J) # ENTRAPMENT EXISTS
uniqueN(diann_pepEntr[Q_VALUE <= 0.01, PCM_J_ID]) # 61364

diann_concatEntr <- readDiann(pathToTsv_concatEntr,
                              proteinsConcatFasta,
                              digestConcatFasta)
any(is.na(diann_concatEntr$ORGANISM)) # FALSE
any(is.na(diann_concatEntr$ORGANISM_FASTA_J)) # FALSE
table(diann_concatEntr$ENTRAPMENT)
table(diann_concatEntr$ENTRAPMENT_FASTA_J)
table(diann_concatEntr$ORGANISM_FASTA_J) # ENTRAPMENT EXISTS
uniqueN(diann_concatEntr[Q_VALUE <= 0.01, PCM_J_ID]) # 60671

# ---- check entrapment decoys ----
any(is.na(diann_concatEntr$DECOY))
any(is.na(diann_concatEntr$ENTRAPMENT))

# ---- combine results from different entrapment analyses ----
diann_entr[,SOFTWARE:='REGULAR']
diann_pepEntr[,SOFTWARE:='PEPTIDE']
diann_concatEntr[,SOFTWARE:='CONCATENATED']
combined <- rbindlist(list(diann_entr,
                           diann_pepEntr,
                           diann_concatEntr),
                      fill = T)

table(combined$ORGANISM)
combined[,ORGANISM:=gsub('(^|;)+ENTRAPMENT', '', ORGANISM)]
combined[,ORGANISM:=gsub('^;', '', ORGANISM)]
combined[ORGANISM=='',ORGANISM:='ENTRAPMENT']
combined[DECOY==1,ORGANISM:='DECOY']
table(combined$ORGANISM)
table(combined$SOFTWARE)

table(combined$ORGANISM_J)
table(combined$ENTRAPMENT_J)
combined[,ORGANISM_J:=gsub('(^|;)+ENTRAPMENT', '', ORGANISM_J)]
combined[,ORGANISM_J:=gsub('^;', '', ORGANISM_J)]
combined[ORGANISM_J=='',ORGANISM_J:='ENTRAPMENT']
combined[DECOY==1,ORGANISM_J:='DECOY']

table(combined$ORGANISM_FASTA_J)
table(combined$ENTRAPMENT_FASTA_J)
combined[,ORGANISM_FASTA_J:=gsub('(^|;)+ENTRAPMENT', '', ORGANISM_FASTA_J)]
combined[,ORGANISM_FASTA_J:=gsub('^;', '', ORGANISM_FASTA_J)]
combined[ORGANISM_FASTA_J=='',ORGANISM_FASTA_J:='ENTRAPMENT']
combined[DECOY==1,ORGANISM_FASTA_J:='DECOY']

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
write.fst(combined, file.path(dataPath, 'data/figure-S8/fst-backup/20241111_figureS8e_diann_bruker_combined_pcms_localPcmEfdr_apexQuan_entr.fst'), compress = 100)

data_efdr <- combined[Q_VALUE<=0.055]
data_efdr[, SOFTWARE := factor(SOFTWARE, efdrLevels, efdrLabels)]
data_efdr[, Q_VALUE_BIN := ceiling(Q_VALUE*10000)/10000]
mean_efdrDiannBruk <- data_efdr[, .(ENTRAPMENT_Q_VALUE = mean(ENTRAPMENT_Q_VALUE),
                                    ENTRAPMENT_Q_VALUE_1 = mean(ENTRAPMENT_Q_VALUE_1)),
                                by=.(SOFTWARE, CONDITION, SAMPLE, Q_VALUE_BIN)]
setnames(mean_efdrDiannBruk, "Q_VALUE_BIN", "Q_VALUE")
mean_efdrDiannBruk[, SAMPLE := gsub("^LFQ_timsTOFPro_diaPASEF_(Condition_._Sample_Alpha_0.)$", "\\1", SAMPLE)]

mean_efdrDiannBruk <-
  rbind(mean_efdrDiannBruk[, .(Q_VALUE = 0, ENTRAPMENT_Q_VALUE = 0, ENTRAPMENT_Q_VALUE_1 = 0),
                           by=.(SOFTWARE, CONDITION, SAMPLE)], mean_efdrDiannBruk)
fwrite(mean_efdrDiannBruk, file.path(figurePath, "figure-E6E-diann-bruker.csv"))

# # ---- re-read intermediately saved data ----
# combined <- read.fst(file.path(dataPath, 'data/figure-S8/fst-backup/20241111_figureS8e_diann_bruker_combined_pcms_localPcmEfdr_apexQuan_entr.fst'), as.data.table = T)
#
# # ---- plot fdr vs efdr ----
# ggplot(data = combined,
#        mapping = aes(x = Q_VALUE, y = ENTRAPMENT_Q_VALUE_1,
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
#         legend.text = element_text(size = 12),
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
