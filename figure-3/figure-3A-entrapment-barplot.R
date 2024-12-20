# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-3")
figurePath <- file.path(dataPath, "data/figure-3")

# ---- pathsToData ----
## chimerys
pathToPdResult <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/20240503_lfq_dia_entrapment_peptides_z1to4_localPcmFdr_True.pdResult")
pathToMs8 <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Chimerys-entrapment/main_search-3.ms8")

## dia-nn
pathToTsv <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN/20240429_lfq_height_entrapment_peptides_report.tsv")

## spectronaut
pathToExport <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Spectronaut/20240925_064248_20240925_SN19_lfq_paper_entrapment_paper_Report_height_noNorm.tsv")

## fastas
pathToFasta_peptides <- file.path(dataPath, "FASTA/CHIMERYS_Benchmark_human-canonical_yeast_ecoli_AND_jpr_2022_contaminants_mimic_peptides.fasta")


# ---- read data including short sanity checks ----
proteinsPeptideFasta <- readProteinsFromFastas(pathToFasta_peptides)

pdresult <- readPdResult_localPcmGrouper(pathToPdResult, pathToMs8, proteinsPeptideFasta)
pdresult[,ORGANISM_J:=ORGANISM]
pdresult[,ENTRAPMENT_J:=ENTRAPMENT]
pdresult[,ORGANISM_FASTA_J:=ORGANISM]
pdresult[,ENTRAPMENT_FASTA_J:=ENTRAPMENT]
setcolorder(pdresult, neworder = c('ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION', 'Q_VALUE', 'SVMSCORE', 'DECOY', 'CHARGE', 'MODIFIED_SEQUENCE', 'QUAN', 'QUANTIFICATION_VALUE', 'IS_IDENTIFIED_BY_MBR', 'SAMPLE', 'FRAGMENTS_USED', 'SPECTRUM_SIMILARITY', 'DPPP', 'APEX_MS8_feature_7', 'ORGANISM', 'ORGANISM_J', 'ORGANISM_FASTA_J', 'ENTRAPMENT', 'ENTRAPMENT_J', 'ENTRAPMENT_FASTA_J', 'PCM_ID', 'PCM_J_ID', 'MODPEP_ID', 'MODPEP_J_ID', 'SOFTWARE'))
any(is.na(pdresult$ORGANISM_FASTA_J)) # FALSE
table(pdresult$ENTRAPMENT_FASTA_J)
table(pdresult$ORGANISM_FASTA_J) # ENTRAPMENT EXISTS
uniqueN(pdresult[Q_VALUE <= 0.01, PCM_J_ID]) # 64048
uniqueN(pdresult[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 62481

diann <- readDiann(pathToTsv, proteinsPeptideFasta)
any(is.na(diann$ORGANISM_FASTA_J)) # FALSE
table(diann$ENTRAPMENT)
table(diann$ENTRAPMENT_FASTA_J)
table(diann$ORGANISM_FASTA_J) # ENTRAPMENT EXISTS
uniqueN(diann[Q_VALUE <= 0.01, PCM_ID]) # 79259
uniqueN(diann[Q_VALUE <= 0.01, PCM_J_ID]) # 79181

sn19 <- readSpectronaut(pathToExport, proteinsPeptideFasta)
any(is.na(sn19$ORGANISM_FASTA_J)) # FALSE
table(sn19$ENTRAPMENT)
table(sn19$ENTRAPMENT_FASTA_J)
table(sn19$ORGANISM_FASTA_J) # ENTRAPMENT EXISTS
uniqueN(sn19[Q_VALUE <= 0.01, PCM_ID]) # 67740
uniqueN(sn19[Q_VALUE <= 0.01, PCM_J_ID]) # 67704
max(sn19[!is.nan(Q_VALUE), Q_VALUE]) # 0.1369829

sn19_flaggedImputation <- readSpectronaut_flagImputation(pathToExport, proteinsPeptideFasta)
any(is.na(sn19_flaggedImputation$ORGANISM_FASTA_J)) # FALSE
table(sn19_flaggedImputation$ENTRAPMENT)
table(sn19_flaggedImputation$ENTRAPMENT_FASTA_J)
table(sn19_flaggedImputation$ORGANISM_FASTA_J) # ENTRAPMENT EXISTS
uniqueN(sn19_flaggedImputation[Q_VALUE <= 0.01, PCM_ID]) # 67740
uniqueN(sn19_flaggedImputation[Q_VALUE <= 0.01, PCM_J_ID]) # 67704

# # ---- intermediate save ----
write.fst(pdresult, file.path(figurePath, 'fst-backup/20241127_pdresult_height_pepEntr_pcms_localPcmGrouper.fst'), compress = 100)
write.fst(diann, file.path(figurePath, 'fst-backup/20241127_diann_height_pepEntr_pcms.fst'), compress = 100)
write.fst(sn19, file.path(figurePath, 'fst-backup/20241127_sn19_height_noNorm_pepEntr_pcms.fst'), compress = 100)
write.fst(sn19_flaggedImputation, file.path(figurePath, 'fst-backup/20241127_sn19_height_noNorm_flaggedImputation_pepEntr_pcms.fst'), compress = 100)

# ---- load diann number of fragments ----
diann_quanFrags <- read.fst(file.path(figurePath, 'fst-backup/20241127_diann_lfq_height_pepEntr_fragCounts.fst'), as.data.table = T)

# ---- diann; add number of fragments for quantification ----
diann <- merge(diann,
               diann_quanFrags[,.(Precursor.Id,
                                  File.Name,
                                  DIANN_PRED_FRAGMENTS,
                                  DIANN_ASSAY_FRAGMENTS,
                                  NotExcludedFromQuanMin,
                                  NotExcludedFromQuanMax)],
               by.x = c('ID', 'SAMPLE'),
               by.y = c('Precursor.Id', 'File.Name'))

# ---- remove ids not fulfilling fragment requirements ----
sn19_filtered <- sn19_flaggedImputation[MIN1_QUAN_FRAGS >= 3]

# ---- change software for filtered spectronaut ----
sn19_filtered[,SOFTWARE:='SPECTRONAUT_FILTERED']

# ---- combine results from different softwares ----
colnames(pdresult)[!colnames(pdresult) %in% colnames(diann)]
colnames(pdresult)[!colnames(pdresult) %in% colnames(sn19)]
colnames(pdresult)[!colnames(pdresult) %in% colnames(sn19_filtered)]
pdresult[,QUAN_FRAGS:=FRAGMENTS_USED]

colnames(diann)[!colnames(diann) %in% colnames(pdresult)]
colnames(diann)[!colnames(diann) %in% colnames(sn19)]
colnames(diann)[!colnames(diann) %in% colnames(sn19_filtered)]
diann[,DIANN_ORGANISM:=NULL]
diann[,QUAN_FRAGS:=NotExcludedFromQuanMax]

colnames(sn19_flaggedImputation)[!colnames(sn19) %in% colnames(pdresult)]
colnames(sn19_flaggedImputation)[!colnames(sn19) %in% colnames(diann)]
colnames(sn19_flaggedImputation)[!colnames(sn19) %in% colnames(sn19_filtered)]
sn19_flaggedImputation[,SN_ORGANISM:=NULL]
sn19_flaggedImputation[,QUAN_FRAGS:=MIN1_QUAN_FRAGS]

colnames(sn19_filtered)[!colnames(sn19_filtered) %in% colnames(pdresult)]
colnames(sn19_filtered)[!colnames(sn19_filtered) %in% colnames(diann)]
colnames(sn19_filtered)[!colnames(sn19_filtered) %in% colnames(sn19)]
sn19_filtered[,SN_ORGANISM:=NULL]
sn19_filtered[,QUAN_FRAGS:=MIN1_QUAN_FRAGS]

combined <- rbindlist(list(pdresult, diann, sn19_flaggedImputation, sn19_filtered), fill = T)
table(combined$ORGANISM)
table(combined$ENTRAPMENT)
combined[,ORGANISM:=gsub('(^|;)+ENTRAPMENT', '', ORGANISM)]
combined[,ORGANISM:=gsub('^;', '', ORGANISM)]
combined[ORGANISM=='',ORGANISM:='ENTRAPMENT']
combined[DECOY==1,ORGANISM:='DECOY']

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

# ---- sanity checks ----
table(combined$ORGANISM)
table(combined$ORGANISM_FASTA_J)
table(combined$SOFTWARE)

# ---- calculate fdr and entrapment fdr ----
combined[,ORGANISM:=ORGANISM_FASTA_J]
combined[,ENTRAPMENT:=ENTRAPMENT_FASTA_J]
combined <- setEntrapmentFDR_local(combined[!is.na(SVMSCORE)])

# ---- remove decoys ----
combined <- combined[DECOY == 0]
uniqueN(combined[SOFTWARE == 'CHIMERYS' & Q_VALUE <= 0.01 & DECOY == 0,PCM_J_ID]) # 62481

# ---- add condition ----
table(combined$SAMPLE)
combined[,CONDITION:=ifelse(grepl("Condition_A", SAMPLE),
                            "A",
                            "B")]

# ---- quantified samples per condition at fdr and entrapment fdr ----
combined[,QUAN_CONDITION_FDR:=uniqueN(SAMPLE[Q_VALUE<=0.01 & QUAN != 0 & !is.na(QUAN)]), by = .(SOFTWARE, CONDITION, PCM_ID)]
combined[SOFTWARE == 'CHIMERYS',QUAN_CONDITION_ENTRAPMENT_FDR:=uniqueN(SAMPLE[Q_VALUE<=0.01 & ENTRAPMENT_Q_VALUE<=0.01 & QUAN != 0 & !is.na(QUAN)]), by = .(SOFTWARE, CONDITION, PCM_ID)]
combined[SOFTWARE != 'CHIMERYS',QUAN_CONDITION_ENTRAPMENT_FDR:=uniqueN(SAMPLE[Q_VALUE<=0.01 & ENTRAPMENT_Q_VALUE_1<=0.01 & QUAN != 0 & !is.na(QUAN)]), by = .(SOFTWARE, CONDITION, PCM_ID)]

# ---- quantified samples per dataset at fdr and entrapment fdr ----
combined[,QUAN_SAMPLES_FDR:=uniqueN(SAMPLE[Q_VALUE<=0.01 & QUAN != 0 & !is.na(QUAN)]), by = .(SOFTWARE, PCM_ID)]
combined[SOFTWARE == 'CHIMERYS',QUAN_SAMPLES_ENTRAPMENT_FDR:=uniqueN(SAMPLE[Q_VALUE<=0.01 & ENTRAPMENT_Q_VALUE<=0.01 & QUAN != 0 & !is.na(QUAN)]), by = .(SOFTWARE,PCM_ID)]
combined[SOFTWARE != 'CHIMERYS',QUAN_SAMPLES_ENTRAPMENT_FDR:=uniqueN(SAMPLE[Q_VALUE<=0.01 & ENTRAPMENT_Q_VALUE_1<=0.01 & QUAN != 0 & !is.na(QUAN)]), by = .(SOFTWARE,PCM_ID)]

# ---- select precursors for cvs ----
combined[,CONDITION_MIN2_QUAN_FDR:=uniqueN(CONDITION[QUAN_CONDITION_FDR >= 2]), by = .(PCM_ID, SOFTWARE)]
combined[,CONDITION_MIN2_QUAN_EFDR:=uniqueN(CONDITION[QUAN_CONDITION_ENTRAPMENT_FDR >= 2]), by = .(PCM_ID, SOFTWARE)]
table(combined$CONDITION_MIN2_QUAN_EFDR)

# ---- flag quan complete ----
combined[,QUAN_COMPLETE_FDR:=QUAN_SAMPLES_FDR==6]
combined[,QUAN_COMPLETE_EFDR:=QUAN_SAMPLES_ENTRAPMENT_FDR==6]

# ---- common pcms ----
common_pcms <- Reduce(intersect, lapply(combined[,unique(SOFTWARE)], function(i) combined[SOFTWARE==i & Q_VALUE<=0.01,unique(PCM_J_ID)]))
length(common_pcms) # 49837

combined[,SHARED_PCM_J_ID := ifelse(PCM_J_ID %in% common_pcms, 1, 0)]

# ---- intermediate save ----
write.fst(combined, file.path(figurePath, 'fst-backup/20241127_figure3a_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst'), compress = 100)
colNames <- c("SOFTWARE", "QUAN", "Q_VALUE", "ENTRAPMENT_Q_VALUE", "ENTRAPMENT_Q_VALUE_1")
fwrite(combined[, .SD, .SDcols = colNames], file.path(figurePath, "figure-3B-quan.csv"))

# precursor readout for E7
dtPrecFrag <- combined[Q_VALUE <= 0.01]
dtPrecFrag[, LOG10QUAN := log10(QUAN)]
softwareLevels <- c("CHIMERYS", "DIA-NN", "SPECTRONAUT", "SPECTRONAUT_FILTERED")
softwareLabels <- c("CHIMERYS", "DIA-NN", "Spectronaut", "Spectronaut\n(curated)")
dtPrecFrag[, SOFTWARE := factor(SOFTWARE, softwareLevels, softwareLabels)]
dtPrecFrag[, MIN1_QUAN_FRAGS := factor(MIN1_QUAN_FRAGS, 0:6)]
dtPrecFrag <- dtPrecFrag[!is.na(LOG10QUAN) & SOFTWARE %in% softwareLabels[3]]
fwrite(dtPrecFrag, file.path(dataPath, "data/figure-E7/figure-E7D-precursors.csv"))

# # ---- re-read intermediately saved data ----
# combined <- read.fst(file.path(figurePath, 'fst-backup/20241127_figure3a_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst'), as.data.table = T)

# ---- count pcms stratified by min 2 quantified samples per condition at fdr and entrapment fdr ----
bardata_fdr <- melt(combined[Q_VALUE<=0.01,
                             list(FDR=uniqueN(PCM_J_ID)),
                             by=.(SOFTWARE, CONDITION_MIN2_QUAN_FDR)],
                    id.vars = c('SOFTWARE', 'CONDITION_MIN2_QUAN_FDR'),
                    value.name = '# of precursors in dataset',
                    variable.name = 'Survive 1%')

bardata_efdr_chimerys <- melt(combined[Q_VALUE<=0.01 &
                                         ENTRAPMENT_Q_VALUE<=0.01 &
                                         SOFTWARE == 'CHIMERYS',
                                       list(eFDR=uniqueN(PCM_J_ID)),
                                       by=.(SOFTWARE, CONDITION_MIN2_QUAN_EFDR)],
                              id.vars = c('SOFTWARE', 'CONDITION_MIN2_QUAN_EFDR'),
                              value.name = '# of precursors in dataset',
                              variable.name = 'Survive 1%')

bardata_efdr <- melt(combined[Q_VALUE<=0.01 &
                                ENTRAPMENT_Q_VALUE_1<=0.01 &
                                SOFTWARE != 'CHIMERYS',
                              list(eFDR=uniqueN(PCM_J_ID)),
                              by=.(SOFTWARE, CONDITION_MIN2_QUAN_EFDR)],
                     id.vars = c('SOFTWARE', 'CONDITION_MIN2_QUAN_EFDR'),
                     value.name = '# of precursors in dataset',
                     variable.name = 'Survive 1%')

setnames(bardata_efdr_chimerys, 'CONDITION_MIN2_QUAN_EFDR', 'CONDITION_MIN2_QUAN_FDR')
setnames(bardata_efdr, 'CONDITION_MIN2_QUAN_EFDR', 'CONDITION_MIN2_QUAN_FDR')
bardata <- rbindlist(list(bardata_fdr, bardata_efdr_chimerys, bardata_efdr))

bardata <- bardata[!SOFTWARE=='SPECTRONAUT_FILTERED']

# ---- intermediate save ----
fwrite(bardata, file.path(figurePath, "figure-3A-entrapment.csv"))

# ---- preliminary plot ----
ggplot(bardata,
       mapping = aes(x = `Survive 1%`,
                     y = `# of precursors in dataset`,
                     fill = factor(CONDITION_MIN2_QUAN_FDR == 2,
                                   levels = c(FALSE, TRUE)))) +
  geom_bar(position = 'stack', stat = 'identity') +
  facet_wrap(.~factor(SOFTWARE,
                      levels = c('CHIMERYS',
                                 'DIA-NN',
                                 'SPECTRONAUT')),
             strip.position = 'bottom') +
  theme(text = element_text(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = .5),
        legend.margin = margin(c(0,0,0,0)),
        legend.position = 'top',
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside",
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5)) +
  scale_y_continuous(labels = label_number(suffix = "k", scale = 1e-3, sep = '')) +
  labs(y = '# of pcms in dataset',
       fill = "Quantified in ≥ 2 replicates per condition")
