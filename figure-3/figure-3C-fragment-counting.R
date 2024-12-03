# ---- setup ----
source(here::here("scripts/data-processing-2.R"))
path <- file.path(here::here(), "figure-3")
figurePath <- file.path(dataPath, "data/figure-3")
#fragmentPath <- file.path(dataPath, "figure-3/fragment-counting")

# ---- pathsToData ----
## chimerys
pathToPdResult <- file.path(figurePath, "fst-backup/20241127_pdresult_height_pepEntr_pcms_localPcmGrouper.fst")

## dia-nn
pathToDiann <- file.path(figurePath, 'fst-backup/20241127_diann_height_pepEntr_pcms.fst')

## spectronaut
pathToSn19 <- file.path(figurePath, 'fst-backup/20241127_sn19_height_noNorm_pepEntr_pcms.fst')
pathToSn19_flaggedImputation <- file.path(figurePath, 'fst-backup/20241127_sn19_height_noNorm_flaggedImputation_pepEntr_pcms.fst')


# ---- read data ----
pdresult <- read.fst(pathToPdResult, as.data.table = T)
diann <- read.fst(pathToDiann, as.data.table = T)
sn19 <- read.fst(pathToSn19, as.data.table = T)
sn19_flaggedImputation <- read.fst(pathToSn19_flaggedImputation, as.data.table = T)

# ---- load diann number of fragments ----
diann_quanFrags <- read.fst(file.path(figurePath, 'fst-backup/20241127_diann_lfq_height_pepEntr_fragCounts.fst'), as.data.table = T)

# ---- diann; add number of fragments for quantification ----
diann <- merge(diann,
               diann_quanFrags[,.(Precursor.Id,
                                  File.Name,
                                  DIANN_PRED_FRAGMENTS,
                                  DIANN_ASSAY_FRAGMENTS,
                                  F.CORRECTED_QUAN_VALUES,
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
pdresult[,ID_FRAGS:=FRAGMENTS_USED]
pdresult[,QUAN_FRAGS:=FRAGMENTS_USED]

colnames(diann)[!colnames(diann) %in% colnames(pdresult)]
colnames(diann)[!colnames(diann) %in% colnames(sn19)]
colnames(diann)[!colnames(diann) %in% colnames(sn19_filtered)]
diann[,DIANN_ORGANISM:=NULL]
diann[,ID_FRAGS:=F.CORRECTED_QUAN_VALUES]
diann[,QUAN_FRAGS:=NotExcludedFromQuanMax]

colnames(sn19_flaggedImputation)[!colnames(sn19) %in% colnames(pdresult)]
colnames(sn19_flaggedImputation)[!colnames(sn19) %in% colnames(diann)]
colnames(sn19_flaggedImputation)[!colnames(sn19) %in% colnames(sn19_filtered)]
sn19_flaggedImputation[,SN_ORGANISM:=NULL]
sn19_flaggedImputation[,ID_FRAGS:=MIN1_ID_FRAGS]
sn19_flaggedImputation[,QUAN_FRAGS:=MIN1_QUAN_FRAGS]

colnames(sn19_filtered)[!colnames(sn19_filtered) %in% colnames(pdresult)]
colnames(sn19_filtered)[!colnames(sn19_filtered) %in% colnames(diann)]
colnames(sn19_filtered)[!colnames(sn19_filtered) %in% colnames(sn19)]
sn19_filtered[,SN_ORGANISM:=NULL]
sn19_filtered[,ID_FRAGS:=MIN1_ID_FRAGS]
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
combined[DECOY==1,ORGANISM:='DECOY']

# ---- sanity checks ----
table(combined$ORGANISM)
table(combined$ORGANISM_J)
table(combined$SOFTWARE)

# ---- calculate fdr and entrapment fdr ----
combined[,ORGANISM:=ORGANISM_J]
combined[,ENTRAPMENT:=ENTRAPMENT_J]
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
write.fst(combined, file.path(figurePath, 'fst-backup/20241127_figure3c_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst'), compress = 100)

# ---- re-read intermediately saved data ----
# combined <- read.fst(file.path(figurePath, 'fst-backup/20241127_figure3c_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst'), as.data.table = T)

combined[SOFTWARE != 'CHIMERYS', ENTRAPMENT_Q_VALUE:=ENTRAPMENT_Q_VALUE_1]

# ---- fdr baed id frag count ----
idFrags_combined_fdr <- combined[Q_VALUE <= 0.01 &
                                   DECOY == 0 &
                                   SOFTWARE != 'SPECTRONAUT_FILTERED',
                                 list(PCM_ID, SAMPLE, CONDITION,
                                      ENTRAPMENT_Q_VALUE,
                                      QUAN,
                                      QUAN_CONDITION_FDR,
                                      CONDITION_MIN2_QUAN_FDR,
                                      ID_FRAGS, ORGANISM, SOFTWARE)]
idFrags_combined_fdr <- unique(idFrags_combined_fdr)
table(idFrags_combined_fdr$ID_FRAGS)

# ---- intermediate save for AH ----
write.fst(idFrags_combined_fdr, file.path(figurePath, '20241127_figure3c_idFrags_noNorm_pepFasta_localPcmFdr_pepFasta_apexQuan.fst'), compress = 100)

# ---- re-read intermediately saved data ----
# idFrags_combined_fdr <- read.fst(file.path(figurePath, '20241127_figure3c_idFrags_noNorm_pepFasta_localPcmFdr_pepFasta_apexQuan.fst'), as.data.table = T)

min(idFrags_combined_fdr$ID_FRAGS)
max(idFrags_combined_fdr$ID_FRAGS)

# ---- preliminary plot ----
ggplot(idFrags_combined_fdr,
       mapping = aes(ID_FRAGS,
                     fill = factor(ENTRAPMENT_Q_VALUE <= 0.01))) +
  geom_bar(stat = 'count', position = 'dodge') +
  facet_wrap(.~SOFTWARE,
             strip.position = 'bottom',
             scales = 'free') +
  theme(text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12 ),
        axis.ticks = element_line(size = .5),
        axis.title = element_text(size = 12),
        legend.position = 'top',
        legend.text=element_text(size = 12),
        legend.margin = margin(c(0,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside") +
  guides(color = guide_legend(
    override.aes=list(shape = 19, size = 2))) +
  labs(x = '# fragments for precursor identification',
       y = '# of file-local precursors',
       fill = 'Entrapment FDR')

# ---- fdr based quan frag count ----
quanFrags_combined_fdr <- combined[Q_VALUE <= 0.01 &
                                     DECOY == 0 &
                                     SOFTWARE != 'SPECTRONAUT_FILTERED',
                                   list(PCM_ID, SAMPLE, CONDITION,
                                        ENTRAPMENT_Q_VALUE,
                                        QUAN,
                                        QUAN_CONDITION_FDR,
                                        CONDITION_MIN2_QUAN_FDR,
                                        QUAN_FRAGS, ORGANISM, SOFTWARE)]
quanFrags_combined_fdr <- unique(quanFrags_combined_fdr)
table(quanFrags_combined_fdr$QUAN_FRAGS)

# ---- intermediate save for AH ----
write.fst(quanFrags_combined_fdr, file.path(figurePath, '20241127_figure3c_quanFrags_noNorm_pepFasta_localPcmFdr_pepFasta_apexQuan.fst'), compress = 100)

# ---- re-read intermediately saved data ----
# quanFrags_combined_fdr <- read.fst(file.path(figurePath, '20241127_figure3c_quanFrags_noNorm_pepFasta_localPcmFdr_pepFasta_apexQuan.fst'), as.data.table = T)

# ---- preliminary plot ----
ggplot(quanFrags_combined_fdr,
       mapping = aes(QUAN_FRAGS,
                     fill = factor(ENTRAPMENT_Q_VALUE <= 0.01))) +
  geom_bar(stat = 'count', position = 'dodge') +
  facet_wrap(.~SOFTWARE,
             strip.position = 'bottom',
             scales = 'free') +
  theme(text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12 ),
        axis.ticks = element_line(size = .5),
        axis.title = element_text(size = 12),
        legend.position = 'top',
        legend.text=element_text(size = 12),
        legend.margin = margin(c(0,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside") +
  guides(color = guide_legend(
    override.aes=list(shape = 19, size = 2))) +
  labs(x = '# fragments for precursor identification',
       y = '# of file-local precursors',
       fill = 'Entrapment FDR')
