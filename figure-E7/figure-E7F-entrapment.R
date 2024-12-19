# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-E7")
figurePath <- file.path(dataPath, "data/figure-E7")

# ---- pathsToData ----
## combined file generated in figure-3/figure-3A-entrapment-barplot.R
pathToCombined <- file.path(dataPath, 'data/figure-3/fst-backup/20241127_figure3a_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst')

# ---- read data ----
combined <- read.fst(pathToCombined, as.data.table = T)

# ---- select precursors for cvs ----
combined[,CONDITION_MIN1_QUAN_FDR:=uniqueN(CONDITION[QUAN_CONDITION_FDR >= 1]), by = .(PCM_ID, SOFTWARE)]
combined[,CONDITION_MIN1_QUAN_EFDR:=uniqueN(CONDITION[QUAN_CONDITION_ENTRAPMENT_FDR >= 1]), by = .(PCM_ID, SOFTWARE)]

# ---- remove spectronaut filtered ----
combined <- combined[SOFTWARE != 'SPECTRONAUT_FILTERED']
sn18_flagged <- combined[SOFTWARE == 'SPECTRONAUT']
table(sn18_flagged$MIN1_QUAN_FRAGS)
sn18_filtered <- sn18_flagged[MIN1_QUAN_FRAGS >= 3]
table(sn18_filtered$MIN1_QUAN_FRAGS)
sn18_filtered[,SOFTWARE:='SPECTRONAUT_FILTERED']
combined <- rbind(combined, sn18_filtered)

# ---- count pcms stratified by full data completeness at fdr and entrapment fdr ----
bardata_fdr <- melt(combined[Q_VALUE<=0.01,
                             list(FDR=uniqueN(PCM_J_ID)),
                             by=.(SOFTWARE, QUAN_COMPLETE_FDR)],
                    id.vars = c('SOFTWARE', 'QUAN_COMPLETE_FDR'),
                    value.name = '# of precursors in dataset',
                    variable.name = 'Survive 1%')

bardata_efdr_chimerys <- melt(combined[Q_VALUE<=0.01 &
                                ENTRAPMENT_Q_VALUE<=0.01 &
                                  SOFTWARE == 'CHIMERYS',
                              list(eFDR=uniqueN(PCM_J_ID)),
                              by=.(SOFTWARE, QUAN_COMPLETE_EFDR)],
                     id.vars = c('SOFTWARE', 'QUAN_COMPLETE_EFDR'),
                     value.name = '# of precursors in dataset',
                     variable.name = 'Survive 1%')

bardata_efdr <- melt(combined[Q_VALUE<=0.01 &
                                ENTRAPMENT_Q_VALUE_1<=0.01 &
                                SOFTWARE != 'CHIMERYS',
                              list(eFDR=uniqueN(PCM_J_ID)),
                              by=.(SOFTWARE, QUAN_COMPLETE_EFDR)],
                     id.vars = c('SOFTWARE', 'QUAN_COMPLETE_EFDR'),
                     value.name = '# of precursors in dataset',
                     variable.name = 'Survive 1%')

setnames(bardata_efdr_chimerys, 'QUAN_COMPLETE_EFDR', 'QUAN_COMPLETE_FDR')
setnames(bardata_efdr, 'QUAN_COMPLETE_EFDR', 'QUAN_COMPLETE_FDR')
bardata <- rbindlist(list(bardata_fdr, bardata_efdr_chimerys, bardata_efdr))

bardata <- bardata[!SOFTWARE %in% c('CHIMERYS', 'DIA-NN')]
bardata <- bardata[!(SOFTWARE == 'SPECTRONAUT_FILTERED' & `Survive 1%` == 'FDR')]
bardata[,LABEL:=paste0(`Survive 1%`, ' - ', SOFTWARE)]
bardata[,LABEL:=factor(LABEL,
                       levels = c('FDR - SPECTRONAUT',
                                  'eFDR - SPECTRONAUT',
                                  'eFDR - SPECTRONAUT_FILTERED'))]

setnames(bardata, c("Survive 1%", "# of precursors in dataset"), c("FDR", "N"))
bardata[, QUAN_COMPLETE_FDR := factor(QUAN_COMPLETE_FDR, c(T, F))]
efdrLevels <- c("FDR - SPECTRONAUT", "eFDR - SPECTRONAUT", "eFDR - SPECTRONAUT_FILTERED")
efdrLabels <- c("FDR", "FDR +\neFDR", "FDR + eFDR +\nquan fragments\nfiltered")
bardata[, LABEL := factor(LABEL, efdrLevels, efdrLabels)]

# ---- intermediate save ----
fwrite(bardata, file.path(figurePath, "figure-E7F-efdr.csv"))
