# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-4")
figurePath <- file.path(dataPath, "data/figure-4")

# ---- pathsToData ----
## combined
pathToCombined <- file.path(figurePath, 'intermediate/20241127_figure4c_combined_pcms_localPcmGrouper_apexQuan_pepEntr.fst')

# ---- read data ----
combined <- read.fst(pathToCombined, as.data.table = T)

combined[CONDITION_MIN2_QUAN_FDR == 2 &
           QUAN_CONDITION_FDR >= 3 & #CHANGED from 2
           Q_VALUE <= 0.01 &
           DECOY == 0 &
           SOFTWARE != 'SPECTRONAUT_FILTERED' &
           #!grepl('M|m', PCM_ID) &
           ORGANISM %in% c('HUMAN', 'ECOLI', 'YEAST')]
combined[, .N, keyby=CONDITION_MIN2_QUAN_FDR]

# ---- fdr based ma plots ----
#ma_combined_fdr[SOFTWARE=="SPECTRONAUT" & PCM_ID=="DNGEFSHHDLAPALDTGTTEEDR_3"]
ma_combined_fdr <- combined[QUAN_SAMPLES_FDR >= 6 & #CHANGED from 2
                              Q_VALUE <= 0.01 &
                              DECOY == 0 &
                              SOFTWARE != 'SPECTRONAUT_FILTERED' &
                              #!grepl('M|m', PCM_ID) &
                              ORGANISM %in% c('HUMAN', 'ECOLI', 'YEAST'),
                            list(PCM_ID, SAMPLE, CONDITION, QUAN, QUAN_FRAGS,
                                 QUAN_COMPLETE_EFDR, QUAN_MIN1_EFDR,
                                 QUAN_SAMPLES_ENTRAPMENT_FDR,
                                 QUAN_MIN1COND_EFDR, ORGANISM, SOFTWARE)]
max(ma_combined_fdr[,table(PCM_ID), by = SOFTWARE][,V1])
ma_combined_fdr[,uniqueN(PCM_ID), by = .(CONDITION, SOFTWARE)]

ma_combined_fdr <- unique(ma_combined_fdr)
any(ma_combined_fdr$QUAN == 0, na.rm=T)
any(is.na(ma_combined_fdr$QUAN))

ma_combined_fdr <- ma_combined_fdr[,.(MEAN_INTENSITY=mean(QUAN),
                                      MIN_INTENSITY=min(QUAN, na.rm = T),
                                      MEAN_INTENSITY_A=mean(QUAN[CONDITION=='A']),
                                      MEAN_INTENSITY_B=mean(QUAN[CONDITION=='B']),
                                      MIN_INTENSITY_A=min(QUAN[CONDITION=='A']),
                                      MIN_INTENSITY_B=min(QUAN[CONDITION=='B']),
                                      SD_INTENSITY_A=sd(QUAN[CONDITION=='A']),
                                      SD_INTENSITY_B=sd(QUAN[CONDITION=='B']),
                                      LOG2RATIO=log2(mean(QUAN[CONDITION=='A'])/mean(QUAN[CONDITION=='B']))),
                                   by=list(SOFTWARE, PCM_ID, ORGANISM, QUAN_COMPLETE_EFDR, QUAN_MIN1_EFDR, QUAN_SAMPLES_ENTRAPMENT_FDR, QUAN_MIN1COND_EFDR)]
any(is.na(ma_combined_fdr$LOG2RATIO))

# ---- intermediate save ----
write.fst(ma_combined_fdr, file.path(figurePath, 'intermediate/20241127_figure4d_ma_noNorm_efdr_localPcmGrouper_pepFasta_min1eFdr_keepMet.fst'), compress = 100)

dtOrg <- copy(ma_combined_fdr)
softwareLevels <- c("CHIMERYS", "DIA-NN", "SPECTRONAUT", "SPECTRONAUT_FILTERED")
softwareLabels <- c("CHIMERYS", "DIA-NN", "Spectronaut", "Spectronaut\n(curated)")
dtOrg[, SOFTWARE := factor(SOFTWARE, softwareLevels, softwareLabels)]
organismLevels <- c("YEAST", "HUMAN", "ECOLI")
organismLabels <- c("Yeast", "Human", "E. coli")
organismRatios <- setNames(log2(c(2, 1, 0.25)), organismLabels)
dtOrg[, ORGANISM := factor(ORGANISM, organismLevels, organismLabels)]
dtOrg[, eFdrLabelComp := ifelse(QUAN_COMPLETE_EFDR, "eFDR all ≤ 1%", "eFDR min 1 > 1%")]
dtOrg[, eFdrLabelComp := factor(eFdrLabelComp, c("eFDR all ≤ 1%", "eFDR min 1 > 1%"))]
dtOrg[, eFdrLabelCond := ifelse(QUAN_MIN1COND_EFDR, "eFDR min 1 per\ncondition ≤ 1%", "eFDR all in any\ncondition > 1%")]
dtOrg[, eFdrLabelCond := factor(eFdrLabelCond, c("eFDR min 1 per\ncondition ≤ 1%", "eFDR all in any\ncondition > 1%"))]
dtMaLines <- data.table(YINTERCEPT = organismRatios, ORGANISM = factor(organismLabels))

fwrite(dtOrg[, .(SOFTWARE, LOG2RATIO, ORGANISM, eFdrLabelCond)],
       file.path(figurePath, "figure-4D-density-keepMet.csv"))


# ---- re-read intermediately saved data ----
# ma_combined_fdr <- read.fst( file.path(figurePath, 'intermediate/20241127_figure4d_ma_noNorm_efdr_localPcmGrouper_pepFasta_min1eFdr_keepMet.fst'), as.data.table = T)

# ---- plot ----

ggplot(ma_combined_fdr,
       mapping = aes(x = log10(MEAN_INTENSITY), y = LOG2RATIO, col = ORGANISM)) +
  geom_point(shape = '.', alpha = 0.5) +
  facet_wrap(.~SOFTWARE,
             strip.position = 'bottom',
             scales = 'free') +
  geom_hline(yintercept = log2(1), lty = 4, linewidth = 0.5) +
  geom_hline(yintercept = log2(2), lty = 4, linewidth = 0.5) +
  geom_hline(yintercept = log2(0.25), lty = 4, linewidth = 0.5) +
  theme(text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12 ),
        axis.ticks = element_line(size = .5),
        axis.title = element_text(size = 12),
        legend.position = 'top',
        legend.text=element_text(size = 12),
        legend.margin = margin(c(0,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5, ),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside") +
  guides(color = guide_legend(
    override.aes=list(shape = 19, size = 2))) +
  scale_x_continuous(expand = c(0.01, 0),
                     limits = c(3, 10)) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(-6, 6)) +
  labs(x = 'log10(Mean intensity, MS2-based)',
       y = 'Foldchange [A/B log2]',
       color = 'Organism')

ggplot(ma_combined_fdr,
       mapping = aes(x = log10(MEAN_INTENSITY), y = LOG2RATIO, col = ORGANISM)) +
  geom_point(shape = '.', alpha = 0.5) +
  facet_wrap(.~SOFTWARE,
             strip.position = 'bottom',
             scales = 'free') +
  geom_hline(yintercept = log2(1), lty = 4, linewidth = 0.5) +
  geom_hline(yintercept = log2(2), lty = 4, linewidth = 0.5) +
  geom_hline(yintercept = log2(0.25), lty = 4,linewidth = 0.5) +
  theme(text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12 ),
        axis.ticks = element_line(size = .5),
        axis.title = element_text(size = 12),
        legend.position = 'top',
        legend.text=element_text(size = 12),
        legend.margin = margin(c(0,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5, ),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside") +
  guides(color = guide_legend(
    override.aes=list(shape = 19, size = 2))) +
  scale_x_continuous(expand = c(0.01, 0),
                     limits = c(0, 10)) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(-15, 15)) +
  labs(x = 'log10(Mean intensity, MS2-based)',
       y = 'Foldchange [A/B log2]',
       color = 'Organism')
