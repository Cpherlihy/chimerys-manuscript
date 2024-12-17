# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-4")
figurePath <- file.path(dataPath, "data/figure-4")
#cvsPath <- file.path(dataPath, "figure-4/cvs")

# ---- pathsToData ----
## combined
pathToCombined <- file.path(dataPath, 'data/figure-3/20241127_figure3a_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst')

# ---- read data ----
combined <- read.fst(pathToCombined, as.data.table = T)

# ---- for peptide-centric search engines, replace ENTRAPMENT_Q_VALUE with ENTRAPMENT_Q_VALUE_1 ----
combined[SOFTWARE=="CHIMERYS", ENTRAPMENT_Q_VALUE_FILTER := ENTRAPMENT_Q_VALUE]
combined[SOFTWARE!="CHIMERYS", ENTRAPMENT_Q_VALUE_FILTER := ENTRAPMENT_Q_VALUE_1]

# ---- flag quan complete ----
combined[,QUAN_MIN1_EFDR:=QUAN_SAMPLES_ENTRAPMENT_FDR>=1]
combined[,QUAN_MIN1COND_EFDR:=all(QUAN_CONDITION_ENTRAPMENT_FDR>=1), by=.(SOFTWARE, PCM_ID)]

# ---- intermediate save ----
write.fst(combined, file.path(figurePath, 'intermediate/20241127_figure4c_combined_pcms_localPcmGrouper_apexQuan_pepEntr.fst'), compress = 100)

# ---- re-read intermediately saved data ----
# combined <- read.fst(file.path(figurePath, 'intermediate/20241127_figure4c_combined_pcms_localPcmGrouper_apexQuan_pepEntr.fst'), as.data.table = T)

# ---- remove MBR for CHIMERYS ----
table(combined$IS_IDENTIFIED_BY_MBR)
min(combined[SOFTWARE == 'CHIMERYS', QUAN])
combined[SOFTWARE=='CHIMERYS' & IS_IDENTIFIED_BY_MBR==1,QUAN:=NA]

# ---- fdr based cvs ----
#combined[SOFTWARE=="CHIMERYS" & PCM_ID=="VLHEAEGHIVTcETNTGEVYR_4"]
cv_combined_fdr <- combined[QUAN_COMPLETE_FDR == TRUE &
                              Q_VALUE <= 0.01 &
                              DECOY == 0 &
                              SOFTWARE != 'SPECTRONAUT_FILTERED' &
                              !grepl('M|m', PCM_ID) &
                              ORGANISM %in% c('HUMAN', 'ECOLI', 'YEAST'),
                            list(PCM_ID, SAMPLE, CONDITION, QUAN_MIN1COND_EFDR, QUAN, QUAN_FRAGS, ORGANISM, SOFTWARE)]
#552,630 - now 657,222
max(cv_combined_fdr[,table(PCM_ID), by = SOFTWARE][,V1])
cv_combined_fdr[,uniqueN(PCM_ID), by = .(CONDITION, SOFTWARE)]

cv_combined_fdr <- unique(cv_combined_fdr)
any(cv_combined_fdr$QUAN == 0)
any(is.na(cv_combined_fdr$QUAN))

cv_combined_fdr <- cv_combined_fdr[,list(CV=sd(QUAN, na.rm=T)/mean(QUAN, na.rm=T)*100,
                                         QUAN_FRAGS=max(QUAN_FRAGS)),
                                   by=list(PCM_ID, CONDITION, QUAN_MIN1COND_EFDR, ORGANISM, SOFTWARE)]
any(is.na(cv_combined_fdr$CV))

# ---- count cvs; fdr ----
cv_combined_fdr[,LABEL:=ifelse(QUAN_FRAGS >= 3, '≥ 3', QUAN_FRAGS)]
cv_combined_fdr[,LABEL:=factor(LABEL,
                               levels = c('0', '1', '2', '≥ 3'))]

cv_combined_fdr[,COUNT:=.N, by = .(SOFTWARE, LABEL)]

# ---- intermediate save ----
write.fst(cv_combined_fdr, file.path(figurePath, 'intermediate/20241127_figure4c_cvs_noNorm_fdr_localPcmGrouper_pepFasta.fst'), compress = 100)

#dtCv[!is.na(CONDITION), .N, keyby=SOFTWARE]
dtCv <- rbind(cbind(TYPE = "no eFDR filter", cv_combined_fdr),
              cbind(TYPE = "eFDR min 1 per\ncondition ≤ 1%", cv_combined_fdr[QUAN_MIN1COND_EFDR==T]))

cvEfdrLabel <- c("no eFDR filter", "eFDR min 1 per\ncondition ≤ 1%")
dtCv[, TYPE := factor(TYPE, cvEfdrLabel)]
dtCv[, COUNT := NULL]
dtCv[, CV := CV/100]
softwareLevels <- c("CHIMERYS", "DIA-NN", "SPECTRONAUT", "SPECTRONAUT_FILTERED")
softwareLabels <- c("CHIMERYS", "DIA-NN", "Spectronaut", "Spectronaut\n(curated)")
dtCv[, SOFTWARE := factor(SOFTWARE, softwareLevels, softwareLabels)]
setkey(dtCv, TYPE, SOFTWARE, LABEL)

fwrite(dtCv, file.path(figurePath, "figure-4C-CVs.csv"))

# ---- re-read intermediately saved data ----
# cv_combined_fdr <- read.fst(file.path(figurePath, 'intermediate/20241127_figure4c_cvs_noNorm_fdr_localPcmGrouper_pepFasta.fst'), as.data.table = T)

# ---- preliminary plot ----
ggplot(cv_combined_fdr,
       aes(x = LABEL,
           y = CV,
           fill = LABEL)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~SOFTWARE,
             strip.position = 'bottom') +
  geom_text(data = cv_combined_fdr,
            aes(LABEL, 180,
                label = COUNT), size = 4) +
  theme(text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(linewidth = .5),
        legend.margin = margin(c(0,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside") +
  scale_x_discrete(drop=FALSE) +
  labs(y = 'CV [%]',
       x = '# fragments for quantification',
       fill = '# fragments for quantification')

# ---- efdr based cvs ----
cv_combined_efdr <- combined[QUAN_COMPLETE_FDR == TRUE &
                               QUAN_MIN1COND_EFDR == TRUE & #CHANGED from QUAN_COMPLETE_EFDR
                               Q_VALUE <= 0.01 &
                               DECOY == 0 &
                               SOFTWARE != 'SPECTRONAUT_FILTERED' &
                               !grepl('M|m', PCM_ID) &
                               ORGANISM %in% c('HUMAN', 'ECOLI', 'YEAST'),
                             list(PCM_ID, SAMPLE, CONDITION, QUAN_MIN1COND_EFDR, QUAN, QUAN_FRAGS, ORGANISM, SOFTWARE)]
#657,222
max(cv_combined_efdr[,table(PCM_ID), by = SOFTWARE][,V1])
cv_combined_efdr[,uniqueN(PCM_ID), by = .(CONDITION, SOFTWARE)]

cv_combined_efdr <- unique(cv_combined_efdr)
any(cv_combined_efdr$QUAN == 0)
any(is.na(cv_combined_efdr$QUAN))

cv_combined_efdr <- cv_combined_efdr[,list(CV=sd(QUAN, na.rm=T)/mean(QUAN, na.rm=T)*100,
                                           QUAN_FRAGS=max(QUAN_FRAGS)),
                                     by=list(PCM_ID, CONDITION, QUAN_MIN1COND_EFDR, ORGANISM, SOFTWARE)]
any(is.na(cv_combined_efdr$CV))

test2 <- combined[PCM_ID == 'HEAAEALGAIASPEVVDVLK_3' & SOFTWARE == 'DIA-NN' &
                    CONDITION == 'A' &
                    CONDITION_MIN2_QUAN_EFDR == 2 &
                    QUAN_CONDITION_ENTRAPMENT_FDR >= 2 &
                    Q_VALUE <= 0.01 &
                    ENTRAPMENT_Q_VALUE <= 0.01 &
                    DECOY == 0 &
                    SOFTWARE != 'SPECTRONAUT_FILTERED' &
                    !grepl('M|m', PCM_ID) &
                    ORGANISM %in% c('HUMAN', 'ECOLI', 'YEAST')]

test2 <- combined[PCM_ID == 'HEAAEALGAIASPEVVDVLK_3' & SOFTWARE == 'DIA-NN']
max(test2$ENTRAPMENT_Q_VALUE)
test2[,list(CV=sd(QUAN, na.rm=T)/mean(QUAN, na.rm=T)*100,
            QUAN_FRAGS=max(QUAN_FRAGS)),
      by=list(PCM_ID, CONDITION, ORGANISM, SOFTWARE)]

# ---- count cvs; fdr ----
cv_combined_efdr[,LABEL:=ifelse(QUAN_FRAGS >= 3, '≥ 3', QUAN_FRAGS)]
cv_combined_efdr[,LABEL:=factor(LABEL,
                                levels = c('0', '1', '2', '≥ 3'))]

cv_combined_efdr[,COUNT:=.N, by = .(SOFTWARE, LABEL)]

# ---- intermediate save ----
write.fst(cv_combined_efdr, file.path(figurePath, 'intermediate/20241127_figure4c_cvs_noNorm_efdr_localPcmGrouper_pepFasta.fst'), compress = 100)

# ---- re-read intermediately saved data ----
# cv_combined_efdr <- read.fst(file.path(figurePath, 'intermediate/20241127_figure4c_cvs_noNorm_efdr_localPcmGrouper_pepFasta.fst'), as.data.table = T)

# ---- plot ----
ggplot(cv_combined_efdr,
       aes(x = LABEL,
           y = CV,
           fill = LABEL)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  # ,scale = "width") +
  facet_wrap(~SOFTWARE,
             strip.position = 'bottom') +
  geom_text(data = cv_combined_efdr,
            aes(LABEL, 180,
                label = COUNT), size = 4) +
  theme(text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(linewidth = .5),
        legend.margin = margin(c(0,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside") +
  scale_x_discrete(drop=FALSE) +
  labs(y = 'CV [%]',
       x = '# fragments for quantification',
       fill = '# fragments for quantification')
