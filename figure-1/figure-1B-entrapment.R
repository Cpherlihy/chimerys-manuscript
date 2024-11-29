#setup
source(here::here("scripts/load-dependencies.R"))
path <- file.path(here::here(), "figure-1")

# ---- pathsToData ----
## chimerys
pathToPdResult_IW1 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW1-(1).pdResult')
pathToPdResult_IW3 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW3-(1).pdResult')
pathToPdResult_IW6 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW6-(1).pdResult')
pathToPdResult_IW8 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW8-(1).pdResult')
pathToPdResult_IW10 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW10-(1).pdResult')
pathToPdResult_IW12 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW12-(1).pdResult')
pathToPdResult_IW15 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW15-(1).pdResult')
pathToPdResult_IW20 <- file.path(dataPath, 'wwDDA/60min-figure2/Chimerys-entrapment/210129_wwDDA_10ug_60min_IW20-(1).pdResult')

## fastas
pathToFasta <- file.path(dataPath, 'FASTA/JoT_uniprot-proteome_musmusculus_review_canonical_20220623_9mimic_-Ie.fasta')


# ---- read data including short sanity checks ----
proteinsFromFasta <- readProteinsFromFastas(pathToFasta)
table(proteinsFromFasta$ORGANISM)
any(grepl('Cont_', proteinsFromFasta$FastaTitleLines))

pdresult_IW1 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW1, proteinsFromFasta)
any(is.na(pdresult_IW1$ORGANISM)) # FALSE
table(pdresult_IW1$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW1[Q_VALUE <= 0.01, PCM_J_ID]) # 59009
uniqueN(pdresult_IW1[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 58467

pdresult_IW3 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW3, proteinsFromFasta)
any(is.na(pdresult_IW3$ORGANISM)) # FALSE
table(pdresult_IW3$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW3[Q_VALUE <= 0.01, PCM_J_ID]) # 63994
uniqueN(pdresult_IW3[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 63408

pdresult_IW6 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW6, proteinsFromFasta)
any(is.na(pdresult_IW6$ORGANISM)) # FALSE
table(pdresult_IW6$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW6[Q_VALUE <= 0.01, PCM_J_ID]) # 58444
uniqueN(pdresult_IW6[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 57908

pdresult_IW8 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW8, proteinsFromFasta)
any(is.na(pdresult_IW8$ORGANISM)) # FALSE
table(pdresult_IW8$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW8[Q_VALUE <= 0.01, PCM_J_ID]) # 53685
uniqueN(pdresult_IW8[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 53191

pdresult_IW10 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW10, proteinsFromFasta)
any(is.na(pdresult_IW10$ORGANISM)) # FALSE
table(pdresult_IW10$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW10[Q_VALUE <= 0.01, PCM_J_ID]) # 49443
uniqueN(pdresult_IW10[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 48994

pdresult_IW12 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW12, proteinsFromFasta)
any(is.na(pdresult_IW12$ORGANISM)) # FALSE
table(pdresult_IW12$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW12[Q_VALUE <= 0.01, PCM_J_ID]) # 45923
uniqueN(pdresult_IW12[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 45495

pdresult_IW15 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW15, proteinsFromFasta)
any(is.na(pdresult_IW15$ORGANISM)) # FALSE
table(pdresult_IW15$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW15[Q_VALUE <= 0.01, PCM_J_ID]) # 41345
uniqueN(pdresult_IW15[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 40962

pdresult_IW20 <- readPdResult_psmGrouper_dda_noConsensusQuan(pathToPdResult_IW20, proteinsFromFasta)
any(is.na(pdresult_IW20$ORGANISM)) # FALSE
table(pdresult_IW20$ORGANISM) # ENTRAPMENT EXISTS
uniqueN(pdresult_IW20[Q_VALUE <= 0.01, PCM_J_ID]) # 35970
uniqueN(pdresult_IW20[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 35629

# ---- add identifier ----
pdresult_IW1[,SOFTWARE:='IW1']
pdresult_IW3[,SOFTWARE:='IW3']
pdresult_IW6[,SOFTWARE:='IW6']
pdresult_IW8[,SOFTWARE:='IW8']
pdresult_IW10[,SOFTWARE:='IW10']
pdresult_IW12[,SOFTWARE:='IW12']
pdresult_IW15[,SOFTWARE:='IW15']
pdresult_IW20[,SOFTWARE:='IW20']

# ---- combine data ----
combined <- rbindlist(list(pdresult_IW1,
                           pdresult_IW3,
                           pdresult_IW6,
                           pdresult_IW8,
                           pdresult_IW10,
                           pdresult_IW12,
                           pdresult_IW15,
                           pdresult_IW20), fill = T)
uniqueN(combined$SAMPLE)
uniqueN(combined$SOFTWARE)

# ---- one row per MODPEP_ID ----
combined[,CHARGE:=NULL]
combined[,PCM_ID:=NULL]
combined[,PCM_J_ID:=NULL]
combined[,MODIFIED_SEQUENCE:=NULL]
combined[,PrecursorAbundance:=NULL]
combined <- unique(combined)

table(combined$ORGANISM)
combined[,ORGANISM:=gsub('(^|;)+ENTRAPMENT', '', ORGANISM)]
combined[,ORGANISM:=gsub('^;', '', ORGANISM)]
combined[ORGANISM=='',ORGANISM:='ENTRAPMENT']
combined[DECOY==1,ORGANISM:='DECOY']
table(combined$ORGANISM)

combined[,ORGANISM_FASTA_J:=ORGANISM]
combined[,ENTRAPMENT_FASTA_J:=ENTRAPMENT]

# ---- calculate fdr and entrapment fdr ----
combined <- setEntrapmentFDR_local(combined[!is.na(SVMSCORE)])

# ---- remove decoys ----
combined <- combined[DECOY == 0]

# ---- plot fdr vs efdr ----
ggplot(data = combined,
       mapping = aes(x = Q_VALUE, y = ENTRAPMENT_Q_VALUE,
                     col = SOFTWARE)) +
  geom_abline(intercept = 0, slope = 1, col = 1, lty = 2) +
  geom_abline(intercept = 0, slope = 1.5, col = 'lightgray', lty = 2) +
  geom_abline(intercept = 0, slope = 2/3, col = 'lightgray', lty = 2) +
  geom_vline(xintercept = 0.01, col = 2, lty = 2) +
  geom_line() +
  theme(text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(linewidth = .5),
        legend.margin = margin(c(0,0,0,0)),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        strip.background = element_blank()) +
  labs(y = 'Entrapment FDR') +
  xlim(c(0,0.05)) +
  ylim(c(0,0.05))

# ---- intermediate save ----
combined <- combined[Q_VALUE<=0.05]
iwLevels <- paste0("IW", c(1, 3, 6, 8, 10, 12, 15, 20))
iwLabels <- as.character(c(1, 3, 6, 8, 10, 12, 15, 20))
combined[, SOFTWARE := factor(SOFTWARE, iwLevels, iwLabels)]
combined[, Q_VALUE_BIN := ceiling(Q_VALUE*10000)/10000]
mean_efdr <- combined[, .(ENTRAPMENT_Q_VALUE = mean(ENTRAPMENT_Q_VALUE),
                          ENTRAPMENT_Q_VALUE_1 = mean(ENTRAPMENT_Q_VALUE_1)),
                      by=.(SOFTWARE, Q_VALUE_BIN)]
setnames(mean_efdr, "Q_VALUE_BIN", "Q_VALUE")

write.fst(combined, file.path(dataPath, 'data/figure-1/figure-1B-entrapment.fst'), compress = 100)
