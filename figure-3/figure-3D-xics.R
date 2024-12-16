# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-3")
figurePath <- file.path(dataPath, "data/figure-3")
#xicPath <- file.path(dataPath, "figure-3/xic")

# ---- pathsToData ----
pathToExport <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Spectronaut/20240925_064248_20240925_SN19_lfq_paper_entrapment_paper_Report_height_noNorm.tsv")

# ---- load xics ----
pathToElutionCurves_empty <- file.path(dataPath, "data/figure-3/intermediate/_GLDDESGPTHGNDSGNHR_.4.tsv")
pathToElutionCurves_good <- file.path(dataPath, "data/figure-3/intermediate/_TIGGGDDSFNTFFSETGAGK_.2.tsv")

# ---- re-load intermediately saved data; including spectronaut and chimerys ----
combined <- read.fst(file.path(dataPath, 'data/figure-3/20241127_figure3a_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst'), as.data.table = T)

# ---- subset per software ----
pdresult <- combined[SOFTWARE == 'CHIMERYS']
sn19_flaggedImputation <- combined[SOFTWARE == 'SPECTRONAUT']
sn19_report <- fread(pathToExport, stringsAsFactors = F, integer64 = 'double')

# ---- get example good / empty meta values
names(sn19_report)[grepl("PEP", names(sn19_report))]
sn19_report[EG.PrecursorId=="_GLDDESGPTHGNDSGNHR_.4" & R.FileName=="LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03",
            .(EG.Qvalue, EG.PEP, EG.Cscore, EG.StartRT, EG.EndRT)]
sn19_report[EG.PrecursorId=="_TIGGGDDSFNTFFSETGAGK_.2" & R.FileName=="LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03",
            .(EG.Qvalue, EG.PEP, EG.Cscore, EG.StartRT, EG.EndRT)]

# --- load xics ----
empty_candidate <- fread(pathToElutionCurves_empty, stringsAsFactors = F, integer64 = 'double')
good_candidate <- fread(pathToElutionCurves_good, stringsAsFactors = F, integer64 = 'double')

# --- prepare xics for plotting ----
ions_empty_xic <- colnames(empty_candidate)[grepl('Intensity', colnames(empty_candidate))]
ions_empty_xic <- strex::str_after_first(ions_empty_xic, '_')
ions_empty_xic <- ions_empty_xic[!grepl('Predicted', ions_empty_xic)]

empty_xic <- data.table()
for(i in 1:length(ions_empty_xic)){
  # i = 1
  ions_empty_xic.i <- ions_empty_xic[i]
  intensities.i <- colnames(empty_candidate)[grepl('Intensity', colnames(empty_candidate)) &
                                                   grepl(ions_empty_xic.i, colnames(empty_candidate), fixed = T)]
  rts.i <- colnames(empty_candidate)[grepl('Retention Time', colnames(empty_candidate)) &
                                           grepl(ions_empty_xic.i, colnames(empty_candidate), fixed = T)]
  cols.i <- c(intensities.i, rts.i)
  ion.i <- empty_candidate[,..cols.i]
  melt(ion.i)
  ion.i[,ion.i:=ions_empty_xic.i]
  setnames(ion.i, c('Intensity', 'Retention time', 'Fragment ion'))
  ion.i <- ion.i[,list(`Fragment ion`, `Retention time`, Intensity)]
  empty_xic <- rbind(empty_xic, ion.i)
}

ions_good_xic <- colnames(good_candidate)[grepl('Intensity', colnames(good_candidate))]
ions_good_xic <- strex::str_after_first(ions_good_xic, '_')
ions_good_xic <- ions_good_xic[!grepl('Predicted', ions_good_xic)]

good_xic <- data.table()
for(i in 1:length(ions_good_xic)){
  # i = 1
  ions_good_xic.i <- ions_good_xic[i]
  intensities.i <- colnames(good_candidate)[grepl('Intensity', colnames(good_candidate)) &
                                               grepl(ions_good_xic.i, colnames(good_candidate), fixed = T)]
  rts.i <- colnames(good_candidate)[grepl('Retention Time', colnames(good_candidate)) &
                                       grepl(ions_good_xic.i, colnames(good_candidate), fixed = T)]
  cols.i <- c(intensities.i, rts.i)
  ion.i <- good_candidate[,..cols.i]
  melt(ion.i)
  ion.i[,ion.i:=ions_good_xic.i]
  setnames(ion.i, c('Intensity', 'Retention time', 'Fragment ion'))
  ion.i <- ion.i[,list(`Fragment ion`, `Retention time`, Intensity)]
  good_xic <- rbind(good_xic, ion.i)
}

# add EG.PrecursorId and R.FileName information
empty_xic[,EG.PrecursorId:="_GLDDESGPTHGNDSGNHR_.4"]
empty_xic[,R.FileName:="LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03"]

good_xic[,EG.PrecursorId:="_TIGGGDDSFNTFFSETGAGK_.2"]
good_xic[,R.FileName:="LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03"]

# ---- add qvalue ----
empty_xic <- merge(empty_xic,
                   sn19_flaggedImputation[,list(ID, SAMPLE,
                                                Q_VALUE, ENTRAPMENT_Q_VALUE,
                                                SVMSCORE, DECOY, QUAN,
                                                MIN1_ID_FRAGS, MIN1_QUAN_FRAGS)],
                   by.x = c('EG.PrecursorId', 'R.FileName'),
                   by.y = c('ID', 'SAMPLE'),
                   all.x = T)
unique(empty_xic$EG.PrecursorId)
unique(empty_xic$R.FileName)
unique(empty_xic$SVMSCORE)
unique(empty_xic$Q_VALUE)
unique(empty_xic$ENTRAPMENT_Q_VALUE)
unique(empty_xic$EG.DatapointsPerPeak)

good_xic <- merge(good_xic,
                   sn19_flaggedImputation[,list(ID, SAMPLE,
                                                Q_VALUE, ENTRAPMENT_Q_VALUE,
                                                SVMSCORE, DECOY, QUAN,
                                                MIN1_ID_FRAGS, MIN1_QUAN_FRAGS)],
                   by.x = c('EG.PrecursorId', 'R.FileName'),
                   by.y = c('ID', 'SAMPLE'),
                   all.x = T)
unique(good_xic$EG.PrecursorId)
unique(good_xic$R.FileName)
unique(good_xic$SVMSCORE)
unique(good_xic$Q_VALUE)
unique(good_xic$ENTRAPMENT_Q_VALUE)
unique(good_xic$EG.DatapointsPerPeak)

# ---- sort by maximum abundance of fragments ----
setnames(empty_xic, c('Retention time', 'Fragment ion'), c('RT', 'IonLabel'))
empty_xic[,APEX_INTENSITY:=max(Intensity), by = IonLabel]
empty_xic <- empty_xic[order(APEX_INTENSITY, decreasing = T)]

setnames(good_xic, c('Retention time', 'Fragment ion'), c('RT', 'IonLabel'))
good_xic[,APEX_INTENSITY:=max(Intensity), by = IonLabel]
good_xic <- good_xic[order(APEX_INTENSITY, decreasing = T)]

# ---- add fragment ion information ----
table(sn19_report$F.Charge)
any(is.nan(sn19_report$F.Charge))
sn19_report[,F.ChargeXicDb:=sapply(F.Charge, function(x) paste0(rep('+', x), collapse = ''))]
sn19_report[,IonLabel:=ifelse(F.FrgLossType == 'noloss',
                              paste0(F.FrgIon, F.ChargeXicDb),
                              paste0(F.FrgIon, F.ChargeXicDb, ' -', F.FrgLossType))]

# ---- information on empty xic and good XIC ----
unique(empty_xic$EG.PrecursorId)
unique(empty_xic$R.FileName)
empty_xic[,IonLabel_reportFormat:=gsub('[', '', IonLabel, fixed = T)]
empty_xic[,IonLabel_reportFormat:=gsub(']', '', IonLabel_reportFormat, fixed = T)]
empty_xic[,IonLabel_reportFormat:=sub("\\ -.*", "", IonLabel_reportFormat)]

empty_xic <- merge(empty_xic,
                   sn19_report[EG.PrecursorId == '_GLDDESGPTHGNDSGNHR_.4' &
                                 R.FileName == 'LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03',
                               list(EG.PrecursorId, R.FileName,
                                    EG.PEP,
                                    FG.PrecMz,
                                    IonLabel, F.FrgIon, F.FrgLossType,
                                    F.FrgMz, F.Charge, F.PeakArea,
                                    EG.StartRT, EG.EndRT)],
                   by.x = c('EG.PrecursorId', 'R.FileName', 'IonLabel_reportFormat'),
                   by.y = c('EG.PrecursorId', 'R.FileName', 'IonLabel'),
                   all.x = T)

empty_xic_info <- empty_xic[,list(EG.PrecursorId, R.FileName,
                                  EG.PEP,
                                  EG.StartRT, EG.EndRT,
                                  FG.PrecMz,
                                  F.FrgIon, F.FrgLossType,
                                  F.FrgMz, F.Charge, F.PeakArea)]
empty_xic_info <- unique(empty_xic_info)

unique(good_xic$EG.PrecursorId)
unique(good_xic$R.FileName)
good_xic[,IonLabel_reportFormat:=gsub('[', '', IonLabel, fixed = T)]
good_xic[,IonLabel_reportFormat:=gsub(']', '', IonLabel_reportFormat, fixed = T)]
good_xic[,IonLabel_reportFormat:=sub("\\ -.*", "", IonLabel_reportFormat)]

good_xic <- merge(good_xic,
                   sn19_report[EG.PrecursorId == '_TIGGGDDSFNTFFSETGAGK_.2' &
                                 R.FileName == 'LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03',
                               list(EG.PrecursorId, R.FileName,
                                    EG.PEP,
                                    FG.PrecMz,
                                    IonLabel, F.FrgIon, F.FrgLossType,
                                    F.FrgMz, F.Charge, F.PeakArea,
                                    EG.StartRT, EG.EndRT)],
                   by.x = c('EG.PrecursorId', 'R.FileName', 'IonLabel_reportFormat'),
                   by.y = c('EG.PrecursorId', 'R.FileName', 'IonLabel'),
                   all.x = T)

# ---- intermediate save for AH and MF ----
write.fst(empty_xic_info, file.path(dataPath, 'data/figure-3/fst-backup/20241127_figure3d_exampleData_emptyXic_apexQuan_pepEntr.fst'), compress = 100)
write.fst(empty_xic, file.path(dataPath, 'data/figure-3/fst-backup/20241127_figure3d_example_emptyXic_apexQuan_pepEntr.fst'), compress = 100)
write.fst(good_xic, file.path(dataPath, 'data/figure-3/fst-backup/20241127_figure3d_example_goodXic_apexQuan_pepEntr.fst'), compress = 100)

fwrite(good_xic, file.path(figurePath, "figure-3D-xic-good.csv"))
fwrite(empty_xic, file.path(figurePath, "figure-3D-xic-bad.csv"))

# ---- re-read intermediately saved data ----
# empty_xic <- read.fst(file.path(dataPath, 'data/figure-3/fst-backup/20241127_figure3d_example_emptyXic_apexQuan_pepEntr.fst'), as.data.table = T)
# good_xic <- read.fst(file.path(dataPath, 'data/figure-3/fst-backup/20241127_figure3d_example_goodXic_apexQuan_pepEntr.fst'), as.data.table = T)

# ---- preliminary plots ----
ggplot(empty_xic, aes(RT, Intensity, group = IonLabel, color = IonLabel)) +
  geom_point() +
  geom_line() +
  theme(text = element_text(size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(size = .5),
        legend.margin = margin(c(0,0,0,0)),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        plot.title=element_text(size=12, face = 'plain')) +
  scale_y_continuous(expand = c(0.02,0)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  labs(color = 'Fragment ions',
       title = paste0('EG.QValue ', unique(empty_xic$Q_VALUE),'\n',
                      'CScore ', unique(empty_xic$SVMSCORE),'\n',
                      'EG.PEP', unique(empty_xic$EG.PEP)))

ggplot(good_xic, aes(RT, Intensity, group = IonLabel, color = IonLabel)) +
  geom_point() +
  geom_line() +
  theme(text = element_text(size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(size = .5),
        legend.margin = margin(c(0,0,0,0)),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        plot.title=element_text(size=12, face = 'plain')) +
  scale_y_continuous(expand = c(0.02,0)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  labs(color = 'Fragment ions',
       title = paste0('EG.QValue ', unique(good_xic$Q_VALUE),'\n',
                      'CScore ', unique(good_xic$SVMSCORE)))
