# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-E7")
figurePath <- file.path(dataPath, "data/figure-E7")

# ---- pathsToData ----
## dia-nn
pathToTsv <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN/20240429_lfq_height_entrapment_peptides_report.tsv")
pathToLib <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN/20240429_lfq_height_entrapment_peptides_report-lib.tsv")

## spectronaut
pathToExport <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/Spectronaut-entrapment/20240829_075057_20240829_SN19_lfq_paper_entrapment_peptides_Report_height_noNorm.tsv")

## combined generated in figure-3/figure-3A-entrapment-barplot.R
pathToCombined <- file.path(dataPath, 'data/figure-3/fst-backup/20241127_figure3a_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst')

# ---- load the data ----
diann <- fread(pathToTsv, stringsAsFactors = F, integer64 = 'double')
diann_lib <- fread(pathToLib, stringsAsFactors = F, integer64 = 'double')

sn <- fread(pathToExport, stringsAsFactors = F, integer64 = 'double')

# ---- load data including entrapment fdr ----
combined <- read.fst(pathToCombined, as.data.table = T)

sn_eFDR <- combined[SOFTWARE == 'SPECTRONAUT']
sn_filtered_eFDR <- combined[SOFTWARE == 'SPECTRONAUT_FILTERED']
diann_eFDR <- combined[SOFTWARE == 'DIA-NN']

# ---- format diann to access fragment quan values ----
strex::str_after_last(diann$File.Name[1], '\\\\')
diann[,File.Name:=strex::str_after_last(File.Name, '\\\\')]
diann[,File.Name:=gsub('.raw', '', File.Name)]

strex::str_after_last(diann_lib$FileName[1], '\\\\')
diann_lib[,FileName:=strex::str_after_last(FileName, '\\\\')]
diann_lib[,FileName:=gsub('.raw', '', FileName)]

diann_lib[,Precursor.Id:=strex::str_before_first(transition_name, '_')]
diann_lib[,ProductMz:=round(ProductMz, 4)]

# ---- fragment id for predictions -----
diann_lib[,Fragment.Id:=seq(1, .N), by = Precursor.Id]

# --- predicted fragments in lib per precursor -----
diann_lib[,DIANN_PRED_FRAGMENTS:=.N, by = Precursor.Id]

# --- ExcludedFromAssay fragments in lib per precursor -----
diann_lib[,DIANN_ASSAY_FRAGMENTS:=uniqueN(Fragment.Id[ExcludeFromAssay == FALSE]),by = Precursor.Id]

# ---- raw quan values per Precursor.Id and File.Name ----
diann[, F.QUAN_VALUES := sapply(Fragment.Quant.Raw, function(x) length(which(strsplit(x, ";")[[1]]>0)))]

# ---- corrected quan values per Precursor.Id and File.Name ----
diann[, F.CORRECTED_QUAN_VALUES := sapply(Fragment.Quant.Corrected, function(x) length(which(strsplit(x, ";")[[1]]>0)))]

# ---- add counts to diann report ----
diann <- merge(diann,
               unique(diann_lib[, list(Precursor.Id, DIANN_PRED_FRAGMENTS, DIANN_ASSAY_FRAGMENTS)]),
               by = 'Precursor.Id')

# ---- add fragment intensities to predictions ----
# NOTE: only max 12 quan values; only added for Precursor.Ids with DIANN_PRED_FRAGMENTS <= 12
fragment_quant_raw <- diann[, list(Fragment.Quant.Raw = unlist(strsplit(Fragment.Quant.Raw, ";"))), by=.(Precursor.Id, File.Name, Precursor.Quantity)]
fragment_quant_raw[,Fragment.Quant.Raw:=as.numeric(Fragment.Quant.Raw)]
fragment_quant_raw[,Fragment.Id:=seq(1, .N), by = .(Precursor.Id, File.Name, Precursor.Quantity)]

fragment_quant_corrected <- diann[, list(Fragment.Quant.Corrected = unlist(strsplit(Fragment.Quant.Corrected, ";"))), by=.(Precursor.Id, File.Name, Precursor.Quantity)]
fragment_quant_corrected[,Fragment.Quant.Corrected:=as.numeric(Fragment.Quant.Corrected)]
fragment_quant_corrected[,Fragment.Id:=seq(1, .N), by = .(Precursor.Id, File.Name, Precursor.Quantity)]

# ---- add entrapment FDR to fragments ----
## diann
table(diann$File.Name)
table(diann_eFDR$SAMPLE)

setkey(diann_eFDR, ID, SAMPLE)
setkey(diann, Precursor.Id, File.Name)
setkey(fragment_quant_corrected, Precursor.Id, File.Name)

fragment_quant_corrected <- fragment_quant_corrected[diann_eFDR[,list(ID, SAMPLE, Q_VALUE, ENTRAPMENT_Q_VALUE_1)]]

## spectronaut
table(sn$R.FileName)
table(sn_eFDR$SAMPLE)

setkey(sn_eFDR, ID, SAMPLE)
setkey(sn, EG.PrecursorId, R.FileName)

sn <- sn[sn_eFDR[,list(ID, SAMPLE, Q_VALUE, ENTRAPMENT_Q_VALUE_1)]]

## spectronaut filtered
table(sn$R.FileName)
table(sn_filtered_eFDR$SAMPLE)

setkey(sn_filtered_eFDR, ID, SAMPLE)
setkey(sn, EG.PrecursorId, R.FileName)

sn_filtered <- sn[sn_filtered_eFDR[,list(ID, SAMPLE, Q_VALUE, ENTRAPMENT_Q_VALUE_1)]]
sn_filtered <- sn_filtered[F.PeakArea > 1]

# ---- combine ----
fragment_quan_sn <- sn[,list(EG.PrecursorId, R.FileName, Q_VALUE, ENTRAPMENT_Q_VALUE_1, F.PeakArea)]
setnames(fragment_quan_sn, c('PRECURSOR', 'SAMPLE', 'Q_VALUE', 'ENTRAPMENT_Q_VALUE_1', 'FRAGMENT_QUAN'))
any(is.na(fragment_quan_sn$Q_VALUE))
str(fragment_quan_sn$Q_VALUE)

fragment_quan_sn_filtered <- sn_filtered[,list(EG.PrecursorId, R.FileName, Q_VALUE, ENTRAPMENT_Q_VALUE_1, F.PeakArea)]
setnames(fragment_quan_sn_filtered, c('PRECURSOR', 'SAMPLE', 'Q_VALUE', 'ENTRAPMENT_Q_VALUE_1', 'FRAGMENT_QUAN'))
any(is.na(fragment_quan_sn_filtered$Q_VALUE))
str(fragment_quan_sn_filtered$Q_VALUE)

fragment_quan_diann <- fragment_quant_corrected[,list(Precursor.Id, File.Name, Q_VALUE, ENTRAPMENT_Q_VALUE_1, Fragment.Quant.Corrected)]
setnames(fragment_quan_diann, c('PRECURSOR', 'SAMPLE', 'Q_VALUE', 'ENTRAPMENT_Q_VALUE_1', 'FRAGMENT_QUAN'))
any(is.na(fragment_quan_diann$Q_VALUE))
str(fragment_quan_diann$Q_VALUE)

fragment_quan_combined <- rbind(fragment_quan_sn[,SOFTWARE:='SPECTRONAUT'],
                                fragment_quan_sn_filtered[,SOFTWARE:='SPECTRONAUT_FILTERED'],
                                fragment_quan_diann[,SOFTWARE:='DIA-NN'])

fragment_quan_combined[FRAGMENT_QUAN>0, LOG10QUAN := log10(FRAGMENT_QUAN)]
softwareLevels <- c("CHIMERYS", "DIA-NN", "SPECTRONAUT", "SPECTRONAUT_FILTERED")
softwareLabels <- c("CHIMERYS", "DIA-NN", "Spectronaut", "Spectronaut\n(curated)")
fragment_quan_combined[, SOFTWARE := factor(SOFTWARE, softwareLevels[2:4], softwareLabels[2:4])]
fragment_quan_combined[, isEfdr001 := factor(ENTRAPMENT_Q_VALUE_1 <= 0.01, c(T, F))]
fragment_quan_combined <- fragment_quan_combined[!is.na(LOG10QUAN) & SOFTWARE %in% softwareLabels[2:3]]

# ---- intermediate save ----
fwrite(fragment_quan_combined[, .SD, .SDcols = c("SOFTWARE", "LOG10QUAN", "isEfdr001")],
       file.path(figurePath, "figure-E7C-fragments.csv"))
