# ---- setup ----
source(here::here("scripts/data-processing-2.R"))
path <- file.path(here::here(), "figure-3")
figurePath <- file.path(dataPath, "data/figure-3")

# ---- pathsToData ----

## dia-nn
pathToTsv <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN/20240429_lfq_height_entrapment_peptides_report.tsv")
pathToLib <- file.path(dataPath, "LFQ_Bench_multispecies/DIA/DIA-NN/20240429_lfq_height_entrapment_peptides_report-lib.tsv")

# ---- load the data ----
diann <- fread(pathToTsv, stringsAsFactors = F, integer64 = 'double')
diann_lib <- fread(pathToLib, stringsAsFactors = F, integer64 = 'double')

# ---- formatting ----
strex::str_after_last(diann$File.Name[1], '\\\\')
diann[,File.Name:=strex::str_after_last(File.Name, '\\\\')]
diann[,File.Name:=gsub('.raw', '', File.Name)]

strex::str_after_last(diann_lib$FileName[1], '\\\\')
diann_lib[,FileName:=strex::str_after_last(FileName, '\\\\')]
diann_lib[,FileName:=gsub('.raw', '', FileName)]

diann_lib[,Precursor.Id:=strex::str_before_first(transition_name, '_')]
diann_lib[,ProductMz:=round(ProductMz, 4)]

# ---- fragment id for predictions -----
diann_lib <- diann_lib[order(Precursor.Id, LibraryIntensity, decreasing=TRUE),]
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

# ggplot(diann, aes(DIANN_PRED_FRAGMENTS)) +
#   geom_bar(stat = 'count')
#
# ggplot(diann, aes(DIANN_ASSAY_FRAGMENTS)) +
#   geom_bar(stat = 'count')

# ---- add fragment intensities to predictions ----
# NOTE: only max 12 quan values
fragment_quant_raw <- diann[, list(Fragment.Quant.Raw = unlist(strsplit(Fragment.Quant.Raw, ";"))), by=.(Precursor.Id, File.Name, Precursor.Quantity)]
fragment_quant_raw[,Fragment.Quant.Raw:=as.numeric(Fragment.Quant.Raw)]
fragment_quant_raw[,Fragment.Id:=seq(1, .N), by = .(Precursor.Id, File.Name, Precursor.Quantity)]

fragment_quant_corrected <- diann[, list(Fragment.Quant.Corrected = unlist(strsplit(Fragment.Quant.Corrected, ";"))), by=.(Precursor.Id, File.Name, Precursor.Quantity)]
fragment_quant_corrected[,Fragment.Quant.Corrected:=as.numeric(Fragment.Quant.Corrected)]
fragment_quant_corrected[,Fragment.Id:=seq(1, .N), by = .(Precursor.Id, File.Name, Precursor.Quantity)]

table(diann_lib$DIANN_PRED_FRAGMENTS)
table(diann_lib$DIANN_PRED_FRAGMENTS > 12)
diann_lib[DIANN_PRED_FRAGMENTS > 12, uniqueN(Precursor.Id)]/uniqueN(diann_lib$Precursor.Id) # 0.01841837

diann_lib_formatted <- merge(diann_lib,
                             fragment_quant_raw,
                             by = c('Precursor.Id',
                                    'Fragment.Id'),
                             all.x = T)

# --- "real" quan fragments per Precursor.Id and File.Name -----
fragment_quant_corrected <- merge(fragment_quant_corrected,
                                  unique(diann_lib[, list(Precursor.Id, DIANN_PRED_FRAGMENTS, DIANN_ASSAY_FRAGMENTS)]),
                                  by = 'Precursor.Id')
fragment_quant_corrected[,QUAN_FRAGS:=uniqueN(Fragment.Id[Fragment.Quant.Corrected > 0]), by = .(Precursor.Id, File.Name)]

isEqual = function(x, ref, rel_tol=1e-7, abs_tol=0){
  stopifnot(length(ref)==1 || length(x)==length(ref))
  isRelEqual = abs(1-x/ref)<=rel_tol
  isAbsEqual = abs(x-ref)<=abs_tol
  return(isRelEqual | isAbsEqual)
}

## one fragment for quan
fragment_quant_corrected_1 <- fragment_quant_corrected[QUAN_FRAGS >= 1 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_1 <- fragment_quant_corrected_1[,list(NotExcludedFromQuan = any(combn(Fragment.Quant.Corrected, 1, FUN = sum) == unique(Precursor.Quantity))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_1 <- fragment_quant_corrected_1[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_1[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_1[,NotExcludedFromQuan:=1]

## two fragment for quan
fragment_quant_corrected_2 <- fragment_quant_corrected[QUAN_FRAGS >= 2 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_2 <- fragment_quant_corrected_2[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 2, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_2 <- fragment_quant_corrected_2[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_2[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_2[,NotExcludedFromQuan:=2]

## three fragment for quan
fragment_quant_corrected_3 <- fragment_quant_corrected[QUAN_FRAGS >= 3 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_3 <- fragment_quant_corrected_3[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 3, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_3 <- fragment_quant_corrected_3[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_3[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_3[,NotExcludedFromQuan:=3]

## four fragment for quan
fragment_quant_corrected_4 <- fragment_quant_corrected[QUAN_FRAGS >= 4 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_4 <- fragment_quant_corrected_4[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 4, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_4 <- fragment_quant_corrected_4[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_4[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_4[,NotExcludedFromQuan:=4]

## five fragment for quan
fragment_quant_corrected_5 <- fragment_quant_corrected[QUAN_FRAGS >= 5 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_5 <- fragment_quant_corrected_5[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 5, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_5 <- fragment_quant_corrected_5[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_5[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_5[,NotExcludedFromQuan:=5]

## six fragment for quan
fragment_quant_corrected_6 <- fragment_quant_corrected[QUAN_FRAGS >= 6 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_6 <- fragment_quant_corrected_6[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 6, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_6 <- fragment_quant_corrected_6[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_6[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_6[,NotExcludedFromQuan:=6]

## seven fragment for quan
fragment_quant_corrected_7 <- fragment_quant_corrected[QUAN_FRAGS >= 7 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_7 <- fragment_quant_corrected_7[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 7, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_7 <- fragment_quant_corrected_7[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_7[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_7[,NotExcludedFromQuan:=7]

## eight fragments for quan
fragment_quant_corrected_8 <- fragment_quant_corrected[QUAN_FRAGS >= 8 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_8 <- fragment_quant_corrected_8[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 8, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_8 <- fragment_quant_corrected_8[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_8[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_8[,NotExcludedFromQuan:=8]

## nine fragments for quan
fragment_quant_corrected_9 <- fragment_quant_corrected[QUAN_FRAGS >= 9 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_9 <- fragment_quant_corrected_9[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 9, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_9 <- fragment_quant_corrected_9[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_9[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_9[,NotExcludedFromQuan:=9]

## ten fragments for quan
fragment_quant_corrected_10 <- fragment_quant_corrected[QUAN_FRAGS >= 10 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_10 <- fragment_quant_corrected_10[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 10, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_10 <- fragment_quant_corrected_10[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_10[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_10[,NotExcludedFromQuan:=10]

## eleven fragments for quan
fragment_quant_corrected_11 <- fragment_quant_corrected[QUAN_FRAGS >= 11 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_11 <- fragment_quant_corrected_11[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 11, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_11 <- fragment_quant_corrected_11[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_11[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_11[,NotExcludedFromQuan:=11]

## twelve fragments for quan
fragment_quant_corrected_12 <- fragment_quant_corrected[QUAN_FRAGS >= 12 &  Fragment.Quant.Corrected > 0]
fragment_quant_corrected_12 <- fragment_quant_corrected_12[,list(NotExcludedFromQuan = any(isEqual(combn(Fragment.Quant.Corrected, 12, FUN = sum), unique(Precursor.Quantity), 1e-5))), by=.(Precursor.Id, File.Name)]
fragment_quant_corrected_12 <- fragment_quant_corrected_12[NotExcludedFromQuan == TRUE]
fragment_quant_corrected_12[,NotExcludedFromQuan:=NULL]
fragment_quant_corrected_12[,NotExcludedFromQuan:=12]

fragment_for_quant <- rbind(fragment_quant_corrected_1,
                            fragment_quant_corrected_2,
                            fragment_quant_corrected_3,
                            fragment_quant_corrected_4,
                            fragment_quant_corrected_5,
                            fragment_quant_corrected_6)

## sanity check
fragment_for_quant[,pasted:=paste0(Precursor.Id, '_', File.Name)]
diann[,pasted:=paste0(Precursor.Id, '_', File.Name)]
length(diann[Precursor.Quantity != 0, pasted])
length(fragment_for_quant$pasted)
length(fragment_for_quant[NotExcludedFromQuan <= 6, pasted])

## min and max fragment combination per Precursor.Id and File.Name
fragment_for_quant_min <- setDT(fragment_for_quant)[, .(NotExcludedFromQuan = min(NotExcludedFromQuan)), by = .(Precursor.Id, File.Name)]
fragment_for_quant_max <- setDT(fragment_for_quant)[, .(NotExcludedFromQuan = max(NotExcludedFromQuan)), by = .(Precursor.Id, File.Name)]

table(fragment_for_quant_min$NotExcludedFromQuan)
table(fragment_for_quant_max$NotExcludedFromQuan)

## sanity check
fragment_for_quant_min[,pasted:=paste0(Precursor.Id, '_', File.Name)]
fragment_for_quant_max[,pasted:=paste0(Precursor.Id, '_', File.Name)]
length(diann[Precursor.Quantity != 0, pasted])
length(fragment_for_quant_min$pasted)
length(fragment_for_quant_max$pasted)
length(fragment_for_quant_max[NotExcludedFromQuan <= 6, pasted])

# ---- save data ----
write.fst(fragment_quant_raw, file.path(figurePath, 'fst-backup/20241127_diann_lfq_height_pepEntr_fragment_quant_raw.fst'), compress = 100)
write.fst(fragment_quant_corrected, file.path(figurePath, 'fst-backup/20241127_diann_lfq_height_pepEntr_fragment_quant_corrected.fst'), compress = 100)
write.fst(fragment_for_quant_min, file.path(figurePath, 'fst-backup/20241127_diann_lfq_height_pepEntr_minQuanFrags.fst'), compress = 100)
write.fst(fragment_for_quant_max, file.path(figurePath, 'fst-backup/20241127_diann_lfq_height_pepEntr_maxQuanFrags.fst'), compress = 100)

# ---- add to diann ----
diann[,pasted:=NULL]
fragment_for_quant_min[,pasted:=NULL]
setnames(fragment_for_quant_min, 'NotExcludedFromQuan', 'NotExcludedFromQuanMin')
diann <- merge(diann, fragment_for_quant_min, by = c('Precursor.Id', 'File.Name'), all.x = T)

fragment_for_quant_max[,pasted:=NULL]
setnames(fragment_for_quant_max, 'NotExcludedFromQuan', 'NotExcludedFromQuanMax')
diann <- merge(diann, fragment_for_quant_max, by = c('Precursor.Id', 'File.Name'), all.x = T)

# ---- save data ----
write.fst(diann, file.path(figurePath, 'fst-backup/20241127_diann_lfq_height_pepEntr_fragCounts.fst'), compress = 100)

ggplot(fragment_for_quant_min, aes(NotExcludedFromQuanMin)) +
  geom_histogram()
