# ---- setup ----
path <- file.path(here::here(), "figure-E3")
figurePath <- file.path(dataPath, "data/figure-E3")
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))

convertBlobAsDoublePd <- function (hexBlob) {
  if(all(is.na(hexBlob)) || length(hexBlob)==0) {
    return(NA)
  } else {
    rawbytes <- hexBlob[rep(c(rep(T, 8), F), length(hexBlob)/9)]
    return(readBin(rawbytes, "double", n = length(rawbytes)/8, size = 8, endian = "little"))
  }
}

# ---- pathsToData ----
## chimerys
pathToPdResult <- file.path(dataPath, "LFQ_Bench_human/Chimerys/LFQ_01_CHIMERYS_v2x7x9_apex_True.pdResult")
pathToMs8 <- file.path(dataPath, "LFQ_Bench_human/Chimerys/main_search-1.ms8")
pathToPercolatorInput <- file.path(dataPath, "LFQ_Bench_human/Chimerys/percolator-noSplit.tab")

## fastas
pathToFasta <- file.path(dataPath, "FASTA/CHIMERYS_Benchmark_human-canonical.fasta")
pathToContaminantsFasta <- file.path(dataPath, "FASTA/jpr_2022_contaminants.fasta")

# ---- read data including short sanity checks ----
proteinsFromFasta <- readProteinsFromFastas(pathToFasta)
proteinsFromContaminantsFasta <- readProteinsFromFastas(pathToContaminantsFasta)
proteinsFromFasta <- rbind(proteinsFromFasta,
                           proteinsFromContaminantsFasta)

pdresult <- readPdResult_psmGrouper_dda_quan_spectralAngles(pathToPdResult,
                                                            pathToMs8,
                                                            proteinsFromFasta,
                                                            pathToPercolatorInput)
any(is.na(pdresult$ORGANISM)) # FALSE
table(pdresult$ORGANISM) # NO ENTRAPMENT
uniqueN(pdresult[Q_VALUE <= 0.01, PCM_J_ID]) # 83598
uniqueN(pdresult[Q_VALUE <= 0.01 & DECOY == 0, PCM_J_ID]) # 82865

# ---- save ----
write.fst(pdresult, file.path(figurePath, 'intermediate/20241119_figureS4_pdresult.fst'), compress = 100)


# Figures A, B ====
data_venn_ptm_all_duell <- fread(file.path(dataPath, "data/figure-E2/figure-E2C-venn-peptides-all.csv"))

dtChim <- data_venn_ptm_all_duell[, .(isChimerys = any(condition %in% "CHIMERYS"),
                                      .N), keyby=ptm_group_J]
dtChim <- dtChim[isChimerys==T, .(isUnique = N==1), keyby=ptm_group_J]
dtChim[, .N, keyby=isUnique][, .(isUnique, N/sum(N))]
dtChim[, pep := gsub("^([A-Z]*)_.*$", "\\1", ptm_group_J)]
dtChim[, MODPEP_J_ID := paste0(pep, "_c", (nchar(pep)-nchar(gsub("C", "", pep))),
                               "_m", gsub("^[A-Z]*_.*_(\\d{1,2})x35_.*$", "\\1", ptm_group_J), "_p0")]
dtChim[, pep := NULL]
dtChim[, ptm_group_J := NULL]
setkey(dtChim, MODPEP_J_ID)


#peptide group data
vecCols <- c("ID", "QUAN", "FRAGMENTS_USED", "SPECTRUM_SIMILARITY", "Q_VALUE", "DECOY", "MODPEP_ID", "MODPEP_J_ID", "SVMSCORE")
dt <- read_fst(file.path(figurePath, 'intermediate/20241119_figureS4_pdresult.fst'), as.data.table = T, columns = vecCols)
# Highest spectral angle per modified peptide and sample.
# Highest number of matched peaks per modified peptide and sample.
# Possible alternative: highest counts per psm with highest MS8_feature_7, i.e. coefficient
# dt <- unique(dt[Q_VALUE<=0.01 & DECOY==0])
dt <- unique(dt)
#dt[, .N, by=ID][order(N)]
dt[ID==204362]

setkey(dt, MODPEP_J_ID)
setkey(dtChim, MODPEP_J_ID)
dt <- dtChim[dt]
# dt <- dt[!is.na(isUnique)] #962

dt[DECOY==1,Category:='Decoy']
dt[isUnique==T & DECOY==0,Category:='CHIMERYS unique target']
dt[isUnique==F & DECOY==0,Category:='Shared target']
dt[is.na(isUnique) & DECOY==0,Category:='FDR >1% target']
dt[,Category:=factor(Category, levels = c('Shared target', 'CHIMERYS unique target', 'FDR >1% target', 'Decoy'), ordered = T)]

fwrite(dt[Q_VALUE<=0.01 & !is.na(QUAN) & !is.na(FRAGMENTS_USED) & !is.na(SPECTRUM_SIMILARITY)],
       file.path(figurePath, "figure-E3A-intensity.csv"))

fwrite(dt[Q_VALUE<=0.01 & !is.na(SVMSCORE) & !is.na(FRAGMENTS_USED) & !is.na(SPECTRUM_SIMILARITY)],
       file.path(figurePath, "figure-E3B-fragments.csv"))

fwrite(dt[Q_VALUE<=0.01 & !is.na(SVMSCORE) & !is.na(FRAGMENTS_USED) & !is.na(SPECTRUM_SIMILARITY)],
       file.path(figurePath, "figure-E3C-spectral-similarity.csv"))

fwrite(dt[!is.na(SVMSCORE)], file.path(figurePath, "figure-E3D-svm.csv"))
