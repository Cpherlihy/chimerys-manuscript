#setup
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

sampleNamesPath <- writeSampleNames(file.path(dataPath, "LFQ_Bench_human"),
                                    outputPath = path, outputName = "sample_names.csv")
filePath_CHIMERYS <- file.path(dataPath, "LFQ_Bench_human/Chimerys/LFQ_01_CHIMERYS_v2x7x9_apex_True.pdResult")
filePath_Sequest <- file.path(dataPath, "LFQ_Bench_human/Sequest-HT/LFQ_01_Sequest.pdResult")
prot <- rbind(readPdProteins(filePath_CHIMERYS, sampleNamesPath, level = "protGroups", loadBackup = F),
              readPdProteins(filePath_Sequest, sampleNamesPath, level = "protGroups", loadBackup = F))

prot_peptides <- dcast(prot, protein~condition, value.var = "peptide_count")
prot_peptides_count <- prot_peptides[, .N, by = c("CHIMERYS", "SequestHT")]
setnames(prot_peptides_count, c('CHIMERYS', 'SequestHT'), c('CHIMERYS', 'Sequest-HT'))

fwrite(prot_peptides_count[`Sequest-HT`>0 & CHIMERYS>0],
       file.path(figurePath, "figure-E3E-proteins.csv"))
