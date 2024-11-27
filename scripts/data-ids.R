require(data.table)
require(foreach)
require(RSQLite)


readIds <- function(dataPath, searchEngine = c("all", "pdResult", "MaxQuant"),
                    writeOutput = F, outputPath, outputName = "ID-counts.csv") {
  #check input arguments
  stopifnot((all(file_test("-f", dataPath)) || (file_test("-d", dataPath) && length(dataPath)==1L)))
  searchEngine <- match.arg(searchEngine)
  stopifnot(is.logical(writeOutput), length(writeOutput)==1L)
  if(missing(outputPath)) outputPath <- getwd()
  stopifnot(file_test("-d", outputPath), length(outputPath)==1L)
  stopifnot(is.character(outputName), length(outputName)==1L)
  
  #identify files to be read
  file_types <- c("all" = "^.*(\\.pdResult|summary\\.txt|psm\\.tsv)$",
                  "pdResult" = "^.*\\.pdResult$",
                  "MaxQuant" = "^.*summary\\.txt$",
                  "MSFragger" = "^.*psm\\.tsv")
  filePaths <- detectFiles(dataPath, file_types[searchEngine])
  progress_bar <- txtProgressBar(min = 0, max = length(filePaths), style = 3)
  
  #load each file via loop
  data_all <- foreach(filePath = filePaths, .combine = c) %do% {
    
    if(grepl(file_types["pdResult"], filePath)) {
      data_count <- readIdsPd(filePath)
    } else if(grepl(file_types["MaxQuant"], filePath)) {
      data_count <- readIdsMaxQuant(filePath)
    } else if(grepl(file_types["MSFragger"], filePath)) {
      data_count <- readIdsMsFragger(filePath)
    }
    
    setTxtProgressBar(progress_bar, which(filePaths %in% filePath))
    return(list(data_count))
  }
  close(progress_bar)
  data_all <- rbindlist(data_all, fill=T)
  
  #reorder columns
  id_columns <-
    c("folder_name", "file_name", "studies",
      "psms_all", "psms_FDR",
      "psms_quantified_all", "psms_quantified_FDR",
      "pcms_all", "pcms_FDR",
      "peptideGroups_all", "peptideGroups_FDR",
      "proteins_all", "proteins_FDR", "proteins_FDR_min2",
      "proteinGroups_all", "proteinGroups_FDR", "proteinGroups_FDR_min2")
  id_columns <- id_columns[id_columns %in% names(data_all)]
  setcolorder(data_all, id_columns)
  
  if(data_all[studies>1, .N]>0) {
    warning("more than 1 study per result file detected")
  }
  if(writeOutput==T) {
    fwrite(data_all, file = file.path(outputPath, outputName))
  }
  return(data_all)
}



readIdsPd <- function(filePath) {
  stopifnot(file_test("-f", filePath))
  
  #read number of studies
  db_connection <- dbConnect(RSQLite::SQLite(), filePath)
  query <- "SELECT COUNT(*) AS studies FROM StudyInformation"
  count_studies <- as.data.table(dbGetQuery(db_connection, query))
  
  #check for quantification columns
  psm_names <- dbListFields(db_connection, "TargetPsms")
  quan_names <- c("QuanValue", "PrecursorAbundance", "Abundances")
  quan_names <- quan_names[quan_names %in% psm_names]
  if(length(quan_names)>1) {
    warning(paste("found", length(quan_names), "quantification columns for",
                  basename(filePath), "using", quan_names[1]))
  }
  
  #read PSMs for all
  if(length(quan_names)==0) {
    query <- "SELECT ModifiedSequence, Charge, FirstScan, StudyFileId, qValue FROM TargetPsms"
    dt <- as.data.table(dbGetQuery(db_connection, query))
    dt[, psm := paste(gsub("(I|L)", "J", ModifiedSequence), Charge, FirstScan, StudyFileId)]
    count_psm <- dt[, .(psms_all = uniqueN(psm),
                        psms_FDR = .SD[qValue<=0.01, uniqueN(psm)])]
  } else {
    query <- sprintf("SELECT ModifiedSequence, Charge, FirstScan, StudyFileId, qValue, %s FROM TargetPsms", quan_names[1])
    dt <- as.data.table(dbGetQuery(db_connection, query))
    dt[, psm := paste(gsub("(I|L)", "J", ModifiedSequence), Charge, FirstScan, StudyFileId)]
    count_psm <- dt[, .(psms_all = uniqueN(psm),
                        psms_FDR = .SD[qValue<=0.01, uniqueN(psm)],
                        psms_quantified_all = .SD[get(quan_names[1])>1, uniqueN(psm)],
                        psms_quantified_FDR = .SD[qValue<=0.01 & get(quan_names[1])>1, uniqueN(psm)])]
  }
  
  #return results as list for wrapping via rbindlist
  if(grepl("^.*\\.pdResult$", filePath)) {
    #read peptides and proteins if pdResult (not msf)
    query <- "SELECT qValue FROM TargetPeptideGroups"
    dt <- as.data.table(dbGetQuery(db_connection, query))
    count_pep <- dt[, .(peptideGroups_all = .N,
                        peptideGroups_FDR = .SD[qValue<=0.01, .N])]
    
    query <- "SELECT Accession, ProteinGroupIDs, Expqvalue, IsMasterProtein, UniquePeptidesCount
    FROM TargetProteins WHERE Expqvalue IS NOT NULL"
    dt <- as.data.table(dbGetQuery(db_connection, query))
    dt[, qValue := sapply(Expqvalue, convertBlobAsDoublePd)]
    
    #extract and map UniquePeptidesCount from proteinGroups (NOT same as from proteins)
    dt_unique <- dt[, .(ProteinGroupID = as.integer(unlist(strsplit(ProteinGroupIDs, ";")))), by = Accession]
    setkey(dt_unique, ProteinGroupID)
    query <- "SELECT * FROM TargetProteinGroups WHERE ExcludedBy IS -1"
    dt_prgrp <- as.data.table(dbGetQuery(db_connection, query))
    setnames(dt_prgrp, "UniquePeptidesCount", "grp_UniquePeptidesCount")
    setkey(dt_prgrp, ProteinGroupID)
    
    dt_merge <- merge(dt_unique, dt_prgrp[, .(ProteinGroupID, grp_UniquePeptidesCount)])
    dt_merge <- dt_merge[, .(grp_UniquePeptidesCount = max(grp_UniquePeptidesCount)), keyby = Accession]
    setkey(dt, Accession)
    dt <- merge(dt, dt_merge, all.x=T, all.y=F)
    
    count_prot <- dt[, .(proteins_all = .N,
                         proteins_FDR = .SD[qValue<=0.01, .N],
                         proteins_FDR_min2 = .SD[qValue<=0.01 & UniquePeptidesCount>1, .N],
                         proteinGroups_all = .SD[IsMasterProtein==0, .N],
                         proteinGroups_FDR = .SD[qValue<=0.01 & IsMasterProtein==0, .N],
                         proteinGroups_FDR_min2 = .SD[qValue<=0.01 & IsMasterProtein==0 & grp_UniquePeptidesCount>1, .N])]
    
    dbDisconnect(db_connection)
    return(cbind(folder_name = basename(dirname(filePath)),
                 file_name = basename(filePath),
                 count_studies, count_psm, count_pep, count_prot))
    
  } else {
    dbDisconnect(db_connection)
    return(cbind(folder_name = basename(dirname(filePath)),
                 file_name = basename(filePath),
                 count_studies, count_psm))
  }
}



readIdsMaxQuant <- function(filePath) {
  stopifnot(file_test("-f", filePath))
  
  data_sum <- fread(filePath)
  count_studies <- data_sum[!`Raw file` %in% "Total", .(studies = .N)]
  
  data_msms <- fread(file.path(dirname(filePath), "msms.txt"),
                     select = c("Reverse", "Precursor Intensity"))
  count_psm <- data_msms[Reverse!="+", .(psms_FDR = .N,
                                         psms_quantified_FDR = .SD[!is.na(`Precursor Intensity`), .N])]
  
  data_ev <- fread(file.path(dirname(filePath), "evidence.txt"),
                   select = c("Reverse", "Modified sequence", "Charge"))
  count_pcm <- data_ev[Reverse!="+", .(pcms_FDR = uniqueN(paste(`Modified sequence`, Charge)))]
  
  data_pep <- fread(file.path(dirname(filePath), "modificationSpecificPeptides.txt"),
                    select = "Reverse")
  count_pep <- data_pep[Reverse!="+", .(peptideGroups_FDR = .N)]
  
  # data_pep <- fread(file.path(dirname(filePath), "peptides.txt"),
  #                   select = "Reverse")
  # count_pep <- data_pep[Reverse!="+", .(peptideGroups_FDR = .N)]
  
  data_prgrp <- fread(file.path(dirname(filePath), "proteinGroups.txt"),
                      select = c("Reverse", "Unique peptides"))
  count_prot <- data_prgrp[Reverse!="+", .(proteinGroups_FDR = .N,
                                           proteinGroups_FDR_min2 = .SD[`Unique peptides`>1, .N])]
  #`Only identified by site`!="+"
  
  return(cbind(folder_name = basename(dirname(filePath)),
               file_name = basename(filePath),
               count_studies, count_psm, count_pcm, count_pep, count_prot))
}



readIdsMsFragger <- function(filePath) {
  stopifnot(file_test("-f", filePath))
  
  data_msms <- fread(filePath, select = c("Spectrum File", "Intensity"))
  count_studies <- data_msms[, .(studies = uniqueN(`Spectrum File`))]
  
  data_msms[Intensity==0, Intensity := NA]
  count_psm <- data_msms[, .(psms_FDR = .N,
                             psms_quantified_FDR = .SD[!is.na(`Intensity`), .N])]
  
  data_ev <- fread(file.path(dirname(filePath), "ion.tsv"), select = c("Intensity"))
  data_ev[Intensity==0, Intensity := NA]
  count_pcm <- data_ev[, .(pcms_FDR = .N)]
  
  data_pep <- fread(file.path(dirname(filePath), "peptide.tsv"), select = c("Intensity"))
  data_pep[Intensity==0, Intensity := NA]
  count_pep <- data_pep[, .(peptides_FDR = .N)]
  
  data_prgrp <- fread(file.path(dirname(filePath), "protein.tsv"), select = c("Unique Peptides"))
  count_prot <- data_prgrp[, .(proteinGroups_FDR = .N,
                               proteinGroups_FDR_min2 = .SD[`Unique Peptides`>1, .N])]
  
  return(cbind(folder_name = basename(dirname(filePath)),
               file_name = basename(filePath),
               count_studies, count_psm, count_pcm, count_pep, count_prot))
}



#BLOB conversion function adapted for Proteome Discoverer (include unwanted separator that needs to be removed)
convertBlobAsDoublePd <- function (hexBlob) {
  if(is.na(hexBlob) || length(hexBlob)==0) {
    return(NA)
  } else {
    rawbytes <- hexBlob[rep(c(rep(T, 8), F), length(hexBlob)/9)]
    return(readBin(rawbytes, "double", n = length(rawbytes)/8, size = 8, endian = "little"))
  }
}


detectFiles <- function(dataPath, dataPattern, recursive = T) {
  #test if dataPath is folder; if yes search for files, if no return unchanged
  if((file_test("-d", dataPath) && length(dataPath)==1L)) {
    files <- list.files(dataPath, pattern = dataPattern, recursive = recursive, full.names = T)
    files <- files[!grepl("_old", files)]
    if(length(files)==0) {
      stop(paste("[error] no", dataPattern, "files detected in folder path"))
    }
    return(files)
  } else {return(dataPath)}
}
