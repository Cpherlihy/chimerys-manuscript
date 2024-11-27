writeSampleNames <- function(dataPath, searchEngine = "all",
                             outputPath = getwd(), outputName = "sample_names.csv", overwrite = F) {
  dataPathIsPath <- all(file_test("-d", dataPath)) && length(dataPath)==1L
  stopifnot((all(file_test("-f", dataPath)) || dataPathIsPath==T))
  searchEngine <- match.arg(searchEngine, choices = names(listFileTypes(level = "study")))
  stopifnot(file_test("-d", outputPath), length(outputPath)==1L)
  stopifnot(is.character(outputName), length(outputName)==1L)
  stopifnot(is.logical(overwrite), length(overwrite)==1L)
  
  if(!overwrite && file.exists(file.path(outputPath, outputName))) {
    message(file.path(outputPath, outputName), " exists already, skipping file export")
    return(file.path(outputPath, outputName))
  }
  
  #identify files to be read
  file_types <- listFileTypes(level = "study")
  filePaths <- detectFiles(dataPath, file_types[searchEngine])
  dataPath <- createSharedFilePath(dataPath)
  
  #initiate progress bar
  n_files <- length(filePaths)
  progress_bar <- txtProgressBar(min = 0, max = n_files, width = n_files, style = 3)
  
  data_samples <- foreach(filePath = filePaths, i = 1:n_files, .errorhandling = "pass") %do% {
    setTxtProgressBar(progress_bar, i-1)
    studyFolderName <- createStudyFolder(filePath, dataPath, dataPathIsPath)
    message(paste(" Loading", file.path(studyFolderName, basename(filePath))))
    
    #check which file type and load samples accordingly
    if(grepl(file_types["pdResult"], basename(filePath))) {
      db_connection <- dbConnect(RSQLite::SQLite(), filePath)
      query_samples <- paste('SELECT SampleIdentifier FROM StudyInformation')
      temp <- as.data.table(dbGetQuery(db_connection, query_samples))
      temp[, raw_file := factor(SampleIdentifier, levels = unique(SampleIdentifier))]
      dbDisconnect(db_connection)
      temp <- temp[, .(order = 1L,
                       study_folder = studyFolderName,
                       study_file = basename(filePath),
                       raw_file)]
      
    } else if(grepl(file_types["DIANN"], basename(filePath))) {
      temp <- fread(filePath, select = "Run")
      temp <- temp[, .(order = 1L,
                       study_folder = studyFolderName,
                       study_file = basename(filePath),
                       raw_file = sort(unique(Run)))]
      
    } else if(grepl(file_types["MaxQuant"], basename(filePath))) {
      temp <- fread(filePath)
      temp <- temp[-.N, .(order = 1L,
                          study_folder = studyFolderName,
                          study_file = basename(filePath),
                          raw_file = sort(`Raw file`))]
      
    } else if(grepl(file_types["FragPipe"], basename(filePath))) {
      temp <- fread(filePath, header = F)
      setnames(temp, c("raw", "experiment", "replicate", "MS"))
      
      #check if experiments and/or replicates are defined
      #if so, they will be stored in subfolders, which need to be appended to the study_folders
      if(temp[!is.na(experiment) | !is.na(replicate), .N]>0) {
        study_folders <- temp[, apply(.SD, 1, function(x) paste(x[!is.na(x)], collapse = "_") ),
                              .SDcols = c("experiment", "replicate")]
        study_folders <- file.path(studyFolderName, study_folders)
      } else {
        study_folders <- studyFolderName
      }
      temp <- temp[, .(order = 1L,
                       study_folder = study_folders,
                       study_file = basename(filePath),
                       raw_file = sort(gsub("^([^.]*?){0,1}\\..*$", "\\1", basename(gsub("\\", "/", raw, fixed=T)))))]
      
    } else if(grepl(file_types["Sage"], basename(filePath))) {
      temp <- fread(filePath, select = "filename")
      temp <- temp[, .(order = 1L,
                       study_folder = studyFolderName,
                       study_file = basename(filePath),
                       raw_file = tools::file_path_sans_ext(sort(unique(filename))))]
      
    } else if(grepl(file_types["Percolator"], basename(filePath))) {
      #Percolator output is file-by-file, not storing raw file information
      temp <- fread(filePath, select = "PSMId")
      temp <- temp[, .(order = 1L,
                       study_folder = studyFolderName,
                       study_file = basename(filePath),
                       raw_file = sort(unique(gsub("^(.*)\\.\\d*\\.\\d*\\.\\d_\\d$", "\\1", PSMId))))]
      
    } else if(grepl(file_types["Metamorpheus"], basename(filePath))) {
      temp <- fread(filePath, select = "File Name")
      temp <- temp[, .(order = 1L,
                       study_folder = studyFolderName,
                       study_file = basename(filePath),
                       raw_file = gsub("-calib", "", sort(unique(`File Name`)), fixed = T))]
      
    } else if(grepl(file_types["MSGFplus"], basename(filePath))) {
      temp <- fread(filePath, select = "#SpecFile")
      temp <- temp[, .(order = 1L,
                       study_folder = studyFolderName,
                       study_file = basename(filePath),
                       raw_file = tools::file_path_sans_ext(sort(unique(`#SpecFile`))))]
    } else {
      stop(paste("no sample parsing defined for file", basename(filePath)))
    }
    
    return(temp)
  }
  setTxtProgressBar(progress_bar, n_files)
  message(" Merging and processing")
  
  #if errors occurred, report files that caused them and return as summarized warnings
  is_error <- sapply(data_samples, inherits, "error")
  if(any(is_error)==T) {
    file_error <- paste0((1:n_files)[is_error], " (", basename(filePaths)[is_error], ")", collapse = ", ")
    warning(paste("file(s) skipped due to error(s):", file_error), call. = F)
    warnErrList(data_samples[is_error])
  }
  
  #return non-error data
  data_samples <- rbindlist(data_samples[!is_error])
  setkey(data_samples, study_folder, study_file, raw_file)
  data_samples[, order := 1L:.N]
  
  if(overwrite==T) {
    data_samples_old <- fread(file.path(outputPath, outputName), na.strings = "")
    by_cols <- c("study_folder", "study_file", "raw_file")
    data_samples <- merge(data_samples,
                          data_samples_old[, .SD, .SDcols = c(by_cols, "sample_name")],
                          by = , all.x=T, all.y=F)
  } else {
    data_samples[, sample_name := character()]
    data_samples[1, sample_name := "condition-treatment_1"]
    if(data_samples[, .N] >= 2) data_samples[2, sample_name := "condition-treatment_2"]
    if(data_samples[, .N] >= 3) data_samples[3, sample_name := "condition-control_1"]
    if(data_samples[, .N] >= 4) data_samples[4, sample_name := "condition-control_2"]
  }
  
  setcolorder(data_samples, c("order", "study_folder", "study_file", "raw_file", "sample_name"))
  setkey(data_samples, order)
  fwrite(data_samples, file.path(outputPath, outputName))
  message("DONE - please edit column 'sample_name' in ", outputName)
  return(file.path(outputPath, outputName))
}


readSampleNames <- function(filePath) {
  stopifnot(file_test("-f", filePath), length(filePath)==1L)
  
  data_samples <- fread(filePath, na.strings = "")
  if(!all(c("order", "study_folder", "study_file", "raw_file", "sample_name") %in% names(data_samples))) {
    stop(paste("'order', 'study_folder", "study_file', 'raw_file' and 'sample_name' need to be part of", basename(filePath)))
  }
  if(!data_samples[, all(grepl("^.*_\\d{1,2}$", sample_name))]) {
    stop("'sample_name' column needs to be of format condition-treatment_1, with replicate as integer after underscore")
  }
  data_samples[, sample_name := gsub("\\n", "\n", sample_name, fixed = T)]
  setkey(data_samples, order, sample_name)
  
  data_samples[, sample_name_old := paste(study_folder, study_file, raw_file, sep = "_")]
  data_samples[, condition_name := gsub("^(.*)_\\d{1,2}$", "\\1", sample_name)]
  
  return(data_samples[])
}


readPsms <- function(dataPath, sampleNamesPath = NULL, searchEngine = "all",
                     parsingRulesPath = here::here("scripts/column_names.csv"), rbindlistFill = T, ...) {
  #check input arguments
  dataPathIsPath <- all(file_test("-d", dataPath)) && length(dataPath)==1L
  stopifnot((all(file_test("-f", dataPath)) || dataPathIsPath==T))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  searchEngine <- match.arg(searchEngine, choices = names(listFileTypes(level = "psms")))
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  stopifnot(is.logical(rbindlistFill), length(rbindlistFill)==1L)
  
  #identify files to be read
  file_types <- listFileTypes(level = "psms")
  filePaths <- detectFiles(dataPath, file_types[searchEngine])
  dataPath <- createSharedFilePath(dataPath)
  
  #initiate progress bar
  n_files <- length(filePaths)
  progress_bar <- txtProgressBar(min = 0, max = n_files, width = n_files, style = 3)
  
  #load each file via loop
  data_all <- foreach(filePath = filePaths, i = 1:n_files, .errorhandling = "pass") %do% {
    #start message line
    if(i>1) {message("")}
    setTxtProgressBar(progress_bar, i-1)
    studyFolderName <- createStudyFolder(filePath, dataPath, dataPathIsPath)
    message(" Loading ", file.path(studyFolderName, basename(filePath)), "...", appendLF = F)
    
    #read and return data
    data_all <- readSearchEngine(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, "psms", ...)
    message(" done", appendLF = F)
    return(data_all)
  }
  message("")
  setTxtProgressBar(progress_bar, n_files)
  message(" Merging and processing")
  
  #if errors occurred, report files that caused them and return as summarized warnings
  is_error <- sapply(data_all, inherits, "error")
  if(any(is_error)==T) {
    file_error <- paste0((1:n_files)[is_error], " (", basename(filePaths)[is_error], ")", collapse = ", ")
    warning(paste("file(s) skipped due to error(s):", file_error), call. = F)
    warnErrList(data_all[is_error])
  }
  #return non-error data
  data_all <- rbindlist(data_all[!is_error], fill=rbindlistFill)
  
  parseReorder(data_all, parsingRulesPath)
  close(progress_bar)
  return(data_all)
}


readPcms <- function(dataPath, sampleNamesPath = NULL, searchEngine = "all",
                     parsingRulesPath = here::here("scripts/column_names.csv"), rbindlistFill = T, ...) {
  #check input arguments
  dataPathIsPath <- all(file_test("-d", dataPath)) && length(dataPath)==1L
  stopifnot((all(file_test("-f", dataPath)) || dataPathIsPath==T))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  searchEngine <- match.arg(searchEngine, choices = names(listFileTypes(level = "pcms")))
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  stopifnot(is.logical(rbindlistFill), length(rbindlistFill)==1L)
  
  #identify files to be read
  file_types <- listFileTypes(level = "pcms")
  filePaths <- detectFiles(dataPath, file_types[searchEngine])
  dataPath <- createSharedFilePath(dataPath)
  
  #initiate progress bar
  n_files <- length(filePaths)
  progress_bar <- txtProgressBar(min = 0, max = n_files, width = n_files, style = 3)
  
  #load each file via loop
  data_all <- foreach(filePath = filePaths, i = 1:n_files, .errorhandling = "pass") %do% {
    #start message line
    if(i>1) {message("")}
    setTxtProgressBar(progress_bar, i-1)
    studyFolderName <- createStudyFolder(filePath, dataPath, dataPathIsPath)
    message(" Loading ", file.path(studyFolderName, basename(filePath)), "...", appendLF = F)
    
    #read and return data
    data_all <- readSearchEngine(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, "pcms", ...)
    message(" done", appendLF = F)
    return(data_all)
  }
  message("")
  setTxtProgressBar(progress_bar, n_files)
  message(" Merging and processing")
  
  #if errors occurred, report files that caused them and return as summarized warnings
  is_error <- sapply(data_all, inherits, "error")
  if(any(is_error)==T) {
    file_error <- paste0((1:n_files)[is_error], " (", basename(filePaths)[is_error], ")", collapse = ", ")
    warning(paste("file(s) skipped due to error(s):", file_error), call. = F)
    warnErrList(data_all[is_error])
  }
  #return non-error data
  data_all <- rbindlist(data_all[!is_error], fill=rbindlistFill)
  
  parseReorder(data_all, parsingRulesPath)
  close(progress_bar)
  return(data_all)
}



readPtmGroups <- function(dataPath, sampleNamesPath = NULL, searchEngine = "all",
                          parsingRulesPath = here::here("scripts/column_names.csv"), rbindlistFill = T, ...) {
  #check input arguments
  dataPathIsPath <- all(file_test("-d", dataPath)) && length(dataPath)==1L
  stopifnot((all(file_test("-f", dataPath)) || dataPathIsPath==T))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  searchEngine <- match.arg(searchEngine, choices = names(listFileTypes(level = "ptmGroups")))
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  stopifnot(is.logical(rbindlistFill), length(rbindlistFill)==1L)
  
  #identify files to be read
  file_types <- listFileTypes(level = "ptmGroups")
  filePaths <- detectFiles(dataPath, file_types[searchEngine])
  dataPath <- createSharedFilePath(dataPath)
  
  #initiate progress bar
  n_files <- length(filePaths)
  progress_bar <- txtProgressBar(min = 0, max = n_files, width = n_files, style = 3)
  
  #load each file via loop
  data_all <- foreach(filePath = filePaths, i = 1:n_files, .errorhandling = "pass") %do% {
    #start message line
    if(i>1) {message("")}
    setTxtProgressBar(progress_bar, i-1)
    studyFolderName <- createStudyFolder(filePath, dataPath, dataPathIsPath)
    message(" Loading ", file.path(studyFolderName, basename(filePath)), "...", appendLF = F)
    
    #read and return data
    data_all <- readSearchEngine(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, "ptmGroups", ...)
    message(" done", appendLF = F)
    return(data_all)
  }
  message("")
  setTxtProgressBar(progress_bar, n_files)
  message(" Merging and processing")
  
  #if errors occurred, report files that caused them and return as summarized warnings
  is_error <- sapply(data_all, inherits, "error")
  if(any(is_error)==T) {
    file_error <- paste0((1:n_files)[is_error], " (", basename(filePaths)[is_error], ")", collapse = ", ")
    warning(paste("file(s) skipped due to error(s):", file_error), call. = F)
    warnErrList(data_all[is_error])
  }
  #return non-error data
  data_all <- rbindlist(data_all[!is_error], fill=rbindlistFill)
  
  parseReorder(data_all, parsingRulesPath)
  close(progress_bar)
  return(data_all)
}



readPdPsms <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                       studyFolderName = NULL, searchEngineName = "PD_psms",
                       filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readPdPsms.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load pdResult data via SQL query
  db_connection <- dbConnect(RSQLite::SQLite(), filePath)
  query_targets <- "SELECT * FROM TargetPsms"
  if(filterFdr==T) {
    query_targets <- paste(query_targets, 'WHERE qValue <= 0.01')
  }
  data_temp <- as.data.table(dbGetQuery(db_connection, query_targets))
  data_temp[, is_decoy := F]
  #load decoys
  if(loadDecoys==T) {
    query_decoys <- "SELECT * FROM DecoyPsms"
    if(filterFdr==T) {
      query_decoys <- paste(query_decoys, 'WHERE qValue <= 0.01')
    }
    data_decoys <- as.data.table(dbGetQuery(db_connection, query_decoys))
    data_decoys[, is_decoy := T]
    data_temp <- rbind(data_temp, data_decoys, fill=T)
  }
  setnames(data_temp, paste0("old_", names(data_temp)))
  
  #load sample naming info via SQL query and create raw file names
  query_samples <- paste('SELECT ID, StudyFileID, FileName, SampleIdentifier FROM StudyInformation')
  data_study <- as.data.table(dbGetQuery(db_connection, query_samples))
  #parse raw file and channels (if applicable)
  data_study[, FileName := tools::file_path_sans_ext(basename(gsub("\\", "/", FileName, fixed=T)))]
  setkey(data_study, ID)
  data_channel <- data_study[, .N, by=.(StudyFileID, FileName)]
  data_temp[, raw_file := factor(old_StudyFileId, levels = data_channel$StudyFileID, labels = data_channel$FileName)]
  dbDisconnect(db_connection)
  
  #check if multi-channel intensities (e.g. TMT/SILAC)
  if("old_Abundances" %in% names(data_temp) && data_channel[, any(N>1)]) {
    #split data and study by StudyFileID
    data_temp[, raw_file := NULL]
    data_list <- split(data_temp, by = "old_StudyFileId")
    study_list <- split(data_study, by = "StudyFileID")
    #parse and melt channel names into raw file
    data_temp <- foreach(data_by = data_list, study_by = study_list, .combine = rbind) %do% {
      temp <- data_by[, (study_by$SampleIdentifier) := transpose(
        lapply(old_Abundances, convertBlobAsDoublePd))]
      melt(temp, measure.vars = study_by$SampleIdentifier,
           variable.name = "raw_file", value.name = "intensity")
    }
  }
  
  #define sample names and reformat columns
  parseSamples(data_temp, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  return(data_temp)
}



readPdPtmGroups <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                            studyFolderName = NULL, searchEngineName = "PD_ptmGroups",
                            filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readPdPtmGroups.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load pdResult data via SQL query
  db_connection <- dbConnect(RSQLite::SQLite(), filePath)
  query_targets <- "SELECT * FROM TargetPeptideGroups"
  if(filterFdr==T) {
    query_targets <- paste(query_targets, 'WHERE qValue <= 0.01')
  }
  data_temp <- as.data.table(dbGetQuery(db_connection, query_targets))
  data_temp[, is_decoy := F]
  #load decoys
  if(loadDecoys==T) {
    query_decoys <- "SELECT * FROM DecoyPeptideGroups"
    if(filterFdr==T) {
      query_decoys <- paste(query_decoys, 'WHERE qValue <= 0.01')
    }
    data_decoys <- as.data.table(dbGetQuery(db_connection, query_decoys))
    data_decoys[, is_decoy := T]
    data_temp <- rbind(data_temp, data_decoys, fill=T)
  }
  setnames(data_temp, paste0("old_", names(data_temp)))
  
  #load sample naming info via SQL query and create raw file names
  query_samples <- paste('SELECT ID, StudyFileID, SampleIdentifier FROM StudyInformation')
  data_study <- as.data.table(dbGetQuery(db_connection, query_samples))
  setkey(data_study, ID)
  dbDisconnect(db_connection)
  
  #explanations of PD columns for future reference
  #lapply(FoundinSamples, convertBlobAsIntegerPd) #3=high, 2=medium, 1=low, 0=NA/peak found (peak = MBR)
  #QuanInfo 1=no quan values, 2=shared (not unique)
  #AbundancesScaled is AbundancesNormalized with total adding up to 100 * n_samples
  #qValuebySearchEngine is min(qValue) from TargetPsms
  
  #melt data to tidy format based on FoundinSamples column
  #CAVE: old_PeptideGroupID is NOT unique across T and D
  data_temp[!is.na(old_FoundinSamples), data_study$SampleIdentifier :=
              transpose(lapply(old_FoundinSamples, convertBlobAsIntegerPd))]
  data_temp[, old_FoundinSamples := NULL]
  data_sample <- melt(data_temp, measure.vars = data_study$SampleIdentifier,
                      variable.name = "raw_file",
                      value.name = "old_FoundinSamples")
  data_sample[, raw_file := factor(raw_file, levels = data_study$SampleIdentifier)]
  setkey(data_sample, raw_file, old_is_decoy, old_PeptideGroupID)
  
  #if available, attach intensity column(s)
  if("old_Abundances" %in% names(data_temp)) {
    data_temp[, data_study$SampleIdentifier :=
                transpose(lapply(old_Abundances, convertBlobAsDoublePd))]
    data_abundance <- melt(data_temp, measure.vars = data_study$SampleIdentifier,
                           id.vars = c("old_is_decoy", "old_PeptideGroupID"),
                           variable.name = "raw_file",
                           value.name = "intensity")
    setkey(data_abundance, raw_file, old_is_decoy, old_PeptideGroupID)
    data_sample <- data_sample[data_abundance]
    data_sample[old_FoundinSamples==0 & intensity>0, old_FoundinSamples := 4]
  }
  
  if("old_AbundancesNormalized" %in% names(data_temp)) {
    data_temp[, data_study$SampleIdentifier :=
                transpose(lapply(old_AbundancesNormalized, convertBlobAsDoublePd))]
    data_abundance <- melt(data_temp, measure.vars = data_study$SampleIdentifier,
                           id.vars = c("old_is_decoy", "old_PeptideGroupID"),
                           variable.name = "raw_file",
                           value.name = "intensity_normalized")
    setkey(data_abundance, raw_file, old_is_decoy, old_PeptideGroupID)
    data_sample <- data_sample[data_abundance]
    data_sample[old_FoundinSamples==0 & intensity>0, old_FoundinSamples := 4]
  }
  
  #filter samples without PSMs and label PSM qvalues
  if(filterFdr==T) {
    data_sample <- data_sample[old_FoundinSamples>0]
  }
  data_sample[, qvalue_psm_label := factor(old_FoundinSamples, levels = c(3:1, 4),
                                           labels = c("<=0.01", "<=0.05", ">0.05", "MBR"))]
  
  #create ptm group column
  modifications <- data_sample[old_Modifications!="", .N, by = old_Modifications
  ][, unlist(strsplit(old_Modifications, ", "))]
  modifications <- unique(gsub("^.*×(.*) \\[.*\\]$", "\\1", modifications))
  parseModificationsColumn(dataTable = data_sample,
                           modifications = modifications,
                           modificationsColumn = "old_Modifications",
                           sequenceColumn = "old_Sequence",
                           unimodColumn = "PD_ptmGroups",
                           unimodExtraction = "gsub",
                           unimodLocation = "^(\\d{1,2})×.*$",
                           unimodSeparator = "; ")
  
  #define sample names and reformat columns
  parseSamples(data_sample, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_sample, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_sample, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_sample, fst_backup_path)
  }
  return(data_sample)
}



readPdProteins <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                           level = c("proteins", "protGroups"), studyFolderName = NULL, searchEngineName = "PD_proteinGroups",
                           filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  level <- match.arg(level)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readPdProteins.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load pdResult data via SQL query
  db_connection <- dbConnect(RSQLite::SQLite(), filePath)
  query_targets <- "SELECT * FROM TargetProteins WHERE PsmCount>0"
  data_temp <- as.data.table(dbGetQuery(db_connection, query_targets))
  data_temp[, Sequence := NULL]
  data_temp[, is_decoy := F]
  #load decoys
  if(loadDecoys==T) {
    query_decoys <- "SELECT * FROM DecoyProteins WHERE PsmCount>0"
    data_decoys <- as.data.table(dbGetQuery(db_connection, query_decoys))
    data_decoys[, is_decoy := T]
    data_temp <- rbind(data_temp, data_decoys, fill=T)
  }
  #filter FDR (requires BLOB decoding)
  data_temp[, qvalue_protein := sapply(Expqvalue, convertBlobAsDoublePd)]
  data_temp[, qvalue_protein_label := sapply(ProteinFDRConfidence, convertBlobAsIntegerPd)]
  data_temp[, qvalue_protein_label := factor(qvalue_protein_label, levels = 3:1,
                                             labels = c("<=0.01", "<=0.05", ">0.05"))]
  if(filterFdr==T) {
    #use label column since old PD versions had qvalue rounding error
    data_temp <- data_temp[qvalue_protein_label=="<=0.01"]
    
  }
  setnames(data_temp, paste0("old_", names(data_temp)))
  
  #load sample naming info via SQL query and create raw file names
  query_samples <- paste('SELECT ID, StudyFileID, SampleIdentifier FROM StudyInformation')
  data_study <- as.data.table(dbGetQuery(db_connection, query_samples))
  setkey(data_study, ID)
  dbDisconnect(db_connection)
  
  #explanations of PD columns for future reference
  #lapply(FoundinSamples, convertBlobAsIntegerPd) #3=high, 2=medium, 1=low, 0=NA/peak found
  #AbundancesScaled is AbundancesNormalized with total adding up to 100 * n_samples
  #UniquePeptidesCount: number of FDR-filtered TargetPeptideGroups that are unique to this protein
  #GroupUniquePeptidesCount: number of FDR-filtered TargetPeptideGroups that are unique to this protein group
  #PeptideGroupCount: number of FDR-filtered TargetPeptideGroups that are shared with this protein group (?)
  
  #melt data to tidy format based on FoundinSamples column
  data_temp[!is.na(old_FoundinSamples), data_study$SampleIdentifier :=
              transpose(lapply(old_FoundinSamples, convertBlobAsIntegerPd))]
  data_temp[, old_FoundinSamples := NULL]
  data_sample <- melt(data_temp, measure.vars = data_study$SampleIdentifier,
                      variable.name = "raw_file",
                      value.name = "old_FoundinSamples")
  data_sample[, raw_file := factor(raw_file, levels = data_study$SampleIdentifier)]
  setkey(data_sample, raw_file, old_UniqueSequenceID)
  
  #if available, attach intensity column(s)
  if("old_Abundances" %in% names(data_temp)) {
    data_temp[, data_study$SampleIdentifier :=
                transpose(lapply(old_Abundances, convertBlobAsDoublePd))]
    data_abundance <- melt(data_temp, measure.vars = data_study$SampleIdentifier,
                           id.vars = c("old_UniqueSequenceID"),
                           variable.name = "raw_file",
                           value.name = "intensity")
    setkey(data_abundance, raw_file, old_UniqueSequenceID)
    data_sample <- data_sample[data_abundance]
    data_sample[old_FoundinSamples==0 & intensity>0, old_FoundinSamples := 4]
  }
  
  if("old_AbundancesNormalized" %in% names(data_temp)) {
    data_temp[, data_study$SampleIdentifier :=
                transpose(lapply(old_AbundancesNormalized, convertBlobAsDoublePd))]
    data_abundance <- melt(data_temp, measure.vars = data_study$SampleIdentifier,
                           id.vars = c("old_UniqueSequenceID"),
                           variable.name = "raw_file",
                           value.name = "intensity_normalized")
    setkey(data_abundance, raw_file, old_UniqueSequenceID)
    data_sample <- data_sample[data_abundance]
    data_sample[old_FoundinSamples==0 & intensity>0, old_FoundinSamples := 4]
  }
  
  #filter samples without PSMs and label PSM qvalues
  if(filterFdr==T) {
    data_sample <- data_sample[old_FoundinSamples>0]
  }
  data_sample[, qvalue_psm_label := factor(old_FoundinSamples, levels = c(3:1, 4),
                                           labels = c("<=0.01", "<=0.05", ">0.05", "MBR"))]
  
  #filter protein groups
  if(level=="protGroups") {
    data_sample <- data_sample[old_IsMasterProtein==0]
  }
  
  #define sample names and reformat columns
  parseSamples(data_sample, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_sample, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_sample, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_sample, fst_backup_path)
  }
  return(data_sample)
}



readDiann <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                      studyFolderName = NULL, searchEngineName = "DIANN",
                      filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  if(filterFdr==F) warning("DIANN report.tsv is always FDR-filtered, ignoring filterFdr==F")
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  if(loadDecoys==T) warning("DIANN report.tsv is always decoy-filtered, ignoring loadDecoys==T")
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readDiann.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load DIANN data via path to summary.txt
  data_temp <- fread(filePath, showProgress = F)
  
  #rename DIANN names
  setnames(data_temp, paste0("old_", names(data_temp)))
  data_temp[, is_decoy := F]
  
  #create raw file names
  data_temp[, raw_file := factor(old_Run, levels = sort(unique(old_Run)))]
  
  
  #define sample names and reformat columns
  parseSamples(data_temp, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  return(data_temp)
}



readMaxQuant <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                         level = c("psms", "pcms", "ptmGroups"), studyFolderName = NULL, searchEngineName = NULL,
                         filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  level <- match.arg(level)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  if(!is.null(searchEngineName)) {
    stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  } else {
    searchEngineName <- paste0("MaxQuant_", level)
  }
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  if(filterFdr==F) warning("MaxQuant is always FDR-filtered, ignoring filterFdr==F")
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readMaxQuant.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load MaxQuant data
  data_temp <- fread(filePath, showProgress = F)
  if(loadDecoys==F) {
    data_temp <- data_temp[!Reverse %in% "+"]
  }
  
  if(level=="ptmGroups") {
    #create ptm group column
    modifications <- data_temp[Modifications!="Unmodified", .N, by = Modifications][, unlist(strsplit(Modifications, ";"))]
    modifications <- unique(gsub("^\\d{0,2} ", "", modifications))
    parseModificationsColumn(dataTable = data_temp,
                             modifications = modifications,
                             modificationsColumn = "Modifications",
                             sequenceColumn = "Sequence",
                             unimodColumn = "MaxQuant",
                             unimodExtraction = "gsub",
                             unimodLocation = "^(\\d{1,2}) .*$",
                             unimodSeparator = ";")
    
    #rename MaxQuant names by replacing spaces / brackets with dots for R compatibility
    setnames(data_temp, paste0("old_", make.names(names(data_temp), unique = T)))
    
    #parse intensity columns to samples
    data_temp[, old_Intensity := NULL]
    intCols <- names(data_temp)[grepl("^old_Intensity\\..*$", names(data_temp))]
    expCols <- names(data_temp)[grepl("^old_Experiment\\..*$", names(data_temp))]
    data_temp_int <- melt(data_temp, measure.vars = intCols,
                          variable.name = "raw_file", value.name = "old_Intensity")
    setkey(data_temp_int, raw_file, old_ptm_group)
    
    #extract PSM count to filter identifications to file-local level
    data_temp_exp <- melt(data_temp, id.vars = c("old_ptm_group"),
                          measure.vars = expCols, variable.name = "raw_file",
                          value.name = "psm_count")
    data_temp_exp[, raw_file := gsub("old_Experiment", "old_Intensity", raw_file, fixed = T)]
    setkey(data_temp_exp, raw_file, old_ptm_group)
    data_temp <- data_temp_int[data_temp_exp]
    data_temp <- data_temp[!is.na(psm_count)]
    data_temp[old_Intensity==0, old_Intensity := NA]
    
    #remap raw_file onto correct study names
    data_summary <- fread(file.path(dirname(filePath), "summary.txt"), showProgress = F)
    data_summary <- data_summary[`Raw file`!="Total"]
    data_temp[, raw_file := factor(gsub("^old_Intensity\\.", "", raw_file),
                                   data_summary$Experiment, data_summary$`Raw file`)]
    
  } else {
    #rename MaxQuant names by replacing spaces / brackets with dots for R compatibility
    setnames(data_temp, paste0("old_", make.names(names(data_temp), unique = T)))
    #create raw file names
    data_temp[, raw_file := factor(old_Raw.file, levels = sort(unique(old_Raw.file)))]
  }
  
  #define sample names and reformat columns
  filePathOverwrite <- file.path(dirname(filePath), "summary.txt")
  parseSamples(data_temp, sampleNamesPath, filePath, studyFolderName, filePathOverwrite = filePathOverwrite)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  return(data_temp)
}



readSage <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                     studyFolderName = NULL, searchEngineName = "Sage",
                     filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readSage.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load Sage data
  data_temp <- fread(filePath, showProgress = F)
  if(filterFdr==T) {
    data_temp <- data_temp[spectrum_q<=0.01]
  }
  if(loadDecoys==F) {
    data_temp <- data_temp[label==1L]
  }
  setnames(data_temp, paste0("old_", names(data_temp)))
  
  #create raw file names
  data_temp[, old_filename := tools::file_path_sans_ext(old_filename)]
  data_temp[, raw_file := factor(old_filename, levels = sort(unique(old_filename)))]
  
  
  #define sample names and reformat columns
  parseSamples(data_temp, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  return(data_temp)
}



readPercolator <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                           studyFolderName = NULL, searchEngineName = "Percolator_psms",
                           filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readPercolator.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load Percolator data
  data_temp <- fread(filePath, showProgress = F)
  setnames(data_temp, paste0("old_", make.names(names(data_temp), unique = T)))
  data_temp[, is_decoy := F]
  
  if(loadDecoys==T) {
    data_temp_decoy <- fread(gsub("^(.*)percolator_target(.*)$", "\\percolator_decoy\\2", filePath), showProgress = F)
    setnames(data_temp_decoy, paste0("old_", make.names(names(data_temp_decoy), unique = T)))
    data_temp_decoy[, is_decoy := T]
    data_temp <- rbind(data_temp, data_temp_decoy)
  }
  if(filterFdr==T) {
    data_temp <- data_temp[old_q.value<=0.01]
  }
  
  #raw file names are parsed from PSMId
  data_temp[, raw_file := gsub("^(.*)\\.\\d*\\.\\d*\\.\\d_\\d$", "\\1", old_PSMId)]
  data_temp[, raw_file := factor(raw_file, levels = sort(unique(raw_file)))]
  
  #define sample names and reformat columns
  parseSamples(data_temp, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  return(data_temp)
}




readFragPipe <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                         level = c("psms", "pcms", "ptmGroups"), studyFolderName = NULL, searchEngineName = "FragPipe_psms",
                         filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  level <- match.arg(level)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  if(filterFdr==F) warning("FragPipe is always FDR-filtered, ignoring filterFdr==F")
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  if(loadDecoys==T) warning("FragPipe is always decoy-filtered, ignoring loadDecoys==T")
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readFragPipe.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load FragPipe data
  data_temp <- fread(filePath, showProgress = F)
  #rename FragPipe names by replacing spaces / brackets with dots for R compatibility
  setnames(data_temp, paste0("old_", make.names(names(data_temp), unique = T)))
  
  #create raw file names
  if(level %in% "psms") {
    old_filenames <- data_temp[, unique(old_Spectrum.File)]
    new_filenames <- gsub("\\", "/", old_filenames, fixed=T)
    new_filenames <- gsub("^interact-([^.]*?)(_rank\\d){0,1}\\..*$", "\\1", basename(new_filenames))
    data_temp[, raw_file := factor(old_Spectrum.File, levels = old_filenames, labels = new_filenames)]
    
    #only psm.txt stores raw file name in file -> need to reconstruct from fp-manifest for others
  } else if(level %in% "pcms") {
    #search for fp-manifest files in same and upper directory
    filePathManifestSame <- list.files(dirname(filePath), pattern = "^fragpipe-files\\.fp-manifest$", full.names = T)
    filePathManifestUp <- list.files(dirname(dirname(filePath)), pattern = "^fragpipe-files\\.fp-manifest$", full.names = T)
    
    #if fp-manifest in same directory, assign if exactly 1 unique raw file detected
    if(length(filePathManifestSame)==1) {
      manifest <- fread(filePathManifestSame)
      setnames(manifest, c("raw", "experiment", "replicate", "MS"))
      if(manifest[, .N, by = raw][, .N]>1) {
        warning("More than 1 raw file detected in FragPipe manifest, cannot separate raw files")
      }
      raw_file_name <- manifest[, unique(raw)][1]
      data_temp[, raw_file := tools::file_path_sans_ext(basename(gsub("\\", "/", raw_file_name, fixed=T)))]
      data_temp[, raw_file := factor(raw_file, levels = sort(unique(raw_file)))]
      
      #if fp-manifest one level higher -> match correct folder
    } else if(length(filePathManifestUp)==1) {
      manifest <- fread(filePathManifestUp)
      setnames(manifest, c("raw", "experiment", "replicate", "MS"))
      manifest[, folder := apply(.SD, 1, function(x) paste(x[!is.na(x)], collapse = "_") ),
               .SDcols = c("experiment", "replicate")]
      raw_file_name <- manifest[folder %in% basename(dirname(filePath)), raw]
      if(length(raw_file_name)>1) {
        warning("More than 1 raw file detected in FragPipe manifest, cannot separate raw files")
      }
      data_temp[, raw_file := tools::file_path_sans_ext(basename(gsub("\\", "/", raw_file_name, fixed=T)))]
      data_temp[, raw_file := factor(raw_file, levels = sort(unique(raw_file)))]
      
      #set to NA if fp-manifest cannot be found
    } else {
      data_temp[, raw_file := NA]
    }
  }
  
  
  #to define sample names and reformat columns, detect if file 1 folder higher than folder with fp-manifest
  if(length(list.files(dirname(dirname(filePath)), pattern = "^fragpipe-files\\.fp-manifest$"))>0) {
    filePathExperiment <- gsub("^(.*)/([^/]*/[^/]*)$", "\\1_\\2", filePath)
    filePathOverwrite <- file.path(dirname(filePath), "fragpipe-files.fp-manifest")
    parseSamples(dataTable = data_temp, sampleNamesPath, filePathExperiment, studyFolderName, filePathOverwrite = filePathOverwrite)
    
    #if not, assume file is is in same folder as fp-manifest
  } else {
    filePathOverwrite <- file.path(dirname(filePath), "fragpipe-files.fp-manifest")
    parseSamples(dataTable = data_temp, sampleNamesPath, filePath, studyFolderName, filePathOverwrite = filePathOverwrite)
  }
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  return(data_temp)
}



readMetamorpheus <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                             studyFolderName = NULL, searchEngineName = "Metamorpheus_psms",
                             filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T, loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readMetamorpheus.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #NOTE: PSMs are not unique, as the notch column (1 = isotope error) allows matching the same PSM to one scan with and without isotope error
  #But, for 102688 PSMs only impacted 1 PSM
  
  #load Metamorpheus data
  data_temp <- fread(filePath, showProgress = F)
  setnames(data_temp, paste0("old_", make.names(names(data_temp), unique = T)))
  data_temp[, is_decoy := F]
  
  if(loadDecoys==F) {
    data_temp <- data_temp[old_Decoy=="N"]
  }
  if(filterFdr==T) {
    data_temp <- data_temp[old_QValue<=0.01]
  }
  
  #raw file names are not parsed, set to NA
  data_temp[, raw_file := gsub("-calib", "", old_File.Name, fixed = T)]
  data_temp[, raw_file := factor(raw_file, levels = sort(unique(raw_file)))]
  
  #define sample names and reformat columns
  parseSamples(data_temp, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  return(data_temp)
}



readMsgfPlus <- function(filePath, sampleNamesPath = NULL, parsingRulesPath = here::here("scripts/column_names.csv"),
                         studyFolderName = NULL, searchEngineName = "MSGFplus",
                         filterFdr = T, loadDecoys = F, loadBackup = T, writeBackup = T,
                         loadSlow = T, loadPtms = F, localizeUnimod = c("unimod:21"), benchmark = F) {
  #check input arguments
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  if(!is.null(studyFolderName)) {
    stopifnot(is.character(studyFolderName), length(studyFolderName)==1L)
  }
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(filterFdr), length(filterFdr)==1L)
  stopifnot(is.logical(loadDecoys), length(loadDecoys)==1L)
  stopifnot(is.logical(loadBackup), length(loadBackup)==1L)
  stopifnot(is.logical(writeBackup), length(writeBackup)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1L)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.character(localizeUnimod), length(localizeUnimod)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  
  #check if fst backup is available and load instead of reprocessing everything
  fst_backup_path <- file.path(dirname(filePath), paste0(tools::file_path_sans_ext(basename(filePath)), ".readMsgfPlus.fst"))
  if(loadBackup==T && file_test("-f", fst_backup_path)) {
    return(readFstBackup(fst_backup_path))
  }
  
  #load MS-GF-plus data
  data_temp <- fread(filePath, showProgress = F)
  #rename MS-GF-plus names by replacing spaces / brackets with dots for R compatibility
  setnames(data_temp, paste0("old_", make.names(names(data_temp), unique = T)))
  if(filterFdr==T) {
    data_temp <- data_temp[old_QValue<=0.01]
  }
  if(loadDecoys==F) {
    data_temp <- data_temp[!grepl("^XXX_.*$", old_Protein)]
  }
  
  #create raw file names
  data_temp[, raw_file := tools::file_path_sans_ext(old_X.SpecFile)]
  data_temp[, raw_file := factor(raw_file, levels = sort(unique(raw_file)))]
  
  
  #define sample names and reformat columns
  parseSamples(data_temp, sampleNamesPath, filePath, studyFolderName)
  warningFileName <- file.path(basename(dirname(filePath)), basename(filePath))
  parseColumns(data_temp, parsingRulesPath, searchEngineName,
               loadSlow, loadPtms, benchmark, warningFileName)
  
  #MSGFplus reports multiple rows per PSM, 1 per protein -> collapse into 1 PSM each
  names_temp <- names(data_temp)[!names(data_temp) %in% c("psm", "protein_group")]
  data_sd <- data_temp[, .SD[1], by = psm, .SDcols = names_temp]
  data_pg <- data_temp[, .(protein_group = paste(sort(unique(protein_group)), collapse = ";")),
                       by = psm][, .(protein_group)]
  data_temp <- cbind(data_sd, data_pg)
  
  if(loadPtms==T) {
    localizePtms(data_temp, localizeUnimod)
  }
  #write fst backup
  if(writeBackup==T) {
    writeFstBackup(data_temp, fst_backup_path)
  }
  return(data_temp)
}



listFileTypes <- function(level = c("study", "psms", "pcms", "ptmGroups", "proteins", "protGroups"), mergeAll = T) {
  level <- match.arg(level)
  stopifnot(is.logical(mergeAll), length(mergeAll)==1L)
  
  if(level == "study") {
    files <- c(
      "pdResult" = ".*\\.pdResult",
      "DIANN" = ".*(report|diann-output)\\.tsv",
      "MaxQuant" = "summary\\.txt",
      "FragPipe" = "fragpipe-files\\.fp-manifest",
      "Sage" = ".*\\.sage\\.tsv",
      "Percolator" = ".*percolator_target_psms\\.tsv",
      "Metamorpheus" = "AllPSMs\\.psmtsv",
      "MSGFplus" = ".*\\.raw\\.tsv"
    )
    
  } else if(level == "psms") {
    files <- c(
      "pdResult" = ".*\\.pdResult",
      "MaxQuant" = "msms\\.txt",
      "FragPipe" = "psm\\.tsv",
      "Sage" = ".*\\.sage\\.tsv",
      "Percolator" = ".*percolator_target_psms\\.tsv",
      "Metamorpheus" = "AllPSMs\\.psmtsv",
      "MSGFplus" = ".*\\.raw\\.tsv"
    )
    
  } else if(level == "pcms") {
    files <- c(
      "DIANN" = ".*(report|diann-output)\\.tsv",
      "MaxQuant" = "evidence\\.txt",
      "FragPipe" = "ion\\.tsv"
    )
    
  } else if(level == "ptmGroups") {
    files <- c(
      "pdResult" = ".*\\.pdResult",
      "MaxQuant" = "modificationSpecificPeptides\\.txt"
    )
    
  }  else if(level == "proteins") {
    files <- c(
      "pdResult" = ".*\\.pdResult"
    )
    
  }  else if(level == "protGroups") {
    files <- c(
      "pdResult" = ".*\\.pdResult"
    )
  }
  
  #add string start and stop, as well as merged file vector
  files_named <- setNames(paste0("^", files, "$"), names(files))
  if(mergeAll==T) {
    files_all <- paste0("^(", paste(files, collapse = "|"), ")$")
    return(c(setNames(files_all, "all"), files_named))
  } else {
    return(files_named)
  }
}



#delete all .fst backup files
deleteFstBackup <- function(dataPath) {
  #detect if input is single folder path or file vector
  isFolderPath <- (file_test("-d", dataPath) && length(dataPath)==1L)
  stopifnot(isFolderPath==T || all(file_test("-f", dataPath)))
  #identify fst backup files
  if(isFolderPath==T) {
    file_list <- detectFiles(dataPath, "\\.fst$", ignore = F)
  } else {
    file_list <- lapply(dataPath, function(path) {
      list.files(dirname(path), full.names = T,
                 pattern = gsub("\\.[^\\.]*$", ".*.fst$", basename(path)))
    })
    file_list <- unique(unlist(file_list))
    if(length(file_list)==0) stop("no .fst files detected")
  }
  message("Found ", length(file_list), " .fst backup files")
  user_prompt <- readline("Delete all files? [y/n]")
  if(user_prompt=="y") {
    file_delete <- file.remove(file_list)
    if(all(file_delete)) {
      message("Successfully deleted all files")
    } else {
      message("Could not delete file(s)")
      print(file_list[!file_delete])
    }
  }
}


#list all files according to pattern
detectFiles <- function(dataPath, dataPattern, recursive = T, ignore = T) {
  stopifnot(is.character(dataPath))
  stopifnot(is.character(dataPattern), length(dataPattern)==1L)
  stopifnot(is.logical(recursive), length(recursive)==1L)
  stopifnot(is.logical(ignore), length(ignore)==1L)
  
  #test if dataPath is folder; if yes search for files, if no return unchanged
  if((all(file_test("-d", dataPath)) && length(dataPath)==1L)) {
    dataPath <- list.files(dataPath, pattern = dataPattern, recursive = recursive, full.names = T)
    if(ignore==T) {
      dataPath <- dataPath[!grepl("_ignore", dataPath)]
    }
    if(length(dataPath)==0) {
      stop(paste("no", dataPattern, "files detected in folder path"))
    }
  }
  
  #remove double slashes and duplicated files
  dataPath <- gsub("//", "/", dataPath, fixed = T)
  dataPath <- dataPath[!duplicated(dataPath)]
  return(dataPath)
}


#create longest shared dataPath from vector of filePaths
createSharedFilePath <- function(filePaths) {
  #if single folder path given, return and exit
  if(all(file_test("-d", filePaths)) && length(filePaths)==1L) {
    return(filePaths)
  }
  filePaths <- normalizePath(dirname(filePaths))
  lsDir <- transpose(strsplit(filePaths, "(?<=.)/", perl=T))
  lsDirN <- sapply(lsDir, uniqueN)
  lsDirShared <- sapply(lsDir[which(cumsum(lsDirN)==seq_along(lsDirN))], unique)
  dataPath <- do.call(file.path, as.list(lsDirShared))
  if(length(dataPath)==0) {
    warning("No shared common folder path detected")
  }
  return(dataPath)
}


#if dataPath is a path, create study_folder as highest common path
createStudyFolder <- function(filePath, dataPath, dataPathIsPath) {
  stopifnot(file_test("-f", filePath), length(filePath)==1L)
  stopifnot(file_test("-d", dataPath), length(dataPath)==1L)
  stopifnot(is.logical(dataPathIsPath), length(dataPathIsPath)==1L)
  
  #if dataPath was provided, remove from filePath to create shorter study_folder
  filePath <- normalizePath(filePath)
  if(dataPathIsPath==T) {
    dataPath <- normalizePath(dataPath)
    study_folder <- gsub("^/(.*)$", "\\1", gsub(dataPath, "", dirname(filePath)))
    
  } else {
    #extract first folder upstream of file
    study_folder <- basename(dirname(filePath))
  }
  return(study_folder)
}


readSearchEngine <- function(filePath, sampleNamesPath = NULL, parsingRulesPath, studyFolderName = NULL, level, ...) {
  stopifnot(file_test("-f", filePath))
  if(!is.null(sampleNamesPath)) {
    stopifnot(file_test("-f", sampleNamesPath), length(sampleNamesPath)==1L)
  }
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  level <- match.arg(level, choices = c("psms", "pcms", "ptmGroups", "proteins", "protGroups"))
  
  #identify the search engine matching the file path
  fileTypes <- listFileTypes(level, mergeAll = F)
  file_match <- sapply(fileTypes, grepl, basename(filePath))
  #error if not exactly 1 match returned
  if(sum(file_match)==0) {
    stop(paste("no sample reading function defined for file", basename(filePath)))
  } else if(sum(file_match)>1) {
    stop(paste("more than one reading function matched for file", basename(filePath)))
  }
  
  #call the respective search engine reader function
  search_engine <- names(which(file_match))
  
  if(search_engine=="pdResult") {
    if(level=="psms") {
      readPdPsms(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, ...)
    } else if(level=="ptmGroups") {
      readPdPtmGroups(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, ...)
    } else if(level=="proteins") {
      readPdProteins(filePath, sampleNamesPath, parsingRulesPath, level, studyFolderName, ...)
    } else if(level=="protGroups") {
      readPdProteins(filePath, sampleNamesPath, parsingRulesPath, level, studyFolderName, ...)
    }
    
  } else if(search_engine=="DIANN") {
    readDiann(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, ...)
    
  } else if(search_engine=="MaxQuant") {
    readMaxQuant(filePath, sampleNamesPath, parsingRulesPath,
                 level, studyFolderName, paste0("MaxQuant_", level), ...)
    
  } else if(search_engine=="FragPipe") {
    readFragPipe(filePath, sampleNamesPath, parsingRulesPath,
                 level, studyFolderName, paste0("FragPipe_", level), ...)
    
  } else if(search_engine=="Sage") {
    readSage(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, ...)
    
  } else if(search_engine=="Percolator") {
    readPercolator(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, ...)
    
  } else if(search_engine=="Metamorpheus") {
    readMetamorpheus(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, ...)
    
  } else if(search_engine=="MSGFplus") {
    readMsgfPlus(filePath, sampleNamesPath, parsingRulesPath, studyFolderName, ...)
  }
}



writeFstBackup <- function(dataTable, fstBackupPath) {
  stopifnot(is.data.table(dataTable))
  stopifnot(file_test("-d", dirname(fstBackupPath)), length(fstBackupPath)==1L)
  write_fst(dataTable, fstBackupPath)
}



readFstBackup <- function(fstBackupPath) {
  stopifnot(file_test("-f", fstBackupPath), length(fstBackupPath)==1L)
  message(" (INFO: existing backup loaded)", appendLF = F)
  return(read_fst(fstBackupPath, as.data.table = T))
}



#BLOB conversion functions adapted for Proteome Discoverer (include unwanted separator that needs to be removed)
convertBlobAsDoublePd <- function (hexBlob) {
  if(all(is.na(hexBlob)) || length(hexBlob)==0) {
    return(NA)
  } else {
    rawbytes <- hexBlob[rep(c(rep(T, 8), F), length(hexBlob)/9)]
    return(readBin(rawbytes, "double", n = length(rawbytes)/8, size = 8, endian = "little"))
  }
}



convertBlobAsFloatPd <- function (hexBlob) {
  rawbytes <- hexBlob[rep(c(rep(T, 4), F), length(hexBlob)/5)]
  return(readBin(rawbytes, "double", n = length(rawbytes)/4, size = 4, endian = "little"))
}



convertBlobAsIntegerPd <- function (hexBlob) {
  rawbytes <- hexBlob[rep(c(rep(T, 4), F), length(hexBlob)/5)]
  return(readBin(rawbytes, "integer", n = length(rawbytes)/4, size = 4, endian = "little"))
}



parseSamples <- function(dataTable, sampleNamesPath = NULL, filePath, studyFolderName = NULL, filePathOverwrite = NULL) {
  #define study columns and column names for reordering
  if(is.null(studyFolderName)) {
    dataTable[, study_folder := factor(basename(dirname(filePath)))]
  } else {
    dataTable[, study_folder := factor(studyFolderName)]
  }
  dataTable[, study_file := factor(basename(filePath))]
  
  if(!is.null(sampleNamesPath)) {
    #define sample, condition and replicate
    data_samples <- readSampleNames(filePath = sampleNamesPath)
    #in case study_file in sample_names was generated by different file name, replace to enable matching
    if(!is.null(filePathOverwrite)) {
      dataTable[, sample := paste(study_folder, basename(filePathOverwrite), raw_file, sep = "_")]
    } else {
      dataTable[, sample := paste(study_folder, study_file, raw_file, sep = "_")]
    }
    
    #check for correct matching of sample names
    sample_match <- dataTable[, unique(sample)] %in% data_samples$sample_name_old
    if(any(!sample_match)) {
      stop(paste("Could not find old sample names",
                 paste(dataTable[, unique(sample)][!sample_match], collapse=", ")))
    }
    
    dataTable[, sample := factor(sample, levels = data_samples$sample_name_old, labels = data_samples$sample_name)]
    dataTable[, condition := factor(gsub("^(.*)_\\d{1,2}$", "\\1", sample), levels = unique(data_samples$condition_name))]
    rep_vec <- as.integer(gsub("^.*_(\\d{1,2})$", "\\1", dataTable$sample))
    dataTable[, replicate := factor(rep_vec, levels = 1L:max(rep_vec))]
  }
  
  #remove TMT/SILAC channels for raw file numbering
  dataTable[, raw_file := factor(gsub("^(.*) - .*$", "\\1", raw_file))]
  dataTable[, raw_number := factor(raw_file, levels = levels(raw_file),
                                   labels = 1L:length(unique(raw_file)))]
}



parseReorder <- function(dataTable, parsingRulesPath) {
  #reorder columns and set key
  col_base <- c("sample", "condition", "replicate", "study_folder", "study_file", "raw_file", "raw_number")
  data_parsing <- fread(parsingRulesPath, na.strings = "")
  col_names <- c(col_base, data_parsing[!is.na(column), column])
  setcolorder(dataTable, col_names[col_names %in% names(dataTable)])
  set_cols <- c("psm", "pcm", "ptm", "peptide", "protein")
  set_col <- set_cols[which(set_cols %in% names(dataTable))[1]]
  if("sample" %in% names(dataTable)) {
    setkeyv(dataTable, c("sample", set_col))
  } else {
    setkeyv(dataTable, c("study_folder", "study_file", "raw_file", set_col))
  }
}



parseColumns <- function(dataTable, parsingRulesPath, searchEngineName,
                         loadSlow=F, loadPtms=F, benchmark=F, warningFileName) {
  #check input arguments
  stopifnot(is.data.table(dataTable))
  stopifnot(file_test("-f", parsingRulesPath), length(parsingRulesPath)==1L)
  stopifnot(is.character(searchEngineName), length(searchEngineName)==1L)
  stopifnot(is.logical(loadSlow), length(loadSlow)==1)
  stopifnot(is.logical(loadPtms), length(loadPtms)==1)
  stopifnot(is.logical(benchmark), length(benchmark)==1)
  stopifnot(is.character(warningFileName), length(warningFileName)==1)
  
  #load and filter column parsing rules file
  data_parsing <- fread(parsingRulesPath, na.strings = "")
  data_parsing <- data_parsing[!is.na(get(searchEngineName))]
  if(data_parsing[, .N]==0) {
    stop(paste("no reformatting entries found for", searchEngineName))
  }
  if(any(!c("column", searchEngineName) %in% names(data_parsing))) {
    stop(paste("column and", searchEngineName, "need to be part of", parsingRulesPath))
  }
  if(loadPtms==F) {
    data_parsing <- data_parsing[PTM==F]
  }
  if(loadSlow==F) {
    data_parsing <- data_parsing[slow==F]
  }
  
  #extract and reformat function calls to format dataTable[, column := expression]
  str_expression <-
    data_parsing[, sprintf("dataTable[, %s := %s][]", column, get(searchEngineName))]
  str_column <- data_parsing[, column]
  #extract required input columns to enable column availability checking
  str_variables <- lapply(str_expression, function(x) all.vars(parse(text = x), unique = F))
  str_variables <- lapply(str_variables, function(x) {
    x <- x[c(-1, -2)]
    unique(x[!x %in% c("T", "F", "x", "convertBlobAsDoublePd", "modSeqPdPsm")])
  })
  names(str_variables) <- data_parsing[, column]
  
  #track missing columns and time
  list_missing <- list()
  if(benchmark==T) {
    time_table <- data.table(column = "start",
                             time = Sys.time())
  }
  
  #execute expressions to create columns 1-by-1
  for(i in 1L:data_parsing[, .N]) {
    #check if variables missing from dataTable
    if(any(!str_variables[[i]] %in% names(dataTable))) {
      list_missing[i] <- list(str_variables[[i]])
    } else {
      eval(parse(text = str_expression[i]))
    }
    if(benchmark==T) {
      time_table <-
        rbind(time_table,
              data.table(column = str_column[i], time = Sys.time()))
    }
  }
  if(benchmark==T) {
    time_table <- time_table[-1, .(column, time_diff = round(time - time_table$time[-.N], 2))]
    print(time_table[order(time_diff)])
  }
  
  #report missing columns as warning
  if(length(list_missing)>0) {
    is_missing <- !sapply(list_missing, is.null)
    input_missing <- unique(unlist(list_missing[is_missing]))
    column_missing <- str_column[which(is_missing)]
    warning(paste("In", warningFileName,
                  paste(input_missing, collapse = ", "),
                  "not available\nSkipping columns",
                  paste(column_missing, collapse = ", ")),
            call. = F)
  }
  
  #delete old column names, reorder columns and set key
  col_delete <- grep("^old_.*$", names(dataTable), value = T)
  dataTable[, (col_delete) := NULL][]
  parseReorder(dataTable, parsingRulesPath)
}


parsePtmUnimod <- function(ptm, searchEngine = c("PD", "Comet", "Thesaurus", "Sage", "MSFragger", "MSFragger_alt", "MaxQuant", "Spectronaut", "DIANN", "Percolator_FragPipe", "Percolator_Metamorpheus", "MSGFplus", "Metamorpheus")){
  stopifnot(is.character(ptm))
  searchEngine <- match.arg(searchEngine)
  
  #Proteome Discoverer
  if(searchEngine=="PD"){
    ptm <- gsub("c", "C", ptm, fixed=T)
    ptm <- gsub("m", "M[unimod:35]", ptm, fixed=T)
    ptm <- gsub("s", "S[unimod:21]", ptm, fixed=T)
    ptm <- gsub("t", "T[unimod:21]", ptm, fixed=T)
    ptm <- gsub("y", "Y[unimod:21]", ptm, fixed=T)
    ptm <- gsub("b", "B", ptm, fixed=T)
    ptm <- gsub("x", "X", ptm, fixed=T)
    ptm <- gsub("z", "Z", ptm, fixed=T)
  }
  
  #Comet
  if(searchEngine=="Comet"){
    ptm <- gsub("n[42.0106]", "[unimod:1]-", ptm, fixed = T)
    ptm <- gsub("[15.9949]", "[unimod:35]", ptm, fixed = T)
  }
  
  #Thesaurus
  if(searchEngine=="Thesaurus"){
    ptm <- gsub("[+79.966331]", "[unimod:21]", ptm, fixed = T)
    ptm <- gsub("(.)\\[\\+42.010565\\]", "[unimod:1]-\\1", ptm)
    ptm <- gsub("[+15.994915]", "[unimod:35]", ptm, fixed = T)
    ptm <- gsub("[+58.00548]", "[unimod:4]", ptm, fixed = T)
    ptm <- gsub("(.)\\[\\+121.976896\\]", "[unimod:1]-\\1[unimod:21]", ptm)
  }
  
  #Sage
  if(searchEngine=="Sage"){
    ptm <- gsub("C[+57.0215]", "C", ptm, fixed=T)
    ptm <- gsub("[+42.0106]-M[-131.0405]", "[unimod:766]-", ptm, fixed=T) #Met loss, then acetylation
    ptm <- gsub("M[-131.0405]", "[unimod:765]-", ptm, fixed=T) #Met loss
    ptm <- gsub("[+42.0106]", "[unimod:1]", ptm, fixed=T)
    ptm <- gsub("[+57.0215]", "[unimod:4]", ptm, fixed=T)
    ptm <- gsub("[+15.9949]", "[unimod:35]", ptm, fixed=T)
    ptm <- gsub("[+79.9663]", "[unimod:21]", ptm, fixed=T)
  }
  
  #MSFragger
  if(searchEngine=="MSFragger"){
    ptm <- gsub("n[43]", "[unimod:1]-", ptm, fixed=T) #n-terminal acetylation
    ptm <- gsub("n[230]", "[unimod:737]-", ptm, fixed=T) #n-terminal TMT10plex
    ptm <- gsub("S[316]", "S[unimod:737]", ptm, fixed=T) #S TMT10plex
    ptm <- gsub("n[305]", "[unimod:2016]-", ptm, fixed=T) #n-terminal TMTpro
    ptm <- gsub("S[391]", "[unimod:2016]", ptm, fixed=T) #S TMTpro
    ptm <- gsub("C[160]", "C[unimod:0]", ptm, fixed=T) #alkylation (bug? always as n[230]C[160], never C anywhere else)
    ptm <- gsub("C[143]", "C[unimod:0]", ptm, fixed=T) #alkylation minus ammonia
    ptm <- gsub("M[147]", "M[unimod:35]", ptm, fixed=T) #Met oxidation
    ptm <- gsub("S[167]", "S[unimod:21]", ptm, fixed=T) #S phospho
    ptm <- gsub("T[181]", "T[unimod:21]", ptm, fixed=T) #T phospho
    ptm <- gsub("Y[243]", "Y[unimod:21]", ptm, fixed=T) #Y phospho
    ptm <- gsub("Q[111]", "Q[unimod:28]", ptm, fixed=T) #Gln->pyro-Glu
    ptm <- gsub("E[111]", "E[unimod:27]", ptm, fixed=T) #Glu->pyro-Glu
  }
  
  #MSFragger alternative
  if(searchEngine=="MSFragger_alt"){
    ptm <- gsub("n[42.0106]", "[unimod:1]-", ptm, fixed=T)
    ptm <- gsub("[57.0215]", "[unimod:4]", ptm, fixed=T)
    ptm <- gsub("[15.9949]", "[unimod:35]", ptm, fixed=T)
  }
  
  #MaxQuant
  if(searchEngine=="MaxQuant"){
    ptm <- gsub("_", "", ptm, fixed=T)
    ptm <- gsub("(Oxidation (M))", "[unimod:35]", ptm, fixed=T)
    ptm <- gsub("(Acetyl (Protein N-term))", "[unimod:1]-", ptm, fixed=T)
    ptm <- gsub("(Acetyl (N-term))", "[unimod:1]-", ptm, fixed=T)
    ptm <- gsub("(Acetyl (K))", "[unimod:1]", ptm, fixed=T)
    ptm <- gsub("(Phospho (STY))", "[unimod:21]", ptm, fixed=T)
  }
  
  #Spectronaut
  if(searchEngine=="Spectronaut"){
    ptm <- gsub("[Oxidation (M)]", "[unimod:35]", ptm, fixed=T)
    ptm <- gsub("[Acetyl (Protein N-term)]", "[unimod:1]-", ptm, fixed=T)
    ptm <- gsub("[Phospho (STY)]", "[unimod:21]", ptm, fixed=T)
    ptm <- gsub("[Gln->pyro-Glu]", "[unimod:28]", ptm, fixed=T)
    ptm <- gsub("[Glu->pyro-Glu]", "[unimod:27]", ptm, fixed=T)
    ptm <- gsub("[Carbamidomethyl (C)]", "", ptm, fixed=T)
  }
  
  #DIA-NN
  if(searchEngine=="DIANN"){
    ptm <- gsub("^(UniMod:1)", "[unimod:1]-", ptm, fixed=T)
    ptm <- gsub("(UniMod:1)", "[unimod:1]", ptm, fixed=T)
    ptm <- gsub("(UniMod:4)", "", ptm, fixed=T) #carbamidomethyl
    ptm <- gsub("(UniMod:21)", "[unimod:21]", ptm, fixed=T)
    ptm <- gsub("(UniMod:26)", "[unimod:26]", ptm, fixed=T) #pyro-carbamidomethyl
    ptm <- gsub("(UniMod:27)", "[unimod:27]", ptm, fixed=T) #Glu->pyro-Glu
    ptm <- gsub("(UniMod:28)", "[unimod:28]", ptm, fixed=T) #Gln->pyro-Glu
    ptm <- gsub("(UniMod:35)", "[unimod:35]", ptm, fixed=T)
  }
  
  #Percolator FragPipe
  if(searchEngine=="Percolator_FragPipe"){
    ptm <- substr(ptm, 3, nchar(ptm)-2) #remove first 2 and last 2 characters
    ptm <- gsub("n", "", ptm, fixed = T)
    ptm <- gsub("[57.0215]", "", ptm, fixed = T)
    ptm <- gsub("[42.0106]", "[unimod:1]", ptm, fixed = T)
    ptm <- gsub("[15.9949]", "[unimod:35]", ptm, fixed = T)
  }
  
  #Percolator Metamorpheus
  if(searchEngine=="Percolator_Metamorpheus"){
    ptm <- substr(ptm, 3, nchar(ptm)-2) #remove first 2 and last 2 characters
    ptm <- gsub("[Common Fixed:Carbamidomethyl on C]", "", ptm, fixed = T)
    ptm <- gsub("[Common Variable:Oxidation on M]", "[unimod:35]", ptm, fixed = T)
  }
  
  #Metamorpheus
  if(searchEngine=="Metamorpheus"){
    ptm <- gsub("[Common Fixed:Carbamidomethyl on C]", "", ptm, fixed = T)
    ptm <- gsub("[Common Variable:Oxidation on M]", "[unimod:35]", ptm, fixed = T)
  }
  
  #MS-gf plus
  if(searchEngine=="MSGFplus") {
    ptm <- substr(ptm, 3, nchar(ptm)-2) #remove first and last 2 characters
    ptm <- gsub("\\+57\\.021", "", ptm, fixed=T) #Carbamidomethyl
    ptm <- gsub("\\+15\\.995", "[unimod:35]", ptm, fixed=T)
  }
  
  return(ptm)
}



#create PTM Group column
parseModificationsColumn <- function(dataTable, modifications,
                                     modificationsColumn = "ptm",
                                     sequenceColumn = "peptide",
                                     unimodColumn = "ptm",
                                     unimodExtraction = "grepl",
                                     unimodLocation,
                                     unimodSeparator = "[A-Z]") {
  unimods <- listUnimods()
  if(missing(modifications)) {
    modifications <- dataTable[, sort(unique(unlist(str_extract_all(ptm, "\\[[^]]*\\]"))))]
  }
  mod_check <- modifications %in% unimods[, get(unimodColumn)]
  if(!all(mod_check)) {
    stop("undefined modification in readPdPeptideGroups, please update function listUnimods() for ",
         paste(modifications[!mod_check], collapse = ", "))
  }
  unimods <- unimods[is_fixed==F]
  #parse number of occurences per modification
  for(i in 1L:unimods[, .N]) {
    if(is.na(unimods[i, get(unimodColumn)])) {
      dataTable[, paste0("n_unimod_", unimods[i, unimod]) := paste0("0x", unimods[i, unimod])][]
    } else {
      dataTable[!like(get(modificationsColumn), unimods[i, get(unimodColumn)], fixed = T),
                paste0("n_unimod_", unimods[i, unimod]) := paste0("0x", unimods[i, unimod])][]
      dataTable[like(get(modificationsColumn), unimods[i, get(unimodColumn)], fixed = T),
                paste0("n_unimod_", unimods[i, unimod]) :=
                  parseModifications(string = get(modificationsColumn)[1],
                                     unimodPattern = unimods[i, get(unimodColumn)],
                                     unimod = unimods[i, unimod],
                                     unimodExtraction = unimodExtraction,
                                     unimodLocation = unimodLocation,
                                     unimodSeparator = unimodSeparator),
                by = modificationsColumn][]
    }
  }
  dataTable[, modifications := apply(.SD[1], 1, paste, collapse = "_"),
            .SDcols = paste0("n_unimod_", unimods[, unimod]),
            by = eval(paste0("n_unimod_", unimods[, unimod]))]
  dataTable[, ptm_group := paste0(get(sequenceColumn), "_", modifications)]
  dataTable[, ptm_group_J := gsub("(I|L)", "J", ptm_group)]
  dataTable[, paste0("n_unimod_", unimods[, unimod]) := NULL]
}



parseModifications <- function(string, unimodPattern, unimod,
                               unimodExtraction = c("gsub", "grepl"),
                               unimodLocation, unimodSeparator) {
  unimodExtraction <- match.arg(unimodExtraction)
  sapply(str_split(string, unimodSeparator), function(x) {
    matched_x <- x[grepl(unimodPattern, x, fixed = T)]
    if(unimodExtraction=="gsub") {
      matched_digit <- sprintf("%sx%s",
                               gsub(unimodLocation, "\\1", matched_x),
                               unimod)
      ifelse(grepl(unimodLocation, matched_x),
             matched_digit, sprintf("1x%s", unimod))
    } else if(unimodExtraction=="grepl") {
      sprintf("%sx%s", length(matched_x), unimod)
    }
  })
}



listUnimods <- function() {
  #add unimods below in following order:
  # unimod integer
  # Monoisotopic mass	numeric
  # TRUE for fixed modifications
  # unimod PSI-MS Name (Interim name if NA)
  # unimod ptm column
  # PD peptide groups modifications column
  # MaxQuant Modifications column
  # FragPipe Assigned Modifications column
  unimods <-
    list(c(1L, 42.010565, F, "Acetyl", "[unimod:1]", "Acetyl", "Acetyl (Protein N-term)", "N-term(42.0106)"),
         c(4L, 57.021464, T, "Carbamidomethyl", "[unimod:4]", "Carbamidomethyl", NA, "C(57.0214)"),
         c(5L, 43.005814, F, "Carbamyl", "[unimod:5]", "Carbamyl", NA, "C(43.0058)"),
         c(21L, 79.966331, F, "Phospho", "[unimod:21]", "Phospho", "Phospho (STY)", NA),
         c(35L, 15.994915, F, "Oxidation", "[unimod:35]", "Oxidation", "Oxidation (M)", "M(15.9949)"),
         c(765L, -131.040485, F, "Met-loss", "[unimod:765]", "Met-loss", NA, NA),
         c(766L, -89.029920, F, "Met-loss+Acetyl", "[unimod:766]", "Met-loss+Acetyl", NA, NA))
  unimods <- as.data.table(transpose(unimods))
  names(unimods) <- c("unimod", "mass", "is_fixed", "name", "ptm", "PD_ptmGroups", "MaxQuant", "FragPipe")
  unimods[, unimod := as.integer(unimod)][]
  unimods[, mass := as.numeric(mass)][]
  return(unimods)
}


modSeqPdPsm <- function(modifications, sequence) {
  test_str <- unlist(strsplit(modifications, "; ", fixed=T))
  test_str <- test_str[!grepl("(Carbamidomethyl|TMT)", test_str)]
  test_str <- gsub("^N-Term\\(.*\\)(\\(.*\\))$", "n1\\1", test_str)
  test_str <- sort(unique(as.integer(gsub("^.(.*)\\(.*\\)$", "\\1", test_str))))
  
  id_vec <- seq(nchar(sequence)) %in% test_str
  id_seq <- unlist(strsplit(sequence, ""))
  paste(ifelse(id_vec, tolower(id_seq), id_seq), collapse = "")
}



mapFastaHeader <- function(fastaHeader, extractor = c("GN=", "OS="), sortUnique = T) {
  extractor <- match.arg(extractor)
  fastaHeader <- lapply(stri_match_all_regex(fastaHeader, "(?=([A-Z][A-Z]=[^=]*)(;| [A-Z][A-Z]=|$))"), function(x) x[,2])
  fastaHeader <- lapply(fastaHeader, function(x) gsub(extractor, "", x[grepl(extractor, x, fixed = T)], fixed = T))
  if(sortUnique) {
    fastaHeader <- sapply(fastaHeader, function(x) paste(sort(unique(x)), collapse = "; "))
  } else {
    fastaHeader <- sapply(fastaHeader, function(x) paste(x, collapse = "; "))
  }
  return(fastaHeader)
}



processIntensity <- function(dataTable, intensityColumns = c("intensity_ms2", "intensity_ms1"),
                             transformLog2 = F, groupColumn = "condition", formatReturn = T) {
  stopifnot(is.data.table(dataTable))
  stopifnot(is.character(intensityColumns))
  stopifnot(is.logical(transformLog2), length(transformLog2)==1)
  stopifnot(is.character(groupColumn), length(groupColumn)==1)
  stopifnot(is.logical(formatReturn), length(formatReturn)==1)
  
  #check which intensity columns are available and log2 transform
  intensityColumns <- intensityColumns[intensityColumns %in% names(dataTable)]
  if(length(intensityColumns)==0) {
    stop("none of the provided intensity columns detected")
  }
  
  #create intensity column if not yet available
  if(!"intensity" %in% names(dataTable)) {
    dataTable[, intensity := numeric()]
  }
  
  #if at least one provided intensity column: consolidate and calculate metrics
  if(length(intensityColumns)>0) {
    #log2 transform intensity columns
    if(transformLog2==T) {
      dataTable[, (intensityColumns) := lapply(.SD, log2), .SDcols = intensityColumns]
    }
    
    #calculate number of valid intensity values before and after intensity column consolidation
    count_pre <- dataTable[, .(intensity_pre = validN(intensity)), by=c(groupColumn)]
    count_col <- dataTable[, lapply(.SD, validN),
                           .SDcols = intensityColumns, by=c(groupColumn)]
    
    #if ALL intensity columns are NA, fill in with other intensity columns, in order of vector
    for(i in 1:length(intensityColumns)) {
      condition_na <- dataTable[, .(all_na = all(is.na(intensity))), by=c(groupColumn)][all_na==T, get(groupColumn)]
      dataTable[get(groupColumn) %in% condition_na, intensity := get(intensityColumns[i])][]
    }
    
    count_post <- dataTable[, .(intensity_post = validN(intensity)), by=c(groupColumn)]
  } else {
    count_pre <- dataTable[, .(intensity = validN(intensity)), by=c(groupColumn)]
    count_col <- NULL
    count_post <- NULL
  }
  
  #calculate metrics on intensity column and report
  count_int <- testIntensity(dataTable = dataTable,
                             intensityColumn = "intensity",
                             groupColumn = groupColumn,
                             formatReturn = F)
  count <- cbind(count_pre, count_col[, -1], count_post[, -1], count_int[, -1])
  
  #set invalid intensities NA
  dataTable[!isValidIntensity(intensity), intensity := NA]
  
  if(formatReturn==T) {
    count <- count[, lapply(.SD, format, big.mark=",", nsmall=2, digits=1)]
  }
  return(count[])
}


testIntensity <- function(dataTable, intensityColumn = c("intensity"),
                          groupColumn = "condition", formatReturn = T) {
  stopifnot(is.data.table(dataTable))
  stopifnot(is.character(intensityColumn), length(intensityColumn)==1L)
  stopifnot(is.character(groupColumn), length(groupColumn)==1)
  stopifnot(is.logical(formatReturn), length(formatReturn)==1)
  
  #calculate metrics on intensity column
  count <- dataTable[, .(N = .N,
                         valid = validIntensityN(get(intensityColumn)),
                         all_valid = .N==validIntensityN(get(intensityColumn)),
                         na = invalidN(get(intensityColumn)),
                         nan = invalidN(get(intensityColumn), is.nan),
                         inf = invalidN(get(intensityColumn), is.infinite),
                         is0 = length(which(get(intensityColumn)==0)),
                         is1 = length(which(get(intensityColumn)==1)),
                         min = round(min(validIntensity(get(intensityColumn))), 2),
                         mean = round(mean(validIntensity(get(intensityColumn))), 2),
                         median = round(median(validIntensity(get(intensityColumn))), 2),
                         max = round(max(validIntensity(get(intensityColumn))), 2)),
                     by = c(groupColumn)]
  if(formatReturn==T) {
    count <- count[, lapply(.SD, format, big.mark=",", nsmall=2, digits=1)]
  }
  return(count[])
}


validN <- function(x, fun = is.na) {
  length(which(!fun(x)))
}

invalidN <- function(x, fun = is.na) {
  length(which(fun(x)))
}

isValidIntensity <- function(x) {
  !is.na(x) & !is.nan(x) & !is.infinite(x) & x!=0 & x!=1
}

validIntensity <- function(x) {
  x[!is.na(x) & !is.nan(x) & !is.infinite(x) & x!=0 & x!=1]
}

validIntensityN <- function(x) {
  length(which(!is.na(x) & !is.nan(x) & !is.infinite(x) & x!=0 & x!=1))
}



testContrasts <- function(dataTable, observationColumn, contrasts,
                          output = c("ratios", "list", "observations"),
                          minQuanReplicates = 3L) {
  #check inputs
  stopifnot(is.data.table(dataTable))
  dtColumns <- c("condition", "sample", "organism", observationColumn, "intensity")
  stopifnot(all(dtColumns %in% names(dataTable)))
  stopifnot(is.character(observationColumn), length(observationColumn)==1L)
  stopifnot(is.character(contrasts))
  if(!all(sapply(str_split(contrasts, "_vs_"), length)==2L)) {
    stop("contrasts needs to be character vector of conditions in dataTable separated by '_vs_'")
  }
  if(!all(unlist(str_split(contrasts, "_vs_")) %in% dataTable[, unique(condition)])) {
    stop("contrasts needs to be character vector of conditions in dataTable separated by '_vs_'")
  }
  output <- match.arg(output)
  stopifnot(is.integer(minQuanReplicates), length(minQuanReplicates)==1L)
  
  #separate contrasts, filter intensity and rename observation column
  lsContrast <- strsplit(contrasts, "_vs_")
  dt <- dataTable[, .SD, .SDcols = dtColumns]
  setnames(dt, observationColumn, "observationColumn")
  
  #check if observations unique
  if(!dt[, .N, by=.(sample, observationColumn)][, all(N==1)]) {
    stop(paste("observations are not unique per sample - aggregate onto", observationColumn, "before trying again"))
  }
  #number of quantified replicates
  dt[, nQuan := sum(!is.na(intensity)), by=.(condition, observationColumn)]
  
  #filter to unique set observations making the quan threshold
  dt <- foreach(contrast = lsContrast, .combine = rbind) %do% {
    dtTemp <- dt[condition %in% contrast]
    dtTemp <- dtTemp[dtTemp[, .(nQuanMin = min(nQuan, na.rm = T)), by=.(observationColumn)],
                     on="observationColumn"][nQuanMin>=minQuanReplicates]
    return(dtTemp[, .(condition, sample, organism, observationColumn, intensity, nQuan)])
  }
  dt <- unique(dt)
  setkey(dt, condition, sample, observationColumn)
  
  #calculate ratios
  if(output %in% c("ratios", "list")) {
    dtRatio <- foreach(contrast = lsContrast, contrastLabel = contrasts, .combine = rbind) %do% {
      #create dt without NAs
      dtInt <- dt[condition %in% contrast & !is.na(intensity),
                  .(condition = factor(condition, contrast),
                    sample, organism, observationColumn, intensity = 2^intensity)]
      setkey(dtInt, condition, sample, observationColumn)
      
      #calculate overall mean intensity
      dtMean <- dtInt[, .(intensityMean = log2(mean(intensity))), keyby=.(observationColumn, organism)]
      #calculate condition level count, mean and sd
      dtCond <- dtInt[, .(.N, meanInt = mean(intensity), sdInt = sd(intensity)),
                      keyby=.(condition, observationColumn)]
      #calculate nQuan, CVs
      dtCvA <- dtCond[condition==contrast[1], .(nQuan1st = N, cv1st = sdInt/meanInt),
                      keyby=observationColumn]
      dtCvB <- dtCond[condition==contrast[2], .(nQuan2nd = N, cv2nd = sdInt/meanInt),
                      keyby=observationColumn]
      #calculate log2 ratios
      dtRatio <- dtCond[, .(contrast = factor(contrastLabel, contrasts),
                            ratio = log2(meanInt[1]/meanInt[2])), keyby=observationColumn]
      #merge all
      dtRatio <- dtCvA[dtCvB][dtRatio][dtMean]
      
      #add NA IDs back if no filtering desired
      if(minQuanReplicates==0) {
        dtFill <- dt[condition %in% contrast,
                     .(contrast = factor(contrastLabel, contrasts),
                       observationColumn = unique(observationColumn))]
        dtRatio <- dtRatio[dtFill, on=c("contrast", "observationColumn")]
        #replace NA counts by 0
        dtRatio[is.na(nQuan1st), nQuan1st := 0L]
        dtRatio[is.na(nQuan2nd), nQuan2nd := 0L]
      }
      #filter on threshold
      dtRatio <- dtRatio[nQuan1st>=minQuanReplicates & nQuan2nd>=minQuanReplicates]
      setcolorder(dtRatio, c("contrast", "organism", "observationColumn", "ratio", "intensityMean",
                             "nQuan1st", "nQuan2nd", "cv1st", "cv2nd"))
    }
    setkey(dtRatio, contrast, observationColumn)
    setnames(dtRatio, "observationColumn", observationColumn)
  }
  
  #return results
  setnames(dt, "observationColumn", observationColumn)
  if(output=="list") {
    return(list("ratios" = dtRatio[], "observations" = dt[]))
  } else if(output=="ratios") {
    return(dtRatio[])
  } else if(output=="observations") {
    return(dt[])
  }
}

mergeSpectra <- function(dtTruth, dtObserved,
                         ppmTolerance = 20, ppmShift = NULL,
                         matchType = c("intensity", "tolerance", "tolerance (matched observed)"),
                         mzTruthName = "mz", mzObservedName = "mz",
                         scanTruthName = "scan_ms2", scanObservedName = "scan_ms2",
                         intensityObservedName = "intensity") {
  stopifnot(is.character(mzTruthName), length(mzTruthName)==1L)
  stopifnot(is.data.table(dtTruth), c(scanTruthName, mzTruthName) %in% names(dtTruth))
  stopifnot(is.character(mzObservedName), length(mzObservedName)==1L)
  stopifnot(is.data.table(dtObserved), c(scanObservedName, mzObservedName) %in% names(dtObserved))
  stopifnot(is.numeric(ppmTolerance), length(ppmTolerance)==1L)
  if(!is.null(ppmShift)) {
    stopifnot(is.numeric(ppmShift), length(ppmShift)==1L)
  } else {
    ppmShift <- 0
  }
  matchType <- match.arg(matchType)
  stopifnot(is.character(scanTruthName), length(scanTruthName)==1L)
  stopifnot(is.character(scanObservedName), length(scanObservedName)==1L)
  stopifnot(is.character(intensityObservedName), length(intensityObservedName)==1L)
  
  #subset and reformat truth
  dtTruth <- dtTruth[, .SD, .SDcols = c(scanTruthName, mzTruthName)]
  setnames(dtTruth, c(scanTruthName, mzTruthName), c("scan", "mz"))
  setkey(dtTruth, scan, mz)
  
  #subset and reformat observed
  if(matchType == "intensity") {
    stopifnot(intensityObservedName %in% names(dtObserved))
    dtObserved <- dtObserved[, .SD, .SDcols = c(scanObservedName, mzObservedName, intensityObservedName)]
    setnames(dtObserved, intensityObservedName, "intensity")
  } else {
    dtObserved <- dtObserved[, .SD, .SDcols = c(scanObservedName, mzObservedName)]
    dtObserved[, intensity := numeric()]
  }
  setnames(dtObserved, c(scanObservedName, mzObservedName), c("scan", "mzObserved"))
  setkey(dtObserved, scan, mzObserved)
  
  
  #calculate tolerances on truth data
  dtTruth[, mzLower := mz + (mz*(ppmShift - ppmTolerance)/1e6)]
  dtTruth[, mzUpper := mz + (mz*(ppmShift + ppmTolerance)/1e6)]
  
  #merge and group by mz
  dtMerge <- dtObserved[dtTruth, .(scan, mzMatch=x.mzObserved, mz, intensity),
                        on=.(scan==scan, mzObserved>=mzLower, mzObserved<=mzUpper)]
  
  if(matchType == "intensity") {
    dtMerge <- dtMerge[!is.na(mz), .(mzMatchedN = uniqueN(mzMatch, na.rm = T),
                                     mzMatch = mzMatch[which.max(intensity)]),
                       by = .(scan, mz)]
    
  } else if(matchType %in% c("tolerance", "tolerance (matched observed)")) {
    if(matchType == "tolerance (matched observed)") {
      setnames(dtMerge, c("mz", "mzMatch"), c("mzMatch", "mz"))
    }
    
    dtMerge <- dtMerge[!is.na(mz), .(mzMatchedN = uniqueN(mzMatch, na.rm = T),
                                     mzMatch = mzMatch[which.min(abs(mz-mzMatch))]),
                       by = .(scan, mz)]
  }
  
  #calculate mz match and difference
  dtMerge[, mzMatched := mzMatchedN>0]
  dtMerge[, mzDif := mz-mzMatch]
  dtCols <- c("scan", "mz", "mzMatched", "mzMatchedN", "mzMatch", "mzDif")
  setcolorder(dtMerge, dtCols)
  setkey(dtMerge, scan, mz)
  
  return(dtMerge[])
}