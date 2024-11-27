#libraries
require(data.table)
require(DBI)
require(fst)
require(ggplot2)
require(scales)
require(pbmcapply)
require(Biostrings)
require(stringr)

#functions
replaceChar <- function(SEQUENCE, POSITION, RESIDUE_REPLACE) {
  POSITION <- unlist(POSITION)
  for (i in seq_along(POSITION)) {
    substr(SEQUENCE, POSITION[i], POSITION[i]) <- RESIDUE_REPLACE
  }
  return(SEQUENCE)
}

convertBlobAsDoublePd <- function (hexBlob) {
  if(is.na(hexBlob) || length(hexBlob)==0) {
    return(NA)
  } else {
    rawbytes <- hexBlob[rep(c(rep(T, 8), F), length(hexBlob)/9)]
    return(readBin(rawbytes, "double", n = length(rawbytes)/8, size = 8, endian = "little"))
  }
}

convertIntensityBlobs <- function(rawbytes) 
{
  n <- length(rawbytes)/4
  return(readBin(rawbytes, "double", n = n, size = 4, endian = "little"))
}

convertMs8Feature7 <- function(double){
  readBin(writeBin(double, raw(), size = 4), "double", 1, size = 4, endian = "little")
}

alldup <- function(x) {duplicated(x) | duplicated(x, fromLast = T)}

readProteinsFromFastas <- function(pathToFasta){
  cat('Reading in data from fasta\n')
  
  fasta <- readAAStringSet(pathToFasta)
  sequencesFromFasta <- as.data.table(fasta)
  setnames(sequencesFromFasta, "x", "Sequence")
  
  proteinsFromFasta <- data.table(names(fasta))
  setnames(proteinsFromFasta, "V1", "FastaTitleLines")
  proteinsFromFasta <- cbind(proteinsFromFasta, sequencesFromFasta)
  
  proteinsFromFasta[grepl('Random_[0-9]', FastaTitleLines),Accession:=FastaTitleLines]
  proteinsFromFasta[!grepl('Random_[0-9]', Accession),Accession:=sapply(FastaTitleLines, function(x) str_split(x, "\\|")[[1]][2])]
  
  proteinsFromFasta[grepl('Random_[0-9]', Accession),ORGANISM:='ENTRAPMENT']
  proteinsFromFasta[!grepl('Random_[0-9]', Accession),ORGANISM:=sapply(strsplit(FastaTitleLines, '(>sp|>tr)|\\|| '), function(i) gsub('[A-Za-z0-9]*_', '', i[3]))]
  proteinsFromFasta[grepl('Cont_', Accession),ORGANISM:='CONTAMINANT']
  
  proteinsFromFasta[,Sequence_J:=gsub('I|L', 'J', Sequence)]
  proteinsFromFasta[,ORGANISM_J:=paste0(sort(unique(unlist(strsplit(ORGANISM, ';')))), collapse = ';'), by = Sequence_J]

  proteinsFromFasta[,c('Sequence', 'Sequence_J'):=NULL]
  setcolorder(proteinsFromFasta, c("Accession", "FastaTitleLines", "ORGANISM", "ORGANISM_J"))
  return(proteinsFromFasta)
}

readDigest <- function(pathToDigest) {
  cat('Reading in data from digest\n')
  digest <- fread(pathToDigest, stringsAsFactors = F, integer64 = 'double')
  digest_firstpep <- digest[,.SD[1,],by=list(Protein_Name)]
  digest_firstpep <- digest_firstpep[grepl('^M', Sequence),]
  digest_firstpep[,Sequence:=gsub('^M', '', Sequence)]
  digest_lastpep <- digest[,.SD[.N,],by=list(Protein_Name)]
  digest_lastpep[,Sequence:=paste0(Sequence, 'K')]
  digest <- rbindlist(list(digest, digest_firstpep, digest_lastpep))
  digest <- digest[,list(Protein_Name, Sequence)]
  digest[grepl('Random_[0-9]', Protein_Name),ORGANISM:='ENTRAPMENT']
  digest[grepl('Cont_', Protein_Name),ORGANISM:='CONTAMINANT']
  digest[!(grepl('Random_[0-9]', Protein_Name) | grepl('Cont_', Protein_Name)),ORGANISM:=sapply(strsplit(Protein_Name, '(>sp|>tr)|\\|| '), function(i) gsub('[A-Za-z0-9]*_', '', i[3]))]
  digest[,Protein_Name:=NULL]
  digest_seq_unique <- digest[!alldup(Sequence),]
  digest_seq_nonunique <- digest[alldup(Sequence),]
  digest_seq_nonunique_concat <- digest_seq_nonunique[,list(ORGANISM=paste0(sort(unique(ORGANISM)), collapse = ';')),by=Sequence]
  digest_seq <- rbindlist(list(digest_seq_unique, digest_seq_nonunique_concat))
  
  digest[,Sequence_J:=gsub('I|L', 'J', Sequence)]
  digest_seq_j_unique <- digest[!alldup(Sequence_J),.(Sequence_J, ORGANISM)]
  setnames(digest_seq_j_unique, 'ORGANISM', 'ORGANISM_J')
  digest_seq_j_nonunique <- digest[alldup(Sequence_J),]
  digest_seq_j_nonunique_concat <- digest_seq_j_nonunique[,list(ORGANISM_J=paste0(sort(unique(ORGANISM)), collapse = ';')),by=Sequence_J]
  digest_seq_j <- rbindlist(list(digest_seq_j_unique, digest_seq_j_nonunique_concat))
  
  digest_seq[,Sequence_J:=gsub('I|L', 'J', Sequence)]
  setkey(digest_seq, Sequence_J)
  setkey(digest_seq_j, Sequence_J)
  digest_seq_j <- digest_seq[digest_seq_j]
  return(digest_seq_j)
}

interpolateQvalues <- function(diann, rule = 1) {
  # Q.Value interpolation for decoys
  # Remove all Q_VALUEs for Decoys
  diann[DECOY==1,Q_VALUE:=NA]
  # Set Q_VALUE for DECOYS which get out-of-distribution SVMSCOREs to 1
  diann[DECOY==1 & SVMSCORE==(-10000000.000000),Q_VALUE:=1]
  # Calculate minimum Q_VALUE per unique SAMPLE,SVMSCORE combination
  qval <- diann[!is.na(Q_VALUE) & SVMSCORE!=(-10000000.000000),list(MIN_Q_VALUE=min(Q_VALUE)),by=list(SAMPLE, SVMSCORE)]
  # Find entries without Q_VALUE
  qnas <- diann[is.na(Q_VALUE)]
  # Sort data.tables by SVMSCORE for subsequent interpolation
  setkey(qval, SVMSCORE)
  setkey(qnas, SVMSCORE)
  # Perform SAMPLE-wise interpolation
  sapply(qnas[,unique(SAMPLE)], function(i) qnas[SAMPLE==i,Q_VALUE:=approx(x = qval[SAMPLE==i,SVMSCORE], y = qval[SAMPLE==i,MIN_Q_VALUE], xout = SVMSCORE, ties = 'ordered', method = 'linear', rule = rule)$y])
  # Concatenate results
  diann_inter <- rbindlist(list(diann[!is.na(Q_VALUE)], qnas))
  return(diann_inter)
}

readPdResult_psmGrouper_dda_noConsensusQuan <- function(pathToPdResult_dda, proteinsFromFasta){
  cat('Reading in data from pdresult\n')
  
  pd <- DBI::dbConnect(RSQLite::SQLite(), pathToPdResult_dda, proteinsFromFasta)
  
  # Read in workflow to get rawfile names
  workflow <- setDT(DBI::dbGetQuery(pd, 'SELECT StudyFileID, FileName From WorkflowInputFiles;'))[!is.na(StudyFileID)]
  workflow[,FileName:=gsub('\\..*', '', sapply(strsplit(FileName, '\\\\'), 'tail', 1))]
  setnames(workflow, 'StudyFileID', 'StudyFileId')
  setkey(workflow, StudyFileId)
  
  # Read in target PSMs and PCMs in order to have access to modified sequences and charges
  pdresult_t <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideGroupID, Sequence, ParentProteinAccessions, qValue, SVMScore from TargetPeptideGroups;'))
  pdresult_mapping_t <- setDT(DBI::dbGetQuery(pd, 'SELECT TargetPsmsPeptideID, TargetPeptideGroupsPeptideGroupID From TargetPsmsTargetPeptideGroups;'))
  pdresult_psms_t <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideID, StudyFileId, Modifications, Charge, PrecursorAbundance From TargetPsms;'))
  setkey(pdresult_psms_t, PeptideID)
  setkey(pdresult_mapping_t, TargetPsmsPeptideID)
  pdresult_psms_t[pdresult_mapping_t,PeptideGroupID:=TargetPeptideGroupsPeptideGroupID]
  pdresult_t <- merge(pdresult_t,
                      pdresult_psms_t[,PeptideID:=NULL], 
                      by = 'PeptideGroupID', 
                      all.x = T)
  # pdresult_t <- pdresult_t[Abundances != 0|is.na(Abundances)]
  pdresult_t <- unique(pdresult_t)
  pdresult_t[,DECOY:=0]
  
  # Read in decoy PSMs and PCMs in order to have access to modified sequences and charges
  pdresult_d <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideGroupID, Sequence, ParentProteinAccessions, qValue, SVMScore from DecoyPeptideGroups;'))
  pdresult_mapping_d <- setDT(DBI::dbGetQuery(pd, 'SELECT DecoyPsmsPeptideID, DecoyPeptideGroupsPeptideGroupID From DecoyPsmsDecoyPeptideGroups;'))
  pdresult_psms_d <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideID, StudyFileId, Modifications, Charge From DecoyPsms;'))
  pdresult_psms_d[,PrecursorAbundance:=NA]
  setkey(pdresult_psms_d, PeptideID)
  setkey(pdresult_mapping_d, DecoyPsmsPeptideID)
  pdresult_psms_d[pdresult_mapping_d,PeptideGroupID:=DecoyPeptideGroupsPeptideGroupID]
  pdresult_d <- merge(pdresult_d,
                      pdresult_psms_d[,PeptideID:=NULL], 
                      by = c('PeptideGroupID'), 
                      all.x = T)
  pdresult_d <- unique(pdresult_d)
  pdresult_d[,PeptideGroupID:=paste0(PeptideGroupID, '_DECOY')]
  pdresult_d[,DECOY:=1]
  setcolorder(pdresult_d, colnames(pdresult_t))
  
  # Concatenate target PCMs and decoy PCMs
  pdresult <- rbindlist(list(pdresult_t, pdresult_d))
  
  # Sample name mapping
  setkey(pdresult, StudyFileId)
  pdresult <- workflow[pdresult]
  
  # Modifications
  pdresult[Modifications%like%'M[0-9]{1,2}',Mpos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Oxidation\\)|M(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'M[0-9]{1,2}',Mpos:=NA]
  pdresult[Modifications%like%'S[0-9]{1,2}',Spos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|S(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'S[0-9]{1,2}',Spos:=NA]
  pdresult[Modifications%like%'T[0-9]{1,2}',Tpos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|T(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'T[0-9]{1,2}',Tpos:=NA]
  pdresult[Modifications%like%'Y[0-9]{1,2}',Ypos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|Y(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'Y[0-9]{1,2}',Ypos:=NA]
  pdresult[,ID:=1:nrow(pdresult)]
  pdresult[!is.na(Spos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 's'),by=ID]
  pdresult[!is.na(Tpos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 't'),by=ID]
  pdresult[!is.na(Ypos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 'y'),by=ID]
  pdresult[!is.na(Mpos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Mpos, 'm'),by=ID]
  pdresult[is.na(MODIFIED_SEQUENCE),MODIFIED_SEQUENCE:=Sequence]
  pdresult[,MODIFIED_SEQUENCE:=gsub('C', 'c', MODIFIED_SEQUENCE)]
  
  # ORGANISM mapping
  proteinsFromFasta[grepl('Random_[0-9]', FastaTitleLines),Accession:=FastaTitleLines]
  proteinsFromFasta[!grepl('Random_[0-9]', Accession),Accession:=sapply(FastaTitleLines, function(x) str_split(x, "\\|")[[1]][2])]
  proteinsFromFasta <- rbindlist(list(proteinsFromFasta,
                                      data.table(Accession='mimic',FastaTitleLines='>sp|mimic|mimic_ENTRAPMENT', 
                                                 ORGANISM='ENTRAPMENT', ORGANISM_J='ENTRAPMENT')))
  setkey(proteinsFromFasta, Accession)
  setkey(pdresult, ParentProteinAccessions)
  pdresult[proteinsFromFasta,ORGANISM:=ORGANISM]
  pdresult[grepl(';', ParentProteinAccessions),ORGANISM:=unlist(pbmclapply(strsplit(ParentProteinAccessions, '; '), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM))], collapse = ';'), mc.cores = 10))]
  
  # ENTRAPMENT
  pdresult[,ENTRAPMENT:=ifelse(ORGANISM=='ENTRAPMENT', 1, 0)]
  
  if(!all(pdresult[ENTRAPMENT == 1,ORGANISM] == 'ENTRAPMENT')){
    print("No entrapments found.")
  }
  
  # DECOYS
  pdresult[DECOY == 1,ORGANISM:='DECOY']
  
  # FORMAT ENTRAPMENT PROTEIN ACCESSIONS
  pdresult[,ParentProteinAccessions:=gsub(' ', '', ParentProteinAccessions)]
  pdresult[,ParentProteinAccessions:=gsub('mimic', 'Random_', ParentProteinAccessions)]
  
  pdresult[,ParentProteinAccessions:=gsub(' ', '', ParentProteinAccessions)]
  pdresult[,ParentProteinAccessions:=gsub('mimic', 'Random_', ParentProteinAccessions)]
  
  # PEPTIDE GROUP IDS
  pdresult[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  pdresult[is.na(m_COUNT),m_COUNT:=0]
  pdresult[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  pdresult[is.na(p_COUNT),p_COUNT:=0]
  pdresult[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  pdresult[is.na(c_COUNT),c_COUNT:=0]
  pdresult[,MODPEP_ID:=paste0(Sequence, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  pdresult[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  
  pdresult[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', Charge)]
  pdresult[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  
  # COLUMN CLEAN-UP AND RENAMING
  pdresult[,c('ID', 'Modifications', 'Mpos', 'Spos', 'Tpos', 'Ypos', 'm_COUNT', 'c_COUNT', 'p_COUNT', 'StudyFileId'):=NULL]
  setnames(pdresult, c('FileName', 'PeptideGroupID', 'Sequence', 'ParentProteinAccessions', 'qValue', 'Charge', 'SVMScore'),
           c('SAMPLE', 'ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION', 'Q_VALUE', 'CHARGE', 'SVMSCORE'))
  pdresult[,SOFTWARE:='CHIMERYS']
  setcolorder(pdresult, neworder = c('ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION', 'Q_VALUE', 'SVMSCORE', 'DECOY', 'CHARGE', 'MODIFIED_SEQUENCE', 'SAMPLE', 'ORGANISM', 'ENTRAPMENT', 'PCM_ID', 'PCM_J_ID', 'MODPEP_ID', 'MODPEP_J_ID', 'SOFTWARE'))
  pdresult <- unique(pdresult)
  
  return(pdresult)
  dbDisconnect(pd)
}

readPdResult_psmGrouper_dda_quan_spectralAngles <- function(pathToPdResult, pathToMs8, proteinsFromFasta, pathToPercolatorInput){
  cat('Reading in data from pdresult\n')
  
  pd <- DBI::dbConnect(RSQLite::SQLite(), pathToPdResult)
  ms8 <- DBI::dbConnect(RSQLite::SQLite(), pathToMs8)
  
  # Read in workflow to get rawfile names
  workflow <- setDT(DBI::dbGetQuery(pd, 'SELECT StudyFileID, FileName From WorkflowInputFiles;'))[!is.na(StudyFileID)]
  workflow[,FileName:=gsub('\\..*', '', sapply(strsplit(FileName, '\\\\'), 'tail', 1))]
  setnames(workflow, 'StudyFileID', 'StudyFileId')
  setkey(workflow, StudyFileId)
  
  # Read in target PSMs and PCMs in order to have access to modified sequences and charges
  pdresult_t <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideGroupID, Sequence, ParentProteinAccessions, qValue, SVMScore, Abundances from TargetPeptideGroups;'))
  pdresult_mapping_t <- setDT(DBI::dbGetQuery(pd, 'SELECT TargetPsmsPeptideID, TargetPeptideGroupsPeptideGroupID From TargetPsmsTargetPeptideGroups;'))
  pdresult_psms_t <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideID, StudyFileId, Modifications, Charge, PrecursorAbundance From TargetPsms;'))
  setkey(pdresult_psms_t, PeptideID)
  setkey(pdresult_mapping_t, TargetPsmsPeptideID)
  pdresult_psms_t[pdresult_mapping_t,PeptideGroupID:=TargetPeptideGroupsPeptideGroupID]
  peptideGroup_quan <- pdresult_psms_t[,list(PeptideGroupID, StudyFileId, PrecursorAbundance)]
  peptideGroup_quan <- unique(peptideGroup_quan)
  peptideGroup_quan[,modPep_PrecursorAbundance:=sum(PrecursorAbundance, na.rm = T), by = list(PeptideGroupID, StudyFileId)]
  peptideGroup_quan[modPep_PrecursorAbundance == 0, modPep_PrecursorAbundance:=NA]
  peptideGroup_quan <- unique(peptideGroup_quan[,list(PeptideGroupID, StudyFileId, modPep_PrecursorAbundance)])
  pdresult_psms_t <- merge(pdresult_psms_t, 
                           peptideGroup_quan,
                           by = c('PeptideGroupID', 'StudyFileId'),
                           all.x = T)
  peptideGroup_abundances <- transpose(pdresult_t[, lapply(Abundances, convertBlobAsDoublePd)])
  setnames(peptideGroup_abundances, workflow$StudyFileId)
  peptideGroup_abundances <- cbind(pdresult_t[,list(PeptideGroupID)],
                                   peptideGroup_abundances)
  peptideGroup_abundances <- melt(peptideGroup_abundances, id.vars = 'PeptideGroupID')
  setnames(peptideGroup_abundances, c('variable', 'value'), c('StudyFileId', 'Abundances'))
  pdresult_t <- merge(pdresult_t[,Abundances:=NULL], 
                      peptideGroup_abundances, 
                      by = 'PeptideGroupID',
                      all.x = T)
  pdresult_t <- merge(pdresult_t,
                      pdresult_psms_t[,PeptideID:=NULL], 
                      by = c('PeptideGroupID', 'StudyFileId'), 
                      all.x = T)
  pdresult_t <- pdresult_t[Abundances != 0|is.na(Abundances)]
  pdresult_t <- unique(pdresult_t)
  pdresult_t[,DECOY:=0]
  
  # Read in decoy PSMs and PCMs in order to have access to modified sequences and charges
  pdresult_d <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideGroupID, Sequence, ParentProteinAccessions, qValue, SVMScore from DecoyPeptideGroups;'))
  pdresult_d[,Abundances:=NA]
  pdresult_mapping_d <- setDT(DBI::dbGetQuery(pd, 'SELECT DecoyPsmsPeptideID, DecoyPeptideGroupsPeptideGroupID From DecoyPsmsDecoyPeptideGroups;'))
  pdresult_psms_d <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideID, StudyFileId, Modifications, Charge From DecoyPsms;'))
  pdresult_psms_d[,PrecursorAbundance:=NA]
  setkey(pdresult_psms_d, PeptideID)
  setkey(pdresult_mapping_d, DecoyPsmsPeptideID)
  pdresult_psms_d[pdresult_mapping_d,PeptideGroupID:=DecoyPeptideGroupsPeptideGroupID]
  peptideGroup_quan <- pdresult_psms_d[,list(PeptideGroupID, StudyFileId, PrecursorAbundance)]
  peptideGroup_quan <- unique(peptideGroup_quan)
  peptideGroup_quan[,modPep_PrecursorAbundance:=NA]
  peptideGroup_quan <- unique(peptideGroup_quan[,list(PeptideGroupID, StudyFileId, modPep_PrecursorAbundance)])
  pdresult_psms_d <- merge(pdresult_psms_d, 
                           peptideGroup_quan,
                           by = c('PeptideGroupID', 'StudyFileId'),
                           all.x = T)
  pdresult_d <- merge(pdresult_d,
                      pdresult_psms_d[,PeptideID:=NULL], 
                      by = c('PeptideGroupID'), 
                      all.x = T)
  pdresult_d <- pdresult_d[Abundances != 0|is.na(Abundances)]
  pdresult_d <- unique(pdresult_d)
  pdresult_d[,PeptideGroupID:=paste0(PeptideGroupID, '_DECOY')]
  pdresult_d[,DECOY:=1]
  setcolorder(pdresult_d, colnames(pdresult_t))
  
  # Concatenate target PCMs and decoy PCMs
  pdresult <- rbindlist(list(pdresult_t, pdresult_d))
  pdresult[,PrecursorAbundance:=NULL]
  pdresult[,modPep_PrecursorAbundance:=NULL]
  
  # Sample name mapping
  setkey(pdresult, StudyFileId)
  pdresult <- workflow[pdresult]
  
  # Modifications
  pdresult[Modifications%like%'M[0-9]{1,2}',Mpos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Oxidation\\)|M(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'M[0-9]{1,2}',Mpos:=NA]
  pdresult[Modifications%like%'S[0-9]{1,2}',Spos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|S(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'S[0-9]{1,2}',Spos:=NA]
  pdresult[Modifications%like%'T[0-9]{1,2}',Tpos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|T(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'T[0-9]{1,2}',Tpos:=NA]
  pdresult[Modifications%like%'Y[0-9]{1,2}',Ypos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|Y(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'Y[0-9]{1,2}',Ypos:=NA]
  pdresult[,ID:=1:nrow(pdresult)]
  pdresult[!is.na(Spos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 's'),by=ID]
  pdresult[!is.na(Tpos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 't'),by=ID]
  pdresult[!is.na(Ypos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 'y'),by=ID]
  pdresult[!is.na(Mpos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Mpos, 'm'),by=ID]
  pdresult[is.na(MODIFIED_SEQUENCE),MODIFIED_SEQUENCE:=Sequence]
  pdresult[,MODIFIED_SEQUENCE:=gsub('C', 'c', MODIFIED_SEQUENCE)]
  
  # ORGANISM mapping
  proteinsFromFasta[grepl('Random_[0-9]', FastaTitleLines),Accession:=FastaTitleLines]
  proteinsFromFasta[!grepl('Random_[0-9]', Accession),Accession:=sapply(FastaTitleLines, function(x) str_split(x, "\\|")[[1]][2])]
  proteinsFromFasta <- rbindlist(list(proteinsFromFasta,
                                      data.table(Accession='mimic',FastaTitleLines='>sp|mimic|mimic_ENTRAPMENT', 
                                                 ORGANISM='ENTRAPMENT', ORGANISM_J='ENTRAPMENT')))
  setkey(proteinsFromFasta, Accession)
  setkey(pdresult, ParentProteinAccessions)
  pdresult[proteinsFromFasta,ORGANISM:=ORGANISM]
  pdresult[grepl(';', ParentProteinAccessions),ORGANISM:=unlist(pbmclapply(strsplit(ParentProteinAccessions, '; '), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM))], collapse = ';'), mc.cores = 10))]
  
  # ENTRAPMENT
  pdresult[,ENTRAPMENT:=ifelse(ORGANISM=='ENTRAPMENT', 1, 0)]
  
  if(!all(pdresult[ENTRAPMENT == 1,ORGANISM] == 'ENTRAPMENT')){
    print("No entrapments found.")
  }
  
  # DECOYS
  pdresult[DECOY == 1,ORGANISM:='DECOY']
  
  # FORMAT ENTRAPMENT PROTEIN ACCESSIONS
  pdresult[,ParentProteinAccessions:=gsub(' ', '', ParentProteinAccessions)]
  pdresult[,ParentProteinAccessions:=gsub('mimic', 'Random_', ParentProteinAccessions)]
  
  pdresult[,ParentProteinAccessions:=gsub(' ', '', ParentProteinAccessions)]
  pdresult[,ParentProteinAccessions:=gsub('mimic', 'Random_', ParentProteinAccessions)]
  
  # PEPTIDE GROUP IDS
  pdresult[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  pdresult[is.na(m_COUNT),m_COUNT:=0]
  pdresult[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  pdresult[is.na(p_COUNT),p_COUNT:=0]
  pdresult[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  pdresult[is.na(c_COUNT),c_COUNT:=0]
  pdresult[,MODPEP_ID:=paste0(Sequence, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  pdresult[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  
  pdresult[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', Charge)]
  pdresult[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  
  pdresult[,c('ID', 'Modifications', 'Mpos', 'Spos', 'Tpos', 'Ypos', 'm_COUNT', 'c_COUNT', 'p_COUNT'):=NULL]
  setnames(pdresult, c('FileName', 'PeptideGroupID', 'Sequence', 'ParentProteinAccessions', 'qValue', 'Charge', 'Abundances', 'SVMScore'),
           c('SAMPLE', 'ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION', 'Q_VALUE', 'CHARGE', 'QUAN', 'SVMSCORE'))
  pdresult[,SOFTWARE:='CHIMERYS']
  pdresult <- unique(pdresult)
  
  # Access spectral angles from percolator input 
  percolatorInput <- fread(pathToPercolatorInput, stringsAsFactors = F, integer64 = 'double')
  percolatorInput <- percolatorInput[,list(Peptide, MS8_feature_99)]
  percolatorInput[,MODIFIED_PEPTIDE_ID:=gsub('_', '', Peptide, fixed = T)]
  percolatorInput[,MODIFIED_PEPTIDE_ID:=gsub('.', '', MODIFIED_PEPTIDE_ID, fixed = T)]
  percolatorInput[,MODIFIED_PEPTIDE_ID:=substr(MODIFIED_PEPTIDE_ID, 1, nchar(MODIFIED_PEPTIDE_ID)-2)]
  percolatorInput[,MODIFIED_PEPTIDE_ID:=as.numeric(MODIFIED_PEPTIDE_ID)]
  percolatorInput <- setDT(percolatorInput)[, .SD[which.max(MS8_feature_99)], by=MODIFIED_PEPTIDE_ID]
  
  # Add quantification data from .ms8 
  ms8_samples <- setDT(DBI::dbReadTable(pd, 'MSnSpectrumInfo'))
  ms8_samples <- ms8_samples[,list(SpectrumFileID, StudyFileId)]
  ms8_samples <- unique(ms8_samples)
  # ms8_samples[,RAWFILE_ID:=SpectrumFileID - min(SpectrumFileID)]
  ms8_samples[,RAWFILE_ID:=seq(0, nrow(ms8_samples)-1, 1)]
  
  ms8_psms <- setDT(DBI::dbGetQuery(ms8, 'SELECT CANDIDATE_ID, MODIFIED_PEPTIDE_ID, MODIFIED_SEQUENCE, PRECURSOR_CHARGE, RAWFILE_ID, DECOY, IS_IDENTIFIED_BY_MBR, MS8_feature_7, MS8_feature_20 from RESULT_MAPPING_THERMO_VIEW;'))
  ms8_psms <- ms8_psms[DECOY == 0]
  ms8_psms <- unique(ms8_psms)
  
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('C[UNIMOD:4]', 'c', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('M[UNIMOD:35]', 'm', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('S[UNIMOD:21]', 's', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('T[UNIMOD:21]', 't', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('Y[UNIMOD:21]', 'y', MODIFIED_SEQUENCE, fixed = TRUE)]
  
  ms8_psms[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  ms8_psms[is.na(m_COUNT),m_COUNT:=0]
  ms8_psms[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  ms8_psms[is.na(p_COUNT),p_COUNT:=0]
  ms8_psms[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  ms8_psms[is.na(c_COUNT),c_COUNT:=0]
  ms8_psms[,MODPEP_ID:=paste0(toupper(MODIFIED_SEQUENCE), '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  ms8_psms[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  ms8_psms[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', PRECURSOR_CHARGE)]
  ms8_psms[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  ms8_psms[,APEX:=MS8_feature_7 == max(MS8_feature_7, na.rm = T), by = .(MODPEP_ID)]
  
  ms8_psms <- merge(ms8_psms, 
                    percolatorInput, 
                    by = 'MODIFIED_PEPTIDE_ID', 
                    all.x = T)
  # any(is.na(ms8_psms$MS8_feature_99))
  ms8_psms <- merge(ms8_psms, 
                    ms8_samples[,list(RAWFILE_ID, StudyFileId)], 
                    by = 'RAWFILE_ID', all.x = TRUE)
  ms8_psms[,SPECTRUM_SIMILARITY:=max(MS8_feature_99, na.rm = T), by = .(MODPEP_ID, StudyFileId)]
  ms8_psms[,FRAGMENTS_USED:=max(MS8_feature_20, na.rm = T), by = .(MODPEP_ID, StudyFileId)]
  # any(is.na(ms8_psms$SPECTRUM_SIMILARITY))
  # any(is.na(ms8_psms$FRAGMENTS_USED))
  
  ms8_merge_psms <- ms8_psms[,list(MODPEP_ID, 
                                   StudyFileId, 
                                   IS_IDENTIFIED_BY_MBR,
                                   SPECTRUM_SIMILARITY,
                                   FRAGMENTS_USED)]
  ms8_merge_psms <- unique(ms8_merge_psms)
  pdresult <- merge(pdresult,
                    ms8_merge_psms, 
                    by = c('MODPEP_ID', 'StudyFileId'),
                    all.x = T)
  # any(is.na(pdresult[DECOY != 1, SPECTRUM_SIMILARITY]))
  # any(is.na(pdresult[DECOY != 1, FRAGMENTS_USED]))
  
  pdresult[,'StudyFileId':=NULL]
  setcolorder(pdresult, neworder = c('ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION', 'Q_VALUE', 'SVMSCORE', 'DECOY', 'CHARGE', 'MODIFIED_SEQUENCE', 'QUAN', 'IS_IDENTIFIED_BY_MBR', 'SAMPLE', 'FRAGMENTS_USED', 'SPECTRUM_SIMILARITY', 'ORGANISM', 'ENTRAPMENT', 'PCM_ID', 'PCM_J_ID', 'MODPEP_ID', 'MODPEP_J_ID', 'SOFTWARE'))      
  return(pdresult)
  
  dbDisconnect(pd)
  dbDisconnect(ms8)
}

interpolateQvalues <- function(diann, rule = 1) {
  # Q.Value interpolation for decoys
  # Remove all Q_VALUEs for Decoys
  diann[DECOY==1,Q_VALUE:=NA]
  # Set Q_VALUE for DECOYS which get out-of-distribution SVMSCOREs to 1
  diann[DECOY==1 & SVMSCORE==(-10000000.000000),Q_VALUE:=1]
  # Calculate minimum Q_VALUE per unique SAMPLE,SVMSCORE combination
  qval <- diann[!is.na(Q_VALUE) & SVMSCORE!=(-10000000.000000),list(MIN_Q_VALUE=min(Q_VALUE)),by=list(SAMPLE, SVMSCORE)]
  # Find entries without Q_VALUE
  qnas <- diann[is.na(Q_VALUE)]
  # Sort data.tables by SVMSCORE for subsequent interpolation
  setkey(qval, SVMSCORE)
  setkey(qnas, SVMSCORE)
  # Perform SAMPLE-wise interpolation
  sapply(qnas[,unique(SAMPLE)], function(i) qnas[SAMPLE==i,Q_VALUE:=approx(x = qval[SAMPLE==i,SVMSCORE], y = qval[SAMPLE==i,MIN_Q_VALUE], xout = SVMSCORE, ties = 'ordered', method = 'linear', rule = rule)$y])
  # Concatenate results
  diann_inter <- rbindlist(list(diann[!is.na(Q_VALUE)], qnas))
  return(diann_inter)
}

readPdResult_localPcmGrouper <- function(pathToPdResult, pathToMs8, proteinsFromFasta, digest = NULL){
  cat('Reading in data from pdresult\n')
  
  pd <- DBI::dbConnect(RSQLite::SQLite(), pathToPdResult)
  ms8 <- DBI::dbConnect(RSQLite::SQLite(), pathToMs8)
  
  # Read in workflow to get rawfile names
  workflow <- setDT(DBI::dbGetQuery(pd, 'SELECT StudyFileID, FileName From WorkflowInputFiles;'))[!is.na(StudyFileID)]
  workflow[,FileName:=gsub('\\..*', '', sapply(strsplit(FileName, '\\\\'), 'tail', 1))]
  setnames(workflow, 'StudyFileID', 'StudyFileId')
  setkey(workflow, StudyFileId)
  
  # Read in target PSMs and PCMs in order to have access to modified sequences and charges
  pdresult_t <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideGroupID, Sequence, ParentProteinAccessions, qValue, SVMScore from TargetPeptideGroups;'))
  pdresult_mapping_t <- setDT(DBI::dbGetQuery(pd, 'SELECT TargetPsmsPeptideID, TargetPeptideGroupsPeptideGroupID From TargetPsmsTargetPeptideGroups;'))
  pdresult_psms_t <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideID, StudyFileId, Modifications, Charge, QuanValue From TargetPsms;'))
  setkey(pdresult_psms_t, PeptideID)
  setkey(pdresult_mapping_t, TargetPsmsPeptideID)
  pdresult_psms_t[pdresult_mapping_t,PeptideGroupID:=TargetPeptideGroupsPeptideGroupID]
  pdresult_t <- merge(pdresult_t, 
                      pdresult_psms_t[,list(PeptideGroupID, StudyFileId, Modifications, Charge, QuanValue)], 
                      by = 'PeptideGroupID', 
                      all.x = T)
  pdresult_t <- unique(pdresult_t)
  pdresult_t[,DECOY:=0]
  
  # Read in decoy PSMs and PCMs in order to have access to modified sequences and charges
  pdresult_d <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideGroupID, Sequence, ParentProteinAccessions, qValue, SVMScore from DecoyPeptideGroups;'))
  pdresult_mapping_d <- setDT(DBI::dbGetQuery(pd, 'SELECT DecoyPsmsPeptideID, DecoyPeptideGroupsPeptideGroupID From DecoyPsmsDecoyPeptideGroups;'))
  pdresult_psms_d <- setDT(DBI::dbGetQuery(pd, 'SELECT PeptideID, StudyFileId, Modifications, Charge From DecoyPsms;'))
  pdresult_psms_d[,QuanValue:=NA]
  setkey(pdresult_psms_d, PeptideID)
  setkey(pdresult_mapping_d, DecoyPsmsPeptideID)
  pdresult_psms_d[pdresult_mapping_d,PeptideGroupID:=DecoyPeptideGroupsPeptideGroupID]
  pdresult_d <- merge(pdresult_d, 
                      pdresult_psms_d[,list(PeptideGroupID, StudyFileId, Modifications, Charge, QuanValue)], 
                      by = 'PeptideGroupID', 
                      all.x = T)
  pdresult_d <- unique(pdresult_d)
  pdresult_d[,PeptideGroupID:=paste0(PeptideGroupID, '_DECOY')]
  pdresult_d[,DECOY:=1]
  
  # Concatenate target PCMs and decoy PCMs
  pdresult <- rbindlist(list(pdresult_t, pdresult_d))
  
  # Sample name mapping
  setkey(pdresult, StudyFileId)
  pdresult <- workflow[pdresult]
  
  # Modifications
  pdresult[Modifications%like%'M[0-9]{1,2}',Mpos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Oxidation\\)|M(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'M[0-9]{1,2}',Mpos:=NA]
  pdresult[Modifications%like%'S[0-9]{1,2}',Spos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|S(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'S[0-9]{1,2}',Spos:=NA]
  pdresult[Modifications%like%'T[0-9]{1,2}',Tpos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|Y[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|T(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'T[0-9]{1,2}',Tpos:=NA]
  pdresult[Modifications%like%'Y[0-9]{1,2}',Ypos:=sapply(strsplit(gsub('C[0-9]{1,2}\\(Carbamidomethyl\\)(; )*|M[0-9]{1,2}\\(Oxidation\\)(; )*|S[0-9]{1,2}\\(Phospho\\)(; )*|T[0-9]{1,2}\\(Phospho\\)(; )*|\\(Phospho\\)|Y(?=[0-9]{1,2})', '', Modifications, perl = T), '; '), as.integer)]
  pdresult[!Modifications%like%'Y[0-9]{1,2}',Ypos:=NA]
  pdresult[,ID:=1:nrow(pdresult)]
  pdresult[!is.na(Spos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 's'),by=ID]
  pdresult[!is.na(Tpos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 't'),by=ID]
  pdresult[!is.na(Ypos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Spos, 'y'),by=ID]
  pdresult[!is.na(Mpos),MODIFIED_SEQUENCE:=replaceChar(Sequence, Mpos, 'm'),by=ID]
  pdresult[is.na(MODIFIED_SEQUENCE),MODIFIED_SEQUENCE:=Sequence]
  pdresult[,MODIFIED_SEQUENCE:=gsub('C', 'c', MODIFIED_SEQUENCE)]
  
  if (!is.null(digest)) {
    warning('proteinsFromFasta ignored! Mapping ORGANISM based on digest instead!')
    setkey(digest, Sequence)
    setkey(pdresult, Sequence)
    pdresult[digest,ORGANISM:=i.ORGANISM]
  } else {
    # ORGANISM mapping
    proteinsFromFasta[grepl('Random_[0-9]', FastaTitleLines),Accession:=FastaTitleLines]
    proteinsFromFasta[!grepl('Random_[0-9]', Accession),Accession:=sapply(FastaTitleLines, function(x) str_split(x, "\\|")[[1]][2])]
    proteinsFromFasta <- rbindlist(list(proteinsFromFasta,
                                        data.table(Accession='mimic',FastaTitleLines='>sp|mimic|mimic_ENTRAPMENT', 
                                                   ORGANISM='ENTRAPMENT', ORGANISM_J='ENTRAPMENT')))
    setkey(proteinsFromFasta, Accession)
    setkey(pdresult, ParentProteinAccessions)
    pdresult[proteinsFromFasta,ORGANISM:=ORGANISM]
    pdresult[grepl(';', ParentProteinAccessions),ORGANISM:=unlist(pbmclapply(strsplit(ParentProteinAccessions, '; '), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM))], collapse = ';'), mc.cores = 10))]
  }
  
  # ENTRAPMENT
  pdresult[,ENTRAPMENT:=ifelse(ORGANISM=='ENTRAPMENT', 1, 0)]
  
  if(!all(pdresult[ENTRAPMENT == 1,ORGANISM] == 'ENTRAPMENT')){
    print("No entrapments found.")
  }
  
  # DECOYS
  pdresult[DECOY == 1,ORGANISM:='DECOY']
  
  pdresult[,ParentProteinAccessions:=gsub(' ', '', ParentProteinAccessions)]
  pdresult[,ParentProteinAccessions:=gsub('mimic', 'Random_', ParentProteinAccessions)]
  
  pdresult[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  pdresult[is.na(m_COUNT),m_COUNT:=0]
  pdresult[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  pdresult[is.na(p_COUNT),p_COUNT:=0]
  pdresult[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  pdresult[is.na(c_COUNT),c_COUNT:=0]
  pdresult[,MODPEP_ID:=paste0(Sequence, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  pdresult[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  
  pdresult[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', Charge)]
  pdresult[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  
  pdresult[,c('ID', 'Modifications', 'Mpos', 'Spos', 'Tpos', 'Ypos', 'm_COUNT', 'c_COUNT', 'p_COUNT'):=NULL]
  setnames(pdresult, c('FileName', 'PeptideGroupID', 'Sequence', 'ParentProteinAccessions', 'qValue', 'Charge', 'QuanValue', 'SVMScore'),
           c('SAMPLE', 'ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION', 'Q_VALUE', 'CHARGE', 'QUAN', 'SVMSCORE'))
  pdresult[,SOFTWARE:='CHIMERYS']
  # pdresult[,ID:=as.character(ID*10+DECOY)]
  pdresult <- unique(pdresult)
  
  # Add quantification data from .ms8 
  ms8_samples <- setDT(DBI::dbReadTable(pd, 'MSnSpectrumInfo'))
  ms8_samples <- ms8_samples[,list(SpectrumFileID, StudyFileId)]
  ms8_samples <- unique(ms8_samples)
  # ms8_samples[,RAWFILE_ID:=SpectrumFileID - min(SpectrumFileID)]
  ms8_samples[,RAWFILE_ID:=seq(0, nrow(ms8_samples)-1, 1)]
  
  ms8_pcm_quan <- setDT(DBI::dbReadTable(ms8, 'PCM_QUANTIFICATION'))
  ms8_quan_intensities <- lapply(ms8_pcm_quan$INTENSITY_BLOB, function(x) convertIntensityBlobs(x))
  ms8_quan_intensities <- lapply(ms8_quan_intensities, unlist)
  names(ms8_quan_intensities) <- unique(ms8_pcm_quan$PCM_ID)
  
  ms8_apex_intensities <- lapply(ms8_quan_intensities, max)
  ms8_apex_intensities <- lapply(ms8_apex_intensities, as.data.table)
  ms8_apex_intensities <- rbindlist(ms8_apex_intensities, idcol = 'name')
  setnames(ms8_apex_intensities, c('PCM_ID', 'APEX_MS8_feature_7'))
  ms8_apex_intensities[,PCM_ID:=bit64::as.integer64.character(PCM_ID)]
  
  ms8_dppp <- lapply(ms8_quan_intensities, function(x) x[x != 0])
  ms8_dppp <- lapply(ms8_dppp, length)
  ms8_dppp <- lapply(ms8_dppp, as.data.table)
  ms8_dppp <- rbindlist(ms8_dppp, idcol = 'name')
  setnames(ms8_dppp, c('PCM_ID', 'DPPP'))
  ms8_dppp[,PCM_ID:=bit64::as.integer64.character(PCM_ID)]
  
  ms8_quan_cols <- c(colnames(ms8_pcm_quan), 'APEX_MS8_feature_7')
  setkey(ms8_pcm_quan, PCM_ID)
  setkey(ms8_apex_intensities, PCM_ID)
  ms8_pcm_quan <- ms8_apex_intensities[ms8_pcm_quan]
  setcolorder(ms8_pcm_quan, ms8_quan_cols)
  
  ms8_quan_cols <- c(colnames(ms8_pcm_quan), 'DPPP')
  setkey(ms8_pcm_quan, PCM_ID)
  setkey(ms8_dppp, PCM_ID)
  ms8_pcm_quan <- ms8_dppp[ms8_pcm_quan]
  setcolorder(ms8_pcm_quan, ms8_quan_cols)
  
  ms8_pcm_cand <- setDT(DBI::dbGetQuery(ms8, 'SELECT CANDIDATE_ID, PCM_ID From CANDIDATE_TO_PCM;'))
  setkey(ms8_pcm_cand, PCM_ID)
  ms8_pcm_quan <- ms8_pcm_cand[ms8_pcm_quan]
  setcolorder(ms8_pcm_quan, c('CANDIDATE_ID', ms8_quan_cols))
  ms8_pcm_quan[,ONTOLOGY_TERM_ID:=NULL]
  ms8_pcm_quan[,RETENTION_TIME_BLOB:=NULL]
  ms8_pcm_quan[,INTENSITY_BLOB:=NULL]
  
  ms8_psms <- setDT(DBI::dbGetQuery(ms8, 'SELECT CANDIDATE_ID, MODIFIED_SEQUENCE, PRECURSOR_CHARGE, RAWFILE_ID, DECOY, IS_IDENTIFIED_BY_MBR, MS8_feature_7, MS8_feature_20 from RESULT_MAPPING_THERMO_VIEW;'))
  ms8_psms <- ms8_psms[DECOY == 0]
  ms8_psms <- merge(ms8_psms, ms8_pcm_quan, by = 'CANDIDATE_ID', all.x = T)
  setnames(ms8_psms, 'PCM_ID', 'MS8_PCM_ID')
  ms8_psms <- unique(ms8_psms)
  ms8_psms[,APEX:=MS8_feature_7 == max(MS8_feature_7, na.rm = T), by = .(MODIFIED_SEQUENCE, PRECURSOR_CHARGE, RAWFILE_ID)]
  # ms8_psms[, APEX:=as.integer(MS8_feature_7) == as.integer(APEX_MS8_feature_7)]
  # ms8_psms[APEX == FALSE, APEX:=MS8_feature_7 == APEX_MS8_feature_7]
  # ms8_psms[APEX == FALSE, APEX:=round(MS8_feature_7,0) == round(APEX_MS8_feature_7,0)]
  
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('C[UNIMOD:4]', 'c', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('M[UNIMOD:35]', 'm', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('S[UNIMOD:21]', 's', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('T[UNIMOD:21]', 't', MODIFIED_SEQUENCE, fixed = TRUE)]
  ms8_psms[,MODIFIED_SEQUENCE:=gsub('Y[UNIMOD:21]', 'y', MODIFIED_SEQUENCE, fixed = TRUE)]
  
  ms8_psms[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  ms8_psms[is.na(m_COUNT),m_COUNT:=0]
  ms8_psms[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  ms8_psms[is.na(p_COUNT),p_COUNT:=0]
  ms8_psms[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  ms8_psms[is.na(c_COUNT),c_COUNT:=0]
  ms8_psms[,MODPEP_ID:=paste0(toupper(MODIFIED_SEQUENCE), '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  ms8_psms[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  ms8_psms[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', PRECURSOR_CHARGE)]
  ms8_psms[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  
  ms8_psms <- merge(ms8_psms, 
                    ms8_samples[,list(RAWFILE_ID, StudyFileId)], 
                    by = 'RAWFILE_ID', all.x = TRUE)
  
  ms8_apex_psms <- ms8_psms[APEX == TRUE]
  ms8_apex_psms <- ms8_apex_psms[,list(PCM_ID, StudyFileId, 
                                       IS_IDENTIFIED_BY_MBR,
                                       QUANTIFICATION_VALUE,
                                       MS8_feature_20,
                                       DPPP,
                                       APEX_MS8_feature_7)]
  ms8_apex_psms <- unique(ms8_apex_psms)
  pdresult <- merge(pdresult,
                    ms8_apex_psms, 
                    by = c('PCM_ID', 'StudyFileId'),
                    all.x = T)
  
  ms8_psms <- ms8_psms[,list(MODPEP_J_ID, PCM_J_ID, QUANTIFICATION_VALUE, StudyFileId)]
  ms8_psms <- ms8_psms[!is.na(QUANTIFICATION_VALUE)]
  ms8_psms <- unique(ms8_psms)
  ms8_psms[,LOCAL_MODPEP_J_ID:=paste0(MODPEP_J_ID, '_', StudyFileId)]
  ms8_psms[,LOCAL_PCM_J_ID:=paste0(PCM_J_ID, '_', StudyFileId)]
  
  pdresult[,LOCAL_MODPEP_J_ID:=paste0(MODPEP_J_ID, '_', StudyFileId)]
  pdresult[,LOCAL_PCM_J_ID:=paste0(PCM_J_ID, '_', StudyFileId)]
  
  ms8_mbr <- ms8_psms[PCM_J_ID %in% pdresult[Q_VALUE <= 0.01 &
                                               DECOY == 0, PCM_J_ID] & 
                        !LOCAL_PCM_J_ID %in% pdresult[!is.na(QUAN), LOCAL_PCM_J_ID],
                      list(PCM_J_ID, StudyFileId, QUANTIFICATION_VALUE, LOCAL_PCM_J_ID)]
  pd_mbr <- pdresult[LOCAL_PCM_J_ID %in% ms8_mbr$LOCAL_PCM_J_ID]
  
  if(!identical(nrow(pd_mbr), nrow(ms8_mbr))){
    cat('Quantification in main_search.ms8 and .pdresult not identical')
    break
  }else{
    pd_mbr[,QUAN:=QUANTIFICATION_VALUE]
    setcolorder(pd_mbr, colnames(pdresult))
    local_pcm_j_id_check <- nrow(pdresult)
    pdresult <- pdresult[!LOCAL_PCM_J_ID %in% pd_mbr$LOCAL_PCM_J_ID]
    pdresult <- rbind(pdresult, 
                      pd_mbr, 
                      fill = T) 
    if(nrow(pdresult) != local_pcm_j_id_check){
      cat('Quantification in main_search.ms8 and .pdresult not identical')
      break
    }else{
      setnames(pdresult, 'MS8_feature_20', 'FRAGMENTS_USED')
      pdresult[,SPECTRUM_SIMILARITY:=NA]
      pdresult[,c('StudyFileId', 'LOCAL_MODPEP_J_ID', 'LOCAL_PCM_J_ID'):=NULL]
      setcolorder(pdresult, neworder = c('ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION', 'Q_VALUE', 'SVMSCORE', 'DECOY', 'CHARGE', 'MODIFIED_SEQUENCE', 'QUAN', 'QUANTIFICATION_VALUE', 'IS_IDENTIFIED_BY_MBR', 'SAMPLE', 'FRAGMENTS_USED', 'SPECTRUM_SIMILARITY', 'DPPP', 'APEX_MS8_feature_7', 'ORGANISM', 'ENTRAPMENT', 'PCM_ID', 'PCM_J_ID', 'MODPEP_ID', 'MODPEP_J_ID', 'SOFTWARE'))      
      return(pdresult)
    }
  }
  dbDisconnect(pd)
  dbDisconnect(ms8)
}

readDiann <- function(pathToTsv, proteinsFromFasta, digest = NULL) {
  cat('Reading in data from dia-nn tsv\n')
  diann <- fread(pathToTsv, stringsAsFactors = F, integer64 = 'double')
  files <- diann[,list(File.Name=unique(File.Name))]
  files[,newname:=gsub('\\..*', '', sapply(strsplit(File.Name, '\\\\'), 'tail', 1))]
  setkey(files, File.Name)
  setkey(diann, File.Name)
  diann[files,File.Name:=newname]
  diann <- rbindlist(list(diann[,list(ID=Precursor.Id, PEPTIDE_SEQUENCE=Stripped.Sequence, PROTEIN_ACCESSION=Protein.Ids, LIB_Q_VALUE=Lib.Q.Value, Q_VALUE=Q.Value, SVMSCORE=CScore, DECOY=0, SAMPLE=File.Name, CHARGE=Precursor.Charge, MODIFIED_SEQUENCE=Modified.Sequence, QUAN=Precursor.Quantity, FRAGMENTS_USED=Fragment.Quant.Raw, SPECTRUM_SIMILARITY=Spectrum.Similarity, DIANN_ORGANISM=Protein.Names, SOFTWARE='DIA-NN')],
                          diann[,list(ID=Precursor.Id, PEPTIDE_SEQUENCE=Stripped.Sequence, PROTEIN_ACCESSION=Protein.Ids, LIB_Q_VALUE=Lib.Q.Value, Q_VALUE=Q.Value, SVMSCORE=Decoy.CScore, DECOY=1, SAMPLE=File.Name, CHARGE=Precursor.Charge, MODIFIED_SEQUENCE=Modified.Sequence, QUAN=Precursor.Quantity, FRAGMENTS_USED=Fragment.Quant.Raw, SPECTRUM_SIMILARITY=Spectrum.Similarity, DIANN_ORGANISM=Protein.Names, SOFTWARE='DIA-NN')]))
  diann <- interpolateQvalues(diann)
  diann[,MODIFIED_SEQUENCE:=gsub('C\\(UniMod\\:4\\)', 'c', MODIFIED_SEQUENCE)]
  diann[,MODIFIED_SEQUENCE:=gsub('M\\(UniMod\\:35\\)', 'm', MODIFIED_SEQUENCE)]
  diann[,FRAGMENTS_USED:=sapply(strsplit(FRAGMENTS_USED, ';'), function(x) length(x[as.numeric(x) > 0]))]
  
  diann[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  diann[is.na(m_COUNT),m_COUNT:=0]
  diann[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  diann[is.na(p_COUNT),p_COUNT:=0]
  diann[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  diann[is.na(c_COUNT),c_COUNT:=0]
  diann[,MODPEP_ID:=paste0(PEPTIDE_SEQUENCE, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  diann[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  
  diann[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', CHARGE)]
  diann[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  
  diann[, c('m_COUNT', 'c_COUNT', 'p_COUNT'):=NULL] 
  
  # ORGANISM mapping
  ## via DIA-NN protein accessions
  diann[,DIANN_ORGANISM:=unlist(pbmclapply(sapply(strsplit(gsub(pattern = '(?<=^|;)[A-Za-z0-9]*_', '', DIANN_ORGANISM, perl = T), ';'), function(j) sort(unique(j))), function(i) paste0(i, collapse = ';'), mc.cores = 10))]
  
  if (!is.null(digest)) {
    ## via digest
    warning('proteinsFromFasta ignored! Mapping ORGANISM based on digest instead!')
    setkey(digest, Sequence)
    setkey(diann, PEPTIDE_SEQUENCE)
    diann[digest,c('ORGANISM', 'ORGANISM_FASTA_J'):=.(ORGANISM, ORGANISM_J)]
  } else {
    ## via fasta
    proteinsFromFasta[grepl('Random_[0-9]', FastaTitleLines),Accession:=FastaTitleLines]
    proteinsFromFasta[!grepl('Random_[0-9]', Accession),Accession:=sapply(FastaTitleLines, function(x) str_split(x, "\\|")[[1]][2])]
    setkey(proteinsFromFasta, Accession)
    setkey(diann, PROTEIN_ACCESSION)
    
    diann[proteinsFromFasta,c('ORGANISM', 'ORGANISM_FASTA_J'):=.(ORGANISM, ORGANISM_J)]
    diann[grepl(';', PROTEIN_ACCESSION),ORGANISM:=unlist(pbmclapply(strsplit(PROTEIN_ACCESSION, ';'), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM))], collapse = ';'), mc.cores = 10))]
    diann[grepl(';', PROTEIN_ACCESSION),ORGANISM_FASTA_J:=unlist(pbmclapply(strsplit(PROTEIN_ACCESSION, ';'), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM_J))], collapse = ';'), mc.cores = 10))]
  }
  
  diann[,ORGANISM_J:=paste0(sort(unique(unlist(strsplit(ORGANISM, ';')))), collapse = ';'), by = PCM_J_ID]
  
  # ENTRAPMENT
  diann[,ENTRAPMENT:=ifelse(ORGANISM=='ENTRAPMENT', 1, 0)]
  diann[,ENTRAPMENT_J:=ifelse(ORGANISM_J=='ENTRAPMENT', 1, 0)]
  diann[,ENTRAPMENT_FASTA_J:=ifelse(ORGANISM_FASTA_J=='ENTRAPMENT', 1, 0)]
  
  # DECOYS
  diann[DECOY == 1,ORGANISM:='DECOY']
  diann[DECOY == 1,ORGANISM_J:='DECOY']
  diann[DECOY == 1,ORGANISM_FASTA_J:='DECOY']
  
  setcolorder(diann, neworder = c('ID', 'PEPTIDE_SEQUENCE',
                                  'PROTEIN_ACCESSION',
                                  'Q_VALUE', 'SVMSCORE', 'DECOY', 'CHARGE',
                                  'MODIFIED_SEQUENCE', 'QUAN',
                                  'SAMPLE', 'FRAGMENTS_USED',
                                  'SPECTRUM_SIMILARITY',
                                  'DIANN_ORGANISM',
                                  'ORGANISM', 'ORGANISM_J', 'ORGANISM_FASTA_J',
                                  'ENTRAPMENT', 'ENTRAPMENT_J', 'ENTRAPMENT_FASTA_J',
                                  'PCM_ID', 'PCM_J_ID', 'MODPEP_ID',
                                  'MODPEP_J_ID', 'SOFTWARE'))
  diann <- unique(diann)
  
  return(diann)
}

readSpectronaut <- function(pathToExport, proteinsFromFasta, digest = NULL){
  cat('Reading in data from Spectronaut export\n')
  sn <- fread(pathToExport, stringsAsFactors = F, integer64 = 'double')
  
  sn <- sn[,list(ID=as.character(EG.PrecursorId), 
                 PEPTIDE_SEQUENCE=PEP.StrippedSequence, 
                 PROTEIN_ACCESSION=PG.ProteinAccessions, 
                 Q_VALUE=EG.Qvalue, 
                 SVMSCORE=EG.Cscore, 
                 DECOY=as.integer(EG.IsDecoy), QUAN=`EG.TotalQuantity (Settings)`,
                 SAMPLE=R.FileName, CHARGE=FG.Charge,
                 MODIFIED_SEQUENCE=EG.ModifiedSequence,
                 FRAGMENTS_USED=FG.FragmentCount,
                 # SN_ORGANISM=PG.Organisms,
                 SOFTWARE='SPECTRONAUT')]
  sn[,MODIFIED_SEQUENCE:=gsub(pattern = 'M\\[Oxidation \\(M\\)\\]', replacement = 'm', x = MODIFIED_SEQUENCE)]
  sn[,MODIFIED_SEQUENCE:=gsub(pattern = 'C\\[Carbamidomethyl \\(C\\)\\]', replacement = 'c', x = MODIFIED_SEQUENCE)]
  sn[,MODIFIED_SEQUENCE:=gsub(pattern = '_', replacement = '', x = MODIFIED_SEQUENCE)]
  
  sn[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  sn[is.na(m_COUNT),m_COUNT:=0]
  sn[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  sn[is.na(p_COUNT),p_COUNT:=0]
  sn[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  sn[is.na(c_COUNT),c_COUNT:=0]
  sn[,MODPEP_ID:=paste0(PEPTIDE_SEQUENCE, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  sn[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  
  sn[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', CHARGE)]
  sn[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  
  sn[, c('m_COUNT', 'c_COUNT', 'p_COUNT'):=NULL]
  
  # ORGANISM mapping
  ## via spectronaut protein accessions
  sn[,SN_ORGANISM:=NA]
  
  if (!is.null(digest)) {
    ## via digest
    warning('proteinsFromFasta ignored! Mapping ORGANISM based on digest instead!')
    setkey(digest, Sequence)
    setkey(sn, PEPTIDE_SEQUENCE)
    sn[digest,c('ORGANISM', 'ORGANISM_FASTA_J'):=.(ORGANISM, ORGANISM_J)]
    # TODO not sure why, but SN identifies this non-tryptic peptide, mapping to
    # Q9UPN3 (in output) and O94854 (in fasta, but not in output)
    sn[PEPTIDE_SEQUENCE=='ILNTVLS',ORGANISM:='HUMAN']
    sn[PEPTIDE_SEQUENCE=='ILNTVLS',ORGANISM_FASTA_J:='HUMAN']
  } else {
    ## via fasta
    proteinsFromFasta[grepl('Random_[0-9]', Accession),Accession:=sapply(FastaTitleLines, function(x) str_split(x, "\\|")[[1]][2])]
    setkey(proteinsFromFasta, Accession)
    setkey(sn, PROTEIN_ACCESSION)
    sn[proteinsFromFasta,c('ORGANISM', 'ORGANISM_FASTA_J'):=.(ORGANISM, ORGANISM_J)]
    sn[grepl(';', PROTEIN_ACCESSION),ORGANISM:=unlist(pbmclapply(strsplit(PROTEIN_ACCESSION, ';'), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM))], collapse = ';'), mc.cores = 10))]
    sn[grepl(';', PROTEIN_ACCESSION),ORGANISM_FASTA_J:=unlist(pbmclapply(strsplit(PROTEIN_ACCESSION, ';'), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM_J))], collapse = ';'), mc.cores = 10))]
    sn[ORGANISM=='',ORGANISM:='ENTRAPMENT']
    sn[ORGANISM_FASTA_J=='',ORGANISM_FASTA_J:='ENTRAPMENT']
  }
  
  sn[,ORGANISM_J:=paste0(sort(unique(unlist(strsplit(ORGANISM, ';')))), collapse = ';'), by = PCM_J_ID]
  
  # ENTRAPMENT
  sn[,ENTRAPMENT:=ifelse(ORGANISM=='ENTRAPMENT', 1, 0)]
  sn[,ENTRAPMENT_J:=ifelse(ORGANISM_J=='ENTRAPMENT', 1, 0)]
  sn[,ENTRAPMENT_FASTA_J:=ifelse(ORGANISM_FASTA_J=='ENTRAPMENT', 1, 0)]
  
  # DECOYS
  sn[DECOY == 1,ORGANISM:='DECOY']
  sn[DECOY == 1,ORGANISM_J:='DECOY']
  sn[DECOY == 1,ORGANISM_FASTA_J:='DECOY']
  
  sn[,SPECTRUM_SIMILARITY:=NA]
  setcolorder(sn, neworder = c('ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION',
                               'Q_VALUE', 'SVMSCORE', 'DECOY', 'CHARGE',
                               'MODIFIED_SEQUENCE',
                               'QUAN', 'SAMPLE',
                               'FRAGMENTS_USED', 'SPECTRUM_SIMILARITY',
                               'SN_ORGANISM', 
                               'ORGANISM', 'ORGANISM_J', 'ORGANISM_FASTA_J',
                               'ENTRAPMENT', 'ENTRAPMENT_J', 'ENTRAPMENT_FASTA_J',
                               'PCM_ID', 'PCM_J_ID',
                               'MODPEP_ID', 'MODPEP_J_ID',
                               'SOFTWARE'))
  sn <- unique(sn)
  return(sn)
}

readSpectronaut_flagImputation <- function(pathToExport, proteinsFromFasta, digest = NULL){
  cat('Reading in data from Spectronaut export\n')
  sn <- fread(pathToExport, stringsAsFactors = F, integer64 = 'double')
  
  sn[,SN_ID_FRAGS:=.N, by = .(EG.PrecursorId, R.FileName)]
  sn[,SN_QUAN_FRAGS:=.SD[F.ExcludedFromQuantification == FALSE, .N], by = .(EG.PrecursorId, R.FileName)]
  sn[,SN_QUAN_SAMPLES_ANY:=uniqueN(R.FileName[SN_QUAN_FRAGS >= 1]), by = .(EG.PrecursorId)]
  sn[,SN_QUAN_SAMPLES:=uniqueN(R.FileName[SN_QUAN_FRAGS >= 3]), by = .(EG.PrecursorId)]
  
  sn[,IMPUTED_MIN1:=ifelse(F.PeakArea <= 1, TRUE, FALSE)]
  sn[,MIN1_ID_FRAGS:=.SD[IMPUTED_MIN1 == FALSE, .N], by = .(EG.PrecursorId, R.FileName)]
  sn[,MIN1_QUAN_FRAGS:=.SD[IMPUTED_MIN1 == FALSE & F.ExcludedFromQuantification == FALSE, .N], by = .(EG.PrecursorId, R.FileName)]
  sn[,MIN1_QUAN_SAMPLES_ANY:=uniqueN(R.FileName[MIN1_QUAN_FRAGS >= 1]), by = .(EG.PrecursorId)]
  sn[,MIN1_QUAN_SAMPLES:=uniqueN(R.FileName[MIN1_QUAN_FRAGS >= 3]), by = .(EG.PrecursorId)]
  
  sn <- sn[,list(ID=as.character(EG.PrecursorId), 
                 PEPTIDE_SEQUENCE=PEP.StrippedSequence, 
                 PROTEIN_ACCESSION=PG.ProteinAccessions, 
                 Q_VALUE=EG.Qvalue, 
                 SVMSCORE=EG.Cscore, 
                 DECOY=as.integer(EG.IsDecoy), QUAN=`EG.TotalQuantity (Settings)`,
                 SAMPLE=R.FileName, CHARGE=FG.Charge,
                 MODIFIED_SEQUENCE=EG.ModifiedSequence,
                 FRAGMENTS_USED=FG.FragmentCount,
                 # SN_ORGANISM=PG.Organisms, 
                 SN_ID_FRAGS,
                 SN_QUAN_FRAGS,
                 SN_QUAN_SAMPLES_ANY,
                 SN_QUAN_SAMPLES, 
                 MIN1_ID_FRAGS,
                 MIN1_QUAN_FRAGS,
                 MIN1_QUAN_SAMPLES_ANY,
                 MIN1_QUAN_SAMPLES, 
                 SOFTWARE='SPECTRONAUT')]
  sn[,MODIFIED_SEQUENCE:=gsub(pattern = 'M\\[Oxidation \\(M\\)\\]', replacement = 'm', x = MODIFIED_SEQUENCE)]
  sn[,MODIFIED_SEQUENCE:=gsub(pattern = 'C\\[Carbamidomethyl \\(C\\)\\]', replacement = 'c', x = MODIFIED_SEQUENCE)]
  sn[,MODIFIED_SEQUENCE:=gsub(pattern = '_', replacement = '', x = MODIFIED_SEQUENCE)]
  
  sn[,m_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'm')]
  sn[is.na(m_COUNT),m_COUNT:=0]
  sn[,p_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 's|t|y')]
  sn[is.na(p_COUNT),p_COUNT:=0]
  sn[,c_COUNT:=stringr::str_count(MODIFIED_SEQUENCE, 'c')]
  sn[is.na(c_COUNT),c_COUNT:=0]
  sn[,MODPEP_ID:=paste0(PEPTIDE_SEQUENCE, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
  sn[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]
  
  sn[,PCM_ID:=paste0(MODIFIED_SEQUENCE, '_', CHARGE)]
  sn[,PCM_J_ID:=gsub('I|L', 'J', PCM_ID)]
  
  sn[, c('m_COUNT', 'c_COUNT', 'p_COUNT'):=NULL]
  
  # ORGANISM mapping
  ## via spectronaut protein accessions
  sn[,SN_ORGANISM:=NA]
  
  if (!is.null(digest)) {
    ## via digest
    warning('proteinsFromFasta ignored! Mapping ORGANISM based on digest instead!')
    setkey(digest, Sequence)
    setkey(sn, PEPTIDE_SEQUENCE)
    sn[digest,c('ORGANISM', 'ORGANISM_FASTA_J'):=.(ORGANISM, ORGANISM_J)]
    # TODO not sure why, but SN identifies this non-tryptic peptide, mapping to
    # Q9UPN3 (in output) and O94854 (in fasta, but not in output)
    sn[PEPTIDE_SEQUENCE=='ILNTVLS',ORGANISM:='HUMAN']
    sn[PEPTIDE_SEQUENCE=='ILNTVLS',ORGANISM_FASTA_J:='HUMAN']
  } else {
    ## via fasta
    proteinsFromFasta[grepl('Random_[0-9]', Accession),Accession:=sapply(FastaTitleLines, function(x) str_split(x, "\\|")[[1]][2])]
    setkey(proteinsFromFasta, Accession)
    setkey(sn, PROTEIN_ACCESSION)
    sn[proteinsFromFasta,c('ORGANISM', 'ORGANISM_FASTA_J'):=.(ORGANISM, ORGANISM_J)]
    sn[grepl(';', PROTEIN_ACCESSION),ORGANISM:=unlist(pbmclapply(strsplit(PROTEIN_ACCESSION, ';'), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM))], collapse = ';'), mc.cores = 10))]
    sn[grepl(';', PROTEIN_ACCESSION),ORGANISM_FASTA_J:=unlist(pbmclapply(strsplit(PROTEIN_ACCESSION, ';'), function(i) paste0(proteinsFromFasta[.(i),sort(unique(ORGANISM_J))], collapse = ';'), mc.cores = 10))]
    sn[ORGANISM=='',ORGANISM:='ENTRAPMENT']
    sn[ORGANISM_FASTA_J=='',ORGANISM_FASTA_J:='ENTRAPMENT']
  }
  
  sn[,ORGANISM_J:=paste0(sort(unique(unlist(strsplit(ORGANISM, ';')))), collapse = ';'), by = PCM_J_ID]
  
  # ENTRAPMENT
  sn[,ENTRAPMENT:=ifelse(ORGANISM=='ENTRAPMENT', 1, 0)]
  sn[,ENTRAPMENT_J:=ifelse(ORGANISM_J=='ENTRAPMENT', 1, 0)]
  sn[,ENTRAPMENT_FASTA_J:=ifelse(ORGANISM_FASTA_J=='ENTRAPMENT', 1, 0)]
  
  # DECOYS
  sn[DECOY == 1,ORGANISM:='DECOY']
  sn[DECOY == 1,ORGANISM_J:='DECOY']
  sn[DECOY == 1,ORGANISM_FASTA_J:='DECOY']
  
  sn[,SPECTRUM_SIMILARITY:=NA]
  setcolorder(sn, neworder = c('ID', 'PEPTIDE_SEQUENCE', 'PROTEIN_ACCESSION',
                               'Q_VALUE', 'SVMSCORE', 'DECOY', 'CHARGE',
                               'MODIFIED_SEQUENCE',
                               'SN_ID_FRAGS', 'SN_QUAN_FRAGS', 
                               'SN_QUAN_SAMPLES_ANY', 'SN_QUAN_SAMPLES', 
                               'MIN1_ID_FRAGS', 'MIN1_QUAN_FRAGS', 
                               'MIN1_QUAN_SAMPLES_ANY', 'MIN1_QUAN_SAMPLES', 
                               'QUAN', 'SAMPLE',
                               'FRAGMENTS_USED', 'SPECTRUM_SIMILARITY',
                               'SN_ORGANISM', 'ORGANISM', 'ORGANISM_J',
                               'ENTRAPMENT', 'ENTRAPMENT_J',
                               'PCM_ID', 'PCM_J_ID',
                               'MODPEP_ID', 'MODPEP_J_ID',
                               'SOFTWARE'))
  sn <- unique(sn)
  return(sn)
}

setEntrapmentFDR_local <- function(localPrecursors) {
  sapply(c('FDR', 'ENTRAPMENT_Q_VALUE', 'ENTRAPMENT_Q_VALUE_1', 'ENTRAPMENT_Q_VALUE_2'), function(i) {
    if(i%in%colnames(localPrecursors)) {
      localPrecursors[,c(i):=NULL]
    }
  })
  # localPrecursors <- copy(ratio1to8)
  setorder(localPrecursors,-SVMSCORE, na.last = F)
  uniqueSVM <- localPrecursors[,list(TARGET=sum(DECOY==0 & ENTRAPMENT_FASTA_J==0),
                                     DECOY=sum(DECOY==1 & ENTRAPMENT_FASTA_J==0),
                                     ENTRAPMENT_TARGET=sum(DECOY==0 & ENTRAPMENT_FASTA_J==1),
                                     ENTRAPMENT_DECOY=sum(DECOY==1 & ENTRAPMENT_FASTA_J==1)),by=list(SOFTWARE, SAMPLE, SVMSCORE)]
  uniqueSVM[,CUMSUM_TARGET:=cumsum(TARGET),by=list(SOFTWARE, SAMPLE)]
  uniqueSVM[,CUMSUM_DECOY:=cumsum(DECOY),by=list(SOFTWARE,SAMPLE)]
  uniqueSVM[,CUMSUM_ENTRAPMENT_TARGET:=cumsum(ENTRAPMENT_TARGET),by=list(SOFTWARE,SAMPLE)]
  uniqueSVM[,CUMSUM_ENTRAPMENT_DECOY:=cumsum(ENTRAPMENT_DECOY),by=list(SOFTWARE,SAMPLE)]
  uniqueSVM[,FDR:=(CUMSUM_DECOY + CUMSUM_ENTRAPMENT_DECOY) / (CUMSUM_TARGET + CUMSUM_ENTRAPMENT_TARGET),by=list(SOFTWARE, SAMPLE)]
  uniqueSVM[,ENTRAPMENT_Q_VALUE:=(CUMSUM_DECOY + CUMSUM_ENTRAPMENT_TARGET) / (CUMSUM_TARGET + CUMSUM_ENTRAPMENT_TARGET),by=list(SOFTWARE, SAMPLE)]
  uniqueSVM[,ENTRAPMENT_Q_VALUE_1:=(CUMSUM_ENTRAPMENT_TARGET) / (CUMSUM_TARGET + CUMSUM_ENTRAPMENT_TARGET),by=list(SOFTWARE, SAMPLE)]
  uniqueSVM[,ENTRAPMENT_Q_VALUE_2:=(CUMSUM_ENTRAPMENT_TARGET) / (CUMSUM_TARGET),by=list(SOFTWARE, SAMPLE)]
  setorder(uniqueSVM, SVMSCORE)
  for (col in c('FDR', 'ENTRAPMENT_Q_VALUE', 'ENTRAPMENT_Q_VALUE_1', 'ENTRAPMENT_Q_VALUE_2')) {
    uniqueSVM[,c(col):=cummin(get(col)),by=list(SOFTWARE, SAMPLE)]
  }
  uniqueSVM[,c('TARGET', 'DECOY', 'ENTRAPMENT_TARGET', 'ENTRAPMENT_DECOY'):=NULL]
  setkey(uniqueSVM, SOFTWARE, SAMPLE, SVMSCORE)
  setkey(localPrecursors, SOFTWARE, SAMPLE, SVMSCORE)
  localPrecursors[uniqueSVM,c('FDR', 'ENTRAPMENT_Q_VALUE', 'ENTRAPMENT_Q_VALUE_1', 'ENTRAPMENT_Q_VALUE_2'):=list(FDR, ENTRAPMENT_Q_VALUE, ENTRAPMENT_Q_VALUE_1, ENTRAPMENT_Q_VALUE_2)]
  return(localPrecursors)
}