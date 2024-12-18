# ---- setup ----
path <- file.path(here::here(), "figure-E3")
figurePath <- file.path(dataPath, "data/figure-E3")
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))


# ---- pathsToData ----
# pdresult
pathToPdResult <- file.path(dataPath, "LFQ_Bench_multispecies/DDA/Chimerys/20240517_lfq_dda_z1to4_True_noNormalization.pdResult")
#pathToPdResult <- '/mnt/paper/01_paper/figures/main/6_figure_DDA_vs_DIA/1_LFQ_Bench_DDA/20231107_lfq_dda_z1to4_True.pdResult'
pathToPdResult_pcms <- file.path(dataPath, "LFQ_Bench_multispecies/DDA/Chimerys/20240524_lfq_dda_z1to4_localPcmFdr_True.pdResult")
#pathToPdResult_pcms <- '/mnt/paper/01_paper/figures/main/6_figure_DDA_vs_DIA/1_LFQ_Bench_DDA/20240524_lfq_dda_z1to4_localPcmFdr_True.pdResult'

# read psms
pd <- dbConnect(RSQLite::SQLite(), pathToPdResult)
psms_t <- dbReadTable(pd, "TargetPsms") %>% as.data.table()
psms_d <- dbReadTable(pd, "DecoyPsms") %>% as.data.table()

# peptide groups
modPeps_t <- dbReadTable(pd, 'TargetPeptideGroups') %>% as.data.table()
modPeps_d <- dbReadTable(pd, 'DecoyPeptideGroups') %>% as.data.table()

# mapping
psms_modPeps_mapping_t <- setDT(DBI::dbGetQuery(pd, 'SELECT TargetPsmsPeptideID, TargetPeptideGroupsPeptideGroupID From TargetPsmsTargetPeptideGroups;'))
psms_modPeps_mapping_d <- setDT(DBI::dbGetQuery(pd, 'SELECT DecoyPsmsPeptideID, DecoyPeptideGroupsPeptideGroupID From DecoyPsmsDecoyPeptideGroups;'))
dbDisconnect(pd)

# pcms
# read corresponding psms
pd <- dbConnect(RSQLite::SQLite(), pathToPdResult_pcms)
pcmGrouper_psms_t <- dbReadTable(pd, "TargetPsms") %>% as.data.table()
pcmGrouper_psms_d <- dbReadTable(pd, "DecoyPsms") %>% as.data.table()

# pcms
pcms_t <- dbReadTable(pd, 'TargetPeptideGroups') %>% as.data.table()
pcms_d <- dbReadTable(pd, 'DecoyPeptideGroups') %>% as.data.table()

# mapping
psms_pcms_mapping_t <- setDT(DBI::dbGetQuery(pd, 'SELECT TargetPsmsPeptideID, TargetPeptideGroupsPeptideGroupID From TargetPsmsTargetPeptideGroups;'))
psms_pcms_mapping_d <- setDT(DBI::dbGetQuery(pd, 'SELECT DecoyPsmsPeptideID, DecoyPeptideGroupsPeptideGroupID From DecoyPsmsDecoyPeptideGroups;'))
dbDisconnect(pd)

# DECOY
psms_t[,DECOY:=0]
psms_d[,DECOY:=1]

modPeps_t[,DECOY:=0]
modPeps_d[,DECOY:=1]

pcmGrouper_psms_t[,DECOY:=0]
pcmGrouper_psms_d[,DECOY:=1]

pcms_t[,DECOY:=0]
pcms_d[,DECOY:=1]

# estimate modified peptides correctly (modified peptide (i.e. peptide group) FDR) and using PSM FDR
colnames(psms_modPeps_mapping_t)
setnames(psms_modPeps_mapping_t, 'TargetPsmsPeptideID', 'PeptideID')
psms_t <- merge(psms_t,
                psms_modPeps_mapping_t,
                by = 'PeptideID',
                all.x = T)
uniqueN(psms_t[qValue <= 0.01, TargetPeptideGroupsPeptideGroupID]) # 121700
uniqueN(modPeps_t[qValue <= 0.01, PeptideGroupID]) # 101662

round((uniqueN(psms_t[qValue <= 0.01, TargetPeptideGroupsPeptideGroupID]) /
         uniqueN(modPeps_t[qValue <= 0.01, PeptideGroupID]) * 100), 2) # 119.71

uniqueN(psms_t[qValue <= 0.01, TargetPeptideGroupsPeptideGroupID]) -
  uniqueN(modPeps_t[qValue <= 0.01, PeptideGroupID]) # 20038

bardata <- data.table('Modified peptides', 'Modified peptides', uniqueN(modPeps_t[qValue <= 0.01, PeptideGroupID]))
setnames(bardata, c('Level estimated', 'By FDR level', '# in dataset'))
bardata <- rbind(bardata,
                 list('Modified peptides', 'PSMs', uniqueN(psms_t[qValue <= 0.01, TargetPeptideGroupsPeptideGroupID])))

# PSM_ID
psms_t[,PSM_ID:=paste0(ModifiedSequence, '_', Charge, '_', StudyFileId, '_', FirstScan)]
psms_t[,PSM_J_ID:=gsub('I|L', 'J', PSM_ID)]

pcmGrouper_psms_t[,PSM_ID:=paste0(ModifiedSequence, '_', Charge, '_', StudyFileId, '_', FirstScan)]
pcmGrouper_psms_t[,PSM_J_ID:=gsub('I|L', 'J', PSM_ID)]

# PCM_ID
psms_t[,PCM_ID:=paste0(ModifiedSequence, '_', Charge)]
psms_t[,PCM_J_ID:=paste0(gsub('I|L', 'J', ModifiedSequence), '_', Charge)]

pcmGrouper_psms_t[,PCM_ID:=paste0(ModifiedSequence, '_', Charge)]
pcmGrouper_psms_t[,PCM_J_ID:=paste0(gsub('I|L', 'J', ModifiedSequence), '_', Charge)]

# MODPEP_ID
psms_t[,m_COUNT:=stringr::str_count(ModifiedSequence, 'm')]
psms_t[is.na(m_COUNT),m_COUNT:=0]
psms_t[,p_COUNT:=stringr::str_count(ModifiedSequence, 's|t|y')]
psms_t[is.na(p_COUNT),p_COUNT:=0]
psms_t[,c_COUNT:=stringr::str_count(ModifiedSequence, 'c')]
psms_t[is.na(c_COUNT),c_COUNT:=0]
psms_t[,MODPEP_ID:=paste0(Sequence, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
psms_t[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]

pcmGrouper_psms_t[,m_COUNT:=stringr::str_count(ModifiedSequence, 'm')]
pcmGrouper_psms_t[is.na(m_COUNT),m_COUNT:=0]
pcmGrouper_psms_t[,p_COUNT:=stringr::str_count(ModifiedSequence, 's|t|y')]
pcmGrouper_psms_t[is.na(p_COUNT),p_COUNT:=0]
pcmGrouper_psms_t[,c_COUNT:=stringr::str_count(ModifiedSequence, 'c')]
pcmGrouper_psms_t[is.na(c_COUNT),c_COUNT:=0]
pcmGrouper_psms_t[,MODPEP_ID:=paste0(Sequence, '_c', c_COUNT, '_m', m_COUNT, '_p', p_COUNT)]
pcmGrouper_psms_t[,MODPEP_J_ID:=gsub('I|L', 'J', MODPEP_ID)]

# add du modified peptides/pcms
colnames(psms_t)
setnames(psms_t, 'TargetPeptideGroupsPeptideGroupID', 'PeptideGroupID')
modPeps_t <- merge(modPeps_t,
                   psms_t[,list(PCM_ID, PCM_J_ID,
                                MODPEP_ID, MODPEP_J_ID,
                                PeptideGroupID)],
                   by = 'PeptideGroupID',
                   all.x = T)

colnames(psms_pcms_mapping_t)
setnames(psms_pcms_mapping_t, 'TargetPsmsPeptideID', 'PeptideID')
pcmGrouper_psms_t <- merge(pcmGrouper_psms_t,
                           psms_pcms_mapping_t,
                           by = 'PeptideID',
                           all.x = T)

colnames(pcmGrouper_psms_t)
setnames(pcmGrouper_psms_t, 'TargetPeptideGroupsPeptideGroupID', 'PeptideGroupID')
pcms_t <- merge(pcms_t,
                pcmGrouper_psms_t[,list(PCM_ID, PCM_J_ID,
                                        MODPEP_ID, MODPEP_J_ID,
                                        PeptideGroupID)],
                by = 'PeptideGroupID',
                all.x = T)

# modified peptides estimated by pcm fdr
uniqueN(modPeps_t[qValue <= 0.01, PeptideGroupID]) # 101662
uniqueN(modPeps_t[qValue <= 0.01, MODPEP_J_ID]) # 101662
uniqueN(pcms_t[qValue <= 0.01, MODPEP_J_ID]) # 113090

round((uniqueN(pcms_t[qValue <= 0.01, MODPEP_J_ID])/
         uniqueN(modPeps_t[qValue <= 0.01, PeptideGroupID]) * 100), 2) # 111.24

bardata <- rbind(bardata,
                 list('Modified peptides', 'PCMs', uniqueN(pcms_t[qValue <= 0.01, MODPEP_J_ID])))

# estimate pcms correctly (pcm (i.e. peptide group) FDR) and using PSM FDR
uniqueN(pcmGrouper_psms_t[qValue <= 0.01, PeptideGroupID])
# 551377
# need to add PCM_ID
uniqueN(pcms_t[qValue <= 0.01, PeptideGroupID]) # 519274

uniqueN(psms_t[qValue <= 0.01, PCM_J_ID]) # 142271
uniqueN(pcmGrouper_psms_t[qValue <= 0.01, PCM_J_ID]) # 142271
# psm based results align

uniqueN(modPeps_t[qValue <= 0.01, unlist(PCM_J_ID)]) # 133979
uniqueN(pcms_t[qValue <= 0.01, unlist(PCM_J_ID)]) # 132534

round((uniqueN(psms_t[qValue <= 0.01, PCM_J_ID])/
         uniqueN(pcms_t[qValue <= 0.01, unlist(PCM_J_ID)]) * 100), 2) # 107.35

round((uniqueN(modPeps_t[qValue <= 0.01, unlist(PCM_J_ID)])/
         uniqueN(pcms_t[qValue <= 0.01, unlist(PCM_J_ID)]) * 100), 2) # 101.09

bardata <- rbind(bardata,
                 list('PCMs', 'PCMs',
                      uniqueN(pcms_t[qValue <= 0.01, unlist(PCM_J_ID)])))
bardata <- rbind(bardata,
                 list('PCMs', 'Modified peptides',
                      uniqueN(modPeps_t[qValue <= 0.01, unlist(PCM_J_ID)])))

bardata <- rbind(bardata,
                 list('PCMs', 'PSMs',
                      uniqueN(psms_t[qValue <= 0.01, PCM_J_ID])))

bardata[`Level estimated` == 'Modified peptides', RELATIVE_INCREASE:= round(`# in dataset`/
                                                                              bardata[`Level estimated` == 'Modified peptides' &
                                                                                        `By FDR level`  == 'Modified peptides', `# in dataset`] * 100 - 100, 2)]

bardata[`Level estimated` == 'PCMs', RELATIVE_INCREASE:= round(`# in dataset`/
                                                                 bardata[`Level estimated` == 'PCMs' &
                                                                           `By FDR level`  == 'PCMs', `# in dataset`] * 100 - 100, 2)]

bardata[,RELATIVE_INCREASE:=paste0('+', RELATIVE_INCREASE, '%')]
bardata[RELATIVE_INCREASE == '+0%', RELATIVE_INCREASE:=NA]

fwrite(bardata, file.path(figurePath, "figure-E3G-fdr.csv"))
