# ---- setup ----
source(here::here("scripts/load-dependencies.R"))
source(here::here("scripts/data-processing-3.R"))
path <- file.path(here::here(), "figure-3")
figurePath <- file.path(dataPath, "data/figure-3")


dt <- fread(file.path(dataPath, "data/figure-E7/figure-E7E-rt.csv"))
mzPrec <- dt[, FG.PrecMz[1]]

dataPath <- "/mnt/paper/01_paper/figures/main/5_figure_DIA"
rawFiles <- c("LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01.raw",
              "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_02.raw",
              "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_02.raw",
              "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_03.raw")
rawPaths <- file.path(dataPath, rawFiles)

raw_scans <- fread(file.path(dataPath, "data/figure-S2/intermediate/raw_scans.csv"))

raw_scans <- raw_scans[StartTime>=23 & StartTime<=26 & MSOrder=="Ms2" &
                         precursorMass-4<=mzPrec & precursorMass+4>=mzPrec]
raw_scans <- raw_scans[file==rawFiles[6]]
scanTypes <- raw_scans[, unique(scanType)]
precMasses <- paste0("Scan filter\n", raw_scans[, unique(gsub("^FTMS \\+ c NSI Full ms2 (.*) \\[.*$", "\\1", scanType))])

lsXic1 <- readChromatogram(rawfile = rawPaths[6], mass = dt$F.FrgMz, tol = 20, filter = scanTypes[1], type = "xic")
names(lsXic1) <- round(dt$F.FrgMz, 4)
dtXic1 <- lapply(lsXic1, function(x) {
  data.table(precMz = precMasses[1], RT = x$times, intensity = x$intensities)
})
dtXic1 <- rbindlist(dtXic1, idcol = "F.FrgMz")


lsXic2 <- readChromatogram(rawfile = rawPaths[6], mass = dt$F.FrgMz, tol = 20, filter = scanTypes[2], type = "xic")
names(lsXic2) <- round(dt$F.FrgMz, 4)
dtXic2 <- lapply(lsXic2, function(x) {
  data.table(precMz = precMasses[2], RT = x$times, intensity = x$intensities)
})
dtXic2 <- rbindlist(dtXic2, idcol = "F.FrgMz")

dtXic <- rbind(dtXic1, dtXic2)
dtXic <- dtXic[RT>=23 & RT<=26]

fwrite(dtXic, file.path(figurePath, "figure-E7E-xic.csv"))
