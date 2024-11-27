#load packages
path <- file.path(here::here(), "figure-2")
suppressMessages(source(here::here("scripts/load-dependencies.R")))
suppressMessages(source(here::here("scripts/data-ids.R")))
msaid_SE <- c("Sequest HT" = msaid_orange,
              "MSFragger" = msaid_green,
              "CHIMERYS" = msaid_blue,
              "CHIMERYS\ntop 15 peaks" = msaid_lightblue)

#load spectra per IW files
filePaths <- list.files(file.path(path, "spectra-per-IW"), full.names = T)
data_spectra <- foreach(i = 1:length(filePaths), .combine = rbind) %do% {
  temp <- fread(filePaths[i])
  return(cbind(file_name = basename(filePaths[i]), temp))
}
# Filter the data to exclude values equal to 0
data_spectra <- data_spectra[`Number of PSMs`>0]

IW_levels <- c(1, 3, 6, 8, 10, 12, 15, 20)
IW_labels <- c("1.4", "3.4", "6.4", "8.4", "10.4", "12.4", "15.4", "20.4")
data_spectra[, condition_IW := factor(gsub("^.*IW(.*)-\\(2\\).*$", "\\1", file_name),
                                      levels = IW_levels, labels = IW_labels)]
count_psms_levels <- c(16:1)
count_psms_labels <- c(rep("≥7", length(count_psms_levels)-6), 6:1)
data_spectra[, counts_PSMs := factor(`Number of PSMs`, count_psms_levels, count_psms_labels)]

count_spectra <- data_spectra[, .N, keyby=.(condition_IW, counts_PSMs)]
fwrite(count_spectra, file.path(path, "spectra-per-IW.csv"))
