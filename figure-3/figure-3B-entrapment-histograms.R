# ---- setup ----
source(here::here("scripts/data-processing-2.R"))
path <- file.path(here::here(), "figure-3")
figurePath <- file.path(dataPath, "data/figure-3")

# ---- pathsToData ----
## combined
pathToCombined <- file.path(figurePath, 'fst-backup/20241127_figure3a_combined_pcms_localPcmGrouper_apexQuan_pepEntr1.fst')

# ---- read data ----
combined <- read.fst(pathToCombined, as.data.table = T)

# ---- change ENTRAPMENT_Q_VALUE for other search engines; not including decoys ----
combined[SOFTWARE != 'CHIMERYS', ENTRAPMENT_Q_VALUE:=ENTRAPMENT_Q_VALUE_1]

# ---- preliminary plot ----
ggplot(combined[!SOFTWARE == 'SPECTRONAUT_FILTERED' &
                  Q_VALUE <= 0.01],
       mapping = aes(log10(QUAN),
                     fill = factor(ENTRAPMENT_Q_VALUE <= 0.01))) +
  geom_histogram(position = 'identity', bins = 200, alpha = 0.5) +
  facet_wrap(.~factor(SOFTWARE,
                      levels = c('CHIMERYS',
                                 'DIA-NN',
                                 'SPECTRONAUT',
                                 'SPECTRONAUT_FILTERED')),
             # scales = 'free',
             strip.position = 'bottom') +
  theme(text = element_text(size = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = .5),
        legend.margin = margin(c(0,0,0,0)),
        legend.position = 'top',
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        strip.placement = "outside",
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5)) +
  scale_y_continuous(labels = label_number(suffix = "k", scale = 1e-3, sep = '')) +
  labs(y = '# precursors in dataset',
       fill = 'Surviving 1% entrapment FDR')
