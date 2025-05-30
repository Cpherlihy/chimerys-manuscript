# Figure E7
MSAID
2025-01-24

- [Setup](#setup)
- [Data](#data)
  - [Runtimes DIA software](#runtimes-dia-software)
  - [IDs LFQ (run-specific)](#ids-lfq-run-specific)
  - [Fragment distribution](#fragment-distribution)
  - [Precursor Fragment distribution](#precursor-fragment-distribution)
  - [XICs](#xics)
  - [eFDR full quan](#efdr-full-quan)
- [Figure](#figure)

# Setup

This document describes how the data analysis and plots for extended
figure 7 were generated. To recreate the figures, make sure to download
all input files (available on
[PRIDE](https://www.ebi.ac.uk/pride/archive?keyword=PXD053241)), place
them under `dataPath` (adjust in `load-dependencies.R` to your own
folder structure) and generate intermediate results in the linked `.R`
scripts.

<details>
<summary>
Details on setup
</summary>

``` r
suppressMessages(source(here::here("scripts/load-dependencies.R")))
path <- file.path(here::here(), "figure-E7")
figurePath <- file.path(dataPath, "data/figure-E7")
msaid_quantified <- c("TRUE" = msaid_darkgray, "FALSE" = msaid_orange)
msaid_eFDR <- c("TRUE" = msaid_darkgray, "FALSE" = msaid_red)
```

</details>

# Data

<details>
<summary>
Details on data processing
</summary>

## Runtimes DIA software

Runtimes were recorded manually as `figure-E7A-runtimes.csv`

``` r
dtRuntimes <- fread(file.path(figurePath, "figure-E7A-runtimes.csv"))
softwareLevels <- c('CHIMERYS 2', 'CHIMERYS 4', 'DIA-NN', 'Spectronaut')
softwareLabels <- c('CHIMERYS\n2.7.9', 'CHIMERYS\n4.0.21', 'DIA-NN\n1.8.1', 'Spectronaut\n19')
dtRuntimes[, Software := factor(Software, softwareLevels, softwareLabels)]

p_runtimes <- ggplot(dtRuntimes, aes(x=Software, y=Runtime)) +
  geom_bar(stat = "identity", fill = msaid_darkgray) +
  geom_text(aes(x=Software, y=Runtime, label = Runtime), vjust=-0.5, family = "Montserrat Light",
            size = 5/.pt, color = msaid_darkgray) +
  theme(axis.title.x = element_blank()) + 
  ylab("Runtime [min]") + ylim(c(0, 1450))
```

## IDs LFQ (run-specific)

[R code to generate input file
`figure-E7B-efdr.csv`](figure-E7B-entrapment.R)

``` r
dtEfdrLfq <- fread(file.path(figurePath, "figure-E7B-efdr.csv"))
softwareLabels <- c("CHIMERYS", "DIA-NN", "Spectronaut")
dtEfdrLfq[, SOFTWARE := factor(SOFTWARE, softwareLabels)]
dtEfdrLfq[, isAll := factor(QUAN_COMPLETE_FDR, c(T, F))]
dtEfdrLfq[, FDR := factor(FDR, c("FDR", "FDR + eFDR"), c("FDR", "FDR +\neFDR"))]

p_efdrLfq <- ggplot(dtEfdrLfq, aes(x=FDR, y=N, fill=isAll)) +
  geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  facet_grid(cols = vars(SOFTWARE)) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual("Quantified in all replicates", values = msaid_quantified) +
  theme(legend.position = "top", legend.location = "plot") +
  xlab("Filtered to ≤1% (run-specific)") + ylab("Dataset global\nprecursors")
```

## Fragment distribution

[R code to generate input file
`figure-E7C-fragments.csv`](figure-E7C-entrapment.R)

``` r
dtQuanFrag <- fread(file.path(figurePath, "figure-E7C-fragments.csv"))
softwareLabels <- c("DIA-NN", "Spectronaut")
dtQuanFrag[, SOFTWARE := factor(SOFTWARE, softwareLabels)]
dtQuanFrag[, isEfdr001 := factor(isEfdr001, c(T, F))]

p_quanFrag <- ggplot(dtQuanFrag, aes(x=LOG10QUAN, fill=isEfdr001)) +
  geom_histogram(binwidth = 0.25) +
  facet_grid(cols = vars(SOFTWARE)) +
  scale_x_continuous(breaks = breaks_extended(6)) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual("eFDR ≤1%", values = msaid_eFDR) +
  theme(legend.position = "top") +
  xlab("Log10 intensity") + ylab("Run-specific\nfragments")
```

## Precursor Fragment distribution

[R code to generate input file
`figure-E7D-precursors.csv`](../figure-3/figure-3A-entrapment-barplot.R)

``` r
dtPrecFrag <- fread(file.path(figurePath, "figure-E7D-precursors.csv"))
softwareLabels <- c("Spectronaut")
dtPrecFrag[, SOFTWARE := factor(SOFTWARE, softwareLabels)]
dtPrecFrag[, MIN1_QUAN_FRAGS := factor(MIN1_QUAN_FRAGS, 0:6)]

p_precFrag <- ggplot(dtPrecFrag, aes(x=LOG10QUAN, fill=MIN1_QUAN_FRAGS)) +
  geom_histogram(binwidth = 0.25) +
  scale_x_continuous(breaks = breaks_extended(6)) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual("Curated fragments", values = rev(c(msaid_gradient_2b(6), msaid_red))) +
  guides(fill = guide_legend(nrow = 1, direction = "vertical")) +
  theme(legend.position = "bottom", legend.location = "plot") +
  xlab("Log10 intensity") + ylab("Run-specific\nprecursors")
```

## XICs

R code to generate input file
[`figure-E7E-rt.csv`](../figure-3/figure-3D-xics.R) and
[`figure-E7E-xic.csv`](figure-E7E-xic.R)

``` r
dtXic <- fread(file.path(figurePath, "figure-E7E-xic.csv"))
dtFrag <- fread(file.path(figurePath, "figure-E7E-rt.csv"))
rtFrag <- dtFrag[, c(unique(EG.StartRT), unique(EG.EndRT))]

p_xic <- ggplot(dtXic,aes(x=RT, y=intensity)) +
  geom_line() +
  geom_vline(xintercept = rtFrag, linetype = "dashed", linewidth = 0.25, color = msaid_darkgray) +
  facet_grid(rows = vars(precMz), cols = vars(F.FrgMz)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  xlab("Retention time [min]") + ylab("Intensity")
```

## eFDR full quan

[R code to generate input file
`figure-E7F-efdr.csv`](figure-E7F-entrapment.R)

``` r
dtEfdrFullLfq <- fread(file.path(figurePath, "figure-E7F-efdr.csv"))
efdrLabels <- c("FDR", "FDR +\neFDR", "FDR + eFDR +\nquan fragments\nfiltered")
dtEfdrFullLfq[, LABEL := factor(LABEL, efdrLabels)]
dtEfdrFullLfq[, QUAN_COMPLETE_FDR := factor(QUAN_COMPLETE_FDR, c(T, F))]

p_efdrLfqFull <- ggplot(dtEfdrFullLfq, aes(x=LABEL, y=N, fill=QUAN_COMPLETE_FDR)) +
  geom_bar(stat = "identity", position = position_stack(reverse = T)) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual("Quantified in all replicates", values = msaid_quantified) +
  theme(legend.position = "top", legend.location = "plot") +
  xlab("Filtered to ≤1% (run-specific)") + ylab("Run-specific\nprecursors")
```

</details>

# Figure

<details>
<summary>
Details on figure generation
</summary>

``` r
#p_design <- c("AAA###\nBBBCCC\nFFGGGG\nHHGGGG")
p_design <- c("AABBCC\nFFGGGG\nHHGGGG")

p_S_dia <- p_runtimes +  p_efdrLfq + p_quanFrag + p_precFrag + p_xic + p_efdrLfqFull +
  plot_layout(design = p_design) +
  plot_annotation(tag_levels = "A")

ggsave2(file.path(path, "figure-E7.pdf"), plot = p_S_dia,
        width = 180, height = 120, units = "mm", device = cairo_pdf)
ggsave2(file.path(path, "figure-E7.png"), plot = p_S_dia,
        width = 180, height = 120, units = "mm")
ggsave2(file.path(path, "figure-E7.eps"), plot = p_S_dia,
        width = 180, height = 120, units = "mm")
```

    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database

    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font family
    'Source Sans 3' not found in PostScript font database

    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Montserrat Light' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Montserrat Light' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Montserrat Light' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Montserrat Light' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Montserrat Light' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    family 'Source Sans 3' not included in postscript() device

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : font
    family 'Source Sans 3' not found in PostScript font database

    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database
    Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    font family 'Source Sans 3' not found in PostScript font database

``` r
ggsave2(file.path(path, "figure-E7.jpeg"), plot = p_S_dia,
        width = 180, height = 120, units = "mm")
ggsave2(file.path(path, "figure-E7.tiff"), plot = p_S_dia,
        width = 180, height = 120, units = "mm")
```

</details>

![figure-E7](figure-E7.png)
