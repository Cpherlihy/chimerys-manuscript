dataPath <- "/mnt/paper/01_paper/submission"

## load libraries
require(data.table)
require(fst)
require(fstcore)
require(foreach)
require(doParallel)
require(RSQLite)
require(ggplot2)
require(cowplot)
require(ggrepel)
require(ggrastr)
require(patchwork)
require(magick)
require(viridis)
require(scales)
require(stringr)
require(stringi)
require(checkmate)
require(rawrr)
require(ggpattern)
require(ggVennDiagram)

## source function scripts
source(here::here("scripts/msaid-theme.R"))
source(here::here("scripts/analysis.R"))
source(here::here("scripts/data-processing.R"))
