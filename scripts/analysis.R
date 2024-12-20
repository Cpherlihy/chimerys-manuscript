
plotUpset <- function(dataTable, groupColumn = "condition",
                      observationColumn = "ptm_group", observationLabel = NULL,
                      maxIntersections = 10L, orderUpset = c("frequency", "degree"),
                      labelType = c("total inside", "both outside", "off"),
                      labelDistance = c(0.1, 0.1), labelSize = c(5L, 5L),
                      plotRelativeWidths = c(0.4, -0.25, 0.6),
                      plotRelativeHeights = c(0.45, -0.125, 0.55),
                      groupFillColumn = NULL, groupFillColors = NULL, returnList = F) {
  stopifnot(is.data.table(dataTable))
  stopifnot(is.character(groupColumn), length(groupColumn)==1L)
  stopifnot(is.character(observationColumn), length(observationColumn)==1L)
  stopifnot(all(c(groupColumn, observationColumn, groupFillColumn) %in% names(dataTable)))
  if(!is.null(observationLabel)) stopifnot(is.character(observationLabel), length(observationLabel)==1L)
  if(is.null(observationLabel)) observationLabel <- observationColumn
  if(!any((is.integer(maxIntersections) && length(maxIntersections)==1L) || is.null(maxIntersections))) {
    stop("maxIntersections needs to be an integer of length 1 or NULL")
  }
  orderUpset <- match.arg(orderUpset)
  labelType <- match.arg(labelType)
  stopifnot(is.numeric(labelDistance), length(labelDistance)==2L)
  stopifnot(is.numeric(labelSize), length(labelSize)==2L)
  stopifnot(is.numeric(plotRelativeWidths), length(plotRelativeWidths)==3L)
  stopifnot(is.numeric(plotRelativeHeights), length(plotRelativeHeights)==3L)
  if(!is.null(groupFillColumn)) stopifnot(is.character(groupFillColumn), length(groupFillColumn)==1L)
  if(!is.null(groupFillColors)) stopifnot(is.character(groupFillColors))
  stopifnot(is.logical(returnList), length(returnList)==1L)

  theme_color <- function(color) {
    theme(plot.background = element_rect(color  = color, linetype = 'dotted', fill = color))
  }


  #prepare data.table
  column_names <- c(groupColumn, observationColumn, groupFillColumn)
  dataTable <- dataTable[, .SD, .SDcols = column_names]
  if(is.null(groupFillColumn)) {
    setnames(dataTable, column_names, c("group", "observation"))
  } else {
    setnames(dataTable, column_names, c("group", "observation", "group_fill"))
  }
  if(is.null(groupFillColors)) {
    groupFillColors <- msaid_col
  }

  #total IDs bar plot
  if(is.null(groupFillColumn)) {
    count_total <- dataTable[, .N, by = group]
  } else {
    count_total <- dataTable[, .(.N, group_fill = as.character(group_fill[1])), by = group]
    count_total[, group_fill := factor(group_fill, names(groupFillColors))]
  }
  order_variable <- count_total[order(-N, group), group]
  count_total[, group := factor(group, levels = order_variable)]
  setkey(count_total, group)

  p_total_bar <- ggplot(count_total, aes(y=group, label=format(N, big.mark=",", trim=T))) +
    scale_y_discrete(position = "right", limits = rev) +
    scale_x_reverse(labels = label_number(scale_cut = cut_short_scale())) +
    xlab(paste("Total unique", observationLabel)) + ylab(NULL) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.text.y.right = element_text(hjust=0.5),
          axis.title.x = element_text(hjust=0),
          plot.background = element_rect(fill = "transparent", color = NA))
  if(is.null(groupFillColumn)) {
    p_total_bar <- p_total_bar +
      geom_bar(aes(x=N), stat="identity")
  } else {
    p_total_bar <- p_total_bar +
      geom_bar(aes(x=N, fill=group_fill), stat="identity") +
      scale_fill_manual("", values = groupFillColors) +
      theme(legend.position = "none")
  }
  if(labelType=="total inside") {
    p_total_bar <- p_total_bar +
      geom_text(aes(x=max(N)*labelDistance[1]),
                color="white", family="Montserrat Light",
                size = labelSize[1]/.pt, hjust = 1)
  } else if(labelType=="both outside") {
    p_total_bar <- p_total_bar +
      geom_text(aes(x=N+max(N)*labelDistance[1]),
                color=msaid_darkgray, family="Montserrat Light",
                size = labelSize[1]/.pt, hjust = 1)
  }


  #upset top bar plot
  dataTable[, is_identified := 1L]
  dataTable <- dcast(dataTable, observation~group, value.var = "is_identified", fill = 0)

  count_wide <- dataTable[, .N, keyby = c(names(dataTable)[-1])]
  count_wide[, upset := apply(.SD, 1, function(x) paste(as.integer(x), collapse="_") ),
             .SDcols = names(count_wide)[!names(count_wide) %in% "N"]]
  count_wide[, display := T]
  #reorder table by frequency or degree
  if(orderUpset=="frequency") {
    upset_order <- count_wide[order(-N, upset), upset]
    count_wide[, upset := factor(upset, levels = upset_order)]
    setkey(count_wide, upset)
  } else if(orderUpset=="degree") {
    group_col <- names(count_wide)[(ncol(count_wide)-3):1]
    count_wide[, group_sum := apply(.SD, 1, sum), .SDcols = c(group_col)]
    setorderv(count_wide, order = 1L, cols = c("group_sum", group_col))
    upset_order <- count_wide[, upset]
    count_wide[, upset := factor(upset, levels = upset_order)]
    setkey(count_wide, upset)
  }
  #filter max intersections
  n_intersections <- count_wide[, .N]
  n_observations <- count_wide[, sum(N, na.rm = T)]
  if(!is.null(maxIntersections) && count_wide[, .N]>maxIntersections) {
    count_wide[!1:maxIntersections, display := F]
  }

  #adjust y axis label margin for better plot alignment
  p_upset_bar <- ggplot(count_wide[display==T],
                        aes(x=upset, label=format(N, big.mark=",", trim=T))) +
    geom_bar(aes(y=N), stat="identity") +
    xlab(NULL) + ylab(paste("Intersections\n", observationLabel)) +
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    theme(axis.text.x = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt"))
  if(labelType %in% c("total inside", "both outside")) {
    p_upset_bar <- p_upset_bar +
      geom_text(aes(y=N+max(N)*labelDistance[2]),
                color=msaid_darkgray, family="Montserrat Light",
                size = labelSize[2]/.pt)
  }


  #upset bottom legend
  count_long <- melt(count_wide[display==T], measure.vars = names(count_wide)[
    !names(count_wide) %in% c("N", "upset", "display")])
  count_long[, variable := factor(variable, levels = rev(order_variable))]
  order_upset <- count_long[order(-N, upset), unique(upset)]
  count_long[, upset := factor(upset, levels = order_upset)]
  upset_line <- count_long[value==1, .(y_min = min(as.integer(variable)),
                                       y_max = max(as.integer(variable))), by = upset]

  p_upset_legend <-
    ggplot(count_long, aes(x=upset, xend=upset)) +
    geom_point(aes(y=variable, color=as.character(value), alpha=as.character(value))) +
    geom_segment(data=upset_line, aes(y=y_min, yend=y_max), color=msaid_darkgray) +
    scale_color_manual("", values = c("1" = msaid_darkgray, "0" = msaid_gray)) +
    scale_alpha_manual("", values = c("1" = 1, "0" = 0.2)) +
    theme(legend.position = "none",
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.text.y = NULL,
          plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent")
          ) +
    xlab(NULL) + ylab(NULL)


  #patchwork plotting instead of cowplot (has spacing issues)
  if(count_wide[display==T, .N] < n_intersections) {
    label_summary <-
      paste0(maxIntersections, " intersections out of ",
             n_intersections,
             "\nordered by ", orderUpset,
             "\n(= ", format(count_wide[display==T, sum(N, na.rm=T)], big.mark=",", trim=T),
             " ", observationLabel, "\nout of ",
             format(n_observations, big.mark=",", trim=T), ")")
  } else {
    label_summary <-
      paste(n_intersections, "intersections\n\n",
            "ordered by", orderUpset, "\n\n",
            format(n_observations, big.mark=",", trim=T), "observations")
  }

  p_description <- ggplot_placeholder(label_summary, 5L)
  p_upset <- wrap_elements(p_description) + plot_spacer() + p_upset_bar +
    plot_spacer() + plot_spacer() + plot_spacer() +
    p_total_bar + plot_spacer() + p_upset_legend +
    plot_layout(widths = plotRelativeWidths, heights = plotRelativeHeights)

  if(returnList==T) {
    return(list("plot" = p_upset,
                "count_total" = count_total,
                "count_intersections" = count_wide))
  } else {
    return(p_upset)
  }
}



plotMirror <- function(dataTable, rawSample, rawPath, rawScans, outputPath, outputName = "mirror_plots.pdf",
                       loadPredictions = T, ppmShift = NULL, ppmTolerance = 20,
                       rawActivation = c("HCD", "CID"), removePtmJ = T, trimPredicted = T,
                       inferysApi = "localhost:8081",
                       inferysModel = "inferys_3.0.0_fragmentation",
                       quiet = F) {
  stopifnot(is.data.table(dataTable))
  col_vec <- c("sample", "scan_ms2", "mz_ratio", "ptm", "charge", "nce", "score_coefficient_normalized")
  if(any(!col_vec %in% names(dataTable))) {
    stop(paste(paste(col_vec, collapse=", "), "need to be part of dataTable"))
  }
  stopifnot(is.character(rawSample), length(rawSample)==1L)
  stopifnot(is.character(rawSample), length(rawSample)==1L)
  stopifnot(file_test("-f", rawPath), length(rawPath)==1L)
  stopifnot(is.integer(rawScans))
  stopifnot(file_test("-d", outputPath), length(outputPath)==1L)
  stopifnot(is.character(outputName), length(outputName)==1L)
  stopifnot(is.logical(loadPredictions), length(loadPredictions)==1L)
  if(!is.null(ppmShift)) {
    stopifnot(is.numeric(ppmShift), length(ppmShift)==1L)
  } else {
    ppmShift <- 0
  }
  stopifnot(is.numeric(ppmTolerance), length(ppmTolerance)==1L)
  rawActivation <- match.arg(rawActivation)
  stopifnot(is.logical(removePtmJ), length(removePtmJ)==1L)
  stopifnot(is.logical(trimPredicted), length(trimPredicted)==1L)
  stopifnot(is.character(inferysApi), length(inferysApi)==1L)
  stopifnot(is.character(inferysModel), length(inferysModel)==1L)


  #format identified peptides
  if(!quiet) print("get predictions")
  rawScans <- unique(rawScans)
  peptides <- dataTable[sample %in% rawSample & scan_ms2 %in% rawScans,
                        .(scan_ms2, mz_ratio, ptm, charge, nce,
                          score_coefficient_normalized)]
  peptides[, scan_ms2 := factor(scan_ms2, levels = rawScans)]
  #remove IL ambiguous peptides and re-normalize normalized CHIMERYS coefficient
  if(removePtmJ==T) {
    peptides[, ptm_J := gsub("(?<!\\[..)(I|L)", "J", ptm, perl=T)]
    peptides <- peptides[, .SD[1], by=.(scan_ms2, ptm_J, charge)]
    peptides[, score_coefficient_normalized := score_coefficient_normalized/sum(score_coefficient_normalized), by=scan_ms2]
  }
  peptides[, ptm := gsub("[unimod:", "[UNIMOD:", ptm, fixed = T)]
  peptides[, score_coefficient_normalized_scaled := score_coefficient_normalized/max(score_coefficient_normalized),
           by = scan_ms2]
  setkey(peptides, scan_ms2, ptm)
  peptides[, id := 1L:.N]
  setkey(peptides, id)

  if(loadPredictions==TRUE) {
    predictions <- fread(file.path(path, "predictions.csv"))
  } else {
    require(chimerysHelpeR)
    #INFERYS predictions
    predictions <- predictSpectrum(sequences = peptides$ptm,
                                   charges = peptides$charge,
                                   collisionEnergies = peptides$nce,
                                   activationType = rawActivation,
                                   model = inferysModel,
                                   spectrumApiIP = inferysApi,
                                   simplify = T)
    predictions <- rbindlist(predictions, idcol = "id")
    setkey(predictions, id)
    predictions <- peptides[predictions]
    predictions[, ion := gsub("^(.).*$", "\\1", annotation)]
    predictions[, fraction := intensity*score_coefficient_normalized_scaled]
    fwrite(predictions, file.path(path, "predictions.csv"))
  }

  #extract peaks from raw file
  if(!quiet) print("get spectra")
  spectra <- readSpectrum(rawPath, scan = rawScans)
  spec_cols <- c("scan", "scanType", "charge", "rtinseconds", "monoisotopicMz", "Monoisotopic M/Z:", "pepmass")
  metas <- rbindlist(lapply(spectra, function(x) x[spec_cols] ))
  metas_mz <- rbindlist(lapply(spectra, function(x) data.table(mz_start = x$massRange[1], mz_stop = x$massRange[2])))
  metas <- cbind(metas, metas_mz)
  peaks <- rbindlist(lapply(spectra, function(x) x[c("scan", "centroid.mZ", "centroid.intensity")] ))
  peaks[, scan := factor(scan, levels = rawScans)]
  setnames(peaks, c("scan", "mz", "intensity"))
  peaks[, fraction := numeric()]


  #merge predictions and peaks
  merged <- rbind(predictions[, .(scan = scan_ms2, origin = "predicted", label = ptm, mz, intensity, fraction, annotation)],
                  peaks[, .(scan, origin = "measured", label = "measured", mz, intensity, fraction, annotation = NA)])

  if(!quiet) print("plot spectra")
  if(!quiet) progress_bar <- txtProgressBar(min = 0, max = nrow(metas), style = 3)
  cairo_pdf(file.path(outputPath, outputName), width = 15, height = 10, onefile = T)
  temp <- foreach(peptide = split(peptides, by = "scan_ms2"),
                  spectrum = split(merged, by = "scan"),
                  meta = split(metas, by = "scan"), .combine = c) %do%
    {
      #trim predicted peaks to measured scan range and at least 1% of base peak intensity
      if(trimPredicted==T) {
        spectrum <- spectrum[origin %in% "measured" | (origin %in% "predicted" & mz >= meta[1, mz_start] & mz <= meta[1, mz_stop])]
        spectrum <- spectrum[origin %in% "measured" | (origin %in% "predicted" & intensity >= 0.01)]
      }

      #factor peptide sequences and labeling
      ptm_orderd <- peptide[order(-score_coefficient_normalized_scaled), ptm]
      ptm_label <- paste0(ptm_orderd, " (", peptide[order(-score_coefficient_normalized_scaled), charge], "+)")
      spectrum[, label := factor(label, levels = c("measured", "precursor", ptm_orderd),
                                 labels = c("measured", "precursor", ptm_label))]
      #re-calibrate mz and calculate mz boundaries
      spectrum[, mz_recal := mz/(1 + ppmShift/1e6)]
      spectrum[, mz_lower := mz_recal - (mz*ppmTolerance/1e6)]
      spectrum[, mz_upper := mz_recal + (mz*ppmTolerance/1e6)]
      #print(spectrum[, .(mz, mz_recal, mz_lower, mz_lower_diff = mz-mz_lower, mz_upper, mz_upper_diff = mz_upper-mz)])

      #identify raw peaks (include precursors)
      mz_predicted <- spectrum[origin=="predicted", mz]
      spectrum[origin=="measured", is_identified := any(mz_lower <= mz_predicted & mz_upper >= mz_predicted), by=mz]
      #identify precursors
      mz_precursor <- c(peptide[, mz_ratio], peptide[, mz_ratio-(18/charge)], peptide[, mz_ratio-(17/charge)])
      spectrum[origin=="measured", is_precursor := any(mz_lower <= mz_precursor & mz_upper >= mz_precursor), by=mz]
      spectrum[is_precursor==T, is_identified := T]
      #spectrum[is_precursor==T, annotation := "prec"]
      spectrum[is_precursor==T, label := "precursor"]
      #identify predictions
      mz_raw <- spectrum[origin=="measured", mz]
      spectrum[origin=="predicted", is_identified := any(mz_lower <= mz_raw & mz_upper >= mz_raw), by=.(label, mz, annotation)]

      #re-scale measured segments to set top identified peak to 100% of top identified predicted peak
      spectrum[, is_truncated := F]
      rel_predicted <- spectrum[origin=="predicted" & is_identified==T][order(-fraction)][1, fraction]
      if(!is.na(rel_predicted)) {
        intensity_raw <- spectrum[is_precursor==F & is_identified==T][order(-intensity)][1, intensity]
        spectrum[origin=="measured", fraction := -intensity/intensity_raw*rel_predicted]
        spectrum[origin=="measured" & intensity>intensity_raw, is_truncated := T]
        spectrum[origin=="measured" & intensity>intensity_raw, fraction := -1]
      } else {
        spectrum[origin=="measured", fraction := -intensity/max(intensity)]
      }

      #stack identical mz on top of each other
      setkey(spectrum, origin, mz, fraction)
      spectrum[, fraction_start := c(0, cumsum(fraction)[-.N]), by=.(mz, origin)]
      spectrum[, fraction_end := cumsum(fraction), by=.(mz, origin)]
      spectrum[!is.na(annotation), annotationMax := ifelse(fraction_end==max(fraction_end), annotation, NA), by=mz]

      #define color scale
      if(nrow(peptide)<=5) {
        spec_color <- c(msaid_gray, msaid_red, msaid_col[1:nrow(peptide)])
      } else {
        spec_color <- c(msaid_gray, msaid_red, msaid_col[-6], viridis((nrow(peptide)-5)))
      }
      names(spec_color) <- c("measured", "precursor", ptm_label)
      #define plot title
      nce_label <- peptide[, .(nce = unique(nce)), keyby=charge][, paste0(nce, "% (", charge, "+)")]
      plot_predicted <- meta[, paste0("predicted NCE: ", paste(nce_label, collapse = ", "))]
      plot_raw <- meta[, paste0("scan ", format(scan, big.mark = ","),
                                " at ", round(rtinseconds/60, 2),
                                " min (", gsub("^(.*) \\+.*(Full .*) \\[.*\\]$", "\\1 \\2", scanType), ")")]

      #mirror plot
      p <- ggplot(spectrum, aes(x=mz, xend=mz, color=label, label=annotationMax)) +
        geom_segment(data = spectrum[is_truncated==F], aes(alpha=is_identified, y=fraction_start, yend=fraction_end)) +
        geom_segment(data = spectrum[is_truncated==T], aes(alpha=is_identified, y=fraction_start, yend=fraction_end),
                     arrow = arrow(length = unit(0.1, "inches")), show.legend = F) +
        geom_text_repel(aes(y=fraction_end), family = "Montserrat Light", size = 5/.pt, show.legend = F,
                        na.rm = T, nudge_y = 0.01, box.padding = 0.05, max.overlaps = 10,
                        segment.size = 0.2, segment.alpha = 0.5,
                        segment.linetype = "dashed", ylim = c(0, 1)) +
        annotate("text", x=max(spectrum$mz)/2, y=1.05, label=plot_predicted,
                 family="Montserrat Light", color=msaid_darkgray, size = 6/.pt) +
        annotate("text", x=max(spectrum$mz)/2, y=-1.05, label=plot_raw,
                 family="Montserrat Light", color=msaid_darkgray, size = 6/.pt) +
        xlab("m/z") + ylab("Relative abundance") +
        scale_color_manual("", values = spec_color) +
        scale_alpha_manual(guide = "none", values = c("TRUE" = 1, "FALSE" = 0.2)) +
        scale_y_continuous(labels = function(x) label_percent()(abs(x))) + #limits = c(-1.05, 1.05)
        theme(legend.position = "top")
      suppressWarnings(print(p))
      if(!quiet) setTxtProgressBar(progress_bar, which(metas$scan %in% meta[1, scan])[1])
      return(list(list(spectrum, p)))
    }
  dev.off()

  #format and return output
  data_spectra <- rbindlist(lapply(temp, function(x) x[[1]] ))
  setkey(data_spectra, scan, origin, label, mz)
  plot_list <- lapply(temp, function(x) x[[2]] )
  names(plot_list) <- rawScans
  return_spectra <- list(data_spectra, plot_list)
  names(return_spectra) <- c("spectra", "plots")
  if(!quiet) close(progress_bar)
  return(return_spectra)
}


#combine two strings
stringCombine <- function(stringList, link = "_", stringPriority=c("first", "last")) {
  if(!is.list(stringList)) stop("'stringList' must be a list of vectors")
  stringPriority <- match.arg(stringPriority)

  if(stringPriority=="first") {
    apply(expand.grid(rev(stringList))[, length(stringList):1], 1, paste, collapse=link)
  } else {
    apply(expand.grid(stringList), 1, paste, collapse=link)
  }
}
