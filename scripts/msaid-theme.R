#info on custom fonts
#http://zevross.com/blog/2014/07/30/tired-of-using-helvetica-in-your-r-graphics-heres-how-to-use-the-fonts-you-like-2/
#https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/

#load themes and fonts first, to ensure fonts registration within ggplot2
#old way to register Montserrat using showtext
#problem: rasterizes all text as paths, cannot be edited in e.g. Illustrator afterwards
#require(showtext)
# font_add_google(name = "Montserrat", family = "mont", regular.wt = 300, bold.wt = 500)
# showtext_auto() #call showtext to print whenever mont is reffered to

#new way to register Montserrat using extrafont -> requires font to be available on Linux
#loading extrafont automatically registers Montserrat, if steps below have been executed at least once
#require(extrafont)
#see https://www.r-project.org/nosvn/pandoc/extrafont.html for more info
#font_import: only execute ONCE after downloading Montserrat to register in R
#then confirm that Motserrat is listed in fonts()
#IMPORTANT: make sure Montserrat font files are not duplicated, otherwise loadfonts() will generate warning:
#"More than one version [...] Skipping setup for this font" -> will not properly function with pdf()
# font_import(paths = "/usr/share/fonts/googlefonts/montserrat/static", prompt = F)
# loadfonts()
# fonts()
# fonttable()
#for normal text use Montserrat Light (= 300)
#for bold text use Montserrat Medium (= 500)
#IMPORTANT: default ggsave() and pdf() throw error "font width unknown for character [...]"
#use ggsave(device = cairo_pdf) and cairo_pdf() instead

#to avoid warning 'font family 'Montserrat' not found in PostScript font database' for cowplot::plot_grid
set_null_device(cairo_pdf)

#msaid theme colors
msaid_blue = "#14A4D9"
msaid_blue_dark <- "#0B5D81"
msaid_orange = "#F06F13"
msaid_yellow = "#F19D1E"
msaid_gray = "#647A84"
msaid_gray_light2 <- "#E3E8EA"
msaid_darkgray = "#3E474C"
msaid_green = "#5DC21F"
msaid_lightblue = "#BCE8FB"
msaid_purple = "#E3449E"
msaid_red = "#FD4A36"
msaid_red_dark <- "#A32011"
msaid_col <- c(msaid_blue, msaid_orange, msaid_green, msaid_yellow, msaid_purple, msaid_red)
msaid_gradient_2 <- colorRampPalette(c(msaid_blue, msaid_orange))
msaid_gradient_2b <- colorRampPalette(c(msaid_darkgray, msaid_orange))


#msaid theme set
theme_set(theme_bw(base_family = "Source Sans 3", base_size = 7) +
            theme(text = element_text(colour = msaid_darkgray),
                  plot.title = element_text(family = "Montserrat Medium"),
                  axis.text = element_text(colour = msaid_darkgray),
                  panel.grid = element_blank(),
                  panel.border = element_rect(colour = msaid_darkgray),
                  axis.ticks = element_line(colour = msaid_darkgray),
                  axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
                  axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
                  legend.key.size = unit(3, "mm")))


## custom ggplot themes and functions

#write face text horizontally
theme_facet_horizontal <- theme(strip.text.y = element_text(angle=0))

#tilt x-axis labels
theme_tilt_xaxis <- function(angle, set_vjust = F) {
  if(set_vjust){
    theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 0.5))
  } else {
    theme(axis.text.x = element_text(angle = angle, hjust = 1))
  }
}

#x-axis dodge
theme_dodge_x_axis <- scale_x_discrete(guide = guide_axis(n.dodge=2))

#ggplot shapes reference plot
show_shape <-
  ggplot(data.table(x = rep(1:5, 5), y = rep(1:5, each = 5), shape = 1:25, label = 1:25),
         aes(x = x, y = y, shape = shape, label = label)) +
  geom_point(size = 9, fill = msaid_orange, color = msaid_blue, stroke = 2) +
  geom_text(aes(x = x+0.2, y = y+0.3)) + scale_shape_identity() +
  theme(axis.text = element_blank())

#empty plot placeholder
ggplot_placeholder <- function(placeholderLabel, size = 11L, family = "Source Sans 3") {
  ggplot() +
    annotate("text", label=placeholderLabel, x=0, y=0, size=size/.pt,
             family=family, color=msaid_darkgray) +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"))
}

#magick image to ggplot without text theme formatting
image_ggplot2 <- function(image, interpolate = FALSE) {
  info <- image_info(image)
  ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    geom_blank() +
    coord_fixed(expand = FALSE, xlim = c(0, info$width),
                ylim = c(0, info$height)) +
    annotation_raster(image, 0, info$width, info$height, 0, interpolate = interpolate) +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"))
}

#cowplot draw ggplot without text theme formatting
ggdraw2 <- function (plot = NULL, xlim = c(0, 1), ylim = c(0, 1), clip = "off") {
  p <- ggplot() + coord_cartesian(xlim = xlim, ylim = ylim,
                                  expand = FALSE, clip = clip) + scale_x_continuous(name = NULL) +
    scale_y_continuous(name = NULL) +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"))
  if (!is.null(plot)) {
    p <- p + draw_plot(plot)
  }
  p
}

#fixes error ""
cut_short_scale2 <- function (space = FALSE) {
  out <- c(0, 1, K = 1000, M = 1e+06, B = 1e+09, T = 1e+12)
  if (space) {
    names(out) <- paste0(" ", names(out))
  }
  out
}
