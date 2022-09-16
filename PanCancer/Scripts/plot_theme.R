## ggplot theme
theme_std <- function(base_size = 11, 
                      base_line_size = base_size/22, 
                      base_rect_size = base_size/22) {
  require(ggplot2)
  theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
      line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"),
      text = element_text(family = 'ArialMT', face = "plain",
                          colour = "black", size = base_size, lineheight = 0.9,
                          hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F),
      axis.text = element_text(colour = "black", family='ArialMT', size=rel(0.8)),
      axis.ticks = element_line(colour = "black", size=rel(1)),
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = rel(1)),
      legend.key = element_blank(),
      strip.background = element_blank())
}