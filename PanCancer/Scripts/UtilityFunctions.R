##----------------+
## Utility functions 
## for Y-chromosome Project
##----------------+
##
## start: 08/11/2021
## revision: 08/24/2022
## revision: 01/10/2023
## chris-kreitzer

CI_z = function (x, ci = 0.95){
  `%>%` <- magrittr::`%>%`
  standard_deviation = sd(x)
  sample_size = length(x)
  Margin_Error = abs(qnorm((1-ci)/2))* standard_deviation/sqrt(sample_size)
  df_out = data.frame(sample_size = length(x), 
                      Mean = mean(x), 
                      sd=sd(x),
                      Margin_Error = Margin_Error,
                      'CI lower limit'= (mean(x) - Margin_Error),
                      'CI Upper limit'= (mean(x) + Margin_Error)) %>%
    tidyr::pivot_longer(names_to = "Measurements", 
                        values_to = "values", 1:6)
  return(df_out)
}


#' calculate the confidence Interval
normConfInt = function(x, alpha = 0.05){
  mean(x, na.rm = T) + qt(1 - alpha / 2, length(x) - 1) * sd(x, na.rm = T) / sqrt(length(x)) * c(-1, 1)
}


#' calculate R^2 (how much variation in Y is explained with X (Y ~ X))
rsq = function (x, y) cor(x, y) ^ 2


##----------------+
## ggsave_golden
##----------------+
ggsave_golden = function(filename, plot, width, ...){
  ggsave(filename = filename, plot = plot, device = cairo_pdf, width = width, height = width / 1.61803398875)
}


##----------------+
## ggplot theme for LOY project
##----------------+
theme_std = function(base_size = 11, 
                     base_line_size = base_size/22, 
                     base_rect_size = base_size/22) {
  require(ggplot2)
  theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
      line = element_line(colour = "black", linewidth = base_line_size, linetype = 1, lineend = "round"),
      text = element_text(family = 'ArialMT', face = "plain",
                          colour = "black", size = base_size, lineheight = 0.9,
                          hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F),
      axis.text = element_text(colour = "black", family='ArialMT', size=rel(0.8)),
      axis.ticks = element_line(colour = "black", linewidth = rel(1)),
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = rel(1)),
      legend.key = element_blank(),
      strip.background = element_blank())
}

theme_Y = theme_bw() +
  theme(
    text = element_text(family = 'ArialMT', size = 14),
    axis.text = element_text(family = 'ArialMT', size = 14, color = 'black'),
    line = element_line(linewidth = .95*(.95/1.6), lineend = 'round'),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = 'black', linewidth = .95*(.95/1.6), lineend = 'round'),
    axis.line.y = element_line(color = 'black', linewidth = .95*(.95/1.6), lineend = 'round'),
    axis.ticks.length = unit(2, 'pt'),
    axis.ticks.x = element_line(color = 'black', linewidth = .95*(.95/1.6), lineend = 'round'),
    axis.ticks.y = element_line(color = 'black', linewidth = .95*(.95/1.6), lineend = 'round'),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.title = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(0.25,0.5,0.25,0.25), 'lines'))

  

#theme_set(theme_Y)
