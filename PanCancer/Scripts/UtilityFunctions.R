##----------------+
## Utility functions 
## for Y-chromosome Project
##----------------+
##
## start: 08/11/2021
## revision: 08/24/2022
## revision: 01/10/2023
## revision: 01/16/2023
## 
## chris-kreitzer


##----------------+
## Required packages
##----------------+
suppressPackageStartupMessages({
  library(plyr)
  library(tidyverse)
  library(patchwork)
  library(ggsignif)
  library(binom)
  library(stringi)
  library(ggrepel)
  library(grid)
  library(readxl)
  library(ggbeeswarm)
  library(here)
  library(survival)
  library(survminer)
})


##----------------+
## Current cohort
##----------------+
cohort = readRDS('~/Documents/MSKCC/10_MasterThesis/Data/00_CohortData/Cohort_071322.rds')
cohort_study = cohort[which(cohort$Study_include == 'yes'), ]

##----------------+
## Cancer-types to include
##----------------+
ctypes = table(cohort_study$CANCER_TYPE)
ctypes_keep = names(ctypes[which((ctypes / length(unique(cohort_study$SAMPLE_ID)) * 100) > 0.2)])
all_ctypes = unique(cohort_study$CANCER_TYPE)


##----------------+ 
## Convenience functions; 
##----------------+
as_perc = function(x) { 100*x }
as_points = function(x, input_unit = 'mm') { as.numeric(convertUnit(unit(x, 'pt'), input_unit)) }
ci_upper = function(n, N) { binom.confint(n, N, methods = 'wilson')[['upper']] }
ci_lower = function(n, N) { binom.confint(n, N, methods = 'wilson')[['lower']] }
mcn_test = function(a, b, c, total) {
  mcnemar.test(matrix(c(a, b, c, total - (a + b + c)), ncol = 2))[['p.value']]
}
f_test = function(a, b, c, d) {
  fisher.test(matrix(c(a, b, c, d), ncol = 2))[['p.value']]
}


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


#' out