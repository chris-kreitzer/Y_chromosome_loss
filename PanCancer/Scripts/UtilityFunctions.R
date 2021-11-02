## Utility functions


CI_z <- function (x, ci = 0.95){
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
