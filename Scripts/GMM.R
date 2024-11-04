##----------------+
## Gaussian mixture model
##----------------+
library(mixtools)
library(cowplot)

GMM = function(data, components){
  ##-- Gaussian mixture model; are there two components?
  mixmdl = mixtools::normalmixEM(data$cnlr, k = components)
  
  plot_mix_comps = function(x, mu, sigma, lam){
    lam * dnorm(x, mu, sigma)
  }
  
  if(components == 2){
    plot = data.frame(x = mixmdl$x) %>%
      ggplot() +
      geom_histogram(aes(x, ..density..), 
                     binwidth = 0.05, 
                     colour = "purple", 
                     fill = "white", 
                     bins = 200,
                     size = 0.3) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                    colour = "red", lwd = 1.2) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                    colour = "blue", lwd = 1.2) +
      geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
                 linetype = 'dashed', 
                 size = 0.3, 
                 color = 'grey15') +
      scale_y_continuous(expand = c(0.01, 0)) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.subtitle = element_text(size = 12),
            plot.title = element_text(size = 16),
            panel.grid = element_blank()) +
      labs(y = "Density", x = 'Copy Number Log Ratio') +
      panel_border(size = 2, color = 'black')
    
  } else if (components == 3){
    
    plot = data.frame(x = mixmdl$x) %>%
      ggplot() +
      geom_histogram(aes(x, ..density..), 
                     binwidth = 0.05, 
                     colour = "purple", 
                     fill = "white", 
                     bins = 200,
                     size = 0.3) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                    colour = "red", lwd = 1.2) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                    colour = "blue", lwd = 1.2) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
                    colour = "black", lwd = 1.2) +
      geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
                 linetype = 'dashed', 
                 size = 0.3, 
                 color = 'grey15') +
      scale_y_continuous(expand = c(0.01, 0)) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.subtitle = element_text(size = 12),
            plot.title = element_text(size = 16),
            panel.grid = element_blank()) +
      labs(y = "Density", x = 'Copy Number Log Ratio') +
      panel_border(size = 2, color = 'black')
    
  } else if (components == 4){
    plot = data.frame(x = mixmdl$x) %>%
      ggplot() +
      geom_histogram(aes(x, ..density..), 
                     binwidth = 0.05, 
                     colour = "purple", 
                     fill = "white", 
                     bins = 200,
                     size = 0.3) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                    colour = "red", lwd = 1.2) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                    colour = "blue", lwd = 1.2) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
                    colour = "black", lwd = 1.2) +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(mixmdl$mu[4], mixmdl$sigma[4], lam = mixmdl$lambda[4]),
                    colour = 'yellow', lwd = 1.2) +
      geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
                 linetype = 'dashed', 
                 size = 0.3, 
                 color = 'grey15') +
      scale_y_continuous(expand = c(0.01, 0)) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.subtitle = element_text(size = 12),
            plot.title = element_text(size = 16),
            panel.grid = element_blank()) +
      labs(y = "Density", x = 'Copy Number Log Ratio') +
      panel_border(size = 2, color = 'black')
  }
  
  return(list(plot = plot,
              model = mixmdl))
}