##----------------+
## Plotting theme for 
## Y_chromosome loss project
##----------------+
## Started: 03/21/2021

## Libraries:
library(ggplot2)
library(sysfonts)
library(patchwork)

theme_Y = function(family = 'ArialMT',
                   panel_border = T,
                   base_size = 12, 
                   base_line_size = base_size / 22) {
  
  if(panel_border){ 
    theme_classic()  %+replace%
      theme(text = element_text(family = family,
                                face = "plain", 
                                colour = "black", 
                                size = base_size, 
                                lineheight = 0.9,
                                hjust = 0.5, 
                                vjust = 0.5, 
                                angle = 0, 
                                margin = margin(), 
                                debug = F),
            # line = element_line(colour = "black", 
            #                     size = base_line_size, 
            #                     linetype = 1, 
            #                     lineend = "round"),
            axis.text = element_text(colour = "black", 
                                     family = family, 
                                     size = rel(0.8)),
            axis.ticks = element_line(colour = "black", 
                                      size = rel(1)),
            panel.grid.major = element_blank(),
            panel.border = element_rect(fill = NA, size = base_line_size, colour = 'black'),
            panel.grid.minor = element_blank(),
            # axis.line = element_line(colour = "black"),
            #                          #size = rel(1)),
            strip.background = element_blank(),
            plot.title.position = 'plot',
            plot.caption.position = 'plot',
            plot.margin = margin(25, 25, 10, 25))
    
  } else {
    theme_classic()  %+replace%
      theme(text = element_text(family = family,
                                face = "plain", 
                                colour = "black", 
                                size = base_size, 
                                lineheight = 0.9,
                                hjust = 0.5, 
                                vjust = 0.5, 
                                angle = 0, 
                                margin = margin(), 
                                debug = F),
            line = element_line(colour = "black", 
                                size = base_line_size, 
                                linetype = 1, 
                                lineend = "round"),
            axis.text = element_text(colour = "black", 
                                     family = family, 
                                     size = rel(0.8)),
            axis.ticks = element_line(colour = "black", 
                                      size = rel(1)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black", 
                                     size = rel(1)),
            legend.key = element_blank(),
            strip.background = element_blank(),
            plot.title.position = 'plot',
            plot.caption.position = 'plot',
            legend.position = 'plot',
            plot.margin = margin(25, 25, 10, 25))
  }
    
}

#' guides for legends
guides(color = guide_colorbar(title.position = 'top', 
                              title.hjust = 0.5, 
                              barwidth = unit(20, 'lines'), 
                              barheight = unit(0.5, 'lines')))


#' save the plot in proper format
ggsave_golden = function(filename, plot, width, ...){
  ggsave(filename = filename, plot = plot, device = cairo_pdf, width = width, height = width / 1.61803398875)
}

# library(sysfonts)
# font_add_google(name = 'Lobster')





# ## broken x-axis
# x_info <- break_axis(dat$nlog10q.gene,maxlower=30,minupper=50,lowerticksize=5,upperticksize=50,ratio_lower_to_upper=0.45)
# dat$nlog10q.gene.new <- x_info$newy
# x_info_gene_labs <- break_axis(gene_labs$nlog10q.gene,maxlower=30,minupper=50,lowerticksize=5,upperticksize=50,ratio_lower_to_upper=0.45)
# gene_labs$nlog10q.gene.new <- x_info_gene_labs$newy
# x_info_residue_labs <- break_axis(residue_labs$nlog10q.gene,maxlower=30,minupper=50,lowerticksize=5,upperticksize=50,ratio_lower_to_upper=0.45)
# residue_labs$nlog10q.gene.new <- x_info_residue_labs$newy
# ## broken y-axis
# y_info <- break_axis(dat$nlog10q.residue,maxlower=10,lowerticksize=2,upperticksize=5,ratio_lower_to_upper=0.25)
# dat$nlog10q.residue.new <- y_info$newy
# y_info_residue_labs <- break_axis(residue_labs$nlog10q.residue,maxlower=10,lowerticksize=2,upperticksize=5,ratio_lower_to_upper=0.25)
# residue_labs$nlog10q.residue.new <- y_info_residue_labs$newy
# ## make plot
# p <- ggplot(dat,aes(x=nlog10q.gene.new,y=nlog10q.residue.new)) +
#   geom_hline(yintercept=sig_val,color='grey',size=0.5,linetype='dashed') +
#   geom_vline(xintercept=sig_val,color='grey',size=0.5,linetype='dashed') +
#   geom_point(aes(color=sig),size=1.2) +
#   geom_text_repel(data=residue_labs,aes(x=nlog10q.gene.new,y=nlog10q.residue.new,label=residue,color=sig),
#                   min.segment.length=0.001,angle=0,hjust=0,vjust=0.5) +
#   geom_text_repel(data=gene_labs,aes(x=nlog10q.gene.new,label=label,color=sig),y=0,angle=90,hjust=1.2,vjust=0.5) +
#   scale_x_continuous(limits=x_info$limits,labels=x_info$labels,breaks=x_info$breaks) +
#   scale_y_continuous(limits=c(-0.5,y_info$limits[2]),labels=y_info$labels,breaks=y_info$breaks) +
#   scale_color_manual(values=cols,name='Leve of enrichment\n(composites):') +
#   theme_classic(base_family='ArialMT',base_size=14) +
#   theme(legend.position='bottom') +
#   labs(x='Gene -log10(q-value)',y='Residue -log10(q-value)')


# break_axis <- function(y, maxlower, minupper=NA, lowerticksize, upperticksize, ratio_lower_to_upper) {
#   if(is.na(minupper)) {
#     breakpos <- maxlower
#     lowerticklabels <- seq(0,breakpos,by=lowerticksize); lowerticklabels
#     upperticklabels <- seq(breakpos+upperticksize,max(y)+upperticksize,by=upperticksize); upperticklabels
#     ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
#     lowertickpos <- lowerticklabels
#     uppertickspacing <- ratio_lower_to_upper * lowerticksize
#     uppertickpos <- breakpos + ((1:length(upperticklabels))*uppertickspacing)
#     tickpos <- c(lowertickpos, uppertickpos)
#     newy <- as.numeric(y)
#     ind <- newy > breakpos
#     newy[ind] <- breakpos + uppertickspacing*((newy[ind]-breakpos) / upperticksize)
#     list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
#   } else {
#     lowerticklabels <- seq(0,maxlower,by=lowerticksize); lowerticklabels
#     upperticklabels <- seq(minupper,max(y)+upperticksize,by=upperticksize); upperticklabels
#     ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
#     lowertickpos <- lowerticklabels
#     uppertickspacing <- ratio_lower_to_upper * lowerticksize
#     uppertickpos <- maxlower + 0.5*lowerticksize + ((1:length(upperticklabels))*uppertickspacing)
#     tickpos <- c(lowertickpos, uppertickpos)
#     newy <- as.numeric(y)
#     ind <- newy > maxlower
#     newy[ind] <- maxlower + 0.5*lowerticksize + 1*uppertickspacing + uppertickspacing*((newy[ind]-minupper) / upperticksize)
#     list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
#   }
# }





