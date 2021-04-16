draw_depth_plot <- function(df, colour = "black", x.title = FALSE, y.title = FALSE, legend = FALSE) {
  # First assign axis limits depending on which enzyme class we're talking about
  get_ymax <- function(x) {
    if(length(unique(x)) !=1) {
      ymax <- 80
    } else {
      ymax <- switch (unique(x)[1],
                      "peptidase" = 40,
                      "glycosylase" = 21,
                      "phosphatase" = 80
      )
    }
    ymax
  }
  max.y <- get_ymax(df$class) # Get the limit for the subpanel
  
  p <- ggplot(df, aes(x=depth.mbsf, y=v0*1000, linetype=treatment, colour = class)) + 
    geom_pointrange(aes(ymin=v0*1000-v0.se*1000, ymax=v0*1000+v0.se*1000, shape = treatment), size = 0.25) +
    geom_line(size = 0.25) +
    geom_vline(xintercept = 51, colour = "gray50") + 
    scale_x_reverse() + 
    scale_linetype_manual(values=c("live"="solid","killed"="dashed")) +
    scale_colour_manual(values = c("peptidase" = "#1b9e77", "glycosylase" = "#d95f02", "phosphatase" = "#7570b3")) +
    scale_shape_manual(values = c(1, 19)) + 
    expand_limits(xmin=0) +
    ylim(c(-2, max.y)) + 
    ylab(expression(paste(v[0], ", ", "nmol ", "substrate ", g^{-1}, " sed ", hr^{-1})))+ 
    xlab("depth, mbsf") + 
    coord_flip() +
    facet_wrap(~enzyme, nrow=1) +
    theme(axis.text.x  = element_text(angle=-45, hjust=0),
          text = element_text(size = 8),
          panel.grid.major = element_line(size = 0.25),
          panel.grid.minor = element_blank())
  if(!legend) {
    p <- p + theme(legend.position = "none")
  }
  if(!x.title) {
    p <- p + theme(axis.title.x = element_blank())
  }
  if(!y.title) {
    p <- p + theme(axis.title.y = element_blank())
  }
  #browser()
  
  # remove depth labels except for a-glucosidase and clostripain
  if(length(unique(df$substrate)) == 1) { # Then we're just looking at one substrate
    if(sum(c("clostripain", "alpha-\nglucosidase") %in% df$enzyme) < 1) { # then it is not clostripain or alpha-glucosidase
      p <- p + theme(axis.text.y = element_blank())
    }
  }
  
  p
}