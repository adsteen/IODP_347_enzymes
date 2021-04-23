###########
# Calculate predicted cell counts at enzyme sampling depths from Joy's cell count data
###########

# Load necessary libraries
#require(readr)
#require(plyr)
#require(ggplot2)


# Read data
cell_d <- read.csv("data/BSB_59_cell_counts.csv")

# Note enzyme depths
enz_depths <- data.frame(depth.m = c(4.5, 11.1, 17.6, 24.3, 30.9, 37.5,43.15, 48.22, 54.95, 61.55, 68.12, 77.92))

# Set up model parameters - span and unconformity depth
unconformity <- 51.68
spn1 <- 0.53
spn2 <- 0.8

# Make two loess models - 1 for above the unconformity, and one below.
m_shallow <- loess(cells ~ depth.m, data = subset(cell_d, depth.m <= unconformity), span=spn1)
m_deep <- loess(cells ~ depth.m, data=subset(cell_d, depth.m > unconformity), span=spn2)

# Make grid of all predictions. THese are depths from 0 - the max integer depth, by quarters, plus depths at which enzymes were measured
pred.depths <- seq(from=0, to=ceiling(max(cell_d$depth.m)), by=0.25)
pred.depths <- c(pred.depths, enz_depths$depth.m)
pred_depths <- data.frame(depth.m = pred.depths) # predict.loess takes a data frame

# Make predictions from the two loess models
pred.shallow <- predict(m_shallow, newdata = pred_depths)
pred_shallow <- data.frame(depth.m = pred_depths$depth.m, pred.cells = pred.shallow)
pred_shallow$section <- "shallow"
pred.deep <- predict(m_deep, newdata = pred_depths)
pred_deep <- data.frame(depth.m = pred_depths$depth.m, pred.cells = pred.deep)
pred_deep$section <- "deep"
pred_cells <- na.omit(rbind(pred_shallow, pred_deep)) %>%
  mutate(core = "predicted")

# Pull out predictions for just the enzyme-containing depths
enz_depth_preds <- pred_cells[pred_cells$depth.m %in% enz_depths$depth.m, ]

# Plot the data
p_cells <- ggplot() + 
  geom_point(data=enz_depth_preds, aes(x=depth.m, y=pred.cells, colour = "interpolated")) +
  geom_point(data=cell_d, aes(x=depth.m, y=cells, colour=core)) + # actual cell counts
  geom_line(data = pred_cells %>% filter(section == "shallow"), aes(x=depth.m, y=pred.cells), colour = "gray75") +
  geom_line(data = pred_cells %>% filter(section == "deep"), aes(x=depth.m, y=pred.cells), colour = "gray75") +
  #geom_smooth(data = cell_d, method="lm", aes(x=depth.m, y=cells, group = 1)) + 
  geom_vline(xintercept = 51.68, colour = "gray50")+
  scale_x_reverse() +
  #scale_colour_brewer(palette="Dark2") +
  xlab("depth, mbsf") +
  ylab(expression(paste("cells ", ml^{-1}, " sed"))) +
  scale_colour_manual(name = "core", values = c("#517C96", "#E65933",
                                 "#000000")) + 
  coord_flip() + 
  theme(legend.position=c(0.81, 0.15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.background = element_blank())
if(print.plots) {
  print(p_cells)
}
if(save.plots) {
  ggsave("plots/cell_counts_smoothed.png", p_cells, height=4, width=3.35, units="in", dpi=300)
}

#ggsave("plots/2016_04_13_cell_smoothing.png", p_cells, height=4, width=3, units="in", dpi=300)



# Export the predictions
#write_csv(enz_depth_preds, "data/2016_04_16_cell_smoothing.csv")

#####
# Merge the cell data into your enzyme activity data
####

#source("R/2016_04_20_LB_Vmax_analysis.R")# Uh-oh

# Note: This script makes predictions (m_shallow and m_deep), saves the predictions to a .csv, and then re-loads them. Deal with it.gif
cell_pred <- read.csv("data/2016_04_16_cell_smoothing.csv")

slopes_and_cells <- merge(samp_slopes, cell_pred, by.x="depth.mbsf", by.y="depth.m")
slopes_and_cells$cell.sp.v0 <- slopes_and_cells$v0 / slopes_and_cells$pred.cells
slopes_and_cells$cell.sp.v0.se <- slopes_and_cells$v0.se / slopes_and_cells$pred.cells

# Change units of cell specific v0
#   First multiply by density (g / cm^3) b/c cell counts are per cm and v0 is per g sed
sed.density <- 1.8 # g/cm^3, I TOTALLY MADE THAT UP
slopes_and_cells$cell.sp.v0 <- slopes_and_cells$cell.sp.v0 * sed.density
slopes_and_cells$cell.sp.v0.se <- slopes_and_cells$cell.sp.v0.se * sed.density

# Now go from umol to amol (umol is 10^-6, fmol is 10^-15, difference is 10^9)
slopes_and_cells$cell.sp.v0 <- slopes_and_cells$cell.sp.v0 * 1e12
slopes_and_cells$cell.sp.v0.se <- slopes_and_cells$cell.sp.v0.se * 1e12
attr(slopes_and_cells$cell.sp.v0, "units") <- "femtomole per cell per hour"

# Copied code from v0.R, but replaced the y, ymin and ymax arguments. 
# Remember kids: if you like writing code, then writing non-extensible code gives you an excuse to write more code!
# draw_depth_plot <- function(df, colour = "black", x.title = FALSE, y.title = FALSE, legend = FALSE) {
#   # First assign axis limits depending on which enzyme class we're talking about
#   get_ymax <- function(x) {
#     if(length(unique(x)) !=1) {
#       ymax <- 80
#     } else {
#       ymax <- switch (unique(x)[1],
#                       "peptidase" = 175,
#                       "glycosylase" = 75,
#                       "phosphatase" = 700
#       )
#     }
#     ymax
#   }
#   max.y <- get_ymax(df$class) # Get the limit for the subpanel
#   
#   p <- ggplot(df, aes(x=depth.mbsf, y=cell.sp.v0, linetype=treatment, colour = class)) + 
#     geom_pointrange(aes(ymin=cell.sp.v0-cell.sp.v0.se, ymax=cell.sp.v0+cell.sp.v0.se, shape = treatment), size = 0.25) +
#     geom_line(size = 0.25) +
#     geom_vline(xintercept = 51, colour = "gray50") + 
#     scale_x_reverse() + 
#     scale_linetype_manual(values=c("live"="solid","killed"="dashed")) +
#     scale_colour_manual(values = c("peptidase" = "#1b9e77", "glycosylase" = "#d95f02", "phosphatase" = "#7570b3")) +
#     scale_shape_manual(values = c(1, 19)) + 
#     expand_limits(xmin=0) +
#     ylim(c(-11, max.y)) + 
#     ylab(expression(paste("cell specific ", v[0], ", amol substrate ", cell^{-1}, " ", hr^{-1})))+ 
#     xlab("depth, mbsf") + 
#     coord_flip() +
#     facet_wrap(~enzyme, nrow=1) +
#     theme(axis.text.x  = element_text(angle=-45, hjust=0),
#           text = element_text(size = 8))
#   if(!legend) {
#     p <- p + theme(legend.position = "none")
#   }
#   if(!x.title) {
#     p <- p + theme(axis.title.x = element_blank())
#   }
#   if(!y.title) {
#     p <- p + theme(axis.title.y = element_blank())
#   }
#   
#   p
# }
# 
# Get plot legend
p_legend <- cowplot::get_legend(draw_depth_plot(slopes_and_cells, legend = TRUE, cell.sp = TRUE))

unique.enzymes <- unique(slopes_and_cells$enzyme)
plot_list <- list()
for(i in unique.enzymes) {
  single_enz <- slopes_and_cells %>% filter(enzyme == i)
  plot_list[[i]] <- draw_depth_plot(single_enz)
}

v0_fig <- cowplot::plot_grid(plot_list[["clostripain"]], plot_list[["gingipain"]], plot_list[["arginyl AP"]], plot_list[["leucyl AP"]], plot_list[["prolyl AP"]], plot_list[["ornithyl AP"]], NULL,
                             plot_list[["alpha-\nglucosidase"]], plot_list[["beta-\nglucosidase"]], plot_list[["cellobiosidase\n"]], plot_list[["beta-\nxylosidase"]], plot_list[["N-acetyl-\nglucosaminidase"]], plot_list[["alkaline\nphosphatase"]], p_legend,
                             ncol = 7)

x.grob <- grid::textGrob(expression(paste(v[0], ", nmol substrate g " , sed^{-1}, " ", hr^{-1})), 
                         gp=grid::gpar(fontface="bold", col="black", fontsize=10))
y.grob <- grid::textGrob("depth, mbsf", 
                         gp=grid::gpar(col="black", fontsize=10), rot = 90)
p_v0_final <- gridExtra::grid.arrange(gridExtra::arrangeGrob(v0_fig, left = y.grob, bottom = x.grob))

if(print.plots) {
  print(p_v0_final)
}
if(save.plots) {
  ggsave("plots/v0_downcore_cell_specific.png", p_v0_final, height = 4, width = 7.08, units = "in", dpi = 300)
}

# Clean out the global environment
unused <- NA
unused <- ls()[!ls() %in% c("print.plots", "print.extra.plots", "save.plots", "samp_slopes", "lm_stats")] # Run this twice, believe it or not, to get the 
rm(list=unused)

