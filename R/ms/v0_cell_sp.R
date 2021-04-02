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
pred_cells <- na.omit(rbind(pred_shallow, pred_deep))

# Pull out predictions for just the enzyme-containing depths
enz_depth_preds <- pred_cells[pred_cells$depth.m %in% enz_depths$depth.m, ]

# Plot the data
p_cells <- ggplot() + 
  geom_point(data=cell_d, aes(x=depth.m, y=cells, colour=core)) + # actual cell counts
  geom_point(data=enz_depth_preds, aes(x=depth.m, y=pred.cells)) +
  scale_x_reverse() +
  scale_colour_brewer(palette="Dark2") +
  xlab("depth, mbsf") +
  ylab(expression(paste("cells ", ml^{-1}, " sed"))) +
  coord_flip() + 
  theme(legend.position="top")
print(p_cells)
#ggsave("plots/2016_04_13_cell_smoothing.png", p_cells, height=4, width=3, units="in", dpi=300)



# Export the predictions
#write_csv(enz_depth_preds, "data/2016_04_16_cell_smoothing.csv")

#####
# Merge the cell data into your enzyme activity data
####

#source("R/2016_04_20_LB_Vmax_analysis.R")# Uh-oh

# Say your data frame of enzyme activities is called 'slopes', and the depth column is called 'depth'
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


# Reordering the factors
slopes_and_cells$substrate <- factor(slopes_and_cells$substrate,
                                levels = c("phe-arg-amc", "phe-val-arg-amc", "arg-amc", "leu-amc", 
                                           "pro-amc", "orn-amc", "mub-b-glu", "mub-a-glu", 
                                           "mub-cell", "mub-xylo", "mub-nag", "mub-po4"))
# Plot the data
p_slopes_and_cells <- ggplot(slopes_and_cells, aes(x=depth.mbsf, y=cell.sp.v0, colour=substrate, linetype=treatment)) + 
  geom_pointrange(aes(ymin=cell.sp.v0-cell.sp.v0.se, ymax=cell.sp.v0+cell.sp.v0.se)) +
  geom_line() +
  scale_x_reverse() + 
  scale_linetype_manual(values=c("live"="solid","killed"="dashed")) +
  scale_colour_manual(values=c("phe-arg-amc"="dark green", "phe-val-arg-amc"="dark green", 
                               "arg-amc"="dark green", "leu-amc"="dark green", "pro-amc"="dark green", 
                               "orn-amc"="dark green", "mub-b-glu"="red", "mub-a-glu"="red", 
                               "mub-cell"="red", "mub-xylo"="red", "mub-nag"="red", 
                               "mub-po4"="blue")) +
  expand_limits(xmin=0) +
  ylab(expression(paste("cell specific ", v[0], ", amol ", "substrate ", cell^{-1}, " ", hr^{-1} )))+ 
  xlab("depth, mbsf") + 
  coord_flip() +
  facet_wrap(~substrate, nrow=4) +
  #theme_set(theme_bw()) +
  #    theme(text=element_text(size=70)) +
  theme(axis.text.x  = element_text(angle=-45, hjust=0)) +
  theme(legend.position="none")
print(p_slopes_and_cells)
ggsave("plots/2016_04_26_cell.sp.V0.png", p_slopes_and_cells, height=8, width=5, units="in", dpi=600)

