#Analysing Vmax assays from Lille Baelt Site
#1-28-16

# Load packages
require(plyr)
require(lubridate)
require(reshape2)
# require(ggplot2)
source("R/ms/draw_depth_plot.R")




###########
# Load the data
############

#Load data for each depth and add a colunm for depth 
#d4 <- read.csv("data/2015_10_06_LB_4.5mbsf.csv")
#d4$depth.mbsf <- 4.5 

#Redo 4.5mbsf 1-26-16
d4_redo <- read.csv("data/2016_01_25_LB_4.5mbsf.csv")

d11 <- read.csv("data/2015_10_19_LB_11.1mbsf.csv")
d11$depth.mbsf <- 11.1 

d17 <- read.csv("data/2015_10_14_LB_17.6mbsf.csv")
d17$depth.mbsf <- 17.6

d24 <-read.csv("data/2015_10_21_LB_24.3mbsf.csv")
d24$depth.mbsf <- 24.3

d30 <-read.csv("data/2015_12_07_LB_30.9mbsf.csv")
d30$depth.mbsf <-30.9

#d37 <- read.csv("data/2015_10_13_LB_37.5mbsf.csv")
#d37$depth.mbsf <- 37.5

d37_redo <- read.csv("data/2016_01_27_LB_37.5mbsf.csv")

d48 <- read.csv("data/2015_12_30_LB_48.22mbsf.csv")

d54 <- read.csv("data/2015_12_10_LB_54.95mbsf.csv")

d61 <-read.csv("data/2015_12_15_LB_61.55mbsf.csv")

d43 <-read.csv("data/2015_12_17_LB_43.15mbsf.csv")

d68 <-read.csv("data/2015_12_22_LB_68.12mbsf.csv")

d77 <-read.csv("data/2015_12_30_LB_77.92mbsf.csv")

#Put all the data into a list
all_data_list <- list(d4_redo=d4_redo, d11=d11, d17=d17, d24=d24, 
                      d30=d30, d37_redo=d37_redo, d43=d43, d48=d48, 
                      d54=d54, d61=d61, d68=d68, d77=d77)

#Put all of the data into a dataframe
all_df <- ldply(all_data_list, identity)

# Calculate elapsed times
some_day <- "2016_01_28"
all_df$Rtime <- ymd_hms(paste(some_day, all_df$elapsed.time))
all_df$elapsed <- as.numeric(all_df$Rtime - min(all_df$Rtime, na.rm=TRUE))/3600
attr(all_df$elapsed, "units") <- "hours"

#Split out sample_data from calibration data
data_allsamp <- subset(all_df, sample.calib=="sample")

########
# Clean out the data
########

# Make raw plots of each substrate and depth
rawplot <- function(x, base.path="plots/small_raw_plots/", print.plot=FALSE, save.plot=TRUE) {
  p <- ggplot(x, aes(x=elapsed, y=RFU, colour=treatment)) +
    geom_point() + 
    geom_smooth(method="lm", se=FALSE) + 
    facet_grid(depth.mbsf ~ substrate, scales="free")
  
  substrate <- x$substrate[1]
  depth <- x$depth.mbsf[1]
  
  fn <- paste0(base.path, substrate,"_", depth, ".png") 
  
  if(print.plot) {
    print(p)
  }
  if(save.plot) { # Omit type="cairo" when you run this on a pc
    ggsave(fn, p, height=3, width=4, units="in", dpi=100)
  }
}

# actually make the small plots
#all_p_list <- d_ply(data_allsamp, c("substrate", "depth.mbsf"), rawplot)

######################
#Edit to omit obvious outliers

# substrate=mub-cell, depth=17.6, treatment==killed, elapsed > 
# subset(data_allsamp, substrate=="mub-cell" & depth.mbsf==17.6 & RFU==38155.35)
d_edit <- subset(data_allsamp, !(substrate=="mub-cell" & depth.mbsf==17.6 & treatment=="killed" & RFU==38155.35))

# substrate=mub=nag, treatment=="killed", RFU > 30000
d_edit <- subset(d_edit, !(substrate=="mub-nag" & depth.mbsf==17.6 & treatment=="killed" & RFU>30000))

#bad live points
d_edit <- subset(d_edit, !(substrate=="phe-val-arg-amc" & depth.mbsf==24.3 & treatment=="live" & RFU==282.87))
#d_edit <- subset(d_edit, !(substrate=="mub-po4" & depth.mbsf==37.5 & treatment=="live" & RFU==1180.24))
#d_edit <- subset(d_edit, !(substrate=="mub-po4" & depth.mbsf==37.5 & treatment=="live" & RFU==1301.44))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==24.3 & treatment=="live" & RFU==1590.29))
#d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==37.5 & treatment=="live" & RFU==39952.89))
#d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==37.5 & treatment=="live" & RFU==39460.46))
d_edit <- subset(d_edit, !(substrate=="mub-po4" & depth.mbsf==4.5 & treatment=="live" & RFU==1747.52))
d_edit <- subset(d_edit, !(substrate=="mub-po4" & depth.mbsf==4.5 & treatment=="live" & RFU==1835))
d_edit <- subset(d_edit, !(substrate=="mub-po4" & depth.mbsf==11.1 & treatment=="live" & RFU==3185.72))
#d_edit <- subset(d_edit, !(substrate=="mub-cell" & depth.mbsf==37.5 & treatment=="live" & RFU==26949.48))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==4.5 & treatment=="live" & RFU==26398.35))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==4.5 & treatment=="live" & RFU==81515.35))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==4.5 & treatment=="live" & RFU==90679.97))
d_edit <- subset(d_edit, !(substrate=="mub-po4" & depth.mbsf==4.5 & treatment=="live" & RFU==614.58))
d_edit <- subset(d_edit, !(substrate=="mub-a-glu" & depth.mbsf==24.3 & treatment=="live" & rep==2))
d_edit <- subset(d_edit, !(substrate=="mub-a-glu" & depth.mbsf==24.3 & treatment=="live" & rep==3))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==4.5 & treatment=="live" & rep==2))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==4.5 & treatment=="live" & rep==1))

#bad killed points
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==24.3 & treatment=="killed" & RFU==750.21))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==11.1 & treatment=="killed" & RFU==10606.48))
d_edit <- subset(d_edit, !(substrate=="mub-po4" & depth.mbsf==17.6 & treatment=="killed" & RFU==977.47))
d_edit <- subset(d_edit, !(substrate=="mub-nag" & depth.mbsf==17.6 & treatment=="killed" & RFU==34521.53))
d_edit <- subset(d_edit, !(substrate=="mub-nag" & depth.mbsf==24.3 & treatment=="killed" & RFU==38.03))
d_edit <- subset(d_edit, !(substrate=="mub-cell" & depth.mbsf==17.6 & treatment=="killed" & RFU==38155.35))
d_edit <- subset(d_edit, !(substrate=="leu-amc" & depth.mbsf==11.1 & treatment=="killed" & RFU==42103.25))
d_edit <- subset(d_edit, !(substrate=="arg-amc" & depth.mbsf==24.3 & treatment=="killed" & RFU==2264.03))
d_edit <- subset(d_edit, !(substrate=="mub-xylo" & depth.mbsf==37.5 & treatment=="killed" & RFU==226.89))
d_edit <- subset(d_edit, !(substrate=="orn-amc" & depth.mbsf==37.5 & treatment=="killed" & RFU==1260.77))

################
# Plot raw data
p_raw_e <- ggplot(d_edit, aes(x=elapsed, y=RFU, colour=treatment)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) + 
  facet_grid(depth.mbsf ~ substrate, scales="free") 

if(print.extra.plots) {
  print(p_raw_e)
}
if(save.plots) {
  ggsave("plots/2016_03_01_LB_raw_data_allsamp.png", p_raw_e, height=8, width=12, units="in", dpi=300)
}
#
d_edit_live <-subset(d_edit, treatment=="live")

p_raw_live <- ggplot(d_edit_live, aes(x=elapsed, y=RFU)) + 
  geom_smooth(method="lm", se=FALSE, size=1) + 
  geom_point() + 
  theme_set(theme_bw()) +
  facet_grid(depth.mbsf ~ substrate, 
             scales="free") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) 
if(print.extra.plots) {
  print(p_raw_live)
}
if(save.plots) {
  ggsave("plots/2016_01_28_LB_raw_data_livesamp.png", p_raw_live, height=20, width=40, units="in", dpi=150)
}
#
# Calculate slopes of fluorescence as a function of time
samp_slopes <- ddply(d_edit, c("substrate", "treatment", "depth.mbsf", "conc.uM"), lm_stats, "elapsed", "RFU")

###########
#Calibrate
###########

calib_data <- subset(all_df, sample.calib=="calib") 
p_calib <- ggplot(calib_data, aes(x=conc.uM, y=RFU, colour=as.factor(depth.mbsf))) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  facet_wrap(~substrate, scales="free_y")
if(print.extra.plots) {
  print(p_calib)
}
#

#Plot to see individual depth plots 
calib_37.5 <- subset(calib_data, depth.mbsf=="37.5")
p_calib_37.5 <-ggplot(calib_37.5, aes(x=conc.uM, y=RFU))+
  geom_point()+
  geom_smooth(method="lm") +
  facet_wrap(~substrate, scales="free_y")
if(print.extra.plots) {
  print(p_calib_37.5) # THe weird thing is that this looks like 3 separate calibration curves (for MUB, at least)
}
#

# Get rid of the bad datapoints
calib_data_e <- subset(calib_data, !(conc.uM==8 & substrate=="amc-std" & depth.mbsf==17.6))
#calib_data_e <- subset(calib_data_e, !(conc.uM==32 & substrate=="amc-std" & depth.mbsf==37.5))
calib_data_e <- subset(calib_data_e, !(conc.uM==80 & substrate=="amc-std" & depth.mbsf==37.5))
calib_data_e <- subset(calib_data_e, !(conc.uM==80 & substrate=="mub-std" & depth.mbsf==37.5))

#edit to take out bad 1st reads with mub-std
calib_data_e <- subset(calib_data_e, !(depth.mbsf=="30.9" & substrate=="mub-std" & elapsed==0.00000))
calib_data_e <- subset(calib_data_e, !(depth.mbsf=="54.95" & substrate=="mub-std" & elapsed==0.00000))
calib_data_e <- subset(calib_data_e, !(depth.mbsf==61.55 & substrate=="mub-std" & elapsed==0.00000))
calib_data_e <- subset(calib_data_e, !(depth.mbsf==43.15 & substrate=="mub-std" & elapsed==0.00000))
calib_data_e <- subset(calib_data_e, !(depth.mbsf==4.5 & substrate=="mub-std" & elapsed==0.00000))
calib_data_e <- subset(calib_data_e, !(depth.mbsf==48.22 & substrate=="mub-std" & elapsed==0.00000))
calib_data_e <- subset(calib_data_e, !(depth.mbsf==68.12 & substrate=="mub-std" & elapsed==0.00000))
calib_data_e <- subset(calib_data_e, !(depth.mbsf==37.5 & substrate=="mub-std" & elapsed==0.00000))

# Check edited calibration data
p_calib <- ggplot(calib_data, aes(x=conc.uM, y=RFU, colour=as.factor(depth.mbsf))) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  # facet_wrap(~substrate, scales="free_y") +
  facet_grid(depth.mbsf~substrate, scales="free_y") +
  #theme_set(theme_bw()) +
  # theme(text=element_text(size=30)) +
  theme(text=element_text(size=10),
        axis.title.x = element_text(face="bold"),
        axis.text.x  = element_text(angle=90, vjust=0.5))
if(print.extra.plots) {
  print(p_calib)
}

###################################
## trying to figure out what's up with calibration curve slopes
###############################

# determine slope of calibration curve for each depth and incubation time and standard

calib_slope_ts <- ddply(calib_data, c("substrate", "depth.mbsf", "elapsed"), 
                        lm_stats, xvar="conc.uM", yvar="RFU")
calib_slope_ts$depth.mbsf <- as.factor(calib_slope_ts$depth.mbsf)

p_calib_slope_ts <-ggplot(calib_slope_ts, aes(x=elapsed, y=slope, colour=depth.mbsf)) + 
  geom_pointrange(aes(ymax=slope+slope.se, ymin=slope-slope.se)) +
  geom_line() + 
  facet_wrap(~substrate)
if(print.extra.plots) {
  print(p_calib_slope_ts)
}

# ggsave("plots/2016_04_20_calib_slope_ts.png", p_calib_slope_ts, height=4, width=6, units="in", dpi=400)

############################################
############################################

# Get calibration slopes (& errors, etc) for each depth and fluorophore
calib_slopes <- ddply(calib_data_e, c("substrate", "depth.mbsf"), lm_stats, "conc.uM", "RFU")

#### Actually calibrate the slopes
# Create empty column for 'calibrate slope'
samp_slopes$calib.slope <- NA

# Define a vector of all the -AMC substrates
AMC.subs <- c("arg-amc", "leu-amc", "orn-amc", "phe-arg-amc","phe-val-arg-amc", "pro-amc")
MUB.subs <- c("mub-a-glu", "mub-b-glu", "mub-cell", "mub-nag", "mub-po4", "mub-xylo")

# Motherfucker, I somehow got rid of the part where I actually create the calib slopes
source("R/ms/clean_v0_calib_units.R")

# Get uM fluorophore per hour by dividing slope of RFU per hour ("slope") by calibration slope ("calib.slope")
samp_slopes$v0.bulk <- samp_slopes$slope / samp_slopes$calib.slope * 1000 # switching from uM to nM units
samp_slopes$v0.bulk.se <- samp_slopes$slope.se / samp_slopes$calib.slope * 1000


# Make a small data frame of sample depths and g sediment per g buffer
mass_df <- data.frame(depth.mbsf=c(4.5, 11.1, 17.6, 24.3, 30.9, 37.5, 43.15, 48.22, 54.95, 61.55, 68.12, 77.92), sed.conc=c(30.2, 30.2, 30.2, 30.2, 30.2, 30.2, 30.2, 30.2, 30.2, 30.2, 30.2, 30.2))

# Merge mass concentrations into big data frame of slopes
samp_slopes <- merge(samp_slopes, mass_df, by="depth.mbsf")

# Then divide by grams of sediment per liter of buffer to get umol fluorophore per g sediemnt per hour
samp_slopes$v0 <- samp_slopes$v0.bulk / samp_slopes$sed.conc
samp_slopes$v0.se <- samp_slopes$v0.bulk.se / samp_slopes$sed.conc

###########
#Plot Down Core Profile of Sample Slopes
###########
# Reordering the factors
samp_slopes$substrate <- factor(samp_slopes$substrate,
                            levels = c("phe-arg-amc", "phe-val-arg-amc", "arg-amc", "leu-amc", 
                                       "pro-amc", "orn-amc", "mub-b-glu", "mub-a-glu", 
                                       "mub-cell", "mub-xylo", "mub-nag", "mub-po4"))

p_samp_slopes <- ggplot(samp_slopes, aes(x=depth.mbsf, y=v0, colour=substrate, linetype=treatment)) + 
  geom_pointrange(aes(ymin=v0-v0.se, ymax=v0+v0.se)) +
  geom_line() +
  scale_x_reverse() + 
  scale_linetype_manual(values=c("live"="solid","killed"="dashed")) +
  scale_colour_manual(values=c("phe-arg-amc"="dark green", "phe-val-arg-amc"="dark green", 
                               "arg-amc"="dark green", "leu-amc"="dark green", "pro-amc"="dark green", 
                               "orn-amc"="dark green", "mub-b-glu"="red", "mub-a-glu"="red", 
                               "mub-cell"="red", "mub-xylo"="red", "mub-nag"="red", 
                               "mub-po4"="blue")) +
  expand_limits(xmin=0) +
  ylab(expression(paste(v[0], ", ", mu, "mol ", "substrate ", g^{-1}, " sed ", hr^{-1})))+ 
  xlab("depth, mbsf") + 
  coord_flip() +
  facet_wrap(~substrate, nrow=4)+
  #theme_set(theme_bw()) +
  #theme(text=element_text(size=60)) +
  theme(axis.text.x  = element_text(angle=-45, hjust=0))+
  theme(legend.position="none")

if(print.plots) {
  print(p_samp_slopes)
}


# Want to write a better approach to make nice plots
# First, replace substrate names with enzyme names
enzyme_names <- function(chr) {
 case_when(chr == "phe-arg-amc" ~ "gingipain",
           chr == "phe-val-arg-amc" ~ "clostripain",
           chr == "arg-amc" ~ "arginyl AP",
           chr == "leu-amc" ~ "leucyl AP",
           chr == "pro-amc" ~ "prolyl AP",
           chr == "orn-amc" ~ "ornithyl AP",
           chr == "mub-b-glu" ~ "beta-\nglucosidase",
           chr == "mub-a-glu" ~ "alpha-\nglucosidase",
           chr == "mub-xylo" ~ "beta-\nxylosidase",
           chr == "mub-nag" ~ "N-acetyl-\nglucosaminidase",
           chr == "mub-cell" ~ "cellobiosidase\n",
           chr == "mub-po4" ~ "alkaline\nphosphatase") 
}

# Classifies enzymes as peptidase, glycosylase, or phosphatase
enzyme_type <- function(chr) {
  type <- rep(NA, length(chr))
  type[str_which(chr, pattern = "amc")] <- "peptidase"
  type[str_detect(chr, pattern = "mub") & !(chr == "mub-po4")] <- "glycosylase" # The fuckup: I'm mixing numeric w/ boolean indexing
  type[chr == "mub-po4"] <- "phosphatase"
  type
}

samp_slopes <- samp_slopes %>%
  mutate(enzyme = enzyme_names(substrate),
         class = enzyme_type(substrate))


# Get a legend for the grouped plots. I won't use the plot, just the legend
p_legend <- cowplot::get_legend(draw_depth_plot(samp_slopes, legend = TRUE))

# Huh: what if I make each plot separately
unique.enzymes <- unique(samp_slopes$enzyme)
plot_list <- list()
for(i in unique.enzymes) {
  single_enz <- samp_slopes %>% filter(enzyme == i)
  plot_list[[i]] <- draw_depth_plot(single_enz)
}

v0_fig <- cowplot::plot_grid(plot_list[["clostripain"]], plot_list[["gingipain"]], plot_list[["arginyl AP"]], plot_list[["leucyl AP"]], plot_list[["prolyl AP"]], plot_list[["ornithyl AP"]], NULL,
                             plot_list[["alpha-\nglucosidase"]], plot_list[["beta-\nglucosidase"]], plot_list[["cellobiosidase\n"]], plot_list[["beta-\nxylosidase"]], plot_list[["N-acetyl-\nglucosaminidase"]], plot_list[["alkaline\nphosphatase"]], p_legend,
                             ncol = 7, 
                             rel_widths = c(1.15, rep(1, 5), 1.3))

x.grob <- grid::textGrob(expression(paste(v[0], ", nmol substrate g " , sed^{-1}, " ", hr^{-1})), 
                  gp=grid::gpar(fontface="bold", col="black", fontsize=10))
y.grob <- grid::textGrob("depth, mbsf", 
                         gp=grid::gpar(col="black", fontsize=10), rot = 90)
p_v0_final <- gridExtra::grid.arrange(gridExtra::arrangeGrob(v0_fig, left = y.grob, bottom = x.grob))

if(print.plots) {
  print(p_v0_final)
}
if(save.plots) {
  ggsave("plots/v0_downcore.tiff", p_v0_final, height = 4, width = 7.08, units = "in", dpi = 300, compression = "lzw")
}

#######
# Delete unused objects
#######

unused <- NA
unused <- ls()[!ls() %in% c("print.plots", "print.extra.plots", "save.plots", "samp_slopes", "lm_stats", "draw_depth_plot")] # Run this twice, believe it or not, to get the 
rm(list=unused)
