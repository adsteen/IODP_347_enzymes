#####################
# Km analysis on clostripain 
######################

# Load packages
library(plyr) # split-apply-combine
library(ggplot2) # plotting
library(reshape2) # Change shape of datasets: I don't think we use this in this script
library(lubridate) # Handle dates

source("R/Jennas_R/lm_stats.R")

# Set graphical theme
theme_set(theme_bw() + 
            theme(text=element_text(size=9)))

#Load the data
d4 <- read.csv("data/2016_02_24_km_4.5mbsf.csv")
d11 <- read.csv("data/2016_02_17_km_11.1mbsf.csv")
d17 <- read.csv("data/2016_02_23_km_17.6mbsf.csv")
d24 <- read.csv("data/2016_02_25_km_24.3mbsf.csv")
#d37 <- read.csv("data/2016_02_27_km_37.5mbsf.csv") #has only two times: 0 hr and 48 hrs. 
d43 <- read.csv("data/2016_02_26_km_43.15mbsf.csv")

# Put all the data into a list
#all_data_list <- list(d4=d4, d11=d11, d17=d17, d43=d43)
all_data_list <- list(d4=d4, d11=d11, d17=d17, d24=d24, d43=d43)

#d37 needs to be processed separately - see R/km_d37.R

# Put all the data into a dataframe
all_df <- ldply(all_data_list, identity)

# Calculate elapsed times
some_day <- "2016_02_24"
all_df$Rtime <- ymd_hm(paste(some_day, all_df$time))
all_df$elapsed <- as.numeric(all_df$Rtime - min(all_df$Rtime, na.rm=TRUE))/3600
attr(all_df$elapsed, "units") <- "hours"

# Add in d37, for which elapsed time was calculated differently because the data were formatted differently
source("R/preprocess_d37_sat_data.R")
all_df <- all_df %>% 

# Calculate slopes of fluorescence as a function of time
slopes <- ddply(all_df, c("depth.mbsf", "treatment", "conc.uM", "sample.calib"), lm_stats, "elapsed", "RFU")

# make a data frame without the killed controls
slopes_live <- subset(slopes, treatment=="live")

#make a data frame with samples only
samp_slopes <- subset(slopes, sample.calib=="sample")

# rename all_df for consistent capitalization
#all_df <- rename(all_df, c("FSU" = "fsu"))

# Plot raw data, look for outliers
p_raw <- ggplot(all_df, aes(x=elapsed, y=RFU, colour=as.factor(conc.uM), linetype = treatment)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~.id, scales="free")
print(p_raw)

# Plot slopes, see how they look
p_uncal_samp_slopes <- ggplot(samp_slopes, aes(x=conc.uM, y=slope, colour=treatment)) +
  geom_pointrange(aes(ymin=slope-slope.se, ymax=slope+slope.se))
#theme(axis.text=element_text(size=16),
      #axis.title=element_text(size=20))
print(p_uncal_samp_slopes)
#ggsave("plots/2016_0_29_Km_4.5_uncal.png", p_uncal_samp_slopes, height=5, width=8, units="in", dpi=400)

############## 
# Calibrate the data
##############

calib_data <- subset(all_df, sample.calib=="calib") # What is going on 
p_calib <- ggplot(calib_data, aes(x=conc.uM, y=RFU)) + 
  geom_point() + 
  geom_smooth(method="lm") 
  #facet_wrap(~substrate, scales="free_y")
print(p_calib)
#ggsave("plots/2016_02_29_km_4.5_calib.png", p_calib, height = 5, width = 8, units = "in", dpi = 400)

# Get the slope of the calibration curve
calib.slope <- summary(lm(RFU~conc.uM, data=calib_data))$coefficients[2, 1]
attr(calib.slope, "units") <- "RFU per uM AMC"

# Divide uncalibrated slope of fluorescence vs time by the calib slope to get Vmax in umol AMC per liter per hour
samp_slopes$v0.raw <- samp_slopes$slope / calib.slope
samp_slopes$v0.raw.se <- samp_slopes$slope.se / calib.slope

#########
# Calculate as umol AMC / hr / g sed
#########

conc.sed <- 3/0.1 # 3 g / 0.1 L
attr(conc.sed, "units") <- "g wet sed / L"

samp_slopes$v0 <- samp_slopes$v0.raw / conc.sed
samp_slopes$v0.se <- samp_slopes$v0.raw.se / conc.sed
attr(samp_slopes$v0, "units") <- "umol AMC / hr / g wet sed"

p_samp_slopes <- ggplot(samp_slopes, aes(x=conc.uM, y=v0, colour=treatment)) + 
  geom_pointrange(aes(ymin=v0-v0.se, ymax=v0+v0.se)) +
#  geom_line(size=2) +
#  expand_limits(xmin=0) +
  ylab(expression(paste(v[0], ", ", mu, "mol ", g^{-1}, " sed ", hr^{-1})))+ 
  xlab("substrate conc.,uM")
print(p_samp_slopes)
#ggsave("plots/2016_03_08_Km_all_depths.png", p_samp_slopes, height=5, width=8, units="in", dpi=150)

samp_slopes_4 <- subset(samp_slopes, depth.mbsf=="4.5")
samp_slopes_11 <- subset(samp_slopes, depth.mbsf=="11.1")
samp_slopes_17 <- subset(samp_slopes, depth.mbsf=="17.6")
#samp_slopes_24 <- subset(samp_slopes, depth.mbsf=="24.3")
samp_slopes_43 <- subset(samp_slopes, depth.mbsf=="43.15")

#########
# Fit MM curves
#########

# Define the michaelis menten formula
mm_form <- formula(v0 ~ (Vmax*conc.uM)/(Km+conc.uM))

# Calculate nls models for each substrate individually
nls_4 <- nls(mm_form, subset(samp_slopes_4, treatment=="live"), start=list(Vmax=0.01, Km=100))
nls_11 <- nls(mm_form, subset(samp_slopes_11, treatment=="live"), start=list(Vmax=0.01, Km=100))
nls_17 <- nls(mm_form, subset(samp_slopes_17, treatment=="live"), start=list(Vmax=0.01, Km=100))
#nls_24 <- nls(mm_form, subset(samp_slopes_24, treatment=="live"), start=list(Vmax=0.01, Km=100))
nls_43 <- nls(mm_form, subset(samp_slopes_43, treatment=="live"), start=list(Vmax=0.01, Km=100))
#nls_depth <- nls(mm_form, subset(slopes, & treatment=="live"), start=list(Vmax=0.01, Km=100))
# continue with all the depths

# Calculate a data frame of predicted values
pred_grid <- data.frame(conc.uM=seq(0, max(samp_slopes$conc.uM), length.out=100))

# Calculate predictions for each nls model separately
preds.4 <- predict(nls_4, newdata=pred_grid) 
preds_4 <- data.frame(conc.uM=pred_grid$conc.uM, 
                        v0=preds.4)
preds.11 <- predict(nls_11, newdata=pred_grid) 
preds_11 <- data.frame(conc.uM=pred_grid$conc.uM, 
                      v0=preds.11)
preds.17 <- predict(nls_17, newdata=pred_grid) 
preds_17 <- data.frame(conc.uM=pred_grid$conc.uM, 
                      v0=preds.17)
#preds.24 <- predict(nls_24, newdata=pred_grid) 
#preds_24 <- data.frame(conc.uM=pred_grid$conc.uM, 
   #                   v0=preds.24)
preds.43 <- predict(nls_43, newdata=pred_grid) 
preds_43 <- data.frame(conc.uM=pred_grid$conc.uM, 
                      v0=preds.43)
#leu.preds <- predict(leu_nls, newdata=pred_grid)
#leu_preds <- data.frame(substrate="Leu-AMC", conc.uM=pred_grid$conc.uM, 
#                        v0=leu.preds)

# nextsubstrate.preds <- predict(nextsubstrate_nls, newdata=pred_grid)
# nextsubstrate_preds <- data.frame(substrate="Nextsubstrate-AMC", conc.uM=pred_grid$conc.uM, 
#                       v0=nextsubstrate.preds)

# Put all the predictions together into one big data frame
preds_4$depth.mbsf <- 4.5
preds_11$depth.mbsf <- 11.1
preds_17$depth.mbsf <- 17.6
#preds_24$depth.mbsf <- 24.3
preds_43$depth.mbsf <- 43.15
preds_list <- list(d4=preds_4, d11=preds_11, d17=preds_17, d43=preds_43)
preds_df <- do.call(rbind, preds_list)
preds_df$treatment <- "live"
preds_df$depth.mbsf <- as.factor(preds_df$depth.mbsf)
samp_slopes$depth.mbsf <- as.factor(samp_slopes$depth.mbsf)

# Complete Plot
# Put all the predictions into one data frame
samp_slopes_live <- subset(samp_slopes, !(treatment=="killed"))
  #d_edit <- subset(d_edit, !(substrate=="phe-val-arg-amc" & depth.mbsf==24.3 & treatment=="live" & RFU==282.87))

# Quick, check the killed rates
p_faceted_sat_curve <- ggplot(samp_slopes, aes(x=conc.uM, y=v0, colour=treatment)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~depth.mbsf, scales="free") # yep, it's fine. Could include as supplemental I suppose.
if(print.plots) {
  print(p_faceted_sat_curve)
}
if(save.plots) {
  ggsave("plots/faceted_sat_curve.png", p_faceted_sat_curve, height = 4, width = 4, units = "in", dpi = 300)
}



p_cal_slopes <- ggplot(samp_slopes_live, aes(x=conc.uM, y=v0, colour=depth.mbsf))+
  geom_pointrange(aes(ymin=v0-v0.se, ymax=v0+v0.se)) +
  geom_line(data=preds_df) +
  scale_colour_brewer(name = "depth, mbsf", palette = "Dark2") +
  ylim(0, 0.025) +
  ylab(expression(paste(v[0], ", ", mu, "M AMC ", g^{-1}, " sed ", hr^{-1}))) +
  xlab(expression(paste("phe-val-arg-AMC conc., ", mu, "M")))
if(print.plots) {
  print(p_cal_slopes)
}
if(save.plots) {
  ggsave("plots/trypsin_sat_curves.png", height = 4, width = 3.5, units = "in", dpi = 300)
}


######
# Create a data frame of Vmax and Km from each prediction
######
get_km <- function(nls_obj) {coef(nls_obj)["Km"]}
get_km_se <- function(nls_obj) {summary(nls_obj)$coefficients["Km", "Std. Error"]}
get_vmax <- function(nls_obj) {coef(nls_obj)["Vmax"]}
get_vmax_se <- function(nls_obj) {summary(nls_obj)$coefficients["Vmax", "Std. Error"]}

depths <- as.numeric(as.character(unique(samp_slopes$depth.mbsf)))
Kms <- c(get_km(nls_4), get_km(nls_11), get_km(nls_17), get_km(nls_43))
Km_ses <- c(get_km_se(nls_4), get_km_se(nls_11), get_km_se(nls_17), get_km_se(nls_43))
Vmaxes <- c(get_vmax(nls_4), get_vmax(nls_11), get_vmax(nls_17), get_vmax(nls_43))
Vmax_ses <- c(get_vmax_se(nls_4), get_vmax_se(nls_11), get_vmax_se(nls_17), get_vmax_se(nls_43))
km_v_depth <- data.frame(depths, Kms, Km_ses, Vmaxes, Vmax_ses)

m_depth_km <- lm(Kms ~ depths, data = km_v_depth)
summary(m_depth_km)

p_km_v_depth <- ggplot(km_v_depth, aes(x=depths, y=Kms, ymin=Kms-Km_ses, ymax=Kms+Km_ses)) + 
  geom_smooth(method = "lm", color = "black") +
  geom_pointrange() + 
  expand_limits(ymin=0, xmin=0) + 
  scale_x_reverse(name = "depth, mbsf") + 
  scale_y_continuous(name = expression(paste(K[m], ", ", mu, "M AMC"))) +
  coord_flip()
if(print.plots) {
  print(p_km_v_depth)
}

p_Vmax_v_depth <- ggplot(km_v_depth, aes(x=depths, y=Vmaxes, ymin=Vmaxes-Vmax_ses, ymax=Vmaxes+Vmax_ses)) + 
  geom_smooth(method = "lm", color = "black") +
  geom_pointrange() + 
  expand_limits(ymin=0, xmin=0) + 
  scale_x_reverse(name = "depth, mbsf") + 
  scale_y_continuous(name = expression(paste(K[m], ", ", mu, "M AMC"))) +
  coord_flip()
if(print.plots) {
  print(p_Vmax_v_depth)
}

p_km_v_depth <- ggplot(km_v_depth, aes(x=depths, y=Kms, ymin=Kms-Km_ses, ymax=Kms+Km_ses)) + 
  geom_smooth(method = "lm", color = "black") +
  geom_pointrange() + 
  expand_limits(ymin=0, xmin=0) + 
  scale_x_reverse(name = "depth, mbsf") + 
  scale_y_continuous(name = expression(paste(K[m], ", ", mu, "M AMC"))) +
  coord_flip()
if(print.plots) {
  print(p_km_v_depth)
}

png("plots/sat_curves.png", height = 3, width = 6, units = "in", res = 300)
  cowplot::plot_grid(p_cal_slopes, p_km_v_depth, labels=c("A", "B"), rel_widths = c(2, 1), label_size = 10)
dev.off()
#ggsave("plots/2016_04_04_Km_four_depths_calib.png", p_cal_slopes, height=5, width=8, units="in", dpi=150)
