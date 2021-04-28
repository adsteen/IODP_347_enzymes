########
#  Calculate saturation curves for substrates other than clostripain, whcih don't seem to be in Jenna's thesis
########
library(tidyverse)

# Quick plot to look at raw data for a single substrate saturation curve
raw_plot <- function(df, substrate, treatments = FALSE) {
  p <- ggplot(df, aes(x=time, y=FSU)) + 
    facet_wrap(~conc.uM) + 
    ggtitle(substrate) # Correctly measured but garbage data
  
  if(treatments) {
    p <- p + geom_point(aes(colour = treatment))
  } else {
    p <- p + geom_point()
  }
}
mub_cell <- read_csv("data/2015_08_27_sc_mub_cell.csv")
p_mub_cell <- raw_plot(mub_cell, substrate = "MUB cell"); print(p_mub_cell) # Correctly obtained but the data are crap

mub_po4 <- read_csv("data/2015_08_27_sc_mub_po4.csv")
p_mub_po4 <- raw_plot(mub_po4, substrate = "MUB po4"); print(p_mub_po4) # Looks sort of OK?

mub_xylo <- read_csv("data/2015_08_27_sc_mub_xylo.csv")
p_mub_xylo <- raw_plot(mub_xylo, substrate = "MUB xylo"); print(p_mub_xylo) # Not very good

mub_a_glu <- read_csv("data/2015_09_08_sc_mub_a_glu.csv")
p_mub_a_glu <- raw_plot(mub_a_glu, substrate = "MUB-a-glu"); print(p_mub_a_glu) # Could work

mub_b_glu <- read_csv("data/2015_09_08_sc_mub_b_glu.csv")
p_mub_b_glu <- raw_plot(mub_b_glu, substrate = "MUB-b-glu"); print(p_mub_b_glu) # looks gorgeous TBH

mub_b_nag <- read_csv("data/2015_09_08_sc_mub_b_nag.csv")
p_mub_b_nag <- raw_plot(mub_b_nag, substrate = "MUB-b-nag"); print(p_mub_b_nag) # positively scrumptious

arg_amc_live <- read_csv("data/2015_09_17_sc_arg_amc_live.csv", skip = 3) %>% mutate(treatment = "live")
arg_amc_killed <- read_csv("data/2015_09_17_sc_arg_amc_kill.csv", skip = 3) %>% mutate(treatment = "killed")
arg_amc <- arg_amc_live %>% rbind(arg_amc_killed)
p_arg_amc <- raw_plot(arg_amc, substrate = "arg-AMC", treatments = TRUE); print(p_arg_amc) # yeah baby

leu_amc_live <- read_csv("data/2015_09_17_sc_leu_amc_live.csv", skip = 3) %>% mutate(treatment = "live") 
leu_amc_killed <- read_csv("data/2015_09_17_sc_leu_amc_kill.csv", skip = 3) %>% mutate(treatment = "killed") 
leu_amc <- leu_amc_live %>% rbind(leu_amc_killed)
p_leu_amc <- raw_plot(leu_amc, substrate = "leu-AMC", treatments = TRUE); print(p_leu_amc) # sure, why not

pro_amc_live <- read_csv("data/2015_09_17_sc_pro_amc_live.csv", skip = 3) %>% mutate(treatment = "live")
pro_amc_killed <- read_csv("data/2015_09_17_sc_pro_amc_kill.csv", skip = 3) %>% mutate(treatment = "killed")
pro_amc <- pro_amc_live %>% rbind(pro_amc_killed)
p_pro_amc <- raw_plot(pro_amc, substrate = "pro_amc", treatments = TRUE); print(p_pro_amc) # yep

orn_amc_live <- read_csv("data/2015_09_24_sc_orn_amc_live.csv", skip = 3) %>% mutate(treatment = "live")
orn_amc_killed <- read_csv("data/2015_09_24_sc_orn_amc_kill.csv", skip = 3) %>% mutate(treatment = "killed")
orn_amc <- orn_amc_live %>% rbind(orn_amc_killed)
p_orn_amc <- raw_plot(orn_amc, substrate = "orn_amc", treatments = TRUE); print(p_orn_amc) # uh-huh
# Calibration curve data look to be in 2015_09_17_sc_cc.csv (not sure which )

ging_live <- read_csv("data/2015_09_24_sc_phe_arg_amc_live.csv", skip = 3) %>% mutate(treatment = "live")
ging_killed <- read_csv("data/2015_09_24_sc_phe_arg_amc_kill.csv", skip = 3) %>% mutate(treatment = "killed")
ging <- ging_live %>% rbind(ging_killed)
p_ging <- raw_plot(ging, substrate = "phe-arg-AMC", treatments = TRUE); print(p_ging) # I guess?

clos_live <- read_csv("data/2015_09_24_sc_phe_val_arg_amc_live.csv", skip = 3) %>% mutate(treatment = "live")
clos_killed <- read_csv("data/2015_09_24_sc_phe_val_arg_amc_kill.csv", skip = 3) %>% mutate(treatment = "killed")
clos <- clos_live %>% rbind(clos_killed)
p_clos <- raw_plot(clos, substrate = "phe-val-arg-AMC", treatments = TRUE); print(p_clos) # looks like substrate inhibition maybe?

# Put it all in a nested data frame for km calculations
sat_curve_list <- list(mub_cell = mub_cell, 
                       mub_po4 = mub_po4, 
                       mub_xylo = mub_xylo, 
                       mub_a_glu = mub_a_glu, 
                       mub_b_glu = mub_b_glu, 
                       mub_b_nag = mub_b_nag, 
                       arg_amc = arg_amc, 
                       leu_amc = leu_amc, 
                       pro_amc = pro_amc, 
                       orn_amc = orn_amc, 
                       ging = ging, 
                       clos = clos)
sat_data <- map_df(sat_curve_list, identity, .id = "substrate") %>%
  group_by(substrate, conc.uM, treatment) %>%
  mutate(Rtime = hms(time),
         elapsed = as.numeric(time - min(time, na.rm = TRUE))/3600,
         treatment = case_when(is.na(treatment) ~ "live",
                               !is.na(treatment) ~ treatment)) %>%
  nest()

get_slope <- function(x) {
  coef(x)[2]
}

#sat_data <- sat_data %>%
#  mutate(lm_col = map(data, lm_stats, "elapsed", "FSU"),
#         uncal.slope = map_dbl(lm_col, get_slope))
# Calc uncalibrated rates (there is calib data, but I don't know whether to believe it; or what fluorophore it uses)
sat_data <- sat_data %>% #
  #group_by(substrate, conc.uM) %>%
  mutate(lms = map(data, ~lm(.x$FSU ~ .x$elapsed)),
         slope = map_dbl(lms, ~coef(.x)[2]),
         slope.se = map_dbl(lms, ~summary(.x)$coefficients[2, 2]))
  #mutate(lms = map(data, lm_stats(.x, "elapsed", "FSU") )) # doesn't work, WTF

# Define the michaelis menten formula
mm_form <- formula(slope ~ (Vmax*conc.uM)/(Km+conc.uM))

safe_nls2 <- function(data, formula, start = list(Km = 100, Vmax = 500)) {
  tryCatch(nls(formula, data, start=start),
           error = function(e) NA)
}
safe_Km <- function(x) {
  tryCatch(summary(x)$coefficient["Km", "Estimate"],
           error = function(x) NA)
}
safe_Km_err <- function(x) {
  tryCatch(summary(x)$coefficient["Km", "Std. Error"],
           error = function(x) NA)
}

safe_predict <- function(mod) {
  grid <- 0:1600
  preds <- tryCatch(predict(mod, newdata = data.frame(conc.uM = grid)),
                    error = function(e) rep(NA, length(grid)))
  data.frame(conc.uM = grid, slope = preds)
}

# Fit saturation curves
v0_data <- sat_data %>%
  group_by(substrate, treatment) %>%
  nest() %>%
  mutate(nls_fit = map(data, safe_nls2, mm_form),
         Km = map_dbl(nls_fit, safe_Km),
         Km.err = map_dbl(nls_fit, safe_Km_err),
         
         preds = map(nls_fit, safe_predict))

# Maybe try again with clostripain and gingipain?
clos_live_bespoke <- nls(mm_form, 
                         sat_data %>% filter(substrate == "clos" & treatment == "live" & conc.uM <= 500),
                         start = list(Km=250, Vmax = 3000),
                         lower = list(Km=0.1, Vmax = 0.1),
                         algorithm = "port")

ging_live_bespoke <- nls(mm_form,
                         sat_data %>% filter(substrate == "ging" & treatment == "live" & conc.uM <= 500),
                         start = list(Km=250, Vmax = 3000),
                         lower = list(Km=0.1, Vmax = 0.1),
                         algorithm = "port")

# Note: I messed around with a bunch of values, but couldn't get the nls to predict Km > Km's lower bound.
# There's strong evidence of substrate inhibition, so I'll cut off everything greatert than 500 uM
v0_data[v0_data$substrate == "clos" & v0_data$treatment == "live", "Km"] <- get_km(clos_live_bespoke)
v0_data[v0_data$substrate == "clos" & v0_data$treatment == "live", "Km.err"] <- get_km_se(clos_live_bespoke)
v0_data[v0_data$substrate == "clos" & v0_data$treatment == "live", "nls_fit"][1] <- clos_live_bespoke


v0_data[v0_data$substrate == "ging" & v0_data$treatment == "live", "Km"] <- get_km(ging_live_bespoke)
v0_data[v0_data$substrate == "ging" & v0_data$treatment == "live", "Km.err"] <- get_km_se(ging_live_bespoke)
v0_data[v0_data$substrate == "ging" & v0_data$treatment == "live", "nls_fit"] <- ging_live_bespoke

# Predict on linear models
preds <- v0_data %>%
  select(substrate, treatment, preds) %>%
  unnest(cols = preds)


ggplot() + 
  geom_pointrange(data = sat_data, aes(x=conc.uM, y=slope, ymin = slope-slope.se, ymax = slope+slope.se, colour = treatment)) + 
  geom_line(data = preds, aes(x=conc.uM, y=slope, colour = treatment)) +
  facet_wrap(~substrate, scales = "free_y")

# Plot of wo  
