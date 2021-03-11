# This is the master script for Jenna's baltic enzymes paper

# Step 1: calculate enzyme activities, including calibrations

# Set teh graphical theme
theme_set(theme_bw() + 
            theme(text = element_text(size=9)))

source("R/Jennas_R/2016_04_20_LB_Vmax_analysis.R")
detach("package:reshape2", unload=TRUE)
detach("package:lubridate", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:plyr", unload=TRUE)



source("R/Jennas_R/2016_04_20_cell_sp_v0.R")

# Generate correlation
library(tidyverse)

# Create v0 column with minimum values
min.v0s <- samp_slopes %>%
  group_by(treatment) %>%
  summarise(min.v0s = min(v0[v0 - v0.se > 0], na.rm=TRUE))
min.v0.killed <- min.v0s$min.v0s[min.v0s$treatment == "killed"]
min.v0.live <- min.v0s$min.v0s[min.v0s$treatment == "live"]


# samp_slopes_2 <- samp_slopes %>%
#   group_by(treatment) %>%
#   mutate(v0.adj = case_when(v0 >= min.v0 ~ v0,
#                             v0 < min(v0, na.rm=TRUE) ~ min(v0, na.rm = TRUE)))


# min live v0


live_killed <- samp_slopes %>%
  select(depth.mbsf, substrate, treatment, v0, v0.se) %>%
  pivot_wider(id_cols = c(depth.mbsf, substrate), # note id_cols defaults to everything not in names_from and values_from
              names_from = treatment,
              values_from = c(v0, v0.se)) %>%
  mutate(v0.adj_killed = case_when(v0_killed < min.v0.killed ~ min.v0.killed,
                                   v0_killed >= min.v0.killed ~ v0_killed),
         v0.adj_live = case_when(v0_live < min.v0.live ~ min.v0.live,
                                 v0_live >= min.v0.live ~ v0_live)) %>%
  mutate(is.neg.x = v0.adj_live - v0.se_live <= 0)

point.color <- "black"
ggplot(live_killed, aes(x=v0.adj_live, y=v0.adj_killed)) + 
  geom_smooth(aes(group=1), method = "lm", color=point.color) +
  geom_crossbar(aes(xmin=v0.adj_live-v0.se_live, xmax=v0.adj_live+v0.se_live), colour=point.color) +
  geom_crossbar(aes(ymin=v0.adj_killed-v0.se_killed, ymax=v0.adj_killed+v0.se_killed), colour=point.color) +
  geom_point(colour="gray10", shape=19) +
  scale_x_log10(name = expression(paste("live ", v[0], ", ", mu, "mol ", g^{-1}, " sed ", hr^{-1}))) + 
  scale_y_log10(name = expression(paste("autoclaved ", v[0], ", ", mu, "mol ", g^{-1}, " sed ", hr^{-1}))) + 
  theme_bw() + 
  theme(text = element_text(size=9))
ggsave("plots/live_vs_killed.png", height=2.5, width=3, units = "in", dpi=300)
