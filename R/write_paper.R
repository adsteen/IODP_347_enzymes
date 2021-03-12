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

# Need to also set min error ranges that are large enough to not screw up the plot



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
  mutate(v0.min_live = case_when(v0.adj_live - v0.se_live <= 1e-5 ~ 1e-5,
                                 TRUE ~ v0.adj_live - v0.se_live),
         v0.min_killed = case_when(v0.adj_killed - v0.se_killed <= 1e-5 ~ 1e-5,
                                   TRUE ~ v0.adj_killed - v0.se_killed),
         se.is.truncated = case_when((v0.min_live == 1e-5 | v0.min_killed == 1e-5)  ~ TRUE,
                                     TRUE ~ FALSE)
         )
 
line.size = 0.25
p_live_killed <- ggplot(live_killed, aes(x=v0.adj_live, y=v0.adj_killed)) + 
  geom_smooth(aes(group=1), method = "lm", color="black") +
  geom_crossbar(aes(xmin=v0.min_live, xmax=v0.adj_live, colour = se.is.truncated), size = line.size) +
  geom_crossbar(aes(xmin=v0.adj_live, xmax=v0.adj_live+v0.se_live), size = line.size) +
  geom_crossbar(aes(ymin=v0.min_killed, ymax=v0.adj_killed, colour=se.is.truncated), size = line.size) +
  geom_crossbar(aes(ymin=v0.adj_killed, ymax=v0.adj_killed+v0.se_killed), size = line.size) +
  geom_point(size=0.5) +
  scale_x_log10(name = expression(paste("live ", v[0], ", ", mu, "mol ", g^{-1}, " sed ", hr^{-1}))) + 
  scale_y_log10(name = expression(paste("autoclaved ", v[0], ", ", mu, "mol ", g^{-1}, " sed ", hr^{-1}))) + 
  scale_colour_manual(values = c("black", "gray75"), guide=FALSE) +
  theme_bw() + 
  theme(text = element_text(size=9))
ggsave("plots/live_vs_killed.png", p_live_killed, height=2.5, width=3, units = "in", dpi=300)


live_killed_mod <- lm(v0.adj_killed ~ v0.adj_live, data = live_killed)
summary(live_killed_mod)

live_killed_mod_subs <- lm(v0.adj_killed ~ v0.adj_live + substrate, data = live_killed)
summary(live_killed_mod_subs)
