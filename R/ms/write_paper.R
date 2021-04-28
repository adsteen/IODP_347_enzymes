# This is the master script for Jenna's baltic enzymes paper

# Step 1: calculate enzyme activities, including calibrations
library(tidyverse)
# Set the graphical theme
theme_set(theme_bw() + 
            theme(text = element_text(size=9)))

print.plots <- FALSE
print.extra.plots <- FALSE
save.plots <- FALSE

# Most of this code was written in pre-tidyverse syntax
# lm_stats performs linear regressions and returns slopes and std errors, a bit like a combination of map, lm and some broom functions
source("R/ms/lm_stats.R")

# Calculate v0 data for figure 1
source("R/ms/v0.R")
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

# Calculate cell-specific v0 data
source("R/ms/v0_cell_sp.R")


# Look at activities of live vs killed enyzmes
source("R/ms/live_vs_killed_analysis.R")
