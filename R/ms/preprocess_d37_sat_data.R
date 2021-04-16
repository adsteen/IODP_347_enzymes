######
# Dealing with saturation curves for the 37 mbsf sample, which apparently took place over 48 hours
#######


d37 <- read.csv("data/2016_02_27_km_37.5mbsf.csv") %>%
  mutate(Rtime = NA,
         elapsed = case_when(time.hr == "0:00" ~ 0,
                             time.hr == "48:00:00" ~ 48)) %>%
  rename(c("time" = "time.hr"))



