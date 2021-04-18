# Clean v0 calib units

# add in calibration slope for depth 4.5, AMC substrates (mub calib depth changed to 61.55 on 4.5, 
#11.1, 17.6, 24.3, and 37.5 just to see 12-18-15)
samp_slopes[samp_slopes$depth==4.5 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==4.5, "slope"]
# add in calibration slope for depth 4.5, MUB substrates
samp_slopes[samp_slopes$depth==4.5 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==4.5, "slope"]

samp_slopes[samp_slopes$depth==11.1 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==11.1, "slope"]
# add in calibration slope for depth 11.1, MUB substrates
samp_slopes[samp_slopes$depth==11.1 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==4.5, "slope"]

# add in calibration slope for depth 17.6, AMC substrates
samp_slopes[samp_slopes$depth==17.6 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==17.6, "slope"]
# add in calibration slope for depth 17.6, MUB substrates
samp_slopes[samp_slopes$depth==17.6 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==4.5, "slope"]

samp_slopes[samp_slopes$depth==24.3 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==24.3, "slope"]
# add in calibration slope for depth 4.5, MUB substrates
samp_slopes[samp_slopes$depth==24.3 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==30.9, "slope"]

samp_slopes[samp_slopes$depth==30.9 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==30.9, "slope"]
# add in calibration slope for depth 4.5, MUB substrates
samp_slopes[samp_slopes$depth==30.9 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==30.9, "slope"]

# add in calibration slope for depth 37.5, AMC substrates
samp_slopes[samp_slopes$depth==37.5 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==37.5, "slope"]
# add in calibration slope for depth 37.5, MUB substrates
samp_slopes[samp_slopes$depth==37.5 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==30.9, "slope"]

samp_slopes[samp_slopes$depth==43.15 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==43.15, "slope"]
# add in calibration slope for depth 43.15, MUB substrates
samp_slopes[samp_slopes$depth==43.15 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==43.15, "slope"]

samp_slopes[samp_slopes$depth==48.22 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==48.22, "slope"]
# add in calibration slope for depth 48.22, MUB substrates
samp_slopes[samp_slopes$depth==48.22 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==48.22, "slope"]

# add in calibration slope for depth 54.95, AMC substrates
samp_slopes[samp_slopes$depth==54.95 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==54.95, "slope"]
# add in calibration slope for depth 54.95, MUB substrates
samp_slopes[samp_slopes$depth==54.95 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==54.95, "slope"]

# add in calibration slope for depth 61.55, AMC substrates
samp_slopes[samp_slopes$depth==61.55 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==61.55, "slope"]
# add in calibration slope for depth 61.55, MUB substrates
samp_slopes[samp_slopes$depth==61.55 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==61.55, "slope"]

# add in calibration slope for depth 68.12, AMC substrates
samp_slopes[samp_slopes$depth==68.12 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==68.12, "slope"]
# add in calibration slope for depth 68.12, MUB substrates
samp_slopes[samp_slopes$depth==68.12 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==68.12, "slope"]

samp_slopes[samp_slopes$depth==77.92 & samp_slopes$substrate %in% AMC.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="amc-std" & calib_slopes$depth.mbsf==77.92, "slope"]
# add in calibration slope for depth 77.92, MUB substrates
samp_slopes[samp_slopes$depth==77.92 & samp_slopes$substrate %in% MUB.subs, "calib.slope"] <- 
  calib_slopes[calib_slopes$substrate=="mub-std" & calib_slopes$depth.mbsf==77.92, "slope"]