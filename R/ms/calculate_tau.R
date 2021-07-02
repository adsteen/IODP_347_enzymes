########
# Calculate tau values for various specific activities found in the literature
########

# Load necessary packages
library(tidyverse)
library(stargazer) # for LaTeX table

# Define constants for model
resp<-1.7e-3 #organic carbon oxidation rate, nmol C g sed-1 hr-1
GE<-0.1 #growth efficiency, unitless
Vmax<-40 #observed clostripain activity, nmol bonds hr-1 g sed-1

# Convert specific activity values in the literature to units of nmol bond hydrolyzed per hr per nmol C in enzyme
SA.wang2005<-8258*60*1000*2*1000*12/181/1e9 #covert 5530 @ 20C U mg-1 into nmol bond hr-1 nmol C in enzyme-1
SA.wang2008<-2351*60*1000*2*1000*12/181/1e9 #covert 2351 @ 20C U mg-1 into nmol bond hr-1 nmol C in enzyme-1; U has units ug Y min-1
SA.szwject1992<-160*60*1000*2*1000*12/1e9 #convert 160 U @ 25C mg-1 into nmol bond hr-1 nmol C in enzyme-1; U has units  mol bond min-1
SA.chevalliers1992<-65*60*1000*2*1000*12/1e9 #convert 65 U @ 37C mg-1 into nmol bond hr-1 nmol C in enzyme-1; U has units  mol bond min-1
SA.ullmann1994<-0.36*3600*1000*2*1000*12/1e9 #convert 0.36 U @ 25C mg-1 into nmol bond hr-1 nmol C in enzyme-1; U has units  umol s-1 mg-1 
SA.mcluskey2016<-1587*60*1e-3*2*1000*12/1e9 #convert 1587 U @ 37C mg-1 into nmol bond hr-1 nmol C in enzyme-1; U has units pmol bond min-1

SA.refs <- c("Wang et al, 2005", "Wang et al, 2008", "Szwajcer-Dey et al, 1992", "Chevalliers et al 1992", "Ullman et al 1994", "McLuskey et al 2016")
SA.vec<-c(SA.wang2005,SA.wang2008,SA.szwject1992,SA.chevalliers1992,SA.ullmann1994,SA.mcluskey2016) #specific activities for different proteases in literature, nmol bond hr-1 nmol C in enzyme-1



# Calculate tau
d <- data.frame(SA = SA.vec,
                refs = SA.refs) %>%
  mutate(tau.hr = Vmax / (SA * GE * resp),
         tau.yr = tau.hr/(24*365)) %>%
  select(SA, tau.hr, tau.yr, refs)

# Print a LaTeX table. I'm not using knitr because I'm a chud.
d %>%
  mutate(SA =signif(SA, digits = 2),
         tau.hr = signif(tau.hr, digits = 2),
         tau.yr = signif(tau.yr, digits = 2)) %>%
  stargazer(summary=FALSE)







#Calculate enzyme concentration based on observed activity and specific activity (supplemental equation 4)
E.calc<-function(A,SA){A/SA}

#Calculate enzyme lifetime based on enzyme concentration, growth efficiency and oxidation rate (supplemental equation 5)
tau.calc<-function(E,R,GE){E/(R*GE)}

#Enzyme concentrations as a function of specific activity, nmol C in enzyme g sed-1
E<-E.calc(A,SA.vec)

#Enyzme lifetimes as a function of enzyme concnetrations, hours
tau<-tau.calc(E,R,GE)

d <- data.frame(SA = SA.vec,
                tau = tau,
                refs = SA.refs)
print(d)

library(stargazer)

