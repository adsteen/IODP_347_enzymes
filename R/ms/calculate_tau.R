R<-1.7e-2 #organic carbon oxidation rate, nmol C g sed-1 hr-1
GE<-0.1 #growth efficiency, unitless
A<-40 #observed clostripain activity, nmol bonds hr-1 g sed-1
SA.vec<-c(1.1e2,2.9e1,2.3e2,9.4e1,3.1e1,2.0e-3) #specific activities for different proteases in literature, nmol bond hr-1 nmol C in enzyme-1

#Calculate enzyme concentration based on observed activity and specific activity (supplemental equation 4)
E.calc<-function(A,SA){A/SA}

#Calculate enzyme lifetime based on enzyme concentration, growth efficiency and oxidation rate (supplemental equation 5)
tau.calc<-function(E,R,GE){E/(R*GE)}


#Enzyme concentrations as a function of specific activity, nmol C in enzyme g sed-1
E<-E.calc(A,SA.vec)

#Enyzme lifetimes as a function of enzyme concnetrations, hours
tau<-tau.calc(E,R,GE)


print(tau)