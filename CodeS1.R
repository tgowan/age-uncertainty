###########################################################################
### Simulation and model code for Gowan et al. age uncertainty model####
############### Modified December 14, 2020 ####################################


# Age model where state (age) is uncertain
# states = live age-0, live age-1, ..., live age-max, or dead
# No dead recoveries, robust design, or time effects

# Simulate data
# Option to export data for E-Surge
# Then fit NIMBLE model with custom likelihood function


####################
# Set simulation Parameters
####################

n.primary <- 8  #number of primary sampling occassions
marked <- rep(100, n.primary) #total number marked each occasion 

age <- 0:4
n.states <- length(age)+1 #number of states (age0, age1, ..., adult, dead)
n.obs <- 4 #number of observation possibilities (seen calf, seen subadult, seen adult, not seen)

s.0 <- rep(0.7)    #survival probability for age-0 (mean)
bet <- 0.5             #slope for how survival changes with age (on logit scale)
s.a <- round(1 / (1 + exp(-1*(qlogis(s.0) + bet*age))), 2) #survival for each age
s <- matrix(rep(s.a, (n.primary-1)), ncol=(n.primary-1), nrow=length(age)) #survival probabilities for each age (row) and each year (column)

pC <- rep(0.5,n.primary)      #recapture probability for age0, age1 for each occasion
pS <- rep(0.6,n.primary)      #recapture probability for age2, age3 for each occasion
pA <- rep(0.7,n.primary)      #recapture probability for age4+ adults for each occasion

pi0 <- rep(0.1,n.primary) #proportion of individuals age0 at first capture
pi1 <- rep(0.2,n.primary) #proportion of individuals age1 at first capture
pi2 <- rep(0.2,n.primary) #proportion of individuals age2 at first capture
pi3 <- rep(0.2,n.primary) #proportion of individuals age3 at first capture
#proportion of individuals age4+ at first capture is calculated by subtraction

d0c <- rep(0.8,n.primary)  #prob recorded as calf, given age0 and detected
d0s <- rep(0.2,n.primary)  #prob recorded as subadult, given age0 and detected
d1c <- rep(0.6,n.primary)  #prob recorded as calf, given age1 and detected
d1s <- rep(0.4,n.primary)  #prob recorded as subadult, given age1 and detected
d2c <- rep(0.2,n.primary)  #prob recorded as calf, given age2 and detected
d2s <- rep(0.5,n.primary)  #prob recorded as subadult, given age2 and detected
d3c <- rep(0.1,n.primary)  #prob recorded as calf, given age3 and detected
d3s <- rep(0.6,n.primary)  #prob recorded as subadult, given age3 and detected
d4c <- rep(0,n.primary)  #prob recorded as calf, given age4+ and detected
d4s <- rep(0.2,n.primary)  #prob recorded as subadult, given age4+ and detected
#prob recorded as adult is calculated by subtraction for all ages

####################


####################
# Set Matrices for simulation
# Adapted from Kery & Schaub 2012 BPA
####################

# a.State at first capture
PSI.INIT <- array(NA, dim=c(1, n.states-1, n.primary))
for (t in 1:n.primary){ #for each occasion...
  PSI.INIT[,,t] <- matrix(c(
    pi0[t], pi1[t], pi2[t], pi3[t], 1-(pi0[t]+pi1[t]+pi2[t]+pi3[t])) #as age0, 1, 2, 3, or 4+
    , nrow=1, byrow=TRUE)
} #t
PSI.INIT[,,t] #view matrix for last occasion

# b.Observation at first capture
totrel <- sum(marked)
PSI.init.obs <- array(NA, dim=c(n.states, n.obs, totrel, n.primary))
for (i in 1:totrel){ #for each individual...
  for (t in 1:n.primary){ #for each occasion...
    PSI.init.obs[,,i,t] <- matrix(c(
      #recorded as calf, subadult, adult, not seen
      d0c[t], d0s[t], (1-(d0c[t]+d0s[t])), 0, #if age0
      d1c[t], d1s[t], (1-(d1c[t]+d1s[t])), 0, #if age1
      d2c[t], d2s[t], (1-(d2c[t]+d2s[t])), 0, #if age2
      d3c[t], d3s[t], (1-(d3c[t]+d3s[t])), 0, #if age3
      d4c[t], d4s[t], (1-(d4c[t]+d4s[t])), 0, #if age4+
      0,      0,      0,                   1) #if dead
      , nrow=n.states, byrow=TRUE)
  } #t
} #i
PSI.init.obs[,,i,t] #view observation matrix for last individual at last occasion

# 1. Survival and State transition matrix
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.primary-1))
for (i in 1:totrel){ #for each individual...
  for (t in 1:(n.primary-1)){ #for each occasion interval...
    PSI.STATE[,,i,t] <- matrix(c(
      0, s[1,t], 0,      0,      0,      1-s[1,t], #if age0, survives to age1 or dies
      0, 0,      s[2,t], 0,      0,      1-s[2,t], #if age1
      0, 0,      0,      s[3,t], 0,      1-s[3,t], #if age2
      0, 0,      0,      0,      s[4,t], 1-s[4,t], #if age3
      0, 0,      0,      0,      s[5,t], 1-s[5,t], #if age4+
      0, 0,      0,      0,      0,      1)        #if dead: must stay dead
      , nrow=n.states, byrow=TRUE)
  } #t
} #i
PSI.STATE[,,i,t] #view transition matrix for last individual at last time-step

# 2.Observation process at recapture matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.primary))
for (i in 1:totrel){ #for each individual...
  for (t in 1:n.primary){ #for each occasion...
    PSI.OBS[,,i,t] <- matrix(c(
      #recorded as calf, subadult, adult, not seen
      pC[t]*d0c[t], pC[t]*d0s[t], pC[t]*(1-(d0c[t]+d0s[t])), 1-pC[t], #if age0
      pC[t]*d1c[t], pC[t]*d1s[t], pC[t]*(1-(d1c[t]+d1s[t])), 1-pC[t], #if age1
      pS[t]*d2c[t], pS[t]*d2s[t], pS[t]*(1-(d2c[t]+d2s[t])), 1-pS[t], #if age2
      pS[t]*d3c[t], pS[t]*d3s[t], pS[t]*(1-(d3c[t]+d3s[t])), 1-pS[t], #if age3
      pA[t]*d4c[t], pA[t]*d4s[t], pA[t]*(1-(d4c[t]+d4s[t])), 1-pA[t], #if age4+
      0,            0,            0,                         1)       #if dead
      , nrow=n.states, byrow=TRUE)
  } #t
} #i
PSI.OBS[,,i,t] #view observation matrix for last individual at last occasion

####################


####################
# Execute simulation
####################

#Empty matrices to store state and captures
CH <- CH.sur <- matrix(0, ncol=n.primary, nrow=sum(marked))
  
# Define a vector with the occasion of first marking
ent.occ <- numeric()
for (t in 1:n.primary){
  ent.occ <- c(ent.occ, rep(t, marked[t]))
}
  
# Simulating survival and state transitions
#states: 1=age0, 2=age1, 3=age2, 4=age3, 5=age4+, 6=dead, 0=not yet marked
for (i in 1:sum(marked)){
  f.state <- which(rmultinom(1, 1, PSI.INIT[,,ent.occ[i]])==1) #simulate state at first capture
  CH.sur[i, ent.occ[i]] <- f.state
  if (ent.occ[i] == n.primary) next # skip if 1st marked on last occasion
  for (t in (ent.occ[i]+1):n.primary){
    state <- which(rmultinom(1, 1, PSI.STATE[CH.sur[i,t-1],,i,t-1])==1)
    CH.sur[i,t] <- state #CH.sur holds true states of individuals at each occasion
  } #t
} #i
  
# Simulating observation at first capture
#observations: 1=calf, 2=subadult, 3=adult, 4=not seen
for (i in 1:sum(marked)){
  t <- ent.occ[i] #occasion of first capture
  event <- which(rmultinom(1, 1, PSI.init.obs[CH.sur[i,t],,i,t])==1)
  CH[i,t] <- event #CH holds observations
} #i
  
# Simulating recaptures
#observations: 1=calf, 2=subadult, 3=adult, 4=not seen
for (i in 1:sum(marked)){
  for (t in 1:n.primary){
    if (CH.sur[i,t] == 0) next # skip (leave as 0) if not marked yet
    if (t == ent.occ[i]) next  # skip (leave as is) if first marked that occasion
    event <- which(rmultinom(1, 1, PSI.OBS[CH.sur[i,t],,i,t])==1)
    CH[i,t] <- event #CH holds observations
  } #t
} #i
  
# Recode CH matrix
# 1=recorded as calf, 2=recorded as subadult, 3=recorded as adult, 4=not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 4

####################


####################
# Option to Export data for E-Surge in MARK .inp format
####################

#eCH <- rCH
#eCH[eCH==4] <- 0 # Replace undetected events with 0
#write.table(eCH, file='age_uncert.inp', sep="", row.names=FALSE, col.names=FALSE, eol=" 1;\n")
####################




############################################################

####################
# Fit model in NIMBLE
# note, code will need to be modified if number of age classes not equal to 5 or number of events not equal to 4
####################

library(nimble)
library(asbio)

# Parameters for MCMC
n.iter = 4000
n.chains = 3
burnin = 2000

####################
# Custom likelihood function in Nimble for analysis
####################

dmsDHMM <- nimbleFunction(
  run = function(
    # argument types
    x = double(1), length = double(),
    pi.b0 = double(), pi.b1 = double(), pi.b2 = double(), pi.b3 = double(), #initial state (mlogit; modify if # age classes <> 5)
    s.0 = double(), s.1 = double(), s.2 = double(), s.3 = double(), s.4 = double(), #survival (modify if # age classes <> 5)
    p.0 = double(), p.1 = double(), p.2 = double(), p.3 = double(), p.4 = double(), #recapture (modify if # age classes <> 5)
    a0c = double(), a1c = double(), a2c = double(), a3c = double(), a4c = double(), #size class assignment (mlogit)
    a0s = double(), a1s = double(), a2s = double(), a3s = double(), a4s = double(), #   (add or remove parameters if # age classes <> 5 or # events <> 4)
    log = double()) { 
    
    # size class assignment; multinomial logit (add or remove parameters if # age classes <> 5 or # events <> 4)
    d0c <- exp(a0c)/(1 + exp(a0c) + exp(a0s) ) #age0 recorded as calf
    d0s <- exp(a0s)/(1 + exp(a0c) + exp(a0s) ) #age0 recorded as subadult
    d0a <- 1 - d0c - d0s                 #age0 recorded as adult
    d1c <- exp(a1c)/(1 + exp(a1c) + exp(a1s) ) #age1
    d1s <- exp(a1s)/(1 + exp(a1c) + exp(a1s) )
    d1a <- 1 - d1c - d1s
    d2c <- exp(a2c)/(1 + exp(a2c) + exp(a2s) ) #age2
    d2s <- exp(a2s)/(1 + exp(a2c) + exp(a2s) )
    d2a <- 1 - d2c - d2s
    d3c <- exp(a3c)/(1 + exp(a3c) + exp(a3s) ) #age3
    d3s <- exp(a3s)/(1 + exp(a3c) + exp(a3s) )
    d3a <- 1 - d3c - d3s
    d4c <- exp(a4c)/(1 + exp(a4c) + exp(a4s) ) #age4+
    d4s <- exp(a4s)/(1 + exp(a4c) + exp(a4s) )
    d4a <- 1 - d4c - d4s
    
    logL <- 0  # initialize log-likelihood
    
    # vector of initial states
    pi <- numeric(6, init = FALSE) # (modify if # age classes <> 5)
    # initial state: age0 through age3 (add or remove parameters if # age classes <> 5)
    pi[1] <- exp(pi.b0) / (1 + exp(pi.b0) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3))
    pi[2] <- exp(pi.b1) / (1 + exp(pi.b1) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3))
    pi[3] <- exp(pi.b2) / (1 + exp(pi.b2) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3))
    pi[4] <- exp(pi.b3) / (1 + exp(pi.b3) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3))
    # initial state: age4+
    pi[5] <- 1 - pi[1] - pi[2] - pi[3] - pi[4]
    # initial state: dead
    pi[6] <- 0
    
    Zpi <- pi
    if(x[1] == 1) {	# if recorded as calf at 1st capture
      #  add or remove parameters if # age classes <> 5 or # events <> 4
      Zpi[1] <- pi[1] * d0c #age 0
      Zpi[2] <- pi[2] * d1c #age 1
      Zpi[3] <- pi[3] * d2c
      Zpi[4] <- pi[4] * d3c		
      Zpi[5] <- pi[5] * d4c 
      Zpi[6] <- 0				 #dead
    } else if(x[1] == 2) {	# if recorded as subadult
      Zpi[1] <- pi[1] * d0s #age 0
      Zpi[2] <- pi[2] * d1s #age 1
      Zpi[3] <- pi[3] * d2s
      Zpi[4] <- pi[4] * d3s		
      Zpi[5] <- pi[5] * d4s
      Zpi[6] <- 0				 #dead
    } else if(x[1] == 3) {	# if recorded as adult
      Zpi[1] <- pi[1] * d0a #age 0
      Zpi[2] <- pi[2] * d1a #age 1
      Zpi[3] <- pi[3] * d2a
      Zpi[4] <- pi[4] * d3a		
      Zpi[5] <- pi[5] * d4a
      Zpi[6] <- 0				 #dead
    } 
    # x[1] cannot be 4 (not seen) at 1st capture
    
    sumZpi <- sum(Zpi)
    logL <- logL + log(sumZpi)
    
    #survival (modify if # age classes <> 5)
    pi[1] <- 0 
    pi[2] <- Zpi[1] * s.0		
    pi[3] <- Zpi[2] * s.1		
    pi[4] <- Zpi[3] * s.2			
    pi[5] <- Zpi[4] * s.3 + Zpi[5] * s.4
    pi[6] <- Zpi[1] * (1 - s.0) + Zpi[2] * (1 - s.1) + Zpi[3] * (1 - s.2) + Zpi[4] * (1 - s.3) + Zpi[5] * (1 - s.4) + Zpi[6]	
    
    pi <- pi/sumZpi
    
    for(t in 2:length) { #for remaining occasions
      Zpi <- pi
      if(x[t] == 1) {	# if recorded as calf
        #  add or remove parameters if # age classes <> 5 or # events <> 4
        Zpi[1] <- pi[1] * p.0 * d0c #age 0
        Zpi[2] <- pi[2] * p.1 * d1c 
        Zpi[3] <- pi[3] * p.2 * d2c
        Zpi[4] <- pi[4] * p.3 * d3c		
        Zpi[5] <- pi[5] * p.4 * d4c
        Zpi[6] <- 0
      } else if(x[t] == 2) {	# if recorded as subadult
        Zpi[1] <- pi[1] * p.0 * d0s #age 0
        Zpi[2] <- pi[2] * p.1 * d1s
        Zpi[3] <- pi[3] * p.2 * d2s
        Zpi[4] <- pi[4] * p.3 * d3s
        Zpi[5] <- pi[5] * p.4 * d4s
        Zpi[6] <- 0
      } else if(x[t] == 3) {	# if recorded as adult
        Zpi[1] <- pi[1] * p.0 * d0a #age 0
        Zpi[2] <- pi[2] * p.1 * d1a
        Zpi[3] <- pi[3] * p.2 * d2a
        Zpi[4] <- pi[4] * p.3 * d3a
        Zpi[5] <- pi[5] * p.4 * d4a
        Zpi[6] <- 0
      } else if(x[t] == 4) {	# if not seen
        Zpi[1] <- pi[1] * (1 - p.0)
        Zpi[2] <- pi[2] * (1 - p.1)
        Zpi[3] <- pi[3] * (1 - p.2)
        Zpi[4] <- pi[4] * (1 - p.3)
        Zpi[5] <- pi[5] * (1 - p.4)
        Zpi[6] <- pi[6]				
      }
      
      sumZpi = sum(Zpi)
      logL <- logL + log(sumZpi)
      
      if(t != length) { 
        # modify if # age classes <> 5
        pi[1] <- 0
        pi[2] <- Zpi[1] * s.0
        pi[3] <- Zpi[2] * s.1
        pi[4] <- Zpi[3] * s.2
        pi[5] <- Zpi[4] * s.3 + Zpi[5] * s.4
        pi[6] <- Zpi[1] * (1 - s.0) + Zpi[2] * (1 - s.1) + Zpi[3] * (1 - s.2) + Zpi[4] * (1 - s.3) + Zpi[5] * (1 - s.4) + Zpi[6]
      }			
      
      pi <- pi/sumZpi
      
    }
    returnType(double())
    if(log){
      return(logL)
    } else{
      return(exp(logL))
    }
  }
)

rmsDHMM <- nimbleFunction( #this function not used but required to register distribution in Nimble
  run = function(
    n = integer(), length = double(),
    pi.b0 = double(), pi.b1 = double(), pi.b2 = double(), pi.b3 = double(), #initial state (modify if # age classes <> 5)
    s.0 = double(), s.1 = double(), s.2 = double(), s.3 = double(),	s.4 = double(), #survival (modify if # age classes <> 5)
    p.0 = double(), p.1 = double(), p.2 = double(), p.3 = double(), p.4 = double(), #recapture (modify if # age classes <> 5)
    a0c = double(), a1c = double(), a2c = double(), a3c = double(), a4c = double(), #size class assignment
    a0s = double(), a1s = double(), a2s = double(), a3s = double(), a4s = double() ) {  #(add or remove parameters if # age classes <> 5 or # events <> 4)
    
    x <- rep(1, length)
    returnType(double(1))
    return(x)
  }
)

registerDistributions(list(
  dmsDHMM = list(
    # add or remove parameters if # age classes <> 5 or # events <> 4
    BUGSdist = 'dmsDHMM(length, pi.b0, pi.b1, pi.b2, pi.b3,
                                s.0, s.1, s.2, s.3, s.4, p.0, p.1, p.2, p.3, p.4,                        
		                            a0c, a1c, a2c, a3c, a4c, a0s, a1s, a2s, a3s, a4s)',
    types = c('value = double(1)', 'length = double()',
              'pi.b0 = double()', 'pi.b1 = double()', 'pi.b2 = double()', 'pi.b3 = double()',
              's.0 = double()', 's.1 = double()', 's.2 = double()', 's.3 = double()',	's.4 = double()',
              'p.0 = double()', 'p.1 = double()', 'p.2 = double()', 'p.3 = double()',	'p.4 = double()',               
              'a0c = double()', 'a1c = double()', 'a2c = double()', 'a3c = double()', 'a4c = double()',
              'a0s = double()', 'a1s = double()', 'a2s = double()', 'a3s = double()', 'a4s = double()'),
    discrete = TRUE
  )
))


####################
# BUGS model in Nimble
####################

MS.model <- nimbleCode({
  
  ## Priors
  # survival as logistic function of age
  mu.s ~ dnorm(0, 0.01) #intercept
  beta.s ~ dnorm(0, 0.01) #slope
    s.est[1:5] <- 1 / (1+exp(-1*(mu.s + beta.s*(0:4)))) #modify if age classes are not 0-4
  
  # initial state (modify if # age classes <> 5)
  #  could alternatively use Dirichlet prior (see Code S2)
  pi.b0 ~ dnorm(0, 0.01)
  pi.b1 ~ dnorm(0, 0.01)
  pi.b2 ~ dnorm(0, 0.01)
  pi.b3 ~ dnorm(0, 0.01)
    pi1 <- exp(pi.b0) / (1 + exp(pi.b0) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3)) #age0
    pi2 <- exp(pi.b1) / (1 + exp(pi.b1) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3)) #age1
    pi3 <- exp(pi.b2) / (1 + exp(pi.b2) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3)) #age2
    pi4 <- exp(pi.b3) / (1 + exp(pi.b3) + exp(pi.b1) + exp(pi.b2) + exp(pi.b3)) #age3
    pi5 <- 1 - pi1 - pi2 - pi3 - pi4                                            #age4+
  
  # recapture
  p.c ~ dunif(0, 1)
    p.0 <- p.c
    p.1 <- p.c
  p.s ~ dunif(0, 1)
    p.2 <- p.s
    p.3 <- p.s
  p.a ~ dunif(0, 1)
    p.4 <- p.a #(modify if # age classes <> 5)
  
  # size class assignments (add or remove parameters if # age classes <> 5 or # events <> 4)
  a0c ~ dnorm(0, 0.01)
  a0s ~ dnorm(0, 0.01)
    d0c <- exp(a0c)/(1 + exp(a0c) + exp(a0s) ) #age0 recorded as calf
    d0s <- exp(a0s)/(1 + exp(a0c) + exp(a0s) ) #age0 recorded as subadult
    d0a <- 1 - d0c - d0s                 #age0 recorded as adult
  a1c ~ dnorm(0, 0.01)
  a1s ~ dnorm(0, 0.01)
    d1c <- exp(a1c)/(1 + exp(a1c) + exp(a1s) ) #age1
    d1s <- exp(a1s)/(1 + exp(a1c) + exp(a1s) )
    d1a <- 1 - d1c - d1s
  a2c ~ dnorm(0, 0.01)
  a2s ~ dnorm(0, 0.01)
    d2c <- exp(a2c)/(1 + exp(a2c) + exp(a2s) ) #age2
    d2s <- exp(a2s)/(1 + exp(a2c) + exp(a2s) )
    d2a <- 1 - d2c - d2s
  a3c ~ dnorm(0, 0.01)
  a3s ~ dnorm(0, 0.01)
    d3c <- exp(a3c)/(1 + exp(a3c) + exp(a3s) ) #age3
    d3s <- exp(a3s)/(1 + exp(a3c) + exp(a3s) )
    d3a <- 1 - d3c - d3s
  a4c ~ dnorm(0, 0.01)
  a4s ~ dnorm(0, 0.01)
    d4c <- exp(a4c)/(1 + exp(a4c) + exp(a4s) ) #age4+
    d4s <- exp(a4s)/(1 + exp(a4c) + exp(a4s) )
    d4a <- 1 - d4c - d4s

  ## Likelihood
  for(i in 1:nind){
    y[i, f[i]:n.occasions] ~ dmsDHMM(length = n.occasions - f[i] + 1,
                                     pi.b0 = pi.b0, pi.b1 = pi.b1, pi.b2 = pi.b2, pi.b3 = pi.b3,
                                     s.0 = s.est[1], s.1 = s.est[2], s.2 = s.est[3], s.3 = s.est[4], s.4 = s.est[5],
                                     p.0 = p.0, p.1 = p.1, p.2 = p.2, p.3 = p.3, p.4 = p.4,                       
                                     a0c = a0c, a1c = a1c, a2c = a2c, a3c = a3c, a4c = a4c,
                                     a0s = a0s, a1s = a1s, a2s = a2s, a3s = a3s, a4s = a4s) #add or remove parameters if # age classes <> 5 or # events <> 4
  }
  
})

#Initial values and parameters to monitor
init <- function(){list(mu.s=rnorm(1), beta.s=rnorm(1), #survival, logit
                        pi.b0=rnorm(1,0,1), pi.b1=rnorm(1,0,1), pi.b2=rnorm(1,0,1), pi.b3=rnorm(1,0,1), #initial state, m-logit
                        p.c=runif(1,0,1), p.s=runif(1,0,1), p.a=runif(1,0,1),  #recapture, real
                        a0c=rnorm(1,0,1), a1c=rnorm(1,0,1), a2c=rnorm(1,0,1), a3c=rnorm(1,0,1), a4c=rnorm(1,0,1),  #size class assignment, m-logit
                        a0s=rnorm(1,0,1), a1s=rnorm(1,0,1), a2s=rnorm(1,0,1), a3s=rnorm(1,0,1), a4s=rnorm(1,0,1)  )}
parameters = list("mu.s", "beta.s", "s.est",
                  "pi1", "pi2", "pi3", "pi4", "pi5",
                  "p.c", "p.s", "p.a",
                  "d0c", "d1c", "d2c", "d3c", "d4c",
                  "d0s", "d1s", "d2s", "d3s", "d4s",
                  "d0a", "d1a", "d2a", "d3a", "d4a")

####################




####################
# Fit model
####################

# Function to get occasion of 1st capture
get.first <- function(x) min(which(x!=4))
f <- apply(rCH, 1, get.first)

# Exclude individuals 1st captured on last occasion
ch <- rCH[!(f>=ncol(rCH)),]
fs <- f[!(f>=ncol(rCH))]

# package data
bugs.data <- list(y = ch)
constants <- list(f=fs, nind=dim(ch)[1], n.occasions=dim(ch)[2])

#Compile model in Nimble
MS.test = nimbleModel(code=MS.model, name="MS.test", constants=constants, data=bugs.data, inits=init() )
MS.testConfig = configureMCMC(MS.test, monitors=parameters)
MS.testMCMC = buildMCMC(MS.testConfig)
Cmodel = compileNimble(MS.test)
MS.testCompile = compileNimble(MS.testMCMC, project=MS.test)

#Run model
out = runMCMC(MS.testCompile, niter=n.iter, nburnin=burnin, nchains=n.chains, inits=init)

#summarize posteriors
parnames <- colnames(out[[1]])
samples = array(NA, dim=c(dim(out[[1]]), n.chains))
for(chain in 1:n.chains){ #save results from each chain
  samples[,,chain] = out[[chain]]
}
smmry = matrix(nrow = dim(samples)[2], ncol = 4)
colnames(smmry) = c("Mean", "Prcntl.025", "Prcntl.975", "Rhat")
rownames(smmry) = colnames(out[[1]])
for(j in 1:dim(samples)[2]){ #for each parameter
  smmry[j,1] = mean(samples[,j,])
  smmry[j,2] = quantile(samples[,j,], 0.025)
  smmry[j,3] = quantile(samples[,j,], 0.975)
  smmry[j,4] = R.hat(samples[,j,], 0)
}
smmry
