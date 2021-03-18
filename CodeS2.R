###########################################################################
### Simulation and model code for Gowan et al. age uncertainty model####
############### Modified Sept 28, 2020 ####################################


# Age model where state (age) is uncertain
# states = live age-0, live age-1, ..., live age-max, or dead
# Robust design
# No dead recoveries or time effects
# State assignment probabilities modeled as trend across age
# Individual random effects on state assignment probabilities

# Simulate data
# Then fit NIMBLE model with custom likelihood function


####################
# Set simulation Parameters
####################

n.primary <- 8  #number of primary sampling occassions
n.secondary <- 2  #number of secondary sampling occassions
n.occasions <- n.primary*n.secondary
open <- rep(c(0,1), n.primary) #which occasions are open for survival/state transitions
marked <- rep(50, n.occasions) #total number marked each occasion 

age <- 0:4
n.states <- length(age)+1 #number of states (age0, age1, ..., adult, dead)
n.obs <- 4 #number of observation possibilities (seen calf, seen subadult, seen adult, not seen)

s.0 <- rep(0.7)    #survival probability for age-0 (mean)
bet <- 0.5             #slope for how survival changes with age (on logit scale)
s.a <- round(1 / (1 + exp(-1*(qlogis(s.0) + bet*age))), 2) #survival for each age
s <- matrix(rep(s.a, (n.primary-1)), ncol=(n.primary-1), nrow=length(age)) #survival probabilities for each age (row) and each year (column)

pC <- 0.5     #recapture probability for age0, age1 for each occasion
pS <- 0.6      #recapture probability for age2, age3 for each occasion
pA <- 0.7      #recapture probability for age4+ adults for each occasion

pi0 <- 0.1 #proportion of individuals age0 at first capture
pi1 <- 0.2 #proportion of individuals age1 at first capture
pi2 <- 0.2 #proportion of individuals age2 at first capture
pi3 <- 0.2 #proportion of individuals age3 at first capture
#proportion of individuals age4+ at first capture is calculated by subtraction

#trend parameters for how assignment probabilities vary with age
bc.int <- 1.6
bc.slope <- -1 #probability of being recorded as calf decreases with age
ba.int <- -3
ba.slope <- 1.2 #probability of being recorded as adult increases with age

#individual variation in assignment probabilities
bc.var <- 0.5 #recorded as calf, variance
ec <- rnorm(sum(marked), 0, bc.var^0.5)
ba.var <- 0.5 #recorded as adult, variance
ea <- rnorm(sum(marked), 0, ba.var^0.5)

#prob (logit scale) recorded as event, given true age and detected
ac.mat <- matrix(nrow=sum(marked), ncol=length(age)) #store values for calf: each individual at each age
aa.mat <- matrix(nrow=sum(marked), ncol=length(age)) #store values for adult: each individual at each age
for(i in 1:sum(marked)){
  ac.mat[i,] <- bc.int + (bc.slope * age) + ec[i]
  aa.mat[i,] <- ba.int + (ba.slope * age) + ea[i]
}

# derive state assignment probabilities on real scale
dc <- matrix(nrow=sum(marked), ncol=length(age)) #store values for recorded as calf
ds <- matrix(nrow=sum(marked), ncol=length(age)) #store values for recorded as subadult
da <- matrix(nrow=sum(marked), ncol=length(age)) #store values for recorded as adult
for(i in 1:sum(marked)){
  dc[i,] <- exp(ac.mat[i,])/(1 + exp(ac.mat[i,]) + exp(aa.mat[i,]) ) #recorded as calf
  da[i,] <- exp(aa.mat[i,])/(1 + exp(ac.mat[i,]) + exp(aa.mat[i,]) ) #recorded as adult
  ds[i,] <-  1 - dc[i,] - da[i,]                             #recorded as subadult
}

####################


####################
# Set Matrices for simulation
# Adapted from Kery & Schaub 2012 BPA
####################

# a.State at first capture
PSI.INIT <- array(NA, dim=c(1, n.states-1, n.occasions))
for (t in 1:n.occasions){ #for each occasion...
  PSI.INIT[,,t] <- matrix(c(
    pi0, pi1, pi2, pi3, 1-(pi0+pi1+pi2+pi3)) #as age0, 1, 2, 3, or 4+
    , nrow=1, byrow=TRUE)
} #t
PSI.INIT[,,t] #view matrix for last occasion

# b.Observation at first capture
totrel <- sum(marked)
PSI.init.obs <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions))
for (i in 1:totrel){ #for each individual...
  for (t in 1:n.occasions){ #for each occasion...
    PSI.init.obs[,,i,t] <- matrix(c(
      #recorded as calf, subadult, adult, not seen
      dc[i,1], ds[i,1], da[i,1], 0, #if age0
      dc[i,2], ds[i,2], da[i,2], 0, #if age1
      dc[i,3], ds[i,3], da[i,3], 0, #if age2
      dc[i,4], ds[i,4], da[i,4], 0, #if age3
      dc[i,5], ds[i,5], da[i,5], 0, #if age4+
      0,      0,      0,                   1) #if dead
      , nrow=n.states, byrow=TRUE)
  } #t
} #i
round(PSI.init.obs[,,i,t], 2) #view observation matrix for last individual at last occasion

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
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions))
for (i in 1:totrel){ #for each individual...
  for (t in 1:n.occasions){ #for each occasion...
    PSI.OBS[,,i,t] <- matrix(c(
      #recorded as calf, subadult, adult, not seen
      pC*dc[i,1], pC*ds[i,1], pC*da[i,1], 1-pC, #if age0
      pC*dc[i,2], pC*ds[i,2], pC*da[i,2], 1-pC, #if age1
      pS*dc[i,3], pS*ds[i,3], pS*da[i,3], 1-pS, #if age2
      pS*dc[i,4], pS*ds[i,4], pS*da[i,4], 1-pS, #if age3
      pA*dc[i,5], pA*ds[i,5], pA*da[i,5], 1-pA, #if age4+
      0,            0,            0,                         1)       #if dead
      , nrow=n.states, byrow=TRUE)
  } #t
} #i
round(PSI.OBS[,,i,t], 2) #view observation matrix for last individual at last occasion

####################


####################
# Execute simulation
####################

#Empty matrices to store state and captures
CH.sur <- matrix(0, ncol=n.primary, nrow=sum(marked))
CH <- CH.sur.aug <- matrix(0, ncol=n.occasions, nrow=sum(marked))

# Define a vector with the occasion of first marking
ent.occ <- numeric()
for (t in 1:n.occasions){
  ent.occ <- c(ent.occ, rep(t, marked[t]))
}
  
# Simulating survival and state transitions
#states: 1=age0, 2=age1, 3=age2, 4=age3, 5=age4+, 6=dead, 0=not yet marked
for (i in 1:sum(marked)){
  f.state <- which(rmultinom(1, 1, PSI.INIT[,,ent.occ[i]])==1) #simulate state at first capture
  p <- ceiling(ent.occ[i]/n.secondary) #primary period of first capture
  CH.sur[i,p] <- f.state #assuming state is constant within primary period
  if (p == n.primary) next # skip if first marked on last primary period
  for (t in (p+1):n.primary){
    state <- which(rmultinom(1, 1, PSI.STATE[CH.sur[i,t-1],,i,t-1])==1)
    CH.sur[i,t] <- state #CH.sur holds true states of individuals at each occasion
  } #t
} #i
  
#fill in state info for secondary occasions
for(i in 1:sum(marked)) {
  x<-numeric(0)
  for(j in 1:n.primary) {
    x<-append(x,rep(CH.sur[i,j],n.secondary))
  }
  CH.sur.aug[i,]<-x
}

# Simulating observations (events)
#observations: 1=calf, 2=subadult, 3=adult, 4=not seen, 0=not yet marked
for (i in 1:sum(marked)){
  CH[i,ent.occ[i]] <- which(rmultinom(1, 1, PSI.init.obs[CH.sur.aug[i,ent.occ[i]],,i,ent.occ[i]])==1) #at first capture
  if (ent.occ[i] == n.occasions) next # skip if first captured on last occasion
  for (t in (ent.occ[i]+1):n.occasions){
    event <- which(rmultinom(1, 1, PSI.OBS[CH.sur.aug[i,t],,i,t])==1)
    CH[i,t] <- event #CH holds observations
  } #t
} #i
  
# Recode CH matrix
# 1=recorded as calf, 2=recorded as subadult, 3=recorded as adult, 4=not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 4

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
burnin = 1000

####################
# Custom likelihood function in Nimble for analysis, Robust design
####################

dmsrdDHMM <- nimbleFunction(
  run = function(
    # argument types
    x = double(1), length = double(), open = double(1),
    pi.0 = double(), pi.1 = double(), pi.2 = double(), pi.3 = double(), pi.4 = double(), #initial state
    s.0 = double(), s.1 = double(), s.2 = double(), s.3 = double(), s.4 = double(), #survival
    p.0 = double(), p.1 = double(), p.2 = double(), p.3 = double(), p.4 = double(), #recapture
    d0c = double(), d1c = double(), d2c = double(), d3c = double(), d4c = double(), #size class assignment
    d0s = double(), d1s = double(), d2s = double(), d3s = double(), d4s = double(),
    d0a = double(), d1a = double(), d2a = double(), d3a = double(), d4a = double(),
    log = double()) { 
    
    logL <- 0  # initialize log-likelihood
    
    # vector of initial states
    pi <- numeric(6, init = FALSE)
    # initial state: age0 through age4+
    pi[1] <- pi.0
    pi[2] <- pi.1
    pi[3] <- pi.2
    pi[4] <- pi.3
    pi[5] <- pi.4
    # initial state: dead
    pi[6] <- 0
    
    Zpi <- pi
    if(x[1] == 1) {	# if recorded as calf at 1st capture
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
    
    #survival
    if(open[1] == 1) { #don't evaulate survival/transitions within a primary
      pi[1] <- 0
      pi[2] <- Zpi[1] * s.0
      pi[3] <- Zpi[2] * s.1
      pi[4] <- Zpi[3] * s.2
      pi[5] <- Zpi[4] * s.3 + Zpi[5] * s.4
      pi[6] <- Zpi[1] * (1 - s.0) + Zpi[2] * (1 - s.1) + Zpi[3] * (1 - s.2) + Zpi[4] * (1 - s.3) + Zpi[5] * (1 - s.4) + Zpi[6]
      
      pi <- pi/sumZpi
    } else{
      pi <- Zpi/sumZpi
    }
    
    for(t in 2:length) { #for remaining occasions
      Zpi <- pi
      if(x[t] == 1) {	# if recorded as calf
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
      
      if(open[t] == 1 & t != length) { #don't evaulate survival/transitions within a primary or after last occasion
        pi[1] <- 0
        pi[2] <- Zpi[1] * s.0
        pi[3] <- Zpi[2] * s.1
        pi[4] <- Zpi[3] * s.2
        pi[5] <- Zpi[4] * s.3 + Zpi[5] * s.4
        pi[6] <- Zpi[1] * (1 - s.0) + Zpi[2] * (1 - s.1) + Zpi[3] * (1 - s.2) + Zpi[4] * (1 - s.3) + Zpi[5] * (1 - s.4) + Zpi[6]
      
        pi <- pi/sumZpi
      } else{
        pi <- Zpi/sumZpi
      }
      
    }
    returnType(double())
    if(log){
      return(logL)
    } else{
      return(exp(logL))
    }
  }
)

rmsrdDHMM <- nimbleFunction( #this function not used but required to register distribution in Nimble
  run = function(
    n = integer(), length = double(), open = double(1),
    pi.0 = double(), pi.1 = double(), pi.2 = double(), pi.3 = double(), pi.4 = double(), #initial state
    s.0 = double(), s.1 = double(), s.2 = double(), s.3 = double(),	s.4 = double(), #survival
    p.0 = double(), p.1 = double(), p.2 = double(), p.3 = double(), p.4 = double(), #recapture
    d0c = double(), d1c = double(), d2c = double(), d3c = double(), d4c = double(), #size class assignment
    d0s = double(), d1s = double(), d2s = double(), d3s = double(), d4s = double(),
    d0a = double(), d1a = double(), d2a = double(), d3a = double(), d4a = double() ) { 
    
    x <- rep(1, length)
    returnType(double(1))
    return(x)
  }
)

registerDistributions(list(
  dmsrdDHMM = list(
    BUGSdist = 'dmsrdDHMM(length, open, pi.0, pi.1, pi.2, pi.3, pi.4,
                                s.0, s.1, s.2, s.3, s.4, p.0, p.1, p.2, p.3, p.4,                        
		                            d0c, d1c, d2c, d3c, d4c, d0s, d1s, d2s, d3s, d4s, d0a, d1a, d2a, d3a, d4a)',
    types = c('value = double(1)', 'length = double()', 'open = double(1)',
              'pi.0 = double()', 'pi.1 = double()', 'pi.2 = double()', 'pi.3 = double()', 'pi.4 = double()',
              's.0 = double()', 's.1 = double()', 's.2 = double()', 's.3 = double()',	's.4 = double()',
              'p.0 = double()', 'p.1 = double()', 'p.2 = double()', 'p.3 = double()',	'p.4 = double()',               
              'd0c = double()', 'd1c = double()', 'd2c = double()', 'd3c = double()', 'd4c = double()',
              'd0s = double()', 'd1s = double()', 'd2s = double()', 'd3s = double()', 'd4s = double()',
              'd0a = double()', 'd1a = double()', 'd2a = double()', 'd3a = double()', 'd4a = double()'),
    discrete = TRUE
  )
))


####################
# BUGS model in Nimble
####################

rdMS.model <- nimbleCode({
  
  ## Priors
  # survival as logistic function of age
  mu.s ~ dnorm(0, 0.01) #intercept
  beta.s ~ dnorm(0, 0.01) #slope
    s.est[1:5] <- 1 / (1+exp(-1*(mu.s + beta.s*(0:4))))
  
  # initial state, Dirichlet prior (could alternatively use m-logit)
  for(n in 1:5) { #for each age
    a[n] ~ dgamma(1, 1)  
    pi[n] <- a[n]/sum(a[1:5])
  }

  # recapture
  p.c ~ dunif(0, 1)
    p.0 <- p.c
    p.1 <- p.c
  p.s ~ dunif(0, 1)
    p.2 <- p.s
    p.3 <- p.s
  p.a ~ dunif(0, 1)
    p.4 <- p.a
  
  # size class assignments as function of age
  bc.int ~ dnorm(0, 0.01) #recorded as calf, intercept
  bc.slope ~ dnorm(0, 0.01) #recorded as calf, slope
  ba.int ~ dnorm(0, 0.01) #recorded as adult, intercept
  ba.slope ~ dnorm(0, 0.01) #recorded as adult, slope
  
  #individual variation in assignment probabilities
  bc.sigma ~ dunif(0, 5) #calf, std dev
  ba.sigma ~ dunif(0, 5) #adult, std dev
  bc.var <- pow(bc.sigma, 2) #variance
  ba.var <- pow(ba.sigma, 2)
  
  # derive assignment probabilities
  for(i in 1:nind){
    ec[i] ~ dnorm(0, sd = bc.sigma)
    ea[i] ~ dnorm(0, sd = ba.sigma)
    #logit scale
    ac.mat[i,1:5] <- bc.int + (bc.slope*(0:4)) + ec[i] #calf at each age
    aa.mat[i,1:5] <- ba.int + (ba.slope*(0:4)) + ea[i] #adult at each age
    #real scale
    dc[i,1:5] <- exp(ac.mat[i,1:5])/(1 + exp(ac.mat[i,1:5]) + exp(aa.mat[i,1:5]) ) #recorded as calf
    da[i,1:5] <- exp(aa.mat[i,1:5])/(1 + exp(ac.mat[i,1:5]) + exp(aa.mat[i,1:5]) ) #recorded as adult
    ds[i,1:5] <-  1 - dc[i,1:5] - da[i,1:5]                             #recorded as subadult
  }

  ## Likelihood
  for(i in 1:nind){
    y[i, f[i]:n.occasions] ~ dmsrdDHMM(length = n.occasions - f[i] + 1, open = open[f[i]:n.occasions],
                                     pi.0 = pi[1], pi.1 = pi[2], pi.2 = pi[3], pi.3 = pi[4], pi.4 = pi[5],
                                     s.0 = s.est[1], s.1 = s.est[2], s.2 = s.est[3], s.3 = s.est[4], s.4 = s.est[5],
                                     p.0 = p.0, p.1 = p.1, p.2 = p.2, p.3 = p.3, p.4 = p.4,
                                     d0c = dc[i,1], d1c = dc[i,2], d2c = dc[i,3], d3c = dc[i,4], d4c = dc[i,5],
                                     d0s = ds[i,1], d1s = ds[i,2], d2s = ds[i,3], d3s = ds[i,4], d4s = ds[i,5],
                                     d0a = da[i,1], d1a = da[i,2], d2a = da[i,3], d3a = da[i,4], d4a = da[i,5])
  }
  
})

#Initial values and parameters to monitor
init <- function(){list(mu.s=rnorm(1), beta.s=rnorm(1), #survival, logit-trend
                        a=rgamma(5,1,1), #initial state, Dirichlet
                        p.c=runif(1,0,1), p.s=runif(1,0,1), p.a=runif(1,0,1),  #recapture, real
                        bc.int=rnorm(1,0,1), bc.slope=rnorm(1,0,1), ba.int=rnorm(1,0,1), ba.slope=rnorm(1,0,1), #size class assignment, logit-trend
                        bc.sigma=runif(1,0,5), ba.sigma=runif(1,0,5) )} #size class assignment, individual variation
parameters = list("mu.s", "beta.s", "s.est",
                  "pi",
                  "p.c", "p.s", "p.a")

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
bugs.data <- list(y=ch, open=open)
constants <- list(f=fs, nind=dim(ch)[1], n.occasions=dim(ch)[2])

#Compile model in Nimble
Rmodel = nimbleModel(code=rdMS.model, constants=constants, data=bugs.data, inits=init() )
conf = configureMCMC(Rmodel, monitors=parameters)
Rmcmc = buildMCMC(conf)
Cmodel = compileNimble(Rmodel)
Cmcmc = compileNimble(Rmcmc, project=Rmodel)

# Run model and save MCMC samples
out = runMCMC(Cmcmc, niter=n.iter, nburnin=burnin, nchains=n.chains, inits=init)

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
