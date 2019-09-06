# Sequential Bayes (indirect effect version) for the simulation study 
# Camiel van Zundert


## Function  ####

Seqbay_indirect <- function(gen_data, eta = 1.5, lambda = 0.08,
                            var = 5,  iterations = 5000, 
                            chains = 3, burn = 1000, thin = 1){
  
  require(rjags)
  
  # Input conversion
  y.i        <- sapply(gen_data, function(X) X[[1]][4,2])
  sigma.sq.i <- sapply(gen_data, function(X) X[[1]][4,3]**2) # Since inverse of variance, not SE
  
  # Model
  indirect_effect_model <- "
  model { 
  for (i in 1:n) { 
  y.i[i]     ~ dnorm(theta.i[i], 1 / sigma.sq.i[i])
  theta.i[i] ~ dnorm(theta, inv.tsq)} # JAGS uses precision not variance 
  
  theta      ~ dnorm(mu, prec) 
  inv.tsq    ~ dgamma(eta, lambda)      # Appropriate prior for precision 
  tsq        <- 1 / inv.tsq 
  }
  "
  
  # Initial values for each chain, random seed selected
  
  inits1 <- list(theta   = -0.92 ,         # Random number
                 inv.tsq =  6.04 ,         # Random number
                 .RNG.name="base::Super-Duper", .RNG.seed=     4163370) # seed
  inits2 <- list(theta   = 1.12 , 
                 inv.tsq = 2.79  , 
                 .RNG.name="base::Wichmann-Hill", .RNG.seed=   4163370)
  inits3 <- list(theta   = -1.32, 
                 inv.tsq =  4.05 , 
                 .RNG.name="base::Mersenne-Twister", .RNG.seed=4163370)
  
  
  # Set-up and run simulation
  BayesMA.jags <- jags.model(file     = textConnection(indirect_effect_model), 
                             n.chains = chains, quiet   = TRUE, 
                             list(y.i = y.i, sigma.sq.i = sigma.sq.i,
                                  n   = length(y.i), 
                                  mu  = 0, eta = eta, lambda = lambda, prec = 1 / var), 
                             inits    = list(inits1, inits2, inits3))
  
  update(BayesMA.jags, burn, progress.bar="none") # Burn-in 
  
  BayesMA.sims <- coda.samples(BayesMA.jags, c("theta", "tsq"), n.iter = iterations,
                               thin = thin, progress.bar = "none")
  
  # Theta HPD credible interval for decision rule 
  all.chains <- as.mcmc(do.call(rbind, BayesMA.sims)) 
  
  # Summary measures 
  theta.median <- median(all.chains[, 1]) 
  theta.HPD.95 <- HPDinterval(all.chains)[1, ] 
  tsq.median   <- median(all.chains[, 2]) 
  tsq.HPD      <- HPDinterval(all.chains)[2, ] 
  
  
  return(list(theta  = theta.median, theta.hpd  = theta.HPD.95, 
              tau.sq = tsq.median  , tau.sq.hpd = tsq.HPD,
              GelmanConverg = gelman.diag(BayesMA.sims)))
}

### Further Setup ##
library(parallel)

## Small Effect .14 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Small5C.RData")

cl                 <- makeCluster(4)
SB_INDI_Small_5C    <- parLapply(cl, Small_5C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Small_5C, file = "SB_INDI_Small5C.RData")

rm(Small_5C   )
rm(SB_INDI_Small_5C)

####
####

# load(file = "~/Simulations/Datageneration/Small10C.RData")

cl               <- makeCluster(4)
SB_INDI_Small_10C <- parLapply(cl, Small_10C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Small_10C, file = "SB_INDI_Small10C.RData")

rm(Small_10C   )
rm(SB_INDI_Small_10C)


####
####

# load(file = "~/Simulations/Datageneration/Small30C.RData")

cl               <- makeCluster(4)
SB_INDI_Small_30C <- parLapply(cl, Small_30C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Small_30C, file = "SB_INDI_Small30C.RData")

rm(Small_30C        )
rm(SB_INDI_Small_30C)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5C.RData")

cl                 <- makeCluster(4)
SB_INDI_Smallmod_5C <- parLapply(cl, Smallmod_5C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Smallmod_5C, file = "SB_INDI_Smallmod5C.RData")

rm(Smallmod_5C        )
rm(SB_INDI_Smallmod_5C)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod10C.RData")

cl                  <- makeCluster(4)
SB_INDI_Smallmod_10C <- parLapply(cl, Smallmod_10C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Smallmod_10C, file = "SB_INDI_Smallmod10C.RData")

rm(Smallmod_10C        )
rm(SB_INDI_Smallmod_10C)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30C.RData")

cl                  <- makeCluster(4)
SB_INDI_Smallmod_30C <- parLapply(cl, Smallmod_30C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Smallmod_30C, file = "SB_INDI_Smallmod30C.RData")

rm(Smallmod_30C        )
rm(SB_INDI_Smallmod_30C)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5C.RData")

cl                  <- makeCluster(4)
SB_INDI_Moder_5C     <- parLapply(cl, Moder_5C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Moder_5C, file = "SB_INDI_Moder5C.RData")

rm(Moder_5C        )
rm(SB_INDI_Moder_5C)

####
####

# load(file = "~/Simulations/Datageneration/Moder10C.RData")

cl                   <- makeCluster(4)
clusterExport(cl,list( "Seqbay_indirect" ))
SB_INDI_Moder_10C  <- parLapply(cl, Moder_10C, function(X) try(Seqbay_indirect(X),T)) 

stopCluster(cl)


save(SB_INDI_Moder_10C, file = "SB_INDI_Moder10C.RData")

rm(Moder_10C        )
rm(SB_INDI_Moder_10C)

####
####

# load(file = "~/Simulations/Datageneration/Moder30C.RData")

cl                   <- makeCluster(4)
clusterExport(cl,list( "Seqbay_indirect" ))
SB_INDI_Moder_30C  <- parLapply(cl, Moder_30C, function(X) try(Seqbay_indirect(X),T)) 
stopCluster(cl)


save(SB_INDI_Moder_30C, file = "SB_INDI_Moder30C.RData")

rm(Moder_30C        )
rm(SB_INDI_Moder_30C)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/NoEff5C.RData")

cl                   <- makeCluster(4)
SB_INDI_NoEff_5C      <- parLapply(cl, NoEff_5C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_NoEff_5C, file = "SB_INDI_NoEff5C.RData")

rm(NoEff_5C        )
rm(SB_INDI_NoEff_5C)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10C.RData")

cl                    <- makeCluster(4)
SB_INDI_NoEff_10C      <- parLapply(cl, NoEff_10C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_NoEff_10C, file = "SB_INDI_NoEff10C.RData")

rm(NoEff_10C        )
rm(SB_INDI_NoEff_10C)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30C.RData")

cl                    <- makeCluster(4)
SB_INDI_NoEff_30C      <- parLapply(cl, NoEff_30C, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_NoEff_30C, file = "SB_INDI_NoEff30C.RData")

rm(NoEff_30C        )
rm(SB_INDI_NoEff_30C)

