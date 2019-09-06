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

# load(file = "~/Simulations/Datageneration/Small5.RData")

cl                 <- makeCluster(4)
SB_INDI_Small_5    <- parLapply(cl, Small_5, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Small_5, file = "SB_INDI_Small5.RData")

rm(Small_5   )
rm(SB_INDI_Small_5)

####
####

# load(file = "~/Simulations/Datageneration/Small10.RData")

cl               <- makeCluster(4)
SB_INDI_Small_10 <- parLapply(cl, Small_10, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Small_10, file = "SB_INDI_Small10.RData")

rm(Small_10   )
rm(SB_INDI_Small_10)


####
####

# load(file = "~/Simulations/Datageneration/Small30.RData")

cl               <- makeCluster(4)
SB_INDI_Small_30 <- parLapply(cl, Small_30, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Small_30, file = "SB_INDI_Small30.RData")

rm(Small_30        )
rm(SB_INDI_Small_30)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5.RData")

cl                 <- makeCluster(4)
SB_INDI_Smallmod_5 <- parLapply(cl, Smallmod_5, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Smallmod_5, file = "SB_INDI_Smallmod5.RData")

rm(Smallmod_5        )
rm(SB_INDI_Smallmod_5)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod10.RData")

cl                  <- makeCluster(4)
SB_INDI_Smallmod_10 <- parLapply(cl, Smallmod_10, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Smallmod_10, file = "SB_INDI_Smallmod10.RData")

rm(Smallmod_10        )
rm(SB_INDI_Smallmod_10)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30.RData")

cl                  <- makeCluster(4)
SB_INDI_Smallmod_30 <- parLapply(cl, Smallmod_30, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Smallmod_30, file = "SB_INDI_Smallmod30.RData")

rm(Smallmod_30        )
rm(SB_INDI_Smallmod_30)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5.RData")

cl                  <- makeCluster(4)
SB_INDI_Moder_5     <- parLapply(cl, Moder_5, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_Moder_5, file = "SB_INDI_Moder5.RData")

rm(Moder_5        )
rm(SB_INDI_Moder_5)

####
####

# load(file = "~/Simulations/Datageneration/Moder10.RData")

cl                   <- makeCluster(4)
clusterExport(cl,list( "Seqbay_indirect" ))
SB_INDI_Moder_10  <- parLapply(cl, Moder_10, function(X) try(Seqbay_indirect(X),T)) 

stopCluster(cl)


save(SB_INDI_Moder_10, file = "SB_INDI_Moder10.RData")

rm(Moder_10        )
rm(SB_INDI_Moder_10)

####
####

# load(file = "~/Simulations/Datageneration/Moder30.RData")

cl                   <- makeCluster(4)
clusterExport(cl,list( "Seqbay_indirect" ))
SB_INDI_Moder_30  <- parLapply(cl, Moder_30, function(X) try(Seqbay_indirect(X),T)) 
stopCluster(cl)


save(SB_INDI_Moder_30, file = "SB_INDI_Moder30.RData")

rm(Moder_30        )
rm(SB_INDI_Moder_30)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/NoEff5.RData")

cl                   <- makeCluster(4)
SB_INDI_NoEff_5      <- parLapply(cl, NoEff_5, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_NoEff_5, file = "SB_INDI_NoEff5.RData")

rm(NoEff_5        )
rm(SB_INDI_NoEff_5)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10.RData")

cl                    <- makeCluster(4)
SB_INDI_NoEff_10      <- parLapply(cl, NoEff_10, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_NoEff_10, file = "SB_INDI_NoEff10.RData")

rm(NoEff_10        )
rm(SB_INDI_NoEff_10)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30.RData")

cl                    <- makeCluster(4)
SB_INDI_NoEff_30      <- parLapply(cl, NoEff_30, Seqbay_indirect)
stopCluster(cl)


save(SB_INDI_NoEff_30, file = "SB_INDI_NoEff30.RData")

rm(NoEff_30        )
rm(SB_INDI_NoEff_30)

