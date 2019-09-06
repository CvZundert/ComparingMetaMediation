# Sequential Bayes (Two - Parameter) for the simulation study 
# Camiel van Zundert


## Function  ####
MultiBayes <- function(gen_data, var = 5,
                       iterations = 5000, chains = 3, 
                       burn = 1000, thin = 1){
  
  require(rjags)
  
  # Input conversion
  yi_a       <- (sapply(gen_data, function(X)  X[[1]][1,2]))
  yi_b       <- (sapply(gen_data, function(X)  X[[1]][2,2]))
  
  sigma_a      <- (sapply(gen_data, function(X) X[[1]][1,3]**2)) 
  # Since inverse of variance, not SE
  sigma_b      <- (sapply(gen_data, function(X) X[[1]][2,3]**2)) 
  # Since inverse of variance, not SE
  
  
  Multi_inv <- "
  model { 
  for (i in 1:n) { 
  y.a.i[i]        ~ dnorm(theta.a.i[i], 1 / sigma.a.sq.i[i])
  y.b.i[i]        ~ dnorm(theta.b.i[i], 1 / sigma.b.sq.i[i])
  
  theta.a.i[i] ~ dnorm(theta.a, inv.tsq.a)
  theta.b.i[i] ~ dnorm(theta.b, inv.tsq.b)
  } # JAGS uses precision not variance 
  
  theta.a      ~ dnorm(mu, prec) 
  theta.b      ~ dnorm(mu, prec) 
  inv.tsq.a    ~ dgamma(1.5, 0.08)      # Appropriate prior for precision 
  inv.tsq.b    ~ dgamma(1.5, 0.08)  
  tsq.a        <- 1 / inv.tsq.a 
  tsq.b        <- 1 / inv.tsq.b
  
  indirect     <- theta.a*theta.b
  }
  "
  
  # Initial values for each chain, random seed selected
  
  inits1 <- list(theta.a   = -0.92 ,         # Random number from   runif(1, -2, 2)
                 inv.tsq.a = 12.19,          # Random number from 1/runif(1,  0, 0.5)
                 theta.b   = -0.33,          # Random number from   runif(1, -2, 2)  
                 inv.tsq.b =  6.04 ,         # Random number from 1/runif(1,  0, 0.5)
                 .RNG.name="base::Super-Duper", .RNG.seed=     4163370) # seed
  inits2 <- list(theta.a   = 0.60 ,  
                 inv.tsq.a = 22.0,        
                 theta.b  =  0.58,         
                 inv.tsq.b =  2.68 ,  
                 .RNG.name="base::Wichmann-Hill", .RNG.seed=   4163370)
  inits3 <- list(theta.a   = 0.34 ,  
                 inv.tsq.a = 11.6,        
                 theta.b   = 0.26 ,         
                 inv.tsq.b = 16.83 ,   
                 .RNG.name="base::Mersenne-Twister", .RNG.seed=4163370)
  
  # Set-up and run simulation
  BayesMA.jags <- jags.model(file     = textConnection(Multi_inv), 
                             n.chains = chains, quiet   = TRUE, 
                             list(y.a.i = yi_a, 
                                  y.b.i = yi_b,
                                  sigma.a.sq.i = sigma_a,
                                  sigma.b.sq.i = sigma_b,
                                  n   = length(yi_a), 
                                  mu  = 0,  prec = 1 / 5), 
                             inits    = list(inits1, inits2, inits3))
  
  update(BayesMA.jags, burn, progress.bar="none") # Burn-in 
  
  BayesMA.sims <- coda.samples(BayesMA.jags, c("indirect", 
                                               "tsq.a", "tsq.b",
                                               "theta.a" , "theta.b")
                               , n.iter = iterations,
                               thin = thin, progress.bar = "none")
  
  # Theta HPD credible interval for decision rule 
  all.chains <- as.mcmc(do.call(rbind, BayesMA.sims)) 
  
  # Summary measures 
  theta.median  <- median(all.chains[, 1]) 
  theta.HPD.95  <- HPDinterval(all.chains)[1, ] 
  
  tsqa.median   <- median(all.chains[, 2]) 
  tsqa.HPD      <- HPDinterval(all.chains)[2, ] 
  
  tsqb.median   <- median(all.chains[, 3]) 
  tsqb.HPD      <- HPDinterval(all.chains)[3, ] 

  a.median      <- median(all.chains[, 4]) 
  a.HPD         <- HPDinterval(all.chains)[4, ] 
  
  b.median      <- median(all.chains[, 5]) 
  b.HPD         <- HPDinterval(all.chains)[5, ] 
  
  return(list(theta      = theta.median , theta.hpd    = theta.HPD.95, 
              tau.a.sq   = tsqa.median  , tau.a.sq.hpd = tsqa.HPD,
              tau.b.sq   = tsqb.median  , tau.b.sq.hpd = tsqb.HPD,
              theta.a.sq = a.median     , theta.a.hpd  = a.HPD,
              theta.b.sq = b.median     , theta.b.hpd  = b.HPD,
              GelmanConverg = gelman.diag(BayesMA.sims)))
}

### Further Setup ##
library(parallel)
setwd("~/Research Archive_ Camiel van Zundert_4163370/Simulations/Sequential Bayesian/Dichotomous/SB Two Para Effect results")

## Small Effect .14 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Small5.RData")

cl                 <- makeCluster(4)
SB_MULT_Small_5    <- parLapply(cl, Small_5, MultiBayes)
stopCluster(cl)


save(SB_MULT_Small_5, file = "SB_MULT_Small5.RData")

rm(Small_5   )
rm(SB_MULT_Small_5)

####
####

# load(file = "~/Simulations/Datageneration/Small10.RData")

cl               <- makeCluster(4)
SB_MULT_Small_10 <- parLapply(cl, Small_10, MultiBayes)
stopCluster(cl)


save(SB_MULT_Small_10, file = "SB_MULT_Small10.RData")

rm(Small_10   )
rm(SB_MULT_Small_10)


####
####

# load(file = "~/Simulations/Datageneration/Small30.RData")

cl               <- makeCluster(4)
SB_MULT_Small_30 <- parLapply(cl, Small_30, MultiBayes)
stopCluster(cl)


save(SB_MULT_Small_30, file = "SB_MULT_Small30.RData")

rm(Small_30        )
rm(SB_MULT_Small_30)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5.RData")

cl                 <- makeCluster(4)
SB_MULT_Smallmod_5 <- parLapply(cl, Smallmod_5, MultiBayes)
stopCluster(cl)


save(SB_MULT_Smallmod_5, file = "SB_MULT_Smallmod5.RData")

rm(Smallmod_5        )
rm(SB_MULT_Smallmod_5)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod10.RData")

cl                  <- makeCluster(4)
SB_MULT_Smallmod_10 <- parLapply(cl, Smallmod_10, MultiBayes)
stopCluster(cl)


save(SB_MULT_Smallmod_10, file = "SB_MULT_Smallmod10.RData")

rm(Smallmod_10        )
rm(SB_MULT_Smallmod_10)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30.RData")

cl                  <- makeCluster(4)
SB_MULT_Smallmod_30 <- parLapply(cl, Smallmod_30, MultiBayes)
stopCluster(cl)


save(SB_MULT_Smallmod_30, file = "SB_MULT_Smallmod30.RData")

rm(Smallmod_30        )
rm(SB_MULT_Smallmod_30)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5.RData")

cl                  <- makeCluster(4)
SB_MULT_Moder_5     <- parLapply(cl, Moder_5, MultiBayes)
stopCluster(cl)


save(SB_MULT_Moder_5, file = "SB_MULT_Moder5.RData")

rm(Moder_5        )
rm(SB_MULT_Moder_5)

####
####

# load(file = "~/Simulations/Datageneration/Moder10.RData")

cl                   <- makeCluster(4)
clusterExport(cl,list( "MultiBayes" ))
SB_MULT_Moder_10  <- parLapply(cl, Moder_10, function(X) try(MultiBayes(X),T)) 

stopCluster(cl)


save(SB_MULT_Moder_10, file = "SB_MULT_Moder10.RData")

rm(Moder_10        )
rm(SB_MULT_Moder_10)

####
####

# load(file = "~/Simulations/Datageneration/Moder30.RData")

cl                   <- makeCluster(4)
clusterExport(cl,list( "MultiBayes" ))
SB_MULT_Moder_30  <- parLapply(cl, Moder_30, function(X) try(MultiBayes(X),T)) 
stopCluster(cl)


save(SB_MULT_Moder_30, file = "SB_MULT_Moder30.RData")

rm(Moder_30        )
rm(SB_MULT_Moder_30)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/NoEff5.RData")

cl                   <- makeCluster(4)
SB_MULT_NoEff_5      <- parLapply(cl, NoEff_5, MultiBayes)
stopCluster(cl)


save(SB_MULT_NoEff_5, file = "SB_MULT_NoEff5.RData")

rm(NoEff_5        )
rm(SB_MULT_NoEff_5)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10.RData")

cl                    <- makeCluster(4)
SB_MULT_NoEff_10      <- parLapply(cl, NoEff_10, MultiBayes)
stopCluster(cl)


save(SB_MULT_NoEff_10, file = "SB_MULT_NoEff10.RData")

rm(NoEff_10        )
rm(SB_MULT_NoEff_10)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30.RData")

cl                    <- makeCluster(4)
SB_MULT_NoEff_30      <- parLapply(cl, NoEff_30, MultiBayes)
stopCluster(cl)


save(SB_MULT_NoEff_30, file = "SB_MULT_NoEff30.RData")

rm(NoEff_30        )
rm(SB_MULT_NoEff_30)


