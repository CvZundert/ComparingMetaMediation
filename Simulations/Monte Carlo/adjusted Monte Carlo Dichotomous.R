# MC_correction function ####

# Note! This method only changes the Confidence Interval Estimation of Monte Carlo (Huang), and is dependent on the
# Estimates of this method. Run the Monte Carlo Huang methods first!


MC_CI_cor <- function(MC_output, Nobs_X, K= 2000){  
  
  meanAB                <- MC_output[[1]]
  sigmaAB               <- MC_output[[3]]   
  setsam                <- array(data = NA, dim = c(Nobs_X, 2, K))
  MeanAtimesB.star      <- rep(0,K)   
  AtimesBsigma.hat.star <- rep(0,K) 

  require(MASS)
  
  for (u in 1:K) {
    setsam[,,u]              <- mvrnorm(n = Nobs_X, mu  = meanAB, Sigma = sigmaAB )
    MeanAtimesB.star[u]      <- mean(setsam[,1,u] * setsam[,2,u])
    AtimesBsigma.hat.star[u] <- sqrt((mean(setsam[,1,u])^2) * (var(setsam[,2,u])/Nobs_X)+
                                     (mean(setsam[,2,u])^2) * (var(setsam[,1,u])/Nobs_X))
  } 
  
  A     <- meanAB[1]
  B     <- meanAB[2]
  Pivot <- (MeanAtimesB.star - (A*B))/ AtimesBsigma.hat.star
  
  ABse  <- sqrt(((A^2)*(sigmaAB[2,2]))+((B^2)*(sigmaAB[1,1]))) 
  #sobel (could also add covariance!!!)
  ## Or t(meanAB)%*% sigmaAB %*% meanAB
  
  CL      <- c((A*B) - ABse * quantile(Pivot, probs = c(.975), names=FALSE),      
               (A*B) - ABse * quantile(Pivot, probs = c(.025), names=FALSE) )
  
  names(CL)<- c("2.5%", "97.5%")
  
  return(list("A_B"      = meanAB, 
              "indirect" = prod(meanAB),
              "vcov"     = sigmaAB,
              "CI"       = CL ))
}

### Further Setup ##
library(parallel)


## Small Effect .14 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Small5.RData")

cl                <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Small_5  <- parLapply(cl, MC_Small_5, function(X) MC_CI_cor(X, 5))

stopCluster(cl)

save(MC_CiCor_Small_5, file = "MC_CiCor_Small5.RData")

rm(MC_Small_5      )
rm(MC_CiCor_Small_5)

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Small10.RData")

cl          <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Small_10 <- parLapply(cl, MC_Small_10[-1031], function(X) MC_CI_cor(X, 10))

stopCluster(cl)

save(MC_CiCor_Small_10, file = "MC_CiCor_Small_10.RData")

rm(MC_Small_10         )
rm(MC_CiCor_Small_10   )

load("MC_CiCor_Small_10.RData")

save(MC_CiCor_Small_10, file = "MC_CiCor_Small10.RData")
####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Small30.RData")

cl           <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Small_30 <- parLapply(cl, MC_Small_30, function(X) MC_CI_cor(X, 30))

stopCluster(cl)

save(MC_CiCor_Small_30, file = "MC_CiCor_Small30.RData")


rm(MC_Small_30         )
rm(MC_CiCor_Small_30   )

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Smallmod5.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Smallmod__5 <- parLapply(cl, MC_Smallmod_5, function(X) MC_CI_cor(X, 5))

stopCluster(cl)

MC_CiCor_Smallmod_5<- MC_CiCor_Smallmod__5
save(MC_CiCor_Smallmod_5, file = "MC_CiCor_Smallmod5.RData")

rm(MC_CiCor_Smallmod_5   )
rm(MC_Smallmod_5         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Smallmod10.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Smallmod_10 <- parLapply(cl, MC_Smallmod_10, function(X) MC_CI_cor(X, 10))

stopCluster(cl)


save(MC_CiCor_Smallmod_10, file = "MC_CiCor_Smallmod10.RData")

rm(MC_CiCor_Smallmod_10   )
rm(MC_Smallmod_10         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Smallmod30.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Smallmod_30 <- parLapply(cl, MC_Smallmod_30, function(X) MC_CI_cor(X, 30))

stopCluster(cl)


save(MC_CiCor_Smallmod_30, file = "MC_CiCor_Smallmod30.RData")

rm(MC_CiCor_Smallmod_30   )
rm(MC_Smallmod_30         )

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Moder5.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Moder_5 <- parLapply(cl, MC_Moder_5[-c(2749, 4597)], function(X) MC_CI_cor(X, 5))

stopCluster(cl)

save(MC_CiCor_Moder_5, file = "MC_CiCor_Moder5.RData")

rm(MC_CiCor_Moder_5   )
rm(MC_Moder_5         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Moder10.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Moder_10 <- parLapply(cl, MC_Moder_10[-5744], function(X) MC_CI_cor(X, 10))

stopCluster(cl)


save(MC_CiCor_Moder_10, file = "MC_CiCor_Moder10.RData")

rm(MC_CiCor_Moder_10   )
rm(MC_Moder_10         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Moder30.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Moder_30 <- parLapply(cl, MC_Moder_30[-c(2117, 4764, 9399)], function(X) MC_CI_cor(X, 30))

stopCluster(cl)


save(MC_CiCor_Moder_30, file = "MC_CiCor_Moder30.RData")

rm(MC_CiCor_Moder_30   )
rm(MC_Moder_30         )

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_NoEff5.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_NoEff_5<- parLapply(cl, MC_NoEff_5[-c(3850, 4399)], function(X) MC_CI_cor(X, 5))

stopCluster(cl)


save(MC_CiCor_NoEff_5, file = "MC_CiCor_NoEff5.RData")

rm(MC_CiCor_NoEff_5   )
rm(MC_NoEff_5)

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_NoEff10.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_NoEff_10 <- parLapply(cl, MC_NoEff_10, function(X) MC_CI_cor(X, 10))

stopCluster(cl)


save(MC_CiCor_NoEff_10, file = "MC_CiCor_NoEff10.RData")

rm(MC_CiCor_NoEff_10   )
rm(MC_NoEff_10         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_NoEff30.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_NoEff_30 <- parLapply(cl, MC_NoEff_30, function(X) MC_CI_cor(X, 30))

stopCluster(cl)


save(MC_CiCor_NoEff_30, file = "MC_CiCor_NoEff30.RData")

rm(MC_CiCor_NoEff_30   )
rm(MC_NoEff_30         )




