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

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Small5C.RData")

cl                <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Small_5C  <- parLapply(cl, MC_Small_5C, function(X) MC_CI_cor(X, 5))

stopCluster(cl)

save(MC_CiCor_Small_5C, file = "MC_CiCor_Small5C.RData")

rm(MC_Small_5C      )
rm(MC_CiCor_Small_5C)

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Small10C.RData")

cl          <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Small_10C <- parLapply(cl, MC_Small_10C, function(X) MC_CI_cor(X, 10))

stopCluster(cl)

save(MC_CiCor_Small_10C, file = "MC_CiCor_Small_10C.RData")

rm(MC_Small_10C         )
rm(MC_CiCor_Small_10C   )

load("MC_CiCor_Small_10C.RData")

save(MC_CiCor_Small_10C, file = "MC_CiCor_Small10C.RData")
####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Small30C.RData")

cl           <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Small_30C <- parLapply(cl, MC_Small_30C, function(X) MC_CI_cor(X, 30))

stopCluster(cl)

save(MC_CiCor_Small_30C, file = "MC_CiCor_Small30C.RData")


rm(MC_Small_30C         )
rm(MC_CiCor_Small_30C   )

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Smallmod5C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Smallmod_5C <- parLapply(cl, MC_Smallmod_5C[-551], function(X) MC_CI_cor(X, 5))

stopCluster(cl)

save(MC_CiCor_Smallmod_5C, file = "MC_CiCor_Smallmod5C.RData")

rm(MC_CiCor_Smallmod_5C   )
rm(MC_Smallmod_5C         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Smallmod10C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Smallmod_10C <- parLapply(cl, MC_Smallmod_10C[-734], function(X) MC_CI_cor(X, 10))

stopCluster(cl)


save(MC_CiCor_Smallmod_10C, file = "MC_CiCor_Smallmod10C.RData")

rm(MC_CiCor_Smallmod_10C   )
rm(MC_Smallmod_10C         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Smallmod30C.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Smallmod_30C <- parLapply(cl, MC_Smallmod_30C, function(X) MC_CI_cor(X, 30))

stopCluster(cl)


save(MC_CiCor_Smallmod_30C, file = "MC_CiCor_Smallmod30C.RData")

rm(MC_CiCor_Smallmod_30C   )
rm(MC_Smallmod_30C         )

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Moder5C.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Moder_5C <- parLapply(cl, MC_Moder_5C, function(X) MC_CI_cor(X, 5))

stopCluster(cl)

save(MC_CiCor_Moder_5C, file = "MC_CiCor_Moder5C.RData")

rm(MC_CiCor_Moder_5C   )
rm(MC_Moder_5C         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Moder10C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Moder_10C <- parLapply(cl, MC_Moder_10C[-c(1882, 2749, 3639, 3953)], function(X) MC_CI_cor(X, 10))

stopCluster(cl)


save(MC_CiCor_Moder_10C, file = "MC_CiCor_Moder10C.RData")

rm(MC_CiCor_Moder_10C   )
rm(MC_Moder_10C         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_Moder30C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_Moder_30C <- parLapply(cl, MC_Moder_30C, function(X) MC_CI_cor(X, 30))

stopCluster(cl)


save(MC_CiCor_Moder_30C, file = "MC_CiCor_Moder30C.RData")

rm(MC_CiCor_Moder_30C   )
rm(MC_Moder_30C         )

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_NoEff5C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_NoEff_5C <- parLapply(cl, MC_NoEff_5C, function(X) MC_CI_cor(X, 5))

stopCluster(cl)


save(MC_CiCor_NoEff_5C, file = "MC_CiCor_NoEff5C.RData")

rm(MC_CiCor_NoEff_5C   )
rm(MC_NoEff_5C         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_NoEff10C.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_NoEff_10C <- parLapply(cl, MC_NoEff_10C, function(X) MC_CI_cor(X, 10))

stopCluster(cl)


save(MC_CiCor_NoEff_10C, file = "MC_CiCor_NoEff10C.RData")

rm(MC_CiCor_NoEff_10C   )
rm(MC_NoEff_10C         )

####
####

# load(file = "~/Simulations/Datageneration/Monte Carlo/MC_NoEff30C.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370     )  # setting Seed
clusterExport      (cl, "MC_CI_cor" )

MC_CiCor_NoEff_30C <- parLapply(cl, MC_NoEff_30C, function(X) MC_CI_cor(X, 30))

stopCluster(cl)


save(MC_CiCor_NoEff_30C, file = "MC_CiCor_NoEff30C.RData")

rm(MC_CiCor_NoEff_30C   )
rm(MC_NoEff_30C         )




