# Data Generation for the simulation study: Continuous predictor 
# Camiel van Zundert

## Function for generating data ####

DataGen_C      <- function(Nstud, 
                           ABvals, 
                           Direct, 
                           Abvari, 
                           Abcor , 
                           Nx, Nsim, seedset){
  
  ## Required packages
  require(mvtnorm)
  require(lavaan)
  
  ## Computing Sigma from variances and correlations given
  R      <- diag(rep(1,2))
  R[1,2] <- Abcor
  R[2,1] <- Abcor
  D      <- diag(sqrt(Abvari))
  Sigma  <- D%*%R%*%D
  
  # Creating lists for parLapply later on (parallel computing)
  outs      <- vector(mode = "list", length = Nsim)
  
  Nstuds    <- replicate(n = Nsim,expr = Nstud, simplify = "array")
  ABvalues  <- replicate(n = Nsim,expr = ABvals, simplify = "array")
  Directs   <- replicate(n = Nsim,expr = Direct, simplify = "array")
  sigs      <- replicate(n = Nsim,expr = Sigma ,simplify = "array")
  Nxs       <- replicate(n = Nsim,expr = Nx, simplify = "array")

  
  for (i in 1:Nsim) {
    outs[[i]][1] <- Nstuds[i]
    outs[[i]][2] <- list(ABvalues[,i])
    outs[[i]][3] <- list(Directs[,i])
    outs[[i]][4] <- list(sigs[,,i])
    outs[[i]][5] <- list(Nxs[,i])
  }
  
  # Study level data generating function using lavaan
  Sim_Dat <- function(Nstudies, ab ,direct, Sigma ,nX){
    
    datagen_effect <- function(a, b, cp, nX){
      
      ntot <- nX                                       # Total participants in data set
      ab   <- a*b                                      # indirect effect
      X    <- rnorm(n = ntot)         
      M    <- a*X + sqrt(abs(1-a^2))*rnorm(n = ntot)   # generated scores of M (prediction times error)
      ey   <- 1-(cp^2 +b^2 + 2*a*cp*b)                 # Simulated error of Y
      Y    <- cp*X + b*M + ey*rnorm(n = ntot)          # Generated scores of Y
      
      Data  <- data.frame(X = X, Y = Y, M = M)
      model <- ' 
      # mediator
      M ~ 1+ a*X
      Y ~ 1+b*M +cp*X
      # indirect effect (a*b)
      ab := a*b
      # total effect
      total := cp + (a*b)
      '
      
      fit   <- sem(model, data = Data, warn = F)
      estim <- parameterestimates(fit)[parameterestimates(fit)[,4] %in% c("cp","a","b","ab", "total"),c(4,5,6)]
      
      return(list("params" = estim, 
                  "cors"   = cor(cbind(X,M,Y)),
                  "n"      = ntot))
    }
    
    ressu <- apply(rmvnorm(n     = Nstudies,
                           mean  = ab, 
                           sigma = Sigma),
                   MARGIN = 1, function(X)datagen_effect(nX = sample(x       = nX,
                                                                     size    = 1,
                                                                     replace = T) ,
                                                         a   = X[1], 
                                                         b   = X[2], 
                                                         cp  = rnorm(n    = 1,
                                                                     mean = direct[1],
                                                                     sd   = direct[2])))
    
    return(ressu)
  }
  
  # Parallel computing generation
  require(parallel)
  
  cl <- makeCluster(4)
  
  clusterSetRNGStream(cl, seedset)                          # Setting seed for reproducibility
  
  clusterExport(cl,list( "rmvnorm", "sem", "parameterestimates", "vcov" ))
  
  Simgen  <- parLapply(cl  ,outs, function(X) Sim_Dat(Nstudies = X[[1]], 
                                                      ab       = X[[2]],
                                                      direct   = X[[3]], 
                                                      Sigma    = X[[4]], 
                                                      nX       = X[[5]]))
  
  stopCluster(cl)
  return(Simgen)
}


## Small Effect .14 (Fritz and MacKinnon) ####


Small_5C <- DataGen_C(Nstud   = 5,                 # number of studies in meta-analysis
                      ABvals  = c(.14,.14),        # a and b true effect
                      Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                      Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                      Abcor   = 0,                 # 0 correlation 
                      Nx      = 50:150,            # sample size range between 50 and 150
                      Nsim    = 10000,             # number of iterations in the Simulation
                      seedset = 4163370)           # Random seed for reproducibility

save(Small_5C ,file = "Small5C.RData")
rm(Small_5C)

####
####


Small_10C <- DataGen_C(Nstud   = 10,                # number of studies in meta-analysis
                       ABvals  = c(.14,.14),        # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility


save(Small_10C,file = "Small10C.RData")
rm(Small_10C)

####
####

Small_30C <- DataGen_C(Nstud   = 30,                # number of studies in meta-analysis
                       ABvals  = c(.14,.14),        # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility


save(Small_30C, file = "Small30C.RData")
rm(Small_30C)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

Smallmod_5C <- DataGen_C(Nstud   = 5,                 # number of studies in meta-analysis
                         ABvals  = c(.26,.26),        # a and b true effect
                         Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                         Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                         Abcor   = 0,                 # 0 correlation 
                         Nx      = 50:150,            # sample size range between 50 and 150
                         Nsim    = 10000,             # number of iterations in the Simulation
                         seedset = 4163370)           # Random seed for reproducibility

save(Smallmod_5C,file = "Smallmod5C.RData")
rm(Smallmod_5C)

####
####

Smallmod_10C <- DataGen_C(Nstud   = 10,               # number of studies in meta-analysis
                          ABvals  = c(.26,.26),       # a and b true effect
                          Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                          Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                          Abcor   = 0,                 # 0 correlation 
                          Nx      = 50:150,            # sample size range between 50 and 150
                          Nsim    = 10000,             # number of iterations in the Simulation
                          seedset = 4163370)           # Random seed for reproducibility

save(Smallmod_10C,file = "Smallmod10C.RData")
rm(Smallmod_10C)

####
####

Smallmod_30C <- DataGen_C(Nstud   = 30,                # number of studies in meta-analysis
                          ABvals  = c(.26,.26),        # a and b true effect
                          Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                          Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                          Abcor   = 0,                 # 0 correlation 
                          Nx      = 50:150,            # sample size range between 50 and 150
                          Nsim    = 10000,             # number of iterations in the Simulation
                          seedset = 4163370)           # Random seed for reproducibility

save(Smallmod_30C,file = "Smallmod30C.RData")
rm(Smallmod_30C)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

Moder_5C <- DataGen_C(Nstud   = 5,                 # number of studies in meta-analysis
                      ABvals  = c(.39,.39),        # a and b true effect
                      Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                      Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                      Abcor   = 0,                 # 0 correlation 
                      Nx      = 50:150,            # sample size range between 50 and 150
                      Nsim    = 10000,             # number of iterations in the Simulation
                      seedset = 4163370)           # Random seed for reproducibility

save(Moder_5C,file = "Moder5C.RData")
rm(Moder_5C)

####
####

Moder_10C <- DataGen_C(Nstud   = 10,                # number of studies in meta-analysis
                       ABvals  = c(.39,.39),        # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility

save(Moder_10C,file = "Moder10C.RData")
rm(Moder_10C)

####
####

Moder_30C <- DataGen_C(Nstud   = 30,                # number of studies in meta-analysis
                       ABvals  = c(.39,.39),        # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility

save(Moder_30C,file = "Moder30C.RData")
rm(Moder_30C)

# No effect, for type 1 error check ####

NoEff_5C <- DataGen_C(Nstud  = 5,                  # number of studies in meta-analysis
                      ABvals = c(0,0),             # a and b true effect
                      Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                      Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                      Abcor   = 0,                 # 0 correlation 
                      Nx      = 50:150,            # sample size range between 50 and 150
                      Nsim    = 10000,             # number of iterations in the Simulation
                      seedset = 4163370)           # Random seed for reproducibility

save(NoEff_5C,file = "NoEff5C.RData")
rm(NoEff_5C)

####
####

NoEff_10C <- DataGen_C(Nstud   = 10,                # number of studies in meta-analysis
                       ABvals  = c(0,0),            # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility

save(NoEff_10C,file = "NoEff10C.RData")
rm(NoEff_10C)

####
####

NoEff_30C <- DataGen_C(Nstud   = 30,                # number of studies in meta-analysis
                       ABvals  = c(0,0),            # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility

save(NoEff_30C,file = "NoEff30C.RData")
rm(NoEff_30C)

