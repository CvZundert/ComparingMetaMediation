# Data Generation for the simulation study 
# Camiel van Zundert, 

## Function for generating data ####

DataGen      <- function(Nstud, 
                        ABvals, 
                        Direct, 
                        Abvari, 
                        Abcor , 
                        Nx, Nprop, Nsim, seedset){
  
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
  Nprops    <- replicate(n = Nsim,expr = Nprop, simplify = "array")
  
  for (i in 1:Nsim) {
    outs[[i]][1] <- Nstuds[i]
    outs[[i]][2] <- list(ABvalues[,i])
    outs[[i]][3] <- list(Directs[,i])
    outs[[i]][4] <- list(sigs[,,i])
    outs[[i]][5] <- list(Nxs[,i])
    outs[[i]][6] <- list(Nprops[,i])
    
  }
  
  # Study level data generating function using lavaan
  Sim_Dat <- function(Nstudies, ab ,direct, Sigma ,nX, Nprop){
    
    datagen_effect <- function(a, b, cp, nX, Nprop){
      
      ntot <- nX                                       # Total participants in data set
      ab   <- a*b                                      # indirect effect
      X    <- c(rep(0, floor(nX*Nprop)),               # Number of participants in each condition
                rep(1, ceiling(nX*(1-Nprop))))       
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
                  "n"      = ntot,
                  "nX"     = table(X)))
    }
    
    ressu <- apply(rmvnorm(n     = Nstudies,
                           mean  = ab, 
                           sigma = Sigma),
                   MARGIN = 1, function(X)datagen_effect(nX = sample(x       = nX,
                                                                     size    = 1,
                                                                     replace = T) , 
                                                         Nprop = runif(n   = 1, 
                                                                       min = Nprop[1],
                                                                       max = Nprop[2]),
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
                                                      nX       = X[[5]],
                                                      Nprop    = X[[6]]))
  
  stopCluster(cl)
  return(Simgen)
}


## Small Effect .14 (Fritz and MacKinnon) ####


Small_5 <- DataGen(Nstud   = 5,                 # number of studies in meta-analysis
                   ABvals  = c(.14,.14),        # a and b true effect
                   Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                   Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                   Abcor   = 0,                 # 0 correlation 
                   Nx      = 50:150,            # sample size range between 50 and 150
                   Nprop   = c(.40,.60),        # distribution of participants across conditions
                   Nsim    = 10000,             # number of iterations in the Simulation
                   seedset = 4163370)           # Random seed for reproducibility

save(Small_5,file = "Small5.RData")
rm(Small_5)

####
####


Small_10 <- DataGen(Nstud   = 10,                # number of studies in meta-analysis
                    ABvals  = c(.14,.14),        # a and b true effect
                    Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                    Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                    Abcor   = 0,                 # 0 correlation 
                    Nx      = 50:150,            # sample size range between 50 and 150
                    Nprop   = c(.40,.60),        # distribution of participants across conditions
                    Nsim    = 10000,             # number of iterations in the Simulation
                    seedset = 4163370)           # Random seed for reproducibility


save(Small_10,file = "Small10.RData")
rm(Small_10)

####
####

Small_30 <- DataGen(Nstud   = 30,                # number of studies in meta-analysis
                    ABvals  = c(.14,.14),        # a and b true effect
                    Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                    Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                    Abcor   = 0,                 # 0 correlation 
                    Nx      = 50:150,            # sample size range between 50 and 150
                    Nprop   = c(.40,.60),        # distribution of participants across conditions
                    Nsim    = 10000,             # number of iterations in the Simulation
                    seedset = 4163370)           # Random seed for reproducibility


save(Small_30,file = "Small30.RData")
rm(Small_30)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

Smallmod_5 <- DataGen(Nstud   = 5,                 # number of studies in meta-analysis
                      ABvals  = c(.26,.26),        # a and b true effect
                      Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                      Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                      Abcor   = 0,                 # 0 correlation 
                      Nx      = 50:150,            # sample size range between 50 and 150
                      Nprop   = c(.40,.60),        # distribution of participants across conditions
                      Nsim    = 10000,             # number of iterations in the Simulation
                      seedset = 4163370)           # Random seed for reproducibility

save(Smallmod_5,file = "Smallmod5.RData")
rm(Smallmod_5)

####
####

Smallmod_10 <- DataGen(Nstud   = 10,               # number of studies in meta-analysis
                       ABvals  = c(.26,.26),       # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nprop   = c(.40,.60),        # distribution of participants across conditions
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility

save(Smallmod_10,file = "Smallmod10.RData")
rm(Smallmod_10)

####
####

Smallmod_30 <- DataGen(Nstud   = 30,                # number of studies in meta-analysis
                       ABvals  = c(.26,.26),        # a and b true effect
                       Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                       Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                       Abcor   = 0,                 # 0 correlation 
                       Nx      = 50:150,            # sample size range between 50 and 150
                       Nprop   = c(.40,.60),        # distribution of participants across conditions
                       Nsim    = 10000,             # number of iterations in the Simulation
                       seedset = 4163370)           # Random seed for reproducibility

save(Smallmod_30,file = "Smallmod30.RData")
rm(Smallmod_30)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

set.seed(999018674)
Moder_5 <- DataGen(Nstud   = 5,                 # number of studies in meta-analysis
                   ABvals  = c(.39,.39),        # a and b true effect
                   Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                   Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                   Abcor   = 0,                 # 0 correlation 
                   Nx      = 50:150,            # sample size range between 50 and 150
                   Nprop   = c(.40,.60),        # distribution of participants across conditions
                   Nsim    = 10000,             # number of iterations in the Simulation
                   seedset = 4163370)           # Random seed for reproducibility

save(Moder_5,file = "Moder5.RData")
rm(Moder_5)

####
####

Moder_10 <- DataGen(Nstud   = 10,                # number of studies in meta-analysis
                    ABvals  = c(.39,.39),        # a and b true effect
                    Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                    Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                    Abcor   = 0,                 # 0 correlation 
                    Nx      = 50:150,            # sample size range between 50 and 150
                    Nprop   = c(.40,.60),        # distribution of participants across conditions
                    Nsim    = 10000,             # number of iterations in the Simulation
                    seedset = 4163370)           # Random seed for reproducibility

save(Moder_10,file = "Moder10.RData")
rm(Moder_10)

####
####

Moder_30 <- DataGen(Nstud   = 30,                # number of studies in meta-analysis
                    ABvals  = c(.39,.39),        # a and b true effect
                    Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                    Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                    Abcor   = 0,                 # 0 correlation 
                    Nx      = 50:150,            # sample size range between 50 and 150
                    Nprop   = c(.40,.60),        # distribution of participants across conditions
                    Nsim    = 10000,             # number of iterations in the Simulation
                    seedset = 4163370)           # Random seed for reproducibility

save(Moder_30,file = "Moder30.RData")
rm(Moder_30)

# No effect, for type 1 error check ####

NoEff_5 <- DataGen(Nstud  = 5,                  # number of studies in meta-analysis
                   ABvals = c(0,0),             # a and b true effect
                   Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                   Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                   Abcor   = 0,                 # 0 correlation 
                   Nx      = 50:150,            # sample size range between 50 and 150
                   Nprop   = c(.40,.60),        # distribution of participants across conditions
                   Nsim    = 10000,             # number of iterations in the Simulation
                   seedset = 4163370)           # Random seed for reproducibility

save(NoEff_5,file = "NoEff5.RData")
rm(NoEff_5)

####
####

NoEff_10 <- DataGen(Nstud   = 10,                # number of studies in meta-analysis
                    ABvals  = c(0,0),            # a and b true effect
                    Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                    Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                    Abcor   = 0,                 # 0 correlation 
                    Nx      = 50:150,            # sample size range between 50 and 150
                    Nprop   = c(.40,.60),        # distribution of participants across conditions
                    Nsim    = 10000,             # number of iterations in the Simulation
                    seedset = 4163370)           # Random seed for reproducibility

save(NoEff_10,file = "NoEff10.RData")
rm(NoEff_10)

####
####

NoEff_30 <- DataGen(Nstud   = 30,                # number of studies in meta-analysis
                    ABvals  = c(0,0),            # a and b true effect
                    Direct  = c(0, sqrt(.0219)), # 0 = true direct effect (sd = 0.15)
                    Abvari  = c(.0219,.0219),    # heterogeneity of a and b
                    Abcor   = 0,                 # 0 correlation 
                    Nx      = 50:150,            # sample size range between 50 and 150
                    Nprop   = c(.40,.60),        # distribution of participants across conditions
                    Nsim    = 10000,             # number of iterations in the Simulation
                    seedset = 4163370)           # Random seed for reproducibility

save(NoEff_30,file = "NoEff30.RData")
rm(NoEff_30)

