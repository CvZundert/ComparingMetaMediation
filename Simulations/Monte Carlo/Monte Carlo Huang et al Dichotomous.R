# Monte Carlo (Huangs version) for the simulation study 
# Camiel van Zundert, 


## Function  ####

MonteCarlo<- function(gen_data, K = 2000, tolerance = 1e-800){
  
  require(mvtnorm)
  require(MASS)
  
  minusLogL.v <-  function(v, Xvals, sigma.aSE, sigma.bSE, correction=T){
    
    L   <- nrow(Xvals)  
    A   <- v[4]
    B   <- v[5]
    V11 <- v[1]
    V12 <- v[2]
    V22 <- v[3]
    V   <- matrix (c(V11, V12,V12, V22), nrow = 2, ncol = 2 )
    
    sigma.asCurrent <- sigma.aSE^2
    sigma.bsCurrent <- sigma.bSE^2
    
    if (correction==T){ 
      f <- (L-1)/L
      }else{
        f <- 1 
        }
    
    rslt2 <- L*log(2*pi)
    
    for (u in 1:L){
      
      M   <- diag(c(sigma.asCurrent[u], sigma.bsCurrent[u])) + V
      
      if (det(M) > 0){
        ldetM <- log (det(M))
        }else{
        cat (" ** det(M) = ", det(M) , v  ,  sigma.asCurrent[u] , " " , sigma.bsCurrent[u]  , "\n")
        ldetM <-  1000   #  make exceptionally large to force det to be positive
        }
      
      rslt2   <- rslt2 + 1/2*f*ldetM + 1/2 * t((Xvals[u,]-c(A,B))) %*% solve (M, tol = tolerance ) %*% (Xvals[u,]-c(A,B))
    }
    
    return (rslt2)
  }
  
  VtoW <- function(v){
    
    w1 <- log(v[1] + 0.001)
    w2 <- v[2] / sqrt(v[1] + 0.001)
    w3 <- log(v[3] - (v[2]^2/(v[1]+0.001)))
    
    if (length(v) == 3)
      return (c(w1, w2, w3))
    
    else 
      return (c(w1, w2, w3 ,v[4:5]))
  }
  
  WtoV <- function(w){
    
    V11 <- exp(w[1])
    V12 <- w[2]*sqrt(V11)
    V22 <- (exp(w[3])) + (V12^2/exp(w[1]))
    
    if(length(w) == 3)
      return(c(V11, V12, V22))
    else
      return(c(V11, V12, V22, w[4:5]))
  }
  
  minusLogL.w <-  function(w, Xvals, sigma.aSE, sigma.bSE){
    
    v   <- WtoV(w)
    ans <- minusLogL.v(v, Xvals, sigma.aSE, sigma.bSE)  
    
    return(ans)
  }
  
  
  gradient.v <- function (v, Xvals, sigma.aSE, sigma.bSE, correction=T){
    
    K <- nrow(Xvals)
    V <- matrix (c(v[1], v[2], v[2], v[3]), ncol = 2)
    A <- v[4]
    B <- v[5]
    
    dsigma.v11 <- matrix (c(1, 0, 0, 0), ncol = 2)
    dsigma.v12 <- matrix (c(0, 1, 1, 0), ncol = 2)
    dsigma.v22 <- matrix (c(0, 0, 0, 1), ncol = 2)
    
    grv11 <- 0
    grv12 <- 0
    grv22 <- 0
    grA   <- 0
    grB   <- 0
    
    if (correction==T)        
    {  f <- (K -1)/K}
    
    else{ f <- 1}
    
    
    for (i in 1:K){
      
      M1         <- matrix(c(sigma.aSE[i]^2, 0, 0, sigma.bSE[i]^2 ), nr=2) + V			
      M1Inv      <- solve(M1,tol = tolerance ) 
      ssq        <- (Xvals[i,] - c(A,B)) %*% t(Xvals [i,] - c(A,B))           ## 2 x 2 SSQ matrix
      DiffMatrix <- ssq %*% M1Inv - f* diag(c(1, 1))
      
      grv11 <- grv11 - 1/2 * sum(diag(DiffMatrix %*% dsigma.v11 %*%  M1Inv))     #formula 11 for d( - R Log L )
      
      grv12 <- grv12 - 1/2 * sum(diag(DiffMatrix %*% dsigma.v12 %*%  M1Inv))    
      
      grv22 <- grv22 - 1/2 * sum(diag(DiffMatrix %*% dsigma.v22 %*%  M1Inv))   
      
      Value <- M1Inv %*% (Xvals[i,] - c(A, B))     #formula 10
      grA   <- grA -  Value[1,1]   
      grB   <- grB -  Value[2,1]   
    }
    
    return(c(grv11, grv12, grv22, grA, grB))
    
  }
  
  
  gradient.w <-function(w, Xvals, sigma.aSE, sigma.bSE){
    
    v       <- WtoV(w)
    grad.v  <- gradient.v(v, Xvals, sigma.aSE, sigma.bSE)
    deriv.w <- t(dVdW(w)) %*% grad.v   
    
    return(deriv.w)
  }
  
  dVdW <- function (w){
    
    # get derivs of V wrt W
    
    dSigma11dw <- c(exp(w[1]), 0, 0, 0, 0)                                   #  = exp ( w[1])
    dSigma12dw <- c( 1/2 * w[2] * exp(1/2 * w[1]), exp(1/2 * w[1]), 0, 0, 0) #  = w[2] * e^ 1/2 w[1] 
    dSigma22dw <- c(0, 2 * w[2], exp(w[3]), 0, 0)                            #  = exp( w[3] ) + w[2]^2 
    
    dvdw <- rbind (dSigma11dw, 
                   dSigma12dw, 
                   dSigma22dw, 
                   c(0, 0, 0, 1, 0), 
                   c(0, 0, 0, 0, 1))
    return( dvdw)
  }
  
  dWdV <- function(v){
    
    v22.1 <- v[3] - v[2]^2/v[1]
    dw1dv <- c(1/v[1], 0, 0, 0, 0)
    dw2dv <- c(- 1/2 * v[2]/v[1]^(1.5), v[1]^(-1/2), 0, 0, 0)
    dw3dv <- c(1/v22.1 * v[2]^2/v[1]^2, - 2/v22.1 * v[2]/v[1], 1/v22.1, 0, 0)
    dw4dv <- c(0, 0, 0, 1, 0)
    dw5dv <- c(0, 0, 0, 0, 1)
    
    dwdv <- rbind(dw1dv, dw2dv, dw3dv, dw4dv, dw5dv)
    
    return ( dwdv)
  }
  
  
  MLsolution.w <-  function(Xvals, sigma.aSE, sigma.bSE, SMALL =0.0001, Bignegative=-10){
    
    L           <- nrow(Xvals)
    sigma.asVar <- sigma.aSE^2
    sigma.bsVar <- sigma.bSE^2
    
    t    <- c(max(var(Xvals[,1]) - mean(sigma.asVar), SMALL), 0, max(var(Xvals[,2]) - mean(sigma.bsVar), SMALL))
    
    w11s <- VtoW(t)[1]
    w12s <-	0
    w22s <- VtoW(t)[3]
    
    a1   <- mean(Xvals[,1])
    b1   <- mean(Xvals[,2])
    
    # use optim function for ML estimates of between trial level variance and covariance, A and B
    
   
    answer <- optim(par       = c(w11s, w12s, w22s, a1, b1), 
                    fn        = minusLogL.w, 
                    gr        = gradient.w, 
                    method    = c("L-BFGS-B"), 
                    lower     = c(Bignegative, -Inf, Bignegative, -Inf, -Inf), 
                    control   = list(trace=6, REPORT=1), 
                    hessian   = TRUE, 
                    Xvals     = Xvals,  
                    sigma.aSE = sigma.aSE, 
                    sigma.bSE = sigma.bSE)
    
    return(answer)
  }
  
  summary.w <- function(rslt.w, Xvals, sigma.aSE, sigma.bSE){
    
    estimate.w <- rslt.w$par  ## on transformed scale (w)
    deriv.w    <- gradient.w(estimate.w, Xvals, sigma.aSE, sigma.bSE)
    
    # convert choleski parameters to var-cov
    
    estimate.v           <- WtoV(estimate.w)
    names(estimate.v)    <- c("v11", "v12", "v22", "A", "B")
    deriv.v              <- gradient.v(estimate.v, Xvals, sigma.aSE, sigma.bSE)
    dvdw                 <- dVdW(estimate.w)
    variancematrix.v     <- dvdw %*% solve(rslt.w$hessian, tol = tolerance) %*% t(dvdw)  # for Sigma and means
    
    
    dimnames(variancematrix.v) <- list(c("v11", "v12", "v22", "A", "B"), 
                                       c("v11", "v12", "v22", "A", "B"))
    
    return(list(convergence      = rslt.w$convergence, 
                minusRlogL       = rslt.w$value, 
                estimate.w       = estimate.w, 
                deriv.w          = deriv.w, 
                estimate.v       = estimate.v, 
                deriv.v          = deriv.v,
                variancematrix.v = variancematrix.v))
  }
  
  
  # Extracting information from generation
  Xvalues <- cbind(a = (sapply(gen_data, function(X) X[[1]][1,2])),
                   b = (sapply(gen_data, function(X) X[[1]][2,2])))
  SE_A    <- sapply(gen_data, function(X) X[[1]][1,3])
  SE_B    <- sapply(gen_data, function(X) X[[1]][2,3])
  
  # Running the above functions
  
  Esties  <- MLsolution.w(Xvals     = Xvalues,
                          sigma.aSE = SE_A,
                          sigma.bSE = SE_B)
  
  results <- summary.w(rslt.w    = Esties,
                       Xvals     = Xvalues,
                       sigma.aSE = SE_A,
                       sigma.bSE = SE_B )
  
  
  setsam                <- array(data = NA, dim = c(nrow(Xvalues), 2, K))
  MeanAtimesB.star      <- rep(0,K)   
  AtimesBsigma.hat.star <- rep(0,K) 
  
  for (u in 1:K) {
    
    setsam[,,u]              <- mvrnorm(n     = nrow(Xvalues), 
                                        mu    = results$estimate.v[4:5], 
                                        Sigma = results$variancematrix.v[4:5,4:5], tol = tolerance)
    
    MeanAtimesB.star[u]      <- mean(setsam[,1,u]*setsam[,2,u])
    
    AtimesBsigma.hat.star[u] <- sqrt((mean(setsam[,1,u])^2)*var(setsam[,2,u])+
                                     (mean(setsam[,2,u])^2)*var(setsam[,1,u]))
    
  } 
  
  A    <- results$estimate.v[4]
  B    <- results$estimate.v[5]
  
  ABse <- sqrt(((A^2)*(results$variancematrix.v[5,5]))+((B^2)*(results$variancematrix.v[4,4]))) #sobel
  
  
  CL       <- c((A*B) - ABse*quantile ((MeanAtimesB.star - (A*B))/AtimesBsigma.hat.star, probs = c(.975), names=FALSE),      
                (A*B) - ABse*quantile ((MeanAtimesB.star - (A*B))/AtimesBsigma.hat.star, probs = c(.025), names=FALSE))
  
  names(CL)<- c("2.5%", "97.5%")
  
  return(list("A_B"      = results$estimate.v[4:5], 
              "indirect" = prod(results$estimate.v[4:5]),
              "vcov"     = results$variancematrix.v[4:5,4:5],
              "CI"       = CL ))
}


### Further Setup ##
library(parallel)

## Small Effect .14 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Small5.RData")

cl                <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  # setting Seed

MC_Small_5        <- parLapply(cl, Small_5, MonteCarlo)

stopCluster(cl)

save(MC_Small_5, file = "MC_Small5.RData")

rm(Small_5   )
rm(MC_Small_5)

####
####

# load(file = "~/Simulations/Datageneration/Small10.RData")

cl          <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  # setting Seed
clusterExport(cl,list( "MonteCarlo" ))
MC_Small_10 <- parLapply(cl, Small_10, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)

save(MC_Small_10, file = "MC_Small10.RData")

rm(Small_10   )
rm(MC_Small_10)

####
####

# load(file = "~/Simulations/Datageneration/Small30.RData")

cl           <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_Small_30   <- parLapply(cl, Small_30, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)

save(MC_Small_30, file = "MC_Small30.RData")

rm(Small_30   )
rm(MC_Small_30)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_Smallmod_5   <- parLapply(cl, Smallmod_5, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_Smallmod_5, file = "MC_Smallmod5.RData")

rm(Smallmod_5   )
rm(MC_Smallmod_5)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod10.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_Smallmod_10   <- parLapply(cl, Smallmod_10, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_Smallmod_10, file = "MC_Smallmod10.RData")

rm(Smallmod_10   )
rm(MC_Smallmod_10)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_Smallmod_30   <- parLapply(cl, Smallmod_30, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_Smallmod_30, file = "MC_Smallmod30.RData")

rm(Smallmod_30   )
rm(MC_Smallmod_30)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_Moder_5     <- parLapply(cl, Moder_5, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_Moder_5, file = "MC_Moder5.RData")

rm(Moder_5   )
rm(MC_Moder_5)

####
####

# load(file = "~/Simulations/Datageneration/Moder10.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_Moder_10     <- parLapply(cl, Moder_10, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_Moder_10, file = "MC_Moder10.RData")

rm(Moder_10   )
rm(MC_Moder_10)

####
####

# load(file = "~/Simulations/Datageneration/Moder30.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_Moder_30     <- parLapply(cl, Moder_30, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_Moder_30, file = "MC_Moder30.RData")

rm(Moder_30   )
rm(MC_Moder_30)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/NoEff5.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_NoEff_5     <- parLapply(cl, NoEff_5, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_NoEff_5, file = "MC_NoEff5.RData")

rm(NoEff_5   )
rm(MC_NoEff_5)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_NoEff_10     <- parLapply(cl, NoEff_10, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_NoEff_10, file = "MC_NoEff10.RData")

rm(NoEff_10   )
rm(MC_NoEff_10)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MonteCarlo" ))
MC_NoEff_30     <- parLapply(cl, NoEff_30, function(X) try(MonteCarlo(X),T)) 

stopCluster(cl)


save(MC_NoEff_30, file = "MC_NoEff30.RData")

rm(NoEff_30   )
rm(MC_NoEff_30)

