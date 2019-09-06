# Two Stage Structural Equation Modelling (TSSEM) for the simulation study 
# Camiel van Zundert


## Function  ####

TSSEM <- function(gen_data, REtype = "Diag"){
  
  require(metaSEM)
  
  # Input Conversion
  X <- lapply(gen_data, function(X) X[[2]])
  n <- sapply(gen_data, function(X) X[[3]])
  
  ## Regression coefficents
  A1 <- create.mxMatrix(c(       0,       0, 0,
                                 "0.1*a" ,       0, 0,
                                 "0.1*cp", "0.1*b", 0),
                        type  ="Full", 
                        byrow = TRUE  , 
                        ncol  = 3     ,
                        nrow  = 3     ,
                        as.mxMatrix = FALSE)
  
  dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("X","M","Y")
  
  ## Covariance matrix
  S1 <- create.mxMatrix(c(1,
                          0, "0.1*var_M",
                          0, 0, "0.1*var_Y"),
                        byrow=TRUE, type="Symm", as.mxMatrix=FALSE)
  
  dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("X","M","Y")
  
  # Step 1: Estimate pooled covariance matrix
  # REM: Random Effects Model
  # Diag: Random effect are independent
  random1 <- tssem1(X, 
                    n,
                    method  ="REM"  , 
                    RE.type = REtype, 
                    acov    = "weighted")
  
  random2 <- tssem2(random1, 
                    Amatrix = A1, 
                    Smatrix = S1, 
                    diag.constraints = TRUE,
                    intervals.type   = "LB",
                    model.name       = "TSSEM2 Own",
                    mx.algebras = list(ab    = mxAlgebra(a*b, name = "ab")))
  
  OG_mxstatus2 <- summary(random2)$Mx.status1
  
  if (!(OG_mxstatus2== 0 |OG_mxstatus2== 1 )){
    random2 <-  rerun(random2, checkHess = F)     # Hessian is not needed when mxAlgebra is used                                     
  }
  
  sumres <- summary(random2)
  
  return(list("indirect"        = sumres$mx.algebras[2],
              "CI"              = sumres$mx.algebras[-2],
              "OpenMX_stage1"   = summary(random1)$Mx.status1,
              "OpenMX_stage2"   = c("final" =sumres$Mx.status1,"original" = OG_mxstatus2),
              "SumInfo"         = sumres,
              "Stage1Info"      = summary(random1)))
  
}

# Further setup
library(parallel)

## Small Effect .14 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Small5.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Small_5     <- parLapply(cl, Small_5, function(X) try(TSSEM(X),TRUE)) # Try is done for catching rare errors

stopCluster(cl)


save(TSSEM_Small_5, file = "TSSEM_Small5.RData")

rm(Small_5   )
rm(TSSEM_Small_5)

####
####

# load(file = "~/Simulations/Datageneration/Small10.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Small_10     <- parLapply(cl, Small_10, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)

save(TSSEM_Small_10, file = "TSSEM_Small10.RData")

rm(Small_10   )
rm(TSSEM_Small_10)


####
####

# load(file = "~/Simulations/Datageneration/Small30.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Small_30     <- parLapply(cl, Small_30, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Small_30, file = "TSSEM_Small30.RData")

rm(Small_30        )
rm(TSSEM_Small_30)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Smallmod_5      <- parLapply(cl, Smallmod_5, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Smallmod_5, file = "TSSEM_Smallmod5.RData")

rm(Smallmod_5        )
rm(TSSEM_Smallmod_5)

####
####

# load(file = "~/Simulations/Datageneration/mallmod10.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Smallmod_10      <- parLapply(cl, Smallmod_10, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Smallmod_10, file = "TSSEM_Smallmod10.RData")

rm(Smallmod_10        )
rm(TSSEM_Smallmod_10)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Smallmod_30      <- parLapply(cl, Smallmod_30, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Smallmod_30, file = "TSSEM_Smallmod30.RData")

rm(Smallmod_30        )
rm(TSSEM_Smallmod_30)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Moder_5          <- parLapply(cl, Moder_5, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Moder_5, file = "TSSEM_Moder5.RData")

rm(Moder_5        )
rm(TSSEM_Moder_5)

####
####

# load(file = "~/Simulations/Datageneration/Moder10.RData")

cl                      <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Moder_10          <- parLapply(cl, Moder_10, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Moder_10, file = "TSSEM_Moder10.RData")

rm(Moder_10        )
rm(TSSEM_Moder_10)

####
####

# load(file = "~/Simulations/Datageneration/Moder30.RData")

cl                      <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Moder_30          <- parLapply(cl, Moder_30, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Moder_30, file = "TSSEM_Moder30.RData")

rm(Moder_30        )
rm(TSSEM_Moder_30)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/NoEff5.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_NoEff_5          <- parLapply(cl, NoEff_5, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_NoEff_5, file = "TSSEM_NoEff5.RData")

rm(NoEff_5        )
rm(TSSEM_NoEff_5)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_NoEff_10          <- parLapply(cl, NoEff_10, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_NoEff_10, file = "TSSEM_NoEff10.RData")

rm(NoEff_10        )
rm(TSSEM_NoEff_10)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_NoEff_30          <- parLapply(cl, NoEff_30, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_NoEff_30, file = "TSSEM_NoEff30.RData")

rm(NoEff_30        )
rm(TSSEM_NoEff_30)


