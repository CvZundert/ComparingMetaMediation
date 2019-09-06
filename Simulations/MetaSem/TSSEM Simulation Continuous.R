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

# load(file = "~/Simulations/Datageneration/Small5C.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Small_5C     <- parLapply(cl, Small_5C, function(X) try(TSSEM(X),TRUE)) # Try is done for catching rare errors

stopCluster(cl)


save(TSSEM_Small_5C, file = "TSSEM_Small5C.RData")

rm(Small_5C   )
rm(TSSEM_Small_5C)

####
####

# load(file = "~/Simulations/Datageneration/Small10C.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Small_10C     <- parLapply(cl, Small_10C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)

save(TSSEM_Small_10C, file = "TSSEM_Small10C.RData")

rm(Small_10C   )
rm(TSSEM_Small_10C)


####
####

# load(file = "~/Simulations/Datageneration/Small30C.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Small_30C     <- parLapply(cl, Small_30C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Small_30C, file = "TSSEM_Small30C.RData")

rm(Small_30C        )
rm(TSSEM_Small_30C)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5C.RData")

cl                 <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Smallmod_5C      <- parLapply(cl, Smallmod_5C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Smallmod_5C, file = "TSSEM_Smallmod5C.RData")

rm(Smallmod_5C        )
rm(TSSEM_Smallmod_5C)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod10C.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Smallmod_10C      <- parLapply(cl, Smallmod_10C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Smallmod_10C, file = "TSSEM_Smallmod10C.RData")

rm(Smallmod_10C        )
rm(TSSEM_Smallmod_10C)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30C.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Smallmod_30C      <- parLapply(cl, Smallmod_30C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Smallmod_30C, file = "TSSEM_Smallmod30C.RData")

rm(Smallmod_30C        )
rm(TSSEM_Smallmod_30C)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5C.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Moder_5C          <- parLapply(cl, Moder_5C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Moder_5C, file = "TSSEM_Moder5C.RData")

rm(Moder_5C        )
rm(TSSEM_Moder_5C)

####
####

# load(file = "~/Simulations/Datageneration/Moder10C.RData")

cl                      <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Moder_10C          <- parLapply(cl, Moder_10C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Moder_10C, file = "TSSEM_Moder10C.RData")

rm(Moder_10C        )
rm(TSSEM_Moder_10C)

####
####

# load(file = "~/Simulations/Datageneration/Moder30C.RData")

cl                      <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_Moder_30C          <- parLapply(cl, Moder_30C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_Moder_30C, file = "TSSEM_Moder30C.RData")

rm(Moder_30C        )
rm(TSSEM_Moder_30C)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/NoEff5C.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_NoEff_5C          <- parLapply(cl, NoEff_5C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_NoEff_5C, file = "TSSEM_NoEff5C.RData")

rm(NoEff_5C        )
rm(TSSEM_NoEff_5C)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10C.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_NoEff_10C          <- parLapply(cl, NoEff_10C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_NoEff_10C, file = "TSSEM_NoEff10C.RData")

rm(NoEff_10C        )
rm(TSSEM_NoEff_10C)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30C.RData")

cl                     <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)        # setting Seed
clusterExport      (cl, list("TSSEM"))
TSSEM_NoEff_30C          <- parLapply(cl, NoEff_30C, function(X) try(TSSEM(X),TRUE)) 

stopCluster(cl)


save(TSSEM_NoEff_30C, file = "TSSEM_NoEff30C.RData")

rm(NoEff_30C        )
rm(TSSEM_NoEff_30C)


