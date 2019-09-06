# Meta Analytic Univariate Effect examination using MASEM for the simulation study 
# Camiel van Zundert


## Function  ####

MASEM <- function(gen_data){
  
  require(metaSEM)
  
  # Input conversion
  y.i        <- sapply(gen_data, function(X) X$`params`[X$`params`[,"label"]=="ab", "est"])    # y.i is study level indirect effect 
  sigma.sq.i <- sapply(gen_data, function(X) X$`params`[X$`params`[,"label"]=="ab", "se"]**2)    # sigma.sq.i is se of indirect effect 
  
  meta_res   <- meta(y = y.i, v = sigma.sq.i, intervals.type = "LB") # Random Effects (variance is not constrained to be equal)
  
  sum_res    <- summary(meta_res)
  
  return(list("indirect"    = sum_res$coefficients[1,1],
              "CI"          = sum_res$coefficients[1,3:4],
              "Mx.status"   = sum_res$Mx.status1,
              "summaryInfo" = sum_res))
}



### Further Setup ##
library(parallel)

## Small Effect .14 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Small5C.RData")

cl                <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  # setting Seed

MASEM_Small_5C        <- parLapply(cl, Small_5C, MASEM)

stopCluster(cl)

save(MASEM_Small_5C, file = "MASEM_Small5C.RData")

rm(Small_5C   )
rm(MASEM_Small_5C)

####
####

# load(file = "~/Simulations/Datageneration/Small10C.RData")

cl          <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  # setting Seed
clusterExport(cl,list( "MASEM" ))
MASEM_Small_10C <- parLapply(cl, Small_10C, function(X) try(MASEM(X),T)) 

stopCluster(cl)

save(MASEM_Small_10C, file = "MASEM_Small10C.RData")

rm(Small_10C   )
rm(MASEM_Small_10C)

####
####

# load(file = "~/Simulations/Datageneration/Small30C.RData")

cl           <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Small_30C   <- parLapply(cl, Small_30C, function(X) try(MASEM(X),T)) 

stopCluster(cl)

save(MASEM_Small_30C, file = "MASEM_Small30C.RData")

rm(Small_30C   )
rm(MASEM_Small_30C)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Smallmod_5C   <- parLapply(cl, Smallmod_5C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Smallmod_5C, file = "MASEM_Smallmod5C.RData")

rm(Smallmod_5C   )
rm(MASEM_Smallmod_5C)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod10C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Smallmod_10C   <- parLapply(cl, Smallmod_10C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Smallmod_10C, file = "MASEM_Smallmod10C.RData")

rm(Smallmod_10C   )
rm(MASEM_Smallmod_10C)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30C.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Smallmod_30C   <- parLapply(cl, Smallmod_30C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Smallmod_30C, file = "MASEM_Smallmod30C.RData")

rm(Smallmod_30C   )
rm(MASEM_Smallmod_30C)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5C.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Moder_5C     <- parLapply(cl, Moder_5C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Moder_5C, file = "MASEM_Moder5C.RData")

rm(Moder_5C   )
rm(MASEM_Moder_5C)

####
####

# load(file = "~/Simulations/Datageneration/Moder10C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Moder_10C     <- parLapply(cl, Moder_10C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Moder_10C, file = "MASEM_Moder10C.RData")

rm(Moder_10C   )
rm(MASEM_Moder_10C)

####
####

# load(file = "~/Simulations/Datageneration/Moder30C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Moder_30C     <- parLapply(cl, Moder_30C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Moder_30C, file = "MASEM_Moder30C.RData")

rm(Moder_30C   )
rm(MASEM_Moder_30C)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/oEff5C.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_NoEff_5C     <- parLapply(cl, NoEff_5C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_NoEff_5C, file = "MASEM_NoEff5C.RData")

rm(NoEff_5C   )
rm(MASEM_NoEff_5C)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10C.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_NoEff_10C     <- parLapply(cl, NoEff_10C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_NoEff_10C, file = "MASEM_NoEff10C.RData")

rm(NoEff_10C   )
rm(MASEM_NoEff_10C)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30C.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_NoEff_30C     <- parLapply(cl, NoEff_30C, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_NoEff_30C, file = "MASEM_NoEff30C.RData")

rm(NoEff_30C   )
rm(MASEM_NoEff_30C)
