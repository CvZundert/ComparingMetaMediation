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

# load(file = "~/Simulations/Datageneration/Small5.RData")

cl                <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  # setting Seed

MASEM_Small_5        <- parLapply(cl, Small_5, MASEM)

stopCluster(cl)

save(MASEM_Small_5, file = "MASEM_Small5.RData")

rm(Small_5   )
rm(MASEM_Small_5)

####
####

# load(file = "~/Simulations/Datageneration/Small10.RData")

cl          <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  # setting Seed
clusterExport(cl,list( "MASEM" ))
MASEM_Small_10 <- parLapply(cl, Small_10, function(X) try(MASEM(X),T)) 

stopCluster(cl)

save(MASEM_Small_10, file = "MASEM_Small10.RData")

rm(Small_10   )
rm(MASEM_Small_10)

####
####

# load(file = "~/Simulations/Datageneration/Small30.RData")

cl           <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Small_30   <- parLapply(cl, Small_30, function(X) try(MASEM(X),T)) 

stopCluster(cl)

save(MASEM_Small_30, file = "MASEM_Small30.RData")

rm(Small_30   )
rm(MASEM_Small_30)

## Small-moderate  Effect .26 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Smallmod5.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Smallmod_5   <- parLapply(cl, Smallmod_5, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Smallmod_5, file = "MASEM_Smallmod5.RData")

rm(Smallmod_5   )
rm(MASEM_Smallmod_5)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod10.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Smallmod_10   <- parLapply(cl, Smallmod_10, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Smallmod_10, file = "MASEM_Smallmod10.RData")

rm(Smallmod_10   )
rm(MASEM_Smallmod_10)

####
####

# load(file = "~/Simulations/Datageneration/Smallmod30.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Smallmod_30   <- parLapply(cl, Smallmod_30, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Smallmod_30, file = "MASEM_Smallmod30.RData")

rm(Smallmod_30   )
rm(MASEM_Smallmod_30)

## Moderate  Effect .39 (Fritz and MacKinnon) ####

# load(file = "~/Simulations/Datageneration/Moder5.RData")

cl             <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Moder_5     <- parLapply(cl, Moder_5, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Moder_5, file = "MASEM_Moder5.RData")

rm(Moder_5   )
rm(MASEM_Moder_5)

####
####

# load(file = "~/Simulations/Datageneration/Moder10.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Moder_10     <- parLapply(cl, Moder_10, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Moder_10, file = "MASEM_Moder10.RData")

rm(Moder_10   )
rm(MASEM_Moder_10)

####
####

# load(file = "~/Simulations/Datageneration/Moder30.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_Moder_30     <- parLapply(cl, Moder_30, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_Moder_30, file = "MASEM_Moder30.RData")

rm(Moder_30   )
rm(MASEM_Moder_30)

# No effect, for type 1 error check ####

# load(file = "~/Simulations/Datageneration/NoEff5.RData")

cl              <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_NoEff_5     <- parLapply(cl, NoEff_5, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_NoEff_5, file = "MASEM_NoEff5.RData")

rm(NoEff_5   )
rm(MASEM_NoEff_5)

####
####

# load(file = "~/Simulations/Datageneration/NoEff10.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_NoEff_10     <- parLapply(cl, NoEff_10, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_NoEff_10, file = "MASEM_NoEff10.RData")

rm(NoEff_10   )
rm(MASEM_NoEff_10)

####
####

# load(file = "~/Simulations/Datageneration/NoEff30.RData")

cl               <- makeCluster(4)

clusterSetRNGStream(cl, 4163370)  
clusterExport(      cl, list( "MASEM" ))
MASEM_NoEff_30     <- parLapply(cl, NoEff_30, function(X) try(MASEM(X),T)) 

stopCluster(cl)


save(MASEM_NoEff_30, file = "MASEM_NoEff30.RData")

rm(NoEff_30   )
rm(MASEM_NoEff_30)
