#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Argument 1 is replicate
Rep <- args[1]

# Argument 2 is burnin vs scenario
burninOrScenario <- args[2]

if (!burninOrScenario %in% c("burnin", "scenario")) {
  stop("Argument 2 must be either burnin or scenario!")
}

# Argument 3 is parameter 1
if (burninOrScenario == "burnin") {
  param1Value <- NA
} else {
  param1Value <- args[3]
}

# Argument 4 is parameter 2
if (burninOrScenario == "burnin") {
  param2Value <- NA
} else {
  param2Value <- args[4]
}
if (burninOrScenario == "burnin") {
#code for burn-in

dir <- paste0("Rep_", Rep)
if (!dir.exists(dir)) {
  stop(paste0("Rep folder is missing for Rep ", Rep, "!"))
}
setwd(dir)
dire <- paste0("stopimp_", param1Value, "_param2_", param2Value)
unlink(dir)
dir.create(path = dire)
tula<- args[3]
setwd("phenoscenario_param1_0.04_param2_NA")
load(paste0("Selection",Rep,"_import0.04.RData"))
setwd("..")
setwd(dire)
param1Value<-as.numeric(tula)
pImport<-param1Value
library(AlphaSimR)
library(ggplot2)
#library(tictoc)
library(R6)
library(nadiv)
library(Matrix)
library(SIMplyBee)
library(dplyr)
#library(ggpubr)
library(tidyr)
getwd()
#This has its own ssh file
year=21
nYear=30
for (year in 21:nYear) {
  print("Starting the cycle")
  #year <- 1 (Use this to check that things are working without setting the whole for loop off )
  #year <- year + 1
  cat(paste0("Year: ", year, "/", nYear, "\n"))
  
  # If this is the first year, create some colonies to start with
  if (year == 1) {
    print("Creating initial colonies")
    age1 <- list(Mel = createMultiColony(x = queens$Mel, n = IrelandSize),
                 Car = createMultiColony(x = queens$Car, n = CarSize))
    
    # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
  } else {
    age2 <- list(Mel = age1$Mel, Car = age1$Car) #, Lig = age1$Lig
    age1 <- list(Mel = age0$Mel, Car = age0$Car) #, Lig = age0$Lig
    age0 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
    age0p1 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
    age0p2 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
  }
  
  #This if we want to set the colonies in a specific location for all reps from getLocations.R
  #we would need to load the x<- load(equis.RData) and the y<-load(y.RData)
  #if (year==11){
  # idmel<-getId(age1$Mel)
  
  #equisde<-x[1:length(idmel)]
  #ygreca<-y[1:length(idmel)]
  #  age1$Mel <-  setLocation(age1$Mel, 
  #                  location = Map(c, equisde, ygreca))
  # idmel2<-getId(age2$Mel)
  # equisde2<-x[length(idmel):length(idmel2)]
  # ygreca2<-y[length(idmel):length(idmel2)]
  #age2$Mel <-  setLocation(age2$Mel, 
  #                        location = Map(c, equisde2, ygreca2))
  #}
  
  
  
  if (year==11){
    #Este codigo coge grupos de colonias entre 2-7 y le asigna una localizacon random a los grupos de colonias.
    idmel<-getId(age1$Mel)
    apiary <- list()
    min_colonies <- 4
    max_colonies <- 10
    while (length(idmel) > 0) {
      apiary_size <- sample(min_colonies:max_colonies, 1)
      
      if (apiary_size > length(idmel)) {
        apiary_size <- length(idmel)
      }
      apiary <- c(apiary, list(idmel[1:apiary_size]))
      idmel <- idmel[-(1:apiary_size)]
    }
    
    a<- pullColonies(age1$Mel, n = length(apiary))
    a<-a$pulled
    a <-  setLocation(a, 
                      location = Map(c, runif(nColonies(a), 0, mapsize), runif(nColonies(a), 0, mapsize)))
    
    k <- unname(sapply(getLocation(a), function(X) X[1]))
    h <- unname(sapply(getLocation(a), function(X) X[2]))
    
    b <- c()
    for (i in 1:length(apiary)) {
      grupo <- apiary[[i]]
      num_ids <- length(grupo)
      b <- c(b, rep(k[i], num_ids))
    }
    d<-c()
    for (i in 1:length(apiary)) {
      grupo <- apiary[[i]]
      num_ids <- length(grupo)
      d <- c(d, rep(h[i], num_ids))
    }
    
    age1$Mel<-setLocation(age1$Mel, location = Map(c, abs((b+runif(length(b), min = -0.01, max = 0.01))), abs((d+runif(length(d), min = -0.01, max = 0.01)))))
    
    #age 2
    idmel2<-getId(age2$Mel)
    apiary2 <- list()
    min_colonies2 <- 4
    max_colonies2 <- 10
    while (length(idmel2) > 0) {
      apiary_size2 <- sample(min_colonies2:max_colonies2, 1)
      
      if (apiary_size2 > length(idmel2)) {
        apiary_size2 <- length(idmel2)
      }
      apiary2 <- c(apiary2, list(idmel2[1:apiary_size2]))
      idmel2 <- idmel2[-(1:apiary_size2)]
    }
    
    a2<- pullColonies(age2$Mel, n = length(apiary2))
    a2<-a2$pulled
    a2 <-  setLocation(a2, 
                       location = Map(c, runif(nColonies(a2), 0, mapsize), runif(nColonies(a2), 0, mapsize)))
    
    k2 <- unname(sapply(getLocation(a2), function(X) X[1]))
    h2 <- unname(sapply(getLocation(a2), function(X) X[2]))
    
    b2 <- c()
    for (i in 1:length(apiary2)) {
      grupo2 <- apiary2[[i]]
      num_ids2 <- length(grupo2)
      b2 <- c(b2, rep(k2[i], num_ids2))
    }
    d2<-c()
    for (i in 1:length(apiary2)) {
      grupo2 <- apiary2[[i]]
      num_ids2 <- length(grupo2)
      d2 <- c(d2, rep(h2[i], num_ids2))
    }
    
    age2$Mel<-setLocation(age2$Mel, location = Map(c, abs((b2+runif(length(b2), min = -0.01, max = 0.01))), abs((d2+runif(length(d2), min = -0.01, max = 0.01)))))
    
    
  }
  # Period1 ------------------------------------------------------------------
  # Build-up the colonies
  print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
  print(Sys.time())
  age1 <- list(Mel = buildUp(age1$Mel),
               Car = buildUp(age1$Car))
  
  if (year > 1) {
    age2 <- list(Mel = buildUp(age2$Mel),
                 Car = buildUp(age2$Car))
    
  }
  
  # Split all age1 colonies
  print("Splitting the colonies")
  print(Sys.time())
  tmp <- list(Mel = split(age1$Mel),
              Car = split(age1$Car))
  
  age1 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  
  # The queens of the splits are 0 years old
  age0p1 <- list(Mel = tmp$Mel$split, Car = tmp$Car$split) 
  x <- sapply(getLocation(age1$Mel), function(X) X[1])
  y <- sapply(getLocation(age1$Mel), function(X) X[2])
  age0p1$Mel<-setLocation(age0p1$Mel, location = Map(c, x,y))
  if (year > 1) {
    
    tmp <- list(Mel = split(age2$Mel),
                Car = split(age2$Car))
    
    age2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant
    )
    #Set the location of splits to a location near where the original colony is
    a <- sapply(getLocation(age2$Mel), function(X) X[1])
    b <- sapply(getLocation(age2$Mel), function(X) X[2])
    tmp$Mel$split<-setLocation(tmp$Mel$split, location = Map(c, a,b))
    # The queens of the splits are 0 years old
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                   Car = c(age0p1$Car, tmp$Car$split))
  }
  
  
  # Create virgin queens
  # Sample colony for the virgin queens
  print("Create virgin queens, period 1")
  print(Sys.time())
  
  # Virgin queens for splits!
  
  
  virginQueens <- list(Mel = createVirginQueens(age1$Mel, collapse=T, nInd=10),
                       Car = createVirginQueens(age1$Car, collapse=T, nInd=10))
  
  virginQueens<-list(Mel = mergePops(virginQueens$Mel),
                     Car = mergePops(virginQueens$Car))
  
  #requeen with carnica and mellifera for Mel and carnica for Car
  nColoniesMel<-nColonies(age0p1$Mel) 
  nColoniesCar<-nColonies(age0p1$Car)
  age0p1 <- list(Mel = reQueen(age0p1$Mel, queen = virginQueens$Mel[sample(1:nInd(virginQueens$Mel), size = nColoniesMel, replace = FALSE)]),
                 Car = reQueen(age0p1$Car, queen = virginQueens$Car[sample(1:nInd(virginQueens$Car), size = nColoniesCar, replace = FALSE)]))
  
  # Swarm a percentage of age1 colonies
  print("Swarm colonies, P1")
  print(Sys.time())
  tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
              Car = pullColonies(age1$Car, p = p1swarm))
  
  age1 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 2),
              Car = swarm(tmp$Car$pulled))
  
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                 Car = c(age0p1$Car, tmp$Car$remnant))
  
  age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
               Car = c(age1$Car, tmp$Car$swarm))
  
  
  
  if (year > 1) {
    # Swarm a percentage of age2 colonies
    tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
                Car = pullColonies(age2$Car, p = p1swarm))
    
    age2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 2),
                Car = swarm(tmp$Car$pulled))
    
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                   Car = c(age0p1$Car, tmp$Car$remnant))
    
    age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                 Car = c(age2$Car, tmp$Car$swarm))
    
  }
  
  # Supersede age1 colonies
  print("Supersede colonies, P1")
  print(Sys.time())
  tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
              Car = pullColonies(age1$Car, p = p1supersede))
  
  age1 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  tmp <- list(Mel = supersede(tmp$Mel$pulled),
              Car = supersede(tmp$Car$pulled))
  
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                 Car = c(age0p1$Car, tmp$Car))
  
  
  if (year > 1) {
    # Supersede age2 colonies
    tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
                Car = pullColonies(age2$Car, p = p1supersede))
    
    age2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                Car = supersede(tmp$Car$pulled))
    
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                   Car = c(age0p1$Car, tmp$Car))
    
  }
  
  # Mate the split colonies
  print("Mate split colonies, P1")
  print(Sys.time())
  if (year == 1) {
    age0p1$Mel <- cross(age0p1$Mel, droneColonies = age1$Mel, crossPlan= "create", spatial= T, radius= 15, nDrones= nDronesPoisson, checkCross = "warning")
    DCACar <- createDCA(age1$Car)
    age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
  } else {
    age0p1$Mel <- cross(age0p1$Mel, droneColonies = c(age1$Mel,age2$Mel), crossPlan= "create", spatial= T, radius= 15, nDrones= nDronesPoisson, checkCross = "warning")
    DCACar <- createDCA(c(age1$Car, age2$Car))
    age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
  }
  dronesped1<-getPed(DCAMel)
  dronespedCar1<-getPed(DCACar)
  # Collapse
  print("Collapse colonies, P1")
  print(Sys.time())
  age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
               Car = selectColonies(age1$Car, p = 1 - p1collapse))
  
  if (year > 1) {
    age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
                 Car = selectColonies(age2$Car, p = 1 - p1collapse))
    
  }
  
  # Period2 ------------------------------------------------------------------
  print("PERIOD 2")
  # Swarm a percentage of age1 colonies
  # Mellifera
  print("Swarm colonies, P2")
  print(Sys.time())
  tmp <- list(Mel = pullColonies(age1$Mel, p = p2swarm),
              Car = pullColonies(age1$Car, p = p2swarm))
  
  age1 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 2),
              Car = swarm(tmp$Car$pulled))
  
  # The queens of the remnant colonies are of age 0
  age0p2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
  
  age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
               Car = c(age1$Car, tmp$Car$swarm))
  
  
  if (year > 1) {
    # Swarm a percentage of age2 colonies
    tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
                Car = pullColonies(age2$Car, p = p2swarm))
    
    age2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 2),
                Car = swarm(tmp$Car$pulled))
    
    # The queens of the remnant colonies are of age 0
    age0p2 <- list(Mel = tmp$Mel$remnant,
                   Car = tmp$Car$remnant)
    
    age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                 Car = c(age2$Car, tmp$Car$swarm))
    
  }
  
  # Supersede a part of age1 colonies
  print("Supersede colonies, P2")
  print(Sys.time())
  
  tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
              Car = pullColonies(age1$Car, p = p2supersede))
  age1 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  tmp <- list(Mel = supersede(tmp$Mel$pulled),
              Car = supersede(tmp$Car$pulled))
  # The queens of superseded colonies are of age 0
  age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                 Car = c(age0p2$Car, tmp$Car))
  
  
  if (year > 1) {
    # Supersede a part of age2 colonies
    tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
                Car = pullColonies(age2$Car, p = p2supersede))
    age2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                Car = supersede(tmp$Car$pulled))
    # The queens of superseded colonies are of age 0
    age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                   Car = c(age0p2$Car, tmp$Car))
    
  }
  
  # Replace all the drones
  print("Replace Drones, P2")
  print(Sys.time())
  
  age1$Mel <- replaceDrones(age1$Mel)
  age1$Car <- replaceDrones(age1$Car)
  
  if (year > 1) {
    age2$Mel <- replaceDrones(age2$Mel)
    age2$Car <- replaceDrones(age2$Car)
  }
  
  # Mate the colonies
  # Import p percentage of carnica colonies into mellifera DCA
  print("Mate colonies, P2")
  print(Sys.time())
  
  
  if (year == 1) {
    age0p2$Mel <- cross(age0p2$Mel, droneColonies = age1$Mel, crossPlan= "create", spatial= T, radius= 15, nDrones= nDronesPoisson, checkCross = "warning")
    DCACar <- createDCA(age1$Car)
    age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
  } else {
    age0p2$Mel <- cross(age0p2$Mel, droneColonies = c(age1$Mel,age2$Mel), crossPlan= "create", spatial= T, radius= 15, nDrones= nDronesPoisson, checkCross = "warning")
    DCACar <- createDCA(c(age1$Car, age2$Car))
    fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
    age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
  }
  
  dronesped2<-getPed(DCAMel)
  dron<-rbind(dronesped1,dronesped2)
  dronesped <- rbind(dronesped,dron)
  dronespedCar2<-getPed(DCACar)
  dronCar<-rbind(dronespedCar1,dronespedCar2)
  dronespedCar <- rbind(dronespedCar,dronCar)
  
  # Collapse
  age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p2collapse),
               Car = selectColonies(age1$Car, p = 1 - p2collapse))
  if (year > 1) {
    age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
                 Car = selectColonies(age2$Car, p = 1 - p2collapse))
  }
  
  # Merge all age 0 colonies (from both periods)
  age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
               Car = c(age0p1$Car, age0p2$Car))
  
  
  # Period3 ------------------------------------------------------------------
  # Collapse age0 queens
  print("PERIOD 3")
  print("Collapse colonies, P3")
  print(Sys.time())
  
  #selection on Fitness queen gv
  failcross<-nVirginQueens(age0$Mel)
  FVQ <- names(which(failcross == 1))
  tmp <- pullColonies(age0$Mel, ID = FVQ)
  age0$Mel <- tmp$remnant
  
  #age0
  #Mellifera
  MelBritFit<-sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,3])#fitness values of mel in Ireland
  NativeFitnessOptima<-mean(MelBritFit) #mean fitness of the native population
  squaredDeviation4Native = function(x) (x - NativeFitnessOptima)^2
  deviations <- sapply(MelBritFit, squaredDeviation4Native) #calculate the deviations
  queensID <- names(sort(deviations)) #sort the deviations from less deviating to most deviating
  Nselectcolon<-round(length(queensID)*(1-p3collapseAge0)) #calculate how many colonies wont collapse
  age0MelqueensID<-queensID[1:Nselectcolon] #select the queens ids that will not collapse
  
  MelBritFit1<-sapply(getGv(age1$Mel, caste = "queen"), function(x) x[1,3])#fitness values of mel in Ireland
  NativeFitnessOptima1<-mean(MelBritFit1) #mean fitness of the native population
  squaredDeviation4Native1 = function(x) (x - NativeFitnessOptima1)^2
  deviations1 <- sapply(MelBritFit1, squaredDeviation4Native1) #calculate the deviations
  queensID1 <- names(sort(deviations1)) #sort the deviations from less deviating to most deviating
  Nselectcolon1<-round(length(queensID1)*(1-p3collapseAge0)) #calculate how many colonies wont collapse
  age1MelqueensID<-queensID1[1:Nselectcolon1] #select the queens ids that will not collapse
  
  CarEuFit<-sapply(getGv(age0$Car, caste = "queen"), function(x) x[1,4])#fitness values of mel in Ireland
  NonNativeFitnessOptima<-mean(CarEuFit) #mean fitness of the native population
  squaredDeviation4NonNative = function(x) (x - NonNativeFitnessOptima)^2
  deviations2 <- sapply(CarEuFit, squaredDeviation4NonNative) #calculate the deviations
  queensID2 <- names(sort(deviations2)) #sort the deviations from less deviating to most deviating
  Nselectcolon2<-round(length(queensID2)*(1-p3collapseAge0)) #calculate how many colonies wont collapse
  age0CarqueensID<-queensID2[1:Nselectcolon2] #select the queens ids that will not collapse
  
  CarEuFit1<-sapply(getGv(age1$Car, caste = "queen"), function(x) x[1,4])#fitness values of mel in Ireland
  NonNativeFitnessOptima1<-mean(CarEuFit1) #mean fitness of the native population
  squaredDeviation4NonNative1 = function(x) (x - NonNativeFitnessOptima1)^2
  deviations3 <- sapply(CarEuFit1, squaredDeviation4NonNative1) #calculate the deviations
  queensID3 <- names(sort(deviations3)) #sort the deviations from less deviating to most deviating
  Nselectcolon3<-round(length(queensID3)*(1-p3collapseAge0)) #calculate how many colonies wont collapse
  age1CarqueensID<-queensID3[1:Nselectcolon3] #select the queens ids that will not collapse
  
  
  #Collapse 
  
  age0 <- list(Mel = selectColonies(age0$Mel, ID = age0MelqueensID),
               Car = selectColonies(age0$Car, ID = age0CarqueensID))
  
  age1 <- list(Mel = selectColonies(age1$Mel, ID = age1MelqueensID),
               Car = selectColonies(age1$Car, ID= age1CarqueensID))
  
  age2 <- list(Mel = NULL, Car = NULL) #We don't need this but just to show the workflow!!!
  
  
  # Maintain the number of colonies ------------------------------------------
  # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
  print("Maintain the number, P2")
  print(Sys.time())
  
  age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
  age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)
  
  
  for (subspecies in c("Mel", "Car")) {     
    if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == IrelandSize
        | (nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == CarSize)
    {
    } 
    else 
    {stop(paste0("The number of colonies for ", subspecies, " does not match the population size!"))}
  }
  
  
  
  
  #option2 of euc distances to be able to put it in the same dataframe
  
  loc<-data.frame(Location= getLocation((age1$Mel), collapse=T))
  x<-75
  y<-75
  has<-c(x,y)
  IdAge <- as.numeric(rownames(loc))
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  EuclAge <- numeric(nrow(loc))
  
  for (i in 1:nrow(loc)){
    EuclAge[i]<-euclidean(has,loc[i,1:2])
  }
  queens = mergePops(getQueen(age1$Mel))
  IBDh<-apply(getIbdHaplo(queens),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
  IBD<-sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
  eucyear<-rep(year, length(IBD))
  age<-rep(1, length(IBD))
  if (year==11){
    Eucdf<-data.frame(IdAge,EuclAge, IBD,eucyear,age)
  } else {
    Eucdfnow<-data.frame(IdAge,EuclAge, IBD,eucyear,age)
    Eucdf<-rbind(Eucdf,Eucdfnow)
  }
  
  loc<-data.frame(Location= getLocation((age0$Mel), collapse=T))
  IdAge <- as.numeric(rownames(loc))
  EuclAge <- numeric(nrow(loc))
  for (i in 1:nrow(loc)){
    EuclAge[i]<-euclidean(has,loc[i,1:2])
  }
  queens = mergePops(getQueen(age0$Mel))
  IBDh<-apply(getIbdHaplo(queens),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
  IBD<-sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
  eucyear<-rep(year, length(IBD))
  age<-rep(0, length(IBD))
  if (year==11){
    Eucdf2<-data.frame(IdAge,EuclAge,IBD,eucyear,age)
  } else {
    Eucdfnow2<-data.frame(IdAge,EuclAge,IBD,eucyear,age)
    Eucdf2<-rbind(Eucdf2,Eucdfnow2)
  }
  
  
  #Alphapart
  correctfathers<-rbind(dronesped,dronespedCar) #We need the pedigree of dornes of pop2 because when I import a honey bee from pop2, I do get ped, and I obtain the mother and the drone, but we need to substitute the drone for the mother of the drone in the father column. And the mother of that drone is in the pedigree of drones of carnica
  E<-getPed(mergePops(getQueen(age0$Mel)))
  catufo<-pullColonies(age0$Mel, ID=IdImportColonies)
  catufo<-catufo$pulled
  lengthcatufo<-length(getQueen(catufo))
  if (lengthcatufo>0){
    import<-getPed(mergePops(getQueen(catufo)))
    Importids<-import$id
    E <- E %>%
      left_join(correctfathers, by = c("father" = "id")) %>%
      select(id, mother = mother.x, father, motheroffather = mother.y) #again I select the mother of drones to be fathers
    E <- E %>%
      distinct(id, .keep_all = TRUE)
    gen<-year+2
    Alphapart4<-data.frame(Generation=rep(gen, length(E$id)),
                           IId=E$id,
                           FId=E$motheroffather,
                           MId=E$mother,
                           Population=rep("Pop1", length(E$id)),
                           BvFi=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,3]),
                           BvFc=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,4]),
                           BvHi=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,1]),
                           BvHc=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,2]))
    
    Alphapart4$Population[Alphapart4$IId %in% Importids] <- "Pop2"
  } else{
    
    E <- E %>%
      left_join(correctfathers, by = c("father" = "id")) %>%
      select(id, mother = mother.x, father, motheroffather = mother.y) #again I select the mother of drones to be fathers
    E <- E %>%
      distinct(id, .keep_all = TRUE)
    gen<-year+2
    Alphapart4<-data.frame(Generation=rep(gen, length(E$id)),
                           IId=E$id,
                           FId=E$motheroffather,
                           MId=E$mother,
                           Population=rep("Pop1", length(E$id)),
                           BvFi=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,3]),
                           BvFc=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,4]),
                           BvHi=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,1]),
                           BvHc=sapply(getGv(age0$Mel, caste = "queen"), function(x) x[1,2]))
    
  }
  Alphapart<-rbind(Alphapart,Alphapart4)
  
  correctfathersCar<-dronespedCar
  Ec<-getPed(mergePops(getQueen(age0$Car)))
  Ec <- Ec %>%
    left_join(correctfathersCar, by = c("father" = "id")) %>%
    select(id, mother = mother.x, father, motheroffather = mother.y) #again I select the mother of drones to be fathers
  gen<-year+2
  AlphapartCar4<-data.frame(Generation=rep(gen, length(Ec$id)),
                            IId=Ec$id,
                            FId=Ec$motheroffather,
                            MId=Ec$mother,
                            Population=rep("Pop2", length(Ec$id)),
                            BvFi=sapply(getGv(age0$Car, caste = "queen"), function(x) x[1,3]),
                            BvFc=sapply(getGv(age0$Car, caste = "queen"), function(x) x[1,4]),
                            BvHi=sapply(getGv(age0$Car, caste = "queen"), function(x) x[1,1]),
                            BvHc=sapply(getGv(age0$Car, caste = "queen"), function(x) x[1,2]))
  AlphapartCar<-rbind(AlphapartCar,AlphapartCar4)
  
  AlphapartCombined<-rbind(AlphapartCombined,Alphapart4,AlphapartCar4)
  
  #this creates two empty dataframes where all the information of each year and each replica will be recorded
  
  
  columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","sdHoneyYieldBrit","HoneyYieldEu","sdHoneyYieldEu",
                   "FitnessBrit","sdFitnessBrit","FitnessEu","sdFitnessEu","pheHoneyYieldBrit","sdpheHoneyYieldBrit","pheHoneyYieldEu","sdpheHoneyYieldEu",
                   "pheFitnessBrit","sdpheFitnessBrit","pheFitnessEu","sdpheFitnessEu","Homocigosity","sdHomocigosity","sdIBD","survivingCar","Eucdist")
  
  if (year==1){ 
    columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","sdHoneyYieldBrit","HoneyYieldEu","sdHoneyYieldEu",
                     "FitnessBrit","sdFitnessBrit","FitnessEu","sdFitnessEu","pheHoneyYieldBrit","sdpheHoneyYieldBrit","pheHoneyYieldEu","sdpheHoneyYieldEu",
                     "pheFitnessBrit","sdpheFitnessBrit","pheFitnessEu","sdpheFitnessEu","Homocigosity","sdHomocigosity","sdIBD","survivingCar","Eucdist")
    MeanVarMel <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for mellifera
    colnames(MeanVarMel)<-columnheaders
    MeanVarCar <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for carnica
    colnames(MeanVarCar)<-columnheaders
  } 
  
  if (year==1){
    uno<-0
  }else{
    uno<-nrow(colonyRecords) #this is to show where to start recording values on the dataframe
  }
  
  
  #record values of Mellifera population
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel", Rep=Rep)
  unoymedio<-nrow(colonyRecords)
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Mel, year = year, population = "Mel", Rep=Rep)
  dos<-nrow(colonyRecords)#this is to show where to end recording values on the dataframe
  
  
  
  #create a dataframe with the mellifera mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
  newrow1<- data.frame(MeanIBD=mean(colonyRecords[(uno+1):dos,"IBD"]), VarIBD=var(colonyRecords[(uno+1):dos,"IBD"]), Year=year,Rep=Rep,Population="Mel",
                       HoneyYieldBrit=mean(colonyRecords[(uno+1):dos,"gvQueens_BritHY"]), sdHoneyYieldBrit=sd(colonyRecords[(uno+1):dos,"gvQueens_BritHY"]),
                       FitnessBrit=mean(colonyRecords[(uno+1):dos,"gvQueens_BritFit"]),sdFitnessBrit=sd(colonyRecords[(uno+1):dos,"gvQueens_BritFit"]),
                       HoneyYieldEu=mean(colonyRecords[(uno+1):dos,"gvQueens_EuHY"]), sdHoneyYieldEu=sd(colonyRecords[(uno+1):dos,"gvQueens_EuHY"]),
                       FitnessEu=mean(colonyRecords[(uno+1):dos,"gvQueens_EuFit"]), sdFitnessEu=sd(colonyRecords[(uno+1):dos,"gvQueens_EuFit"]),
                       pheHoneyYieldBrit=mean(colonyRecords[(uno+1):dos,"pheQueens_BritHY"]), sdpheHoneyYieldBrit=sd(colonyRecords[(uno+1):dos,"pheQueens_BritHY"]),
                       pheFitnessBrit=mean(colonyRecords[(uno+1):dos,"pheQueens_BritFit"]),sdpheFitnessBrit=sd(colonyRecords[(uno+1):dos,"pheQueens_BritFit"]),
                       pheHoneyYieldEu=mean(colonyRecords[(uno+1):dos,"pheQueens_EuHY"]), sdpheHoneyYieldEu=sd(colonyRecords[(uno+1):dos,"pheQueens_EuHY"]),
                       pheFitnessEu=mean(colonyRecords[(uno+1):dos,"pheQueens_EuFit"]), sdpheFitnessEu=sd(colonyRecords[(uno+1):dos,"pheQueens_EuFit"]),
                       Homocigosity=mean(colonyRecords[(uno+1):dos,"pHomBrood"]),sdHomocigosity=sd(colonyRecords[(uno+1):dos,"pHomBrood"]),
                       sdIBD=sd(colonyRecords[(uno+1):dos,"IBD"]), survivingCar=sum(colonyRecords[(uno+1):unoymedio,"IBD"]==0), Eucdist=1)
  
  
  
  #record values of Carnica population
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car", Rep=Rep)
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car", Rep=Rep)
  tres<-nrow(colonyRecords)
  
  #create a dataframe with the carnica mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
  newrow2<- data.frame(MeanIBD=mean(colonyRecords[(dos+1):tres,"IBD"]), VarIBD=var(colonyRecords[(dos+1):tres,"IBD"]), Year=year,Rep=Rep,Population="Car",
                       HoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]), sdHoneyYieldBrit=sd(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]),
                       FitnessBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritFit"]),sdFitnessBrit=sd(colonyRecords[(dos+1):tres,"gvQueens_BritFit"]),
                       HoneyYieldEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]), sdHoneyYieldEu=sd(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]),
                       FitnessEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuFit"]), sdFitnessEu=sd(colonyRecords[(dos+1):tres,"gvQueens_EuFit"]),
                       pheHoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"pheQueens_BritHY"]), sdpheHoneyYieldBrit=sd(colonyRecords[(dos+1):tres,"pheQueens_BritHY"]),
                       pheFitnessBrit=mean(colonyRecords[(dos+1):tres,"pheQueens_BritFit"]),sdpheFitnessBrit=sd(colonyRecords[(dos+1):tres,"pheQueens_BritFit"]),
                       pheHoneyYieldEu=mean(colonyRecords[(dos+1):tres,"pheQueens_EuHY"]), sdpheHoneyYieldEu=sd(colonyRecords[(dos+1):tres,"pheQueens_EuHY"]),
                       pheFitnessEu=mean(colonyRecords[(dos+1):tres,"pheQueens_EuFit"]), sdpheFitnessEu=sd(colonyRecords[(dos+1):tres,"pheQueens_EuFit"]),
                       Homocigosity=mean(colonyRecords[(dos+1):tres,"pHomBrood"]),sdHomocigosity=sd(colonyRecords[(dos+1):tres,"pHomBrood"]),
                       sdIBD=sd(colonyRecords[(dos+1):tres,"IBD"]), survivingCar=sd(colonyRecords[(dos+1):tres,"IBD"]),Eucdist=sd(colonyRecords[(dos+1):tres,"IBD"]))
  
  #Combine what we had in each dataframe with the new info, so each year the dataframe updates with new values
  MeanVarMel<-rbind(MeanVarMel,newrow1)
  MeanVarCar<-rbind(MeanVarCar,newrow2)
  save.image("Prueba.RData")
  print("prueba saved")
    
  
} #end of year loop
ahaha<-rbind(Eucdf2,Eucdf)
CombinedDf<-rbind(MeanVarCar,MeanVarMel)
save.image(paste0("Selection",Rep,"_import",param1Value,".RData"))

sink(file = paste0("Selection",Rep,"_import",param1Value,".RData")) # saving data from the burnin
cat(save.image(paste0("Selection",Rep,"_import",param1Value,".RData")))
sink()
write.csv(CombinedDf, paste0("scenarioData_rep",Rep,"_import",param1Value,".csv"))
sink(file = paste0("scenarioData_rep",Rep,"_import",param1Value,".csv")) # TODO: convert this to a CSV export
cat(write.csv(CombinedDf, paste0("scenarioData_rep",Rep,"_import",param1Value,".csv")))
sink()
write.csv(AlphapartCar, paste0("AlphapartCar",Rep,"_import",param1Value,".csv"))
sink(file = paste0("AlphapartCar",Rep,"_import",param1Value,".csv"))
cat(write.csv(AlphapartCar, paste0("AlphapartCar",Rep,"_import",param1Value,".csv")))
sink()
write.csv(Alphapart, paste0("Alphapart",Rep,"_import",param1Value,".csv"))
sink(file = paste0("Alphapart",Rep,"_import",param1Value,".csv"))
cat(write.csv(Alphapart, paste0("Alphapart",Rep,"_import",param1Value,".csv")))
sink()
write.csv(AlphapartCombined, paste0("AlphapartCombined",Rep,"_import",param1Value,".csv"))
sink(file = paste0("AlphapartCombined",Rep,"_import",param1Value,".csv"))
cat(write.csv(AlphapartCombined, paste0("AlphapartCombined",Rep,"_import",param1Value,".csv")))
sink()
write.csv(ahaha, paste0("Eucdist_rep",Rep,"_import",param1Value,".csv"))
sink(file = paste0("Eucdist_rep",Rep,"_import",param1Value,".csv")) # TODO: convert this to a CSV export
cat(write.csv(ahaha, paste0("Eucdist_rep",Rep,"_import",param1Value,".csv")))
sink()
remove(MeanVarCar)
remove(MeanVarMel)
setwd("..") # scenario folder
setwd("..") # rep folder
}
