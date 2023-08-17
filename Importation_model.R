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

cat("We are doing:\n")
cat("  * replicate: ", Rep, "\n")
cat("  *", burninOrScenario, "\n")
cat("  * parmeter1: ", param1Value, "\n")
cat("  * parmeter2: ", param2Value, "\n")

if (burninOrScenario == "burnin") {
  
  dir <- paste0("Rep_", Rep) #name of the directory
  unlink(dir, recursive = TRUE) # delete the folder if it exists 
  dir.create(path = dir) # create the folder
  setwd(dir) # go into the folder
  
  dir <- "burnin" 
  unlink(dir) #delete the burnin folder if there is a previous one
  dir.create(path = dir)
  setwd(dir) #create a burnin folder where we will store the data for the burnin
  
  # code for burnin goes here
  maintainCarSize <- function(age0 = NULL, age1 = NULL) {
    if ((nColonies(age0) + nColonies(age1)) > CarSize) { # check if the sum of all colonies is greater than apiary size
      IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
      splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
      age0split <- splits0$pulled # create an object for age 0 splits
      age0swarm <- splits0$remnant # create an object for swarms and superseded colonies
      age0needed <- CarSize - nColonies(age1) # calculate the number of age 0 colonies that are needed to fill up the apiary
      splitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
      if (age0needed <= nColonies(age0swarm)) { # check if the number of age 0 colonies needed is lower or equal to age 0 swarms
        swarmID <- sample(getId(age0swarm), age0needed) # if yes, select the ids of swarms that will stay in apiary
        swarmTMP <- pullColonies(age0swarm, ID = swarmID) # pull out those selected age0 swarms
        age0 <- swarmTMP$pulled # put selected swarms to age 0 object
      } else if (age0needed > nColonies(age0swarm)) { # in case when age 0 needed is grater than number of swarm select splits
        nSplitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
        gvCarQueens_EuHY <- sapply(getPheno(age0split, caste = "queen"), function(x) x[1,2]) #get pheno for production of age 0 carnica queens
        queensID<-names(sort(gvCarQueens_EuHY,decreasing=T))
        age0CarqueensID<-queensID[1:nSplitsNeeded]
        splitTmp <- pullColonies(age0split, ID = age0CarqueensID) # pull the splits
        splits <- splitTmp$pulled # select pulled splits
        age0 <- c(age0swarm, splits) # combine splits and swarms in age 0 object
      }
      return(age0)
    }
  }
  
  maintainIrelandSize <- function(age0 = NULL, age1 = NULL) {
    if ((nColonies(age0) + nColonies(age1)) > CarSize) { # check if the sum of all colonies is greater than apiary size
      IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
      splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
      age0split <- splits0$pulled # create an object for age 0 splits
      age0swarm <- splits0$remnant # create an object for swarms and superseded colonies
      age0needed <- CarSize - nColonies(age1) # calculate the number of age 0 colonies that are needed to fill up the apiary
      splitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
      if (age0needed <= nColonies(age0swarm)) { # check if the number of age 0 colonies needed is lower or equal to age 0 swarms
        swarmID <- sample(getId(age0swarm), age0needed) # if yes, select the ids of swarms that will stay in apiary
        swarmTMP <- pullColonies(age0swarm, ID = swarmID) # pull out those selected age0 swarms
        age0 <- swarmTMP$pulled # put selected swarms to age 0 object
      } else if (age0needed > nColonies(age0swarm)) { # in case when age 0 needed is grater than number of swarm select splits
        nSplitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
        splitId <- sample(getId(age0split), nSplitsNeeded) # select ids of splits
        splitTmp <- pullColonies(age0split, ID = splitId) # pull the splits
        splits <- splitTmp$pulled # select pulled splits
        age0 <- c(age0swarm, splits) # combine splits and swarms in age 0 object
      }
      return(age0)
    }
  }
  
  # Load packages
  library(AlphaSimR)
  library(ggplot2)
  #library(tictoc)
  library(R6)
  library(nadiv)
  library(Matrix)
  library(SIMplyBee)
  library(dplyr)
  library(tidyr)

  # TODO: replace with devtools installation from Github once the package is operational
  # Source the development version of AlphaSimR
  
  
  # Founder population parameters -------------------------------------------------------------------
  nMelN = 200                  # Number of Mellifera
  nCar =200                   # Number of Carnica
  nMelnI = 200
  nChr = 1                     # Number of chromomsome
  nDronesPerQueen = 100
  nSegSites = 100              # Number of segregating sites
  
  # Population parameters -------------------------------------------------------------------
  # Number of repeats
  nYear <- 10                    # Number of years
  IrelandSize<-100              #Ireland population size
  IrelandnISize<-100
  CarSize<-100                  #Carnica pop size
  #LigSize<-100                  #Ligustica pop size
  nWorkers <- 10                # Number of workers in a full colony
  nDrones <- 100                 # Number of drones in a full colony (typically nWorkers * 0.2 (not in the example))
  pFathers <- nFathersPoisson   # Number of drones the queen mates with (could also be a function)
  nVirginQueens <- 1            # Number of created virgin queens
  
  # Period parameters -------------------------------------------------------------------
  # Period 1 (spring)
  p1swarm <- 0.05              # Percentage of colonies that swarm in period 1
  p1supersede <- 0.05          # Percentage of colonies that supersede in period 1
  p1collapse <- 0.10           # Percentage of colonies that  collapse in period 1
  
  # Period2 (summer)
  p2swarm <- 0.01              # Percentage of colonies that swarm in period 2
  p2supersede <- p1supersede   # Percentage of colonies that supersede in period 2
  p2collapse <- p1collapse     # Percentage of colonies that collapse in period 2
  
  # Period3 (winter)
  p3collapseAge0 <- 0.25      # Percentage of age 0 colonies that collapse in period 3
  p3collapseAge1 <- 0.3       # Percentage of age 2 colonies that collapse in period 3
  
  #Import parameters -------------------------------------------------------------------
  pImport <- param1Value              # Percentage import from carnica to mellifera
  
  # Create data frames for recording the number of age0 and age1 colonies, csd variability and for recording cpu time
  loopTime <- data.frame(Rep = NA, tic = NA, toc = NA, msg = NA, time = NA)
  
  # Prepare recording function
  data_rec <- function(datafile, colonies, year, population, Rep) {
    queens = mergePops(getQueen(colonies))
    IBDh = apply(getIbdHaplo(queens), MARGIN = 1, FUN =  function(X) sum(X %in% 1:((nMelN+nMelnI)*2)/length(X)))
    datafile = rbind(datafile,
                     data.frame(colonies             = deparse(substitute(colonies)),
                                population           = population,
                                year                 = year,
                                Rep                  = Rep,
                                Id                   = queens@id,
                                MId                  = queens@mother,
                                FId                  = queens@father,
                                nFathers             = nFathers(queens),
                                nDPQ                 = sapply(getFathers(queens), function(x) length(unique(x@mother))),
                                nCsdAlColony         = sapply(colonies@colonies, function(x) nCsdAlleles(x, collapse = TRUE)),
                                nCsdApiary           = rep(nCsdAlleles(colonies, collapse = TRUE), queens@nInd),
                                pHomBrood            = calcQueensPHomBrood(queens),
                                gvQueens_BritHY  = sapply(getGv(colonies, caste = "queen"), function(x) x[1,1]),
                                gvQueens_EuHY    = sapply(getGv(colonies, caste = "queen"), function(x) x[1,2]),
                                gvQueens_BritFit  = sapply(getGv(colonies, caste = "queen"), function(x) x[1,3]),
                                gvQueens_EuFit    = sapply(getGv(colonies, caste = "queen"), function(x) x[1,4]),
                                pheQueens_BritHY  = sapply(getPheno(colonies, caste = "queen"), function(x) x[1,1]),
                                pheQueens_EuHY    = sapply(getPheno(colonies, caste = "queen"), function(x) x[1,2]),
                                pheQueens_BritFit  = sapply(getPheno(colonies, caste = "queen"), function(x) x[1,3]),
                                pheQueens_EuFit    = sapply(getPheno(colonies, caste = "queen"), function(x) x[1,4]),
                                IBD = sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
                                
                                
                     ))}
  colonyRecords = NULL
  
  
  pImport <- param1Value #parameter 1
  
  load(file="/exports/cmvm/eddie/eb/groups/HighlanderLab/visitors/icarlos_honeybee_introgression/FounderGenomes_ThreePop_16chr.RData")
  
  
  # quick haplo to get the founder genomes for now

  
  # STEP 2: Create SP object and write in the global simulation/population parameters
  SP <- SimParamBee$new(founderGenomes, csdChr = ifelse(nChr >= 3, 3, 1), nCsdAlleles = 128)
  SP$nWorkers <- nWorkers
  SP$nDrones <- nDrones
  SP$nFathers <- pFathers
  SP$nVirginQueens <- nVirginQueens
  SP$swarmP <- 0.5                # Probability of swarming? (TODO: double check this)
  SP$splitP <- 0.3
  SP$setTrackPed(TRUE)            # Track the pedigree
  SP$setTrackRec(TRUE)            # Track the recombination
  
  SP$addSnpChip(nSnpPerChr = 10)   # Add a SNP chip with 3 SNPs per chromosome
  csdChr <- SP$csdChr             # define csd chromomsome
  
  # Add traits - taken from the QuantGen vignette 
  mean <- c(0,0,0,0)
  varA <- c(0.25,0.25,0.1,0.1)
  corA <- matrix(data = c(  1.0, 0.9,  0.0, 0.0, 
                            0.9, 1.0,  0.0, 0.0, 
                            0.0, 0.0,  1.0, -0.5, 
                            0.0, 0.0, -0.5, 1.0), nrow = 4, byrow = TRUE)
  
  SP$addTraitA(nQtlPerChr = 100, mean = mean, var = varA, corA = corA,
               name = c("QueenHYBrit", "QueenHYEu", "FitnessBrit", "FitnessEu"))
  
  varE <- c(0.75,0.75,0.9,0.9)
  
  # TODO: what is a reasonable environmental correlation between queen and worker effects?
  corE <- matrix(data =     c(  1.0, 0.0,  0.0, 0.0, 
                                0.0, 1.0,  0.0, 0.0, 
                                0.0, 0.0,  1.0, 0.0, 
                                0.0, 0.0,  0.0, 1.0 ), nrow = 4, byrow = TRUE)
  
  SP$setVarE(varE = varE, corE = corE)
  
  # STEP 3: Set up your base population
  # Create a base population for A. m. mellifera, A. m. mellifera cross, and A. m. carnica (400 of each)
  #virginQueens <- list(Mel = createVirginQueens(x = founderGenomes[1:(nMelN)]),
  #MelnI = createVirginQueens(x = founderGenomes[(nMelN+1):(nMelnI+nMelN)]),
  #Car = createVirginQueens(x = founderGenomes[(nMelN+nMelnI+1):(nMelN + nMelnI + nCar)]))
  
  basePop <- newPop(founderGenomes)
  basePopMel <- basePop[1:100]
  virginQueensMel <- randCross(basePopMel, nProgeny = 1, nCrosses = nMelN)
  virginQueensMel <- SIMplyBee:::editCsdLocus(virginQueensMel, simParamBee = SP)
  virginQueensMel@sex[] <- "F"
  SP$changeCaste(id = virginQueensMel@id, caste = "V")
  
  basePopMelnI <- basePop[101:200]
  virginQueensMelnI <- randCross(basePopMelnI, nProgeny = 1, nCrosses = nMelnI)
  virginQueensMelnI <- SIMplyBee:::editCsdLocus(virginQueensMelnI, simParamBee = SP)
  virginQueensMelnI@sex[] <- "F"
  SP$changeCaste(id = virginQueensMelnI@id, caste = "V")
  
  basePopCar <- basePop[801:900]
  virginQueensCar <- randCross(basePopCar, nProgeny = 1, nCrosses = nCar)
  virginQueensCar <- SIMplyBee:::editCsdLocus(virginQueensCar, simParamBee = SP)
  virginQueensCar@sex[] <- "F"
  SP$changeCaste(id = virginQueensCar@id, caste = "V")
  
  virginQueens <- list(Mel =  virginQueensMel, MelnI =  virginQueensMelnI, Car =  virginQueensCar)
  # Create drones for A. m. mellifera, A. m. mellifera cross, and A. m. carnica
  drones <- list(Mel = createDrones(x = virginQueens$Mel[(IrelandSize+1):(nMelN)], nInd = nDronesPerQueen),
                 MelnI = createDrones(x = virginQueens$MelnI[(IrelandnISize+1):(nMelnI)], nInd = nDronesPerQueen),
                 Car = createDrones(x = virginQueens$Car[(CarSize+1):nCar], nInd = nDronesPerQueen))

  # Get fathers for Mel, MelCross and Car
  fathersMel <- pullDroneGroupsFromDCA(drones$Mel, n = nInd(virginQueens$Mel[1:IrelandSize]), nDrones = nFathersPoisson)
  fathersMelnI <- pullDroneGroupsFromDCA(drones$MelnI, n = nInd(virginQueens$MelnI[1:IrelandnISize]), nDrones = nFathersPoisson)
  fathersCar <- pullDroneGroupsFromDCA(drones$Car, n = nInd(virginQueens$Car[1:CarSize]), nDrones = nFathersPoisson)
  
  # Mate virgin queens with fathers to make them queens
  queens <- list(Mel = SIMplyBee::cross(x = virginQueens$Mel[1:IrelandSize], drones = fathersMel),
                 MelnI = SIMplyBee::cross(x = virginQueens$MelnI[1:IrelandnISize], drones = fathersMelnI),
                 Car = SIMplyBee::cross(x = virginQueens$Car[1:CarSize], drones = fathersCar))
  
  
  year=1
  nYear=10
  
  # Start the year-loop ------------------------------------------------------------------
  for (year in 1:nYear) {
    print("Starting the cycle")
    #year <- 1 (Use this to check that things are working without setting the whole for loop off )
    #year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    
    # If this is the first year, create some colonies to start with
    if (year == 1) {
      print("Creating initial colonies")
      age1 <- list(Mel = createMultiColony(x = queens$Mel, n = IrelandSize),
                   MelnI = createMultiColony(x = queens$MelnI, n = IrelandnISize),
                   Car = createMultiColony(x = queens$Car, n = CarSize))
      
      # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
    } else {
      age2 <- list(Mel = age1$Mel, MelnI = age1$MelnI, Car = age1$Car) #, Lig = age1$Lig
      age1 <- list(Mel = age0$Mel, MelnI = age0$MelnI, Car = age0$Car) #, Lig = age0$Lig
      age0 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #, Lig = NULL
      age0p1 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #, Lig = NULL
      age0p2 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #, Lig = NULL
    }
    
    # Period1 ------------------------------------------------------------------
    # Build-up the colonies
    print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
    print(Sys.time())
    age1 <- list(Mel = buildUp(age1$Mel),
                 MelnI = buildUp(age1$MelnI),
                 Car = buildUp(age1$Car))
    
    if (year > 1) {
      age2 <- list(Mel = buildUp(age2$Mel),
                   MelnI = buildUp(age2$MelnI),
                   Car = buildUp(age2$Car))
      
    }
    
    # Split all age1 colonies
    print("Splitting the colonies")
    print(Sys.time())
    tmp <- list(Mel = split(age1$Mel),
                MelnI = split(age1$MelnI),
                Car = split(age1$Car))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    
    # The queens of the splits are 0 years old
    age0p1 <- list(Mel = tmp$Mel$split, MelnI = tmp$MelnI$split, Car = tmp$Car$split) #,Lig = tmp$Lig$split
    
    if (year > 1) {
      # Split all age2 colonies
      tmp <- list(Mel = split(age2$Mel),
                  MelnI = split(age2$MelnI),
                  Car = split(age2$Car))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant
      ) 
      # The queens of the splits are 0 years old
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                     MelnI = c(age0p1$MelnI, tmp$MelnI$split),
                     Car = c(age0p1$Car, tmp$Car$split))
      
    }
    
    
    # Create virgin queens
    # Sample colony for the virgin queens
    print("Create virgin queens, period 1")
    print(Sys.time())
   # virginDonor <- list(Mel = sample.int(n = nColonies(age1$Mel), size = 1),
   #                    MelnI = sample.int(n = nColonies(age1$MelnI), size = 1),
   #                   Car = sample.int(n = nColonies(age1$Car), size = 1))
    
    # Virgin queens for splits!
    
    
    pImport<-0
    
    virginQueens <- list(Mel = createVirginQueens(age1$Mel, collapse=T, nInd=10),
                         MelnI = createVirginQueens(age1$MelnI, collapse=T, nInd=10),
                         Car = createVirginQueens(age1$Car, collapse=T, nInd=10))
    
    virginQueens<-list(Mel = mergePops(virginQueens$Mel),
                       MelnI = mergePops(virginQueens$MelnI),
                       Car = mergePops(virginQueens$Car))
               
    #requeen with carnica and mellifera for Mel and carnica for Car
    nColoniesMel<-nColonies(age0p1$Mel) 
    nColoniesCar<-nColonies(age0p1$Car)
    nColoniesMelnI<-nColonies(age0p1$MelnI)
    age0p1 <- list(Mel = reQueen(age0p1$Mel, queen = virginQueens$Mel[sample(1:nInd(virginQueens$Mel), size = nColoniesMel, replace = FALSE)]),
                   MelnI = reQueen(age0p1$MelnI, queen = virginQueens$MelnI[sample(1:nInd(virginQueens$MelnI), size = nColoniesMelnI, replace = FALSE)]),
                   Car = reQueen(age0p1$Car, queen = virginQueens$Car[sample(1:nInd(virginQueens$Car), size = nColoniesCar, replace = FALSE)]))
    
    
    
    
    # Swarm a percentage of age1 colonies
    print("Swarm colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
                MelnI = pullColonies(age1$MelnI, p = p1swarm),
                Car = pullColonies(age1$Car, p = p1swarm))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                MelnI = swarm(tmp$MelnI$pulled),
                Car = swarm(tmp$Car$pulled))
    
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                   MelnI = c(age0p1$MelnI, tmp$MelnI$remnant),
                   Car = c(age0p1$Car, tmp$Car$remnant))
    
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 MelnI = c(age1$MelnI, tmp$MelnI$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))
    
    
    
    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
                  MelnI = pullColonies(age2$MelnI, p = p1swarm),
                  Car = pullColonies(age2$Car, p = p1swarm))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
                  MelnI = swarm(tmp$MelnI$pulled),
                  Car = swarm(tmp$Car$pulled))
      
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                     MelnI = c(age0p1$MelnI, tmp$MelnI$remnant),
                     Car = c(age0p1$Car, tmp$Car$remnant))
      
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                   MelnI = c(age2$MelnI, tmp$MelnI$swarm),
                   Car = c(age2$Car, tmp$Car$swarm))
      
    }
    
    # Supersede age1 colonies
    print("Supersede colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
                MelnI = pullColonies(age1$MelnI, p = p1supersede),
                Car = pullColonies(age1$Car, p = p1supersede))
    #,Lig = pullColonies(age1$Lig, p = p1supersede))
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    #,Lig = tmp$Lig$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                MelnI = supersede(tmp$MelnI$pulled),
                Car = supersede(tmp$Car$pulled))
    #,Lig = supersede(tmp$Lig$pulled))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                   MelnI = c(age0p1$MelnI, tmp$MelnI),
                   Car = c(age0p1$Car, tmp$Car))
    
    
    if (year > 1) {
      # Supersede age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
                  MelnI = pullColonies(age2$MelnI, p = p1supersede),
                  Car = pullColonies(age2$Car, p = p1supersede))
      #,Lig = pullColonies(age2$Lig, p = p1supersede))
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      #,Lig = tmp$Lig$remnant)
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  MelnI = supersede(tmp$MelnI$pulled),
                  Car = supersede(tmp$Car$pulled))
      #,Lig = supersede(tmp$Lig$pulled))
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                     MelnI = c(age0p1$MelnI, tmp$MelnI),
                     Car = c(age0p1$Car, tmp$Car))
      #,Lig = c(age0p1$Lig, tmp$Lig))
    }
    
    # Mate the split colonies
    print("Mate split colonies, P1")
    print(Sys.time())
    if (year == 1) {
      DCAMel <- createDCA(age1$Mel)
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
      DCAMelnI <- createDCA(age1$MelnI)
      age0p1$MelnI <- cross(age0p1$MelnI, drones = pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p1$MelnI), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
      DCAMelnI <- createDCA(c(age1$MelnI, age2$MelnI)) 
      age0p1$MelnI <- cross(age0p1$MelnI, drones = pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p1$MelnI), nDrones = nFathersPoisson))
    }
    
    # Collapse
    print("Collapse colonies, P1")
    print(Sys.time())
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
                 MelnI = selectColonies(age1$MelnI, p = 1 - p1collapse),
                 Car = selectColonies(age1$Car, p = 1 - p1collapse))
    #,Lig = selectColonies(age1$Lig, p = 1 - p1collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
                   MelnI = selectColonies(age2$MelnI, p = 1 - p1collapse),
                   Car = selectColonies(age2$Car, p = 1 - p1collapse))
      #,Lig = selectColonies(age2$Lig, p = 1 - p1collapse))
    }
    
    # Period2 ------------------------------------------------------------------
    print("PERIOD 2")
    # Swarm a percentage of age1 colonies
    # Mellifera
    print("Swarm colonies, P2")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p2swarm),
                MelnI = pullColonies(age1$MelnI, p = p2swarm),
                Car = pullColonies(age1$Car, p = p2swarm))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                MelnI = swarm(tmp$MelnI$pulled),
                Car = swarm(tmp$Car$pulled))
    
    # The queens of the remnant colonies are of age 0
    age0p2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
    
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 MelnI = c(age1$MelnI, tmp$MelnI$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))
    
    
    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
                  MelnI = pullColonies(age2$MelnI, p = p2swarm),
                  Car = pullColonies(age2$Car, p = p2swarm))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
                  MelnI = swarm(tmp$MelnI$pulled),
                  Car = swarm(tmp$Car$pulled))
      
      # The queens of the remnant colonies are of age 0
      age0p2 <- list(Mel = tmp$Mel$remnant,
                     MelnI = tmp$MelnI$remnant,
                     Car = tmp$Car$remnant)
      
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                   MelnI = c(age2$MelnI, tmp$MelnI$swarm),
                   Car = c(age2$Car, tmp$Car$swarm))
      
    }
    
    # Supersede a part of age1 colonies
    print("Supersede colonies, P2")
    print(Sys.time())
    
    tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
                MelnI = pullColonies(age1$MelnI, p = p2supersede),
                Car = pullColonies(age1$Car, p = p2supersede))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                MelnI = supersede(tmp$MelnI$pulled),
                Car = supersede(tmp$Car$pulled))
    
    # The queens of superseded colonies are of age 0
    age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                   MelnI = c(age0p2$MelnI, tmp$MelnI),
                   Car = c(age0p2$Car, tmp$Car))
    
    
    if (year > 1) {
      # Supersede a part of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
                  MelnI = pullColonies(age2$MelnI, p = p2supersede),
                  Car = pullColonies(age2$Car, p = p2supersede))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  MelnI = supersede(tmp$MelnI$pulled),
                  Car = supersede(tmp$Car$pulled))
      
      # The queens of superseded colonies are of age 0
      age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                     MelnI = c(age0p2$MelnI, tmp$MelnI),
                     Car = c(age0p2$Car, tmp$Car))
      
    }
    
    # Replace all the drones
    print("Replace Drones, P2")
    print(Sys.time())
    
    age1$Mel <- replaceDrones(age1$Mel)
    age1$MelnI <- replaceDrones(age1$MelnI)
    age1$Car <- replaceDrones(age1$Car)
    
    if (year > 1) {
      age2$Mel <- replaceDrones(age2$Mel)
      age2$MelnI <- replaceDrones(age2$MelnI)
      age2$Car <- replaceDrones(age2$Car)
      
    }
    
    # Mate the colonies
    # Import p percentage of carnica colonies into mellifera DCA
    print("Mate colonies, P2")
    print(Sys.time())
    
    
    if (year == 1) {
      DCAMel <- createDCA(age1$Mel)
      age0p2$Mel <- cross(age0p2$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
      DCAMelnI <- createDCA(age1$MelnI)
      age0p2$MelnI <- cross(age0p2$MelnI, drones = pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p2$MelnI), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      fathersMel <- pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson)
      age0p2$Mel <- cross(age0p2$Mel, drones = fathersMel)
      DCACar <- createDCA(c(age1$Car, age2$Car))
      fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
      age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
      DCAMelnI <- createDCA(c(age1$MelnI, age2$MelnI))
      fathersMelnI <-  pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p2$MelnI), nDrones = nFathersPoisson)
      age0p2$MelnI <- cross(age0p2$MelnI, drones = fathersMelnI)
    }
    
    # Collapse
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p2collapse),
                 MelnI = selectColonies(age1$MelnI, p = 1 - p2collapse),
                 Car = selectColonies(age1$Car, p = 1 - p2collapse))
    
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
                   MelnI = selectColonies(age2$MelnI, p = 1 - p2collapse),
                   Car = selectColonies(age2$Car, p = 1 - p2collapse))
      
    }
    
    # Merge all age 0 colonies (from both periods)
    age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
                 MelnI = c(age0p1$MelnI, age0p2$MelnI),
                 Car = c(age0p1$Car, age0p2$Car))
    
    
    
    # Period3 ------------------------------------------------------------------
    # Collapse age0 queens
    print("PERIOD 3")
    print("Collapse colonies, P3")
    print(Sys.time())
    
    #selection on Fitness queen gv
    
    #age0
    #Mellifera
    
    MelBritFit<-sapply(getPheno(age0$Mel, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID<-names(sort(MelBritFit,T))
    Nselectcolon<-round(length(queensID)*(1-p3collapseAge0)) #calculate how many colonies will collapse
    age0MelqueensID<-queensID[1:Nselectcolon] #select the queens ids that will not collapse
    
    MelBritFit1<-sapply(getPheno(age1$Mel, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID1<-names(sort(MelBritFit1,T))
    Nselectcolon1<-round(length(queensID1)*(1-p3collapseAge1)) #calculate how many colonies will collapse
    age1MelqueensID<-queensID1[1:Nselectcolon1] #select the queens ids that will not collapse
    
    MelnIBritFit<-sapply(getPheno(age0$MelnI, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID2<-names(sort(MelnIBritFit,T))
    Nselectcolon2<-round(length(queensID2)*(1-p3collapseAge0)) #calculate how many colonies will collapse
    age0MelnIqueensID<-queensID2[1:Nselectcolon2] #select the queens ids that will not collapse
    
    MelnIBritFit1<-sapply(getPheno(age1$MelnI, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID3<-names(sort(MelnIBritFit1,T))
    Nselectcolon3<-round(length(queensID3)*(1-p3collapseAge1)) #calculate how many colonies will collapse
    age1MelnIqueensID<-queensID3[1:Nselectcolon3] #select the queens ids that will not collapse
    
    CarEuFit<-sapply(getPheno(age0$Car, caste = "queen"), function(x) x[1,4])#select on Mel Brit fitness
    queensID<-names(sort(CarEuFit,T))
    Nselectcolon<-round(length(queensID)*(1-p3collapseAge0)) #calculate how many colonies will collapse
    age0CarqueensID<-queensID[1:Nselectcolon] #select the queens ids that will not collapse
    
    CarEuFit1<-sapply(getPheno(age1$Car, caste = "queen"), function(x) x[1,4])#select on Mel Brit fitness
    queensID1<-names(sort(CarEuFit1,T))
    Nselectcolon1<-round(length(queensID1)*(1-p3collapseAge1)) #calculate how many colonies will collapse
    age1CarqueensID<-queensID1[1:Nselectcolon1] #select the queens ids that will not collapse
    
    #Collapse 
    age0 <- list(Mel = selectColonies(age0$Mel, ID = age0MelqueensID),
                 MelnI = selectColonies(age0$MelnI, ID = age0MelnIqueensID),
                 Car = selectColonies(age0$Car, ID = age0CarqueensID))
    
    age1 <- list(Mel = selectColonies(age1$Mel, ID = age1MelqueensID),
                 MelnI = selectColonies(age1$MelnI, ID = age1MelnIqueensID),
                 Car = selectColonies(age1$Car, ID= age1CarqueensID))
    
    age2 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #We don't need this but just to show the workflow!!!
    
    # Maintain the number of colonies ------------------------------------------
    # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
    print("Maintain the number, P2")
    print(Sys.time())
    
    age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
    age0$MelnI <- maintainIrelandSize(age0 = age0$MelnI, age1 = age1$MelnI)
    age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)
    
    
    for (subspecies in c("Mel", "MelnI", "Car")) {     
      if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == IrelandSize
          | (nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == CarSize) 
      {
      } 
      else 
      {stop(paste0("The number of colonies for ", subspecies, " does not match the population size!"))}
    }
    
    
    #track mean IBD and variance of IBD
    
    #this creates two empty dataframes where all the information of each year and each replica will be recorded
    
    columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","sdHoneyYieldBrit","HoneyYieldEu","sdHoneyYieldEu",
                     "FitnessBrit","sdFitnessBrit","FitnessEu","sdFitnessEu","pheHoneyYieldBrit","sdpheHoneyYieldBrit","pheHoneyYieldEu","sdpheHoneyYieldEu",
                     "pheFitnessBrit","sdpheFitnessBrit","pheFitnessEu","sdpheFitnessEu","Homocigosity","sdHomocigosity","sdIBD")
    
    if (year==1){ 
      columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","sdHoneyYieldBrit","HoneyYieldEu","sdHoneyYieldEu",
                       "FitnessBrit","sdFitnessBrit","FitnessEu","sdFitnessEu","pheHoneyYieldBrit","sdpheHoneyYieldBrit","pheHoneyYieldEu","sdpheHoneyYieldEu",
                       "pheFitnessBrit","sdpheFitnessBrit","pheFitnessEu","sdpheFitnessEu","Homocigosity","sdHomocigosity","sdIBD")
      MeanVarMel <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for mellifera
      colnames(MeanVarMel)<-columnheaders
      MeanVarCar <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for carnica
      colnames(MeanVarCar)<-columnheaders
      MeanVarMelnI <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for carnica
      colnames(MeanVarMelnI)<-columnheaders
    } 
    
    if (year==1){
      uno<-0
    }else{
      uno<-nrow(colonyRecords) #this is to show where to start recording values on the dataframe
    }
    
    
    #record values of Mellifera population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel", Rep=Rep)
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
                         sdIBD=sd(colonyRecords[(uno+1):dos,"IBD"]))
    
    #record values for Mellifera without imports population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$MelnI, year = year, population = "MelnI", Rep=Rep)
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$MelnI, year = year, population = "MelnI", Rep=Rep)
    tres<-nrow(colonyRecords)#this is to show where to end recording values on the dataframe
    
    
    #create a dataframe with the mellifera mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
    newrow2<- data.frame(MeanIBD=mean(colonyRecords[(dos+1):tres,"IBD"]), VarIBD=var(colonyRecords[(dos+1):tres,"IBD"]), Year=year,Rep=Rep,Population="MelnI",
                         HoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]), sdHoneyYieldBrit=sd(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]),
                         FitnessBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritFit"]),sdFitnessBrit=sd(colonyRecords[(dos+1):tres,"gvQueens_BritFit"]),
                         HoneyYieldEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]), sdHoneyYieldEu=sd(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]),
                         FitnessEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuFit"]), sdFitnessEu=sd(colonyRecords[(dos+1):tres,"gvQueens_EuFit"]),
                         pheHoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"pheQueens_BritHY"]), sdpheHoneyYieldBrit=sd(colonyRecords[(dos+1):tres,"pheQueens_BritHY"]),
                         pheFitnessBrit=mean(colonyRecords[(dos+1):tres,"pheQueens_BritFit"]),sdpheFitnessBrit=sd(colonyRecords[(dos+1):tres,"pheQueens_BritFit"]),
                         pheHoneyYieldEu=mean(colonyRecords[(dos+1):tres,"pheQueens_EuHY"]), sdpheHoneyYieldEu=sd(colonyRecords[(dos+1):tres,"pheQueens_EuHY"]),
                         pheFitnessEu=mean(colonyRecords[(dos+1):tres,"pheQueens_EuFit"]), sdpheFitnessEu=sd(colonyRecords[(dos+1):tres,"pheQueens_EuFit"]),
                         Homocigosity=mean(colonyRecords[(dos+1):tres,"pHomBrood"]),sdHomocigosity=sd(colonyRecords[(dos+1):tres,"pHomBrood"]),
                         sdIBD=sd(colonyRecords[(dos+1):tres,"IBD"]))
    
    #record values of Carnica population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car", Rep=Rep)
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car", Rep=Rep)
    cua<-nrow(colonyRecords)
    
    #create a dataframe with the carnica mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
    newrow3<- data.frame(MeanIBD=mean(colonyRecords[(tres+1):cua,"IBD"]), VarIBD=var(colonyRecords[(tres+1):cua,"IBD"]), Year=year,Rep=Rep,Population="Car",
                         HoneyYieldBrit=mean(colonyRecords[(tres+1):cua,"gvQueens_BritHY"]), sdHoneyYieldBrit=sd(colonyRecords[(tres+1):cua,"gvQueens_BritHY"]),
                         FitnessBrit=mean(colonyRecords[(tres+1):cua,"gvQueens_BritFit"]),sdFitnessBrit=sd(colonyRecords[(tres+1):cua,"gvQueens_BritFit"]),
                         HoneyYieldEu=mean(colonyRecords[(tres+1):cua,"gvQueens_EuHY"]), sdHoneyYieldEu=sd(colonyRecords[(tres+1):cua,"gvQueens_EuHY"]),
                         FitnessEu=mean(colonyRecords[(tres+1):cua,"gvQueens_EuFit"]), sdFitnessEu=sd(colonyRecords[(tres+1):cua,"gvQueens_EuFit"]),
                         pheHoneyYieldBrit=mean(colonyRecords[(tres+1):cua,"pheQueens_BritHY"]), sdpheHoneyYieldBrit=sd(colonyRecords[(tres+1):cua,"pheQueens_BritHY"]),
                         pheFitnessBrit=mean(colonyRecords[(tres+1):cua,"pheQueens_BritFit"]),sdpheFitnessBrit=sd(colonyRecords[(tres+1):cua,"pheQueens_BritFit"]),
                         pheHoneyYieldEu=mean(colonyRecords[(tres+1):cua,"pheQueens_EuHY"]), sdpheHoneyYieldEu=sd(colonyRecords[(tres+1):cua,"pheQueens_EuHY"]),
                         pheFitnessEu=mean(colonyRecords[(tres+1):cua,"pheQueens_EuFit"]), sdpheFitnessEu=sd(colonyRecords[(tres+1):cua,"pheQueens_EuFit"]),
                         Homocigosity=mean(colonyRecords[(tres+1):cua,"pHomBrood"]),sdHomocigosity=sd(colonyRecords[(tres+1):cua,"pHomBrood"]),
                         sdIBD=sd(colonyRecords[(tres+1):cua,"IBD"]))
    
    #Combine what we had in each dataframe with the new info, so each year the dataframe updates with new values
    MeanVarMel<-rbind(MeanVarMel,newrow1)
    MeanVarMelnI<-rbind(MeanVarMelnI,newrow2)
    MeanVarCar<-rbind(MeanVarCar,newrow3)
    
  } # Year-loop
  
  print("Saving image data")
  save.image(paste0("Burnin",Rep,".RData"))
  
  
  sink(file = paste0("Burnin",Rep,".RData")) # saving data from the burnin
  cat(save.image(paste0("Burnin",Rep,".RData")))
  sink()
  
  
  setwd("..") # burnin folder
  setwd("..") # rep folder
}

if (burninOrScenario == "scenario") {
  dir <- paste0("Rep_", Rep)
  if (!dir.exists(dir)) {
    stop(paste0("Rep folder is missing for Rep ", Rep, "!"))
  }
  setwd(dir)
  
  if (!dir.exists("burnin")) {
    stop(paste0("Burnin folder is missing for Rep ", Rep, "!"))
  }
  
  # code for loading burnin goes here
  dire <- paste0("scenario_param1_", param1Value, "_param2_", param2Value)
  unlink(dir)
  dir.create(path = dire)
  tula<- args[3]
  setwd("burnin")
  load(paste0("Burnin",Rep,".RData"))
  setwd("..")
  setwd(dire)
  param1Value<-as.numeric(tula)
  pImport<-param1Value
  
  # code for scenario goes here
  # Load packages
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
  year=11
  nYear=20
  for (year in 11:nYear) {
    print("Starting the cycle")
    #year <- 1 (Use this to check that things are working without setting the whole for loop off )
    #year <- year + 1
    cat(paste0("Year: ", year, "/", nYear, "\n"))
    
    # If this is the first year, create some colonies to start with
    if (year == 1) {
      print("Creating initial colonies")
      age1 <- list(Mel = createMultiColony(x = queens$Mel, n = IrelandSize),
                   MelnI = createMultiColony(x = queens$MelnI, n = IrelandSize),
                   Car = createMultiColony(x = queens$Car, n = CarSize))
      
      # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
    } else {
      age2 <- list(Mel = age1$Mel, MelnI = age1$MelnI, Car = age1$Car) #, Lig = age1$Lig
      age1 <- list(Mel = age0$Mel, MelnI = age0$MelnI, Car = age0$Car) #, Lig = age0$Lig
      age0 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #, Lig = NULL
      age0p1 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #, Lig = NULL
      age0p2 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #, Lig = NULL
    }
    
    # Period1 ------------------------------------------------------------------
    # Build-up the colonies
    print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
    print(Sys.time())
    age1 <- list(Mel = buildUp(age1$Mel),
                 MelnI = buildUp(age1$MelnI),
                 Car = buildUp(age1$Car))
    
    if (year > 1) {
      age2 <- list(Mel = buildUp(age2$Mel),
                   MelnI = buildUp(age2$MelnI),
                   Car = buildUp(age2$Car))
      
    }
    
    # Split all age1 colonies
    print("Splitting the colonies")
    print(Sys.time())
    tmp <- list(Mel = split(age1$Mel),
                MelnI = split(age1$MelnI),
                Car = split(age1$Car))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    
    # The queens of the splits are 0 years old
    age0p1 <- list(Mel = tmp$Mel$split, MelnI = tmp$MelnI$split, Car = tmp$Car$split) #,Lig = tmp$Lig$split
    
    if (year > 1) {
      
      tmp <- list(Mel = split(age2$Mel),
                  MelnI = split(age2$MelnI),
                  Car = split(age2$Car))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant
      )
      # The queens of the splits are 0 years old
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                     MelnI = c(age0p1$MelnI, tmp$MelnI$split),
                     Car = c(age0p1$Car, tmp$Car$split))
      #,Lig = c(age0p1$Lig, tmp$Lig$split))
    }
    
    
    # Create virgin queens
    # Sample colony for the virgin queens
    print("Create virgin queens, period 1")
    print(Sys.time())
  
    # Virgin queens for splits!
    
    #requeen with carnica and mellifera for Mel and carnica for Car
    
    pImport<-param1Value
    tmp <- (Mel = pullColonies(age0p1$Mel, p=pImport)) #pull colonies to requeen with imports
    IdImportColonies<-getId(tmp$pulled) #get the ids of the imported colonies
    age0p1 <- list(Mel = tmp$remnant,
                   MelImport = tmp$pulled,
                   MelnI = c(age0p1$MelnI, tmp$MelnI$split),
                   Car = c(age0p1$Car, tmp$Car$split))
    
   
    virginQueens <- list(Mel = createVirginQueens(age1$Mel, collapse=T, nInd=10),
                         MelnI = createVirginQueens(age1$MelnI, collapse=T, nInd=10),
                         Car = createVirginQueens(age1$Car, collapse=T, nInd=10))
    
    virginQueens<-list(Mel = mergePops(virginQueens$Mel),
                       MelnI = mergePops(virginQueens$MelnI),
                       Car = mergePops(virginQueens$Car))
    
    carqueens<-nColonies(age0p1$Car)+nColonies(age0p1$MelImport)
    
    virginQueens<-list(Mel = virginQueens$Mel[sample(1:nInd(virginQueens$Mel), size = nColonies(age0p1$Mel), replace = FALSE)],
                       MelnI =  virginQueens$MelnI[sample(1:nInd(virginQueens$MelnI), size = nColonies(age0p1$MelnI), replace = FALSE)],
                       Car = virginQueens$Car[sample(1:nInd(virginQueens$Car), size = carqueens, replace = FALSE)])
                       
    
    
    #requeen with carnica and mellifera for Mel and carnica for Car
    nColoniesMelImport<-nColonies(age0p1$MelImport)
    nColoniesCar<-nColonies(age0p1$Car)+nColonies(age0p1$MelImport)
    nColoniesMelnI<-nColonies(age0p1$MelnI)
    age0p1 <- list(Mel = c(reQueen(age0p1$Mel, queen = (virginQueens$Mel)) ,
                           reQueen(age0p1$MelImport, queen = c((virginQueens$Car)[1:nColoniesMelImport]))),
                   MelnI = reQueen(age0p1$MelnI, queen = virginQueens$MelnI[1:nColoniesMelnI]),
                   Car = reQueen(age0p1$Car, queen = virginQueens$Car[(nColoniesMelImport+1):nColoniesCar]))
    
    
    # Swarm a percentage of age1 colonies
    print("Swarm colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
                MelnI = pullColonies(age1$MelnI, p = p1swarm),
                Car = pullColonies(age1$Car, p = p1swarm))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                MelnI = swarm(tmp$MelnI$pulled),
                Car = swarm(tmp$Car$pulled))
    
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                   MelnI = c(age0p1$MelnI, tmp$MelnI$remnant),
                   Car = c(age0p1$Car, tmp$Car$remnant))
    
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 MelnI = c(age1$MelnI, tmp$MelnI$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))
    
    
    
    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
                  MelnI = pullColonies(age2$MelnI, p = p1swarm),
                  Car = pullColonies(age2$Car, p = p1swarm))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
                  MelnI = swarm(tmp$MelnI$pulled),
                  Car = swarm(tmp$Car$pulled))
      
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                     MelnI = c(age0p1$MelnI, tmp$MelnI$remnant),
                     Car = c(age0p1$Car, tmp$Car$remnant))
      
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                   MelnI = c(age2$MelnI, tmp$MelnI$swarm),
                   Car = c(age2$Car, tmp$Car$swarm))
      
    }
    
    # Supersede age1 colonies
    print("Supersede colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
                MelnI = pullColonies(age1$MelnI, p = p1supersede),
                Car = pullColonies(age1$Car, p = p1supersede))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                MelnI = supersede(tmp$MelnI$pulled),
                Car = supersede(tmp$Car$pulled))
    
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                   MelnI = c(age0p1$MelnI, tmp$MelnI),
                   Car = c(age0p1$Car, tmp$Car))
    
    
    if (year > 1) {
      # Supersede age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
                  MelnI = pullColonies(age2$MelnI, p = p1supersede),
                  Car = pullColonies(age2$Car, p = p1supersede))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  MelnI = supersede(tmp$MelnI$pulled),
                  Car = supersede(tmp$Car$pulled))
      
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                     MelnI = c(age0p1$MelnI, tmp$MelnI),
                     Car = c(age0p1$Car, tmp$Car))
      
    }
    
    # Mate the split colonies
    print("Mate split colonies, P1")
    print(Sys.time())
    if (year == 1) {
      DCAMel <- createDCA(age1$Mel)
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
      DCAMelnI <- createDCA(age1$MelnI)
      age0p1$MelnI <- cross(age0p1$MelnI, drones = pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p1$MelnI), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
      DCAMelnI <- createDCA(c(age1$MelnI, age2$MelnI)) 
      age0p1$MelnI <- cross(age0p1$MelnI, drones = pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p1$MelnI), nDrones = nFathersPoisson))
    }
    
    # Collapse
    print("Collapse colonies, P1")
    print(Sys.time())
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
                 MelnI = selectColonies(age1$MelnI, p = 1 - p1collapse),
                 Car = selectColonies(age1$Car, p = 1 - p1collapse))
    
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
                   MelnI = selectColonies(age2$MelnI, p = 1 - p1collapse),
                   Car = selectColonies(age2$Car, p = 1 - p1collapse))
      
    }
    
    # Period2 ------------------------------------------------------------------
    print("PERIOD 2")
    # Swarm a percentage of age1 colonies
    # Mellifera
    print("Swarm colonies, P2")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p2swarm),
                MelnI = pullColonies(age1$MelnI, p = p2swarm),
                Car = pullColonies(age1$Car, p = p2swarm))
    
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                MelnI = swarm(tmp$MelnI$pulled),
                Car = swarm(tmp$Car$pulled))
    
    # The queens of the remnant colonies are of age 0
    age0p2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
    
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 MelnI = c(age1$MelnI, tmp$MelnI$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))
    
    
    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
                  MelnI = pullColonies(age2$MelnI, p = p2swarm),
                  Car = pullColonies(age2$Car, p = p2swarm))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
                  MelnI = swarm(tmp$MelnI$pulled),
                  Car = swarm(tmp$Car$pulled))
      
      # The queens of the remnant colonies are of age 0
      age0p2 <- list(Mel = tmp$Mel$remnant,
                     MelnI = tmp$MelnI$remnant,
                     Car = tmp$Car$remnant)
      
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                   MelnI = c(age2$MelnI, tmp$MelnI$swarm),
                   Car = c(age2$Car, tmp$Car$swarm))
      
    }
    
    # Supersede a part of age1 colonies
    print("Supersede colonies, P2")
    print(Sys.time())
    
    tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
                MelnI = pullColonies(age1$MelnI, p = p2supersede),
                Car = pullColonies(age1$Car, p = p2supersede))
    age1 <- list(Mel = tmp$Mel$remnant,
                 MelnI = tmp$MelnI$remnant,
                 Car = tmp$Car$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                MelnI = supersede(tmp$MelnI$pulled),
                Car = supersede(tmp$Car$pulled))
    # The queens of superseded colonies are of age 0
    age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                   MelnI = c(age0p2$MelnI, tmp$MelnI),
                   Car = c(age0p2$Car, tmp$Car))
    
    
    if (year > 1) {
      # Supersede a part of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
                  MelnI = pullColonies(age2$MelnI, p = p2supersede),
                  Car = pullColonies(age2$Car, p = p2supersede))
      age2 <- list(Mel = tmp$Mel$remnant,
                   MelnI = tmp$MelnI$remnant,
                   Car = tmp$Car$remnant)
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  MelnI = supersede(tmp$MelnI$pulled),
                  Car = supersede(tmp$Car$pulled))
      # The queens of superseded colonies are of age 0
      age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                     MelnI = c(age0p2$MelnI, tmp$MelnI),
                     Car = c(age0p2$Car, tmp$Car))
      
    }
    
    # Replace all the drones
    print("Replace Drones, P2")
    print(Sys.time())
    
    age1$Mel <- replaceDrones(age1$Mel)
    age1$MelnI <- replaceDrones(age1$MelnI)
    age1$Car <- replaceDrones(age1$Car)
    
    if (year > 1) {
      age2$Mel <- replaceDrones(age2$Mel)
      age2$MelnI <- replaceDrones(age2$MelnI)
      age2$Car <- replaceDrones(age2$Car)
    }
    
    # Mate the colonies
    # Import p percentage of carnica colonies into mellifera DCA
    print("Mate colonies, P2")
    print(Sys.time())
    
    
    if (year == 1) {
      DCAMel <- createDCA(age1$Mel)
      age0p2$Mel <- cross(age0p2$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
      DCAMelnI <- createDCA(age1$MelnI)
      age0p2$MelnI <- cross(age0p2$MelnI, drones = pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p2$MelnI), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      fathersMel <- pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson)
      age0p2$Mel <- cross(age0p2$Mel, drones = fathersMel)
      DCACar <- createDCA(c(age1$Car, age2$Car))
      fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
      age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
      DCAMelnI <- createDCA(c(age1$MelnI, age2$MelnI))
      fathersMelnI <-  pullDroneGroupsFromDCA(DCA = DCAMelnI, n = nColonies(age0p2$MelnI), nDrones = nFathersPoisson)
      age0p2$MelnI <- cross(age0p2$MelnI, drones = fathersMelnI)
    }
    
    # Collapse
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p2collapse),
                 MelnI = selectColonies(age1$MelnI, p = 1 - p2collapse),
                 Car = selectColonies(age1$Car, p = 1 - p2collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
                   MelnI = selectColonies(age2$MelnI, p = 1 - p2collapse),
                   Car = selectColonies(age2$Car, p = 1 - p2collapse))
    }
    
    # Merge all age 0 colonies (from both periods)
    age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
                 MelnI = c(age0p1$MelnI, age0p2$MelnI),
                 Car = c(age0p1$Car, age0p2$Car))
    
    
    # Period3 ------------------------------------------------------------------
    # Collapse age0 queens
    print("PERIOD 3")
    print("Collapse colonies, P3")
    print(Sys.time())
    
    #selection on Fitness queen gv
    
    #age0
    #Mellifera
    MelBritFit<-sapply(getPheno(age0$Mel, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID<-names(sort(MelBritFit,T))
    Nselectcolon<-round(length(queensID)*(1-p3collapseAge0)) #calculate how many colonies will collapse
    age0MelqueensID<-queensID[1:Nselectcolon] #select the queens ids that will not collapse
    
    MelBritFit1<-sapply(getPheno(age1$Mel, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID1<-names(sort(MelBritFit1,T))
    Nselectcolon1<-round(length(queensID1)*(1-p3collapseAge1)) #calculate how many colonies will collapse
    age1MelqueensID<-queensID1[1:Nselectcolon1] #select the queens ids that will not collapse
    
    MelnIBritFit<-sapply(getPheno(age0$MelnI, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID2<-names(sort(MelnIBritFit,T))
    Nselectcolon2<-round(length(queensID2)*(1-p3collapseAge0)) #calculate how many colonies will collapse
    age0MelnIqueensID<-queensID2[1:Nselectcolon2] #select the queens ids that will not collapse
    
    MelnIBritFit1<-sapply(getPheno(age1$MelnI, caste = "queen"), function(x) x[1,3])#select on Mel Brit fitness
    queensID3<-names(sort(MelnIBritFit1,T))
    Nselectcolon3<-round(length(queensID3)*(1-p3collapseAge1)) #calculate how many colonies will collapse
    age1MelnIqueensID<-queensID3[1:Nselectcolon3] #select the queens ids that will not collapse
    
    CarEuFit<-sapply(getPheno(age0$Car, caste = "queen"), function(x) x[1,4])#select on Mel Brit fitness
    queensID<-names(sort(CarEuFit,T))
    Nselectcolon<-round(length(queensID)*(1-p3collapseAge0)) #calculate how many colonies will collapse
    age0CarqueensID<-queensID[1:Nselectcolon] #select the queens ids that will not collapse
    
    CarEuFit1<-sapply(getPheno(age1$Car, caste = "queen"), function(x) x[1,4])#select on Mel Brit fitness
    queensID1<-names(sort(CarEuFit1,T))
    Nselectcolon1<-round(length(queensID1)*(1-p3collapseAge1)) #calculate how many colonies will collapse
    age1CarqueensID<-queensID1[1:Nselectcolon1] #select the queens ids that will not collapse
    
    #Collapse 
    Carqueenpre<-nColonies(selectColonies(age0$Mel,ID= IdImportColonies))
    Carqueenpre1<-nColonies(selectColonies(age1$Mel,ID= IdImportColonies))
    age0 <- list(Mel = selectColonies(age0$Mel, ID = age0MelqueensID),
                 MelnI = selectColonies(age0$MelnI, ID = age0MelnIqueensID),
                 Car = selectColonies(age0$Car, ID = age0CarqueensID))
    
    age1 <- list(Mel = selectColonies(age1$Mel, ID = age1MelqueensID),
                 MelnI = selectColonies(age1$MelnI, ID = age1MelnIqueensID),
                 Car = selectColonies(age1$Car, ID= age1CarqueensID))
    
    age2 <- list(Mel = NULL, MelnI = NULL, Car = NULL) #We don't need this but just to show the workflow!!!
    Carqueenpost<-nColonies(selectColonies(age0$Mel,ID= IdImportColonies))
    Carqueenpost1<-nColonies(selectColonies(age1$Mel,ID= IdImportColonies))
    deadqueens<-(Carqueenpre-Carqueenpost)+(Carqueenpre1-Carqueenpost1)
    print(deadqueens)
    
    
    # Maintain the number of colonies ------------------------------------------
    # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
    print("Maintain the number, P2")
    print(Sys.time())
    
    age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
    age0$MelnI <- maintainIrelandSize(age0 = age0$MelnI, age1 = age1$MelnI)
    age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)
    
    
    for (subspecies in c("Mel", "MelnI", "Car")) {     
      if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == IrelandSize
          | (nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == CarSize)
      {
      } 
      else 
      {stop(paste0("The number of colonies for ", subspecies, " does not match the population size!"))}
    }
    
    #track mean IBD and variance of IBD
    
    #this creates two empty dataframes where all the information of each year and each replica will be recorded
    
    
    columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","sdHoneyYieldBrit","HoneyYieldEu","sdHoneyYieldEu",
                     "FitnessBrit","sdFitnessBrit","FitnessEu","sdFitnessEu","pheHoneyYieldBrit","sdpheHoneyYieldBrit","pheHoneyYieldEu","sdpheHoneyYieldEu",
                     "pheFitnessBrit","sdpheFitnessBrit","pheFitnessEu","sdpheFitnessEu","Homocigosity","sdHomocigosity","sdIBD")
    
    if (year==1){ 
      columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","sdHoneyYieldBrit","HoneyYieldEu","sdHoneyYieldEu",
                       "FitnessBrit","sdFitnessBrit","FitnessEu","sdFitnessEu","pheHoneyYieldBrit","sdpheHoneyYieldBrit","pheHoneyYieldEu","sdpheHoneyYieldEu",
                       "pheFitnessBrit","sdpheFitnessBrit","pheFitnessEu","sdpheFitnessEu","Homocigosity","sdHomocigosity","sdIBD")
      MeanVarMel <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for mellifera
      colnames(MeanVarMel)<-columnheaders
      MeanVarCar <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for carnica
      colnames(MeanVarCar)<-columnheaders
      MeanVarMelnI <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for carnica
      colnames(MeanVarMelnI)<-columnheaders
    } 
    
    if (year==1){
      uno<-0
    }else{
      uno<-nrow(colonyRecords) #this is to show where to start recording values on the dataframe
    }
    
    
    #record values of Mellifera population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel", Rep=Rep)
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
                         sdIBD=sd(colonyRecords[(uno+1):dos,"IBD"]))
    
    #record values for Mellifera without imports population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$MelnI, year = year, population = "MelnI", Rep=Rep)
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$MelnI, year = year, population = "MelnI", Rep=Rep)
    tres<-nrow(colonyRecords)#this is to show where to end recording values on the dataframe
    
    
    #create a dataframe with the mellifera mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
    newrow2<- data.frame(MeanIBD=mean(colonyRecords[(dos+1):tres,"IBD"]), VarIBD=var(colonyRecords[(dos+1):tres,"IBD"]), Year=year,Rep=Rep,Population="MelnI",
                         HoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]), sdHoneyYieldBrit=sd(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]),
                         FitnessBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritFit"]),sdFitnessBrit=sd(colonyRecords[(dos+1):tres,"gvQueens_BritFit"]),
                         HoneyYieldEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]), sdHoneyYieldEu=sd(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]),
                         FitnessEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuFit"]), sdFitnessEu=sd(colonyRecords[(dos+1):tres,"gvQueens_EuFit"]),
                         pheHoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"pheQueens_BritHY"]), sdpheHoneyYieldBrit=sd(colonyRecords[(dos+1):tres,"pheQueens_BritHY"]),
                         pheFitnessBrit=mean(colonyRecords[(dos+1):tres,"pheQueens_BritFit"]),sdpheFitnessBrit=sd(colonyRecords[(dos+1):tres,"pheQueens_BritFit"]),
                         pheHoneyYieldEu=mean(colonyRecords[(dos+1):tres,"pheQueens_EuHY"]), sdpheHoneyYieldEu=sd(colonyRecords[(dos+1):tres,"pheQueens_EuHY"]),
                         pheFitnessEu=mean(colonyRecords[(dos+1):tres,"pheQueens_EuFit"]), sdpheFitnessEu=sd(colonyRecords[(dos+1):tres,"pheQueens_EuFit"]),
                         Homocigosity=mean(colonyRecords[(dos+1):tres,"pHomBrood"]),sdHomocigosity=sd(colonyRecords[(dos+1):tres,"pHomBrood"]),
                         sdIBD=sd(colonyRecords[(dos+1):tres,"IBD"]))
    
    #record values of Carnica population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car", Rep=Rep)
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car", Rep=Rep)
    cua<-nrow(colonyRecords)
    
    #create a dataframe with the carnica mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
    newrow3<- data.frame(MeanIBD=mean(colonyRecords[(tres+1):cua,"IBD"]), VarIBD=var(colonyRecords[(tres+1):cua,"IBD"]), Year=year,Rep=Rep,Population="Car",
                         HoneyYieldBrit=mean(colonyRecords[(tres+1):cua,"gvQueens_BritHY"]), sdHoneyYieldBrit=sd(colonyRecords[(tres+1):cua,"gvQueens_BritHY"]),
                         FitnessBrit=mean(colonyRecords[(tres+1):cua,"gvQueens_BritFit"]),sdFitnessBrit=sd(colonyRecords[(tres+1):cua,"gvQueens_BritFit"]),
                         HoneyYieldEu=mean(colonyRecords[(tres+1):cua,"gvQueens_EuHY"]), sdHoneyYieldEu=sd(colonyRecords[(tres+1):cua,"gvQueens_EuHY"]),
                         FitnessEu=mean(colonyRecords[(tres+1):cua,"gvQueens_EuFit"]), sdFitnessEu=sd(colonyRecords[(tres+1):cua,"gvQueens_EuFit"]),
                         pheHoneyYieldBrit=mean(colonyRecords[(tres+1):cua,"pheQueens_BritHY"]), sdpheHoneyYieldBrit=sd(colonyRecords[(tres+1):cua,"pheQueens_BritHY"]),
                         pheFitnessBrit=mean(colonyRecords[(tres+1):cua,"pheQueens_BritFit"]),sdpheFitnessBrit=sd(colonyRecords[(tres+1):cua,"pheQueens_BritFit"]),
                         pheHoneyYieldEu=mean(colonyRecords[(tres+1):cua,"pheQueens_EuHY"]), sdpheHoneyYieldEu=sd(colonyRecords[(tres+1):cua,"pheQueens_EuHY"]),
                         pheFitnessEu=mean(colonyRecords[(tres+1):cua,"pheQueens_EuFit"]), sdpheFitnessEu=sd(colonyRecords[(tres+1):cua,"pheQueens_EuFit"]),
                         Homocigosity=mean(colonyRecords[(tres+1):cua,"pHomBrood"]),sdHomocigosity=sd(colonyRecords[(tres+1):cua,"pHomBrood"]),
                         sdIBD=sd(colonyRecords[(tres+1):cua,"IBD"]))

    #Combine what we had in each dataframe with the new info, so each year the dataframe updates with new values
    MeanVarMel<-rbind(MeanVarMel,newrow1)
    MeanVarMelnI<-rbind(MeanVarMelnI,newrow2)
    MeanVarCar<-rbind(MeanVarCar,newrow3)
  } #end of year loop
  CombinedDf<-rbind(MeanVarCar,MeanVarMel,MeanVarMelnI)
  save.image(paste0("Selection",Rep,"_import",param1Value,".RData"))
  
  sink(file = paste0("Selection",Rep,"_import",param1Value,".RData")) # saving data from the burnin
  cat(save.image(paste0("Selection",Rep,"_import",param1Value,".RData")))
  sink()
  write.csv(CombinedDf, paste0("scenarioData_rep",Rep,"_import",param1Value,".csv"))
  sink(file = paste0("scenarioData_rep",Rep,"_import",param1Value,".csv")) # TODO: convert this to a CSV export
  cat(write.csv(CombinedDf, paste0("scenarioData_rep",Rep,"_import",param1Value,".csv")))
  sink()
  remove(MeanVarCar)
  remove(MeanVarMelnI) 
  remove(MeanVarMel)
  setwd("..") # scenario folder
  setwd("..") # rep folder
  
}
