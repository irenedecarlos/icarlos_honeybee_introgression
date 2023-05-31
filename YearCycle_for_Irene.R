#setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Clean workspace
rm(list = ls())

# Define functions
maintainIrelandSize <- function(age0 = NULL, age1 = NULL) {
  if ((nColonies(age0) + nColonies(age1)) > IrelandSize) { # check if the sum of all colonies is greater than apiary size
    IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
    splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
    age0split <- splits0$pulled # create an object for age 0 splits
    age0swarm <- splits0$remnant # create an object for swarms and superseded colonies
    age0needed <- IrelandSize - nColonies(age1) # calculate the number of age 0 colonies that are needed to fill up the apiary
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
library(tictoc)
library(R6)
library(nadiv)
library(Matrix)
library(SIMplyBee)
library(dplyr)
library(tidyr)
# TODO: replace with devtools installation from Github once the package is operational
# Source the development version of AlphaSimR
getwd()

# Founder population parameters -------------------------------------------------------------------
nMelN = 450                   # Number of Mellifera
nCar = 150                    # Number of Carnica
nChr = 1                     # Number of chromomsome
nDronesPerQueen = 50
nSegSites = 100              # Number of segregating sites

# Population parameters -------------------------------------------------------------------
nRep <- 1                     # Number of repeats
nYear <- 10                   # Number of years
#apiarySize <- 300             # Number of colonies in the apiary
IrelandSize<-300            #remove apiary size from code
CarSize<-100
nWorkers <- 10                # Number of workers in a full colony
nDrones <- 50                 # Number of drones in a full colony (typically nWorkers * 0.2 (not in the example))
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
pImport <- 0.5              # Percentage import from carnica to mellifera

# Create data frames for recording the number of age0 and age1 colonies, csd variability and for recording cpu time
loopTime <- data.frame(Rep = NA, tic = NA, toc = NA, msg = NA, time = NA)

# Prepare recording function
data_rec <- function(datafile, colonies, year, population) {
  queens = mergePops(getQueen(colonies))
  datafile = rbind(datafile,
                   data.frame(colonies             = deparse(substitute(colonies)),
                              population           = population,
                              year                 = year,
                              Id                   = queens@id,
                              MId                  = queens@mother,
                              FId                  = queens@father,
                              nFathers             = nFathers(queens),
                              nDPQ                 = sapply(getFathers(queens), function(x) length(unique(x@mother))),
                              nCsdAlColony         = sapply(colonies@colonies, function(x) nCsdAlleles(x, collapse = TRUE)),
                              nCsdApiary           = rep(nCsdAlleles(colonies, collapse = TRUE), queens@nInd),
                              pHomBrood            = calcQueensPHomBrood(queens),
                              gvQueens_QueenTrait  = sapply(getGv(colonies, caste = "queen"), function(x) x[1,1]),
                              gvQueens_WorkerTrait = sapply(getGv(colonies, caste = "queen"), function(x) x[1,2])
                   ))}
colonyRecords = NULL

# Start of the rep-loop ---------------------------------------------------------------------
for (Rep in 1:nRep) {
  # Rep <- 1 (you can use this to check if your code is running alright for one whole loop)
  cat(paste0("Rep: ", Rep, "/", nRep, "\n"))
  tic(paste0(nYear, 'y loop'))         # Measure cpu time
  Rprof()                              # Start profiling


  # Founder population ---------------------------------------------------------
  # STEP 1:  Create a founder population of A. m. mellifera and A. m. carnica bees (un-# the one you want to use)

  # Using simulateHoneyBeeGenomes is a great way to get a basic founder gene pool
  # founderGenomes <- simulateHoneyBeeGenomes(nMelN = nMelN,
  #                                           nCar = nCar,
  #                                           nChr = nChr,
  #                                           nSegSites = nSegSites)
  #
  # save(founderGenomes, file="founderGenomes_ThreePop.RData")

  # Or you can load your previously made founder population
  # print("Loading in the founderData")
  # load("FounderGenomes_ThreePop_16chr.RData")
  # load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/FounderGenomes_ThreePop_16chr.RData")

# quick haplo to get the founder genomes for now.
founderGenomes<- quickHaplo(sum(nMelN,nCar),1,segSites = 100)
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
  # Skip this for now
  # Add traits - taken from the QuantGen vignette 
  mean <- c(20, 0)
  varA <- c(1, 1 / SP$nWorkers)
  corA <- matrix(data = c( 1.0, -0.5,
                           -0.5,  1.0), nrow = 2, byrow = TRUE)
  SP$addTraitA(nQtlPerChr = 100, mean = mean, var = varA, corA = corA,
               name = c("queenTrait", "workersTrait"))

  varE <- c(3, 3 / SP$nWorkers)

  # TODO: what is a reasonable environmental correlation between queen and worker effects?
  corE <- matrix(data = c(1.0, 0.3,
                          0.3, 1.0), nrow = 2, byrow = T)
  
  SP$setVarE(varE = varE, corE = corE)

  
  # STEP 3: Set up your base population
  # Create a base population for A. m. mellifera, A. m. mellifera cross, and A. m. carnica (400 of each)
  virginQueens <- list(Mel = createVirginQueens(x = founderGenomes[1:(nMelN)]),
                       Car = createVirginQueens(x = founderGenomes[(nMelN +1):(nMelN + nCar)]))
  # Create drones for A. m. mellifera, A. m. mellifera cross, and A. m. carnica
  drones <- list(Mel = createDrones(x = virginQueens$Mel[(IrelandSize+1):(nMelN)], nInd = nDronesPerQueen),
                 Car = createDrones(x = virginQueens$Car[(CarSize+1):nCar], nInd = nDronesPerQueen))
  # Get fathers for Mel, MelCross and Car
  fathersMel <- pullDroneGroupsFromDCA(drones$Mel, n = nInd(virginQueens$Mel[1:IrelandSize]), nDrones = nFathersPoisson)
  fathersCar <- pullDroneGroupsFromDCA(drones$Car, n = nInd(virginQueens$Car[1:CarSize]), nDrones = nFathersPoisson)

  # Mate virgin queens with fathers to make them queens
  queens <- list(Mel = SIMplyBee::cross(x = virginQueens$Mel[1:IrelandSize], drones = fathersMel),
                 Car = SIMplyBee::cross(x = virginQueens$Car[1:CarSize], drones = fathersCar))

  #skip this
  #Set allele frequency for queens
  tmp <- c(virginQueens$Mel, virginQueens$Car)
  alleleFreqBaseQueens <- calcBeeAlleleFreq(x = getSegSiteGeno(tmp),
                                            sex = tmp@sex)

  alleleFreqBaseQueensCar <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Car),
                                               sex = virginQueens$Car@sex)

  alleleFreqBaseQueensMel <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Mel),
                                               sex = virginQueens$Mel@sex)
  #Get allele freq for csd locus
  csdLocus <- paste0(SP$csdChr, "_", SP$csdPosStart:SP$csdPosStop)
  alleleFreqCsdLocusBaseQueens <- alleleFreqBaseQueens[csdLocus]
  alleleFreqCsdLocusBaseCar <- alleleFreqBaseQueensCar[csdLocus]
  alleleFreqCsdLocusBaseMel <- alleleFreqBaseQueensMel[csdLocus]

  #Get allele freq for csd Chromosome - this pulls out only the 3rd chromosome
  alleleFreqCsdChrBaseQueens <- t(as.data.frame(alleleFreqBaseQueens))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueens))))] %>% t()
  alleleFreqCsdChrBaseCar <- t(as.data.frame(alleleFreqBaseQueensCar))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensCar))))] %>% t()
  alleleFreqCsdChrBaseMel <- t(as.data.frame(alleleFreqBaseQueensMel))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensMel))))] %>% t()

year=1
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
                   Car = createMultiColony(x = queens$Car, n = CarSize))
      print("Record initial colonies")
      colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Mel, year = year, population = "Mel")
      colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car")

      # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
    } else {
      age2 <- list(Mel = age1$Mel, Car = age1$Car)
      age1 <- list(Mel = age0$Mel, Car = age0$Car)
      age0 <- list(Mel = NULL, Car = NULL)
      age0p1 <- list(Mel = NULL, Car = NULL)
      age0p2 <- list(Mel = NULL, Car = NULL)
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

    if (year > 1) {
      # Split all age2 colonies
      tmp <- list(Mel = split(age2$Mel),
                  Car = split(age2$Car))
      age2 <- list(Mel = tmp$Mel$remnant,
                   Car = tmp$Car$remnant)
      # The queens of the splits are 0 years old
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                     Car = c(age0p1$Car, tmp$Car$split))
    }

    # Create virgin queens
    # Sample colony for the virgin queens
    print("Create virgin queens, period 1")
    print(Sys.time())
    virginDonor <- list(Mel = sample.int(n = nColonies(age1$Mel), size = 1),
                        Car = sample.int(n = nColonies(age1$Car), size = 1))
    # Virgin queens for splits!
    virginQueens <- list(Mel = createVirginQueens(age1$Mel[[virginDonor$Mel]], nInd = nColonies(age0p1$Mel)),
                         Car = createVirginQueens(age1$Car[[virginDonor$Car]], nInd = nColonies(age0p1$Car)))

    # Requeen the splits --> queens are now 0 years old
    age0p1 <- list(Mel = reQueen(age0p1$Mel, queen = virginQueens$Mel),
                   Car = reQueen(age0p1$Car, queen = virginQueens$Car))

    # Swarm a percentage of age1 colonies
    print("Swarm colonies, P1")
    print(Sys.time())
    tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
                Car = pullColonies(age1$Car, p = p1swarm))
    age1 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                Car = swarm(tmp$Car$pulled))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                   Car = c(age0p1$Car, tmp$Car$remnant))
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))


    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
                  Car = pullColonies(age2$Car, p = p1swarm))
      age2 <- list(Mel = tmp$Mel$remainingColonies,
                   Car = tmp$Car$remnant)
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
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
      DCAMel <- createDCA(age1$Mel)
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      #DCAMelCross <- createDCA(c(age1$MelCross,
                                 #selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)))) #Is this 0 the P argument?
      #age0p1$MelCross <- cross(age0p1$MelCross, drones = pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p1$MelCross), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      #DCAMelCross <- createDCA(c(age1$MelCross,
                                # selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 #age2$MelCross,
                                 #selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      #age0p1$MelCross <- cross(age0p1$MelCross, drones = pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p1$MelCross), nDrones = nFathersPoisson))
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
    }

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
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
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
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
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
      DCAMel <- createDCA(age1$Mel)
      age0p2$Mel <- cross(age0p2$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson))
      #DCAMelCross <- createDCA(c(age1$MelCross,
                                 #selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0))))
      #age0p2$MelCross <- cross(age0p2$MelCross, drones = pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p2$MelCross), nDrones = nFathersPoisson))
      DCACar <- createDCA(age1$Car)
      age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      fathersMel <- pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson)
      fathersMel[[1]] <- c(fathersMel[[1]], createDrones(age1$Mel[[1]], nInd = 2))
      age0p2$Mel <- cross(age0p2$Mel, drones = fathersMel)
      #DCAMelCross <- createDCA(c(age1$MelCross,
                                 #selectColonies(age1$Car, n = round(nColonies(age1$MelCross) * pImport, 0)),
                                 #age2$MelCross,
                                 #selectColonies(age2$Car, n = round(nColonies(age2$MelCross) * pImport, 0))))
      #fathersMelCross <- pullDroneGroupsFromDCA(DCA = DCAMelCross, n = nColonies(age0p2$MelCross), nDrones = nFathersPoisson)
      #fathersMelCross[[1]] <-  c(fathersMelCross[[1]], createDrones(age1$MelCross[[1]], nInd = 2))
      #age0p2$MelCross <- cross(age0p2$MelCross, drones = fathersMelCross)
      DCACar <- createDCA(c(age1$Car, age2$Car))
      fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
      fathersCar[[1]] <- c(fathersCar[[1]], createDrones(age1$Car[[1]], nInd = 2))
      age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
    }

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
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel")
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car")

    # Period3 ------------------------------------------------------------------
    # Collapse age0 queens
    print("PERIOD 3")
    print("Collapse colonies, P3")
    print(Sys.time())

    age0 <- list(Mel = selectColonies(age0$Mel, p = (1 - p3collapseAge0)),
                 Car = selectColonies(age0$Car, p = (1 - p3collapseAge0)))
    age1 <- list(Mel = selectColonies(age1$Mel, p = (1 - p3collapseAge1)),
                 Car = selectColonies(age1$Car, p = (1 - p3collapseAge1)))
    age2 <- list(Mel = NULL, MelCross = NULL, Car = NULL) #We don't need this but just to show the workflow!!!


    # Maintain the number of colonies ------------------------------------------
    # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
    print("Maintain the number, P2")
    print(Sys.time())

    age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
    age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)

    for (subspecies in c("Mel", "MelCross", "Car")) {
      if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) != apiarySize) {
        stop(paste0("The number of colonies for ", subspecies, " does not match the apiary size!"))
      }
    }

  } # Year-loop

  a <- toc()
  loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))

} # Rep-loop

print("Saving image data")
save.image("SpringerSimulation_import.RData")








