#setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Clean workspace

rm(list = ls())
getwd()

# Define functions
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
      gvCarQueens_EuHY <- sapply(getGv(age0split, caste = "queen"), function(x) x[1,3]) #get gv for fitness of age 0 carnica queens
      queensID<-names(sort(gvCarQueens_EuHY,decreasing=T))
      age0CarqueensID<-queensID[1:nSplitsNeeded]
      splitTmp <- pullColonies(age0split, ID = age0CarqueensID) # pull the splits
      splits <- splitTmp$pulled # select pulled splits
      age0 <- c(age0swarm, splits) # combine splits and swarms in age 0 object
    }
    return(age0)
  }
}


maintainIrelandSizes <- function(age0 = NULL, age1 = NULL) {
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
      queens <- mergePops(getQueen(age0split))#get the queens from the split
      IBDh <- apply(getIbdHaplo(queens),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))#get IBD for each haplotype
      IBD = sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)#IBd for every queen
      MelHY<-sapply(getGv(age0split, caste = "queen"), function(x) x[1,1])# get the gvHoneyYield for mellifera
      MelHY<-MelHY*IBD #Multiply the honey yield by the IBD
      MelImpHY<-sapply(getGv(age0split, caste = "queen"), function(x) x[1,2])#get the honey yield for imported carnica
      MelImpHY<-MelImpHY*(1-IBD)#Multiply the honey yield by the IBD of carnica
      MelPopHY<-MelHY+MelImpHY #sum both so that I have the hybridised mellifera, and imported carnica HY 
      queensID<-names(sort(MelPopHY,decreasing=T))#order the IDs of the queens with higher HY
      age0MelqueensID<-queensID[1:nSplitsNeeded] # select ids of splits
      splitTmp <- pullColonies(age0split, ID = age0MelqueensID) # pull the splits
      splits <- splitTmp$pulled # select pulled splits
      age0 <- c(age0swarm, splits) #Colonies that will be used to mantain size
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
library(tictoc)
library(R6)
library(nadiv)
library(Matrix)
library(SIMplyBee)
library(dplyr)
library(ggpubr)

#library(tidyr)
# TODO: replace with devtools installation from Github once the package is operational
# Source the development version of AlphaSimR
getwd()

# Founder population parameters -------------------------------------------------------------------
nMelN = 300                  # Number of Mellifera
nCar =300                   # Number of Carnica
#nLig = 150
nChr = 1                     # Number of chromomsome
nDronesPerQueen = 50
nSegSites = 100              # Number of segregating sites

# Population parameters -------------------------------------------------------------------
nRep <- 1                     # Number of repeats
nYear <- 10                    # Number of years
IrelandSize<-200              #Ireland population size
CarSize<-200                  #Carnica pop size
#LigSize<-100                  #Ligustica pop size
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
pImport <- 0.3              # Percentage import from carnica to mellifera

# Create data frames for recording the number of age0 and age1 colonies, csd variability and for recording cpu time
loopTime <- data.frame(Rep = NA, tic = NA, toc = NA, msg = NA, time = NA)

# Prepare recording function
data_rec <- function(datafile, colonies, year, population, Rep) {
  queens = mergePops(getQueen(colonies))
  IBDh = apply(getIbdHaplo(queens),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
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
                              gvQueens_BritFit  = sapply(getGv(colonies, caste = "queen"), function(x) x[1,2]),
                              gvQueens_EuHY    = sapply(getGv(colonies, caste = "queen"), function(x) x[1,3]),
                              gvQueens_EuFit    = sapply(getGv(colonies, caste = "queen"), function(x) x[1,4]),
                              IBD = sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
                              
                              
                   ))}
colonyRecords = NULL

Rep=1
nRep=5
# Start of the rep-loop ---------------------------------------------------------------------
for (Rep in 1:nRep) {
 
  #dir.create(paste0("C:/Users/Usuario/Honeybeeimports/Replica",Rep))
  #setwd(paste0("C:/Users/Usuario/Honeybeeimports/Replica",Rep))
 
  
  dir.create(paste0("C:/Users/Usuario/Desktop/Máster/Bees/Testeando/Replica",Rep))
  setwd(paste0("C:/Users/Usuario/Desktop/Máster/Bees/Testeando/Replica",Rep))
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
  #FounderGenomes<-load("FounderGenomes_ThreePop_16chr.RData")
  # load("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation/PlottingData/FounderGenomes_ThreePop_16chr.RData")
  
  # quick haplo to get the founder genomes for now.
  founderGenomes<- quickHaplo(sum(nMelN,nCar),4,segSites = 1000)
  nMelN<-300
  
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
  corA <- matrix(data = c(  1.0, 0.0,  0.0, 0.0, 
                            0.0, 1.0,  0.0, 0.0, 
                            0.0, 0.0,  1.0, 0.0, 
                            0.0, 0.0,  0.0, 1.0), nrow = 4, byrow = TRUE)
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
  virginQueens <- list(Mel = createVirginQueens(x = founderGenomes[1:(nMelN)]),
                       Car = createVirginQueens(x = founderGenomes[(nMelN+1):(nMelN + nCar)]))
  SP$nTraits
  
  #,Lig = createVirginQueens(x = founderGenomes[(nMelN+nCar+1):(nMelN+nCar+nLig)])
  
  # Create drones for A. m. mellifera, A. m. mellifera cross, and A. m. carnica
  drones <- list(Mel = createDrones(x = virginQueens$Mel[(IrelandSize+1):(nMelN)], nInd = nDronesPerQueen),
                 Car = createDrones(x = virginQueens$Car[(CarSize+1):nCar], nInd = nDronesPerQueen))
  #,Lig = createDrones(x = virginQueens$Car[(LigSize+1):nLig], nInd = nDronesPerQueen)
  # Get fathers for Mel, MelCross and Car
  fathersMel <- pullDroneGroupsFromDCA(drones$Mel, n = nInd(virginQueens$Mel[1:IrelandSize]), nDrones = nFathersPoisson)
  fathersCar <- pullDroneGroupsFromDCA(drones$Car, n = nInd(virginQueens$Car[1:CarSize]), nDrones = nFathersPoisson)
  #fathersLig <- pullDroneGroupsFromDCA(drones$Lig, n = nInd(virginQueens$Lig[1:LigSize]), nDrones = nFathersPoisson)
  # Mate virgin queens with fathers to make them queens
  queens <- list(Mel = SIMplyBee::cross(x = virginQueens$Mel[1:IrelandSize], drones = fathersMel),
                 Car = SIMplyBee::cross(x = virginQueens$Car[1:CarSize], drones = fathersCar))
  #,Lig = SIMplyBee::cross(x = virginQueens$Lig[1:LigSize], drones = fathersLig)
  
  #QueenHaploCodes<-unique(getIbdHaplo(queens$Mel))
  #hapl<-as.data.frame.array(QueenHaploCodes,row.names=NULL)
  #HaploCodes<-hapl[,1]
  
  #Set allele frequency for queens
  tmp <- c(virginQueens$Mel, virginQueens$Car) #, virginQueens$Lig
  
  alleleFreqBaseQueens <- calcBeeAlleleFreq(x = getSegSiteGeno(tmp),
                                            sex = tmp@sex)
  
  alleleFreqBaseQueensCar <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Car),
                                               sex = virginQueens$Car@sex)
  
  #alleleFreqBaseQueensLig <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Lig),
  #sex = virginQueens$Lig@sex)
  
  alleleFreqBaseQueensMel <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Mel),
                                               sex = virginQueens$Mel@sex)
  
  #Get allele freq for csd locus
  csdLocus <- paste0(SP$csdChr, "_", SP$csdPosStart:SP$csdPosStop)
  alleleFreqCsdLocusBaseQueens <- alleleFreqBaseQueens[csdLocus]
  alleleFreqCsdLocusBaseCar <- alleleFreqBaseQueensCar[csdLocus]
  #alleleFreqCsdLocusBaseLig <- alleleFreqBaseQueensLig[csdLocus]
  alleleFreqCsdLocusBaseMel <- alleleFreqBaseQueensMel[csdLocus]
  
  #Get allele freq for csd Chromosome - this pulls out only the 3rd chromosome
  alleleFreqCsdChrBaseQueens <- t(as.data.frame(alleleFreqBaseQueens))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueens))))] %>% t()
  alleleFreqCsdChrBaseCar <- t(as.data.frame(alleleFreqBaseQueensCar))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensCar))))] %>% t()
  #alleleFreqCsdChrBaseLig <- t(as.data.frame(alleleFreqBaseQueensLig))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensLig))))] %>% t()
  alleleFreqCsdChrBaseMel <- t(as.data.frame(alleleFreqBaseQueensMel))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensMel))))] %>% t()
  
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
                   Car = createMultiColony(x = queens$Car, n = CarSize))
     
      # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
    } else {
      age2 <- list(Mel = age1$Mel, Car = age1$Car) #, Lig = age1$Lig
      age1 <- list(Mel = age0$Mel, Car = age0$Car) #, Lig = age0$Lig
      age0 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
      age0p1 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
      age0p2 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
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
    age0p1 <- list(Mel = tmp$Mel$split, Car = tmp$Car$split) #,Lig = tmp$Lig$split
    
    if (year > 1) {
      # Split all age2 colonies
      tmp <- list(Mel = split(age2$Mel),
                  Car = split(age2$Car))
      
      age2 <- list(Mel = tmp$Mel$remnant,
                   Car = tmp$Car$remnant
      ) 
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
    
    if (year<11){
      pImport<-0
    }else {
      pImport<-0.3
    }
   
    tmp <- (Mel = pullColonies(age0p1$Mel, p=pImport)) #pull colonies to requeen with imports
    IdImportColonies<-getId(tmp$pulled) #get the ids of the imported colonies
    age0p1 <- list(Mel = tmp$remnant,
                   MelImport = tmp$pulled,
                   Car = c(age0p1$Car, tmp$Car$split))
   
    virginQueens <- list(Mel = createVirginQueens(age1$Mel[[virginDonor$Mel]], nInd = nColonies(age0p1$Mel)),
                         Car = createVirginQueens(age1$Car[[virginDonor$Car]], nInd = nColonies(age0p1$Car)+nColonies(age0p1$MelImport)))
   
  
    #requeen with carnica and mellifera for Mel and carnica for Car
    nColoniesMelImport<-nColonies(age0p1$MelImport)
    nColoniesCar<-nColonies(age0p1$Car)+nColonies(age0p1$MelImport)
    
    age0p1 <- list(Mel = c(reQueen(age0p1$Mel, queen = (virginQueens$Mel)) ,
                           reQueen(age0p1$MelImport, queen = c((virginQueens$Car)[1:nColoniesMelImport]))),      
                   Car = reQueen(age0p1$Car, queen = virginQueens$Car[(nColoniesMelImport+1):nColoniesCar]))
   
    
    
    
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
      
      age2 <- list(Mel = tmp$Mel$remnant,     
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
    #,Lig = pullColonies(age1$Lig, p = p1supersede))
    age1 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    #,Lig = tmp$Lig$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                Car = supersede(tmp$Car$pulled))
    #,Lig = supersede(tmp$Lig$pulled))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                   Car = c(age0p1$Car, tmp$Car))
    
    
    if (year > 1) {
      # Supersede age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
                  Car = pullColonies(age2$Car, p = p1supersede))
      #,Lig = pullColonies(age2$Lig, p = p1supersede))
      age2 <- list(Mel = tmp$Mel$remnant,
                   Car = tmp$Car$remnant)
      #,Lig = tmp$Lig$remnant)
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  Car = supersede(tmp$Car$pulled))
      #,Lig = supersede(tmp$Lig$pulled))
      age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
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
      #DCALig <- createDCA(age1$Lig)
      #age0p1$Lig <- cross(age0p1$Lig, drones = pullDroneGroupsFromDCA(DCA = DCALig, n = nColonies(age0p1$Lig), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
      DCACar <- createDCA(c(age1$Car, age2$Car))
      age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
      #DCALig <- createDCA(c(age1$Lig, age2$Lig)) 
      #age0p1$Lig <- cross(age0p1$Lig, drones = pullDroneGroupsFromDCA(DCA = DCALig, n = nColonies(age0p1$Lig), nDrones = nFathersPoisson))
    }
    
    # Collapse
    print("Collapse colonies, P1")
    print(Sys.time())
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
                 Car = selectColonies(age1$Car, p = 1 - p1collapse))
    #,Lig = selectColonies(age1$Lig, p = 1 - p1collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
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
                Car = pullColonies(age1$Car, p = p2swarm))
    #,Lig = pullColonies(age1$Lig, p = p2swarm))
    age1 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    #,Lig = tmp$Lig$remnant)
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                Car = swarm(tmp$Car$pulled))
    #,Lig = swarm(tmp$Lig$pulled))
    # The queens of the remnant colonies are of age 0
    age0p2 <- list(Mel = tmp$Mel$remnant,
                   Car = tmp$Car$remnant)
    #,Lig = tmp$Lig$remnant)
    age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
                 Car = c(age1$Car, tmp$Car$swarm))
    #,Lig = c(age1$Lig, tmp$Lig$swarm))
    
    if (year > 1) {
      # Swarm a percentage of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
                  Car = pullColonies(age2$Car, p = p2swarm))
      #,Lig = pullColonies(age2$Lig, p = p2swarm))
      age2 <- list(Mel = tmp$Mel$remnant,
                   Car = tmp$Car$remnant)
      #,Lig = tmp$Lig$remnant)
      tmp <- list(Mel = swarm(tmp$Mel$pulled),
                  Car = swarm(tmp$Car$pulled))
      #,Lig = swarm(tmp$Lig$pulled))
      # The queens of the remnant colonies are of age 0
      age0p2 <- list(Mel = tmp$Mel$remnant,
                     Car = tmp$Car$remnant)
      #,Lig = tmp$Lig$remnant)
      age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                   Car = c(age2$Car, tmp$Car$swarm))
      #,Lig = c(age2$Lig, tmp$Lig$swarm))
    }
    
    # Supersede a part of age1 colonies
    print("Supersede colonies, P2")
    print(Sys.time())
    
    tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
                Car = pullColonies(age1$Car, p = p2supersede))
    #,Lig = pullColonies(age1$Lig, p = p2supersede))
    age1 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    #,Lig = tmp$Lig$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                Car = supersede(tmp$Car$pulled))
    #,Lig = supersede(tmp$Lig$pulled))
    # The queens of superseded colonies are of age 0
    age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                   Car = c(age0p2$Car, tmp$Car))
    #,Lig = c(age0p2$Lig, tmp$Lig))
    
    if (year > 1) {
      # Supersede a part of age2 colonies
      tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
                  Car = pullColonies(age2$Car, p = p2supersede))
      #,Lig = pullColonies(age2$Lig, p = p2supersede))
      age2 <- list(Mel = tmp$Mel$remnant,
                   Car = tmp$Car$remnant)
      #Lig = tmp$Lig$remnant)
      tmp <- list(Mel = supersede(tmp$Mel$pulled),
                  Car = supersede(tmp$Car$pulled))
      #Lig = supersede(tmp$Lig$pulled))
      # The queens of superseded colonies are of age 0
      age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                     Car = c(age0p2$Car, tmp$Car))
      #, Lig = c(age0p2$Lig, tmp$Lig))
    }
    
    # Replace all the drones
    print("Replace Drones, P2")
    print(Sys.time())
    
    age1$Mel <- replaceDrones(age1$Mel)
    age1$Car <- replaceDrones(age1$Car)
    #age1$Lig <- replaceDrones(age1$Lig)
    if (year > 1) {
      age2$Mel <- replaceDrones(age2$Mel)
      age2$Car <- replaceDrones(age2$Car)
      #age2$Lig <- replaceDrones(age2$Lig)
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
      #DCALig <- createDCA(age1$Lig)
      #age0p2$Lig <- cross(age0p2$Lig, drones = pullDroneGroupsFromDCA(DCA = DCALig, n = nColonies(age0p2$Lig), nDrones = nFathersPoisson))
    } else {
      DCAMel <- createDCA(c(age1$Mel, age2$Mel))
      fathersMel <- pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson)
      fathersMel[[1]] <- c(fathersMel[[1]], createDrones(age1$Mel[[1]], nInd = 2))
      age0p2$Mel <- cross(age0p2$Mel, drones = fathersMel)
      DCACar <- createDCA(c(age1$Car, age2$Car))
      fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
      fathersCar[[1]] <- c(fathersCar[[1]], createDrones(age1$Car[[1]], nInd = 2))
      age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
      
    }
    
    # Collapse
    age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p2collapse),
                 Car = selectColonies(age1$Car, p = 1 - p2collapse))
    #,Lig = selectColonies(age1$Lig, p = 1 - p2collapse))
    if (year > 1) {
      age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
                   Car = selectColonies(age2$Car, p = 1 - p2collapse))
      #,Lig = selectColonies(age2$Lig, p = 1 - p2collapse))
    }
    
    # Merge all age 0 colonies (from both periods)
    age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
                 Car = c(age0p1$Car, age0p2$Car))
    #,Lig = c(age0p1$Lig, age0p2$Lig))
    
   
    # Period3 ------------------------------------------------------------------
    # Collapse age0 queens
    print("PERIOD 3")
    print("Collapse colonies, P3")
    print(Sys.time())
    
    #selection on Fitness queen gv
    
    #age0
    #Mellifera
    
      MelBritFit<-sapply(getPheno(age0$Mel, caste = "queen"), function(x) x[1,2])#select on Mel Brit fitness
      queensID<-names(sort(MelBritFit,T))
      Nselectcolon<-round(length(queensID)*(1-p3collapseAge0)) #calculate how many colonies will collapse
      age0MelqueensID<-queensID[1:Nselectcolon] #select the queens ids that will not collapse
      
      MelBritFit1<-sapply(getPheno(age1$Mel, caste = "queen"), function(x) x[1,2])#select on Mel Brit fitness
      queensID1<-names(sort(MelBritFit1,T))
      Nselectcolon1<-round(length(queensID1)*(1-p3collapseAge1)) #calculate how many colonies will collapse
      age1MelqueensID<-queensID1[1:Nselectcolon1] #select the queens ids that will not collapse
    
    #Collapse 
    age0 <- list(Mel = selectColonies(age0$Mel, ID=age0MelqueensID),
                 Car = selectColonies(age0$Car, p = 1 - p3collapseAge0))

    age1 <- list(Mel = selectColonies(age1$Mel, ID=age1MelqueensID),
                 Car = selectColonies(age1$Car, p = 1 - p3collapseAge1))
 
    age2 <- list(Mel = NULL, Car = NULL) #We don't need this but just to show the workflow!!!
    
    # Maintain the number of colonies ------------------------------------------
    # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
    print("Maintain the number, P2")
    print(Sys.time())
    
    age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
    age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)
    #age0$Lig <- maintainCarSize(age0 = age0$Lig, age1 = age1$Lig)
    
    for (subspecies in c("Mel", "Car")) {     #,"Lig"
      if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == IrelandSize
          | (nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == CarSize)
      {
      } 
      else 
      {stop(paste0("The number of colonies for ", subspecies, " does not match the population size!"))}
    }
    
    
    #track mean IBD and variance of IBD
    
    #this creates two empty dataframes where all the information of each year and each replica will be recorded
    if (year==1){ 
      columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","HoneyYieldEu","FitnessBrit","FitnessEu","Homocigosity")
      MeanVarMel <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for mellifera
      colnames(MeanVarMel)<-columnheaders
      MeanVarCar <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for carnica
      colnames(MeanVarCar)<-columnheaders
      
    } 
    
    if (year==1 & Rep==1){
      uno<-0
    }else{
      uno<-nrow(colonyRecords) #this is to show where to start recording values on the dataframe
    }
    
    
    #record values of Mellifera population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel", Rep=Rep)
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Mel, year = year, population = "Mel", Rep=Rep)
    dos<-nrow(colonyRecords)#this is to show where to end recording values on the dataframe
  
    #create a dataframe with the mellifera mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
    newrow1<- data.frame(MeanIBD=mean(colonyRecords[(uno+1):dos,"IBD"]), VarIBD=var(colonyRecords[(uno+1):dos,"IBD"]), Year=year,Rep=Rep,Population="Mel"
                         ,HoneyYieldBrit=mean(colonyRecords[(uno+1):dos,"gvQueens_BritHY"]), FitnessBrit=mean(colonyRecords[(uno+1):dos,"gvQueens_BritFit"])
                         ,HoneyYieldEu=mean(colonyRecords[(uno+1):dos,"gvQueens_EuHY"]), FitnessEu=mean(colonyRecords[(uno+1):dos,"gvQueens_EuFit"])
                         ,Homocigosity=mean(colonyRecords[(uno+1):dos,"pHomBrood"]))
   
    #record values of Carnica population
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car", Rep=Rep)
    colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car", Rep=Rep)
    tres<-nrow(colonyRecords)
    
    #create a dataframe with the carnica mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
    newrow2<- data.frame(MeanIBD=mean(colonyRecords[(dos+1):tres,"IBD"]), VarIBD=var(colonyRecords[(dos+1):tres,"IBD"]), Year=year,Rep=Rep, Population="Car"
                         ,HoneyYieldEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]), FitnessEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuFit"])
                         ,HoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]), FitnessBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritFit"])
                         ,Homocigosity=mean(colonyRecords[(dos+1):tres,"pHomBrood"]))
    
    #Combine what we had in each dataframe with the new info, so each year the dataframe updates with new values
    MeanVarMel<-rbind(MeanVarMel,newrow1)
    MeanVarCar<-rbind(MeanVarCar,newrow2)
    
  } # Year-loop
  
  a <- toc()
  loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))
  print("Saving image data")
  save.image(paste0("Burnin",Rep,".RData"))
remove(MeanVarMel)
remove(MeanVarCar)
} # Rep-loop






Rep<-1
nRep<-5
for (Rep in 1:nRep) {
  cat(paste0("Rep: ", Rep, "/", nRep, "\n"))
  tic(paste0(nYear, 'y loop'))         # Measure cpu time
  Rprof()
 
#dir.create(paste0("C:/Users/Usuario/Honeybeeimports/Replica",Rep,"/ScenarioRep",Rep))
#load(paste0("C:/Users/Usuario/Honeybeeimports/Replica",Rep,"/Burnin",Rep,".RData")) #load the burn in
#setwd(paste0("C:/Users/Usuario/Honeybeeimports/Replica",Rep,"/ScenarioRep",Rep))


  
  dir.create(paste0("C:/Users/Usuario/Desktop/Máster/Bees/Testeando/Replica",Rep,"/ScenarioRep",Rep))
  load(paste0("C:/Users/Usuario/Desktop/Máster/Bees/Testeando/Replica",Rep,"/Burnin",Rep,".RData")) #load the burn in
  setwd(paste0("C:/Users/Usuario/Desktop/Máster/Bees/Testeando/Replica",Rep,"/ScenarioRep",Rep))
  
  
year=11
nYear=30
for (year in 11:nYear) {
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
  
  # Period1 ------------------------------------------------------------------
  # Build-up the colonies
  print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
  print(Sys.time())
  age1 <- list(Mel = buildUp(age1$Mel),
               Car = buildUp(age1$Car))
  #,Lig = buildUp(age1$Lig))
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
  age0p1 <- list(Mel = tmp$Mel$split, Car = tmp$Car$split) #,Lig = tmp$Lig$split
  
  if (year > 1) {
    # Split all age2 colonies
    tmp <- list(Mel = split(age2$Mel),
                Car = split(age2$Car))
    #,Lig = split(age2$Lig))
    age2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant
    ) #,Lig = tmp$Lig$remnant
    # The queens of the splits are 0 years old
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                   Car = c(age0p1$Car, tmp$Car$split))
    #,Lig = c(age0p1$Lig, tmp$Lig$split))
  }
  
  
  # Create virgin queens
  # Sample colony for the virgin queens
  print("Create virgin queens, period 1")
  print(Sys.time())
  virginDonor <- list(Mel = sample.int(n = nColonies(age1$Mel), size = 1),
                      Car = sample.int(n = nColonies(age1$Car), size = 1))
  #Lig = sample.int(n = nColonies(age1$Lig), size = 1))
  # Virgin queens for splits!
  
  if (year<11){
    pImport<-0
  }else {
    pImport<-0.3
  }
  
  tmp <- (Mel = pullColonies(age0p1$Mel, p=pImport)) #pull colonies to requeen with imports
  IdImportColonies<-getId(tmp$pulled) #get the ids of the imported colonies
  age0p1 <- list(Mel = tmp$remnant,
                 MelImport = tmp$pulled,
                 Car = c(age0p1$Car, tmp$Car$split))
  #,Lig = c(age0p1$Lig, tmp$Lig$split))
  virginQueens <- list(Mel = createVirginQueens(age1$Mel[[virginDonor$Mel]], nInd = nColonies(age0p1$Mel)),
                       Car = createVirginQueens(age1$Car[[virginDonor$Car]], nInd = nColonies(age0p1$Car)+nColonies(age0p1$MelImport)))
  #,Lig = createVirginQueens(age1$Lig[[virginDonor$Lig]], nInd = nColonies(age0p1$Lig)+(nColonies(age0p1$MelImport)/2)))
  
  #requeen with carnica and mellifera for Mel and carnica for Car
  nColoniesMelImport<-nColonies(age0p1$MelImport)
  nColoniesCar<-nColonies(age0p1$Car)+nColonies(age0p1$MelImport)
  #nColoniesLig<-nColonies(age0p1$Lig)+(nColonies(age0p1$MelImport)/2)
  age0p1 <- list(Mel = c(reQueen(age0p1$Mel, queen = (virginQueens$Mel)) ,
                         reQueen(age0p1$MelImport, queen = c((virginQueens$Car)[1:nColoniesMelImport]))),       #,(virginQueens$Lig)[1:(nColoniesMelImport/2)]))),
                 Car = reQueen(age0p1$Car, queen = virginQueens$Car[(nColoniesMelImport+1):nColoniesCar]))
  #,Lig = reQueen(age0p1$Lig, queen = virginQueens$Lig[((nColoniesMelImport/2)+1):nColoniesLig]))
  
 
  # Swarm a percentage of age1 colonies
  print("Swarm colonies, P1")
  print(Sys.time())
  tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
              Car = pullColonies(age1$Car, p = p1swarm))
  #,Lig = pullColonies(age1$Lig, p = p1swarm))
  age1 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  #,Lig = tmp$Lig$remnant)
  tmp <- list(Mel = swarm(tmp$Mel$pulled),
              Car = swarm(tmp$Car$pulled))
  #,Lig = swarm(tmp$Lig$pulled))
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                 Car = c(age0p1$Car, tmp$Car$remnant))
  #,Lig = c(age0p1$Lig, tmp$Lig$remnant))
  age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
               Car = c(age1$Car, tmp$Car$swarm))
  #,Lig = c(age1$Lig, tmp$Lig$swarm))
  
  
  if (year > 1) {
    # Swarm a percentage of age2 colonies
    tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
                Car = pullColonies(age2$Car, p = p1swarm))
    #, Lig = pullColonies(age2$Lig, p = p1swarm))
    age2 <- list(Mel = tmp$Mel$remnant,     
                 Car = tmp$Car$remnant)
    #,Lig = tmp$Lig$remnant)
    tmp <- list(Mel = swarm(tmp$Mel$pulled),
                Car = swarm(tmp$Car$pulled))
    #,Lig = swarm(tmp$Lig$pulled))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                   Car = c(age0p1$Car, tmp$Car$remnant))
    #,Lig = c(age0p1$Lig, tmp$Lig$remnant))
    age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
                 Car = c(age2$Car, tmp$Car$swarm))
    #,Lig = c(age2$Lig, tmp$Lig$swarm))
  }
  
  # Supersede age1 colonies
  print("Supersede colonies, P1")
  print(Sys.time())
  tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
              Car = pullColonies(age1$Car, p = p1supersede))
  #,Lig = pullColonies(age1$Lig, p = p1supersede))
  age1 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  #,Lig = tmp$Lig$remnant)
  tmp <- list(Mel = supersede(tmp$Mel$pulled),
              Car = supersede(tmp$Car$pulled))
  #,Lig = supersede(tmp$Lig$pulled))
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                 Car = c(age0p1$Car, tmp$Car))
  #, Lig = c(age0p1$Lig, tmp$Lig))
  
  if (year > 1) {
    # Supersede age2 colonies
    tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
                Car = pullColonies(age2$Car, p = p1supersede))
    #,Lig = pullColonies(age2$Lig, p = p1supersede))
    age2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
    #,Lig = tmp$Lig$remnant)
    tmp <- list(Mel = supersede(tmp$Mel$pulled),
                Car = supersede(tmp$Car$pulled))
    #,Lig = supersede(tmp$Lig$pulled))
    age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
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
    
  } else {
    DCAMel <- createDCA(c(age1$Mel, age2$Mel))
    age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
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
  print("Mate colonies, P2")
  print(Sys.time())
  
  
  if (year == 1) {
    DCAMel <- createDCA(age1$Mel)
    age0p2$Mel <- cross(age0p2$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson))
    DCACar <- createDCA(age1$Car)
    age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
    
  } else {
    DCAMel <- createDCA(c(age1$Mel, age2$Mel))
    fathersMel <- pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson)
    fathersMel[[1]] <- c(fathersMel[[1]], createDrones(age1$Mel[[1]], nInd = 2))
    age0p2$Mel <- cross(age0p2$Mel, drones = fathersMel)
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
 
  
  
  # Period3 ------------------------------------------------------------------
  # Collapse age0 queens
  print("PERIOD 3")
  print("Collapse colonies, P3")
  print(Sys.time())
  
  #selection on Fitness queen gv
  
  #age0
  #Mellifera
  age0$Mel
  coco<-selectColonies(age0$Mel,ID= IdImportColonies)
  
  MelBritFit<-sapply(getPheno(age0$Mel, caste = "queen"), function(x) x[1,2])#select on Mel Brit fitness
  queensID<-names(sort(MelBritFit,T))
  Nselectcolon<-round(length(queensID)*(1-p3collapseAge0)) #calculate how many colonies will collapse
  age0MelqueensID<-queensID[1:Nselectcolon] #select the queens ids that will not collapse
  
  MelBritFit1<-sapply(getPheno(age1$Mel, caste = "queen"), function(x) x[1,2])#select on Mel Brit fitness
  queensID1<-names(sort(MelBritFit1,T))
  Nselectcolon1<-round(length(queensID1)*(1-p3collapseAge1)) #calculate how many colonies will collapse
  age1MelqueensID<-queensID1[1:Nselectcolon1] #select the queens ids that will not collapse
  
  #Collapse 
  age0 <- list(Mel = selectColonies(age0$Mel, ID=age0MelqueensID),
               Car = selectColonies(age0$Car, p = 1 - p3collapseAge0))
  
  age1 <- list(Mel = selectColonies(age1$Mel, ID=age1MelqueensID),
               Car = selectColonies(age1$Car, p = 1 - p3collapseAge1))
  
  age2 <- list(Mel = NULL, Car = NULL) #We don't need this but just to show the workflow!!!
  
  # Maintain the number of colonies ------------------------------------------
  # Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
  print("Maintain the number, P2")
  print(Sys.time())
  
  age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
  age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)
  #age0$Lig <- maintainCarSize(age0 = age0$Lig, age1 = age1$Lig)
  
  for (subspecies in c("Mel", "Car")) {     #,"Lig"
    if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == IrelandSize
        | (nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == CarSize)
    {
    } 
    else 
    {stop(paste0("The number of colonies for ", subspecies, " does not match the population size!"))}
  }
  
  #track mean IBD and variance of IBD
  
  #this creates two empty dataframes where all the information of each year and each replica will be recorded
  if (year==1 & Rep==1){ 
    columnheaders<-c("MeanIBD","VarIBD","Year","Rep","Population","HoneyYieldBrit","HoneyYieldEu","FitnessBrit","FitnessEu","Homocigosity")
    MeanVarMel <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for mellifera
    colnames(MeanVarMel)<-columnheaders
    MeanVarCar <- data.frame(matrix(ncol = length(columnheaders), nrow = 0)) #dataframe for carnica
    colnames(MeanVarCar)<-columnheaders
    uno<-0
  } else{
    uno<-nrow(colonyRecords) #this is to show where to start recording values on the dataframe
  }
  
  
  #record values of Mellifera population
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel", Rep=Rep)
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Mel, year = year, population = "Mel", Rep=Rep)
  dos<-nrow(colonyRecords)#this is to show where to end recording values on the dataframe
  
  #create a dataframe with the mellifera mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
  newrow1<- data.frame(MeanIBD=mean(colonyRecords[(uno+1):dos,"IBD"]), VarIBD=var(colonyRecords[(uno+1):dos,"IBD"]), Year=year,Rep=Rep,Population="Mel"
                       ,HoneyYieldBrit=mean(colonyRecords[(uno+1):dos,"gvQueens_BritHY"]), FitnessBrit=mean(colonyRecords[(uno+1):dos,"gvQueens_BritFit"])
                       ,HoneyYieldEu=mean(colonyRecords[(uno+1):dos,"gvQueens_EuHY"]), FitnessEu=mean(colonyRecords[(uno+1):dos,"gvQueens_EuFit"])
                       ,Homocigosity=mean(colonyRecords[(uno+1):dos,"pHomBrood"]))
  
  #record values of Carnica population
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car", Rep=Rep)
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car", Rep=Rep)
  tres<-nrow(colonyRecords)
  
  #create a dataframe with the carnica mean IBD, variance of IBD, mean Honey yield, mean fitness and mean homocigosity
  newrow2<- data.frame(MeanIBD=mean(colonyRecords[(dos+1):tres,"IBD"]), VarIBD=var(colonyRecords[(dos+1):tres,"IBD"]), Year=year,Rep=Rep, Population="Car"
                       ,HoneyYieldEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuHY"]), FitnessEu=mean(colonyRecords[(dos+1):tres,"gvQueens_EuFit"])
                       ,HoneyYieldBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritHY"]), FitnessBrit=mean(colonyRecords[(dos+1):tres,"gvQueens_BritFit"])
                       ,Homocigosity=mean(colonyRecords[(dos+1):tres,"pHomBrood"]))
  
  #Combine what we had in each dataframe with the new info, so each year the dataframe updates with new values
  MeanVarMel<-rbind(MeanVarMel,newrow1)
  MeanVarCar<-rbind(MeanVarCar,newrow2)
 } #end of year loop
save.image(paste0("Selection",Rep,".RData"))
CombinedDf<-rbind(MeanVarCar,MeanVarMel)
write.csv(CombinedDf, paste0("Repi",Rep,".csv"))
remove(MeanVarCar)
remove(MeanVarMel)

}#end of rep loop


Rep<-5

setwd(paste0("C:/Users/Usuario/Honeybeeimports/Replica",Rep,".RData/ScenarioRep",Rep))
load((paste0("Selection",Rep,".RData")))
Rep1M<-MeanVarMel
Rep1C<-MeanVarCar
Rep2M<-MeanVarMel
Rep2C<-MeanVarCar
Rep3M<-MeanVarMel
Rep3C<-MeanVarCar
Rep4M<-MeanVarMel
Rep4C<-MeanVarCar
Rep5M<-MeanVarMel
Rep5C<-MeanVarCar
remove(MeanVarMel)
remove(MeanVarCar)
Melif<-rbind(Rep1M,Rep2M,Rep3M,Rep4M,Rep5M)
Carnic<-rbind(Rep1C,Rep2C,Rep3C,Rep4C,Rep5C)
#Dataframes with mean replicates

Mother<-read.csv("mother.csv")

#This calculates the mean and sd through different replicates of the mean IBD, mean Honey Yield, fitness and homocigosity
CombinedDf<-CombinedDf %>% group_by(Population)
Carni<-patata2[1:20,]
Meli<-patata2[21:40,]
groupedCar<-Carnic %>% group_by(Year) %>% summarise(meanIBD = mean(MeanIBD),
                                                        sdIBD = sd(MeanIBD),
                                                    quanIBDl= quantile(MeanIBD,p=0.025),
                                                    quanIBDh= quantile(MeanIBD,p=0.975),
                                                        meanHoneyYield = mean(HoneyYieldEu),
                                                        sdHoneyYield = sd(HoneyYieldEu),
                                                    quanHYl= quantile(HoneyYieldEu,p=0.025),
                                                    quanHYh= quantile(HoneyYieldEu,p=0.975),
                                                        meanFitness = mean(FitnessEu),
                                                        sdFitness = sd(FitnessEu),
                                                    quanFl= quantile(FitnessEu,p=0.025),
                                                    quanFh= quantile(FitnessEu,p=0.975),
                                                        meanHomocigosity = mean(Homocigosity),
                                                        sdHomocigosity = sd(Homocigosity),
                                                    quanHol= quantile(Homocigosity,p=0.025),
                                                    quanHoh= quantile(Homocigosity,p=0.975),
                                                        Pop= "Car"
)

groupedMel<-Melif  %>% group_by(Year) %>% summarise(meanIBD = mean(MeanIBD),
                                                     sdIBD = sd(MeanIBD),
                                                    quanIBDl= quantile(MeanIBD,p=0.025),
                                                    quanIBDh= quantile(MeanIBD,p=0.975),
                                                     meanHoneyYield = mean(HoneyYieldBrit),
                                                     sdHoneyYield = sd(HoneyYieldBrit),
                                                    quanHYl= quantile(HoneyYieldBrit,p=0.025),
                                                    quanHYh= quantile(HoneyYieldBrit,p=0.975),
                                                     meanFitness = mean(FitnessBrit),
                                                     sdFitness = sd(FitnessBrit),
                                                    quanFl= quantile(FitnessBrit,p=0.025),
                                                    quanFh= quantile(FitnessBrit,p=0.975),
                                                     meanHomocigosity = mean(Homocigosity),
                                                     sdHomocigosity = sd(Homocigosity),
                                                    quanHol= quantile(Homocigosity,p=0.025),
                                                    quanHoh= quantile(Homocigosity,p=0.975),
                                                     Pop= "Mel"
                                                    
)


#Plot the mean of the replicas

#plot for IBD
grouped<-rbind(groupedMel,groupedCar)

a<-ggplot(data = grouped, aes(x = Year, group=Pop)) + geom_line(aes(y= meanIBD, colour=Pop)) +
  geom_ribbon(aes(ymin = quanIBDl, ymax = quanIBDh), alpha = 0.1)+
  geom_vline(xintercept= 10, linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))

#plot for Fitness
#(ymin = meanFitness - sdFitness, ymax = meanFitness + sdFitness)
b<-ggplot(data = grouped, aes(x = Year, y = meanFitness, group=Pop)) + geom_line(aes(colour=Pop)) +
  geom_ribbon(aes(ymin = quanFl, ymax = quanFh), alpha = 0.1) +
  geom_vline(xintercept= 10, linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))

#plot for HoneyYield
c<-ggplot(data = grouped, aes(x = Year, y = meanHoneyYield, group=Pop)) + geom_line(aes(colour=Pop)) +
  geom_ribbon(aes(ymin = quanHYl, ymax = quanHYh), alpha = 0.1)+
  geom_vline(xintercept= 10, linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) 
#plot for Homocigosity
d<-ggplot(data = grouped, aes(x = Year, y = meanHomocigosity, group=Pop)) + geom_line(aes(colour=Pop)) +
  geom_ribbon(aes(ymin = quanHol, ymax = quanHoh), alpha = 0.1)+
  geom_vline(xintercept= 10, linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
ggarrange(a, b, c, d, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2, common.legend = T)





#plot of mellifera traits

groupedMel2<-Melif %>% group_by(Year) %>% summarise(
                                                        meanHoneyYieldBrit = mean(HoneyYieldBrit),
                                                        sdHoneyYieldBrit = sd(HoneyYieldBrit),
                                                        meanFitnessBrit = mean(FitnessBrit),
                                                        sdFitnessBrit = sd(FitnessBrit),
                                                        meanHoneyYieldEu = mean(HoneyYieldEu),
                                                        sdHoneyYieldEu = sd(HoneyYieldEu),
                                                        meanFitnessEu = mean(FitnessEu),
                                                        sdFitnessEu = sd(FitnessBrit),
                                                        Pop= "Mel"
)

hola<-ggplot(data=groupedMel2, aes(x=Year)) + geom_line(aes(y=meanHoneyYieldBrit),colour="orange") +
  geom_ribbon(aes(ymin = meanHoneyYieldBrit - sdHoneyYieldBrit, ymax = meanHoneyYieldBrit + sdHoneyYieldBrit), alpha = 0.1)+
 geom_line(aes(y=meanFitnessBrit), colour="#484AF5") +
  geom_ribbon(aes(ymin = meanFitnessBrit - sdFitnessBrit, ymax = meanFitnessBrit + sdFitnessBrit), alpha = 0.1)+
 geom_line(aes(y=meanHoneyYieldEu), colour="red") +
  geom_ribbon(aes(ymin = meanHoneyYieldEu - sdHoneyYieldEu, ymax = meanHoneyYieldEu + sdHoneyYieldEu), alpha = 0.1)+
 geom_line(aes(y=meanFitnessEu),colour="#1CABF7") +
  geom_ribbon(aes(ymin = meanFitnessEu - sdFitnessEu, ymax = meanFitnessEu + sdFitnessEu), alpha = 0.1)+
  geom_vline(xintercept= 10, linetype="dotted")+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  labs(y="traits")

library(tidyr)
df <- groupedMel2 %>% pivot_longer( cols= c("meanHoneyYieldBrit","meanFitnessBrit", "meanHoneyYieldEu","meanFitnessEu"), names_to="Meanval",values_to= "Mean")
ggplot(data=df, aes(x=Year, y=Mean)) + geom_line(aes(colour=Meanval))+
  geom_ribbon(aes(group=Sd, ymin = Mean- Sd, ymax = Mean + Sd), alpha = 0.1)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  labs(y="traits")
df


#plot of carnica traits
groupedCar2<-Carnic %>% group_by(Year) %>% summarise(
  meanHoneyYieldBrit = mean(HoneyYieldBrit),
  sdHoneyYieldBrit = sd(HoneyYieldBrit),
  meanFitnessBrit = mean(FitnessBrit),
  sdFitnessBrit = sd(FitnessBrit),
  meanHoneyYieldEu = mean(HoneyYieldEu),
  sdHoneyYieldEu = sd(HoneyYieldEu),
  meanFitnessEu = mean(FitnessEu),
  sdFitnessEu = sd(FitnessBrit),
  Pop= "Car")

adiso<-ggplot(data=groupedCar2, aes(x=Year)) + geom_line(aes(y=meanHoneyYieldBrit),colour="orange") +
  geom_ribbon(aes(ymin = meanHoneyYieldBrit - sdHoneyYieldBrit, ymax = meanHoneyYieldBrit + sdHoneyYieldBrit), alpha = 0.1)+
  geom_line(aes(y=meanFitnessBrit), colour="#484AF5") +
  geom_ribbon(aes(ymin = meanFitnessBrit - sdFitnessBrit, ymax = meanFitnessBrit + sdFitnessBrit), alpha = 0.1)+
  geom_line(aes(y=meanHoneyYieldEu), colour="red") +
  geom_ribbon(aes(ymin = meanHoneyYieldEu - sdHoneyYieldEu, ymax = meanHoneyYieldEu + sdHoneyYieldEu), alpha = 0.1)+
  geom_line(aes(y=meanFitnessEu),colour="#1CABF7") +
  geom_ribbon(aes(ymin = meanFitnessEu - sdFitnessEu, ymax = meanFitnessEu + sdFitnessEu), alpha = 0.1)+
  geom_vline(xintercept= 10, linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  labs(y="traits")



ggarrange(hola, adios, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)




























#Plot of each of the replicas

df <- bind_rows(
  MeanVarMel %>% mutate(Population = "Mel"),
  MeanVarCar %>% mutate(Population="Car"))
df
#plot of Mean IBD
ggplot(df, aes(x=Year, y=MeanIBD, group=Population)) + 
  geom_smooth(aes(colour=Population),se=F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))

#plot Fitness for different reps
ggplot(df, aes(x=Year, y=Fitness, group=Population)) + 
  geom_smooth(aes(colour=Population),se=F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))


#plot Honey Yield
ggplot(df, aes(x=Year, y=HoneyYield, group=Population)) + 
  geom_smooth(aes(colour=Population),se=F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))

#Plot Homocigosity 
ggplot(df, aes(x=Year, y=Homocigosity, group=Population)) + 
  geom_smooth(aes(colour=Population),se=F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))


Rep<-5
load(paste0("C:/Users/Usuario/Honeybeeimports/Replica",Rep,".RData/Burnin",Rep,".RData"))

MeanVarMel[1,1]<-1
MeanVarMel[1,2]<-0
Rep1<-rbind(MeanVarMel,MeanVarCar)
Rep2<-rbind(MeanVarMel,MeanVarCar)
Rep3<-rbind(MeanVarMel,MeanVarCar)
Rep4<-rbind(MeanVarMel,MeanVarCar)
Rep5<-rbind(MeanVarMel,MeanVarCar)
remove(MeanVarCar)
remove(MeanVarMel)
motherburnin<-rbind(Rep1,Rep2,Rep3,Rep4,Rep5)
motherburnin

load("C:/Users/Usuario/Honeybeeimports/5-burnin-5-selection-run-by-my-computer.RData")