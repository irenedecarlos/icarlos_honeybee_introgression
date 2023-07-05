#!/usr/bin/env Rscript
 
args <- commandArgs(trailingOnly = TRUE)
 
# Argument 1 is number of replicates we want to do
nRep <- args[1]
 
 
# ---- Specify factors and their levels and combinations -----------------------
#nRep<-1
reps <- 1:nRep
 
combinations <- data.frame(Rep = NA,
                            setting = c("burnin", "scenario", "scenario"),
                            selection = c(NA, 0.3, 0.2))
 
for (Rep in 1:nRep) {
  tmp <- combinations
   tmp$Rep <- Rep
   if (Rep == 1) {
     combinationsAll <- tmp
   } else {
     combinationsAll <- rbind(combinationsAll, tmp)
   }
}
 
 # ---- Run combinations --------------------------------------------------------
 
for (combination in 1:nrow(combinationsAll)) {
   # combination <- 1
   # ... write the script
 
 Rep <- combinationsAll[combination, "Rep"]
 setting <- combinationsAll[combination, "setting"]
 selection <- combinationsAll[combination, "selection"]
 scriptFilename <- paste0("script_Rep_", Rep,
                          "_setting_", setting,
                          "_selection_", selection,
                          ".sh")

 cat("Preparing & running ", scriptFilename, "\n")

 sink(file = scriptFilename)
 cat("#!/bin/sh\n")
 cat("#\n")
 cat("# ---- Grid Engine options (lines prefixed with #$) ----\n")
 cat("#\n")
 cat("#$ -P roslin_HighlanderLab\n")
 cat("#$ -N HBeeImport_", Rep, "_", setting, "_", selection, "\n", sep = "")
 cat("#$ -cwd\n")
 cat("#$ -l h_vmem=40G\n")
 cat("#$ -pe sharedmem 1\n")
 cat("#$ -l h_rt=20:00:00\n")
 cat("#$ -j yes\n")
 cat("#$ -o HBeeImport_", Rep, "_", setting, "_", selection, ".out\n", sep = "")
 cat("#\n")
 cat("# ---- The commands ----\n")
 cat("#\n")
 cat("# Initialise the environment modules\n")
 cat(". /etc/profile.d/modules.sh\n")
 cat("module load roslin/R/4.1.0\n")
 cat("#\n")
 cat("# Run the program\n")
 cat("./simulation.R", Rep, setting, selection,"\n")
 cat("#\n")
 sink()

 # ... launch the script
 if (setting=="burnin"){
   system(command = paste0("qsub ", scriptFilename))
 } else {
   system(command = paste0("qsub -hold_jid HBeeImport_", Rep,"_burnin_NA ", scriptFilename))
   
 }
 
}

 
 

