#!/usr/bin/env Rscript

# Gather reps
reps <- dir()
reps <- reps[grepl(x = reps, pattern = "rep_")]

results <- NULL
for (rep in reps) {
  setwd(rep)

  # Gather scenarios
  scenarios <- dir()
  scenarios <- scenarios[grepl(x = scenarios, pattern = "scenario_")]
  for (scenario in scenarios) {
    resultsFilename <- paste0(scenario, "/scenarioData.csv")
    if (file.exists(resultsFilename)) {
      result <- read.csv(file = resultsFilename)
      result$scenario <- scenario
      is (is.null(results)) {
        results <- result
      } else {
        results <- rbind(results, result)
      }
    }
  }

  setwd("..")
}

# Parse parameters/factors
# TODO: strsplit(results$scenario)

# Analyse results
# TODO: whatever plotting or summarisation we need
