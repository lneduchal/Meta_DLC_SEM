
library(glue)

options(scipen = 999)
rm(list = ls())
graphics.off()

set.seed(24)
setwd("~/Meta_DLC_SEM")

# Information about the run
source("configuration.R")

invisible(
  
  {
    
    cat("Model evaluation:", "\n")
    cat("----------------------\n\n")
        
    time_start <- Sys.time()
    cat(glue("Started running {candidate_model} at {time_start}.\n"))
    
    # Running the candidate model with the generated dataset as targets
    source(glue("candidate_models/model_evaluation/evaluation.R"))
    
    time_end <- Sys.time()
    cat("\n")
    cat(glue("Finished running at {time_end}."), "\n")
      
    cat("\n----------------------\n")
    cat("Evaluation successful.")
    
  }
  
)
