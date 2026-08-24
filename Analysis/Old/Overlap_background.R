set.seed(123)  # Set seed for reproducibility
`%nin%` = Negate('%in%')
library(tidyverse)
library(viridis)
library(SIBER)
library(factoextra)
library(lme4)
library(future.apply) ## Package for parallelizing lapply
library(progressr)
system.time({
## Overlap calculation is computationally intensive across the many post draws and species pairs
### Parallelize the process to help it run faster
#### Note - 1000 draws appears to take 24+ hrs on laptop
plan(multisession, workers = parallel::detectCores() - 1) ## # cores computer has subtract cores you want to keep free
## Create a little progress bar in the consol so you know how it's working
handlers(global = TRUE)
handlers("txtprogressbar")

with_progress({
  
  p = progressor(along = overlap_subset) ## progress bar
  
  ## List with the overlap data to compress after loop
  overlap_list = future_lapply(overlap_subset, function(h) {
    
    p()  # update progress bar
    
    overlap(
      comm = h,
      dr = 1000,
      posterior = siber.setup[[2]],
      pair_df = pair_df
      
    )
  })
})
})
save(overlap_list, file = "Data/RData/overlap.list.RData")

