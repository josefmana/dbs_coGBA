
rm( list = ls() ) # clear environment
s = 87542 # seed for reproducibility
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

library(here) # directory management
library(tidyverse) # data wrangling
library(brms) # model fitting
library(tidybayes) # posterior manipulation
library(gt) # tables

# read outcome data
d0 <-
  
  lapply(
    
    setNames( c("motor","psych"), c("motor","psych") ),
    function(i)
      read.csv( here( "_raw", paste0(i,"_longit.csv") ), sep = "," ) %>%
      mutate( gba = ifelse( id %in% c( read.table( here("_raw","gba.txt"), header = F ) )$V1, 1, 0 ), .after = sex )
  
  )
