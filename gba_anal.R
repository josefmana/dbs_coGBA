
rm( list = ls() ) # clear environment
s = 87542 # seed for reproducibility
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

library(here) # directory management
library(tidyverse) # data wrangling
library(brms) # model fitting
library(tidybayes) # posterior manipulation
library(patchwork) # plot grids
#library(gt) # tables
theme_set( theme_minimal() )

# read outcome data
d0 <-
  
  lapply(
    
    setNames( c("motor","psych"), c("motor","psych") ),
    function(i)
      read.csv( here( "_raw", paste0(i,"_longit.csv") ), sep = "," ) %>%
      mutate( gba = ifelse( id %in% c( read.table( here("_raw","gba.txt"), header = F ) )$V1, 1, 0 ), .after = sex )
  
  )


# REGRESSION MODELS OF DRS-2 ----

# prepare folder for models
if (!dir.exists("mods") ) dir.create("mods")

# compute shifted time from surgery
md_time <- median( d0$psych[ d0$psych$event == "screening" & complete.cases(d0$psych$drsii), "stimtime_years"], na.rm = T )

# prepare data for analysis
d1 <-
  d0$psych %>%
  filter( complete.cases(drsii) ) %>% # keep rows with DRS-2 only
  mutate( time = stimtime_years - md_time ) %>% # shift time such that median pre-surgery assessment is t0
  mutate( GBA = as.character(gba) ) %>% # for brms to not read GBA as a number
  select( id, GBA, time, drsii, faq, pdaq, bdi, staix1, staix2 ) # select variables of interest

# set-up model formulas (i.e., linear models)
f <-
  
  list(
    replicate = bf( drsii ~ 1 + GBA * time + (1 + time | id) ) + gaussian(), # replication of previous studies
    heterosce = bf( drsii ~ 1 + GBA * time + (1 + time | id), sigma ~ 1 + time ) + gaussian(), # adding time-varying heteroscedasticity,
    betabinom = bf( drsii | trials(144) ~ 1 + GBA * time + (1 + time | id), phi ~ 1 + time ) + beta_binomial() # beta-binomial with time varying overdispersion
  )

# prepare names for them
n <-
  
  list(
    replicate = "mu ~ 1 + GBA * Time + (1 + Time | ID),\nsigma ~ 1, Gaussian",
    heterosce = "mu ~ 1 + GBA * Time + (1 + Time | ID),\nsigma ~ 1 + Time, Gaussian",
    betabinom = "mu ~ 1 + GBA * Time + (1 + Time | ID),\nphi ~ 1 + Time, Beta-Binomial"
  )

# use default brms priors for this exploration
p <- NULL

# fit them
m <-
  
  lapply(
    
    setNames( names(f), names(f) ),
    function(i)
      brm( formula = f[[i]], prior = p, data = d1, seed = s, file = here( "mods", paste0(i,".rds") ) )
    
  )

# POSTERIOR PREDICTIVE CHECKS ----

# prepare folder for figures
if (!dir.exists("figs") ) dir.create("figs")

## statistics prediction ----

ppc_stat <-
  
  # loop through models
  lapply(
    
    setNames( names(m), names(m) ),
    function(i)
      
      # loop through stats
      lapply(
        
        setNames( c("mean","sd","median","IQR"), c("mean","sd","median","IQR") ),
        function(j)
          
          pp_check( m[[i]], type = "stat", stat = j ) +
          labs( x = j ) +
          theme( legend.position = "none" )
        
      )
    
  )

# add 2D mean/SD posterior prediction
for ( i in names(ppc_stat) ) ppc_stat[[i]]$`2d` <- pp_check( m[[i]], type = "stat_2d", stat = c("mean","sd") ) + labs( title = n[[i]] ) + theme( legend.position = "none", plot.title = element_text( hjust = 0.5 ) )

# pull mean/SD figures together
with(
  ppc_stat,
  ( replicate$`2d` | heterosce$`2d` | betabinom$`2d` ) /
    ( replicate$mean | heterosce$mean | betabinom$mean ) /
    ( replicate$sd | heterosce$sd | betabinom$sd ) /
    ( replicate$median | heterosce$median | betabinom$median ) /
    ( replicate$IQR | heterosce$IQR | betabinom$IQR ) +
    plot_layout( heights = c(1.5,1,1,1,1) )
)

# save it
ggsave( plot = last_plot(), filename = here("figs","ppc_stats.jpg"), dpi = 300, width = 12.8, height = 13.3 )


## shape prediction ----

# start with PITs and interval predictions
ppc_shape <-
  
  lapply(
    
    setNames( names(m), names(m) ),
    function(i)
      
      # loop through stats
      lapply(
        
        setNames( c("pit_ecdf","intervals"), c("pit_ecdf","intervals") ),
        function(j)
          
          pp_check( m[[i]], type = j, ndraws = NULL ) + theme( legend.position = "none" )
        
      )
 
  )

# add dens overlay/bar plots
for ( i in names(ppc_shape) ) ppc_shape[[i]]$dens <- pp_check( m[[i]], type = "dens_overlay", ndraws = 100 ) + labs( title = n[[i]] ) + theme( legend.position = "none", plot.title = element_text( hjust = 0.5 ) )

# pull shape PPC together
with(
  ppc_shape,
  ( replicate$dens | heterosce$dens | betabinom$dens ) /
    ( replicate$pit_ecdf | heterosce$pit_ecdf | betabinom$pit_ecdf ) /
    ( replicate$intervals | heterosce$intervals | betabinom$intervals ) +
    plot_layout( heights = c(1.3,1,1) )
)

# save it
ggsave( plot = last_plot(), filename = here("figs","ppc_shape.jpg"), dpi = 300, width = 12.8, height = 13.3 )


# RESULTS ----

# prepare a list for results
post_comp <- list()

# fill-in three-layers of the conditional (i.e., marginal) effects
for ( i in names(m) ) {
  
  post_comp[[i]] <- list()
  post_comp[[i]]$sam <- plot( conditional_effects( m[[i]], method = "posterior_epred", re_formula = NA, effects = "time:GBA" ), points = T )
  post_comp[[i]]$pop <- plot( conditional_effects( m[[i]], method = "posterior_epred", re_formula = NULL, effects = "time:GBA" ), points = T )
  post_comp[[i]]$obs <- plot( conditional_effects( m[[i]], method = "posterior_predict", effects = "time:GBA" ), points = T )
  
}

# prepare a 3x3 grid
with(
  post_comp,
  ( replicate$sam$`time:GBA` | heterosce$sam$`time:GBA` | betabinom$sam$`time:GBA` ) /
    ( replicate$pop$`time:GBA` | heterosce$pop$`time:GBA` | betabinom$pop$`time:GBA` ) /
    ( replicate$obs$`time:GBA` | heterosce$obs$`time:GBA` | betabinom$obs$`time:GBA` ) +
    plot_layout( guides = "collect" ) +
    plot_annotation(
      tag_levels = "A",
      title = "GBA vs non-GBA post-surgery cognitive decline",
      subtitle = "Plot shows marginalized post-surgery cognitive decline of GBA vs non-GBA patients estimated\nby three statistical models (columns) and three levels of prediction (rows)",
      theme = theme( plot.title = element_text( hjust = 0.5, face = "bold" ),  plot.subtitle =  element_text(hjust = 0.5) )
    )
)

# save it
ggsave( plot = last_plot(), filename = here("figs","marginal_effects.jpg"), dpi = 300, width = 12.8, height = 8.32 )
