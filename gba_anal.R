
rm( list = ls() ) # clear environment
s = 1 # seed for reproducibility
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

library(here) # directory management
library(tidyverse) # data wrangling
library(readODS) # reading .ods data
library(openxlsx) # reading .xlsx data
library(brms) # model fitting
library(tidybayes) # posterior manipulation
library(patchwork) # plot grids
#library(gt) # tables

theme_set( theme_minimal() )

# prepare folders
sapply( c("figs","mods"), function(i) if (!dir.exists(i) ) dir.create(i) )

# READ DATA ----

d0 <- read.csv( here("_raw","psych_longit.csv"), sep = "," )
d.roma <- read_ods( here("_raw","ROMA-ITEMPOVyeteni_DATA_2024-05-02.ods") )
d.gba <- read.xlsx( here("_raw","GBA_iTEMPO.xlsx") )

# pre-process
d1 <-
  
  d0 %>%
  filter( id %in% substr(d.roma$ipn, 1, 6) ) %>% # keep only patients with genetic data
  filter( stn_dbs == 1 ) %>% # keep only STN DBS patients
  filter( complete.cases(drsii) ) %>% # keep rows with DRS-2 only
  mutate( GBA = factor( ifelse( id %in% d.gba$IPN, 1, 0 ) ), .after = sex ) %>% # add GBA indicator
  mutate( md_time = median(stimtime_years[event == "screening"], na.rm = T), time = stimtime_years - md_time ) %>% # shift time such that median pre-surgery assessment is t0
  select( id, event, GBA, time, md_time, drsii ) # select variables of interest

# plot assessment distributions
table( d1[ , c("id","event","GBA") ] ) %>%
  as.data.frame() %>%
  mutate(
    `GBA: ` = case_when(
      GBA == 1 & Freq == 1 ~ "GBA+",
      GBA == 0 & Freq == 1 ~ "GBA-",
      .default = NA
    ),
    event = factor(event, levels = c("screening", paste0( "y",seq(1,19,2) ) ), ordered = T )
  ) %>%
  filter( complete.cases(`GBA: `) ) %>%
  
  ggplot() +
  aes( y = id, x = event, fill = `GBA: ` ) +
  geom_tile( colour = "white" ) +
  scale_fill_manual( values = c("grey82", "#0072B2" ) ) +
  labs(
    title = "DRS-2 assessment distributions across patients",
    subtitle = "coloured fields indicate cases present in the data",
    y = NULL,
    x = NULL
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text( size = 16, face = "bold", hjust = .5 ),
    plot.subtitle = element_text( size = 14, hjust = .5 )
  )

# save it
ggsave( plot = last_plot(), filename = here("figs","assessment_distributions.jpg"), dpi = 300, width = 8, height = 16.2 )


# REGRESSION MODELS OF DRS-2 ----

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
    replicate = "Normal(mu, sigma) likelihood\n\nmu ~ 1 + GBA * Time + (1 + Time | ID),\nsigma ~ 1",
    heterosce = "Normal(mu, sigma) likelihood\n\nmu ~ 1 + GBA * Time + (1 + Time | ID),\nsigma ~ 1 + Time",
    betabinom = "Beta-Binomial(mu, phi) likelihood\n\nmu ~ 1 + GBA * Time + (1 + Time | ID),\nphi ~ 1 + Time"
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
    plot_layout( heights = c(1.6,1,1,1,1) )
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
