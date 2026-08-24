## Source setup for use across scripts
set.seed(123)  # Set seed for reproducibility
`%nin%` = Negate('%in%')
# ---------------------- Library Loads ---------------------------------------
library(tidyverse)
library(viridis)
library(SIBER)
library(factoextra)
library(lme4)
library(future.apply) ## Package for parallelizing lapply
library(progressr)
#--------------------- Predefined Values ----------------------------------- 
d18Water = -37.18 ## Little Moose d18O
## Jags setup
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors for jags
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3


Nsamples=1000
n.posts <- 1000;
p.ell <- 0.90 # How much data to include? Standard ellipses --> p.ell = .9
n.points = 1000


#--------------------- Data frames ---------------------------------------------

df = read.csv("Data/Metadata and Measurement data/isotopes_lengths.csv") %>%
  mutate(length_s = LENGTH - ave(LENGTH, SPECIES))

## Legend names
common.names = c("Brook trout", "Common shiner", 
                 "Pumpkinseed", "Smallmouth bass",
                 "White sucker", "Creek chub", 
                 "C. Mudminnow", "A. Salmon", 
                 "Slimy sculpin", "Rainbow smelt", 
                 "Lake trout") %>% rev()

## Assign colors ------------
col = df %>%
  group_by(CODE) %>%
  summarize(d180_mean = mean(D18O_VPDB)) %>%
  arrange(-d180_mean) %>%
  mutate(color = viridis(n = 11)) 

## Assign legend for summary graph -----------------
legend = data.frame(CODE = col$CODE, 
                    common = common.names) %>%
  left_join(col) %>%
  arrange(CODE) %>%
  mutate(group = as.numeric(as.factor(CODE))) %>%
  arrange(-d180_mean) %>%
  mutate(guild = case_when(CODE %in% c("LT", "SS") ~ "cold", CODE %in% c("RS", "LLS", "CC", "WS") ~ "mixed",
                           CODE %in% c("MM", "SMB", "PS","CS", "ST")~ "cool"))

## Join legend with data
df = df %>% 
  left_join(legend) %>%
  arrange(-d180_mean) 



#---------------------- Functions ----------------------------------------------

## Plotting function for 95% quantiles
bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, .95)) {
  r <- quantile(x, probs = probs, na.rm = TRUE)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


## Experienced temperature function -------
ex.temp = function(d18Water,d18Otolith){
  
  a = (d18Otolith + 1000) / (d18Water + 1000)
  a.ln = 1000*log(a)
  temp = 20690 / (a.ln + 41.69) - 273
  return(temp)
}

## Overlap function ---------
overlap = function(data_input, comm, dr, posterior, pair_df, names_modified = NULL) {

  # Subset pairs
  name_matrix = pair_df[pair_df$C1 == comm & pair_df$C2 == comm, c("Var1","Var2")]
  
  n_pairs = nrow(name_matrix)
  
  # Preallocate
  S1_first = matrix(NA_real_, nrow = dr, ncol = n_pairs)
  S2_first = matrix(NA_real_, nrow = dr, ncol = n_pairs)
  # Compute overlap
  for (i in seq_len(n_pairs)) {
    
    overlap_res = bayesianOverlap(
      ellipse1 = name_matrix[i, 1],
      ellipse2 = name_matrix[i, 2],
      posterior,
      draws = dr,## Number of draws from the posterior used
      p.interval = 0.95, 
      n = 100 ## Number of points the ellipse is estimated around (higher = more smooth sides)
    )
    
    
    S1_first[,i] = overlap_res[,3] / overlap_res[,1]
    S2_first[,i] = overlap_res[,3] / overlap_res[,2]
  
  }
  # Data frame
  S1_mat = t(S1_first) %>%
    as.data.frame() %>% 
    mutate(Sp1 = name_matrix$Var1,
           Sp2 = name_matrix$Var2)
  S2_mat = t(S2_first) %>%
    as.data.frame() %>% 
    mutate(Sp1 = name_matrix$Var2,
           Sp2 = name_matrix$Var1)
  
  full_olap_mat = rbind(S1_mat, S2_mat)
  
 
   
  # Return the data frame
  return((full_olap_mat))
}




## -------------- Data function ----------------------------

## Removed PD because not sure identification 
## Must put in Siber formatted data_input

data_setup = function(data_input, com_num){
  dat_codes = data_input$group
  data = data_input %>% 
    filter(community == com_num)
  data = data[order(data$group),] %>% as.data.frame() 
  siber.example <- createSiberObject(data)
  posterior <- siberMVN(siber.example, parms, priors)
  
  both = list(siber.example, posterior, data,dat_codes)
  return(both)
}


## ------------------------ Ellipse data function ---------------------



ellip_data = function(numb_species, numb_posts, posterior){
  for(i in 1:numb_species){
    dat = vector()
    for(j in 1:numb_posts){
      ellipse_data =  ellipse::ellipse(x = matrix(c(posterior[[i]][j,1],
                                                    median(posterior[[i]][j,2]),
                                                    median(posterior[[i]][j,3]),
                                                    median(posterior[[i]][j,4])),
                                                  2,2), 
                                       centre = c(median(posterior[[i]][j,5]),
                                                  median(posterior[[i]][j,6])), level = .95, 
                                       npoints = n.points) %>% as.data.frame() %>%
        summarise_all(list("min"=min, "max"=max)) %>% as.matrix()
      
      dat = rbind(dat, ellipse_data)
    }
    
    ellip[,,i] = dat
    
  }
  
  return(ellip)
}

