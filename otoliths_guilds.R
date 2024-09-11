library(tidyverse)
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/LML_SMB_removal/Data/ALC/")
library(viridis)

## Guilds ----------------

cold = c("LT","CS", "RS","SS")
salmonid = c("LT","CS", "ST", "LLS")
warm = c("LT","CS", "CC", "PS", "MM", "SMB")
common.names = c("Brook trout", "Common shiner", 
                 "Pumpkinseed", "Smallmouth bass",
                 "White sucker", "Creek chub", 
                 "Mudminnow", "Atlantic salmon", 
                 "Slimy sculpin", "Rainbow smelt", 
                 "Lake trout") %>% rev()
## Read in data --------------

oto = read.csv("Airey_SIA_database_final_measurement.csv") %>%
  separate(ISO_YSAMP_N, into = c("SIC", "WATER", "YEAR", "YSAMP")) %>%
  filter(WATER == "LML",
         YEAR > 2018, 
         MATERIAL == "OTOLITH") %>%
  rename("CODE" = "TAXA")

## Assign colors ------------
col = oto %>%
  group_by(CODE) %>%
  summarize(d180_mean = mean(D18O)) %>%
  arrange(-d180_mean) %>%
  mutate(color = viridis(n = 11)) 

## Assign legend for summary graph -----------------
legend = data.frame(CODE = col$CODE, 
                    common = common.names) %>%
  left_join(col) %>%
  arrange(CODE) %>%
  mutate(group = as.numeric(as.factor(CODE))) %>%
  arrange(-d180_mean)

## Join legend with data
dat = oto %>% 
  left_join(legend) %>%
  arrange(-d180_mean)

## Visualize Data

dat %>% 
  ggplot(aes(x = D13C, y = D18O, col =reorder(common, -d180_mean)))  + 
  geom_point(size = 2) +
  scale_color_manual(values = legend$col )  +
  theme_minimal(base_size = 14) +
  stat_ellipse(lwd = .8, level = .4, key_glyph = "rect") +
  labs(col = "Species") +
  xlim(-18, -8.5) +
  ylim(-10, -4)+
  theme(plot.margin = margin(t = 10, r = 2, b = 10, l = 20))


## SIBER
library(SIBER)

dat = oto %>%
  arrange(CODE) %>%
  rename("group" = "CODE", 
         "iso1" = "D18O",
         "iso2" = "D13C") %>%
  mutate(community = 1)  %>% 
  select(iso1, iso2, group, community)

siber.dat = createSiberObject(dat)

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3


posterior <- siberMVN(siber.dat, parms, priors)

setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/thermal-guild-otolith")
source("isotope_functions_update.R")
## Overlap

overlap_list = list() 
for(h in 1:2){
  
  overlap_list[[h]] = overlap(dat,c(1:2),100, posterior) %>% as.data.frame()
}



# Code to reduce list to matrix 
overlap_data = Reduce(full_join, overlap_list) %>%
  select('SppPair', everything())

save(file = "Data/overlap.RData", overlap_data)

overlap_coocur = overlap_data %>% 
  pivot_longer(2:length(.[1,]),
               names_to = "Community", 
               values_to = "Values") %>% 
  na.omit() %>%
  rename("Species_Pair" = `SppPair`) %>% 
  separate(Community, into = c("c", "Community", "post")) %>%
  select(-c) %>%
  group_by(Species_Pair, post) %>%
  mutate(Values = as.numeric(Values))%>%
  summarize(med_overlap = mean(Values, na.rm = T)) %>%
  ungroup() %>%
  group_by(Species_Pair) %>%
  summarize(med_overlap = mean(med_overlap, na.rm = T)) %>%
  #summarize(mean = mean(as.numeric(Values), na.rm = T)) %>%
  separate(Species_Pair, into = c("s1", "v", "s2"), sep = " ") %>% 
  select(-v) %>%
  filter(s1 != s2) %>% rename(Species = "s1") %>% 
  mutate(Species = str_replace(Species, "NPD", "PD")) %>% 
  mutate(Species = str_replace(Species, "BNM", "BM")) %>%
  mutate(s2 = str_replace(s2, "NPD", "PD")) %>% 
  mutate(s2 = str_replace(s2, "BNM", "BM")) %>%
  #left_join(prop_comm) %>% 
  full_join(coccur_matrix) %>%
  left_join(dat_spec) %>% 
  left_join(dat_spec2) %>%
  unite("species_pair",family, familys2, remove = F) %>%
  filter(med_overlap != "NA")
