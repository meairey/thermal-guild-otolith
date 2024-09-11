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

## Experienced temperature
ex.temp = function(d18Water,d18Otolith){
  
  a = (d18Otolith + 1000) / (d18Water + 1000)
  a.ln = 1000*log(a)
  temp = 20690 / (a.ln + 41.69) - 273
  return(temp)
}

# ex.temp(d18Water, d18Otolith) # example

## Little Moose d18O

d18Water = -37.18

## Read in data --------------

oto = read.csv("Airey_SIA_database_final_measurement.csv") %>%
  separate(ISO_YSAMP_N, into = c("SIC", "WATER", "YEAR", "YSAMP")) %>%
  filter(WATER == "LML",
         YEAR > 2018, 
         MATERIAL == "OTOLITH") %>%
  rename("CODE" = "TAXA") %>%
  mutate(experienced_temp = ex.temp(d18Water = d18Water , d18Otolith = D18O))

## Plot experienced temperatures 

## Remember that this is mostly genus specific... so this isnt a good thing to do

oto %>% 
  ggplot(aes(x = CODE, y = experienced_temp)) +
  geom_boxplot()


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

overlap.dat = overlap(dat,1,5, posterior) %>% as.data.frame() %>% ## Overlap function
  pivot_longer(1:length(.[1,])-1, names_to = "post", values_to = "overlap" ) %>%
  separate(post, into = c("text", "comm", "post")) %>%
  mutate(post = as.numeric(post)) %>%
  separate(SppPair, into = c("s1", "v", "s2")) %>%
  select(-text,-v, -comm)
  

save(file = "Data/overlap.RData", overlap_data)

overlap.dat %>%
  mutate(overlap = as.numeric(overlap)) %>%
  group_by(s1, s2) %>%
  summarize(mean_overlap = mean(overlap)) %>%
  ggplot(aes(x = s1, y = s2, fill = mean_overlap)) + 
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  labs(fill = "Mean Ovelap") 



## Niche area


## Ellipse area just use full posterior with full length of species/ellipses  
ellipse.area = siberEllipses(posterior) %>% 
  as.data.frame() %>%
  rename_with(~ names(posterior)) %>% 
  mutate(post_n = seq(1:length(.[,1]))) %>%
  pivot_longer(1:length(posterior), 
               names_to = "comm",
               values_to = "area") %>%
  separate(comm, into = c("community", "group")) %>%
  
  group_by(post_n, community) %>%
  mutate(total_area = sum(area)) %>%
  mutate(relative_area = area / total_area) %>% 
  group_by(community, group) 
# Function to calculate boxplot percentiles

bp.pctiles = function (x, probs = c(0.05, 0.25, 0.5, 0.75, .95)) {
  r <- quantile(x, probs = probs, na.rm = TRUE)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


ellipse.area %>% 
  mutate(group = (as.factor(group))) %>%
  mutate(group = as.numeric(group)) %>%
  left_join(legend) %>%
  ggplot(aes(x = common, y = area)) + 
  stat_summary(fun.data=bp.pctiles, geom="boxplot", width= .5) +
  theme_minimal(base_size = 13) + 
  scale_x_discrete(labels = label_wrap_gen(width = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), 
        axis.title.x = element_blank()) + 
  ylab("Niche Area") 
  




## Guild clustering --------------------------------------
## K means clustering on the d18O 

# Sample data
df = oto %>% select(CODE, D18O, D13C)

# Normalize the temperature column (optional but helpful for clustering)
df$temperature_scaled <- scale(df$D18O)

# Run k-means clustering
set.seed(123)  # Set seed for reproducibility
k <- 3  # Choose the number of clusters (guilds)
clusters <- kmeans(df$temperature_scaled, centers = 3)

# Add the cluster labels to the data frame
df$guild <- as.factor(clusters$cluster)


species.guild =df %>% group_by(CODE, guild) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CODE) %>%
  filter(n == max(n)) %>%
  ungroup() %>% 
  select(-n)

# View the clustered data
dat %>% 
  left_join(species.guild) %>%
  select(guild, common, d180_mean) %>%
  unique() %>% 
  ggplot(aes(x = guild, fill = reorder(common, -d180_mean))) +
  geom_bar()  +
  theme_minimal(base_size = 14) +
  scale_fill_manual("Species", values = legend$col)+
  scale_x_discrete( labels = c("1" = "Cold", "2" = "Cool", "3" = "Warm")) + 
  theme(axis.title.x = element_blank())
