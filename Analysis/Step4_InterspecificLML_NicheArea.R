## Setup ------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/thermal-guild-otolith")
source("Analysis/isotope_functions_thermalguildotolith.R") ## data frames, functions, and libraries get loaded in

load(file = "Data/RData/SIBERSetup.RData")
#### Niche area ----------------------
## Ellipse area just use full posterior with full length of species/ellipses  
ellipse.area = siberEllipses(siber.setup[[2]]) %>% 
  as.data.frame() %>%
  rename_with(~ names(siber.setup[[2]])) %>% 
  mutate(post_n = seq(1:length(.[,1]))) %>%
  pivot_longer(1:length(siber.setup[[2]]), 
               names_to = "comm",
               values_to = "area") %>%
  separate(comm, into = c("community", "group")) %>%
  group_by(post_n, community) %>%
  mutate(total_area = sum(area)) %>%
  mutate(relative_area = area / total_area) %>% 
  group_by(community, group) 


#save(ellipse.area, file = "Data/ellipse.area.RData")
#load( file = "Data/ellipse.area.RData")
##### Figure 2B: Niche Area ---------------------------
ellipse.area %>% 
  mutate(group = as.numeric(group)) %>%
  left_join(legend) %>%
  ggplot(aes(y = reorder(common,  d180_mean), x = area,
             fill = reorder(common,  d180_mean)
             ), 
        ) + 
  stat_summary(fun.data=bp.pctiles, geom="boxplot", width= .5) +
  theme_minimal(base_size = 11) + 
  theme(axis.title.y = element_blank(), 
        legend.position = "none") + 
  scale_fill_manual("Taxon",values = (legend %>% arrange(d180_mean))$color, 
                    labels =(legend %>% arrange(d180_mean))$common) +
  
  xlab("SEAc (95% CI)") -> Figure_2B

#ggsave(Figure_2B, file = "Graphics/Figure_2B_AreaSummary.pdf", width = 2.5, height = 4, units = "in", dpi = 500)


### SD of d18O as metric for thermal breadth use because of metabolic complication with d13C

## I don't like this as a permanent fix because without SEAc the whole thing becomes very dependent on the sample size again. So, I'm thinking of a supplemental figure that shows that sd of d18O roughly scales with the niche area. plus d13C isnt totally independent of temperature
mean_ellipse_area = ellipse.area %>% 
  group_by(group) %>%
  summarize(mean_area = mean(area)) %>% 
  mutate(group = as.numeric(group)) %>%
  left_join(legend) %>%
  rename(SPECIES = CODE)
sd_O = df %>% group_by(SPECIES) %>% 
  summarize(mean_d18O = mean(D18O_VPDB, na.rm = T), 
            sd_d18O = sd(D18O_VPDB, na.rm = T)) %>% 
  left_join(mean_ellipse_area)

lm(mean_area ~ sd_d18O, data = sd_O) %>% summary()
  
eq_label = paste0("italic(y) == ", 0.03, " + ", 3.82, " * italic(x)  ~ (italic(R)^2 == ", 0.59, ")")
sd_O %>%
  ggplot(aes(x = sd_d18O, y = mean_area, label = SPECIES)) +
  geom_smooth(method = "lm", lty = "dashed", alpha = .5, col = "black") + 
  geom_text() + 
  theme_minimal() + 
  xlab(expression(SD ~ delta^18 * O)) + 
  ylab("Mean SEAc") + 
  annotate("text", x = .5, y = 6, label = eq_label, parse = TRUE, size = 3)

ggsave(file = "Graphics/Figure_S3_SEAcvSD.pdf", height = 4, width = 4)
 
## Guild clustering --------------------------------------
## K means clustering on the d18O 

# Sample data
cluster.dat = df %>% select(CODE, D18O_VPDB, D13C.otolith)

# Normalize the temperature column (optional but helpful for clustering)
cluster.dat$temperature_scaled = scale(df$D18O)

# Run k-means clustering

#Compute WCSS for different values of k
wcss = numeric(10)  # Store WCSS values

# Compute the elbow method and visualize optimal k
fviz_nbclust(cluster.dat$temperature_scaled, kmeans, method = "wss")
k = 3  # Choose the number of clusters (guilds)
# Split out the clusters
clusters = kmeans(cluster.dat$temperature_scaled, centers =2)

# Add the cluster labels to the data frame
cluster.dat$guild = as.factor(clusters$cluster)

## Format the guild cluster data frame
species.guild = cluster.dat %>% group_by(CODE, guild) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CODE) %>%
  filter(n == max(n)) %>%
  ungroup() %>% 
  select(-n)

## Format data for graphing guilds
g = cluster.dat %>%
  left_join(legend %>% 
              select(-guild), by = c("CODE")) %>%
  group_by(common, d180_mean) %>%
  mutate(sum =n()) %>% 
  ungroup() %>% 
  group_by(common, sum, color,d180_mean, guild) %>%
  summarize(guild_sum = n()) %>%
  mutate(proportion = guild_sum / sum) %>%
  arrange(guild, proportion)%>% 
  mutate(order = case_when(guild == 2 ~ proportion, guild == 1 ~ 1 - proportion)) %>%
  arrange(order, -d180_mean) %>%
  ungroup() %>%
  mutate(order = c(1:15))

##### Figure 3: Guilds ----------------------
  
g %>%
  ggplot(aes(x = reorder(common, order), y = proportion, fill = guild)) + 
  geom_bar(stat = "identity") + 
  theme_minimal(base_size = 11) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) + 
  scale_fill_manual("Thermal Guild", values = col$color[c(2,6)], 
                    labels = c("1" = "Cold", "2" = "Cool")) + 
  ylab("Proportion \nof Samples") -> Figure_3

ggsave(Figure_3, file = "Graphics/Figure_3_Guilds.pdf", width = 6, height = 2, dpi = 500 )
