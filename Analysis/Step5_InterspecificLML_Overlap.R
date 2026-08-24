## Setup ------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/thermal-guild-otolith")
source("Analysis/isotope_functions_thermalguildotolith.R") ## data frames, functions, and libraries get loaded in

load(file = "Data/RData/SIBERSetup.RData")
#### Overlap -------------------------


# Get unique species and their families
species_order = legend %>%
  select(CODE, common)%>%
  arrange(common) %>%
  pull(CODE)
## Now expand the grid to get a subsettable matrix
pair_df = expand.grid(species_order, species_order, stringsAsFactors = FALSE)
pair_df = pair_df[pair_df$Var1 != pair_df$Var2, ] ## Remove same species overlap pairs

pair_df = pair_df %>%
  mutate(Var1 = factor(Var1, levels = species_order ), 
         Var2 = factor(Var2, levels = species_order)) %>%
  filter(as.numeric((Var1)) > (as.numeric(Var2))) %>% 
  mutate(C1 = 1, C2 = 1) 
overlap_subset = 1
# --------------------- Overlap Parallelized ------------------------

## Run the overlap calculations using the custom overlap function that interfaces w/ SIBER

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

overlap_list
## Load in data so you don't have to run this every time
load(file = "Data/RData/overlap.list.RData")
## Join together all components of the overlap_list for one wide dataframe
overlap_data = Reduce(full_join, overlap_list) %>%
  select(Sp1, Sp2, everything()) %>%
  pivot_longer(3:length(.[1,]),
               names_to = "post", 
               values_to = "Values") %>%
  mutate(post = parse_number(post)) %>%
  na.omit() %>%
  rename(overlap = Values, 
         s1 = Sp1, 
         s2 = Sp2)

# save(file = "Data/overlap.RData", overlap_data) # Save if desired
# load(file = "Data/overlap.RData") # load if desired

##### Figure 2A: Overlap --------------
## Visualizing the overlap data for the thermal niches
slate_copper_cb <- c(
  "#2F3E46",
  "#4A6C6F",
  "#9DB8A0",
  "#E9C46A",
  "#BC6C25"
)

overlap.graph = overlap_data %>%
  mutate(overlap = as.numeric(overlap)) %>%
  group_by(s1, s2) %>%
  summarize(mean_overlap = mean(overlap)) %>%
    left_join(legend, by = c("s1" = "CODE")) %>% 
  rename(s1_common = common) %>%
  left_join(legend, by = c("s2" = "CODE")) %>%
  rename(s2_common= common) %>%
  #ggplot(aes(x = s1_common, y = s2_common, fill = mean_overlap)) + 
  ggplot(aes(x = reorder(s1_common, d180_mean.x), y = reorder(s2_common, d180_mean.y), fill = mean_overlap)) + 
  geom_tile() +
  scale_fill_gradientn(colors = slate_copper_cb) +
  theme_minimal(base_size = 11) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(fill = "Mean Overlap") 
overlap.graph

ggsave(overlap.graph, file = "Graphics/Figure_2A_Overlap.Summary.pdf", width = 4, height = 4)

## Credible intervals of overlap data

overlap.table = overlap.dat %>%
  mutate(overlap = as.numeric(overlap)) %>%
  group_by(s1, s2) %>%
  summarize(mean_overlap = round(mean(overlap), digits = 2), 
            upper = round(quantile(overlap, .975), digits = 2),
            lower = round(quantile(overlap, .025), digits = 2)) %>%
  ungroup() %>% 
  left_join(legend, by = c("s1" = "CODE")) %>% 
  rename(s1_common =common) %>%
  left_join(legend, by = c("s2" = "CODE")) %>%
  rename(s2_common= common) %>%
  mutate(CI = paste(mean_overlap, " [", lower, ", ",upper, "]" ) ) %>%
  select(s1_common, s2_common, CI ) %>%
  pivot_wider(values_from = CI, names_from = s2_common)

#write.csv(overlap.table, file = "Data/overlap.table.csv", row.names = F)