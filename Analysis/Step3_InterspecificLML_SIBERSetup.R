## Setup ------
setwd("C:/Users/monta/OneDrive - Airey Family/GitHub/thermal-guild-otolith")
source("Analysis/isotope_functions_thermalguildotolith.R") ## data frames, functions, and libraries get loaded in

## SIBER -----------------------------
#### SIBER Niches --------------

# SIBER data setup
siber.data = df %>%
  arrange(CODE) %>%
  rename(
         "iso1" = "D18O_VPDB",
         "iso2" = "D13C.otolith") %>%
  mutate(community = 1)  %>% 
  select(iso1, iso2, group, community)
# siber object and posterior generation
siber.setup = data_setup(siber.data, 1) # custom function generates and formats posterior
spp=length(names(siber.setup[[2]])) # total number of species
save(siber.setup, file = "Data/RData/SIBERSetup.RData")
##### Figure 1: LML Fish -----------------
# Plot setup (ellipses get added to this)
p = ggplot() + geom_point(siber.setup[[3]] %>%
                              as.data.frame() %>%
                              left_join(legend),
                          mapping =  aes(x = iso2, y =iso1 , 
                                         col =  reorder(CODE, d180_mean),
                                         alpha = guild), 
                          show.legend = FALSE) + 
    theme_minimal()

ellip = array(0, dim=c((n.points),length(c(0,1)),spp))


# Create an array that includes each species' ellipse
for(i in 1:spp){
  ellipse_data =  ellipse::ellipse(x = matrix(c(mean(siber.setup[[2]][[i]][,1]),
                                                mean(siber.setup[[2]][[i]][,2]),
                                                mean(siber.setup[[2]][[i]][,3]),
                                                mean(siber.setup[[2]][[i]][,4])),
                                                2,2),
                                   centre = c(mean(siber.setup[[2]][[i]][,5]),
                                              mean(siber.setup[[2]][[i]][,6])),
                                   level = .4,
                                   npoints = n.points)
  ellip[,,i] = ellipse_data
}


# Format ellipse data for plotting
ellip.graph = data.frame(xax = as.vector((ellip[,2,])), yax = as.vector((ellip[,1,])), CODE = (rep(as.factor(unique((siber.setup[[3]]$group))), each = 1000))) %>%
  rename(group = CODE) %>%
    left_join(legend %>%
                mutate(group = as.factor(group))) %>%
    arrange(d180_mean) ## line is important for legend to work



## Final graph with isotope niches
p = p + geom_path(data = ellip.graph, aes(x = xax,y = yax,
                                          color = CODE,
                                          alpha = guild),
                  key_glyph = "rect", 
                  lwd = 1) +
  ylab(expression(paste(delta ^18, "O"))) +
  xlab(expression(paste(delta ^13, "C"))) + 
  scale_color_manual(values = (legend %>% arrange(d180_mean))$color, 
                       labels =(legend %>% arrange(d180_mean))$common,
                       name = "Species") +
  theme(text = element_text(size = 13)) + 
  scale_alpha_manual("Guilds", values = c(1,1,1))  + 
  guides(alpha = "none")
  
print(p)

#ggsave(p, file = "Graphics/Figure_1_LML.Web.pdf", height = 5, width = 6, dpi = 500)
#ggsave(p, file = "Graphics/Figure_1_LML_cold.Web.png", height = 5, width = 6, dpi = 500)
  




## Supplemental Table S1 ----------- 
# Morphological table of LML fish

df = read.csv("Data/Metadata and Measurement data/isotopes_lengths.csv") %>%
  mutate(length_s = LENGTH - ave(LENGTH, SPECIES))

df %>% 
  ggplot(aes(x = CODE, y= LENGTH)) +
  geom_violin() +
  geom_jitter(width = .1) + 
  theme_minimal(base_size = 10) + 
  scale_y_log10() 

df %>% 
  ggplot(aes(x = LENGTH,y = D18O_VPDB)) +
  geom_point() + 
  theme_minimal(base_size = 14) + 
  geom_smooth(method = lm, se = F, col = "black", lty = 2) +
  facet_wrap(~SPECIES, scales = "free") + 
  xlab("Length (mm)") +
  ylab("d18O (otolith)")

library(lmerTest)
m1 <- lmerTest::lmer(
  D18O_VPDB ~ length_s + (1| SPECIES),
  data = df,
  REML = TRUE
)

summary(m1)
