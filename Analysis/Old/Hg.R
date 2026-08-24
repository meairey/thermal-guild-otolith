library(tidyverse)
`%nin%` = Negate(`%in%`)
library(wesanderson)
library(lme4)
source("isotope_functions_thermalguildotolith.R")
## Read in data --------------
oto = read.csv("../1.clean_isotope/SI_MEASUREMENT.csv") %>%
  left_join(read.csv("../1.clean_isotope/SI_SAMPLE.csv"), by = "ISO_YSAMP_N") %>%
  filter(SAMPLE_TYPE == "OTOLITH",
         #WATER == "LML" , ## commented out for steve's fish
         YEAR > 2018)  %>%
  rename("CODE" = "TAXON") %>%
  mutate(experienced_temp = ex.temp(d18Water = d18Water , d18Otolith = D18O_VPDB))%>% ## Function found above

  filter(ITEM_N != "22-39") %>% ## Removing smaller LLS that may have hatchery influence
  select(-COMMENT.x, -COMMENT.y)



## Load in new data that has CN for those taxa

oto.samps = (oto %>% filter(SAMPLE_TYPE == "OTOLITH"))$ISO_FISH_N ## IDs for fish w/ otolith measurements

tissue = read.csv("../1.clean_isotope/SI_MEASUREMENT.csv") %>% ## Tissue samples for fish that had otoliths measured
  left_join(read.csv("../1.clean_isotope/SI_SAMPLE.csv"), by = "ISO_YSAMP_N") %>%
  filter(SAMPLE_TYPE == "TISSUE", 
         #WATER == "LML" , # commented out for steve's fish
         YEAR > 2018, 
         GROUP == "FISH", 
         corrected != "duplicate delete", 
        # ISO_FISH_N %in% oto.samps ## commented out for steve's fish
        )  %>% 
  rename("CODE" = "TAXON") %>%
  mutate(experienced_temp = ex.temp(d18Water = d18Water , d18Otolith = D18O_VPDB)) %>%
  select(-COMMENT.x, -COMMENT.y)

new = read.csv("../1.clean_isotope/iso_measurement.csv") %>%
  #filter(ISO_FISH_N %in% oto.samps) %>%  ## commented out for steve's fish
  left_join(read.csv("../1.clean_isotope/isotope_sample.csv"), by = "ISO_YSAMP_N") %>%
  select(-COMMENT.x, -COMMENT.y)


data.full = rbind(oto, tissue) %>%
  select(ISO_FISH_N, DATE_COL,  ITEM_N, CATEGORY, GROUP, CODE, SAMPLE_TYPE, YEAR, D13C, D18O_VPDB, D15N) %>% 
 
  rbind(., new %>% 
          select(ISO_FISH_N, DATE_COL, ITEM_N, CATEGORY, GROUP, TAXON, SAMPLE_TYPE, YEAR, D13C, D18O_VPDB, D15N )  %>%
         rename("CODE" = "TAXON")  ) 


x = data.full %>%
 select(ISO_FISH_N, DATE_COL, SAMPLE_TYPE, CODE , D13C, D18O_VPDB ) %>% 
  unique() %>%
  filter(SAMPLE_TYPE == "OTOLITH") %>%
  rename(D13C.oto = D13C) %>%
  select(ISO_FISH_N, CODE, D13C.oto, D18O_VPDB)

y = data.full %>%
 select(ISO_FISH_N, DATE_COL,  SAMPLE_TYPE,CODE,   D15N, D13C ) %>% 
  unique() %>%
  filter(SAMPLE_TYPE == "TISSUE") %>%
  rename(D13C.tissue = D13C) %>%
  select(ISO_FISH_N, CODE, D13C.tissue, D15N)





## Steve's isotope data
cat = read.csv("Data/lake_categories.csv")



hg = read.csv("Data/Data_Steve.csv") %>%
  select(AFRP_fish_id, length_mm, SeLengthCorrected_ug_g_dry_weight, HgLengthCorrected_ug_g_dry_weight)

combined = data.full %>%
  mutate(experienced_temp = ex.temp(d18Water = d18Water , d18Otolith = D18O_VPDB)) %>% 
  left_join(hg, by = c("ISO_FISH_N" = "AFRP_fish_id"))  %>%
  filter(ISO_FISH_N %in% hg$AFRP_fish_id) %>%
  separate(ISO_FISH_N, into = c("sp","WATER","date","gear","num")) %>%
  left_join(cat) %>%
  filter(WATER != "POL")

combined %>%
  ggplot(aes(x = experienced_temp, y = HgLengthCorrected_ug_g_dry_weight, col = OXY_CAT)) + 
  geom_point() + 
  geom_smooth(method = lm, se = F) +
  theme_minimal(base_size = 13) +
  xlab("Growing Temperature (C)") + ylab("Hg (length corrected)") +
  scale_color_manual("Oxythermal\nCategory", values = wes_palette("Moonrise2", n = 4))


install.packages("MuMIn")
combined %>% 
  filter(OXY_CAT == "Constrained", 
         experienced_temp > 10
        ) %>% 
  lmer(data = ., HgLengthCorrected_ug_g_dry_weight ~ experienced_temp + 
       SeLengthCorrected_ug_g_dry_weight + (1|WATER)
     ) %>%
  AIC()
  #MuMIn::r.squaredGLMM()
  summary()

combined %>% 
  filter(OXY_CAT == "Constrained", 
         #experienced_temp > 10
         ) %>% 
  ggplot(aes(x = experienced_temp, y = HgLengthCorrected_ug_g_dry_weight)) + 
  geom_point() + 
  geom_smooth(method = lm)


## Variation of temperature across lake categories

combined %>%
  filter(sp == "ST") %>%
  ggplot(aes(x = OXY_CAT, y = experienced_temp)) + 
  geom_boxplot() +
  geom_point()

library(mgcv)

# fit GAM
## The GAM suggests that the relationship is pretty much linear so I went back and just used the lmer above with lake as a random effect
gam_dat = combined %>% 
  filter(OXY_CAT == "Constrained"#, 
        # experienced_temp > 10
         ) %>%
  mutate(WATER = as.factor(WATER))
  
gam_fit <- gam(
  HgLengthCorrected_ug_g_dry_weight ~ s(experienced_temp, k = 5) + 
   SeLengthCorrected_ug_g_dry_weight  + s(WATER, bs="re"),
  data = gam_dat,
  method = "REML"
)

# quick model check
summary(gam_fit)
gam.check(gam_fit)



newdat <- data.frame(
  experienced_temp = seq(
    min(gam_dat$experienced_temp, na.rm = TRUE),
    max(gam_dat$experienced_temp, na.rm = TRUE),
    length.out = 200
  ),
  SeLengthCorrected_ug_g_dry_weight = mean(
    gam_dat$SeLengthCorrected_ug_g_dry_weight,
    na.rm = TRUE
  )
)
# predictions + SE
pred <- predict(gam_fit, newdat, se.fit = TRUE)


newdat$fit <- pred$fit
newdat$se  <- pred$se.fit

ggplot() +
  geom_point(data = gam_dat, aes(experienced_temp, HgLengthCorrected_ug_g_dry_weight), size = 2) +
  geom_line(data = newdat, aes(experienced_temp, fit), linewidth = 1) +
  geom_ribbon(
    data = newdat,
    aes(
      experienced_temp,
      ymin = fit - 2 * se,
      ymax = fit + 2 * se
    ),
    alpha = 0.2
  ) +
  labs(
    x = "Experienced temperature",
    y = "Hg (length-corrected)"
  ) +
  theme_classic()



### Isotope investigations


tiss_oto = left_join(x, y) %>%
  filter(CODE == "ST") %>%
  na.omit() %>%
  mutate(experienced_temp = ex.temp(d18Water = d18Water , d18Otolith = D18O_VPDB),
        D13C.tissue = as.numeric(D13C.tissue)) %>%
  left_join(hg %>%
              select("AFRP_fish_id", "length_mm"), by = c("ISO_FISH_N" = "AFRP_fish_id")) %>%
  separate(ISO_FISH_N, into = c("sp","WATER","date","gear","num")) 


library(lmerTest)

tiss_oto %>% 
  filter(experienced_temp > 10) %>%
  ggplot(aes(x = experienced_temp, y= as.numeric(D13C.tissue), col = WATER)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F)



tiss_oto %>% 
  filter(experienced_temp >10) %>%
  ggplot(aes(x = experienced_temp ,y= as.numeric(D15N), col = WATER)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F)

tiss_oto %>% 
 # filter(experienced_temp >10) %>%
  lmer(data = ., experienced_temp ~ D15N + length_mm + (1|WATER) ) %>%
  summary()


tiss_oto %>% 
 # filter(experienced_temp >10) %>%
  lmer(data = ., experienced_temp ~ D13C.tissue + (1|WATER)) %>%
  summary()

tiss_oto %>%
  ggplot(aes(x =D13C.tissue, y = D15N, col = WATER)) +
  geom_point() +
  stat_ellipse(level = .4)

tiss_oto %>%
  ggplot(aes(x = experienced_temp))

## Standardizing for length
tiss_oto = tiss_oto %>%
  na.omit() %>%
  mutate(D15N_std = resid(lm(D15N ~ length_mm, data = .))) 

tiss_oto %>% 
  filter(experienced_temp > 10) %>%
  lmer(data = ., experienced_temp ~ D15N_std +(1|WATER) ) %>%
  summary()


