library(tidyverse)
## Analysis across gear types

`%nin%` = Negate(`%in%`)
fish = read.csv("Data/FISH_MEASUREMENT_LML.csv")
sample = read.csv("Data/FISH_SAMPLE_2024.csv")


catch.data = fish %>% left_join(sample) %>%
  filter(GEAR %nin% c("RANG"))


cpue.data = catch.data %>% group_by(YEAR, WATER, GEAR, SITE_N, 
                        DATE_COL,DAY_N, GEAR_CODE, SPECIES,EFFORT) %>%
  summarize(catch = n()) %>%
  ungroup() %>%
  filter(EFFORT != "") %>%
  mutate(EFFORT = parse_number(EFFORT)) %>%
  mutate(EFFORT = case_when(GEAR == "TPN" ~ EFFORT*24, GEAR != "TPN" ~ EFFORT)) %>%
  #mutate(cpue = catch / (EFFORT/60/60)) ## transform to individuals per hour
  mutate(cpue = catch / (EFFORT)) 
  

cpue.data %>%
  filter(SPECIES == "LT") %>%
  ggplot(aes(x = YEAR, y = cpue, col = GEAR)) + 
  geom_smooth(se=F)#+
#  scale_y_log10()

catch.data %>% filter(YEAR< 2000) %>%
  select(YEAR, WATER, GEAR, SITE_N, 
           DATE_COL,DAY_N, GEAR_CODE, SPECIES,EFFORT) %>%
  select(SPECIES) %>% unique()



catch.data %>% 
  filter(GEAR == "TPN") %>% 
  select(TIME_START, TIME_END,EFFORT, EFFORT_UNIT) %>%
  mutate(EFFORT_hours = parse_number(EFFORT) * 24)
