library(tidyverse)
library(lubridate)
source("isotope_functions_thermalguildotolith.R")
##Replicating Demetra's graphs


lake_characteristics = read.csv("Data/TempDO Data/lake_metadata.csv") %>%
  group_by(afrp_abbrev) %>% 
  slice(1)
otolith_dem = read.csv("Data/Metadata and Measurement data/ThesisOtolith_Demetra_Processed.csv")

## 2022 LML ST data 
## Water chemistry
d18Water = -37.18 ## Little Moose d18O

df_lml = read.csv("Data/Metadata and Measurement data/isotopes_lengths.csv") %>%
  mutate(length_s = LENGTH - ave(LENGTH, SPECIES)) %>%
  filter(CODE == "ST", 
         grepl("BEF", ISO_FISH_N)) %>%
  separate(ISO_FISH_N, into = c("st", "afrp_abbrev", "SAMPLE_DATE", "gear", "num"), remove = F) %>%
  mutate(SAMPLE_DATE = mdy(SAMPLE_DATE)) %>%
  mutate(SAMPLE_DATE = format(SAMPLE_DATE, "%m/%d/%Y")) %>%
  mutate(CalcWaterTemp = ex.temp(d18Otolith = D18O_VPDB, d18Water = d18Water)) %>%
  mutate(lake_type = "Buffered",
         ORIGIN = "NRP",
         Lake = "Little Moose Lake") %>%
  rename(FISH_N = ISO_FISH_N, 
         FISH_LENGTH_mm = LENGTH) %>%
  select(-X, -SPECIES)


## Maximum temperature in lake type 

df_cross %>% 
  group_by(lake_type) %>%
  summarize(max_temp.min = min(max_temp), 
            max_temp.max = max(max_temp))


##
df_cross = left_join(otolith_dem, lake_characteristics) %>%
  filter(Lake != "Pico Lake") %>%
  mutate(hatchery_cont = .8/ OTOLITH_WEIGHT_mg) %>%
  filter(!(ORIGIN == "STK" & hatchery_cont >= .15)) %>%
  select(hatchery_cont, OTOLITH_WEIGHT_mg, ORIGIN, everything()) %>%
  arrange(ORIGIN, hatchery_cont) %>% 
  bind_rows(df_lml)

### Lake by category
df_cross %>% ggplot(aes(y = CalcWaterTemp, x = Lake)) +
  geom_hline(yintercept = 15, lty = "dashed") +
  facet_wrap(~lake_type,strip.position = "bottom",
             ncol = 4, 
             scales = "free_x") + 
  geom_boxplot(fill = NA, outliers = F) +
  geom_jitter(aes(shape= ORIGIN , col = ORIGIN ),
             size = 3, width = .2, height = 0) + 
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = c(0.89, 0.8),   # x, y inside plot (0 to 1)
  legend.background = element_rect(fill = "white", color = "black"),
  legend.box.background = element_rect(color = "black")) +
  labs(col = "Fish Origin", shape = "Fish Origin") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_color_manual("Fish Origin", values = c("#705653","#536D70"), 
                     labels = c("Wild-spawned", "Hatchery-reared")) + 
  scale_shape_manual("Fish Origin", values = c("circle", "triangle"),
                     labels = c("Wild-spawned", "Hatchery-reared"))+
  ylab("Average Growth Temperature (C)")  -> Lake_Temp_Graph

#ggsave(Lake_Temp_Graph, file = "Graphics/Lake_Temp_Graph.png", width = 7, height = 4, dpi = 500)


## LMER actually is AGT ~ TDO5 + TL + (1|Water) ## So I think I need to replot the fit line
mod.data = df_cross %>% 
  filter(afrp_abbrev %nin% c("LML", "FBL", "HAL")) %>% 
  mutate(lake_type = factor(lake_type, levels = c("Constrained", "Squeezed", "Overheated"))) %>%
  mutate(tdo5_avg = scale(tdo5_avg), 
         FISH_LENGTH_mm = scale(FISH_LENGTH_mm))
model.lmer = lmerTest::lmer(CalcWaterTemp ~ tdo5_avg + FISH_LENGTH_mm + (1|lake_name), data = mod.data)
summary(model.lmer)
pred_dat <- data.frame(
  TDO5 = seq(min(mod.data$tdo5_avg, na.rm = TRUE),
             max(mod.data$tdo5_avg, na.rm = TRUE),
             length.out = 100),
  TL = mean(mod.data$FISH_LENGTH_mm, na.rm = TRUE)
)

pred_dat <- data.frame(
  tdo5_avg = seq(
    min(mod.data$tdo5_avg, na.rm = TRUE),
    max(mod.data$tdo5_avg, na.rm = TRUE),
    length.out = 100
  ),
  FISH_LENGTH_mm = mean(mod.data$FISH_LENGTH_mm, na.rm = TRUE)
)
library(emmeans)
pred_dat <- emmeans(
  model.lmer,
  ~ tdo5_avg,
  at = list(
    tdo5_avg = seq(
      min(mod.data$tdo5_avg, na.rm = TRUE),
      max(mod.data$tdo5_avg, na.rm = TRUE),
      length.out = 100
    ),
    FISH_LENGTH_mm = mean(mod.data$FISH_LENGTH_mm, na.rm = TRUE)
  )
) %>%
  as.data.frame()

## Lake by TDO5



ggplot() + 
 
  geom_ribbon(aes(x = tdo5_avg, y = emmean, ymin =lower.CL, ymax = upper.CL), 
              data = pred_dat, alpha = .1) + 
  geom_line(aes(x = pred_dat$tdo5_avg, y = pred_dat$emmean), lty = "dashed") +
  geom_point(aes(x = tdo5_avg, y = CalcWaterTemp, 
                 col = lake_type, shape = lake_type), size = 3,
             data = mod.data) +
  xlab("TDO5 (C)") + ylab("Average Growing Temperature (C)") +
  scale_color_manual("Lake Category", values = viridis(n = 3))  + #c("#839755", "#976955","#558397")
  scale_shape_manual("Lake Category", values = c("circle", "triangle", "square")) +
  theme_minimal(base_size = 11) + 
  theme(legend.position = "top")
ggsave(file = "Graphics/Figure_4A.pdf", width = 3.5, height = 3, dpi = 500)

## Trying to derive the cumulative frequency graph 

df_cross %>% 
  filter(lake_type != "Buffered") %>% 
  filter(CalcWaterTemp <= 15) %>%
  dim() ## 42 / 62 fish had calcwater temp below or = to 15

df_cross %>% 
  group_by(lake_type) %>% 
  summarize(mean_temp = mean(CalcWaterTemp))
df_cross %>%
  filter(lake_type != "Buffered") %>%
  filter(!is.na(CalcWaterTemp)) %>%
  group_by(lake_type) %>%
  arrange(CalcWaterTemp) %>%
  mutate(
    cumulative_frequency = cume_dist(CalcWaterTemp)
  ) %>% 
  select(lake_type, CalcWaterTemp, cumulative_frequency) %>%
  filter(cumulative_frequency > .75) %>%
  slice_head()

ggplot( aes(x = CalcWaterTemp, color = lake_type), data = df_cross %>% 
          filter(lake_type != "Buffered") %>% 
          mutate(lake_type = factor(lake_type, 
                                    levels = c("Constrained", "Squeezed", "Overheated")))) +
  stat_ecdf(linewidth = 1.2) +
  labs(
    x = "Lifelong Experienced Temperature (°C)",
    y = "Cumulative Frequency"
  ) +
  theme_minimal(base_size = 11) + 
  scale_color_manual("Lake Category", values = viridis(n = 3)) + 
  theme(legend.position = "top")
ggsave(file = "Graphics/Figure_4B.pdf", width = 3.5, height = 3)

df_cross %>%
  filter(afrp_abbrev == "LML") %>%
  summarize(range(CalcWaterTemp))


## Table of ST characteristics across all lakes to mimic other supplementary graph 

df_cross %>% 
  separate(FISH_N, into = c("species", "lake", "date", "gear","num")) %>%
  group_by(Lake) %>%
  mutate(gear = paste(unique(gear), collapse = ", ")) %>%
  mutate(SAMPLE_DATE = paste(unique(SAMPLE_DATE), collapse = ", ")) %>%
  mutate(total_fish = n()) %>%
  group_by(Lake, total_fish, gear, SAMPLE_DATE) %>%
  summarize(mean_length = mean(FISH_LENGTH_mm, na.rm = T),
            SD_length = sd(FISH_LENGTH_mm, na.rm = T),
            min_length = min(FISH_LENGTH_mm, na.rm = T), 
            max_length = max(FISH_LENGTH_mm, na.rm = T),
            gear = unique(gear)) %>%
  select(Lake, total_fish, mean_length, SD_length, min_length, max_length, gear, SAMPLE_DATE) -> brook.trout.supp

write.csv(brook.trout.supp, file = "Data/Tables/TableS3A_ST_Summary.csv", row.names = F)  


df_cross %>% 
  filter(afrp_abbrev == "LML") %>% 
  select(FISH_N, FINAL_AGE, FISH_LENGTH_mm, FISH_WEIGHT_g) %>% 
  write.csv(., file = "Data/Metadata and Measurement data/LML_ST.csv", row.names = F)

## ST length by lake visualization

df_cross %>% 
  ggplot(aes(y = afrp_abbrev, x = FISH_LENGTH_mm)) + 
  geom_point() 


## Kruskal -Wallist Test and posthoc comparisons for visualization


# Kruskal-Wallis test
kruskal.test(scale(CalcWaterTemp) ~ lake_name, data = df_cross %>% 
               filter(is.na(YSAMP_N)))

# Post hoc pairwise Wilcoxon tests
## None of these are showing up as significant?
pairwise = pairwise.wilcox.test(
  df_cross$CalcWaterTemp,
  df_cross$lake_name,
  p.adjust.method = "holm"
) 

pairwise$p.value %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "L1") %>% 
  pivot_longer(2:ncol(.)) %>% 
  filter(value < .05)

