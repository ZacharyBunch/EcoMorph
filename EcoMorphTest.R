#### Libraries ####
library(readxl)

#### File read ####

file_path <- "C:/Users/bunch/OneDrive - The Pennsylvania State University/EcoMorph Project/EcoMorph Testing/EcoMorph_Callibration_Data.xlsx"

ecocal <- read_excel(file_path)


#### Data organization ####

library(dplyr)
library(stringr)
library(ggplot2)

mean_itd_year_site <- ecocal %>%
  mutate(
    Year = str_extract(Sample, "\\d{4}") %>% as.integer(),
    Site = str_extract(Sample, "^[^_]+_[^_]+")
  ) %>%
  filter(!is.na(Year), !is.na(Site)) %>%
  group_by(Site, Year) %>%
  summarise(
    mean_ITD_mm = mean(`ITD (mm)`, na.rm = TRUE),
    sd_ITD_mm   = sd(`ITD (mm)`, na.rm = TRUE),
    n           = n(),
    .groups = "drop"
  ) %>%
  arrange(Site, Year)

mean_itd_year_site


#### Bar Graph ####

library(ggplot2)

ggplot(mean_itd_year_site,
       aes(x = factor(Year), y = mean_ITD_mm)) +
  geom_col(fill = "grey60", width = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_ITD_mm - sd_ITD_mm / sqrt(n),
      ymax = mean_ITD_mm + sd_ITD_mm / sqrt(n)
    ),
    width = 0.2
  ) +
  geom_text(
    aes(label = paste0("n = ", n)),
    vjust = -0.6,
    size = 3.5
  ) +
  facet_wrap(~ Site, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Mean ITD (mm)"
  ) +
  coord_cartesian(clip = "off")


#### Stats ####
library(dplyr)
library(stringr)
library(emmeans)

ecocal2 <- ecocal %>%
  mutate(
    Year = str_extract(Sample, "\\d{4}") %>% factor(),
    Site = str_extract(Sample, "^[^_]+_[^_]+")
  ) %>%
  filter(!is.na(Year), !is.na(Site))

pairwise_years <- ecocal2 %>%
  group_by(Site) %>%
  group_modify(~{
    m <- lm(`ITD (mm)` ~ Year, data = .x)
    
    as.data.frame(
      emmeans(m, pairwise ~ Year, adjust = "tukey")$contrasts
    )
  }) %>%
  ungroup()

pairwise_years

