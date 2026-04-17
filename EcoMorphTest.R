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


#### Weather ####

library(terra)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# 1) Find all GeoTIFFs (nested OK)
files <- list.files(
  getwd(),
  pattern = "\\.(tif|tiff)$",
  full.names = TRUE,
  ignore.case = TRUE,
  recursive = TRUE
)
stopifnot(length(files) > 0)

# 2) Read each file safely, and name layers correctly (handles multi-band TIFFs)
ras_list <- lapply(files, function(f) {
  x <- rast(f)
  base <- tools::file_path_sans_ext(basename(f))
  if (nlyr(x) == 1) {
    names(x) <- base
  } else {
    names(x) <- paste0(base, "_band", seq_len(nlyr(x)))
  }
  x
})

r <- do.call(c, ras_list)  # combine into one SpatRaster

# 3) County-wide mean per layer
m <- terra::global(r, fun = "mean", na.rm = TRUE)
m$layer <- rownames(m)
rownames(m) <- NULL

# 4) Parse var + year from layer names
df_long <- m %>%
  rename(value = mean) %>%
  mutate(
    var  = str_extract(layer, "^(pr|tmax|tmin)"),
    year = as.integer(str_extract(layer, "\\d{4}"))
  ) %>%
  filter(!is.na(var), !is.na(year)) %>%
  select(var, year, value) %>%
  group_by(year, var) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop")

# 5) Wide + compute annual average temp and rainfall
annual <- df_long %>%
  pivot_wider(names_from = var, values_from = value) %>%
  mutate(
    tavg = (tmax + tmin) / 2
  ) %>%
  select(year, tavg, pr) %>%
  arrange(year)

print(annual)

# 6) Plot
ggplot(annual, aes(year, tavg)) +
  geom_line() + geom_point() +
  labs(x = "Year", y = "Average temperature (dataset units)")

ggplot(annual, aes(year, pr)) +
  geom_line() + geom_point() +
  labs(x = "Year", y = "Precipitation (dataset units)")


library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# assumes you already have: annual with columns year, tavg, pr

# If you want prettier axis labels (and optional unit conversions)
annual2 <- annual %>%
  mutate(
    # uncomment if Kelvin
    # tavg = tavg - 273.15,
    # uncomment if precipitation is meters and you want mm
    # pr = pr * 1000
  )

# Long format for a clean faceted plot
pdat <- annual2 %>%
  pivot_longer(c(tavg, pr), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    tavg = "Average temperature",
                    pr   = "Total precipitation"
    )
  )

ggplot(pdat, aes(x = year, y = value)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    x = "Year",
    y = NULL,
    title = "Adams County climate summary",
    subtitle = "County mean from BeeSpatial rasters"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    panel.spacing = unit(1.1, "lines")
  )


#### Combine ####
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

# climate: Year, tavg, pr
climate <- annual %>%
  rename(Year = year) %>%
  select(Year, tavg, pr)

# ITD summary: Site, Year, mean_ITD_mm, se_ITD_mm, n
file_path <- "C:/Users/bunch/OneDrive - The Pennsylvania State University/EcoMorph Project/EcoMorph Testing/EcoMorph_Callibration_Data.xlsx"
ecocal <- read_excel(file_path)

mean_itd_year_site <- ecocal %>%
  mutate(
    Year = as.integer(str_extract(Sample, "\\d{4}")),
    Site = str_extract(Sample, "^[^_]+_[^_]+")
  ) %>%
  filter(!is.na(Year), !is.na(Site)) %>%
  group_by(Site, Year) %>%
  summarise(
    mean_ITD_mm = mean(`ITD (mm)`, na.rm = TRUE),
    se_ITD_mm   = sd(`ITD (mm)`, na.rm = TRUE) / sqrt(dplyr::n()),
    n           = dplyr::n(),
    .groups = "drop"
  )

# 1) Make a climate grid for EVERY site x EVERY climate year
sites <- sort(unique(mean_itd_year_site$Site))

site_year_climate <- tidyr::crossing(
  Site = sites,
  climate
)

# 2) Add ITD onto that grid (this keeps climate-only years too)
itd_climate_full <- site_year_climate %>%
  left_join(mean_itd_year_site, by = c("Site", "Year")) %>%
  arrange(Site, Year)

itd_climate_full


#### Climate Graph ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Use your full table that includes climate-only years
# (from your screenshot it looks like: Site, Year, tavg, pr, mean_ITD_mm, se_ITD_mm, n)
dat <- itd_climate_full %>%
  mutate(
    Site = str_replace_all(Site, "_", " "),
    Year = as.integer(Year),
    ymin = mean_ITD_mm - se_ITD_mm,
    ymax = mean_ITD_mm + se_ITD_mm
  )

# Keep the year axis consistent across all panels
year_breaks <- sort(unique(dat$Year))

# Build long-format data for climate-only panels (so they always show full time series)
clim_long <- dat %>%
  select(Site, Year, tavg, pr) %>%
  pivot_longer(c(tavg, pr), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    tavg = "Mean temperature",
                    pr   = "Precipitation"
    )
  )

# ITD panel data (only where ITD exists)
itd_only <- dat %>%
  filter(!is.na(mean_ITD_mm)) %>%
  transmute(
    Site, Year,
    metric = "ITD (mm)",
    value = mean_ITD_mm,
    ymin, ymax, n
  )

# Combine for faceting
plot_dat <- bind_rows(
  itd_only %>% select(Site, Year, metric, value),
  clim_long
) %>%
  mutate(metric = factor(metric, levels = c("ITD (mm)", "Mean temperature", "Precipitation")))

# Plot: climate shows all years, ITD shows only sampled years (no fake points)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

dat <- itd_climate_full %>%
  mutate(
    Site = str_replace_all(Site, "_", " "),
    Year = as.integer(Year),
    ymin = mean_ITD_mm - se_ITD_mm,
    ymax = mean_ITD_mm + se_ITD_mm
  )

year_breaks <- sort(unique(dat$Year))

# Climate in long format
clim_long <- dat %>%
  select(Site, Year, tavg, pr) %>%
  pivot_longer(c(tavg, pr), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
                    tavg = "Mean temperature",
                    pr   = "Precipitation")
  )

# ITD only where sampled
itd_only <- dat %>%
  filter(!is.na(mean_ITD_mm)) %>%
  transmute(
    Site, Year,
    metric = "ITD (mm)",
    value = mean_ITD_mm,
    ymin, ymax, n
  )

ggplot() +
  # Temperature line
  geom_line(
    data = clim_long %>% filter(metric == "Mean temperature"),
    aes(x = Year, y = value, group = 1),
    linewidth = 0.9
  ) +
  geom_point(
    data = clim_long %>% filter(metric == "Mean temperature"),
    aes(x = Year, y = value),
    size = 2.3
  ) +
  
  # Precipitation line (NEW)
  geom_line(
    data = clim_long %>% filter(metric == "Precipitation"),
    aes(x = Year, y = value, group = 1),
    linewidth = 0.9,
    linetype = "solid"
  ) +
  geom_point(
    data = clim_long %>% filter(metric == "Precipitation"),
    aes(x = Year, y = value),
    size = 2.3
  ) +
  
  # ITD with SE
  geom_line(
    data = itd_only,
    aes(x = Year, y = value, group = 1),
    linewidth = 1.0
  ) +
  geom_point(
    data = itd_only,
    aes(x = Year, y = value),
    size = 2.6
  ) +
  geom_errorbar(
    data = itd_only,
    aes(x = Year, ymin = ymin, ymax = ymax),
    width = 0.15,
    linewidth = 0.6
  ) +
  geom_text(
    data = itd_only,
    aes(x = Year, y = value, label = paste0("n=", n)),
    vjust = 1.7,
    size = 3
  ) +
  
  facet_grid(metric ~ Site, scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = year_breaks) +
  labs(
    title = "ITD and annual climate through time",
    subtitle = "Climate shown for all years; ITD shown only where sampled (± SE)",
    x = "Year",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(face = "bold", size = 13),
    strip.text.y.left = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    axis.text = element_text(size = 11),
    panel.spacing = unit(1.1, "lines")
  )

#### PRISM Data ####
library(readr)
library(dplyr)
library(lubridate)

prism <- read_csv(
  "PRISM_ppt_tmin_tmean_tmax_stable_4km_200905_202008_39.9340_-77.2553.csv",
  skip = 11,
  col_names = c("Date", "ppt_in", "tmin_F", "tmean_F", "tmax_F"),
  show_col_types = FALSE
) %>%
  mutate(
    Date  = ym(Date),
    Year  = year(Date),
    Month = month(Date)
  )

prism


library(dplyr)

summer_means <- prism %>%
  filter(Month %in% 5:8) %>%   # May–August
  group_by(Year) %>%
  summarise(
    summer_ppt_in  = mean(ppt_in,  na.rm = TRUE),
    summer_tmin_F  = mean(tmin_F,  na.rm = TRUE),
    summer_tmean_F = mean(tmean_F, na.rm = TRUE),
    summer_tmax_F  = mean(tmax_F,  na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  ) %>%
  arrange(Year)

summer_means


library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

#### 1) Summer climate (May–August means) ####
summer_means <- prism %>%
  filter(Month %in% 5:8) %>%
  group_by(Year) %>%
  summarise(
    summer_ppt_in  = mean(ppt_in,  na.rm = TRUE),
    summer_tmean_F = mean(tmean_F, na.rm = TRUE),
    .groups = "drop"
  )

#### 2) Merge with ITD summary ####
clim_itd_summer <- mean_itd_year_site %>%
  left_join(summer_means, by = "Year")

#### 3) Make full Site × Year grid (keep climate-only years) ####
sites <- unique(mean_itd_year_site$Site)

full_grid <- tidyr::crossing(
  Site = sites,
  summer_means
) %>%
  left_join(mean_itd_year_site, by = c("Site", "Year")) %>%
  arrange(Site, Year)

#### 4) Prepare plotting data ####
dat <- full_grid %>%
  mutate(
    Site = str_replace_all(Site, "_", " "),
    ymin = mean_ITD_mm - se_ITD_mm,
    ymax = mean_ITD_mm + se_ITD_mm
  )

clim_long <- dat %>%
  select(Site, Year, summer_tmean_F, summer_ppt_in) %>%
  pivot_longer(
    cols = c(summer_tmean_F, summer_ppt_in),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    summer_tmean_F = "Summer mean temperature (°F)",
                    summer_ppt_in  = "Summer mean precipitation (in)")
  )

itd_only <- dat %>%
  filter(!is.na(mean_ITD_mm)) %>%
  transmute(
    Site, Year,
    metric = "ITD (mm)",
    value = mean_ITD_mm,
    ymin, ymax, n
  )

year_breaks <- sort(unique(dat$Year))

#### 5) Plot ####
ggplot() +
  
  # Climate lines
  geom_line(
    data = clim_long,
    aes(x = Year, y = value, group = 1),
    linewidth = 0.9
  ) +
  geom_point(
    data = clim_long,
    aes(x = Year, y = value),
    size = 2.3
  ) +
  
  # ITD with SE
  geom_line(
    data = itd_only,
    aes(x = Year, y = value, group = 1),
    linewidth = 1.0
  ) +
  geom_point(
    data = itd_only,
    aes(x = Year, y = value),
    size = 2.6
  ) +
  geom_errorbar(
    data = itd_only,
    aes(x = Year, ymin = ymin, ymax = ymax),
    width = 0.15,
    linewidth = 0.6
  ) +
  geom_text(
    data = itd_only,
    aes(x = Year, y = value, label = paste0("n=", n)),
    vjust = 1.7,
    size = 3
  ) +
  
  facet_grid(metric ~ Site, scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = year_breaks) +
  labs(
    title = "Summer climate (May–Aug) and ITD through time",
    subtitle = "Climate shown for all years; ITD only where sampled (± SE)",
    x = "Year",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    strip.text.x = element_text(face = "bold"),
    strip.text.y.left = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Helper to build the full grid + plot for a chosen summer temp variable
make_itd_summer_plot <- function(temp_var = c("tmin_F", "tmax_F"),
                                 temp_label = c("Summer mean Tmin (°F)", "Summer mean Tmax (°F)")) {
  
  temp_var <- match.arg(temp_var)
  temp_label <- temp_label[match(temp_var, c("tmin_F", "tmax_F"))]
  
  # 1) Summer (May–Aug) means from PRISM
  summer_means <- prism %>%
    filter(Month %in% 5:8) %>%
    group_by(Year) %>%
    summarise(
      summer_ppt_in = mean(ppt_in, na.rm = TRUE),
      summer_temp   = mean(.data[[temp_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2) Full Site × Year climate grid and join ITD (keeps climate-only years)
  sites <- sort(unique(mean_itd_year_site$Site))
  
  full_grid <- tidyr::crossing(
    Site = sites,
    summer_means
  ) %>%
    left_join(mean_itd_year_site, by = c("Site", "Year")) %>%
    arrange(Site, Year) %>%
    mutate(
      Site = str_replace_all(Site, "_", " "),
      ymin = mean_ITD_mm - se_ITD_mm,
      ymax = mean_ITD_mm + se_ITD_mm
    )
  
  # 3) Climate long format
  clim_long <- full_grid %>%
    select(Site, Year, summer_temp, summer_ppt_in) %>%
    pivot_longer(c(summer_temp, summer_ppt_in), names_to = "metric", values_to = "value") %>%
    mutate(
      metric = recode(metric,
                      summer_temp   = temp_label,
                      summer_ppt_in = "Summer mean precipitation (in)")
    )
  
  # 4) ITD only where sampled
  itd_only <- full_grid %>%
    filter(!is.na(mean_ITD_mm)) %>%
    transmute(
      Site, Year,
      metric = "ITD (mm)",
      value = mean_ITD_mm,
      ymin, ymax, n
    )
  
  year_breaks <- sort(unique(full_grid$Year))
  
  # 5) Plot
  ggplot() +
    geom_line(
      data = clim_long,
      aes(x = Year, y = value, group = 1),
      linewidth = 0.9
    ) +
    geom_point(
      data = clim_long,
      aes(x = Year, y = value),
      size = 2.3
    ) +
    geom_line(
      data = itd_only,
      aes(x = Year, y = value, group = 1),
      linewidth = 1.0
    ) +
    geom_point(
      data = itd_only,
      aes(x = Year, y = value),
      size = 2.6
    ) +
    geom_errorbar(
      data = itd_only,
      aes(x = Year, ymin = ymin, ymax = ymax),
      width = 0.15,
      linewidth = 0.6
    ) +
    geom_text(
      data = itd_only,
      aes(x = Year, y = value, label = paste0("n=", n)),
      vjust = 1.7,
      size = 3
    ) +
    facet_grid(metric ~ Site, scales = "free_y", switch = "y") +
    scale_x_continuous(breaks = year_breaks) +
    labs(
      title = paste0("Summer climate (May–Aug) and ITD through time: ", gsub("Summer mean ", "", temp_label)),
      subtitle = "Climate shown for all years; ITD only where sampled (± SE)",
      x = "Year",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      strip.text.x = element_text(face = "bold"),
      strip.text.y.left = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# Plot 1: Tmin
p_tmin <- make_itd_summer_plot("tmin_F", c("Summer mean Tmin (°F)", "Summer mean Tmax (°F)"))
print(p_tmin)

# Plot 2: Tmax
p_tmax <- make_itd_summer_plot("tmax_F", c("Summer mean Tmin (°F)", "Summer mean Tmax (°F)"))
print(p_tmax)


