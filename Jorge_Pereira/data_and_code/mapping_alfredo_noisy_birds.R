# =============================================================================
# Making your meta-analysis spatial - practical
# AEET course, 22 October 2026
# Jorge Pereira (University of Coimbra)
#
# Data: effects of anthropogenic noise on birds (meta-analysis kindly shared by Alfredo)
# Each row is one effect size. Effect sizes come from studies, and one study
# can give many effect sizes.
#
# What we will do in this (apparently long) but simple script.
# We have little time to practice coding today, so the session is mostly
# about running the code step by step (Ctrl + Enter on Windows, or Cmd + Enter on Mac)
# and trying to understand the results.
#
#   Part 1 - put the studies on a map and look at where the evidence comes from
#   Part 2 - use the coordinates to build moderators (latitude, longitude, and other
#            variables that you can calculate for your data)
#   Part 3 - test the influence of those moderators on the effect size
# =============================================================================


# 0. Start clean ---------------------------------------------------------------

# Remove everything from the environment, so we start from zero
rm(list = ls())

# Set the working directory to the folder where this script is saved.
# This only works in RStudio and only if the script is saved in the same folder as the dataset.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()


# 0b. Packages -----------------------------------------------------------------

# Install once (remove the # and run). After that, keep it commented,
# otherwise R reinstalls everything every time you run the script
 install.packages(c("sf", "ggplot2", "dplyr", "readr", "patchwork",
                    "rnaturalearth", "rnaturalearthdata", "glmmTMB", "emmeans",
                    "rstudioapi"))

library(readr)          # read csv files
library(dplyr)          # work with tables (filter, group_by, summarise...)
library(ggplot2)        # plots and maps
library(patchwork)      # put several ggplots side by side
library(sf)             # spatial data: points, polygons, projections
library(rnaturalearth)  # world map (countries) ready to use
library(glmmTMB)        # mixed models
library(emmeans)        # compare multiple groups after a model

theme_set(theme_bw(base_size = 12))   # same clean look for all plots


# =============================================================================
# PART 1. MAPPING THE EVIDENCE
# =============================================================================

# 1. Read and look at the data -------------------------------------------------

# Study_ID is a code, not a number, so we read it as text
es_raw <- read_csv("map_dataset_noise_and_birds.csv", col_types = cols(Study_ID = col_character()))

glimpse(es_raw)
summary(es_raw$yi.SMD.H.signed)

# The original names are long, so we give shorter ones
es <- es_raw |>
  rename(yi      = yi.SMD.H.signed,
         lat     = Latitude_N,
         lon     = Longitude_E,
         species = species.updated)

# How much evidence do we have?
n_effects <- nrow(es)
n_studies <- n_distinct(es$Study_ID)
n_species <- n_distinct(es$species)
n_sites   <- nrow(distinct(es, lat, lon))

n_effects
n_studies
n_species
n_sites

# Number of effect sizes per study
es_per_study <- es |>
  count(Study_ID, sort = TRUE)
es_per_study


# 2. First map: where are the studies? -----------------------------------------

# Turn the table into a spatial object.
# coords = which columns are x and y (longitude first!)
# crs = 4326 means plain longitude/latitude (WGS84, the GPS system)
# remove = FALSE keeps lon and lat as normal columns too
es_sf <- st_as_sf(es, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Country polygons from Natural Earth (scale 50 = medium detail)
world <- ne_countries(scale = 50, returnclass = "sf")

# A quick map, without any projection on purpose (compare with the next one)
map_first <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.1) +
  geom_sf(data = es_sf, colour = "red", size = 2) +
  labs(subtitle = "Where are the studies?")
map_first

# Many effect sizes share the same coordinates, so points are on top of each
# other. We make a table with one row per site (unique lat/lon)
sites <- es |>
  group_by(lat, lon) |>
  summarise(n_es      = n(),                    # effect sizes at the site
            n_studies = n_distinct(Study_ID),   # studies at the site
            mean_yi   = mean(yi),
            .groups   = "drop")

# Give each site a number, we will use it later to join tables
sites$site_id <- 1:nrow(sites)

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Now with a proper projection. Equal Earth keeps the areas right, so Africa
# is not smaller than Greenland. EPSG:8857 is the code
crs_world <- "EPSG:8857"

map_sites <- ggplot() +
  geom_sf(data = world, fill = "grey92", colour = "grey70", linewidth = 0.1) +
  geom_sf(data = sites_sf, aes(size = n_es), shape = 21,
          fill = "#2C7BB6", colour = "white", alpha = 0.8) +
  scale_size_area(max_size = 9, name = "Effect sizes") +
  coord_sf(crs = crs_world) +
  labs(subtitle = paste(n_effects, "effect sizes,", n_studies, "studies,",
                        nrow(sites), "sites"))
map_sites


# 3. What does a coordinate represent? -----------------------------------------

# Step 1: coordinates as text
lat_text <- as.character(sites$lat)
lon_text <- as.character(sites$lon)

# Step 2: keep only what comes after the dot and count the characters.
# sub(".*\\.", "", x) means "delete everything up to the dot".
# If there is no dot (e.g. "30") there are 0 decimals
sites$dec_lat <- ifelse(grepl(".", lat_text, fixed = TRUE),
                        nchar(sub(".*\\.", "", lat_text)), 0)
sites$dec_lon <- ifelse(grepl(".", lon_text, fixed = TRUE),
                        nchar(sub(".*\\.", "", lon_text)), 0)

# Step 3: the less precise of the two, and the approximate error in km
sites$decimals     <- pmin(sites$dec_lat, sites$dec_lon)
sites$precision_km <- 111 / 10^sites$decimals

# Sites ordered from the least to the most precise
sites |>
  select(lat, lon, decimals, precision_km) |>
  arrange(decimals)

table(sites$decimals)


# 4. Countries and continents --------------------------------------------------

# We only need the name and the continent of each country
world_simple <- world |>
  select(country = name, continent)

# Shorter names for the plots
world_simple$country[world_simple$country == "United States of America"] <- "USA"
world_simple$country[world_simple$country == "United Kingdom"] <- "UK"

# In which country is each site?
# st_nearest_feature gives the closest country. We use it instead of a normal
# "point inside polygon" join because coastal sites can fall in the sea
nearest <- st_nearest_feature(sites_sf, world_simple)

sites$country   <- world_simple$country[nearest]
sites$continent <- world_simple$continent[nearest]

# Natural Earth puts the whole of Russia in Europe. Our Russian site is in
# Siberia (east of the Urals, ~60 E), so we move it to Asia by hand
sites$continent[sites$country == "Russia" & sites$lon > 60] <- "Asia"

# Update the spatial version with the new columns
sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Bring site_id, country and continent back to every effect size
es <- es |>
  left_join(select(sites, lat, lon, site_id, country, continent),
            by = c("lat", "lon"))


# 4a. How many studies per country and continent?

studies_country <- es |>
  group_by(country, continent) |>
  summarise(n_studies = n_distinct(Study_ID),
            n_es      = n(),
            .groups   = "drop") |>
  arrange(desc(n_studies))
studies_country

# reorder() sorts the bars from the largest to the smallest
bar_studies_country <- ggplot(studies_country,
                              aes(x = reorder(country, -n_studies),
                                  y = n_studies, fill = continent)) +
  geom_col() +
  labs(x = NULL, y = "Number of studies", fill = "Continent",
       title = "Studies per country") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

studies_continent <- es |>
  group_by(continent) |>
  summarise(n_studies = n_distinct(Study_ID), .groups = "drop")
studies_continent

bar_studies_continent <- ggplot(studies_continent,
                                aes(x = reorder(continent, -n_studies),
                                    y = n_studies, fill = continent)) +
  geom_col(show.legend = FALSE) +
  labs(x = NULL, y = "Number of studies", title = "Studies per continent")


# 4b. Mean effect size per country and continent

effect_country <- es |>
  group_by(country, continent) |>
  summarise(mean_yi = mean(yi),
            n_es    = n(),
            .groups = "drop")
effect_country

# The label goes above positive bars and below negative bars
effect_country$label_y <- ifelse(effect_country$mean_yi >= 0,
                                 effect_country$mean_yi + 0.03,
                                 effect_country$mean_yi - 0.03)

bar_effect_country <- ggplot(effect_country,
                             aes(x = reorder(country, mean_yi), y = mean_yi,
                                 fill = continent)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  geom_text(aes(y = label_y, label = paste0("n=", n_es)), size = 3) +
  labs(x = NULL, y = "Mean effect size (SMD)", fill = "Continent",
       title = "Mean effect per country") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

effect_continent <- es |>
  group_by(continent) |>
  summarise(mean_yi = mean(yi),
            n_es    = n(),
            .groups = "drop")
effect_continent

effect_continent$label_y <- ifelse(effect_continent$mean_yi >= 0,
                                   effect_continent$mean_yi + 0.03,
                                   effect_continent$mean_yi - 0.03)

bar_effect_continent <- ggplot(effect_continent,
                               aes(x = reorder(continent, mean_yi), y = mean_yi,
                                   fill = continent)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_text(aes(y = label_y, label = paste0("n=", n_es)), size = 3) +
  labs(x = NULL, y = "Mean effect size (SMD)", title = "Mean effect per continent")

# The four bar plots together (patchwork: | side by side, / one below the other)
# Note: n in the bar labels is the number of effect sizes, not studies
fig_bars <- (bar_studies_country | bar_studies_continent) /
            (bar_effect_country  | bar_effect_continent) +
  plot_layout(guides = "collect", widths = c(2, 1))
fig_bars


# 4c. Map of the number of studies per country (choropleth)

world_studies <- world_simple |>
  left_join(select(studies_country, country, n_studies), by = "country")

# Countries without studies get NA and are drawn in very light grey
map_choropleth <- ggplot(world_studies) +
  geom_sf(aes(fill = n_studies), colour = "grey70", linewidth = 0.1) +
  scale_fill_viridis_c(option = "mako", direction = -1,
                       na.value = "grey95", name = "Studies") +
  coord_sf(crs = crs_world) +
  labs(subtitle = "Number of studies per country")
map_choropleth

# Question: where is there NO evidence?
# Can we say anything about tropical birds? About Africa or Asia?


# 5. Map the effect sizes: one dot per study -----------------------------------

# The study is our unit of replication, so we want one dot per study.
# A few studies worked at more than one site. We keep one dot per study and
# site, because averaging their coordinates would put a dot where nobody worked
studies_sites <- es |>
  group_by(Study_ID, lat, lon, continent) |>
  summarise(n_es    = n(),
            mean_yi = mean(yi),
            .groups = "drop")

studies_sites_sf <- st_as_sf(studies_sites, coords = c("lon", "lat"),
                             crs = 4326, remove = FALSE)

# We fix the colour and size limits here, so all the maps below use the
# same scale and can be compared
lim_yi <- max(abs(studies_sites$mean_yi))
max_es <- max(studies_sites$n_es)

# Several studies were done at the same site and their dots overlap.
# st_jitter moves each dot a little (up to 60 km) so we can see all of them.
# set.seed makes the random move the same every time we run the script
set.seed(1)
studies_world <- studies_sites_sf |>
  st_transform(crs_world) |>
  st_jitter(amount = 60000)          # metres, because Equal Earth is in metres

# Colours: red = negative, blue = positive, white = around zero
map_studies <- ggplot() +
  geom_sf(data = world, fill = "grey92", colour = "grey70", linewidth = 0.1) +
  geom_sf(data = studies_world, aes(size = n_es, fill = mean_yi),
          shape = 21, colour = "grey20", stroke = 0.3, alpha = 0.9) +
  scale_fill_distiller(palette = "RdBu", direction = 1,
                       limits = c(-lim_yi, lim_yi),
                       name = "Mean effect\nper study (SMD)") +
  scale_size_area(max_size = 8, limits = c(1, max_es),
                  name = "Effect sizes\nper study") +
  coord_sf(crs = crs_world) +
  labs(subtitle = paste("One dot per study:", n_studies, "studies,",
                        n_effects, "effect sizes"))
map_studies


# 6. Regional maps with local projections --------------------------------------

# At regional scale we use projections made for each region
crs_namerica <- "ESRI:102008"   # North America Albers Equal Area
crs_europe   <- "EPSG:3035"     # Europe Lambert Azimuthal Equal Area

namerica <- world |> filter(continent == "North America")
europe   <- world |> filter(continent == "Europe")

# Smaller jitter (20 km) because we are zoomed in
set.seed(1)
studies_namerica <- studies_sites_sf |>
  filter(continent == "North America") |>
  st_transform(crs_namerica) |>
  st_jitter(amount = 20000)

studies_europe <- studies_sites_sf |>
  filter(continent == "Europe") |>
  st_transform(crs_europe) |>
  st_jitter(amount = 20000)

# xlim and ylim cut the map to the area we want (in metres of each projection)
map_na <- ggplot() +
  geom_sf(data = namerica, fill = "grey92", colour = "grey70", linewidth = 0.1) +
  geom_sf(data = studies_namerica, aes(fill = mean_yi, size = n_es),
          shape = 21, colour = "grey20", alpha = 0.9) +
  scale_fill_distiller(palette = "RdBu", direction = 1,
                       limits = c(-lim_yi, lim_yi), name = "Mean effect\nper study") +
  scale_size_area(max_size = 8, limits = c(1, max_es),
                  name = "Effect sizes\nper study") +
  coord_sf(crs = crs_namerica, xlim = c(-2.5e6, 3e6), ylim = c(-2e6, 1.8e6)) +
  labs(subtitle = "North America")

map_eu <- ggplot() +
  geom_sf(data = europe, fill = "grey92", colour = "grey70", linewidth = 0.1) +
  geom_sf(data = studies_europe, aes(fill = mean_yi, size = n_es),
          shape = 21, colour = "grey20", alpha = 0.9) +
  scale_fill_distiller(palette = "RdBu", direction = 1,
                       limits = c(-lim_yi, lim_yi), name = "Mean effect\nper study") +
  scale_size_area(max_size = 8, limits = c(1, max_es),
                  name = "Effect sizes\nper study") +
  coord_sf(crs = crs_europe, xlim = c(2.6e6, 5.6e6), ylim = c(1.5e6, 4.4e6)) +
  labs(subtitle = "Europe")

# Same scales in both maps, so patchwork can show a single legend
map_na + map_eu + plot_layout(guides = "collect")


# =============================================================================
# PART 2. BUILDING MODERATORS FROM THE COORDINATES
# =============================================================================

# 7. Geographic moderators -----------------------------------------------------

# We keep the moderators in their own table, one row per site.
# Any new moderator you calculate for your sites (elevation, land cover,
# distance to roads...) can be added here as a new column. In section 8 we
# join this table to the effect sizes
site_mods <- sites |>
  select(site_id, lat, lon)

# Absolute latitude = distance from the equator, no matter the hemisphere.
# It is a proxy for climate, seasonality and day length
site_mods$abs_lat <- abs(site_mods$lat)

# Hemisphere
site_mods$hemisphere <- ifelse(site_mods$lat >= 0, "North", "South")
table(site_mods$hemisphere)

# Question: with so few sites in the South, can we test the hemisphere?


# =============================================================================
# PART 3. TESTING THE MODERATORS WITH MIXED MODELS
# =============================================================================

# 8. Why a mixed model? --------------------------------------------------------

# Effect sizes from the same study are not independent: same birds, same team,
# same methods. If we ignore this, the study with 56 effect sizes counts 56
# times and the study with 1 effect size counts once.
#
# A mixed model uses all the effect sizes but adds a random effect for study:
#   (1 | Study_ID) = each study can have its own average
# so the model knows which effect sizes belong together.
#
# Note: this file has no sampling variances (vi), so effect sizes are not
# weighted by their precision as in a full meta-analysis (metafor::rma.mv).
# For the same reason we cannot add a random effect for each effect size
# (nested in study): it would be the same as the residual.

# Add the site moderators to every effect size (lat and lon are already in es)
es_mod <- es |>
  left_join(select(site_mods, -lat, -lon), by = "site_id")


# 9. Is there an overall effect of noise? --------------------------------------

# Only an intercept (1) and the random effect of study
m0 <- glmmTMB(yi ~ 1 + (1 | Study_ID), data = es_mod)
summary(m0)

# The intercept is the overall mean effect. Its p-value tests if it differs from 0
mean_effect <- fixef(m0)$cond[1]
mean_effect

plot_overall <- ggplot(es_mod, aes(x = yi)) +
  geom_histogram(binwidth = 0.1, fill = "grey70", colour = "white") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = mean_effect, colour = "red", linewidth = 1) +
  labs(x = "Effect size (SMD)", y = "Number of effect sizes",
       title = "Overall effect of noise",
       subtitle = "Red line = mean from the mixed model; dashed line = no effect (0)")
plot_overall

# Where is the variation? Between studies, or within studies (residual)?
VarCorr(m0)

# The within-study SD (~0.48) is more than twice the between-study SD (~0.20).
# Latitude or continent only change BETWEEN sites, so they can explain only a
# small part of the variation.
# Question: what changes WITHIN a study?


# 10. Does the effect change with the moderators? ------------------------------

# For each moderator: fit the model, read the test, then plot it

# 10a. Latitude

m_lat <- glmmTMB(yi ~ abs_lat + (1 | Study_ID), data = es_mod)
summary(m_lat)       # look at the abs_lat row: estimate and p-value

# Intercept and slope of the model, to draw the line
coef_lat <- fixef(m_lat)$cond
coef_lat

# geom_jitter moves the points a little sideways so they don't overlap
plot_lat <- ggplot(es_mod, aes(x = abs_lat, y = yi)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_jitter(width = 0.3, alpha = 0.3) +
  geom_abline(intercept = coef_lat[1], slope = coef_lat[2],
              colour = "red", linewidth = 1) +
  labs(x = "Absolute latitude (degrees)", y = "Effect size (SMD)",
       title = "Latitude")
plot_lat

# The line goes slightly down, but p is far from 0.05.
# "No evidence of an effect" is not the same as "evidence of no effect".


# 10b. Longitude

# Same model as latitude, now with longitude
m_lon <- glmmTMB(yi ~ lon + (1 | Study_ID), data = es_mod)
summary(m_lon)

coef_lon <- fixef(m_lon)$cond

plot_lon <- ggplot(es_mod, aes(x = lon, y = yi, colour = continent)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_jitter(width = 0.3, alpha = 0.4) +
  geom_abline(intercept = coef_lon[1], slope = coef_lon[2],
              colour = "black", linewidth = 1) +
  labs(x = "Longitude (degrees)", y = "Effect size (SMD)", colour = "Continent",
       title = "Longitude")
plot_lon

# The points form one group per continent.
# Question: what does a change "per degree of longitude" mean for a bird?


# 10c. Continent

# Asia has 1 study, Oceania and South America have 2 each: be careful
m_cont <- glmmTMB(yi ~ continent + (1 | Study_ID), data = es_mod)
summary(m_cont)

# Does continent improve the model? We compare it with the model without it
# (likelihood ratio test)
anova(m0, m_cont)

# Mean effect of each continent from the model
emm_cont <- emmeans(m_cont, ~ continent)
emm_cont

# All pairs of continents. With 5 continents there are 10 comparisons, so the
# p-values are corrected (Tukey) to avoid false positives
pairs(emm_cont)


n_per_continent <- es_mod |>
  group_by(continent) |>
  summarise(n_studies = n_distinct(Study_ID), .groups = "drop")
n_per_continent

# Height where we write the labels: a bit above the highest point
label_y <- max(es_mod$yi) + 0.3

plot_continent <- ggplot(es_mod, aes(x = continent, y = yi)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +   # no outliers: points are shown below
  geom_jitter(width = 0.15, alpha = 0.3) +
  geom_text(data = n_per_continent,
            aes(x = continent, y = label_y, label = paste0("n = ", n_studies)),
            size = 3.5) +
  labs(x = NULL, y = "Effect size (SMD)", title = "Continent",
       subtitle = "n = number of studies")
plot_continent


# 10d. Species

# Are some species more affected by noise? First, how many studies per species?
sp_summary <- es_mod |>
  group_by(species) |>
  summarise(n_studies = n_distinct(Study_ID),
            n_es      = n(),
            .groups   = "drop") |>
  arrange(desc(n_studies))
sp_summary

table(sp_summary$n_studies)

# Most species come from a single study. For those, the species effect and
# the study effect are the same thing, and we cannot separate them.
# So we don't compare species one by one. We add species as a second random
# effect (crossed with study) and ask how much variation it explains
m_sp <- glmmTMB(yi ~ 1 + (1 | Study_ID) + (1 | species), data = es_mod)
summary(m_sp)

VarCorr(m_sp)       # SD between studies, between species, and residual
anova(m0, m_sp)     # is the model better with species?

# Also look at the intercept of m_sp and compare it with m0.
# Is the overall effect still significant? How robust is it?

# Plot only the species with 3 or more studies (descriptive, not a test)
sp_common <- sp_summary |>
  filter(n_studies >= 3)

es_common <- es_mod |>
  filter(species %in% sp_common$species)

label_y_sp <- max(es_common$yi) + 0.4

# reorder() sorts the species by their median effect size
plot_species <- ggplot(es_common,
                       aes(x = reorder(species, yi, FUN = mean), y = yi)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, alpha = 0.3) +
  geom_text(data = sp_common,
            aes(x = species, y = label_y_sp, label = paste0("n = ", n_studies)),
            size = 3.5) +
  coord_flip() +                         # species names are long: horizontal
  labs(x = NULL, y = "Effect size (SMD)",
       title = "Species with 3 or more studies",
       subtitle = "n = number of studies")
plot_species
