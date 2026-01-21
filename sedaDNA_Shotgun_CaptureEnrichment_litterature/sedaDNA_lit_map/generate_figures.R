# ===============================================================
# Ancient Metagenomic Publication Data Visualization Pipeline
# ===============================================================
# Description:
# This script processes and visualizes geospatial and temporal 
# trends in ancient metagenomic studies using data from a Google Sheet.
# It includes:
#   - Setup and data import
#   - Spatial reprojecting and site mapping
#   - Yearly and cumulative publication bar plots
#   - Analysis of reference databases used
#
# Output:
# Several figures saved in .png or .pdf format for publication or presentation
# ===============================================================


# -----------------------------
# Load Required Libraries (Install if Missing)
# -----------------------------
# These packages are needed for geospatial handling, plotting, 
# reading from Google Sheets, and data wrangling.

required_packages <- c(
  "sf", "ggplot2", "ggrepel", "rnaturalearth", "rnaturalearthdata",
  "readr", "readxl", "tidyverse", "googlesheets4", "dplyr", 
  "scales", "viridis"
)

# Install any missing packages
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load all packages
lapply(required_packages, library, character.only = TRUE)


# -----------------------------
# Authenticate and Read Data
# -----------------------------
# Authenticate with Google Sheets and read the dataset ("Table1" sheet)
gs4_auth()
coordinates <- read_sheet(
  "https://docs.google.com/spreadsheets/d/13cmBUi4cigUaTKtQeFLFvS0gXT8AeWxWKzHv2UcOBCI/edit?gid=0#gid=0", 
  sheet = "Eukaryotes"
)


# -----------------------------
# Data Inspection & Cleaning
# -----------------------------
# Check column names and distinct TargetGroup values.
# Filter out NA Latitude values and "Microorganisms | NA" group.
colnames(coordinates)
unique(coordinates$TargetGroup)
coordinates_noNA <- coordinates %>%
  filter(Latitude != "NA", TargetGroup != "Microorganisms | NA")


# -----------------------------
# Geospatial Preparation
# -----------------------------
# Load and reproject world map to Robinson projection for better visualization.
# Convert coordinates to sf object and project to match world map.
world <- ne_countries(scale = "medium", returnclass = "sf")
world_proj <- st_transform(world, crs = st_crs("+proj=robin"))

coordinates_sf <- st_as_sf(coordinates_noNA, coords = c("Longitude", "Latitude"), crs = 4326)
coordinates_proj <- st_transform(coordinates_sf, crs = st_crs("+proj=robin"))
coordinates_proj <- cbind(coordinates_noNA, st_coordinates(coordinates_proj))


# -----------------------------
# Filter & Plot Spatial Data
# -----------------------------
# Filter dataset for valid entries and plot the world map with study sites.
# Label sites and color by MolecularMethod shape.
coordinates_proj_SG_TE <- coordinates_proj %>%
  mutate(Lab = as.factor(unlist(Lab))) %>%
  filter(
    year_published != "NA",
    TargetGroup != "Microorganisms",
    TargetTaxa != "Microorganisms",
    TargetTaxa != "Prokaryotes"
  )

ggplot() +
  geom_sf(data = world_proj, fill = "lightgrey", color = "black") +
  geom_point(data = coordinates_proj_SG_TE, 
             aes(x = X, y = Y, shape = MolecularMethod), 
             color = "black", size = 3) +
  geom_text_repel(data = coordinates_proj_SG_TE, 
                  aes(x = X, y = Y, label = SiteName),
                  color = "black", size = 3, fontface = "bold", box.padding = 0.3) +
  coord_sf(crs = st_crs("+proj=robin")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Ancient Metagenomic Study Sites",
    x = "Longitude",
    y = "Latitude",
    shape = "Molecular Method"
  )

ggsave("../../figures/SG_TE_map_method.png", width = 10, height = 7, dpi = 300)

colnames(coordinates)

# -----------------------------
# Yearly Publication Barplot
# -----------------------------
# Count publications per year (after filtering) and plot as bar chart.
coordinates %>%
  filter(
    year_published != "NA",
    TargetGroup != "Microorganisms",
    TargetTaxa != "Microorganisms"
  ) %>%
  select(year_published, MolecularMethod, DOI) %>%
  unique() %>%
  group_by(year_published, MolecularMethod) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(x = year_published, y = n)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Ancient Metagenome Publication Counts by Year",
    x = "Year Published",
    y = "Number of Publications",
    fill = "Molecular Method"
  )

ggsave("../../figures/barplot_no_publications.png", width = 6, height = 4, dpi = 300)


# -----------------------------
# Cumulative Publication Trend
# -----------------------------
# Compute cumulative counts over time for each method.
coordinates %>%
  filter(
    year_published != "NA",
    TargetGroup != "Microorganisms",
    TargetTaxa != "Microorganisms"
  ) %>%
  select(year_published, MolecularMethod, DOI) %>%
  unique() %>%
  group_by(year_published, MolecularMethod) %>%
  count() %>%
  ungroup() %>%
  arrange(year_published) %>%
  group_by(MolecularMethod) %>%
  mutate(cumulative_n = cumsum(n)) %>%
  ggplot(aes(x = year_published, y = cumulative_n, fill = MolecularMethod)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  theme_minimal(base_size = 14) +
  theme_test() +
  labs(
    title = "Cumulative Publications Over Time by Molecular Method",
    x = "Year Published",
    y = "Cumulative Number of Publications",
    fill = "Molecular Method"
  )

ggsave("../../figures/barplot_cumsum_no_publications_methods.png", width = 6, height = 4, dpi = 300)



# -----------------------------
# Yearly Stacked Barplot by Method
# -----------------------------
# Stacked bars for number of publications each year, broken down by method.
coordinates %>%
  filter(
    year_published != "NA",
    TargetGroup != "Microorganisms",
    TargetTaxa != "Microorganisms"
  ) %>%
  select(year_published, MolecularMethod, DOI) %>%
  unique() %>%
  group_by(year_published, MolecularMethod) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(x = year_published, y = n, fill = MolecularMethod)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  labs(
    title = "Ancient Metagenome Publication Counts by Year and Molecular Method",
    x = "Year Published",
    y = "Number of Publications",
    fill = "Molecular Method"
  )

ggsave("../../figures/barplot_no_publications_methods.png", width = 8, height = 5, dpi = 300)


# -----------------------------
# Publications by Reference Database
# -----------------------------
# Count number of publications using each reference database and visualize.
coordinates %>%
  filter(
    year_published != "NA",
    TargetGroup != "Microorganisms",
    TargetTaxa != "Microorganisms",
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(reference_target, DOI) %>%
  unique() %>%
  group_by(reference_target) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(reference_target, n)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(
    title = "Publications by Reference Databases",
    x = "Reference Target",
    y = "Number of Publications"
  )

ggsave("../../figures/barplot_N_databases.png", width = 8, height = 5, dpi = 300)

# -----------------------------
# Publications by Mappers
# -----------------------------
# Count number of publications using each reference database and visualize.
coordinates %>%
  filter(
    !is.na(year_published),
    TargetGroup != "Microorganisms",
    TargetTaxa != "Microorganisms",
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(mapper, DOI) %>%
  unique() %>%
  group_by(mapper) %>%
  count() %>%
  ungroup() %>%
  filter(!is.na(mapper), !is.na(n)) %>%              # Remove NAs
  mutate(mapper = as.character(mapper)) %>%          # Ensure mapper is character
  ggplot(aes(mapper, n)) +
  geom_bar(stat = "identity", position = "dodge", na.rm = TRUE) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(
    title = "Publications by Mapper",
    x = "Mapper",
    y = "Number of Publications"
  )

ggsave("../../figures/barplot_N_mappers.png", width = 8, height = 5, dpi = 300)



# -----------------------------
# Color by SampleType and shape by TargetGroup
# -----------------------------
# Filter dataset for valid entries and plot the world map with study sites.
# Label sites and color by MolecularMethod shape.
# Load libraries
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

# Prepare data
coordinates_proj_SG_TE <- coordinates_proj %>%
  mutate(Lab = as.factor(unlist(Lab))) %>%
  filter(
    year_published != "NA",
    TargetGroup != "Microorganisms",
    TargetTaxa != "Microorganisms",
    TargetTaxa != "Prokaryotes"
  )

# Refined Okabe–Ito style palette (13 distinct, colorblind-safe tones)
okabe_ito_mod <- c(
  "#E69F00", # orange-yellow
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#8A2BE2", # violet (replaces lemon yellow)
  "#0072B2", # blue
  "#D55E00", # vermilion
  "#CC79A7", # reddish purple
  "#999999", # grey
  "#A52A2A", # brown
  "#00CED1", # turquoise
  "#FFD700", # golden yellow
  "#228B22", # forest green
  "#DA70D6"  # orchid pink
)

# Main plot
ggplot() +
  # World map background
  geom_sf(data = world_proj, fill = "lightgrey", color = "black") +
  
  # Sample points
  geom_point(
    data = coordinates_proj_SG_TE,
    aes(x = X, y = Y, shape = TargetGroup, color = SampleType),
    size = 3
  ) +
  
  # Site labels
#  geom_text_repel(
#    data = coordinates_proj_SG_TE,
#    aes(x = X, y = Y, label = SiteName),
#    color = "black", size = 3, fontface = "bold", box.padding = 0.3
#  ) +
  
  # Apply refined color palette
  scale_color_manual(values = okabe_ito_mod) +
  
  # Coordinate system
  coord_sf(crs = st_crs("+proj=robin")) +
  
  # Clean minimal theme
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  
  # Labels and titles
  labs(
    title = "Ancient Metagenomic Study Sites",
    subtitle = "Global distribution of sampling locations by Sample Type and Studied Organism",
    x = "Longitude",
    y = "Latitude",
    shape = "Target Organism",
    color = "Sample Type"
  )

ggsave("../../figures/SG_TE_map_SampleType_Target.png", width = 10, height = 7, dpi = 300)

colnames(coordinates)




#### END OF SCRIPT

